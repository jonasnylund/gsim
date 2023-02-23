#include "octtree.h"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <fstream>
#include <memory>
#include <stack>
#include <vector>

#include <omp.h>

#include "particle.h"
#include "timers.h"


namespace model {

Node::Node(
    Tree* tree,
    Node* parent,
    int depth,
    const numerical_types::ndarray& center,
    numerical_types::real width)
    : tree(tree),
      parent(parent),
      depth(depth),
      center(center),
      width(width) {
  assert(depth == (parent == nullptr ? 0: parent->depth + 1));
  this->center_of_mass.fill(0.0);
  this->total_mass = 0.0;
  this->particles.fill(nullptr);
}

Node* Node::getSubnode(bool indices[numerical_types::num_dimensions]) {
  int linear_index = Node::SubnodeIndex(indices);

  // If the node does not yet exist, create it.
  if (!this->hasChildren()) {
    this->allocateChildren();
  }

  return this->child(linear_index);
}

void Node::addMass(numerical_types::real mass, const numerical_types::ndarray& position) {
  this->total_mass += mass;
  const numerical_types::real com_scale = mass / this->total_mass;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->center_of_mass[i] *= (1.0 - com_scale);
    this->center_of_mass[i] += position[i] * com_scale;
  }
}

void Node::add(Particle* particle) {
  this->num_particles_contained++;

  // If this one of the first particle added, we are currently a leaf node.
  if (this->num_particles_contained <= this->particles.size()) {
    assert(this->num_particles_local + 1 == this->num_particles_contained);
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i] == nullptr) {
        this->particles[i] = particle;
        this->num_particles_local++;
        particle->containing_node = this;
        return;
      }
    }
    assert(false);
  }

  // If we have local particles, but no more of them will fit, we must
  // first move them to the correct node.
  bool indices[numerical_types::num_dimensions];
  if (this->num_particles_local > 0) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i] == nullptr){
        assert(false);
        break;
      }
      this->indexOf(this->particles[i]->position, indices);
      this->getSubnode(indices)->add(this->particles[i]);
      this->particles[i] = nullptr;
    }
    this->num_particles_local = 0;
  }

  // Add the current particle to the respective subnode.
  this->indexOf(particle->position, indices);
  this->getSubnode(indices)->add(particle);
}

void Node::remove(Particle* particle) {
  for (int i = 0; i < max_num_particles; i++) {
    if (this->particles[i] == particle) {
      this->particles[i] = nullptr;
      this->num_particles_local--;
      break;
    }
  }
  assert(this->num_particles_local == this->num_particles_contained);
}

void Node::gatherChildParticles() {
  assert(this->num_particles_contained <= max_num_particles);
  if (!this->hasChildren()) {
    return;
  }

  this->num_particles_contained = 0;
  this->num_particles_local = 0;
  for (int i = 0; i < num_subnodes; i++) {
    if (this->child(i)->num_particles_contained > 0) {
      if (this->child(i)->num_particles_local == 0) {
        this->child(i)->gatherChildParticles();
      }
      for (int j = 0; j < max_num_particles; j++) {
        if (this->child(i)->particles[j] == nullptr){
          continue;
        }
        this->add(this->child(i)->particles[j]);
      }
      this->child(i)->particles.fill(nullptr);
      this->child(i)->num_particles_local = 0;
      this->child(i)->num_particles_contained = 0;
    }
  }
}

void Node::allocateChildren() {
  assert(!this->hasChildren());
  const numerical_types::real half_width = this->width / 2;

  this->children.reserve(num_subnodes);
  for (int i = 0; i < num_subnodes; i++) {
    numerical_types::ndarray center;
    for (int j = 0; j < numerical_types::num_dimensions; j++) {
      int index = (i >> j) & 0x1;
      center[j] = this->center[j] - half_width + this->width * index;
    }
    this->children.emplace_back(this->tree, this, this->depth + 1, center, half_width);
  }
  assert(this->children.size() == num_subnodes);
}

void Node::aggregateQuantities() {
  this->dirty = false;
  this->total_mass = 0.0;

  if (this->hasParticles()) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i] != nullptr) {
        this->addMass(this->particles[i]->mass,
                      this->particles[i]->position);
      }
    }
  }
  else if (this->hasChildren()) {
    for (int i = 0; i < num_subnodes; i++) {
      if (this->child(i)->num_particles_contained > 0){
        this->child(i)->aggregateQuantities();
        this->addMass(this->child(i)->total_mass,
                      this->child(i)->center_of_mass);
      }
    }
  }
}

void Node::clear() {
  if (this->hasParticles()) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i] != nullptr) {
        this->particles[i]->containing_node = nullptr;
      }
      this->particles[i] = nullptr;
    }
  }

  this->num_particles_local = 0;
  this->num_particles_contained = 0;

  if (this->hasChildren()) {
    for (int i = 0; i < num_subnodes; i++) {
      this->child(i)->clear();
    }
  }
}

bool Node::contains(const numerical_types::ndarray& point) const {
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    if (point[i] < this->center[i] - width)
      return false;
    if (point[i] > this->center[i] + width)
      return false;
  }
  return true;
}

void Node::prune() {
  if (this->parent != nullptr &&
      this->parent->num_particles_contained == 0) {
    this->children.clear();
  }
  else if (this->hasChildren()) {
    for (int i = 0; i < num_subnodes; i++) {
      this->child(i)->prune();
    }
  }
}

int Node::computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const {
  int num_calculations = 0;

  // If this node has particles, compute the force against them.
  if (this->hasParticles()) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i] == nullptr || this->particles[i] == &particle) {
        continue;
      }
      // Compute the accelleration against this particle.
      this->computeAccelleration(
        this->particles[i]->mass,
        this->particles[i]->position,
        particle,
        theta,
        epsilon,
        result);
      num_calculations++;
    }
    return num_calculations;
  }
  // Check if we should go deeper into the tree to compute the force.
  else if (this->hasChildren()) {
    // Calculate the distance between the the object and the cell.
    numerical_types::real distance_sq = 0;
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      const numerical_types::real d = this->center[i] - particle.position[i];
      distance_sq += d * d;
    }

    // The closest a particle in this node can be to the given particle is
    // the distance between the given particle and the node's center, subtract
    // the hypothenuse from the node center to its corners. In N dimensions,
    // this is sqrt(N * width ^ 2).
    // We keep all values squared though to avoid computing the sqrt.
    const numerical_types::real width_sq = this->width * this->width;
    distance_sq -= numerical_types::num_dimensions * width_sq;
    // If we have subnodes and the particle is close.
    if (theta * distance_sq < width_sq) {
      for (int i = 0; i < num_subnodes; i++) {
        if (this->constChild(i)->num_particles_contained > 0)
          num_calculations += this->constChild(i)->computeAccelleration(
            particle, theta, epsilon, result);
      }
      return num_calculations;
    }
  }
    // Compute the force against the local combined mass and position.
  this->computeAccelleration(this->total_mass, this->center_of_mass, particle, theta, epsilon, result);
  
  return num_calculations++;
}

void Node::computeAccelleration(
    numerical_types::real mass,
    const numerical_types::ndarray& position,
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const {
  numerical_types::ndarray distance_array;
  numerical_types::real distance_sq = epsilon;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    distance_array[i] = particle.position[i] - position[i];
    distance_sq += distance_array[i] * distance_array[i];
  }
  const numerical_types::real distance_sq_inv = 1.0 / distance_sq;
  // Compute the accelleration due to gravity.
  const numerical_types::real distance_inv = std::sqrt(distance_sq_inv);
  const numerical_types::real accelleration = mass * distance_sq_inv;

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] -= accelleration * distance_array[i] * distance_inv;
  }
}

// Tree

void Tree::rebuild(std::vector<Particle>& particles) {
  // Calculate the extent of the boundingbox of all particles.
  std::array<numerical_types::ndarray, 3> stats = Particle::positionExtent(particles);
  
  numerical_types::ndarray average = stats[0];
  numerical_types::ndarray min = stats[1];
  numerical_types::ndarray max = stats[2];

  // Calculate the largest width from average to either end of the bounding-box.
  numerical_types::real width = 0.0;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    if (std::abs(average[i] - min[i]) > width)
      width = std::abs(average[i] - min[i]);
    if (std::abs(max[i] - average[i]) > width)
      width = std::abs(max[i] - average[i]);
  }
  
  width *= 1.5;	// Have some margin in the size.
  // Center the root node on the average position of all the particles.
  // This should balance the tree somewhat when a few particles are far away. 
  this->root_node.reset(new Node(this, nullptr, 0, average, width));
  // Add all particles to the tree.
  for (Particle& particle: particles) {
    this->add(&particle, this->root_node.get());
  }
  this->root_node->aggregateQuantities();
}

bool Tree::update(std::vector<Particle>& particles) {
  assert(this->root_node.get() != nullptr);
  bool rebuild_required = false;

  // Attempt to update the tree, rebuilding if required.
  Timer::byName("Relocate")->set();
  for (Particle& particle: particles) {
    if (!this->relocate(&particle)) {
      rebuild_required = true;
      break;
    }
  }
  Timer::byName("Relocate")->reset();
  if (rebuild_required) {
    this->rebuild(particles);
  }
  else{
    Timer::byName("Aggregate")->set();
    this->root_node->aggregateQuantities();
    Timer::byName("Aggregate")->reset();
  }
  return !rebuild_required;
}

void Tree::add(Particle* particle, Node* node) {
  // Walks the tree until a leaf node is found and adds the particle to it.
  assert(node != nullptr);

  bool indices[numerical_types::num_dimensions];
  while (node->num_particles_contained > max_num_particles) {
    node->num_particles_contained++;

    node->indexOf(particle->position, indices);
    node = node->getSubnode(indices);
  }

  node->add(particle);
}

int Tree::computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const {
  int computations = this->root_node->computeAccelleration(particle, theta, epsilon, result);

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] *= numerical_types::G;
  }
  return computations;
}

bool Tree::relocate(Particle* particle) {
  Node* node = particle->containing_node;
  assert(node != nullptr);

  int i = 0;
  // Move the particle up in the tree if the parent node only
  // contains particles in this node.
  while (node->parent != nullptr &&
         node->parent->num_particles_contained <= max_num_particles) {
    i++;
    node->num_particles_contained--;
    node = node->parent;
  }
  // Mark the current node as the new leaf for all its contained
  // particles.
  Node* next_left_node = node;
  // Move the particle up in the tree if the current node does not
  // contain the particle position.
  while (node != nullptr && !node->contains(particle->position)) {
    i++;
    node->num_particles_contained--;
    node = node->parent;
  }
  // If the root node did not contain the particle, we cannot add it.
  if (node == nullptr){
    return false;
  }
  // Don't add the particle again if the node does not
  // change from any of the above.
  if (node == particle->containing_node) {
    return true;
  }
  // Subtract the particle from the node since adding it
  // will add one particle again.
  node->num_particles_contained--;

  // Remove the particle from its previous node.
  particle->containing_node->remove(particle);

  // If we can move particles up the tree to a new leaf node,
  // do it.
  if (next_left_node != particle->containing_node) {
    next_left_node->gatherChildParticles();
  }
  particle->containing_node = nullptr;

  // node->add(particle);
  this->add(particle, node);
  return true;
}

std::vector<Node*> Tree::getNodesAtDepth(int depth) {
  assert(depth >= 0);

  if (depth == 0) {
    return {this->root_node.get()};
  }
  std::vector<Node*> parents = this->getNodesAtDepth(depth - 1);
  std::vector<Node*> children;
  // Reserv capacity for the largest possible number of children.
  children.reserve(parents.size() * num_subnodes);

  for (Node* parent : parents) {
    if (parent->hasChildren()) {
      for (int i = 0; i < num_subnodes; i++) {
        children.push_back(parent->child(i));
      }
    }
  }
  // Shrink the vector to the actual number of children.
  children.shrink_to_fit();
  return children;
}

void Tree::zero() {
  this->zero(this->root_node.get());
}

void Tree::zero(Node* node) {
  assert(node != nullptr);
  std::stack<Node*> stack;
  stack.push(node);

  while(!stack.empty()) {
    Node* current = stack.top();
    stack.pop();
    current->total_mass = 0.0;
    current->dirty = true;

    if (!current->hasChildren()) {
      continue;
    }
    for (int i = 0; i < num_subnodes; i++) {
      if (!current->constChild(i)->dirty)
        stack.push(current->child(i));
    }
  }
}

void Tree::write(std::ofstream& file) const {
  std::stack<Node*> stack;
  stack.push(this->root_node.get());

  while(!stack.empty()) {
    Node* current = stack.top();
    stack.pop();

    if (current->hasChildren()) {
      for (int i = 0; i < num_subnodes; i++) {
        stack.push(current->child(i));
      }
    }

    file << current->width << ", ";
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      file << current->center[i] << ", ";
    }
  }
}

}  // namespace model
