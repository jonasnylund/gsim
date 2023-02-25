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
  this->center_of_mass.fill(0.0);
  this->total_mass = 0.0;
  this->particles.fill({});
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
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particle(i) == nullptr) {
        this->particle(i) = particle;
        this->num_particles_local++;
        particle->containing_node = this;
        return;
      }
    }
  }

  // If we have local particles, but no more of them will fit, we must
  // first move them to the correct node.
  bool indices[numerical_types::num_dimensions];
  if (this->num_particles_local > 0) {
    assert(this->num_particles_local == max_num_particles);
    for (int i = 0; i < max_num_particles; i++) {
      this->indexOf(this->particle(i)->position, indices);
      this->getSubnode(indices)->add(this->particle(i));
      this->particle(i) = nullptr;
    }
    this->num_particles_local = 0;
  }

  // Add the current particle to the respective subnode.
  this->indexOf(particle->position, indices);
  this->getSubnode(indices)->add(particle);
}

void Node::remove(Particle* particle) {
  for (int i = 0; i < max_num_particles; i++) {
    if (this->particle(i) == particle) {
      this->particle(i) = {};
      this->num_particles_local--;
      break;
    }
  }
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
        if (this->child(i)->particle(j) == nullptr){
          continue;
        }
        this->add(this->child(i)->particle(j));
      }
      this->child(i)->particles.fill({});
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
}

void Node::aggregateQuantities() {
  this->dirty = false;
  this->total_mass = 0.0;

  if (this->hasParticles()) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particle(i) != nullptr) {
        this->particles[i].mass = this->particle(i)->mass;
        this->particles[i].position = this->particle(i)->position;
        this->addMass(this->particles[i].mass,
                      this->particles[i].position);
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
      if (this->particle(i) != nullptr) {
        this->particle(i)->containing_node = nullptr;
      }
      this->particle(i) = nullptr;
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

  if (theta * distance_sq > width_sq) {
    // If the particle is sufficiently far away, approximate the accelleration
    // from the cells mass and center of mass.
    this->computeAccelleration(this->total_mass, this->center_of_mass, particle, theta, epsilon, result);
    num_calculations++;
  }
  else if (this->hasParticles()) {
    // If this node is a leaf node, compute the accelleration against each
    // particle individually.
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particles[i].particle == nullptr ||
          this->particles[i].particle == &particle) {
        continue;
      }
      // Compute the accelleration against this particle.
      this->computeAccelleration(
        this->particles[i].mass,
        this->particles[i].position,
        particle,
        theta,
        epsilon,
        result);
      num_calculations++;
    }
  }
  else if (this->hasChildren()) {
    // If the node has children, iterate over each and compute the accelleration.
    for (int i = 0; i < num_subnodes; i++) {
      if (this->constChild(i)->num_particles_contained > 0)
        num_calculations += this->constChild(i)->computeAccelleration(
          particle, theta, epsilon, result);
    }
  }
  return num_calculations;
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
  Timer::byName("Tree: rebuild")->set();
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
  Timer::byName("Tree: rebuild")->reset();
}

bool Tree::update(std::vector<Particle>& particles) {
  assert(this->root_node.get() != nullptr);
  bool rebuild_required = false;

  // Attempt to update the tree, rebuilding if required.
  Timer::byName("Tree: relocate")->set();
  for (Particle& particle: particles) {
    if (!this->relocate(&particle)) {
      rebuild_required = true;
      break;
    }
  }
  Timer::byName("Tree: relocate")->reset();
  if (rebuild_required) {
    this->rebuild(particles);
  }
  else{
    Timer::byName("Tree: aggregate")->set();
    this->root_node->aggregateQuantities();
    Timer::byName("Tree: aggregate")->reset();
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
  Node* next_leaf_node = node;
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
  if (next_leaf_node != particle->containing_node) {
    next_leaf_node->gatherChildParticles();
  }
  particle->containing_node = nullptr;

  // node->add(particle);
  this->add(particle, node);
  return true;
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
