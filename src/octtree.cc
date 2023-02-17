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
    numerical_types::ndarray center,
    numerical_types::real width,
    Node* parent) {
  memcpy(&this->center[0], &center[0], numerical_types::ndarray_bytes_size);
  this->width = width;
  this->parent = parent;
  
  memset(&this->center_of_mass[0], 0.0, numerical_types::ndarray_bytes_size);
  this->total_mass = 0.0;
}

Node* Node::getSubnode(bool indices[numerical_types::num_dimensions]) {
  int linear_index = Node::SubnodeIndex(indices);

  // If the node does not yet exist, create it.
  if (!this->hasChildren()) {
    this->allocateChildren();
  }

  return &this->children[linear_index];
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
  // If this is the first particle added, we are currently a leaf node.
  if (this->num_particles_contained == 1) {
    this->particle = particle;
    this->particle->containing_node = this;
    return;
  }

  bool indices[numerical_types::num_dimensions];
  // If this is the second particle, we must move the first into a subnode as well.
  if (this->particle != nullptr) {
    this->indexOf(this->particle->position, indices);
    this->getSubnode(indices)->add(this->particle);
    this->particle = nullptr;
  }

  // Add the current particle to the respective subnode.
  this->indexOf(particle->position, indices);
  this->getSubnode(indices)->add(particle);
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
    this->children.emplace_back(center, half_width, this);
  }
}

void Node::aggregateQuantities() {
  this->dirty = false;
  if (this->particle != nullptr) {
    this->total_mass = this->particle->mass;
    this->center_of_mass = this->particle->position;
    return;
  }

  this->total_mass = 0.0;
  if (!this->hasChildren()) {
    return;
  }
  for (int i = 0; i < num_subnodes; i++) {
    if (this->child(i)->num_particles_contained > 0){
      this->child(i)->aggregateQuantities();
      this->addMass(this->child(i)->total_mass,
                    this->child(i)->center_of_mass);
    }
  }
}

void Node::clear() {
  if (this->particle != nullptr)
    this->particle->containing_node = nullptr;

  this->particle = nullptr;
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
  if (this->particle == &particle) {
    return num_calculations;
  }
  else if (this->particle == nullptr) {
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
      if (this->hasChildren()) {
        for (int i = 0; i < num_subnodes; i++) {
          if (this->constChild(i)->num_particles_contained > 0)
            num_calculations += this->constChild(i)->computeAccelleration(
              particle, theta, epsilon, result);
        }
      }
      return num_calculations;
    }
  }

  // Recalculate the distance to center of mass, not center of the cell.
  numerical_types::ndarray distance_array;
  numerical_types::real distance_sq = epsilon;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    distance_array[i] = particle.position[i] - this->center_of_mass[i];
    distance_sq += distance_array[i] * distance_array[i];
  }
  const numerical_types::real distance_sq_inv = 1.0 / distance_sq;
  // Compute the accelleration due to gravity.
  const numerical_types::real distance_inv = std::sqrt(distance_sq_inv);
  const numerical_types::real accelleration = this->total_mass * distance_sq_inv;

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] -= accelleration * distance_array[i] * distance_inv;
  }

  return ++num_calculations;
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
  this->root_node.reset(new Node(average, width, nullptr));
  // Add all particles to the tree.
  for (Particle& particle: particles) {
    this->add(&particle, this->root_node.get());
  }
  // this->aggregateQuantities();
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
    // this->aggregateQuantities();
    Timer::byName("Aggregate")->reset();
  }
  return !rebuild_required;
}

void Tree::add(Particle* particle, Node* node) {
  // Walks the tree until a leaf node is found and adds the particle to it.
  assert(node != nullptr);

  bool indices[numerical_types::num_dimensions];
  while (node->num_particles_contained > 1) {
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
         node->parent->num_particles_contained == 1) {
    i++;
    node->num_particles_contained--;
    node = node->parent;
  }
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
  particle->containing_node->particle = nullptr;
  particle->containing_node = nullptr;

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
