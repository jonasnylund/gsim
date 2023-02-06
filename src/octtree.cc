#include "octtree.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <memory>
#include <vector>

#include "particle.h"

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
  if (this->children[linear_index] == nullptr) {
    numerical_types::ndarray center;
    numerical_types::real half_width = this->width / 2;
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      center[i] = this->center[i] - half_width + this->width * indices[i];
    }
    this->children[linear_index] = std::make_unique<Node>(center, half_width, this);
  }

  return this->children[linear_index].get();
}

void Node::addMass(numerical_types::real mass, const numerical_types::ndarray& position) {
  this->total_mass += mass;
  numerical_types::real com_scale = mass / this->total_mass;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->center_of_mass[i] += position[i] * com_scale;
  }
  this->num_particles_contained++;
}

void Node::add(const Particle* particle) {
  // Update internal state
  this->addMass(particle->mass, particle->position);
  
  // If this is the first particle added, we are currently a leaf node.
  if (this->num_particles_contained == 1) {
    this->particle = particle;
    return;
  }

  bool indices[numerical_types::num_dimensions];
  // If this is the second particle, we must move the first into a subnode as well.
  if (this->particle != nullptr) {
    this->indexOf(this->particle, indices);
    this->getSubnode(indices)->add(this->particle);
    this->particle = nullptr;
  }

  // Add the current particle to the respective subnode.
  this->indexOf(particle, indices);
  this->getSubnode(indices)->add(particle);
}

void Node::clear() {
  this->particle = nullptr;
  this->num_particles_contained = 0;
  this->total_mass = 0;
  memset(&this->center_of_mass[0], 0.0, numerical_types::ndarray_bytes_size);

  for (int i = 0; i < num_subnodes; i++) {
    if (this->children[i] != nullptr && this->children[i]->num_particles_contained > 0) {
      this->children[i]->clear();
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

void Node::prune(int depth) {
  for (int i = 0; i < num_subnodes; i++) {
    if (depth <= 0 && this->children[i] != nullptr) {
      if (this->children[i]->num_particles_contained == 0) {
        this->children[i].release();
      }
    }
    if (this->children[i] != nullptr) {
      this->children[i]->prune(depth - 1);
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
  // Calculate the distance to the object.
  numerical_types::real distance_sq = 0;
  numerical_types::ndarray distance_array;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    distance_array[i] = particle.position[i] - this->center_of_mass[i];
    distance_sq += distance_array[i] * distance_array[i];
  }

  // If we have subnodes and the particle is close.
  if (this->particle == nullptr && theta * distance_sq < this->width * this->width) {
    for (int i = 0; i < num_subnodes; i++) {
      if (this->children[i] != nullptr &&
      this->children[i]->num_particles_contained > 0)
        num_calculations += this->children[i]->computeAccelleration(
          particle, theta, epsilon, result);
    }
    return num_calculations;
  }
  numerical_types::real distance_sq_inv = 1 / (distance_sq + epsilon);
  // Compute the accelleration due to gravity.
  numerical_types::real distance_inv = std::sqrt(distance_sq_inv);
  numerical_types::real force = numerical_types::G * this->total_mass * distance_sq_inv;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] -= force * distance_array[i] * distance_inv;
  }
  return ++num_calculations;
}

void Tree::rebuild(const std::vector<Particle>& particles) {
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
  // Add all particles to the tree recursively.
  for (const Particle& particle : particles) {
    this->root_node->add(&particle);
  }
}

void Tree::update(const std::vector<Particle>& particles) {
  this->root_node->clear();
  for (const Particle& particle: particles) {
    this->add(&particle);
    // this->root_node->add(&particle);
  }
}

Node* Tree::getRoot() const {
  return this->root_node.get();
}

void Tree::add(const Particle* particle) {
  // Walks the tree until a leaf node is found and adds the particle to it.

  Node* current_node = this->root_node.get();

  bool indices[numerical_types::num_dimensions];
  while (current_node->getNumContainedParticles() > 1) {
    current_node->addMass(particle->mass, particle->position);

    current_node->indexOf(particle, indices);
    current_node = current_node->getSubnode(indices);
  }

  current_node->add(particle);
}

void Tree::clear() {
  this->root_node->clear();
}

}  // namespace model
