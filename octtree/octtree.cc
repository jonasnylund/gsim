#include "gsim/octtree/octtree.h"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <fstream>
#include <memory>
#include <stack>
#include <vector>

#include <omp.h>

#include "gsim/common/numerical_types.h"
#include "gsim/common/timers.h"
#include "gsim/octtree/particle.h"


namespace gsim {

Tree::Node::Node(
    Tree* tree,
    Tree::Node* parent,
    numerical_types::NodeKey id,
    const numerical_types::ndarray& center,
    numerical_types::real width)
    : id(id),
      parent_id(parent != nullptr ? parent->id : numerical_types::emptykey),
      tree(tree),
      center(center),
      width(width),
      depth(id.depth) {
  this->center_of_mass.fill(0.0);
  this->total_mass = 0.0;
  this->particles.fill({});
  this->children.fill(numerical_types::emptykey);
}

void Tree::Node::addMass(numerical_types::real mass, const numerical_types::ndarray& position) {
  this->total_mass += mass;
  const numerical_types::real com_scale = mass / this->total_mass;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    this->center_of_mass[i] *= (1.0 - com_scale);
    this->center_of_mass[i] += position[i] * com_scale;
  }
}

void Tree::Node::add(Particle* particle) {
  this->num_particles_contained++;

  // If this one of the first particle added, we are currently a leaf node.
  if (this->num_particles_contained <= this->particles.size()) {
    this->particle(this->num_particles_local++) = particle;
    particle->containing_node = this->id;
    return;
  }

  // If we don't have child nodes yet, create them.
  if (!this->hasChildren()) {
    this->allocateChildren();
  }

  // If we have local particles, but no more of them will fit, we must
  // first move them to the correct node.
  int indices[numerical_types::num_dimensions];
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

void Tree::add(Particle* particle, Tree::Node* node) {
  // Walks the tree until a leaf node is found and adds the particle to it.
  // assert(node != nullptr);

  int indices[numerical_types::num_dimensions];
  // while (node->num_particles_contained > max_num_particles) {
  while (node->hasChildren()) {
    node->num_particles_contained++;

    node->indexOf(particle->position, indices);
    node = node->getSubnode(indices);
  }

  node->add(particle);
}

void Tree::clear() {
  for (auto& nodes : this->nodes) {
    nodes.clear();
  }
  assert(this->empty());
}

void Tree::Node::remove(Particle* particle) {
  for (int i = 0; i < max_num_particles; i++) {
    if (this->particle(i) == particle) {
      this->particle(i) = {};
      this->num_particles_local--;
      this->num_particles_contained--;
      break;
    }
  }
  particle->containing_node = numerical_types::emptykey;
}

void Tree::Node::updateKey(numerical_types::NodeKey key) {
  Node* parent = this->parent();
  if (parent != nullptr) {
    for (int i = 0; i < num_subnodes; i++) {
      if (parent->children[i] == this->id) {
        parent->children[i] = key;
        break;
      }
    }
  }
  if (this->hasParticles()) {
    for (int i = 0; i < max_num_particles; i++) {
      if (this->particle(i) != nullptr) {
        this->particle(i)->containing_node = key;
      }
    }
  }
  if (this->hasChildren()) {
    for (int i = 0; i < num_subnodes; i++) {
      this->child(i)->parent_id = key;
    }
  }

  this->id = key;
}

void Tree::Node::gatherChildParticles() {
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
      for (int j = 0; j < this->child(i)->num_particles_local; j++) {
        this->add(this->child(i)->particle(j));
      }
      this->child(i)->particles.fill({});
      this->child(i)->num_particles_local = 0;
      this->child(i)->num_particles_contained = 0;
    }
  }
  assert(this->num_particles_contained <= max_num_particles);
}

void Tree::Node::allocateChildren() {
  assert(!this->hasChildren());
  const numerical_types::real half_width = this->width / 2.0;
  int indices[numerical_types::num_dimensions];

  for (int i = 0; i < num_subnodes; i++) {
    numerical_types::ndarray center;
    for (int j = 0; j < numerical_types::num_dimensions; j++) {
      int index = (i >> j) & 0x1;
      center[j] = this->center[j] - half_width + this->width * index;
    }
    this->indexOf(center, indices);
    assert(this->SubnodeIndex(indices) == i);

    const numerical_types::NodeKey key = this->tree->allocateNode(this, center, half_width);
    this->children[i] = key;
  }
}

void Tree::Node::clear() {
  if (this->hasChildren()) {
    for (int i = 0; i < num_subnodes; i++) {
      if (this->num_particles_contained > 0)
        this->child(i)->clear();
    }
  }
  this->particles.fill({});
  this->num_particles_local = 0;
  this->num_particles_contained = 0;
}

bool Tree::Node::contains(const numerical_types::ndarray& point) const {
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    if (point[i] < this->center[i] - width)
      return false;
    if (point[i] > this->center[i] + width)
      return false;
  }
  return true;
}

int Tree::Node::computeAccelleration(
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
  constexpr numerical_types::real root = std::sqrt(numerical_types::num_dimensions);
  const numerical_types::real distance_to_edge = std::sqrt(distance_sq) - root * this->width;

  if (theta * distance_to_edge > this->width) {
    // If the particle is sufficiently far away, approximate the accelleration
    // from the cells mass and center of mass.
    this->computeAccelleration(
      this->total_mass, this->center_of_mass, particle, theta, epsilon, result
    );
    num_calculations++;
  }
  else if (this->hasParticles()) {
    // If this node is a leaf node, compute the accelleration against each
    // particle individually.
    for (int i = 0; i < this->num_particles_local; i++) {
      if (this->particles[i].particle == &particle) {
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
      if (this->child(i)->num_particles_contained > 0) {
        num_calculations += this->child(i)->computeAccelleration(
          particle, theta, epsilon, result
        );
      }
    }
  }
  return num_calculations;
}

void Tree::Node::computeAccelleration(
    numerical_types::real mass,
    const numerical_types::ndarray& position,
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const {
  numerical_types::ndarray distance_array;
  numerical_types::real distance_sq = epsilon;
  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    distance_array[i] = position[i] - particle.position[i];
    distance_sq += distance_array[i] * distance_array[i];
  }
  const numerical_types::real distance_sq_inv = 1.0 / distance_sq;
  // Compute the accelleration due to gravity.
  const numerical_types::real distance_inv = std::sqrt(distance_sq_inv);
  const numerical_types::real accelleration = mass * distance_sq_inv;

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] += accelleration * distance_array[i] * distance_inv;
  }
}

int Tree::computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const {
  int computations = this->rootNode()->computeAccelleration(particle, theta, epsilon, result);

  for (int i = 0; i < numerical_types::num_dimensions; i++) {
    result[i] *= numerical_types::G;
  }
  return computations;
}

// Tree
void Tree::rebuild(std::vector<Particle>& particles) {
  // Clear the tree of nodes.
  this->clear();

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
  
  Timer::byName("Tree: allocate")->set();
  width *= 1.5;	// Have some margin in the size.
  // Center the root node on the average position of all the particles.
  // This should balance the tree somewhat when a few particles are far away.
  this->allocateNode(nullptr, average, width);
  // Add all particles to the tree.
  for (Particle& particle: particles) {
    particle.containing_node = numerical_types::emptykey;
    this->add(&particle, this->rootNode());
  }
  Timer::byName("Tree: allocate")->reset();
  this->update();
  this->rebuilds++;
}

bool Tree::relocate(std::vector<Particle>& particles) {
  assert(this->rootNode() != nullptr);

  Timer::byName("Tree: relocate")->set();
  #pragma omp parallel for schedule(guided)
  for (auto& nodes : this->nodes) {
    for (auto& node : nodes) {
      node.num_particles_contained = 0;
      node.num_particles_local = 0;
      node.particles.fill({});
      node.total_mass = 0.0;
    }
  }

  // Attempt to update the tree, rebuilding if required.
  bool rebuild_required = false;
  #pragma omp parallel for reduction(|| : rebuild_required)
  for (Particle& particle: particles) {
    // Find the position of each particle in the tree in parallel,
    // record the position with the particle.
    if (this->rootNode()->contains(particle.position)) {
      const numerical_types::NodeKey key = this->getContainingNode(particle.position);
      particle.containing_node = key;
    } else {
      rebuild_required = true;
    }
  }

  if (!rebuild_required) {
    // Add all particles to the respective nodes serially, avoiding
    // races if a node needs to be split into subnodes.
    for (Particle& particle: particles) {
      this->getNode(particle.containing_node)->add(&particle);
    }
  } else {
    this->rebuild(particles);
  }
  Timer::byName("Tree: relocate")->reset();
  return !rebuild_required;
}

void Tree::update() {
  Timer::byName("Tree: aggregate")->set();
  const int deepest_index = this->height() - 1;

  // Compute node quantities in parallel at each depth. Begin with
  // the leaf nodes and iterate over each depth in the tree in reverse.
  for (int depth = deepest_index; depth >= 0; depth--) {
    const int num_iterations = this->nodes[depth].size();
    #pragma omp parallel for \
            if (num_iterations > 32)
    for (int i = 0; i < num_iterations; i++) {
      Node& node = this->nodes[depth][i];
      node.total_mass = 0.0;
      node.num_particles_contained = 0;

      if (node.hasChildren()) {
        for (int j = 0; j < num_subnodes; j++) {
          if (node.child(j)->num_particles_contained > 0) {
            node.num_particles_contained += node.child(j)->num_particles_contained;
            node.addMass(node.child(j)->total_mass,
                         node.child(j)->center_of_mass);
          }
        }
        // Move particles up the tree if they can be.
        if (node.num_particles_contained > 0 &&
            node.num_particles_contained <= max_num_particles) {
          node.gatherChildParticles();
        }
      }
      // (Re)Compute center of mass.
      if (node.hasParticles()) {
        node.total_mass = 0.0;
        for (int j = 0; j < node.num_particles_local; j++) {
          // if (node.particle(j) != nullptr) {
            node.particles[j].mass = node.particle(j)->mass;
            node.particles[j].position = node.particle(j)->position;
            node.addMass(node.particles[j].mass,
                         node.particles[j].position);
          // }
        }
        node.num_particles_contained = node.num_particles_local;
      }
    }
  }

  Timer::byName("Tree: aggregate")->reset();
}

numerical_types::NodeKey Tree::allocateNode(
    Tree::Node* parent,
    const numerical_types::ndarray& center,
    numerical_types::real width) {
  numerical_types::depth_key_t depth;
  if (parent != nullptr) {
    depth = parent->depth + 1;
  }
  else {
    depth = 0;
  }

  assert(this->height() + 1 >= depth);
  if (this->height() <= depth) {
    this->nodes.resize(depth + 1);
  }

  std::vector<Node>* node_vector = &this->nodes[depth];
  const numerical_types::index_key_t index = node_vector->size();

  numerical_types::NodeKey key = {depth, index};
  node_vector->emplace_back(this, parent, key, center, width);

  return key;
}

void Tree::prune() {
  Timer::byName("Tree: pruning")->set();

  // Start pruning at the outermost leaf nodes.
  // Don't prune all the way to the root node.
  for (int depth = this->height() - 1; depth > 0; depth--) {
    std::vector<Node>& node_vector = this->nodes[depth];

    int index = 0;
    while (index < node_vector.size()) {
      Node& node = node_vector[index];
      // Parent is never null since depth > 0.
      Node* parent = node.parent();

      // Don't prune nodes whos' parent node has particles.
      if (parent->num_particles_contained > 0) {
        // Skip all nodes that share this parent, as they should
        // all be pruned together.
        index += num_subnodes;
        continue;
      }
      // Get iterators to the range to prune.
      std::vector<Node>::iterator erase_begin = node_vector.begin() + index;
      std::vector<Node>::iterator erase_end = node_vector.begin() + index + num_subnodes;

      // Erase the nodes, and update the parent to be child-less.
      node_vector.erase(erase_begin, erase_end);
      parent->children.fill(numerical_types::emptykey);
    }

    // Recalculate the ids for the nodes that remain.
    for (index = 0; index < node_vector.size(); index++) {
      numerical_types::NodeKey key;
      key.depth = depth;
      key.index = index;
      Node& node = node_vector[index];

      if (node.id == key) {
        // Nothing to update.
        continue;
      }
      node.updateKey(key);
    } 
  }

  Timer::byName("Tree: pruning")->reset();
}

numerical_types::NodeKey Tree::getContainingNode(
  numerical_types::ndarray position) const {
  const Tree::Node* node = this->rootNode();

  int indices[numerical_types::num_dimensions];
  while (node->hasChildren()) {
    node->indexOf(position, indices);
    node = node->getSubnode(indices);
  }
  return node->id;
}

int Tree::numNodes() const {
  int num_nodes = 0;
  for (int i = 0; i < this->height(); i++) {
    num_nodes += this->nodes[i].size();
  }
  return num_nodes;
}

std::pair<int, int> Tree::countLeafNodes() const {
  int num_nodes = 0;
  int num_particles = 0;
  for (const auto& nodes : this->nodes) {
    for (const auto& node : nodes) {
      int local_particles = 0;
      for (int i = 0; i < max_num_particles; i++) {
        if (node.particles[i].particle != nullptr) {
          local_particles++;
        }
      }
      assert(local_particles == node.num_particles_local);
      if (local_particles > 0) {
        assert(local_particles == node.num_particles_contained);
        num_nodes++;
        num_particles += local_particles;
      }
    }
  }
  return {num_nodes, num_particles};
}

void Tree::write(std::ofstream& file) const {
  const size_t real_bytes = sizeof(numerical_types::real);
  std::stack<const Tree::Node*> stack;
  stack.push(this->rootNode());
  while(!stack.empty()) {
    const Tree::Node* current = stack.top();
    stack.pop();

    if (current->hasChildren()) {
      for (int i = 0; i < num_subnodes; i++) {
        stack.push(current->child(i));
      }
    }

    file.write(reinterpret_cast<const char*>(&current->width), real_bytes);
    file.write(reinterpret_cast<const char*>(current->center.data()),
               numerical_types::ndarray_bytes_size);
  }
}

}  // namespace gsim
