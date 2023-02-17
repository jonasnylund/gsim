#pragma once

#include <array>
#include <fstream>
#include <memory>
#include <vector>

#include "numerical_types.h"
#include "particle.h"

namespace model {

constexpr int num_subnodes = (1 << numerical_types::num_dimensions);

class Node {
public:

  // Return the linear index corresponding to the multidimensional index
  // given.
  static inline int SubnodeIndex(const bool indices[numerical_types::num_dimensions]) {
    int index = 0;
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      index += indices[i] * (1 << i);
    }
    return index;
  }

  Node(
    numerical_types::ndarray center,
    numerical_types::real width,
    Node* parent);

  // Checks whether a particle is in the bounding box of this node.
  bool contains(const numerical_types::ndarray& point) const;

  // Recursively computes the total mass and center of mass
  // of this node and all subnodes.
  void aggregateQuantities();

  // Compute the accelleration due to gravity for the given position. Counts and returns
  // the number of interactions used to produce the result.
  int computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const;

  int getNumContainedParticles() const { return this->num_particles_contained; }

  Node* getSubnode(bool indices[numerical_types::num_dimensions]);

  // Removes nodes in the tree with no contiained particles.
  void prune();

  inline void indexOf(
      const Particle* particle,
      bool indices[numerical_types::num_dimensions]) const {
    for (int i = 0; i < numerical_types::num_dimensions; i++) {
      indices[i] = this->center[i] < particle->position[i];
    }
  }

  inline bool hasChildren() const { return !this->children.empty(); }
  inline Node* child(int index) { return &this->children[index]; };
  inline const Node* constChild(int index) const { return &this->children[index]; }

  inline numerical_types::real mass() const { return this->total_mass; }
  inline numerical_types::ndarray centerOfMass() const { return this->center_of_mass; }

 protected:

  // Add a particle to this node.
  void add(Particle* particle);

  void addMass(numerical_types::real mass, const numerical_types::ndarray& position);

  void allocateChildren();

  // Clear all particles from this node and all child nodes.
  void clear();

  numerical_types::ndarray center;
  numerical_types::real width;

  Node* parent = nullptr;
  std::vector<Node> children;
  Particle* particle = nullptr;
  int num_particles_contained = 0;
  bool dirty = true;
  
  numerical_types::ndarray center_of_mass = {0.0};
  numerical_types::real total_mass = 0.0;

  friend class Tree;
};


class Tree {
 public:
  // Builds and populates a new tree from a list of particles.
  void rebuild(std::vector<Particle>& particles);
  // Attempts to update the tree without reallocating. Returns
  // true on success, and false if the tree required rebuilding.
  bool update(std::vector<Particle>& particles);

  int computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const;

  Node* getRoot() const { return this->root_node.get(); }

  void write(std::ofstream& file) const;

 private:
  // Add a particle to the tree.
  void add(Particle* particle, Node* root);

  // Updates a particles position in the tree, moving it up until
  // either a cell containing it is found, or the root node
  // is reached.
  bool relocate(Particle* particle);

  // Returns a vector with all nodes at a certain tree depth.
  std::vector<Node*> getNodesAtDepth(int depth);

  // Clear the state of each node.
  void zero();
  void zero(Node* root);
  std::unique_ptr<Node> root_node;
};

}  // namespace model
