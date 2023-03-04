#pragma once

#include <array>
#include <fstream>
#include <memory>
#include <vector>

#include "numerical_types.h"
#include "particle.h"

namespace model {

constexpr int num_subnodes = (1 << numerical_types::num_dimensions);
constexpr int max_num_particles = 8;


class Tree {
 public:
  // Builds and populates a new tree from a list of particles.
  void rebuild(std::vector<Particle>& particles);
  // Attempts to update the tree without reallocating. Returns
  // true on success, and false if the tree required rebuilding.
  bool update(std::vector<Particle>& particles);

  // Partially updates the nodes of the tree without relocating
  // the particles.
  void update(float fraction = 1.0f);

  // Compute the accelleration for a single particle.
  int computeAccelleration(
    const Particle& particle,
    numerical_types::real theta,
    numerical_types::real epsilon,
    numerical_types::ndarray& result) const;

  void clear();

  // Prune empty nodes from the tree.
  void prune();

  // Write the tree structure to a text file.
  void write(std::ofstream& file) const;

  // // Returns true if the tree is empty, false otherwise.
  inline bool empty() const {
    return this->rootNode() == nullptr || this->rootNode()->num_particles_contained == 0;
  }

  // Returns the maximum height of the tree.
  inline int height() const { return this->nodes.size(); }

  // Returns the total number of nodes allocated in the tree.
  int numNodes() const;

 protected:
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
      Tree* tree,
      Node* parent,
      numerical_types::NodeKey id,
      const numerical_types::ndarray& center,
      numerical_types::real width);

    Node(Node&&) = default;
    Node& operator=(Node&&) = default; 	

    // Checks whether a point is in the bounding box of this node.
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

    // Returns the total number of particles in this node and all its children.
    int getNumContainedParticles() const { return this->num_particles_contained; }

    // Computes the subnode containing the given point.
    inline void indexOf(
        const numerical_types::ndarray& position,
        bool indices[numerical_types::num_dimensions]) const {
      for (int i = 0; i < numerical_types::num_dimensions; i++) {
        indices[i] = this->center[i] < position[i];
      }
    }

    // Returns true if this node has child nodes, false otherwise.
    inline bool hasChildren() const { return this->children[0] != numerical_types::emptykey; }
    inline Node* child(int index) {
      return this->tree->getNode(this->children[index]);
    }
    inline const Node* child(int index) const {
      return this->tree->getNode(this->children[index]);
    }
    inline Node* parent() const { return this->tree->getNode(this->parent_id); }

    // Returns true if this node has particles, false otherwise.
    inline bool hasParticles() const { return this->num_particles_local > 0; }
    inline Particle *& particle(int index) { return this->particles[index].particle; }

  protected:
    struct PseudoParticle {
      Particle* particle = nullptr;
      numerical_types::real mass = 0.0;
      numerical_types::ndarray position = {};
    };

    // Add a particle to this node.
    void add(Particle* particle);

    void addMass(numerical_types::real mass, const numerical_types::ndarray& position);

    void allocateChildren();

    void computeAccelleration(
      numerical_types::real mass,
      const numerical_types::ndarray& position,
      const Particle& particle,
      numerical_types::real theta,
      numerical_types::real epsilon,
      numerical_types::ndarray& result) const;

    // Remove a particle from this node.
    void remove(Particle* particle);

    // Moves all contained particles into this node, if possible.
    void gatherChildParticles();

    void updateKey(numerical_types::NodeKey key);

    // Clear all particles from this node and all child nodes.
    void clear();

    Node* getSubnode(bool indices[numerical_types::num_dimensions]);

    numerical_types::ndarray center;
    numerical_types::real width;

    numerical_types::NodeKey id;
    numerical_types::NodeKey parent_id;

    Tree* tree;
    numerical_types::depth_key_t depth;
    std::array<numerical_types::NodeKey, num_subnodes> children;
    std::array<PseudoParticle, max_num_particles> particles;
    int num_particles_local = 0;
    int num_particles_contained = 0;
    bool dirty = true;
    
    numerical_types::ndarray center_of_mass;
    numerical_types::real total_mass = 0.0;

    friend class Tree;
  };


  inline Node* rootNode() {
    if (this->nodes.size() > 0 && this->nodes[0].size() > 0) {
      return &this->nodes[0][0];
    }
    return nullptr;
  }
  inline const Node* rootNode() const {
    if (this->nodes.size() > 0 && this->nodes[0].size() > 0) {
      return &this->nodes[0][0];
    }
    return nullptr;
  }

  inline Node* getNode(numerical_types::NodeKey key) {
    if (key == numerical_types::emptykey)
      return nullptr;
    return &this->nodes[key.depth][key.index];
  }
  inline const Node* getNode(numerical_types::NodeKey key) const {
    if (key == numerical_types::emptykey)
      return nullptr;
    return &this->nodes[key.depth][key.index];
  }

  numerical_types::NodeKey allocateNode(
    Node* parent, const numerical_types::ndarray& center, numerical_types::real width);

 private:
  // Add a particle to the tree.
  void add(Particle* particle, Node* node);

  // Updates a particles position in the tree, moving it up until
  // either a cell containing it is found, or the root node
  // is reached.
  bool relocate(Particle* particle);

  std::vector<std::vector<Node>> nodes;
  int rebuilds = 0;

  friend class Node;
};

}  // namespace model
