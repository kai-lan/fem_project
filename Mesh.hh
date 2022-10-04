/**
 * @file Mesh.hh
 * @author Kai Lan (kai.weixian.lan@gmail.com)
 * @brief 
 * 
 * @date 2022-05-10
 */
#ifndef MESH
#define MESH
#include <iostream>
#include <vector>
#include "Energy.hh"

constexpr double g = -9.80665; // gravity acceleration
constexpr int g_direction = 1; // Set y axis as the gravity direction
// TODO: add nonzero Neumann condition(external force)
class Mesh
{
    using MXd  = Eigen::MatrixXd;
    using MXdRM = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
    using MXi  = Eigen::MatrixXi;
    using VXi  = Eigen::VectorXi;
    using VXd  = Eigen::VectorXd;
    //using AXb  = Eigen::Array<bool, -1, 1>;
private:
    std::unique_ptr<Energy> energy_;
    std::vector<MXd> elementDms_; // Cache this
    VXd elementVolumes_; // Cache this
    VXd massLumpingMatrix_; // Cache this
    VXd weight_;
    void initializeElementDms();
    void initializeElementVolumes();
    void initializeMassMatrix();
    void initializeWeight();
    void enforceDirichletCondition(); // "Eliminate" corresponding rows and cols in mass matrix and weight vector
public:
    const int DIM;
    const double R0; // Assuming R0 constant in X for now
    const MXdRM V0; // Original positions
    const MXi E; // Indices of vertices of one element stored in one row. **Vertices should be stores in positive orientation.
    const VXi dirichletNodes; // Assuming those nodes have zero displacement
    const VXi neumannNodes; 
    const MXd externalForces; // Forces on Neumann BC

    int numNodes() const { return V0.rows(); }
    int numElements() const { return E.rows(); }
    const std::vector<MXd>& Dm() const { return elementDms_; }
    const VXd& volumes() const { return elementVolumes_; }
    const VXd& massMatrix() const { return massLumpingMatrix_; }
    const VXd& weight() const { return weight_; }

    Mesh(const MXd& V0, const MXi& E, double rho, const VXi& Dirichlet, const VXi& Neumann, const MXd& forces) 
    : DIM(V0.cols()), R0(rho), V0(V0), E(E), dirichletNodes(Dirichlet), neumannNodes(Neumann), externalForces(forces) {
        if (neumannNodes.size() != externalForces.rows()) throw std::runtime_error("Number of Neumann nodes and nodes under external forces unmatch!");
        energy_ = std::make_unique<Energy>(DIM);
        initializeElementDms();
        initializeElementVolumes();
        initializeMassMatrix();
        initializeWeight();
        enforceDirichletCondition();
    }
    Mesh(const MXd& V0, const MXi& E, double rho) : DIM(V0.cols()), R0(rho), V0(V0), E(E) {
        energy_ = std::make_unique<Energy>(DIM);
        initializeElementDms();
        initializeElementVolumes();
        initializeMassMatrix();
        initializeWeight();
    }

    //VXd flatten(const MXd& M) const { return VXd::Map(M.data(), M.size()); }
    VXd flatten(const MXdRM& M) const { return VXd::Map(M.data(), M.size()); }
    MXdRM unflatten(const VXd& V) const { return MXdRM::Map(V.data(), V.size()/DIM, DIM); }

    // phi_e(zeta) = X_0 + D_e * zeta
    // zeta in R^d, d = 2 for triangle mesh and d = 3 for tetrahedron mesh
    VXd X_from_zeta(int e, const VXd& zeta) const;

    // Inverse of phi_e = D_e^-1 * (X - X_0)
    VXd zeta_from_X(int e, const VXd& X) const;

    // N_e: X -> [Nhat_0(zeta), Nhat_1(zeta), Nhat_2(zeta)]^T
    VXd N_e(int e, const VXd& X) const;

    // Mass of an element
    double element_mass(int e) const;

    // Mass for each node of an element: element_mass / 3 for triangle and element_mass / 4 for tetrahedron
    double node_element_mass(int e) const;   
    
    // Deformation gradient F for element E, given phi_j(t) (same shape as V0) function
    MXd element_deformation(int e, const MXd &x) const;

    // Total potential energy
    double potential_energy(const MXd &x) const;

    // P used to assemble f
    MXd element_pressure(int e, const MXd &x) const;

    // f force
    VXd force(const MXd &x) const;

    // Eliminate corresponding rows and cols in force to enforce Dirichlet conditions
    void enforce_dirichlet_condition(VXd &forces) const;

    ~Mesh(){};
};
#endif /* MESH */
