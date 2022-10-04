#include "Mesh.hh"

void Mesh::initializeElementDms() {
    for (int e = 0; e < numElements(); ++e) {
        MXd D = MXd::Zero(DIM, DIM);
        for (int i = 0; i < DIM; ++i) 
            D.col(i) = V0.row(E(e, i+1)) - V0.row(E(e, 0)); // [X_1 - X_0, X_2 - X0] for 2D, [X_1 - X_0, X_2 - X_0, X_3 - X_0] for 3D
        elementDms_.push_back(D);
    }
}
void Mesh::initializeElementVolumes() {
    elementVolumes_.setZero(numElements());
    for (int e = 0; e < numElements(); ++e) 
        elementVolumes_(e) = elementDms_[e].determinant() / (DIM * (DIM - 1));
}
// Mass lumping splits mass of an element evenly to its vertices.
// We will loop over each element to distribute its weight to vertices at appropriate locations
// It is a diagonal matrix, so we store it as a long vector.
void Mesh::initializeMassMatrix() {
    massLumpingMatrix_.setZero(DIM * numNodes());
    const int numVerticesPerEle = DIM + 1;
    for (int e = 0; e < numElements(); ++e) {
        double m = node_element_mass(e);
        for (int ie = 0; ie < numVerticesPerEle; ++ie) {
            int mesh_e = E(e, ie); // global index of ie th node of element e
            massLumpingMatrix_.segment(DIM * mesh_e, DIM) += m * VXd::Ones(DIM);
        }
    }
}
void Mesh::enforceDirichletCondition() {
    for (int nd = 0; nd < dirichletNodes.size(); ++nd) {
        massLumpingMatrix_.segment(DIM * dirichletNodes(nd), DIM) = VXd::Ones(DIM);
        weight_.segment(DIM * dirichletNodes(nd), DIM) = VXd::Zero(DIM);
    }
}
// Build the weight vector
void Mesh::initializeWeight() {
    weight_.setZero(DIM * numNodes());
    const int numVerticesPerEle = DIM + 1;
    VXd gravity = VXd::Zero(DIM);
    gravity(g_direction) = g;
    for (int e = 0; e < numElements(); ++e) {
        double m = node_element_mass(e);
        for (int ie = 0; ie < numVerticesPerEle; ++ie) {
            int mesh_e = E(e, ie);
            weight_.segment(DIM * mesh_e, DIM) += m * gravity;
        }     
    }
}

// phi_e(zeta) = X0 + D_e * zeta
// zeta in R^d, d = 2 for triangle mesh and d = 3 for tetrahedron mesh
Mesh::VXd Mesh::X_from_zeta(int e, const VXd& zeta) const {
    return V0.row(E(e, 0)).transpose() + elementDms_[e] * zeta;
}
// Inverse of phi_e = D_e^-1 * (X - X0)
Mesh::VXd Mesh::zeta_from_X(int e, const VXd& X) const {
    return elementDms_[e].inverse() * (X - V0.row(E(e, 0)).transpose());
}
// N_e: X -> [Nhat_0(zeta), Nhat_1(zeta), Nhat_2(zeta)]^T
Mesh::VXd Mesh::N_e(int e, const VXd& X) const {
    VXd zeta = zeta_from_X(e, X);
    VXd N = VXd::Zero(DIM + 1);
    N(0) = 1;
    for (int i = 1; i < DIM + 1; ++i) {
        N(0) -= zeta(i - 1);
        N(i) = zeta(i - 1); 
    }
    return N;
}

// Mass of an element
double Mesh::element_mass(int e) const {
    return R0 * elementVolumes_[e];
}

// Mass for each node of an element: element_mass / 3 for triangle and element_mass / 4 for tetrahedron
double Mesh::node_element_mass(int e) const {
    return element_mass(e) / (DIM + 1);
}

// Deformation gradient F for element E, given phi_j(t) (same shape as V0) function
Mesh::MXd Mesh::element_deformation(int e, const MXd &x) const {
    MXd edges = MXd::Zero(DIM, DIM);
    for (int i = 0; i < DIM; ++i)
        edges.col(i) = x.row(E(e, i+1)) - x.row(E(e, 0));
    return edges * elementDms_[e].inverse();
}

// Total potential energy
double Mesh::potential_energy(const MXd &x) const {
    double pe = 0;
    for (int e = 0; e < numElements(); ++e)
        pe += energy_->energy(element_deformation(e, x)) * elementVolumes_(e);
    return pe;
}
// Pressure for element e
Mesh::MXd Mesh::element_pressure(int e, const MXd &x) const {
    return energy_->pressure(element_deformation(e, x));
}
// f force
Mesh::VXd Mesh::force(const MXd &x) const {
    const int numVerticesPerEle = DIM + 1;
    VXd f = VXd::Zero(DIM * numNodes());
    MXd f_e = MXd::Zero(DIM, numVerticesPerEle);
    for (int e = 0; e < numElements(); ++e) {
        f_e.rightCols(DIM) = - elementVolumes_(e) * element_pressure(e, x) * elementDms_[e].inverse().transpose();
        f_e.col(0) = - f_e.rightCols(DIM).rowwise().sum();
        for (int ie = 0; ie < numVerticesPerEle; ++ie) {
            int mesh_i = E(e, ie);
            f.segment(DIM * mesh_i, DIM) += f_e.col(ie);
        }
    }
    enforce_dirichlet_condition(f);
    return f;
}

void Mesh::enforce_dirichlet_condition(VXd &force) const {
    for (int nd = 0; nd < dirichletNodes.size(); ++nd)
        force.segment(DIM * dirichletNodes(nd), DIM) = VXd::Zero(DIM);
}