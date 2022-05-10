/**
 * @file Mesh.hh
 * @author Kai Lan (kai.weixian.lan@gmail.com)
 * @brief 
 * 
 * @date 2022-05-10
 */
#ifndef MESH
#define MESH

#include <Eigen/Dense>

class Mesh
{
    using Vertices = Eigen::Matrix<double, Eigen::Dynamic, 3>; // Eigen::ColMajor or Eigen::RowMajor
private:
    int DIM;
    Vertices m_V;
    /* data */
public:
    Mesh(/* args */);
    ~Mesh();
};

Mesh::Mesh(/* args */)
{
}

Mesh::~Mesh()
{
}


#endif /* MESH */
