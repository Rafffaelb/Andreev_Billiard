#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

using namespace std;
using namespace Eigen;

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B);

MatrixXcd product_vector_matrix(Vector3d vector, MatrixXcd paulimatrix_x, MatrixXcd paulimatrix_y, MatrixXcd paulimatrix_z);

#endif
