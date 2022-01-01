#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/Auxiliary_Functions.h"
#include <cmath>
#include <complex>
#include <chrono>

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B){

	MatrixXcd C(A.rows() * B.rows(), A.cols() * B.rows());

	for (int i = 0; i < A.rows(); i++){
		for (int j = 0; j < A.cols(); j++){
			C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
		}
	}

	return C;
}

MatrixXcd product_vector_matrix(Vector3d vector, MatrixXcd paulimatrix_x, MatrixXcd paulimatrix_y, MatrixXcd paulimatrix_z){

	return vector[0]*paulimatrix_x + vector[1]*paulimatrix_y + vector[2]*paulimatrix_z;
}
