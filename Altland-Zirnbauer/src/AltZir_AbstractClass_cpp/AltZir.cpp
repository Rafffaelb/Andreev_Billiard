#include <iostream>
#include "../../include/AltZir_AbstractClass_h/AltZir.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

AltZir::AltZir() {};

AltZir::~AltZir() {};

void AltZir::Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) {};

void AltZir::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) {};

void AltZir::Create_H(MatrixXcd* H_pointer, int _ress, double _lambda) {};

void AltZir::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) {};

void AltZir::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1) {};

void AltZir::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps, int N1) {};
