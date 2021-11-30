#ifndef QUANTUM_CHAOTIC_BILLIARD_H
#define QUANTUM_CHAOTIC_BILLIARD_H

using namespace std;
using namespace Eigen;

class Quantum_chaotic_billiard{

	private:
		int _N1;
		int _N2;
		int _electron_hole_deg;

		complex<double> _G;
		complex<double> _P;

		MatrixXcd _H;
		MatrixXcd _W;
		MatrixXcd _C1;
		MatrixXcd _C2;
		MatrixXcd _S;

	public:
		Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2, int electron_hole_deg);

		void Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2, int electron_hole_deg);

		void Calculate_Smatrix();

		void Calculate_G_and_P();

		complex<double> getG();

		complex<double> getP();
};

#endif
