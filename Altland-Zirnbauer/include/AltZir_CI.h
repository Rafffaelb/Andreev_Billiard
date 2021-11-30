#ifndef ALTZIR_CI_H
#define ALTZIR_CI_H

#include "AltZir_AbstractClass_h/AltZir.h"

class AltZir_CI: public AltZir{

	public:
		AltZir_CI(double lambda, int num_steps, int spin_deg, int electron_hole_deg);

		~AltZir_CI();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda);

		void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps);
};
#endif
