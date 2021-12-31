#ifndef ALTZIR_DIII_H
#define ALTZIR_DIII_H

#include "AltZir_AbstractClass_h/AltZir.h"

class AltZir_DIII: public AltZir{

	public:
		AltZir_DIII(double lambda, int num_steps, int spin_deg, int electron_hole_deg);

		~AltZir_DIII();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda);

		void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps);
		void Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1);
		void Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps, int N1);
		void Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps);
};
#endif
