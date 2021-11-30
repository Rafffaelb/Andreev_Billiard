#include <iostream>
#include <cstring>
#include "include/AltZir_AbstractClass_h/AltZir.h"
#include <cmath>
#include "include/AltZir_D.h"
#include "include/AltZir_DIII.h"
#include "include/AltZir_C.h"
#include "include/AltZir_CI.h"

using namespace std;

int main(int argc, char **argv){

	double lambda;
	int num_steps, spin_deg, electron_hole_deg;
	
	lambda = 0.5;
	num_steps = 100000;

	for (int i = 1; i < argc; i++){

		if (strcmp(argv[i],"AltZir_D") == 0){
			
			electron_hole_deg = 1;		
			spin_deg = 1;

			AltZir_D altzir_d(lambda, num_steps, spin_deg, electron_hole_deg);
		
			for (int j = 1; j < argc; j++){
				
				if (strcmp(argv[j],"Channel") == 0){
				
					cout << "\n ###### Running Class D Altland-Zirnbauer Ensemble (variable: Channel) ##### \n" << endl;

					altzir_d.Run_Simulation_Conductance_Channels();
				}
			}
			altzir_d.~AltZir_D();
		}
		else{
			if (strcmp(argv[i],"AltZir_DIII") == 0){
				
				electron_hole_deg = 1;
				spin_deg = 2;

				AltZir_DIII altzir_diii(lambda, num_steps, spin_deg, electron_hole_deg);

				for (int j = 1; j < argc; j++){
				
					if (strcmp(argv[j],"Channel") == 0){
				
						cout << "\n ###### Running Class DIII Altland-Zirnbauer Ensemble (variable: Channel) ##### \n" << endl;
						
						altzir_diii.Run_Simulation_Conductance_Channels();
					}
				}
				altzir_diii.~AltZir_DIII();	
			}
			else{

				if (strcmp(argv[i],"AltZir_C") == 0){

					electron_hole_deg = 2;
					spin_deg = 1;

					AltZir_C altzir_c(lambda, num_steps, spin_deg, electron_hole_deg);

					for (int j = 1; j < argc; j++){
				
						if (strcmp(argv[j],"Channel") == 0){
				
							cout << "\n ###### Running Class C Altland-Zirnbauer Ensemble (variable: Channel) ##### \n" << endl;
						
							altzir_c.Run_Simulation_Conductance_Channels();
						}
					}
					altzir_c.~AltZir_C();
				}
				else{
	
					if (strcmp(argv[i],"AltZir_CI") == 0){

						electron_hole_deg = 2;
						spin_deg = 1;

						AltZir_CI altzir_ci(lambda, num_steps, spin_deg, electron_hole_deg);

						for (int j = 1; j < argc; j++){
				
							if (strcmp(argv[j],"Channel") == 0){
				
								cout << "\n ###### Running Class CI Altland-Zirnbauer Ensemble (variable: Channel) ##### \n" << endl;
						
								altzir_ci.Run_Simulation_Conductance_Channels();
							}
						}
						altzir_ci.~AltZir_CI();	
					}
				}
			}
		}
	}
	return 0;
}
