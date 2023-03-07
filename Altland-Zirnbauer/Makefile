CC = g++

all:
	$(CC) -mfma -fopenmp -O3 main.cpp src/AltZir_AbstractClass_cpp/AltZir.cpp src/AltZir_D.cpp src/AltZir_DIII.cpp src/AltZir_C.cpp src/AltZir_CI.cpp src/Quantum_chaotic_billiard.cpp src/Auxiliary_Functions.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Conductance_Channels.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Conductance_Gamma.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Concurrence_Gamma.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Ress.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Gamma.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Fixed_Base.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Conductance_Energy.cpp src/AltZir_AbstractClass_cpp/Run_Simulation_Conductance_Energy_Gamma.cpp -o AltZir

clean:
	rm AltZir
