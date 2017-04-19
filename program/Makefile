mean : energy
	./mean.sh

energy : initial_configuration
	gfortran pathintegral2.f90 mtfort90.f90 -o pathintegral2
	./pathintegral2
initial_configuration : create_input
	gfortran initial_conf_interaction.f90 mtfort90.f90 -o initial_conf_interaction
	./initial_conf_interaction

create_input	: input_parameters.txt
	./create_input.sh
