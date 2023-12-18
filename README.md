# Hydrosolver1D
This is a code submitted as an exercise for the numerical hydrodynamics course at uni tuebingen. It was meant as an one week assignment and therefore should be seen as a showcase of the numerical method and the physics involved and not as a ready-to-use code. 
It's a Numerical Hydrodynamics code to solve the 1DEuler Equations using a Finite Volume approach. One can choose between gaussian or shock tube initial conditions, different EOSs and numerical viscosity. It uses reflecting boundary conditions.

## Initial conditions
1. Gaussian
2. Shock Tube


## Code Manual
To use the code one must first compile the main.cpp file and run the resulting binary. The parameters of the simulations can be set in the main function as constructer attributes of the simulation class. I chose this approach since the exercise required us to run the simulation a few times with different parameters. There are no dependencies on any library except the std lib. 

The simulation results can be plotted via a Jupyter NB script and can be subsequently animated using ffmpeg in a sh script. Therefore downloading ffmpeg is required to get working animations.
