# Basic-Molecular-Dynamics-Simulation                             By Sudheer Kumar Peddathimmareddy.
My first attempt at writing a code to perform Molecular Dynamics simulations using basic Newton Laws for Argon ideal gas.
All the modules and functions used for this MD simulation is in one single file.

How to run:
1.  g++ MD_simulation.cpp
2.  ./a.out

This is the very basic way of running this code.

What to expect? (NOTE: Each step generates a files with containing relevant data.)
1.  The code is pre-configured for 64 Argon atoms.
2.  It generates:
    a.  Random intial positions in x, y, z format along a file containing redundant data about distances between atoms (for debugging purpose).
    b.  Random intial velocities in x, y, z format i.e. velocity of a given atom in x direction, y direction and z direction.
4.  After intialization the Molecular Dynamics simulation starts based on the intial velocities and positions according to Velocity - Verlet algorithm.
5.  Results from each step is recorded in the files generated (with self-explanatory names) for example, Kinetic Energies.txt contains information about how KE is changing for all the atoms with each step, etc.
6.  The results, i.e. changes in positions of Ar atoms for each step, from the simulation can be viewed by opening "Final Positions.xyz" in any molecular visualization program supporting .xyz format.

Future work:
1.  Make the code much more user friendly by sub-dividing the Code into modules and such.
2.  Make the code faster and more accurate.
