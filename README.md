# CFD
This is a DNS code for a flat channel with cylindrical shaped roughness developed by Dr Ricardo Garcia-Mayoral and his research team. 

The Navier-Stokes equations are advanced in time using the method in Kim, J. & Moin, P. 1985 Application of a Fractional-Step Method to Incompressible Navier-Stokes Equations. J. Comp. Phys. 59, 308–323.

This code uses MPI so a parallel computing architecture is needed. 

The main program is inside littleharsh.f90 and the parameters for the code are set inside input.in. 
