
RESONANCES BETWEEN MASSIVE BODIES
http://www.fisica.edu.uy/~gallardo/atlas/plares.html

1) TWO PLANETS. If you are interested in mutual resonances between two planets in a one star + 2 planets system, the most
the simple version is

version 3: plaresPoincare.f

It works in the Poincaré system, but in the case of planets similar to Earth or super-Earths, the Poincaré system is almost the same as
an astrocentric system, so you can get good estimates of equilibrium points, libration periods, and widths.


2) SEVERAL PLANETS. If you are interested in several mutual resonances between two planets in a one star + multi-planet system,
the best option is versions 1 or 2 (in Jacobi system).

They are almost the same code in two versions:
version1: you need the input file "atlaspr.inp" with the planetary system data.
version2: planetary system data incorporated in the code.
Both codes calculate an "atlas" of resonances, that is, resonance widths in space (a,e) for a given orbital inclination.


3) BINARY. If you are interested in studying a planet around a binary star system, you must use version 1 or 2 (the Jacobi system is required).


4) SATELLITES. If you are interested in inter-satellite resonances around a planet with negligible J2, you can test the codes at your own risk.


5) HAMILTONIAN. If you want to study the level curves of the Hamiltonian the code is hamiltplares.f.





For very weak resonances the codes can produce erroneous outputs, you can contact us and we will try to help you:
 
gallardo@fisica.edu.uy
cbeauge@unc.edu.ar
cristian.giuppone@unc.edu.ar

Reference: Gallardo, Beauge, Giuppone, 2021, 
A&A 646, A148 (2021)
"Semianalytical model for planetary resonances
Application to planets around single and binary stars"
