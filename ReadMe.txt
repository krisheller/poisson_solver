This project is a discretized solver for poisson's equation for gravity.
It seeks to achieve a second-order accurate solution to this equation by
starting on a test problem and converging to this solution through convergence
testing.  Lapack is used to solve the linear systems created by the main program problem2d.f08

The Main folder holds files meant to be edited while the lapack-3.8.0 folder
just contains lapack and should not be edited at all.

Author:   Kris Heller
Email:	  kristoferheller@gmail.com
