# Quick start:

gvpp needs a single configuration file called <boatname>.geom and describing
the boat in terms of its geometric parameters. An example is provided in
asso99.geom, and contains the description of each parameter. The vpp can then
be run from the command line as

'python vpp.py <boatname> <windSpeed> [<windSpeed> <windSpeed> ...]'


where windSpeed is the true wind speed at which the boat performance are
computed. At least one, and as many as you want, with speeds must be provided.

# GVPP 

gvpp is a Velocity Prediction Program for sailing boats. It is designed to
calculate the velocity of the boat under different wind conditions.  gvpp is
based on the DSYHS for the hull's hydrodynamic and the 1980 Hazen's model for
the sail's aerodynamics. The solution is obtained by solving a constrained
maximization problem whose variable are the boat speed *V*, the heel angle
*phi*}, the crew arm *b*, and a factor *F* taking account of the eventual
flattening of the sails. The variable to be maximized is *V*.  Constraints are
applied on the four variables, which are required to assume a value between
user specified minimum and maximum, the equilibrium equations for the forces
along the longitudinal direction and for the heeling and righting moment.


