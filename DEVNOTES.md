DEVELOPMENT NOTES
-----------------

Using CELLS, I could make solar sail classes that automatically
propogate the trajectory.

By specifying various up front parameters, compose the equations of
motion to propagate, then simply use RKA from BLD-ODE to solve them.
The composition process includes specifying things like sail
properties, force models, coordinates, gravitating bodies, and
luminous bodies.  Seems like this would work best as a closure with
accessible properties.  Thinking how to do this with CLOS...  CLOS can
store a closure of the equations of motion, and a wrapper function
around RKA can bind the sail to be accessible to it.

Implementation of the Kustaanheimo-Stiefel equations of motion is a
goal. The work I need to do is to convert the prior quaternion
representations to geometric algebra. Hestenes went only so far,
falling short of a perturbed orbital element model. I want to get this
working before starting with other architectural improvements.

For the time being, I'm defining a *sc* parameter that stores problem
information in a SAIL class object, which is accessed by slot within
functions for equations of motion, acceleration, and pointing. *sc*
can be redefined to specific instances so that different problems can
be solved.
