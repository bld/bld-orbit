DEVELOPMENT NOTES
=================

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

More Notes 3/11/2013
--------------------

BLD-ODE allows a state to be an arbitrary collection of objects so
long as they have some form of state arithmetic defined on them. By
default, they have hash tables & vectors defined. I just added
lists. Problem is, I need to reach bach & access the parameters of the
spacecraft that the state variable is referencing. One way is with an
additional slot in a state class that access the body calling
it. Trouble is, I'm creating these classes all over the place, and it
would be laborious to repopulate it. Another way is to make it a class
allocated slot, in which case I need to define a new state class for
each and every state variable. Makes the *SC* variable look downright
clean.

3/12/2013
---------

So, the pieces I have are a geometric algebra library and Runge Kutta
adaptive stepsize ODE solver. Now, I need a way to build up a state
and equations of motion and solve them. The ODE solver will solve
states of a variety of types.

Suggestion: for the state variable and equation of motion, define a
closure that contains the problem data, like spacecraft,
bodies, and other parameters of the problem.