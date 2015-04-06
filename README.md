bld-orbit
=========

bld-orbit is an orbital mechanics in Common Lisp. This library uses
the orbit equations of Kustaanheimo and Stiefel instead of the
traditional Kepler orbital elements and their singularities. Orbits
are represented by the geometry of the orbit plane (using geometric
algebra), and time is scaled by the distance from the orbit center to
make integration more uniform for highly elliptical orbits.

<TBD> Better explanation on what all this is. The included document
"Perturbed Kepler Orbital Element Problem in Geometric Algebra" goes
into a mathematical derivation, but is hardly an introduction.

Comparison with General Mission Analysis Tool (GMAT)
----------------------------------------------------

This library is tested against NASA's open source GMAT program,
version R2013a, to demonstrate that it is capable of the accuracy
needed for real mission analysis. GMAT R2013a was the staging release
for R2013b (an unreleased version) that was deemed ready for
operational use on the Advanced Composition Explorer (ACE) mission.

### Test scripts

Load the file "gmat/Test01.script" to load the first test, which is to
propagate a spacecraft for 1 year, starting from Earth's position and
velocity, with only solar gravity.

The report files require an absolute path, and Test01.script is saved
with the relative path. To set the absolute path:

1. Run GMAT R2013a
2. Open "gmat/Test01.script"
3. Open the "Output" folder in the "Resources" tab.
4. Open "ReportFile01Sun" (double-click)
5. Click "Browse" and select the "gmat/ReportFile01Sun.txt"
6. Repeat process for "ReportFile01SSB".
7. Click "Run" to regenerate report files.
8. Click "Output" tab, and double-click report files to look at the results in GMAT.

Other libraries
---------------

### bld-ga & bld-e3

The orbits are described by geometric algebra objects that represent
the plane and planar velocity of the orbit. These are implemented in
bld-ga (geometric algebra) and bld-e3 (geometric algebra of Euclidean
three-dimensional space).

### bld-ode

bld-orbit uses bld-ode's adaptive step-size Runge-Kutta 7-8th order
integrators. There are versions of the integrator that end at a
specified integrator time, or other stopping condition & value.

### cl-spice

cl-spice is a library that accesses precision orbit and other data of
celestial bodies and spacecraft using JPL's SPICE system.

Usage
-----

See "examples.lisp" for some Kepler orbit test cases.

<TBD>

