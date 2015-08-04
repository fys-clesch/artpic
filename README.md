# artpic
A ray tracing program in C

Copyright 2013 Clemens Schäfermeier, clemens ( at ) fh-muenster.de

    This file is part of artpic.

    artpic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    artpic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with artpic.  If not, see <http://www.gnu.org/licenses/>.

artpic is an open source project for ray tracing.
It can be used to simulated scattering effects where geomtrical effects are governing the system. Rays are treated as parts of plane eletric waves with a complex amplitude and a polarisation.
So far, the only primitive is a sphere, but complex system can be composed of intersecting spheres.

The program relys on  freeglut, an OpenGL context handler to make it portable over a variety of operating systems, for 3D data visualisation.
OpenMP is used to trace rays throgh the system in parallel.
It is convenient to have gnuplot installed for exporting results to plots.
