The code as packaged here contains a Makefile that has been tested
using gcc 4.5.2 on a 64-bit Linux host. The C++ code also compiles
as a Visual Studio 2010 console application project on Windows 7
64-bit. This archive contains the following files:

00readme.txt    - This file
Makefile        - UNIX/Linux makefile
stdafx.h        - C/C++ header file (mainly for Visual Studio)
parameters3.h   - C/C++ header file with configurable parameters
bec3p.c++       - C/C++ program file
bec3p.cl        - OPENCL code (experimental)
animate.sh      - Shell (bash) script to create animated GIFs
video.sh        - Shell (bash) script to convert animations into MP4

NB: Support for the FORTRAN version has been discontinued, due to the
superior performance of the C++ version.

WARNING: This program produces a large number of output files. It may
be a good idea to run this program in an empty folder to avoid
clobbering existing data.

Compilation using Makefile provided
-----------------------------------

Adjust run-time parameters in parameters.h as desired,
then compile the code on Linux using:

make all    - compiles both the C++ and Fortran versions
make bec3pc - compiles the C++ version

Test run
--------
To run the program enter,
  ./bec3pc  - C++ executable

The sample parameter file (parameters.h) is configured to simulate
a BEC star or stellar core using 100000 iterations and a grid size
of 80x80x80. With these settings, the program takes several hours
to run on a typical computer.

OPENCL
------

An experimental OPENCL version is also provided. To enable OPENCL,
set the OPENCL macro in parameters.h. OPENCL has only been tested
on a Windows platform, and the performance gain of this implementation
was modest.

Visualisation
-------------

The output files from the program can be converted to GIF and MPEG-4
formats using the scripts provided; more details below.

./animate.sh -o bec3p_out.gif

The bec3_out.gif file can be viewed in a browser window or using the
'animate' utility of the ImageMagick suite, which is preinstalled in
many Linux distributions.

./video.sh -i bec2p_out.gif -o bec3p_out.mp4 

The bec30_out.mp4 can be viewed in Windows Media Player or VLC media
player.



Detailed description
--------------------

Our self-gravitating BEC code has several components. A stand-alone
executable compiled from C++ (bec3p.c++) source performs the simulation.
User-definable parameters are in a header file (parameters3.h). A change
in parameters necessitates recompiling the executables, but this is
quickly accomplished using the supplied make file (Makefile). The
simulation code creates data files at set intervals, specifically
formatted to facilitate plotting using the gnuplot utility. A UNIX style
shell script (animate.sh) is provided that uses gnuplot to create
animations in the form of animated GIF files. Finally, yet another shell
script (video.sh) can be used to convert these GIF files into MPEG-4
video animations.

For additional details, see arXiv:1207.5249.

The structure of the main source file (bec3p.c++) is as follows:

* The main program (_tmain function in the C++ imlementation)
  initializes parameters and executes the main simulation loop.

* The function init creates the initial distribution. This is the
  only user-modifiable subroutine in the main program file; in the
  test implementation, it creates a naive initial distribution
  inspired by the Lane-Emden equation.

* The function get_U computes the potential term (without rotation) in
  the Gross-Pitaevskii equation, as given in Eq. (3).

* The functions calc_rhs_x, calc_rhs_y and calc_rhs_z compute the
  right-hand side of the Gross-Pitaevskii equation (1).

* The function solve_x, solve_y and solve_z solve for the left-hand side
  of the Gross-Pitaevskii equation according to the scheme represented
  by Eq. (4).

* The function thomas solves a tridiagonal system of linear equations
  using the Thomas algorithm.

* The function get_normsimp obtains a numerical approximation of the
  integral of the mod square of the wave function using Simpson's
  formula.

* The functions movieX, movieY and movieZ generate output after a set
  number of iterations (controlled by the variable nstep2) that is set
  in parameters3.h or parameters3.f90 in a format that is suitable for
  use with gnuplot.

* The function get_phi computes the gravitational potential
  corresponding to the current values of \psi using the relaxation
  method.

* The function get_density computes the density |\psi|^2.


Output files produced by the program include the following:

* The files densZnnnnnnn.dat contain values of density |\psi|^2, sorted
  and grouped to facilitate plotting with gnuplot in the plane
  perpendicular to the Z-axis (i.e., the XY plane). Similarly, the files
  densXnnnnnnn.dat and densYnnnnnnn.dat contain density values sorted
  and grouped for plotting in the plane perpendicular to the X and
  Y-axis.

* The files gravZnnnnnnn.dat, gravXnnnnnnn.dat and gravYnnnnnnn.dat
  contain values of the gravitational potential, arranged similarly.

* The files phasZnnnnnnn.dat, phasXnnnnnnn.dat and phasYnnnnnnn.dat
  contain values of the complex phase of \psi, again arranged as before.


The shell script animate.sh provides simple visualization. By default,
this script create an animated GIF of density values, in the plane
perpendicular to the Z-axis and across the origin, with the color range
normalized to the 1000th iteration step. Running

  animate.sh -h

prints a brief help message that lists all the options of this script;
in particular, the option -p that allows the user to change the plane
for visualization and the option -t that allows the user to change the
plot type. The script can also output an animated plot of the estimated
Keplerian rotational velocity (option -t vrot), for circular orbits
centered around the origin in the XY plane.

Finally, the shell script video.sh can be used to convert animated GIF
output into MPEG-4 video. For this script to work, the utilities ffmpeg
and mplayer must be installed on the workstation where the script is
being used.

