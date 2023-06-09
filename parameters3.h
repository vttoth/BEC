// BEC3P.cpp : Self-gravitating Bose-Einstein condensate
//
// This version dated 2013/08/03.


// Set this to enable OPENCL (experimental;
//                            GPU needed, negligible gain in performance)
//#define USECL

#ifdef USECL
#define Float float
#else
#define Float double
#endif
#ifndef KERNEL
const Float pi = 4 * atan((Float)1);
#endif

#define SQ(x) ((x)*(x))

// User-configurable parameters begin here

// Self-gravitating or trapped condensate?
#define GRAV

// Grid size
#define Nx 80
#define Ny 80
#define Nz 80

// A comment on units: The Gross-Pitaevskii equation is solved in a
// "dimensionless" form where the particle mass (M) and the reduced Planck
// constant (hbar) are both set to 1. This leaves us free to define one of
// the three fundamental units.
//
// For example, consider a Bose-star with particles weighing 5.6e-11 eV,
// or 1e-46 kg. This is our unit of mass. The condition that hbar = 1 then
// defines the ratio between area ([L]^2) and time ([T]) such that that 1
// m^2/s = 1e12 [L]^2/[T]. Choosing the unit of length to be 1 km = 1000 m
// implies that time must be measured in microseconds, [T] = 1 us = 1e-6 s.
//
// Angular velocity (omega) is measured in radians/[T]. Therefore, with the
// above unit of time, omega = 1e-2 amounts to 10,000 radians or about 1591
// revolutions per second. For a star with a 50 km radius, this implies a
// rotational velocity of 500,000 km/s, which is superluminal. This value is
// used therefore only for program testing and not meant to represent a "real"
// Bose star.
//
// The combined Gross-Pitaevskii-Poisson system remains invariant under the
// following transformation: t -> l^2 t, r -> lr, V -> V / l^2, mu -> mu / l^2,
// c -> c / k^2 / l^2, G -> G / k^2 / l^4, psi -> k psi, where k and l are
// arbitrary scaling constants. In particular, this allows us to scale
// simultaneously the particle number N ~ psi^2, and the gravitational constant
// G as follows: N -> k^2 N, G -> G / k^2.
//
// By way of example, if the Bose star above was of 1 solar mass (~ 2e30 kg),
// it would contain N = 2e76 particles. The gravitational constant in the units
// chosen above is 6.67e-11 m^3/kg/s^2 = 6.67e-78 [L]^3/[M][T]^2. Choosing
// k = 1e-38 allows us to rescale these quantities such that N = 2, G = 0.0667.
// This is especially useful if single-precision arithmetic is used (e.g., when
// using OPENCL on a GPU.)
//
// For a galactic dark matter cloud consisting of ultralight bosons, very
// different values are needed. Assuming a boson mass of 1e-23 eV ~ 1.78e-59 kg
// and using the [L] = 1 kpc ~ 3.09e19 m as our unit of length results in
// [M] = 1.78e-59 kg and [T] ~ 1.66e14 s ~ 5.26 Myr as our unit of time. In
// these units, the gravitational constant is G = 1.11e-99 [L]^3/[M]/[T]^2,
// and a 1e11 solar mass halo contains N = 1.12e100 particles. Rescaling using
// k = 1e-50 yields G = 11.1 and N = 1.12. In this case, angular velocity is
// measured in radians / [T], so omega = 1 corresponds to 0.03 revolutions per
// million years. For a cloud with a radius of 50 kpc, this is equivalent to
// a nonrelativistic surface velocity of roughly 3% of the speed of light.
//
// It is also possible to simulate a microscopic 87Rb condensate in a trapping
// potential. For this condensate, the scattering length is a = 5.82 nm, and
// the particle mass is 1.445e-25 kg. Choosing [T] = 0.1/2pi s ~ 0.016 s yields
// a unit of length [L] = 3.42e-6 m. In these units, a = 0.0017.
// 
// For a meaningful simulation, a trap with a radius of 50 L ~ 1.71e-4 m can 
// used, containing N = 21000 particles. The condensate radius was 6.5 L, and
// the simulation step is 0.01 T. A simulation lasting 100,000 [T] therefore
// amounts to ~16 seconds.
//
// Lastly, the value of the scattering parameter, a, given by the formula
// a = G * SQ(R/pi) is consistent with the Thomas-Fermi approximation (i.e.,
// this is the formula used by Boehmer and Harko that relates the values of a
// and R). It is not necessary to use this value; other values may be
// physically justified. In restored units, the actual equation is
// hbar^2 a / m = Gm^2 (R/pi)^2, or R = pi sqrt(hbar^2 a / Gm^3).

// Physical size of simulation volume in units of [L]
const Float xl = -160.0f, yl = -160.0f, zl = -160.0f;
const Float xr = 160.0f, yr = 160.0f, zr = 160.0f;

#ifndef KERNEL
// Simulation parameters
const Float tau = 10;					// Time step (units of [T])
const int time_n = 100000;				// Number of iterations to run
const Float G = 0.0667;					// Newton's constant (may be scaled)
const Float N = 4;						// Particle number (may be scaled)
const Float R = 80.0;					// Size of initial condensate (in [L])
const Float a = 0.5 * G * SQ(R/pi);			// Scattering length (TF default)
const Float c = 4 * pi * a;				// BEC interaction coupling strength
const Float ex = 0.0;					// Softening parameters
const Float ey = 0.0;
const Float ez = 0.0;
const Float omega0 = 1e-4;				// Initial angular velocity (in rad/[T])
const Float gamma0 = 0.0;				// Softening parameter
const int despin_n = 10000;				// When to stop spinning the condensate

// Iteration tolerances
const Float tolGPE = 1e-4;				// GPE nonlinear term iteration
const Float tolPSN = 1e-4;				// Poisson relaxation method iteration
const Float tolREL = 1e-6;				// Imaginary time system relaxation

// Output control
const int nstep0 = 1;	// number of steps of initial transient without output
const int nstep1 = 2;	// every how many steps energy is output
const int nstep2 = 20;	// every how many steps contour plot is output
#endif
