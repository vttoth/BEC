// BEC3P.cpp : Self-gravitating Bose-Einstein condensate
//
// Copyright (c) 2023 Eniko J. M. Madarassy and Viktor T. Toth
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// This version dated 2013/08/05.

#include "stdafx.h"

#include "parameters3.h"

#ifdef USECL
#include <CL/opencl.h>
#include <SDKCommon.hpp>
#include <SDKFile.hpp>
#endif

#pragma warning(disable:4996)

using namespace std;

const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1);

#define NRMN 100
#define NRMC 1000
#define NRMP 10000
#define NTOL 1e-6

complex<Float> eye, dt;
Float t, dx, dy, dz, idx2, idy2, idz2;

complex<Float> psi[Nn];
complex<Float> psi_n[Nn];
Float dens[Nn];
Float phi[Nn];
Float UU[Nn];
Float res[Nn];

#define ijk(i,j,k) ((((i) * (Ny + 1)) + (j)) * (Nz + 1) + (k))

#define psi(i,j,k) (psi[ijk(i,j,k)])
#define psi_n(i,j,k) (psi_n[ijk(i,j,k)])
#define density(i,j,k) (dens[ijk(i,j,k)])
#define phi(i,j,k) (phi[ijk(i,j,k)])
#define U(i,j,k) (UU[ijk(i,j,k)])
#define res(i,j,k) (res[ijk(i,j,k)])

// Forward declarations
Float init(int, int, int);
Float fermi(Float, int, int, int);
Float get_normsimp();
void get_cond(Float &, Float &, Float &, Float, Float &);

void get_phi();
void get_Vtr();
void get_U(Float);
void calc_rhs_x(complex<Float>, complex<Float>,
				complex<Float>, complex<Float>);
void solve_x(complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>);
void calc_rhs_y(complex<Float>, complex<Float>,
				complex<Float>, complex<Float>);
void solve_y(complex<Float>, complex<Float>,
			 complex<Float>, complex<Float>);
void calc_rhs_z(complex<Float>, complex<Float>, complex<Float>);
void solve_z(complex<Float>, complex<Float>, complex<Float>);
Float energy(Float, FILE*);
void movie(int);
void thomas(complex<Float> *, complex<Float> *, complex<Float> *,
			complex<Float> *, int m);
void get_density();

#ifdef USECL

cl_mem clpsi, clpsi_n, clU, cldensity, clphi, clres;
cl_kernel kernelNorm_psi, kernelGet_U, kernelGet_res;
cl_kernel kernelCalc_rhs_x, kernelCalc_rhs_y, kernelCalc_rhs_z;
cl_kernel kernelSolve_x, kernelSolve_y, kernelSolve_z;
cl_command_queue command_queue;
#define GROUP_SIZE 256
size_t globalThreads = (((Nn + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE);
size_t localThreads = GROUP_SIZE;
cl_context context = 0;
cl_program program = 0;

bool initializeCL()
{
	// OpenCL structures
	cl_int status;
	cl_platform_id platform;
	cl_uint num_platforms;
	cl_device_id device_id;
	cl_uint num_devices;
	cl_context_properties cps[] = {CL_CONTEXT_PLATFORM, 0, 0};
	streamsdk::SDKFile includeFile;
	streamsdk::SDKFile kernelFile;	// create CL program using the kernel source
	std::string includePath = "parameters3.h";
	std::string kernelPath = "BEC3P.cl";
	const char * source;
	size_t sourceSize;
	const char szCLflags[] = "";
	//const char szCLflags[] = "-cl-mad-enable -cl-no-signed-zeros "
	//						   "-cl-unsafe-math-optimizations "
	//						   "-cl-finite-math-only -cl-fast-relaxed-math";
	char szCLopts[sizeof(szCLflags) + MAX_PATH + 10];
	char szCLpath[MAX_PATH];

	if (clGetPlatformIDs(1, &platform, &num_platforms) != CL_SUCCESS ||
		num_platforms != 1 ||
		clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device_id,
					   &num_devices) != CL_SUCCESS || num_devices != 1)
	{
		printf("Failed to initialize OpenCL.\n");
		return false;
	}

	cps[1] = (cl_context_properties)platform;
    context = clCreateContextFromType(cps, CL_DEVICE_TYPE_ALL, NULL, NULL,
									  &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to create CL context.\n"); context = 0; return false; }
	command_queue = clCreateCommandQueue(context, device_id, 0, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to create CL command queue.\n"); return false; }

	GetShortPathNameA(_getcwd(NULL, 0), szCLpath, sizeof(szCLpath));
	for (char *p = szCLpath; *p != 0; p++) if (*p == '\\') *p = '/';
	strcpy(szCLopts, "-I ");
	strcat(szCLopts, szCLpath);
	strcat(szCLopts, " ");
	strcat(szCLopts, szCLflags);

	source = "#include \"BEC3P.cl\"\n";
	sourceSize = strlen(source);
	program = clCreateProgramWithSource(context, 1, &source, &sourceSize,
										&status);

	if (status != CL_SUCCESS)
		{ printf("Failed to create CL program.\n"); program = 0; return false; }

	if ((status = clBuildProgram(program, 0, NULL, szCLopts, NULL, NULL))
		!= CL_SUCCESS)
	{
		printf("Failed to build CL program: %d.\n", status);
		char buffer[4096];
		size_t length;
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
							  sizeof(buffer), buffer, &length);
		printf("%s\n", buffer);
		return false;
	}

	kernelGet_res = clCreateKernel(program, "get_res", &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to create CL kernel get_res.\n"); return false; }

	kernelCalc_rhs_x = clCreateKernel(program, "calc_rhs_x", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel Calc_rhs_x.\n"); return false; }

	kernelCalc_rhs_y = clCreateKernel(program, "calc_rhs_y", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel Calc_rhs_y.\n"); return false; }

	kernelCalc_rhs_z = clCreateKernel(program, "calc_rhs_z", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel Calc_rhs_z.\n"); return false; }

	kernelSolve_x = clCreateKernel(program, "solve_x", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel solve_x.\n"); return false; }

	kernelSolve_y = clCreateKernel(program, "solve_y", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel solve_y.\n"); return false; }

	kernelSolve_z = clCreateKernel(program, "solve_z", &status);
    if(status != CL_SUCCESS)
		{ printf("Failed to create CL kernel solve_z.\n"); return false; }

	clpsi = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
						   sizeof(psi), psi, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate clpsi\n"); return false; }

	clpsi_n = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
							 sizeof(psi_n), psi_n, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate clpsi_n\n"); return false; }

	cldensity = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
							   sizeof(dens), dens, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate cldensity\n"); return false; }

	clphi = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
						   sizeof(phi), phi, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate clphi\n"); return false; }

	clU = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
						 sizeof(UU), UU, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate clU\n"); return false; }

	clres = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
						   sizeof(res), res, &status);
	if (status != CL_SUCCESS)
		{ printf("Failed to allocate clres\n"); return false; }

	size_t sz;
	char inf[10000];
	memset(inf, 0, sizeof(inf));

	clGetPlatformInfo(platform, CL_PLATFORM_PROFILE, sizeof(inf) - 1, inf,
					  NULL);
	printf("\n\n\nPROFILE: %s\n", inf);
	clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(inf) - 1, inf, NULL);
	printf("NAME: %s\n", inf);
	clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(inf) - 1, inf, NULL);
	printf("VENDOR: %s\n", inf);

	clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS, sizeof(inf) - 1, inf,
					  NULL);
	printf("PLATFORM EXTENSIONS: %s\n", inf);

	clGetDeviceInfo(device_id, CL_DEVICE_PROFILE, sizeof(inf), inf, NULL);
	printf("DEVICE PROFILE: %s\n", inf);

	clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(inf), inf, NULL);
	printf("DEVICE NAME: %s\n", inf);

	clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(inf), inf, NULL);
	printf("DEVICE VENDOR: %s\n", inf);

	clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, sizeof(inf), inf, NULL);
	printf("DEVICE EXTENSIONS: %s\n", inf);

	clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(sz), &sz,
					NULL);
	printf("WORKGROUP SIZE: %d\n", sz);

	printf("====================================================\n\n\n");

	return true;
}

void finalizeCL()
{
	if (program) clReleaseProgram(program);
	if (context) clReleaseContext(context);
}
#endif

//**********************************************************************
void loopdetect(Float *nrma, Float norm, const char *pdir, int &nrmc)
{
	int i, j;
	int nrmn = nrmc++;

	if (nrmc > NRMC)
	{
		fprintf(stderr, "No convergence after %d "
				 "iterations in the %s direction.\n", nrmc, pdir);
		exit(1);
	}

	if (nrmn >= NRMN)
	{
		memmove(nrma, &(nrma[1]), sizeof(nrma) - sizeof(nrma[0]));
		nrmn = NRMN - 1;
	}
	nrma[nrmn] = norm;

	for (i = nrmn - 1; i >= 0; i--)
	{
		if (fabs(norm - nrma[i]) < fabs(NTOL * norm))
		{
			fprintf(stderr, "Oscillatory loop in the %s direction is "
						 "detected, solution is not convergent.\n", pdir);
			for (j = i; j <= nrmn; j++)
				fprintf(stderr, "norm[%d]=%lg\n", j, nrma[j]);
			exit(1);
		}
	}
}

//**********************************************************************
int _tmain(int argc, _TCHAR* argv[])
{
	complex<Float> xi;
	Float omega, gamma, mu;
	Float norm, norm0;
	Float E0, E1 = 0, E2 = 0;
	int i, j, k, itime, ktime;
	FILE *file21, *file22, *file23, *file24, *file33, *file34, *file41;
	Float nrma[NRMN];
	int nrmc;
	complex<Float> foo1X, foo1XS, foo1Y, foo1YS, foo1ZS;
	complex<Float> foo2XS, foo2YS;
	complex<Float> foo3X, foo3XS, foo3Y, foo3YS, foo3Z, foo3ZS;
	complex<Float> foo4;
	complex<Float> foo5X, foo5Y, foo5Z;

#ifdef USECL
	initializeCL();
#endif

	// Fixed parameters

	eye = complex<Float>(0, 1);
	dx = (xr - xl) / Nx;
	dy = (yr - yl) / Ny;
	dz = (zr - zl) / Nz;

	idx2 = 1 / dx / dx;
	idy2 = 1 / dy / dy;
	idz2 = 1 / dz / dz;

	norm0 = N;
     
	// Initial conditions

	// Zero arrays
	memset(phi, 0, sizeof(phi));
	memset(psi, 0, sizeof(psi));
	memset(psi_n, 0, sizeof(psi));
	memset(UU, 0, sizeof(UU));
	memset(res, 0, sizeof(res));
	memset(dens, 0, sizeof(dens));

	dt = complex<Float>(0, -tau);
	omega = 0;
	gamma = 0;
#ifdef GRAV
	mu = 0;
#else
	mu = 0.5 / SQ(R);
#endif

	xi = (eye + gamma) / (gamma * gamma + 1);

	// Helper variables that need to be computed only once
	foo1X = eye * omega * xi * dt / (4 * dx);
	foo1XS = -xi * dt * idx2 / (Float)4;
	foo1Y = eye * omega * xi * dt / (4 * dy);
	foo1YS = -xi * dt * idy2 / (Float)4;
	foo1ZS = -xi * dt * idz2 / (Float)4;
	foo2XS = eye * omega * xi * dt / (4 * dx);
	foo2YS = eye * omega * xi * dt / (4 * dy);
	foo3X = xi * dt * idx2 / (Float)4;
	foo3XS = (Float)1 + xi * dt * idx2 / (Float)2;
	foo3Y = xi * dt * idy2 / (Float)4;
	foo3YS = (Float)1 + xi * dt * idy2 / (Float)2;
	foo3Z = xi * dt * idz2 / (Float)4;
	foo3ZS = (Float)1 + xi * dt * idz2 / (Float)2;
	foo4 = xi * dt / (Float)6;
	foo5X = (Float)1 - xi * dt * idx2 / (Float)2;
	foo5Y = (Float)1 - xi * dt * idy2 / (Float)2;
	foo5Z = (Float)1 - xi * dt * idz2 / (Float)2;

	// Initial state

	file34 = fopen("psi34.dat", "w");
	t = 0.0;
	itime = 0;
	ktime = 0;
	for (i = 0; i <= Nx; i++)
	{
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
		{
#ifdef GRAV
			Float rho = init(i, j, k);
			Float x2 = SQ(xl + i * dx);
			Float y2 = SQ(yl + j * dy);
			Float z2 = SQ(zl + k * dz);
			Float r = sqrt(x2 + y2 + z2);

			psi(i, j, k) = sqrt(rho);
			phi(i, j, k) = (Float)(-G * N / (r > (.25 * dx) ? r : .5 * dx));
#else
			psi(i, j, k) = fermi(mu, i, j, k);
#endif
			fprintf(file34, "%lg %lg %lg %lg\n", xl + i * dx, yl + j * dy,
											zl + k * dz, real(psi(i, j, k)));
		}
		fprintf(file34, "\n");	// For Gnuplot
	}
	fclose(file34);

	// Boundary conditions psi=0
	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
	{
		psi(i, j, 0) = 0;
		psi(i, j, Nz) = 0;
	}
	for (j = 0; j <= Ny; j++)
		for (k = 0; k <= Nz; k++)
	{
		psi(0, j, k) = 0;
		psi(Nx, j, k) = 0;
	}
	for (k = 0; k <= Nz; k++)
		for (i = 0; i <= Nx; i++)
	{
		psi(i, 0, k) = 0;
		psi(i, Ny, k) = 0;
	}

	get_density();

#ifdef GRAV
	get_phi();
#else
	get_Vtr();
#endif

	// get initial U 
	file41 = fopen("erg41.dat", "w");
	get_U(mu);

	norm =  get_normsimp();
	printf("N=%6d, t=%11.4lg, E=%11.4lg, P=%11.4lg\n", 0, 0.0,
			energy(mu, file41), norm);
	fflush(stdout);

	movie(itime);	// Output data for contour plots

	file21 = fopen("psi21.dat", "w");
	file22 = fopen("psi22.dat", "w");
	file23 = fopen("psi23.dat", "w");
	file24 = fopen("psi24.dat", "w");

	// Time loop

	for (itime = 1, ktime = 0; itime - ktime <= time_n; itime++)
	{
		bool bLoop;
		t += real(dt);
		if ((itime - ktime) == despin_n) omega = 0;

		// Find psi'=Rx(psi_old)
		calc_rhs_x(foo1X, foo3X, foo4, foo5X);

		for (bLoop = true, nrmc = 0; bLoop;)// Solve Lx(psi'')=psi' and iterate
		{									// due to nonlinear term in U
			solve_x(foo1XS, foo2XS, foo3XS, foo4);
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "X", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_x required %d interations.\n", nrmc);
	fflush(stderr);
#endif
		// Find psi'''=Ry(psi'')
		calc_rhs_y(foo1Y, foo3Y, foo4, foo5Y);
		for (bLoop = true, nrmc = 0; bLoop;)
		{	// Solve Ly(psi'''')=psi'''
			solve_y(foo1YS, foo2YS, foo3YS, foo4);
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "Y", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_y required %d interations.\n", nrmc);
	fflush(stderr);
#endif
		calc_rhs_z(foo3Z, foo4, foo5Z);		// Find psi'''''=Rz(psi'''')
		for (bLoop = true, nrmc = 0; bLoop;)
		{
			solve_z(foo1ZS, foo3ZS, foo4);	// Solve Lz(psi_new)=psi'''''
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
			Float norm1 = get_normsimp();
			if (fabs(norm1 - norm) < fabs(tolGPE * norm)) bLoop = false;
			else loopdetect(nrma, norm1, "Z", nrmc);
			norm = norm1;
		}
#ifdef SHOW_LOOPS
	fprintf(stderr, "H_z required %d interations.\n", nrmc);
	fflush(stderr);
#endif

		E0 = energy(mu, file41);
		printf("N=%6d, t=%11.4lg, E=%11.4lg, P=%11.4lg\n", itime, t, E0, norm);
		fflush(stdout);

		if (itime > 10 && itime % nstep1 == 0)
		{
			file33 = fopen("psi33.dat", "w");

			for (i = 0; i <= Nx; i++)
			{
				for (j = 0; j <= Ny; j++)
					for (k = 0; k <= Nz; k++)
						fprintf(file33, "%lg %lg %lg\n", xl + i * dx,
												yl + j * dy, density(i, j, k));
				fprintf(file33, "\n");  // For Gnuplot
			}
			fclose(file33);
		}
		if (itime > nstep0 && itime % nstep2 == 0)
		{
			movie(itime);			// Output data for contour plots
		}

		{
			Float x2, y2, alpha1, elz;

			get_cond(x2, y2, alpha1, omega, elz);

			fprintf(file21, "%lg %lg\n", t, alpha1);
			fprintf(file22, "%lg %lg\n", t, elz);
			fprintf(file23, "%lg %lg\n", t, E0);
			fprintf(file24, "%lg %lg\n", t, norm);

			fflush(file21);
			fflush(file22);
			fflush(file23);
			fflush(file24);
		}

		if (real(dt) < 1e-15)
		{

			if (fabs(E0 - E1) < tolREL * fabs(E0) &&
				fabs(2 * E1 - E0 - E2) < tolREL * fabs(E0))
			{
				dt = tau;
				omega = omega0;
				gamma = gamma0;
#ifndef GRAV
				mu = 0;
#endif
				xi = (eye + gamma) / (gamma * gamma + 1);

				// Recompute helper variables
				foo1X = eye * omega * xi * dt / (4 * dx);
				foo1XS = -xi * dt * idx2 / (Float)4;
				foo1Y = eye * omega * xi * dt / (4 * dy);
				foo1YS = -xi * dt * idy2 / (Float)4;
				foo1ZS = -xi * dt * idz2 / (Float)4;
				foo2XS = eye * omega * xi * dt / (4 * dx);
				foo2YS = eye * omega * xi * dt / (4 * dy);
				foo3X = xi * dt * idx2 / (Float)4;
				foo3XS = (Float)1 + xi * dt * idx2 / (Float)2;
				foo3Y = xi * dt * idy2 / (Float)4;
				foo3YS = (Float)1 + xi * dt * idy2 / (Float)2;
				foo3Z = xi * dt * idz2 / (Float)4;
				foo3ZS = (Float)1 + xi * dt * idz2 / (Float)2;
				foo4 = xi * dt / (Float)6;
				foo5X = (Float)1 - xi * dt * idx2 / (Float)2;
				foo5Y = (Float)1 - xi * dt * idy2 / (Float)2;
				foo5Z = (Float)1 - xi * dt * idz2 / (Float)2;
			}
			else
			{
				E2 = E1;
				E1 = E0;
				ktime++;
			}
#if 1
			Float renorm = sqrt(N / norm);
			for (i = 0; i <= Nx; i++)
				for (j = 0; j <= Ny; j++)
					for (k = 0; k <= Nz; k++)
				psi(i, j, k) *= renorm;
			get_density();
#ifdef GRAV
			get_phi();
#endif
			get_U(mu);
#endif
		}

	}	// Close time loop

	fclose(file21);
	fclose(file22);
	fclose(file23);
	fclose(file24);
	fclose(file41);

#ifdef USECL
	finalizeCL();
#endif

	return true;
}

//*********************************************************************
Float init(int i, int j, int k)		// Initial state
{
	Float F, x, y, z, r;

	x = xl + i * dx;
	y = yl + j * dy;
	z = zl + k * dz;
	r = sqrt((1 + ex) * x * x + (1 + ey) * y * y + (1 + ez) * z * z);
	F = 0.0;
//	if (r < R) F = (Float)(N / 4 / pi / (r > 0 ? SQ(r) : SQ(.5 * dx)) / R);

//	if (r < R) F = (Float)N * 27 / 8 / pi * SQ(log(r / R + 1e-29))

//	if (r < R) F = (Float)N * 27 / 8 / pi * SQ(log((r > 0 ? r : .5 * dx) / R))
//					/ (4 * pi * R * R * R / 3);

	if (r > 0) F = sin(r / R) / (r / R);
	else if (r == 0) F = 1;
	if (F <= 0) F = 0;
	else F = (Float)N * pow(F, 6) / (4 * pi * R * R * R / 3);
	return F;
}

//**********************************************************************
Float energy(Float mu, FILE *file41)
{
	Float EK, EU, EI, mdpsisq;
	int i, j, k;
	const static Float dV = dx * dy * dz;

	EK = EU = EI = 0;
	for (k = 1; k < Nz; k++)
		for (j = 1; j < Ny; j++)
			for (i = 1; i < Nx; i++)
	{
		mdpsisq = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
		EI += SQ(mdpsisq);
		EU += (
#ifdef GRAV
			0.5 *		// See Wang, PRD64, 124009
#endif
			phi(i, j, k) - mu) * mdpsisq;
		EK += SQ(real(psi(i + 1, j, k) - psi(i - 1, j, k))) * idx2;
		EK += SQ(imag(psi(i + 1, j, k) - psi(i - 1, j, k))) * idx2;
		EK += SQ(real(psi(i, j + 1, k) - psi(i, j - 1, k))) * idy2;
		EK += SQ(imag(psi(i, j + 1, k) - psi(i, j - 1, k))) * idy2;
		EK += SQ(real(psi(i, j, k + 1) - psi(i, j, k - 1))) * idz2;
		EK += SQ(imag(psi(i, j, k + 1) - psi(i, j, k - 1))) * idz2;
	}
	EI *= (Float)0.5 * c * dV;
	EU *= dV;
	EK *= (Float)0.125 * dV;
	fprintf(file41, "%lg\t%lg\t%lg\t%lg\n", t, EK, EU, EI);
	fflush(file41);
	return EK + EU + EI;
}

//**********************************************************************
void get_U(Float mu)	// Find U
{
	int i, j, k;

	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
		U(i, j, k) = c * density(i, j, k) + phi(i, j, k) - mu;
}

//*********************************************************************
void get_Vtr()	// Trapping potential
{
	int i, j, k;
	const Float d = (Float)0.25 * (SQ(xr - xl) + SQ(yr - yl) + SQ(zr - zl));

	for (k = 0; k <= Nz; k++)
		for (j = 0; j <= Ny; j++)
			for (i = 0; i <= Nx; i++)
	{
		phi(i, j, k) = (Float)0.5 * ((1 + ex) * SQ(xl + i * dx) +
								(1 + ey) * SQ(yl + j * dy) +
								(1 + ez) * SQ(zl + k * dz) - d) / SQ(SQ(R));
	}
}

//*********************************************************************
Float fermi(Float mu, int i, int j, int k)	// Fermi-Thomas initial state
{
	Float x, y, z, r2, R2;
	const Float norm = 15 * N * sqrt(2 * mu) * SQ(mu) / pi;

	x = xl + i * dx;
	y = yl + j * dy;
	z = zl + k * dz;
	r2 = (1 + ex) * SQ(x) + (1 + ey) * SQ(y) + (1 + ez) * SQ(z);
	R2 = SQ(R);
	if (r2 < R2) return (Float)sqrt((0.5 * (R2 - r2)) * norm);
	else return 0.0;
}

//**********************************************************************
void calc_rhs_x(complex<Float> foo1,	// Find psi_n= Rx(psi)
			complex<Float> foo3, complex<Float> foo4, complex<Float> foo5)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelCalc_rhs_x, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelCalc_rhs_x, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelCalc_rhs_x, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelCalc_rhs_x, a++, sizeof(Float), (void *)&dy);
	clSetKernelArg(kernelCalc_rhs_x, a++, 2 * sizeof(Float), (void *)&foo1);
	clSetKernelArg(kernelCalc_rhs_x, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelCalc_rhs_x, a++, 2 * sizeof(Float), (void *)&foo4);
	clSetKernelArg(kernelCalc_rhs_x, a++, 2 * sizeof(Float), (void *)&foo5);

	clEnqueueWriteBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						 NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelCalc_rhs_x, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						0, NULL, NULL);
#else
	Float y;
	complex<Float> foo2;
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		y = yl + j * dy;
		foo2 = foo1 * y;
		for (k = 1; k < Nz; k++)
		{
			psi_n(i, j, k) = (foo3 - foo2) * psi(i - 1, j, k) +
							 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
							 (foo3 + foo2) * psi(i + 1, j, k);
		}
	}
#endif
}

//**********************************************************************
void calc_rhs_y(complex<Float> foo1,	// Find psi_n=Ry(psi)
			complex<Float> foo3, complex<Float> foo4, complex<Float> foo5)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelCalc_rhs_y, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelCalc_rhs_y, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelCalc_rhs_y, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelCalc_rhs_y, a++, sizeof(Float), (void *)&dx);
	clSetKernelArg(kernelCalc_rhs_y, a++, 2 * sizeof(Float), (void *)&foo1);
	clSetKernelArg(kernelCalc_rhs_y, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelCalc_rhs_y, a++, 2 * sizeof(Float), (void *)&foo4);
	clSetKernelArg(kernelCalc_rhs_y, a++, 2 * sizeof(Float), (void *)&foo5);

	clEnqueueWriteBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						 NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelCalc_rhs_y, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						0, NULL, NULL);
#else
	Float x;
	complex<Float> foo2;
	int i, j, k;

	for (i = 1; i < Nx; i++)
	{
		x = xl + i * dx;
		foo2 = foo1 * x;
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
		{
			psi_n(i, j, k) = (foo3 + foo2) * psi(i, j - 1, k) +
							 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
							 (foo3 - foo2) * psi(i, j + 1, k);
		}
	}
#endif
}

//**********************************************************************
void calc_rhs_z(complex<Float> foo3,	// Find psi_n=Ry(psi)
				complex<Float> foo4, complex<Float> foo5)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelCalc_rhs_z, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelCalc_rhs_z, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelCalc_rhs_z, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelCalc_rhs_z, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelCalc_rhs_z, a++, 2 * sizeof(Float), (void *)&foo4);
	clSetKernelArg(kernelCalc_rhs_z, a++, 2 * sizeof(Float), (void *)&foo5);

	clEnqueueWriteBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						 NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelCalc_rhs_z, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						0, NULL, NULL);
#else
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
			for (k = 1; k < Nz; k++)
		psi_n(i, j, k) = foo3 * psi(i, j, k - 1) +
						 (foo5 - foo4 * U(i, j, k)) * psi(i, j, k) +
						 foo3 * psi(i, j, k + 1);
#endif
}

//**********************************************************************
void solve_x(complex<Float> foo1,	// Solve Lx(psi)=psi_n for psi
			 complex<Float> foo2, complex<Float> foo3, complex<Float> foo4)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelSolve_x, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelSolve_x, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelSolve_x, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelSolve_x, a++, 2 * sizeof(Float), (void *)&foo1);
	clSetKernelArg(kernelSolve_x, a++, 2 * sizeof(Float), (void *)&foo2);
	clSetKernelArg(kernelSolve_x, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelSolve_x, a++, 2 * sizeof(Float), (void *)&foo4);
	clSetKernelArg(kernelSolve_x, a++, sizeof(Float), (void *)&dy);

	clEnqueueWriteBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelSolve_x, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						NULL, NULL);
#else
	Float y;
	complex<Float> l[Nx], d[Nx], u[Nx], r[Nx];
	int i, j, k;

	for (k = 1; k < Nz; k++)
		for (j = 1; j < Ny; j++)
	{
		y = yl + j * dy;
		l[1] = foo1 + foo2 * y;
		d[1] = foo3;
		u[1] = foo1 - foo2 * y;
		r[1] = psi_n(1, j, k);
		for (i = 2; i < Nx; i++)
		{
			l[i] = l[1];
			d[i] = d[1] + foo4 * U(i, j, k);
			u[i] = u[1];
			r[i] = psi_n(i, j, k);
		}
		d[1] += foo4 * U(1, j, k);
		thomas(l, d, u, r, Nx - 1);
		for (i = 1; i < Nx; i++) psi(i, j, k) = r[i];
	}
#endif
}

//**********************************************************************
void solve_y(complex<Float> foo1,	// Solve Ly(psi)=psi_n for psi
			 complex<Float> foo2, complex<Float> foo3, complex<Float> foo4)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelSolve_y, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelSolve_y, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelSolve_y, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelSolve_y, a++, 2 * sizeof(Float), (void *)&foo1);
	clSetKernelArg(kernelSolve_y, a++, 2 * sizeof(Float), (void *)&foo2);
	clSetKernelArg(kernelSolve_y, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelSolve_y, a++, 2 * sizeof(Float), (void *)&foo4);
	clSetKernelArg(kernelSolve_y, a++, sizeof(Float), (void *)&dx);

	clEnqueueWriteBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelSolve_y, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						NULL, NULL);
#else
	Float x;
	complex<Float> l[Ny], d[Ny], u[Ny], r[Ny];
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (k = 1; k < Nz; k++)
	{
		x = xl + i * dx;
		l[1] = foo1 - foo2 * x;
		d[1] = foo3;
		u[1] = foo1 + foo2 * x;
		r[1] = psi_n(i, 1, k);
		for (j = 2; j < Ny; j++)
		{
			l[j] = l[1];
			d[j] = d[1] + foo4 * U(i, j, k);
			u[j] = u[1];
			r[j] = psi_n(i, j, k);
		}
		d[1] += foo4 * U(i, 1, k);
		thomas(l, d, u, r, Ny - 1);
		for (j = 1; j < Ny; j++) psi(i, j, k) = r[j];
	}
#endif
}

//**********************************************************************
void solve_z(complex<Float> foo1,	// Solve Ly(psi)=psi_n for psi
			 complex<Float> foo3, complex<Float> foo4)
{
#ifdef USECL
	int a = 0;

	clSetKernelArg(kernelSolve_z, a++, sizeof(cl_mem), (void *)&clpsi);
	clSetKernelArg(kernelSolve_z, a++, sizeof(cl_mem), (void *)&clpsi_n);
	clSetKernelArg(kernelSolve_z, a++, sizeof(cl_mem), (void *)&clU);
	clSetKernelArg(kernelSolve_z, a++, 2 * sizeof(Float), (void *)&foo1);
	clSetKernelArg(kernelSolve_z, a++, 2 * sizeof(Float), (void *)&foo3);
	clSetKernelArg(kernelSolve_z, a++, 2 * sizeof(Float), (void *)&foo4);

	clEnqueueWriteBuffer(command_queue, clpsi_n, CL_TRUE, 0, sizeof(psi_n), psi_n,
						 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, clU, CL_TRUE, 0, sizeof(UU), UU, 0, NULL,
						 NULL);
	clEnqueueNDRangeKernel(command_queue, kernelSolve_z, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, clpsi, CL_TRUE, 0, sizeof(psi), psi, 0,
						NULL, NULL);
#else
	complex<Float> l[Nz], d[Nz], u[Nz], r[Nz];
	int i, j, k;

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++)
	{
		l[1] = u[1] = foo1;
		d[1] = foo3;
		r[1] = psi_n(i, j, 1);
		for (k = 2; k < Nz; k++)
		{
			l[k] = u[k] = foo1;
			d[k] = d[1] + foo4 * U(i, j, k);
			r[k] = psi_n(i, j, k);
		}
		d[1] += foo4 * U(i, j, 1);
		thomas(l, d, u, r, Nz - 1);
		for (k = 1; k < Nz; k++) psi(i, j, k) = r[k];
	}
#endif
}

//**********************************************************************
void get_cond(Float &x2, Float &y2, Float &alpha1,	// Find the mean values
			 Float omega, Float &elz)				// <x^2>, <y^2>, <elz>
{													// and the distortion
													// alpha of the condensate
	complex<Float> px, py, l;
	Float x, y;
	int i, j, k;

	x2 = 0.0;
	y2 = 0.0;
	elz = 0.0;
	l = complex<Float>(0, 0);

	for (i = 0; i <= Nx; i++)
	{
		x = xl + i * dx;
		for (j = 0; j <= Ny; j++)
		{
			y = yl + j * dy;
			for (k = 0; k <= Nz; k++)
			{
				x2 += x * x * (SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))));
				y2 += y * y * (SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k))));
				if(i > 0 && i < Nx)
					px = (psi(i + 1, j, k) - psi(i - 1, j, k)) / (2 * dx);
				else if (i == 0) px = (psi(1, j, k) - psi(0, j, k)) / dx;
				else if (i == Nx) px = (psi(Nx, j, k) - psi(Nx - 1, j, k)) / dx;
				if (j > 0 && j < Ny)
					py = (psi(i, j + 1, k) - psi(i, j - 1, k)) / (2 * dy);
				else if (j == 0) py = (psi(i, 1, k) - psi(i, 0, k)) / dy;
				else if (j == Ny) py = (psi(i, Ny, k) - psi(i, Ny - 1, k)) / dy;
				l += conj(psi(i, j, k)) * eye * (x * py - y * px);
			}
		}
	}
	elz = sqrt(SQ(real(l)) + SQ(imag(l))) * dx * dy * dz;
	alpha1 = omega * ((x2 - y2) / (x2 + y2));
}

//**********************************************************************
void thomas(complex<Float> *l,	// Tridiagonal matrix solver with Thomas
			complex<Float> *d,	// algorithm. The arrays l,d,u and r contain
			complex<Float> *u,	// respectively the subdiagonal, the diagonal,
			complex<Float> *r,	// the superdiagonal and the right hand side.
			int m)				// m is the size of the system. At output the
{								// solution is put in r, while the original
								// right hand side is destroyed.
	int j;
	complex<Float> a;
    
	for (j = 2; j <= m; j++)
	{
		a = -l[j] / d[j - 1];
		d[j] = d[j] + a * u[j - 1];
		r[j] = r[j] + a * r[j - 1];
	}
	r[m] = r[m] / d[m];
	for (j = m - 1; j >= 1; j--)
		r[j] = (r[j] - u[j] * r[j + 1]) / d[j];
}

//***************************************************************************
Float get_normsimp()	// Find the norm using Simpson's rule
{
	static Float normz[Nx + 1][Ny + 1], normy[Nx + 1], norm;
	int i, j, k;

	norm = 0.0;
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
	{
		normz[i][j] = SQ(real(psi(i, j, 0))) + SQ(imag(psi(i, j, 0))) +
					SQ(real(psi(i, j, Nz - 1))) + SQ(imag(psi(i, j, Nz - 1)));
		for (k = 1; k <= Nz - 3; k += 2)
			normz[i][j] += 4 * (SQ(real(psi(i, j, k))) +
								SQ(imag(psi(i, j, k)))) +
						   2 * (SQ(real(psi(i, j, k + 1))) +
								SQ(imag(psi(i, j, k + 1))));
	}
	for (i = 0; i < Nx; i++)
	{
		normy[i] = normz[i][0] + normz[i][Ny - 1];
		for (j = 1; j <= Ny - 3; j += 2)
			normy[i] += 4 * normz[i][j] + 2 * normz[i][j + 1];
	}
	norm = normy[0] + normy[Nx - 1];
	for (i = 1; i <= Nx - 3; i += 2)
		norm += 4 * normy[i] + 2 * normy[i + 1];
	return norm * dx * dy * dz / 27;
}
    
//**********************************************************************
void movieZ(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densZ%07d.dat", itime);
	sprintf(cs, "phasZ%07d.dat", itime);
	sprintf(cu, "gravZ%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (k = 0; k <= Nz; k += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (i = 0; i <= Nx; i += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movieX(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densX%07d.dat", itime);
	sprintf(cs, "phasX%07d.dat", itime);
	sprintf(cu, "gravX%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (i = 0; i <= Nx; i += 2)
		for (j = 0; j <= Ny; j += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movieY(int itime)	// Outputs files
{
	Float x, y, z, phpsi, depsi;
	char ch[17], cs[17], cu[17];
	int i, j, k;
	FILE *file8, *file10, *file11;

	sprintf(ch, "densY%07d.dat", itime);
	sprintf(cs, "phasY%07d.dat", itime);
	sprintf(cu, "gravY%07d.dat", itime);
	file8 = fopen(ch, "w");
	file10 = fopen(cs, "w");
	file11 = fopen(cu, "w");

	// Output files for contour plots

	for (j = 0; j <= Ny; j += 2)
		for (i = 0; i <= Nx; i += 2)
	{
		for (k = 0; k <= Nz; k += 2)
		{
			x = xl + i * dx;
			y = yl + j * dy;
			z = zl + k * dz;
			depsi = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
			phpsi = atan2(imag(psi(i, j, k)) + 1e-15f, real(psi(i, j, k)));
			fprintf(file8, "%lf %lf %lf %lg\n", x, y, z, depsi);
			fprintf(file11, "%lf %lf %lf %lg\n", x, y, z, phi(i, j, k));
			fprintf(file10, "%lf %lf %lf %lg\n", x, y, z, phpsi);
		}
		fprintf(file8, "\n");	// Leave blank line for Gnuplot before next j
		fprintf(file11, "\n");
		fprintf(file10, "\n");
	}

	fclose(file8);
	fclose(file10);
	fclose(file11);
}

//**********************************************************************
void movie(int itime)
{
	movieZ(itime);
	movieX(itime);
	movieY(itime);
}

//**********************************************************************
#ifdef USECL
inline void get_phi_kernel(bool bDir, Float G4pi, Float idx2, Float idy2,
						   Float idz2, Float h2)
{
	int a = 0;

	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), &cldensity);
	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), bDir ? &clphi : &clres);
	clSetKernelArg(kernelGet_res, a++, sizeof(cl_mem), bDir ? &clres : &clphi);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &G4pi);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idx2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idy2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &idz2);
	clSetKernelArg(kernelGet_res, a++, sizeof(Float), &h2);

	clEnqueueWriteBuffer(command_queue, cldensity, CL_TRUE, 0, sizeof(dens), dens,
						 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, bDir ? clphi : clres, CL_TRUE, 0,
			bDir ? sizeof(phi) : sizeof(res), bDir ? phi : res, 0, NULL, NULL);
	clEnqueueNDRangeKernel(command_queue, kernelGet_res, 1, NULL, &globalThreads,
						   &localThreads, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, bDir ? clres : clphi, CL_TRUE, 0,
			bDir ? sizeof(res) : sizeof(phi), bDir ? res : phi, 0, NULL, NULL);
}
#endif

void get_phi()	// Grav. potential via Poisson's Eq.
{
	Float rtot1, rtot2;
	const Float h2 = (Float)0.5 / (idx2 + idy2 + idz2);
	int i, j, k;
	const Float G4pi = 4 * pi * G;
	int nrmp = 0;

	rtot1 = -1;
	rtot2 = 1;
	for (nrmp = 0; fabs(rtot2 - rtot1) > tolPSN * fabs(rtot1); nrmp++)
	{
		if (nrmp > NRMP)
		{
			fprintf(stderr, "Poisson not converging after %d iterations "
							"(last: %lg -> %lg)\n", nrmp, rtot2, rtot1);
			fflush(stderr);
			return;
		}
		rtot2 = rtot1;
		rtot1 = 0;
#ifndef USECL
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
		{
			res(i, j, k) = ((phi(i - 1, j, k) + phi(i + 1, j, k)) * idx2 +
							(phi(i, j - 1, k) + phi(i, j + 1, k)) * idy2 +
							(phi(i, j, k - 1) + phi(i, j, k + 1)) * idz2 -
							G4pi * density(i, j, k)) * h2;
		}
#else
		get_phi_kernel(true, G4pi, idx2, idy2, idz2, h2);
		get_phi_kernel(false, G4pi, idx2, idy2, idz2, h2);
#endif
		for (k = 1; k < Nz; k++)
			for (j = 1; j < Ny; j++)
				for (i = 1; i < Nx; i++)
		{
#ifndef USECL
			phi(i, j, k) = ((res(i - 1, j, k) + res(i + 1, j, k)) * idx2 +
							(res(i, j - 1, k) + res(i, j + 1, k)) * idy2 +
							(res(i, j, k - 1) + res(i, j, k + 1)) * idz2 -
							G4pi * density(i, j, k)) * h2;
#endif
			rtot1 += phi(i, j, k);
		}
	}
#ifdef SHOW_LOOPS
	fprintf(stderr, "Poisson required %d interations.\n", nrmp);
	fflush(stderr);
#endif
}

//**********************************************************************
void get_density()
{
	int i, j, k;

	for (i = 0; i <= Nx; i++)
		for (j = 0; j <= Ny; j++)
			for (k = 0; k <= Nz; k++)
		density(i, j, k) = SQ(real(psi(i, j, k))) + SQ(imag(psi(i, j, k)));
}
