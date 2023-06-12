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

#pragma OPENCL EXTENSION cl_amd_printf:enable
//#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics:enable

#define KERNEL
#include "parameters3.h"

#define NR ((Nx>Ny?Nx:Ny)>Nz?(Nx>Ny?Nx:Ny):Nz)
#define NN ((Nx + 1) * (Ny + 1) * (Nz + 1))
#define DX ((Ny + 1) * (Nz + 1))
#define DY (Nz + 1)
#define DZ (1)
#define NX ((Ny - 1) * (Nz - 1))
#define NY ((Nx - 1) * (Nz - 1))
#define NZ ((Nx - 1) * (Ny - 1))

void thomas(float2 *l,	// Tridiagonal matrix solver with Thomas
			float2 *d,	// algorithm. The arrays l,d,u and r contain
			float2 *u,	// respectively the subdiagonal, the diagonal,
			float2 *r,	// the superdiagonal and the right hand side.
			int m)		// m is the size of the system. At output the
{						// solution is put in r, while the original
						// right hand side is destroyed.
	int j;
	float2 a;
    
	for (j = 2; j <= m; j++)
	{
		// a = -l[j] / d[j - 1];
		a.x = -(l[j].x * d[j - 1].x + l[j].y * d[j - 1].y) / (SQ(d[j - 1].x) + SQ(d[j - 1].y));
		a.y = -(l[j].y * d[j - 1].x - l[j].x * d[j - 1].y) / (SQ(d[j - 1].x) + SQ(d[j - 1].y));
		d[j].x += a.x * u[j - 1].x - a.y * u[j - 1].y;
		d[j].y += a.x * u[j - 1].y + a.y * u[j - 1].x;
		r[j].x += a.x * r[j - 1].x - a.y * r[j - 1].y;
		r[j].y += a.x * r[j - 1].y + a.y * r[j - 1].x;
	}
	//	r[m] = r[m] / d[m];
	r[m].x = (r[m].x * d[m].x + r[m].y * d[m].y) / (SQ(d[m].x) + SQ(d[m].y));
	r[m].y = (r[m].y * d[m].x - r[m].x * d[m].y) / (SQ(d[m].x) + SQ(d[m].y));
	for (j = m - 1; j >= 1; j--)
	{
		float2 rr;
		//	r[j] = (r[j] - u[j] * r[j + 1]) / d[j];
		rr.x = r[j].x - u[j].x * r[j + 1].x + u[j].y * r[j + 1].y;
		rr.y = r[j].y - u[j].x * r[j + 1].y - u[j].y * r[j + 1].x;
		r[j].x = (rr.x * d[j].x + rr.y * d[j].y) / (SQ(d[j].x) + SQ(d[j].y));
		r[j].y = (rr.y * d[j].x - rr.x * d[j].y) / (SQ(d[j].x) + SQ(d[j].y));
	}
}

__kernel void solve_x(__global float2 *psi, __global float2 *psi_n,
					  __global float *U, float2 foo1, float2 foo2, float2 foo3,
					  float2 foo4, float dy)
{
	float2 l[NR], d[NR], u[NR], r[NR];
	float y;
	int i, j, k, n;
	uint tsz = get_global_size(0);
	uint tid = get_global_id(0);

	for (; tid < NX; tid += tsz)
	{
		k = tid / (Ny - 1) + 1;
		j = tid % (Ny - 1) + 1;

		y = yl + j * dy;
		l[1] = foo1 + foo2 * y;
		d[1] = foo3;
		u[1] = foo1 - foo2 * y;
		for (i = 2; i < Nx; i++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			l[i] = l[1];
			d[i] = d[1] + foo4 * U[n];
			u[i] = u[1];
			r[i] = psi_n[n];
		}
		n = (Ny + 1 + j) * (Nz + 1) + k;
		d[1] += foo4 * U[n];
		r[1] = psi_n[n];
		thomas(l, d, u, r, Nx - 1);
		for (i = 1; i < Nx; i++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			psi[n] = r[i];
		}
		psi[n = j * (Nz + 1) + k] = 0;
		psi[(Nx * (Ny + 1) + j) * (Nz + 1) + k] = 0;
	}
}

__kernel void solve_y(__global float2 *psi, __global float2 *psi_n,
					  __global float *U, float2 foo1, float2 foo2, float2 foo3,
					  float2 foo4, float dx)
{
	float2 l[NR], d[NR], u[NR], r[NR];
	float x;
	int i, j, k, n;
	uint tsz = get_global_size(0);
	uint tid = get_global_id(0);

	for (; tid < NY; tid += tsz)
	{
		i = tid / (Nz - 1) + 1;
		k = tid % (Nz - 1) + 1;
		x = xl + i * dx;
		l[1] = foo1 - foo2 * x;
		d[1] = foo3;
		u[1] = foo1 + foo2 * x;
		for (j = 2; j < Ny; j++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			l[j] = l[1];
			d[j] = d[1] + foo4 * U[n];
			u[j] = u[1];
			r[j] = psi_n[n];
		}
		n = (i * (Ny + 1) + 1) * (Nz + 1) + k;
		d[1] += foo4 * U[n];
		r[1] = psi_n[n];
		thomas(l, d, u, r, Ny - 1);
		for (j = 1; j < Ny; j++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			psi[n] = r[j];
		}
		psi[i * (Ny + 1) * (Nz + 1) + k] = 0;
		psi[(i * (Ny + 1) + Ny) * (Nz + 1) + k] = 0;
	}
}

__kernel void solve_z(__global float2 *psi, __global float2 *psi_n,
					  __global float *U, float2 foo1, float2 foo3, float2 foo4)
{
	float2 l[NR], d[NR], u[NR], r[NR];
	int i, j, k, n;
	uint tsz = get_global_size(0);
	uint tid = get_global_id(0);

	for (; tid < NZ; tid += tsz)
	{
		i = tid / (Ny - 1) + 1;
		j = tid % (Ny - 1) + 1;

		l[1] = u[1] = foo1;
		d[1] = foo3;
		for (k = 2; k < Nz; k++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			l[k] = u[k] = l[1];
			d[k] = d[1] + foo4 * U[n];
			r[k] = psi_n[n];
		}
		n = (i * (Ny + 1) + j) * (Nz + 1) + 1;
		r[1] = psi_n[n];
		d[1] += foo4 * U[n];
		thomas(l, d, u, r, Nz - 1);
		for (k = 1; k < Nz; k++)
		{
			n = (i * (Ny + 1) + j) * (Nz + 1) + k;
			psi[n] = r[k];
		}
		psi[(i * (Ny + 1) + j) * (Nz + 1)] = 0;
		psi[(i * (Ny + 1) + j) * (Nz + 1) + Nz] = 0;
	}
}

__kernel void calc_rhs_x(__global float2 *psi, __global float2 *psi_n,
						 __global float *U, float dy, float2 foo1,
						 float2 foo3, float2 foo4, float2 foo5)
{
	int i, j, k;
	float y;
	float2 foo2, foo6, foo7, foo8;
	uint tid = get_global_id(0);
	uint tsz = get_global_size(0);

	for (; tid < NN; tid += tsz)
	{
		k = tid % (Nz + 1);
		if (k > 0 && k < Nz)
		{
			j = tid / (Nz + 1);
			i = j / (Ny + 1);
			if (i > 0 && i < Nx)
			{
				j %= (Ny + 1);
				if (j > 0 && j < Ny)
				{
					y = yl + j * dy;
					foo2 = foo1 * y;
					foo6 = foo3 - foo2;
					foo7.x = foo5.x - foo4.x * U[tid];
					foo7.y = foo5.y - foo4.y * U[tid];
					foo8 = foo3 + foo2;
					psi_n[tid].x = foo6.x * psi[tid - DX].x - foo6.y * psi[tid - DX].y +
							foo7.x * psi[tid].x - foo7.y * psi[tid].y +
							foo8.x * psi[tid + DX].x - foo8.y * psi[tid + DX].y;
					psi_n[tid].y = foo6.y * psi[tid - DX].x + foo6.x * psi[tid - DX].y +
							foo7.y * psi[tid].x + foo7.x * psi[tid].y +
							foo8.y * psi[tid + DX].x + foo8.x * psi[tid + DX].y;
				}
				else psi_n[tid] = 0;
			}
			else psi_n[tid] = 0;
		}
		else psi_n[tid] = 0;
	}
}

__kernel void calc_rhs_y(__global float2 *psi, __global float2 *psi_n,
						 __global float *U, float dx, float2 foo1,
						 float2 foo3, float2 foo4, float2 foo5)
{
	int i, j, k;
	float x;
	float2 foo2, foo6, foo7, foo8;
	uint tid = get_global_id(0);
	uint tsz = get_global_size(0);

	for (; tid < NN; tid += tsz)
	{
		k = tid % (Nz + 1);
		if (k > 0 && k < Nz)
		{
			j = tid / (Nz + 1);
			i = j / (Ny + 1);
			if (i > 0 && i < Nx)
			{
				j %= (Ny + 1);
				if (j > 0 && j < Ny)
				{
					x = xl + i * dx;
					foo2 = foo1 * x;
					foo6 = foo3 + foo2;
					foo7.x = foo5.x - foo4.x * U[tid];
					foo7.y = foo5.y - foo4.y * U[tid];
					foo8 = foo3 - foo2;
					psi_n[tid].x = foo6.x * psi[tid - DY].x - foo6.y * psi[tid - DY].y +
							foo7.x * psi[tid].x - foo7.y * psi[tid].y +
							foo8.x * psi[tid + DY].x - foo8.y * psi[tid + DY].y;
					psi_n[tid].y = foo6.y * psi[tid - DY].x + foo6.x * psi[tid - DY].y +
							foo7.y * psi[tid].x + foo7.x * psi[tid].y +
							foo8.y * psi[tid + DY].x + foo8.x * psi[tid + DY].y;
				}
				else psi_n[tid] = 0;
			}
			else psi_n[tid] = 0;
		}
		else psi_n[tid] = 0;
	}
}

__kernel void calc_rhs_z(__global float2 *psi, __global float2 *psi_n,
						 __global float *U, float2 foo3, float2 foo4,
						 float2 foo5)
{
	int i, j, k;
	float2 foo7;
	uint tid = get_global_id(0);
	uint tsz = get_global_size(0);

	for (; tid < NN; tid += tsz)
	{
		k = tid % (Nz + 1);
		if (k > 0 && k < Nz)
		{
			j = tid / (Nz + 1);
			i = j / (Ny + 1);
			if (i > 0 && i < Nx)
			{
				j %= (Ny + 1);
				if (j > 0 && j < Ny)
				{
					foo7.x = foo5.x - foo4.x * U[tid];
					foo7.y = foo5.y - foo4.y * U[tid];
					psi_n[tid].x = foo3.x * psi[tid - DZ].x - foo3.y * psi[tid - DZ].y +
							foo7.x * psi[tid].x - foo7.y * psi[tid].y +
							foo3.x * psi[tid + DZ].x - foo3.y * psi[tid + DZ].y;
					psi_n[tid].y = foo3.y * psi[tid - DZ].x + foo3.x * psi[tid - DZ].y +
							foo7.y * psi[tid].x + foo7.x * psi[tid].y +
							foo3.y * psi[tid + DZ].x + foo3.x * psi[tid + DZ].y;
				}
				else psi_n[tid] = 0;
			}
			else psi_n[tid] = 0;
		}
		else psi_n[tid] = 0;
	}
}

__kernel void get_res(__global float *density, __global float *phi,
					  __global float *res, float G4pi, float idx2,
					  float idy2, float idz2, float h2)
{
	int i, j, k;
	uint tid = get_global_id(0);
	uint tsz = get_global_size(0);
	for (; tid < NN; tid += tsz)
	{
		k = tid % (Nz + 1);
		if (k > 0 && k < Nz)
		{
			j = tid / (Nz + 1);
			i = j / (Ny + 1);
			if (i > 0 && i < Nx)
			{
				j %= (Ny + 1);
				if (j > 0 && j < Ny)
				{
					res[tid] = ((phi[tid - DX] + phi[tid + DX]) * idx2 +
								(phi[tid - DY] + phi[tid + DY]) * idy2 +
								(phi[tid - DZ] + phi[tid + DZ]) * idz2 -
								G4pi * density[tid]) * h2;
				}
				else res[tid] = phi[tid];
			}
			else res[tid] = phi[tid];
		}
		else res[tid] = phi[tid];
	}
}
