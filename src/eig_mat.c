/*
 * eig_mat.c
 *
 *  Created on: Sep 10, 2011
 *      Author: sebastien
 */



/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   domain Java Matrix library JAMA. */

#include <eig_mat.h>

static void compute_roots(matrix mat, rvec roots)
{
	real c0, c1, c2;
	real one_over_3 = (real) 1 / (real) 3;
	real sqrt_3 = sqrt(3);
	real c2_over_3, a_over_3;
	real half_mb, q, magnitude, angle, cs, sn, root0, root1, root2;
	real mat_XX = mat[XX][XX];
	real mat_XY = mat[XX][YY];
	real mat_XZ = mat[XX][ZZ];
	real mat_YY = mat[YY][YY];
	real mat_YZ = mat[YY][ZZ];
	real mat_ZZ = mat[ZZ][ZZ];

	/*
	 * The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
	 * eigenvalues are the roots to this equation, all guaranteed to be
	 * real-valued, because the matrix is symmetric.
	 */
	c0 = mat_XX * mat_YY * mat_ZZ + 2.0 * mat_XY * mat_XZ * mat_YZ - mat_XX * mat_YZ * mat_YZ -
			mat_YY * mat_XZ * mat_XZ - mat_ZZ * mat_XY * mat_XY;

	c1 = mat_XX * mat_YY - mat_XY * mat_XY + mat_XX * mat_ZZ - mat_XZ * mat_XZ +
			mat_YY * mat_ZZ - mat_YZ * mat_YZ;

	c2 = mat_XX + mat_YY + mat_ZZ;

	/*
	 * Construct the parameters used in classifying the roots of the equation
	 * and in solving the equation for the roots in closed form.
	 */
	c2_over_3 = c2 * one_over_3;
	a_over_3 = (c1 - c2 * c2_over_3) * one_over_3;

	if (a_over_3 > 0.0)
	{
		a_over_3 = 0.0;
	}

	half_mb = 0.5 * (c0 + c2_over_3 * (2.0 * c2_over_3 * c2_over_3 - c1));

	q = half_mb * half_mb + a_over_3 * a_over_3 * a_over_3;

	if (q > 0.0)
	{
		q = 0.0;
	}

	/*
	 * Compute the eigenvalues by solving for the roots of the polynomial.
	 */

	magnitude = sqrt( - a_over_3);
	angle = atan2(sqrt(-q), half_mb) * one_over_3;
	cs = cos(angle);
	sn = sin(angle);

	root0 = c2_over_3 + 2.0 * magnitude * cs;
	root1 = c2_over_3 - magnitude * (cs + sqrt_3 * sn);
	root2 = c2_over_3 - magnitude * (cs - sqrt_3 * sn);


	/*
	 * Sort roots so root0 < root1 < root2
	 */
	if (root1 >= root0)
	{
		roots[0] = root0;
		roots[1] = root1;
	}
	else
	{
		roots[0] = root1;
		roots[1] = root0;
	}

	if (root2 >= roots[1])
	{
		roots[2] = root2;
	}
	else
	{
		roots[2] = roots[1];
		if (root2 >= roots[0])
		{
			roots[1] = root2;
		}
		else
		{
			roots[1] = roots[0];
			roots[0] = root2;
		}
	}

}

static int compute_rank(matrix m)
{
	// Compute the maximum magnitude matrix entry.
	real abs, save, max = -1, inv_max;
	int row, col, maxRow = -1, maxCol = -1;

	for (row = 0; row < DIM; row++)
	{
		for (col = row; col < DIM; col++)
		{
			abs = fabs(m[row][col]);
			if (abs > max)
			{
				max = abs;
				maxRow = row;
				maxCol = col;
			}
		}
	}

	if (max == 0)
	{
		/*
		 * The rank is 0. The eigenvalue has multiplicity 3.
		 */
		return 0;
	}

	/*
	 * The rank is at least 1. Swap the row containing the maximum-magnitude
	 * entry with row 0.
	 */

	if (maxRow != 0)
	{
		for (col = 0; col < 3; col++)
		{
			save = m[0][col];
			m[XX][col] = m[maxRow][col];
			m[maxRow][col] = save;
		}
	}

	/*
	 * Row-reduce the matrix...
	 * Scale row 0 to generate a 1-valued pivot.
	 */
	inv_max = 1/m[XX][maxCol];
	m[XX][XX] *= inv_max;
	m[XX][YY] *= inv_max;
	m[XX][ZZ] *= inv_max;

	/*
	 * Eliminate the maxCol column entries in rows 1 and 2.
	 */
	if (maxCol == XX)
	{
		m[YY][YY] -= m[YY][XX] * m[XX][YY];
		m[YY][ZZ] -= m[YY][XX] * m[XX][ZZ];
		m[ZZ][YY] -= m[ZZ][XX] * m[XX][YY];
		m[ZZ][ZZ] -= m[ZZ][XX] * m[XX][ZZ];
		m[YY][XX] = 0;
		m[ZZ][XX] = 0;
	}
	else if (maxCol == YY)
	{
		m[YY][XX] -= m[YY][YY] * m[XX][XX];
		m[YY][ZZ] -= m[YY][YY] * m[XX][ZZ];
		m[ZZ][XX] -= m[ZZ][YY] * m[XX][XX];
		m[ZZ][ZZ] -= m[ZZ][YY] * m[XX][ZZ];
		m[YY][YY] = 0;
		m[ZZ][YY] = 0;
	}
	else
	{
		m[YY][XX] -= m[YY][ZZ] * m[XX][XX];
		m[YY][YY] -= m[YY][ZZ] * m[XX][YY];
		m[ZZ][XX] -= m[ZZ][ZZ] * m[XX][XX];
		m[ZZ][YY] -= m[ZZ][ZZ] * m[XX][YY];
		m[YY][ZZ] = 0;
		m[ZZ][ZZ] = 0;
	}

	/*
	 * Compute the maximum-magnitude entry of the last two rows of the
	 * row-reduced matrix.
	 */
	max = -1;
	maxRow = -1;
	maxCol = -1;
	for (row = 1; row < DIM; row++)
	{
		for (col = 0; col < DIM; col++)
		{
			abs = fabs(m[row][col]);
			if (abs > max)
			{
				max = abs;
				maxRow = row;
				maxCol = col;
			}
		}
	}

	if (max == XX)
	{
		/*
		 * The rank is 1. The eigenvalue has multiplicity 2.
		 */
		return 1;
	}

	/*
	 * If row 2 has the maximum-magnitude entry, swap it with row 1.
	 */
	if (maxRow == ZZ)
	{
		for (col = 0; col < DIM; col++)
		{
			save = m[YY][col];
			m[YY][col] = m[ZZ][col];
			m[ZZ][col] = save;
		}
	}

	/*
	 * Scale row 1 to generate a 1-valued pivot.
	 */
	inv_max = 1 / m[YY][maxCol];
	m[YY][XX] *= inv_max;
	m[YY][YY] *= inv_max;
	m[YY][ZZ] *= inv_max;

	/*
	 * The rank is 2. The eigenvalue has multiplicity 1.
	 */
	return 2;

}

void complete_basis(rvec z_axis, rvec x_axis, rvec y_axis)
{
	rvec_normalize(z_axis);

	if ( fabs(z_axis[XX]) >= fabs(z_axis[YY]) )
	{
		x_axis[XX] = -z_axis[ZZ];
		x_axis[YY] = 0;
		x_axis[ZZ] = z_axis[XX];
	}
	else
	{
		x_axis[XX] = 0;
		x_axis[YY] = z_axis[ZZ];
		x_axis[ZZ] = -z_axis[YY];
	}

	rvec_normalize(x_axis);
	rvec_cprod(z_axis, x_axis, y_axis);
}

void eigen_33_sym(matrix a, matrix eig_vec, rvec eig_val)
{
	matrix scaled_mat, m0, m1;
	real max_val;
	rvec roots;
	int i,j, rank0, rank1;

	/* scale the matrix */
	max_val = a[XX][XX];
	for (i=0; i<DIM; i++)
		for (j=i; j<DIM; j++)
			if (fabs(a[i][j]) > max_val)
				max_val = fabs(a[i][j]);

	for (i=0; i<DIM; i++)
	  for (j=0; j<DIM; j++)
	  {
		  scaled_mat[i][j] = a[i][j] / max_val;
	  }
	

	/* Get the eigenvalues */
	compute_roots(scaled_mat, roots);

	eig_val[XX] = roots[XX];
	eig_val[YY] = roots[YY];
	eig_val[ZZ] = roots[ZZ];

	/*
	 * Compute A - eig_val[i] * I
	 */
	mat_copy(scaled_mat, m0);
	m0[XX][XX] -= eig_val[0];
	m0[YY][ZZ] -= eig_val[0];
	m0[ZZ][ZZ] -= eig_val[0];

	rank0 = compute_rank(m0);
	if (rank0 == 0)
	{
		/*
		 * eigenvalues are identical
		 */
		rvec_clear(eig_vec[XX]);
		eig_vec[XX][XX] = 1.0;
		rvec_clear(eig_vec[YY]);
		eig_vec[YY][YY] = 1.0;
		rvec_clear(eig_vec[ZZ]);
		eig_vec[ZZ][ZZ] = 1.0;
		return;
	}
	if (rank0 == 1)
	{
		/*
		 * eig_val[0] = eig_val[1] < eig_val[2]
		 */
		complete_basis(m0[0], eig_vec[XX],eig_vec[YY]);
		rvec_cprod(eig_vec[XX], eig_vec[YY], eig_vec[ZZ]);
		return;
	}

	/*
	 * rank0 == 2
	 */
	rvec_cprod_norm(m0[0], m0[1], eig_vec[XX]);

	/*
	 * Compute A - eig_val[1] * I
	 */
	mat_copy(scaled_mat, m1);
	m1[XX][XX] -= eig_val[1];
	m1[YY][ZZ] -= eig_val[1];
	m1[ZZ][ZZ] -= eig_val[1];

	rank1 = compute_rank(m1); // zero rank detected earlier, rank1 must be positive
	if (rank1 == 1)
	{
		/*
		 * eig_val[0] < eig_val[1] = eig_val[2]
		 */
		complete_basis(eig_vec[XX], eig_vec[YY], eig_vec[ZZ]);
		return;
	}
	/*
	 * rank1 == 2
	 */
	rvec_cprod_norm(m1[0], m1[1], eig_vec[YY]);

	/*
	 * eigenvalues must be distinct at this point, rank2 must be 2
	 */
	rvec_cprod(eig_vec[XX], eig_vec[YY], eig_vec[ZZ]);

}
