#include <typedefs.h>

/* Useful functions, also ripped from GROMACS */
real real_min(real a, real b)
{
    return min(a,b);
}

real real_max(real a, real b)
{
    return max(a,b);
}

real real_abs(real a)
{
    return (a < 0) ? -a : a;
}

real rvec_norm2(const rvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

real rvec_norm(const rvec a)
{
  return sqrt(rvec_norm2(a));
}

void rvec_normalize(rvec a)
{
	real vec_norm = rvec_norm(a);

	a[XX] /= vec_norm;
	a[YY] /= vec_norm;
	a[ZZ] /= vec_norm;
}

void rvec_add(const rvec a,const rvec b,rvec c)
{
  real x,y,z;

  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];

  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

void rvec_inc(rvec a,const rvec b)
{
  real x,y,z;

  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];

  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

void rvec_dec(rvec a,const rvec b)
{
  real x,y,z;

  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];

  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

void rvec_sub(const rvec a,const rvec b,rvec c)
{
  real x,y,z;

  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];

  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

void rvec_smul(real a,const rvec v1,rvec v2)
{
	v2[XX]=a*v1[XX];
	v2[YY]=a*v1[YY];
 	v2[ZZ]=a*v1[ZZ];
}

void rvec_cprod(const rvec a, const rvec b, rvec c)
{
	c[XX] = a[YY] * b[ZZ] - a[ZZ] * b[YY];
	c[YY] = a[ZZ] * b[XX] - a[XX] * b[ZZ];
	c[ZZ] = a[XX] * b[YY] - a[YY] * b[XX];
}

void rvec_cprod_norm(const rvec a, const rvec b, rvec c)
{
	rvec_cprod(a, b, c);
	rvec_normalize(c);
}

real rvec_dprod(const rvec a, const rvec b)
{
	return a[XX] * b[XX] + a[YY] * b[YY] + a[ZZ] * b[ZZ];
}

void rvec_copy(const rvec src,rvec dest)
{
  dest[XX]=src[XX];
  dest[YY]=src[YY];
  dest[ZZ]=src[ZZ];
}

void rvec_clear(rvec a)
{
  a[XX]=0.0;
  a[YY]=0.0;
  a[ZZ]=0.0;
}

void mat_clear(matrix a)
{
	rvec_clear(a[XX]);
	rvec_clear(a[YY]);
	rvec_clear(a[ZZ]);
}

void mat_copy(matrix src,matrix dest)
{
  rvec_copy(src[XX],dest[XX]);
  rvec_copy(src[YY],dest[YY]);
  rvec_copy(src[ZZ],dest[ZZ]);
}



/* Caution: matrix's determinant is not tested to be != 0 */
void invert_mat(matrix mat, matrix inverted_mat)
{
	real det;
	real A, B, C, D, E, F, G, H, K;
	real a, b, c, d, e, f, g, h, k;

	a = mat[XX][XX];
	b = mat[XX][YY];
	c = mat[XX][ZZ];
	d = mat[YY][XX];
	e = mat[YY][YY];
	f = mat[YY][ZZ];
	g = mat[ZZ][XX];
	h = mat[ZZ][YY];
	k = mat[ZZ][ZZ];

	A = e * k - f * h; /* ek-fh */
	B = f * g - d * k; /* fg-dk */
	C = d * h - e * g; /* dh-eg */
	D = c * h - b * k; /* ch-bk */
	E = a * k - c * g; /* ak-cg */
	F = g * b - a * h; /* gb-ah */
	G = b * f - c * e; /* bf-ce */
	H = c * d - a * f; /* cd-af */
	K = a * e - b * d; /* ae-bd */

	det = a * A + b * B + c * C;
	det = 1 / det; /* Beware: determinant MUST not be 0! */

	inverted_mat[XX][XX] = det * A;
	inverted_mat[XX][YY] = det * D;
	inverted_mat[XX][ZZ] = det * G;
	inverted_mat[YY][XX] = det * B;
	inverted_mat[YY][YY] = det * E;
	inverted_mat[YY][ZZ] = det * H;
	inverted_mat[ZZ][XX] = det * C;
	inverted_mat[ZZ][YY] = det * F;
	inverted_mat[ZZ][ZZ] = det * K;
}

void mat_from_rvec(rvec v1, rvec v2, rvec v3, matrix mat)
{
	mat[XX][XX] = v1[XX];
	mat[XX][YY] = v2[XX];
	mat[XX][ZZ] = v3[XX];
	mat[YY][XX] = v1[YY];
	mat[YY][YY] = v2[YY];
	mat[YY][ZZ] = v3[YY];
	mat[ZZ][XX] = v1[ZZ];
	mat[ZZ][YY] = v2[ZZ];
	mat[ZZ][ZZ] = v3[ZZ];
}
