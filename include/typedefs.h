/*
 * typedefs.h
 *
 *  Created on: 22 janv. 2010
 *      Author: sebastien
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

/*
 * make sure Python is included first
 */
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>


/*
 * Constants
 */
#define XX	    0			/* Defines for indexing in	*/
#define	YY	    1			/* vectors			*/
#define ZZ   	2
#define DIM   	3			/* Dimension of vectors		*/
#define XXXX    0           /* defines to index matrices */
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b) )
#endif
#ifndef even
#define even(a) ( ( (a+1) / 2) == (a / 2) )
#endif

/*
 * Types used by FATSLiM
 */
typedef npy_int64 fsl_int;

typedef npy_uint64 fsl_uint;

typedef npy_float64 real;

typedef real rvec[DIM];

typedef real matrix[DIM][DIM];

/*
 * Useful functions
 */
real real_min(real a, real b);
real real_max(real a, real b);
real real_abs(real a);
real rvec_norm2(const rvec a);
real rvec_norm(const rvec a);
void rvec_normalize(rvec a);
void rvec_add(const rvec a,const rvec b,rvec c);
void rvec_inc(rvec a,const rvec b);
void rvec_dec(rvec a,const rvec b);
void rvec_sub(const rvec a,const rvec b,rvec c);
void rvec_smul(real a,const rvec v1,rvec v2);
void rvec_cprod(const rvec a, const rvec b, rvec c);
void rvec_cprod_norm(const rvec a, const rvec b, rvec c);
real rvec_dprod(const rvec a, const rvec b);
void rvec_copy(const rvec src,rvec dest);
void rvec_clear(rvec a);
void mat_clear(matrix a);
void mat_copy(matrix src,matrix dest);
void invert_mat(matrix mat, matrix inverted_mat);
void mat_from_rvec(rvec v1, rvec v2, rvec v3, matrix mat);

#endif /* _TYPEDEFS_H_ */
