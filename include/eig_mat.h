/*
 * eig_mat.h
 *
 *  Created on: Sep 10, 2011
 *      Author: sebastien
 */

#ifndef EIG_MAT_H_
#define EIG_MAT_H_

#include <Python.h>
#include <math.h>
#include "typedefs.h"

/*void eigen_decomposition(matrix a, matrix eig_vec, rvec eig_val);*/

void complete_basis(rvec u, rvec v, rvec w);
void eigen_33_sym(matrix a, matrix eig_vec, rvec eig_val);

#endif /* EIG_MAT_H_ */
