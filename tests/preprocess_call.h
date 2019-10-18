#include <Rinternals.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparse.h"
#include "eigen.h"

SEXP preprocess(SEXP dim_in, SEXP val_in, SEXP row_ind_in, SEXP col_ptr_in);


