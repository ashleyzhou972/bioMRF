#include "add.h"
#include "dm_call.h"
#include <R_ext/Rdynload.h>


void R_init_bioMRF(DllInfo *info) {
	R_RegisterCCallable("bioMRF", "add",  (DL_FUNC) &add_);
  R_RegisterCCallable("bioMRF", "pmrf",  (DL_FUNC) &double_metropolis_cont);
}
