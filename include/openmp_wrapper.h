#ifndef _OPENMPWRAPPER_
#define _OPENMPWRAPPER_

#ifdef _OPENMP
#include <omp.h>
#else

int omp_get_max_threads(void)
{
    return 1;
}

int omp_get_num_procs(void)
{
    return 1;
}

int omp_get_thread_num(void)
{
    return 0;
}

#endif /* _OPENMP */

#endif /* _OPENMPWRAPPER_ */