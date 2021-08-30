
#include "MyTest.H"
#include "initProb_K.H"

using namespace amrex;

void
MyTest::initProbPoisson ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        const auto prob_lo = geom.ProbLoArray();
        const auto dx      = geom.CellSizeArray();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto rhsfab = rhs.array(mfi);
            auto exactfab = exact_solution.array(mfi);
            AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                actual_init_poisson(i,j,k,rhsfab,exactfab,prob_lo,dx);
            });
        }

        solution.setVal(0.0);
    }
}

