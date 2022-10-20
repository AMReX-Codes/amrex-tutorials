#ifndef AMR_WAVE_H_
#define AMR_WAVE_H_

#include <AMReX_Amr.H>

class AmrWave
    : public amrex::Amr
{
public:
    using amrex::Amr::Amr;

    virtual ~AmrWave ();

    // If we need to override any virtual functions in amrex::Amr, we can do
    // it here.
};

#endif
