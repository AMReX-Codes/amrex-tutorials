#include "AmrWave.H"
#include "AmrLevelWave.H"

using namespace amrex;

AmrWave::~AmrWave ()
{
    MultiFab const& S = this->getLevel(0).get_new_data(State_Type);
    amrex::Print() << "At the end of simulation, the min and max of the wave are "
                   << S.min(0) << " and " << S.max(0) << "\n\n";
}
