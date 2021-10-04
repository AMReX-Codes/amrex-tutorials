#include "particles/ElectrostaticParticleContainer.H"

void ElectrostaticParticleContainer::writeParticles(int n) {
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
