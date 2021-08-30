#ifndef _INDEXING_H_
#define _INDEXING_H_

// Indexes into the grid data: velocity, processor id, and density phi
namespace Idx {
    enum GridIndexes {AMREX_D_DECL(xvel=0, yvel, zvel), proc, phi};
}

#endif
