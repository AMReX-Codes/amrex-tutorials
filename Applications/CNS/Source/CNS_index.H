#ifndef CNS_INDEX_H_
#define CNS_INDEX_H_

// Components of the conserved state.  All three momentum components are kept
// in 2D as well; the out-of-plane one is simply advected.
inline constexpr int URHO  = 0;
inline constexpr int UMX   = 1;
inline constexpr int UMY   = 2;
inline constexpr int UMZ   = 3;
inline constexpr int UEDEN = 4;  // rho*E (total energy)
inline constexpr int UEINT = 5;  // rho*e (internal energy)
inline constexpr int UTEMP = 6;  // temperature
inline constexpr int NCONS = 7;

// The first NEQNS components are evolved by fluxes.  The others are
// recomputed from them after each stage.
inline constexpr int NEQNS = 5;

// Primitive variables used by the hydro kernels
inline constexpr int QRHO  = 0;
inline constexpr int QU    = 1;
inline constexpr int QV    = 2;
inline constexpr int QW    = 3;
inline constexpr int QPRES = 4;
inline constexpr int QCS   = 5;  // sound speed
inline constexpr int NPRIM = 6;

#endif
