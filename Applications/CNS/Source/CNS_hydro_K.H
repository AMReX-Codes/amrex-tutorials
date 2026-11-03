#ifndef CNS_HYDRO_K_H_
#define CNS_HYDRO_K_H_

// Second-order Godunov kernels: conserved -> primitive conversion,
// characteristic slope limiting, and an approximate Riemann solver.  The
// direction-dependent kernels take the direction as an argument.

#include "CNS_index.H"
#include "CNS_parm.H"

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include <cmath>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_ctoprim (int i, int j, int k,
             amrex::Array4<amrex::Real const> const& u,
             amrex::Array4<amrex::Real> const& q,
             Parm const& parm) noexcept
{
    using amrex::Real;

    Real rho = amrex::max(u(i,j,k,URHO),parm.smallr);
    Real rhoinv = Real(1.0)/rho;
    Real ux = u(i,j,k,UMX)*rhoinv;
    Real uy = u(i,j,k,UMY)*rhoinv;
    Real uz = u(i,j,k,UMZ)*rhoinv;
    Real kineng = Real(0.5)*rho*(ux*ux+uy*uy+uz*uz);
    Real ei = u(i,j,k,UEDEN) - kineng;
    if (ei <= Real(0.0)) { ei = u(i,j,k,UEINT); }
    Real p = amrex::max((parm.eos_gamma-Real(1.0))*ei,parm.smallp);

    q(i,j,k,QRHO) = rho;
    q(i,j,k,QU) = ux;
    q(i,j,k,QV) = uy;
    q(i,j,k,QW) = uz;
    q(i,j,k,QPRES) = p;
    q(i,j,k,QCS) = std::sqrt(parm.eos_gamma*p*rhoinv);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
cns_limiter (amrex::Real dlft, amrex::Real drgt) noexcept
{
    using amrex::Real;

    Real dcen = Real(0.5)*(dlft+drgt);
    Real dsgn = std::copysign(Real(1.0), dcen);
    Real slop = Real(2.0) * amrex::min(std::abs(dlft),std::abs(drgt));
    Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * amrex::min(dlim,std::abs(dcen));
}

// Velocity components normal (qn) and tangential (qt1, qt2) to direction dir
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_vel_index (int dir, int& qn, int& qt1, int& qt2) noexcept
{
    qn  = QU + dir;
    qt1 = (dir == 0) ? QV : QU;
    qt2 = (dir == 2) ? QV : QW;
}

// Limited slopes of the characteristic variables in direction dir
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope (int dir, int i, int j, int k,
           amrex::Array4<amrex::Real> const& dq,
           amrex::Array4<amrex::Real const> const& q) noexcept
{
    using amrex::Real;

    int qn, qt1, qt2;
    cns_vel_index(dir, qn, qt1, qt2);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);

    auto lft = [&] (int n) { return q(i,j,k,n) - q(i-di,j-dj,k-dk,n); };
    auto rgt = [&] (int n) { return q(i+di,j+dj,k+dk,n) - q(i,j,k,n); };

    const Real rho = q(i,j,k,QRHO);
    const Real cs = q(i,j,k,QCS);
    const Real cs2 = cs*cs;

    Real d0 = cns_limiter(Real(0.5)*lft(QPRES)/cs - Real(0.5)*rho*lft(qn),
                          Real(0.5)*rgt(QPRES)/cs - Real(0.5)*rho*rgt(qn));
    Real d1 = cns_limiter(lft(QRHO) - lft(QPRES)/cs2,
                          rgt(QRHO) - rgt(QPRES)/cs2);
    Real d2 = cns_limiter(Real(0.5)*lft(QPRES)/cs + Real(0.5)*rho*lft(qn),
                          Real(0.5)*rgt(QPRES)/cs + Real(0.5)*rho*rgt(qn));
    Real d3 = cns_limiter(lft(qt1), rgt(qt1));
    Real d4 = cns_limiter(lft(qt2), rgt(qt2));

    dq(i,j,k,0) = d0;
    dq(i,j,k,1) = d1;
    dq(i,j,k,2) = d2;
    dq(i,j,k,3) = d3;
    dq(i,j,k,4) = d4;
}

// Approximate (two-shock) Riemann solver returning the Godunov fluxes
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann (const amrex::Real gamma, const amrex::Real smallp,
             const amrex::Real rl, const amrex::Real ul, const amrex::Real pl,
             const amrex::Real ut1l, const amrex::Real ut2l,
             const amrex::Real rr, const amrex::Real ur, const amrex::Real pr,
             const amrex::Real ut1r, const amrex::Real ut2r,
             amrex::Real& flxrho, amrex::Real& flxu, amrex::Real& flxut,
             amrex::Real& flxutt, amrex::Real& flxe) noexcept
{
    using amrex::Real;

    constexpr Real weakwv = Real(1.e-3);
    constexpr Real small = Real(1.e-6);

    Real clsql = gamma*pl*rl;
    Real clsqr = gamma*pr*rr;
    Real wl = std::sqrt(clsql);
    Real wr = std::sqrt(clsqr);
    Real cleft = wl/rl;
    Real cright = wr/rr;
    Real ccsmall = small*(cleft+cright);

    Real pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);
    Real pstnm1 = pstar;

    Real wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
    Real wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    Real ustarp = ul - (pstar-pl)/wl;
    Real ustarm = ur + (pstar-pr)/wr;

    pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = amrex::max(pstar,smallp);

    Real ustar = Real(0.);
    for (int iter = 0; iter < 3; ++iter)
    {
        wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
        wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

        wl = Real(1.)/std::sqrt(wlsq);
        wr = Real(1.)/std::sqrt(wrsq);

        Real ustnm1 = ustarm;
        Real ustnp1 = ustarp;

        ustarm = ur - (pr - pstar)*wr;
        ustarp = ul + (pl - pstar)*wl;

        Real dpditer = std::abs(pstnm1-pstar);
        Real zp = std::abs(ustarp-ustnp1);
        if (zp-weakwv*cleft < Real(0.0) ) {
            zp = dpditer*wl;
        }
        Real zm = std::abs(ustarm-ustnm1);
        if (zm-weakwv*cright < Real(0.0) ) {
            zm = dpditer*wr;
        }

        Real zz = zp+zm;
        Real denom = dpditer/ amrex::max(zz,ccsmall);
        pstnm1 = pstar;
        pstar = pstar - denom*(ustarm-ustarp);
        pstar = amrex::max(pstar,smallp);
        ustar = Real(0.5)*(ustarm+ustarp);
    }

    Real ro, uo, po, sgnm, utrans1, utrans2;
    if (ustar > Real(0.)) {
        ro = rl;
        uo = ul;
        po = pl;
        sgnm = Real(1.);
        utrans1 = ut1l;
        utrans2 = ut2l;
    } else if (ustar < Real(0.)) {
        ro = rr;
        uo = ur;
        po = pr;
        sgnm = Real(-1.);
        utrans1 = ut1r;
        utrans2 = ut2r;
    } else {
        uo = Real(0.5)*(ur+ul);
        po = Real(0.5)*(pr+pl);
        ro = Real(2.)*(rl*rr)/(rl+rr);
        sgnm = Real(1.);
        utrans1 = Real(0.5)*(ut1l+ut1r);
        utrans2 = Real(0.5)*(ut2l+ut2r);
    }
    Real wosq = (Real(0.5)*(gamma-Real(1.))*(pstar+po)+pstar)*ro;
    Real co = std::sqrt(gamma * po / ro);
    Real wo = std::sqrt(wosq);
    Real dpjmp = pstar-po;
    Real rstar = ro/(Real(1.)-ro*dpjmp/wosq);
    Real cstar = std::sqrt(gamma * pstar / rstar);
    Real spout = co-sgnm*uo;
    Real spin = cstar - sgnm*uo;
    if(pstar >= po) {
        spin = wo/ro-sgnm*uo;
        spout = spin;
    }
    Real ss = amrex::max(spout-spin, spout+spin);
    Real frac = Real(0.5)*(Real(1.)+(spin+spout)/amrex::max(ss,ccsmall));

    Real rgdnv, ugdnv, pgdnv;
    if (spout < Real(0.)) {
        rgdnv = ro;
        ugdnv = uo;
        pgdnv = po;
    } else if(spin >= Real(0.)) {
        rgdnv = rstar;
        ugdnv = ustar;
        pgdnv = pstar;
    } else {
        rgdnv = frac*rstar + (Real(1.) - frac)* ro;
        ugdnv = frac*ustar + (Real(1.) - frac)* uo;
        pgdnv = frac*pstar + (Real(1.) - frac)* po;
    }

    flxrho = rgdnv*ugdnv;
    flxu = rgdnv*ugdnv*ugdnv+pgdnv;
    flxut = rgdnv*ugdnv*utrans1;
    flxutt = rgdnv*ugdnv*utrans2;
    flxe = ugdnv*(Real(0.5)*rgdnv*(ugdnv*ugdnv+utrans1*utrans1+utrans2*utrans2) + pgdnv/(gamma -Real(1.)) + pgdnv);
}

// Flux through the face between cells (i,j,k)-e_dir and (i,j,k).  Only the
// first NEQNS components have fluxes.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_flux (int dir, int i, int j, int k,
          amrex::Array4<amrex::Real> const& f,
          amrex::Array4<amrex::Real const> const& dq,
          amrex::Array4<amrex::Real const> const& q,
          Parm const& parm) noexcept
{
    using amrex::Real;

    int qn, qt1, qt2;
    cns_vel_index(dir, qn, qt1, qt2);
    const int il = i - (dir == 0), jl = j - (dir == 1), kl = k - (dir == 2);

    // Left state, extrapolated from cell (il,jl,kl)
    Real cspeed = q(il,jl,kl,QCS);
    Real rl = q(il,jl,kl,QRHO) + Real(0.5) * ( (dq(il,jl,kl,0)+dq(il,jl,kl,2))/cspeed + dq(il,jl,kl,1));
    rl = amrex::max(rl, parm.smallr);
    Real ul = q(il,jl,kl,qn) + Real(0.5) * ( (dq(il,jl,kl,2)-dq(il,jl,kl,0))/q(il,jl,kl,QRHO));
    Real pl = q(il,jl,kl,QPRES) + Real(0.5) *  (dq(il,jl,kl,0)+dq(il,jl,kl,2))*cspeed;
    pl = amrex::max(pl, parm.smallp);
    Real ut1l = q(il,jl,kl,qt1) + Real(0.5) * dq(il,jl,kl,3);
    Real ut2l = q(il,jl,kl,qt2) + Real(0.5) * dq(il,jl,kl,4);

    // Right state, extrapolated from cell (i,j,k)
    cspeed = q(i,j,k,QCS);
    Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
    rr = amrex::max(rr, parm.smallr);
    Real ur = q(i,j,k,qn) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
    Real pr = q(i,j,k,QPRES) - Real(0.5) * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
    pr = amrex::max(pr, parm.smallp);
    Real ut1r = q(i,j,k,qt1) - Real(0.5) * dq(i,j,k,3);
    Real ut2r = q(i,j,k,qt2) - Real(0.5) * dq(i,j,k,4);

    // Momentum components normal and tangential to the face
    const int umn = UMX + (qn - QU);
    const int umt1 = UMX + (qt1 - QU);
    const int umt2 = UMX + (qt2 - QU);

    cns_riemann(parm.eos_gamma, parm.smallp,
                rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r,
                f(i,j,k,URHO), f(i,j,k,umn), f(i,j,k,umt1), f(i,j,k,umt2), f(i,j,k,UEDEN));

    for (int n = NEQNS; n < NCONS; ++n) { f(i,j,k,n) = Real(0.0); }
}

#endif
