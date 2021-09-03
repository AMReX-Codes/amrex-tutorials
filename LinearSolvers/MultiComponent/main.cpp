//
// Tutorial:    MultiComponent Linear Solve
//
// File:        main.cpp
//
// Author:      Brandon Runnels
//              University of Colorado Colorado Springs
//              brunnels@uccs.edu
//              solids.uccs.edu
//
// Date:        September 3, 2019
//
// Description: This tutorial demonstrates how to implement a multi-component
//              nodal linear operator. This tutorial demonstrates the
//              "CFStrategy::ghostnodes" method for computing the reflux at
//              the coarse/fine boundary.
//
// See also:    ./MCNodalLinOp.H, ./MCNodalLinOp.cpp
//              for implementation of the MC linear operator.
//

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Machine.H>
#include <AMReX_MLMG.H>
#include <AMReX_PlotFileUtil.H>
#include "MCNodalLinOp.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    Initialize(argc, argv);
    {

    //
    // Read in mesh structure (# amr levels, # nodes)
    // Geometry will be a nnodes x nnodes (x nnodes) grid
    // refined in the center up to nlevels.
    //
    struct {
        int nlevels = 3;
        int nnodes = 32;
        int max_grid_size = 10000000;
        Vector<int> ref_ratio = {4};
    } mesh;
    {
        ParmParse pp("mesh");
        pp.query("nlevels",mesh.nlevels);
        pp.query("nnodes",mesh.nnodes);
        pp.query("max_grid_size",mesh.max_grid_size);
        pp.queryarr("ref_ratio",mesh.ref_ratio);
        //if (mesh.ref_ratio.size() == 1) mesh.ref_ratio.reset(mesh.nlevels-1, mesh.ref_ratio[0]);
    }

    //
    // Read in linear operator parameters:
    //    ncomp = number of components for a MC solve
    //    coeff = ncomp x ncomp list of coefficients
    //
    struct {
        int ncomp=1;
        Vector<Real> coeff = {1.0};
    } op;
    {
        ParmParse pp("op");
        pp.query("ncomp",op.ncomp);
        pp.queryarr("coeff",op.coeff);
    }

    //
    // Read in MLMG solver parameters
    //
    struct {
        int verbose = -1;
        int bottom_verbose = -1;
        int max_iter = -1;
        int fixed_iter = -1;
        int max_fmg_iter = -1;
        int linop_maxorder = -1;
        int agglomeration = -1;
        int consolidation = -1;
        int max_coarsening_level = -1 ;
        int pre_smooth = -1;
        int post_smooth = -1;
    } mlmg;
    {
        ParmParse pp("mlmg");
        pp.query("verbose",mlmg.verbose);
        pp.query("bottom_verbose",mlmg.bottom_verbose );
        pp.query("max_iter",mlmg.max_iter);
        pp.query("max_fmg_iter",mlmg.max_fmg_iter);
        pp.query("agglomeration",mlmg.agglomeration);
        pp.query("consolidation",mlmg.consolidation);
        pp.query("max_coarsening_level",mlmg.max_coarsening_level);
        pp.query("fixed_iter",mlmg.fixed_iter);
        pp.query("pre_smooth",mlmg.pre_smooth);
        pp.query("post_smooth",mlmg.post_smooth);
    }


    //
    // Initialize geometry and grids
    //
    Vector<Geometry> geom;
      Vector<BoxArray> cgrids, ngrids;
     Vector<DistributionMapping> dmap;
      Vector<MultiFab> solution, solgn, rhs, rhsgn, res, resgn, b, bgn;
     geom.resize(mesh.nlevels);
     cgrids.resize(mesh.nlevels);
     ngrids.resize(mesh.nlevels);
     dmap.resize(mesh.nlevels);
     solution.resize(mesh.nlevels);
     solgn.resize(mesh.nlevels);
     rhs.resize(mesh.nlevels);
     rhsgn.resize(mesh.nlevels);
     res.resize(mesh.nlevels);
     resgn.resize(mesh.nlevels);
     b.resize(mesh.nlevels);
     bgn.resize(mesh.nlevels);
    RealBox rb({AMREX_D_DECL(-0.5,-0.5,-0.5)},
              {AMREX_D_DECL(0.5,0.5,0.5)});
    Geometry::Setup(&rb, 0);
    Box NDomain(IntVect{AMREX_D_DECL(0,0,0)},
                IntVect{AMREX_D_DECL(mesh.nnodes,mesh.nnodes,mesh.nnodes)},
                IntVect::TheNodeVector());
    Box CDomain = convert(NDomain, IntVect::TheCellVector());

    //
    // Refine the grid
    //
    Box domain = CDomain;
     for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
         {
             geom[ilev].define(domain);
             if (ilev < mesh.nlevels-1) domain.refine(mesh.ref_ratio[ilev]);
         }
    Box cdomain = CDomain;
    int fac = 1;
     for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
    {
        cgrids[ilev].define(cdomain);
        cgrids[ilev].maxSize(mesh.max_grid_size); // TODO

        
        if (ilev > 0) fac *= (mesh.ref_ratio[ilev-1]/2);
        cdomain.grow(-(mesh.nnodes/4) * fac);
        if (ilev < mesh.nlevels-1) cdomain.refine(mesh.ref_ratio[ilev]);
        ngrids[ilev] = cgrids[ilev];
        ngrids[ilev].convert(IntVect::TheNodeVector());
    }

    //
    // Initialize the solution and rhs fabs.
    // Initialize the RHS fab to the function:
    //    RHS[0] = x1*(1-x1) * x2(1-x2) * x3(1-x3)
    //    RHS[1] = 0
    //    RHS[2] = 0 ... etc
    //
    int nghost = 2;
    for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
    {
        if (ilev > 0) nghost = mesh.ref_ratio[ilev-1];
        dmap   [ilev].define(cgrids[ilev]);
        solution[ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        solution[ilev].setVal(0.0);
        solution[ilev].setMultiGhost(true);
        b[ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        b[ilev].setMultiGhost(true);
        rhs     [ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        rhs     [ilev].setVal(0.0);
        rhs     [ilev].setMultiGhost(true);
        res     [ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        res     [ilev].setVal(0.0);
        res     [ilev].setMultiGhost(true);
        {
            BoxArray grids = ngrids[ilev];
            grids.grow(nghost);
            solgn[ilev].define(grids,dmap[ilev],op.ncomp,0); solgn[ilev].setMultiGhost(true);
            rhsgn[ilev].define(grids,dmap[ilev],op.ncomp,0); rhsgn[ilev].setMultiGhost(true);
            resgn[ilev].define(grids,dmap[ilev],op.ncomp,0); resgn[ilev].setMultiGhost(true);
            bgn[ilev].define(grids,dmap[ilev],op.ncomp,0); bgn[ilev].setMultiGhost(true);
        }

        Box dom(geom[ilev].Domain());
        const Real AMREX_D_DECL( dx = geom[ilev].CellSize()[0],
                                 dy = geom[ilev].CellSize()[1],
                                 dz = geom[ilev].CellSize()[2]);
        const Real AMREX_D_DECL( minx = geom[ilev].ProbLo()[0],
                                 miny = geom[ilev].ProbLo()[1],
                                 minz = geom[ilev].ProbLo()[2]);
        dom.convert(IntVect::TheNodeVector());
        dom.grow(-1); // Shrink domain so we don't operate on any boundaries
        for (MFIter mfi(solution[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            bx.grow(nghost);        // Expand to cover first layer of ghost nodes
            bx = bx & dom;  // Take intersection of box and the problem domain

            Array4<Real> const& RHS  = rhs[ilev].array(mfi);
            for (int n = 0; n < op.ncomp; n++)
                ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    Real AMREX_D_DECL(x1 = i*dx + minx,
                                      x2 = j*dy + miny,
                                      x3 = k*dz + minz);

                    if (n==0) RHS(i,j,k,n) = AMREX_D_TERM(   (x1-0.5)*(x1+0.5),
                                                           * (x2-0.5)*(x2+0.5),
                                                           * (x3-0.5)*(x3+0.5));
                    else RHS(i,j,k,n) = 0.0;
                });
         }
        rhs[ilev].FillBoundary(false,true);
    }

    //
    // Set params to be passed to MLMG solver
    //
    LPInfo info;
    if (mlmg.agglomeration >= 0)        info.setAgglomeration(mlmg.agglomeration);
    if (mlmg.consolidation >= 0)        info.setConsolidation(mlmg.consolidation);
//    std::cout << "MAX CRS LEVEL = " << mlmg.max_coarsening_level << std::endl;
    if (mlmg.max_coarsening_level >= 0) info.setMaxCoarseningLevel(mlmg.max_coarsening_level);

    //
    // Initialize the MCNodalLinOp linear operator
    // (see ./MCNodalLinOp.cpp, ./MCNodalLinOp.H for implementation)
    //
    MCNodalLinOp linop;
    linop.setNComp(op.ncomp);
    linop.setCoeff(op.coeff);
    linop.define(geom,cgrids,dmap,mesh.ref_ratio,info);
    linop.setDomainBC({AMREX_D_DECL(amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet)},
                      {AMREX_D_DECL(amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet)});
    for (int ilev = 0; ilev < mesh.nlevels; ++ilev) linop.setLevelBC(ilev,&solution[ilev]);

    //
    // Initialize the MLMG solver
    //
    MLMG solver(linop);
    if (mlmg.verbose >= 0)     solver.setVerbose(mlmg.verbose);
    if (mlmg.bottom_verbose >= 0)  solver.setBottomVerbose(mlmg.bottom_verbose);
    if (mlmg.fixed_iter >= 0)  solver.setFixedIter(mlmg.fixed_iter);
    if (mlmg.max_iter >= 0)    solver.setMaxIter(mlmg.max_iter);
    if (mlmg.max_fmg_iter >= 0)solver.setMaxFmgIter(mlmg.max_fmg_iter);
    if (mlmg.pre_smooth >= 0) solver.setPreSmooth(mlmg.pre_smooth);
    if (mlmg.post_smooth >= 0) solver.setPostSmooth(mlmg.post_smooth);
    // IMPORTANT! Use the "CFStrategy::ghostnodes" strategy to avoid
    // having to implement a complicated "reflux" routine!
    solver.setCFStrategy(MLMG::CFStrategy::ghostnodes);

    //
    // Perform the solve
    //
    Real tol_rel = 1E-8, tol_abs = 1E-8;
    solver.solve(GetVecOfPtrs(solution),GetVecOfConstPtrs(rhs),tol_rel,tol_abs);
    solver.compResidual(GetVecOfPtrs(res),GetVecOfPtrs(solution),GetVecOfConstPtrs(rhs));
    solver.apply(GetVecOfPtrs(b),GetVecOfPtrs(solution));

    //
    // Write the output to ./solution
    //
    WriteMLMF ("solution",GetVecOfConstPtrs(solution),geom);
    WriteMLMF ("residual",GetVecOfConstPtrs(res),geom);

    for (int i = 0; i < resgn.size(); i++)
    {
        amrex::MultiFab::Copy(solgn[i],solution[i],0,0,op.ncomp,0); // Dx = x
        amrex::MultiFab::Copy(rhsgn[i],rhs[i],0,0,op.ncomp,0); // Dx = x
        amrex::MultiFab::Copy(resgn[i],res[i],0,0,op.ncomp,0); // Dx = x
        amrex::MultiFab::Copy(bgn[i],b[i],0,0,op.ncomp,0); // Dx = x
        
        //resgn[i].copy();
    }
    WriteMLMF ("solgn",GetVecOfConstPtrs(solgn),geom);
    WriteMLMF ("rhsgn",GetVecOfConstPtrs(rhsgn),geom);
    WriteMLMF ("residualgn",GetVecOfConstPtrs(resgn),geom);
    WriteMLMF ("bgn",GetVecOfConstPtrs(bgn),geom);    

    }
    Finalize();
}

