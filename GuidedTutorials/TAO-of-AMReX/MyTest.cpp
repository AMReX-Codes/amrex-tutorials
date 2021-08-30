#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>
#include <AMReX_BCUtil.H>
#include <cassert>

using namespace amrex;

namespace ExtTaoBC
{
    std::map<BoxCornerTuple, BoxBoundaryValues> ext_dir_bcs;
}

BoxCornerTuple make_box_corner_tuple(const amrex::Box& bx)
{
    const auto bx_lo = amrex::lbound(bx);
    const auto bx_hi = amrex::ubound(bx);
    
    return std::make_tuple(bx_lo.x, bx_lo.y, bx_lo.z,
                           bx_hi.x, bx_hi.y, bx_hi.z);
}

MyTest::MyTest()
{
    readParameters();
}

void MyTest::solve()
{
    solvePoisson(solution, rhs);
}

void MyTest::update_counter()
{
    iteration_counter++;
}

std::string MyTest::get_iteration_filename(std::string filename)
{
    std::stringstream sstream;
    sstream << filename << std::setw(5) << std::setfill('0') << iteration_counter;
    std::string file_with_counter = sstream.str();
    return file_with_counter;
}

void MyTest::write_plotfile(bool minimal_plotfile)
{
    if (minimal_plotfile)
        writeMinimalPlotfile();
    else
        writePlotfile();
}

void MyTest::get_number_global_bcs(int& num_lower, int& num_left, int& num_upper)
{
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    num_lower = (domain_hi.x - domain_lo.x + 1); // lower, upper edges

    num_upper = (domain_hi.x - domain_lo.x + 1); // lower, upper edges

    num_left = domain_hi.y - domain_lo.y + 1; // left edge
    num_left += 2; // left/lower and left/upper corners
}

void MyTest::get_number_local_bcs(int& local_num_lower, int& local_num_left, int& local_num_upper)
{
    // Get number of boundary values local to this MPI rank
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    local_num_lower = 0;
    local_num_left = 0;
    local_num_upper = 0;

    int im = 0;
    for (MFIter mfi(solution, false); mfi.isValid(); ++mfi)
        {
            int num_lower = 0;
            int num_left = 0;
            int num_upper = 0;

            const Box& bx = mfi.validbox();
            const auto bx_lo = lbound(bx);
            const auto bx_hi = ubound(bx);

            im++;

            const Box& gbx = mfi.growntilebox();

            auto key = make_box_corner_tuple(gbx); 
            BoxBoundaryValues vvr;

            bool aligned_lower = false;
            bool aligned_left  = false;
            bool aligned_upper = false;

            // check lower boundary
            if (bx_lo.y == domain_lo.y) {
                aligned_lower = true;
                num_lower += (bx_hi.x - bx_lo.x + 1);
            }

            // check left boundary
            if (bx_lo.x == domain_lo.x) {
                aligned_left = true;
                num_left += (bx_hi.y - bx_lo.y + 1);
            }

            // check upper boundary
            if (bx_hi.y == domain_hi.y) {
                aligned_upper = true;
                num_upper += (bx_hi.x - bx_lo.x + 1);
            }

            // add corners to left edge
            if (aligned_left && aligned_lower)
                num_left++;

            if (aligned_left && aligned_upper)
                num_left++;

            // add vectors to ext_dir_bcs mapped to the grown box
            EdgeBoundaryValues xbc_lower; xbc_lower.resize(num_lower);
            vvr.push_back(xbc_lower);

            EdgeBoundaryValues xbc_left; xbc_left.resize(num_left);
            vvr.push_back(xbc_left);

            EdgeBoundaryValues xbc_upper; xbc_upper.resize(num_upper);
            vvr.push_back(xbc_upper);

            auto key_val = std::make_pair(key, vvr);
            ExtTaoBC::ext_dir_bcs.insert(key_val); 

            local_num_lower += num_lower;
            local_num_left += num_left;
            local_num_upper += num_upper;
        }
}


void MyTest::update_boundary_values(int nlower, const Real *xlower,
                                    int nleft, const Real *xleft,
                                    int nupper, const Real *xupper)
{
    int ilower = 0;
    int ileft = 0;
    int iupper = 0;

    auto iter = ExtTaoBC::ext_dir_bcs.begin();
    int iiter = 0;
    while (iter != ExtTaoBC::ext_dir_bcs.end())
    {
        auto& vvr = iter->second;

        const int size_lower = vvr(ExtTaoBC::lower_boundary).size();
        const int size_left = vvr(ExtTaoBC::left_boundary).size();
        const int size_upper = vvr(ExtTaoBC::upper_boundary).size();

        iiter++;

        for (int i = 0; i < size_lower; ++i)
        {
            vvr(ExtTaoBC::lower_boundary, i) = xlower[ilower + i];
            assert(ilower + i <= nlower);
        }

        for (int i = 0; i < size_left; ++i)
        {
            vvr(ExtTaoBC::left_boundary, i) = xleft[ileft + i];
            assert(ileft + i <= nleft);
        }

        for (int i = 0; i < size_upper; ++i)
        {
            vvr(ExtTaoBC::upper_boundary, i) = xupper[iupper + i];
            assert(iupper + i <= nupper);
        }

        ilower += size_lower;
        ileft += size_left;
        iupper += size_upper;

        iter++;
    }

    solution.FillBoundary(geom.periodicity());
    FillDomainBoundary(solution, geom, {bcs});
}

void MyTest::setup_adjoint_system()
{
    // setup the (dR/du)^T * lambda = - \partial f/\partial u linear system
    // by evaluating \partial f/\partial u = \int_V (u(p) - u_t) dV

    adjoint_rhs = 0.0;
    adjoint = 0.0;

    // for right boundary, adjoint_rhs(cell) = -dfdu = target solution(cell) - poisson solution(cell)
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    const int k = 0;

    for (MFIter mfi(solution, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto bx_lo = lbound(bx);
        const auto bx_hi = ubound(bx);

        const auto sol_arr = solution[mfi].array();
        const auto exact_sol_arr = exact_solution[mfi].array();
        auto adj_arr = adjoint_rhs[mfi].array();

        // check if we have part of the right boundary
        if (bx_hi.x == domain_hi.x) {
            const int i = bx_hi.x;
            for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                adj_arr(i, j, k) = exact_sol_arr(i, j, k) - sol_arr(i, j, k);
                adj_arr(i, j, k) *= AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
            }
        }
    }
}

void MyTest::solve_adjoint_system()
{
    // solve (dR/du)^T * lambda = - \partial f/\partial u

    solvePoisson(adjoint, adjoint_rhs);
}

void MyTest::set_target_solution(Real (*ftarget)(Real* coords))
{
    target_function = ftarget;

    update_target_solution();
}

void MyTest::update_target_solution()
{
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    Copy(exact_solution, solution, 0, 0, 1, 0);

    const int k = 0;

    for (MFIter mfi(exact_solution, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto bx_lo = lbound(bx);
        const auto bx_hi = ubound(bx);

        auto exact_sol_arr = exact_solution[mfi].array();

        if (bx_hi.x == domain_hi.x) {
            const int i = bx_hi.x;
            for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                IntVect cell_indices;
                AMREX_D_TERM(cell_indices[0] = i;, cell_indices[1] = j;, cell_indices[2] = k;)
                Vector<Real> cell_location(AMREX_SPACEDIM);
                geom.CellCenter(cell_indices, cell_location);
                exact_sol_arr(i, j, k) = target_function(cell_location.dataPtr());
            }
        }
    }
}

Real MyTest::calculate_obj_val()
{
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    Real fobj_local = 0.0;

    const int k = 0;

    for (MFIter mfi(solution, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto bx_lo = lbound(bx);
        const auto bx_hi = ubound(bx);

        const auto sol_arr = solution[mfi].array();
        const auto exact_sol_arr = exact_solution[mfi].array();

        // check if we have part of the right boundary
        if (bx_hi.x == domain_hi.x) {
            const int i = bx_hi.x;
            // loop over right edge, excluding corners
            for (int j = std::max(bx_lo.y, domain_lo.y+1); j <= std::min(bx_hi.y, domain_hi.y-1); ++j) {
                fobj_local += 0.5 * std::pow(sol_arr(i, j, k) - exact_sol_arr(i, j, k), 2.0);
            }
        }
    }

    fobj_local *= AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));

    Real fobj_global = fobj_local;
    ParallelDescriptor::ReduceRealSum(fobj_global);

    return fobj_global;
}

void MyTest::calculate_opt_gradient(int nlower, Real* glower,
                                    int nleft, Real* gleft,
                                    int nupper, Real* gupper)
{
    // calculate df/dp - \partial f/\partial p = (dR/dp)^T * lambda
    const Box& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    for (int i = 0; i < nlower; ++i) glower[i] = 0.0;
    for (int i = 0; i < nleft; ++i) gleft[i] = 0.0;
    for (int i = 0; i < nupper; ++i) gupper[i] = 0.0;

    // indices into dfdp arrays from TAO
    int ilower = 0;
    int ileft = 0;
    int iupper = 0;

    const int k = 0;

    for (MFIter mfi(adjoint, false); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            const auto bx_lo = lbound(bx);
            const auto bx_hi = ubound(bx);

            bool aligned_lower = false;
            bool aligned_left  = false;
            bool aligned_upper = false;

            const auto adjoint_arr = adjoint[mfi].array();

            // For each iteration of the MFIter, we fill TAO arrays in the same
            // order as we set number of entries in get_number_local_bcs().

            // check lower
            if (bx_lo.y == domain_lo.y) {
                aligned_lower = true;
                const int j = bx_lo.y;
                for (int i = bx_lo.x; i <= bx_hi.x; ++i) {
                    glower[ilower] += adjoint_arr(i, j, k);
                    ilower++;
                }
            }

            // check left (not including corners)
            if (bx_lo.x == domain_lo.x) {
                aligned_left = true;
                const int i = bx_lo.x;
                for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                    gleft[ileft] += adjoint_arr(i, j, k);
                    ileft++;
                }
            }

            // check upper
            if (bx_hi.y == domain_hi.y) {
                aligned_upper = true;
                const int j = bx_hi.y;
                for (int i = bx_lo.x; i <= bx_hi.x; ++i) {
                    gupper[iupper] += adjoint_arr(i, j, k);
                    iupper++;
                }
            }

            // lower left corner
            if (aligned_left && aligned_lower) {
                gleft[ileft] += adjoint_arr(bx_lo.x, bx_lo.y, k);
                ileft++;
            }

            // upper left corner
            if (aligned_left && aligned_upper) {
                gleft[ileft] += adjoint_arr(bx_lo.x, bx_hi.y, k);
                ileft++;
            }
        }
}

void MyTest::solvePoisson(amrex::MultiFab &solution,
                          amrex::MultiFab &rhs)
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = mlmg_tol_rel;
    const Real tol_abs = mlmg_tol_abs;

    const int nlevels = 1;

    if (composite_solve)
    {

        MLPoisson mlpoisson({geom}, {grids}, {dmap}, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Neumann,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});

        const Real* dx = geom.CellSize();
        mlpoisson.setDomainBCLoc({AMREX_D_DECL(0.5*dx[0],0.5*dx[1],0.5*dx[2])},
                                 {AMREX_D_DECL(0.5*dx[0],0.5*dx[1],0.5*dx[2])});

        mlpoisson.setLevelBC(0, &solution);

        MLMG mlmg(mlpoisson);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre)
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc)
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);
    }
    else
    {
        MLPoisson mlpoisson({geom}, {grids}, {dmap}, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
            {AMREX_D_DECL(LinOpBCType::Neumann,
                          LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet)});

        mlpoisson.setLevelBC(0, &solution);

        MLMG mlmg(mlpoisson);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre)
            {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc)
            {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

        mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);
    }
}

void MyTest::readParameters()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("composite_solve", composite_solve);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("mlmg_tol_rel", mlmg_tol_rel);
    pp.query("mlmg_tol_abs", mlmg_tol_abs);

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    pp.query("hypre_interface", hypre_interface_i);
    if (hypre_interface_i == 1)
    {
        hypre_interface = Hypre::Interface::structed;
    }
    else if (hypre_interface_i == 2)
    {
        hypre_interface = Hypre::Interface::semi_structed;
    }
    else
    {
        hypre_interface = Hypre::Interface::ij;
    }
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc", use_petsc);
#endif

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(use_hypre && use_petsc),
                                     "use_hypre & use_petsc cannot be both true");
}

void MyTest::initData()
{
    int nlevels = 1;

    RealBox rb({AMREX_D_DECL(0., 0., 0.)}, {AMREX_D_DECL(1., 1., 1.)});
    Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0, 0, 0)}, IntVect{AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1)});
    Box domain = domain0;
    geom.define(domain);

    grids.define(domain);
    grids.maxSize(max_grid_size);
    dmap.define(grids);
    solution.define(grids, dmap, 1, ngrow);
    adjoint.define(grids, dmap, 1, ngrow);
    adjoint_rhs.define(grids, dmap, 1, ngrow);
    rhs.define(grids, dmap, 1, ngrow);
    exact_solution.define(grids, dmap, 1, 0);

    solution = 0.0;

    // set up boundary conditions:
    // - first set all boundary conditions to external Dirichlet
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        bcs.setLo(i, BCType::ext_dir);
        bcs.setHi(i, BCType::ext_dir);
    }

    // - then set outflow to right = first order extrapolation
    bcs.setHi(0, BCType::foextrap);

    initProbPoisson();
}
