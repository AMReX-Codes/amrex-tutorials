#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>
#include <AMReX_Dim3.H>
#include <map>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

struct EdgeBoundaryValues
{
    amrex::Vector<amrex::Real> x;

    EdgeBoundaryValues() { x.resize(0); }

    amrex::Real& operator()(int i) { return x[static_cast<long>(i)]; }

    void resize(int n) { x.resize(n); }

    int size() { return x.size(); }

    void push_back(amrex::Real xi) { x.push_back(xi); }
};

struct BoxBoundaryValues
{
    amrex::Vector<EdgeBoundaryValues> x;

    BoxBoundaryValues() { x.resize(0); }

    EdgeBoundaryValues& operator()(int i) { return x[i]; }

    amrex::Real& operator()(int i, int j) { return x[i](j); }

    void resize(int n) { x.resize(n); }

    int size() { return x.size(); }

    void push_back(EdgeBoundaryValues xi) { x.push_back(xi); }
};

typedef std::tuple<int, int, int, int, int, int> BoxCornerTuple;

BoxCornerTuple make_box_corner_tuple(const amrex::Box& bx);

namespace ExtTaoBC
{
    extern std::map<BoxCornerTuple, BoxBoundaryValues> ext_dir_bcs;

    const int lower_boundary = 0;
    const int left_boundary = 1;
    const int upper_boundary = 2;
}

class MyTest
{
    public:

        MyTest ();
        void initData ();
        void solve ();
        void writePlotfile ();
        void writeMinimalPlotfile ();
        int n_cell = 128;
        int nb, nl, nt;
        int Nb, Nl, Nt;
        bool fd_grad = false;
	bool plot = false;

    public: // make these public for cuda
        void initProbPoisson ();
        void initProbABecLaplacian ();
        void get_number_global_bcs(int& num_lower, int& num_left, int& num_upper);
        void get_number_local_bcs(int& num_lower, int& num_left, int& num_upper);
        void set_target_solution(amrex::Real (*ftarget)(amrex::Real* coords));
        void update_boundary_values(int nb, const amrex::Real *xb, int nl, const amrex::Real *xl, int nt, const amrex::Real *xt);

        amrex::Real calculate_obj_val();
        void calculate_opt_gradient(int nlower, amrex::Real* glower,
                int nleft, amrex::Real* gleft,
                int nupper, amrex::Real* gupper);
        void setup_adjoint_system();
        void solve_adjoint_system();
        void write_plotfile(bool minimal=false);
        void update_counter();

    private:

        void readParameters ();
        void solvePoisson(amrex::MultiFab &solution, amrex::MultiFab &rhs);
        void solveABecLaplacian ();
        std::string get_iteration_filename(std::string filename);
        void update_target_solution();

        int max_level = 1;
        int ref_ratio = 2;
        int max_grid_size = 64;
        const int ngrow = 1; // number of ghost cells for solution

        bool composite_solve = true;

        int prob_type = 1;  // 1. Poisson,  2. ABecLaplacian

        int iteration_counter = 0;

        // For MLMG solver
        int verbose = 2;
        int bottom_verbose = 0;
        int max_iter = 100;
        int max_fmg_iter = 0;
        int linop_maxorder = 2;
        bool agglomeration = true;
        bool consolidation = true;
        int max_coarsening_level = 30;
        bool use_hypre = false;
        bool use_petsc = false;
        amrex::Real mlmg_tol_rel = 1.0e-12;
        amrex::Real mlmg_tol_abs = 1.0e-6;

#ifdef AMREX_USE_HYPRE
        int hypre_interface_i = 1;  // 1. structed, 2. semi-structed, 3. ij
        amrex::Hypre::Interface hypre_interface = amrex::Hypre::Interface::structed;
#endif

        amrex::BCRec bcs; // BCs for 1 component

        amrex::Geometry geom;
        amrex::BoxArray grids;
        amrex::DistributionMapping dmap;

        amrex::MultiFab solution;
        amrex::MultiFab rhs;
        amrex::MultiFab adjoint;
        amrex::MultiFab adjoint_rhs;
        amrex::MultiFab exact_solution;
        amrex::MultiFab acoef;
        amrex::MultiFab bcoef;

        amrex::Real ascalar = 1.e-3;
        amrex::Real bscalar = 1.0;

        amrex::Real (*target_function)(amrex::Real* physical_coordinates);
};

#endif
