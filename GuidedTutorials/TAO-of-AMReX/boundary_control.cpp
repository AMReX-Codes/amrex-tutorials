#include <chrono>
using namespace std::chrono;

#include <AMReX.H>
#include "MyTest.H"

#include <petsc.h>

amrex::Real TargetSolution(amrex::Real* coords);
PetscErrorCode FormFunction(Tao tao, Vec P, PetscReal *f, void *ptr);
PetscErrorCode FormFunctionGradient(Tao tao, Vec P, PetscReal *f, Vec G, void *ptr);

int main (int argc, char* argv[])
{
    PetscErrorCode     ierr;
    PetscReal          one = 1.0;
    int                nb, nl, nt, n;
    int                Nb, Nl, Nt, N;
    int                no_args = 2;
    Vec                P;
    Tao                tao;
    PetscBool          flg, fd=PETSC_FALSE, plot=PETSC_FALSE;
    PetscMPIInt        rank;

    amrex::Initialize(no_args, argv);
    ierr = PetscInitialize(&argc, &argv, (char*)0, (char*)0); if (ierr) return ierr;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    
    {
    BL_PROFILE("main");
    MyTest mytest;

    // Get domain sizing and gradient type from command line options
    ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&mytest.n_cell,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL,NULL,"-fd",&fd,&flg);CHKERRQ(ierr);
    if (fd) mytest.fd_grad = true;
    ierr = PetscOptionsGetBool(NULL,NULL,"-plot",&plot,&flg);CHKERRQ(ierr);
    if (plot) mytest.plot = true;

    // Prep the solver and set the target solution we want to recover on the right edge of the domain
    // u_target = 4.0*(y-0.5)^2;
    mytest.initData();
    mytest.set_target_solution(TargetSolution);

    // Create PETSc vector for the bottom, left and top Dirichlet boundaries
    mytest.get_number_local_bcs(mytest.nb, mytest.nl, mytest.nt);
    n = mytest.nb + mytest.nl + mytest.nt; // total local size
    mytest.get_number_global_bcs(mytest.Nb, mytest.Nl, mytest.Nt);
    N = mytest.Nb + mytest.Nl + mytest.Nt; // total global size
    ierr = VecCreateMPI(PETSC_COMM_WORLD, n, N, &P);CHKERRQ(ierr);
    ierr = VecSet(P, one); CHKERRQ(ierr);

    // Create and setup the TAO optimization algorithm
    ierr = TaoCreate(PETSC_COMM_WORLD, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOBQNLS); CHKERRQ(ierr); // TAOBQNLS is a bound-constrained quasi-Newton linesearch alg,
    ierr = TaoSetInitialVector(tao, P); CHKERRQ(ierr);
    if (!mytest.fd_grad) {
        // Sets the callback for the gradient computed with the adjoint method
        ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient, &mytest); CHKERRQ(ierr);
    } else {
        // TAO provides a default gradient function for the central finite-difference formula
        ierr = TaoSetObjectiveRoutine(tao, FormFunction, &mytest);CHKERRQ(ierr);
        ierr = TaoSetGradientRoutine(tao, TaoDefaultComputeGradient, &mytest); CHKERRQ(ierr);
    }
    
    // This lets TAO read command line options
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

    // Start the optimization solution
    auto start = high_resolution_clock::now();
    ierr = TaoSolve(tao); CHKERRQ(ierr);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    if (rank == 0) std::cout << "TaoSolve() duration: " << duration.count() << " microseconds" << std::endl;
    }

    // Cleanup and exit
    ierr = TaoDestroy(&tao);CHKERRQ(ierr);
    ierr = VecDestroy(&P);CHKERRQ(ierr);
    ierr = PetscFinalize();
    amrex::Finalize();
}

/* -------------------------------------------------------------------- */

amrex::Real TargetSolution(amrex::Real* coords)
{
    amrex::Real y = coords[1];
    amrex::Real utarg = 4.*(y-0.5)*(y-0.5) - 0.5;
    return utarg;
}

/* -------------------------------------------------------------------- */

// This function computes only the objective value
PetscErrorCode FormFunction(Tao tao, Vec P, PetscReal *f, void *ptr)
{
    MyTest            *mytest = (MyTest *) ptr;
    PetscErrorCode    ierr;
    PetscReal         ff=0;
    const PetscScalar *pp;
    amrex::Real       *pb, *pl, *pt;

    PetscFunctionBeginUser;

    // Allocate arrays for bottom, left and top Dirichlet boundaries
    ierr = PetscMalloc1(mytest->nb * sizeof(amrex::Real), &pb);CHKERRQ(ierr);
    ierr = PetscMalloc1(mytest->nl * sizeof(amrex::Real), &pl);CHKERRQ(ierr);
    ierr = PetscMalloc1(mytest->nt * sizeof(amrex::Real), &pt);CHKERRQ(ierr);

    // Split the PETSc vector of optimization variables into its components
    ierr = VecGetArrayRead(P, &pp);CHKERRQ(ierr);
    for (int i=0; i<mytest->nb; i++) pb[i] = pp[i];
    for (int i=0; i<mytest->nl; i++) pl[i] = pp[mytest->nb + i];
    for (int i=0; i<mytest->nt; i++) pt[i] = pp[mytest->nb + mytest->nl + i];
    ierr = VecRestoreArrayRead(P, &pp);CHKERRQ(ierr);

    // Set the boundary values from TAO into the AMReX solver
    mytest->update_boundary_values(mytest->nb, pb, mytest->nl, pl, mytest->nt, pt);

    // Clean-up
    ierr = PetscFree(pb);CHKERRQ(ierr);
    ierr = PetscFree(pl);CHKERRQ(ierr);
    ierr = PetscFree(pt);CHKERRQ(ierr);

    /// Solve the Laplace equations with the prescribed boundary values
    mytest->solve();

    // Compute the objective function using the AMReX solution
    // f = 0.5 * \int_0^1 (u(1, y) - u_targ)^2 dy
    ff = mytest->calculate_obj_val();
    *f = ff;

    if (mytest->fd_grad && mytest->plot) {
        // Perform other misc operations like visualization and I/O
        mytest->write_plotfile();
        mytest->update_counter();
    }

    PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------- */

PetscErrorCode FormFunctionGradient(Tao tao, Vec P, PetscReal *f, Vec G, void *ptr)
{
    MyTest            *mytest = (MyTest *) ptr;
    PetscErrorCode    ierr;
    amrex::Real       *gb, *gl, *gt;
    PetscScalar       *gg;

    PetscFunctionBeginUser;

    // Evaluate the objective function
    ierr = FormFunction(tao, P, f, ptr);CHKERRQ(ierr);

    // Allocate arrays for bottom, left and top Dirichlet boundaries
    ierr = PetscMalloc1(mytest->nb * sizeof(amrex::Real), &gb);CHKERRQ(ierr);
    ierr = PetscMalloc1(mytest->nl * sizeof(amrex::Real), &gl);CHKERRQ(ierr);
    ierr = PetscMalloc1(mytest->nt * sizeof(amrex::Real), &gt);CHKERRQ(ierr);

    // Solve the adjoint problem and assemble the gradient for each boundary
    mytest->setup_adjoint_system(); 
    mytest->solve_adjoint_system();
    mytest->calculate_opt_gradient(mytest->nb, gb, mytest->nl, gl, mytest->nt, gt);

    // Write the AMReX gradient data into the PETSc vector
    ierr = VecGetArray(G, &gg);CHKERRQ(ierr);
    for (int i=0; i<mytest->nb; i++) gg[i] = gb[i];
    for (int i=0; i<mytest->nl; i++) gg[mytest->nb + i] = gl[i];
    for (int i=0; i<mytest->nt; i++) gg[mytest->nb + mytest->nl + i] = gt[i];
    ierr = VecRestoreArray(G, &gg);CHKERRQ(ierr);
    ierr = VecScale(G, mytest->n_cell*mytest->n_cell);CHKERRQ(ierr);

    // Clean-up
    ierr = PetscFree(gb);CHKERRQ(ierr);
    ierr = PetscFree(gl);CHKERRQ(ierr);
    ierr = PetscFree(gt);CHKERRQ(ierr);

    // Perform other misc operations like visualization and I/O
    if (mytest->plot) {
        mytest->write_plotfile();
        mytest->update_counter();
    }

    PetscFunctionReturn(0);
}
