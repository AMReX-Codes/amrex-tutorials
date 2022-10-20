This solves the scalar wave equation $\frac{\partial^2 u}{\partial t^2} = \nabla u$.
If we define $v = \frac{\partial u}{\partial t}$, then we have two
first-order partial differential equations, $\frac{\partial u}{\partial t} = v$
and $\frac{\partial v}{\partial t} = \nabla u$.  The Laplacian operator is
discretized dimension by dimension with a fourth-order stencil,

$$\frac{\partial^2 u}{\partial x^2} = \left(-\frac{5}{2} u_i + \frac{4}{3} (u_{i-1} + u_{i+1}) - \frac{1}{12} (u_{i-2} + u_{i+2})\right) / \Delta x^2$$

The time stepping is done with a Runge-Kutta method (RK2, RK3 or RK4).  In
this test, the displacement at the x-direction boundaries is zero, and the
it's periodic in the y-direction.  Note that refluxing is not implemented in
this test code.
