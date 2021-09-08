# Direct_Shooting_Simulation-_Single_and_Multiple_Shooting

This is an optimal control solving algorithm  using MATLAB's function handle and optimization toolbox. Basic direct shooting algorithm is reformulated, its main idea to discrete control (DSS), or control and state (DMS) to get global optimal input of dynamical system.  These code could work on low-dimension simple OCP, with free terminal time $t_f$, state-control constrained $C(x(t),u(t),t)\leq 0$.

 [Git link](git@github.com:OdellDotson/Direct_Shooting_Simulation-_Single_and_Multiple_Shooting.git)

## problem

A simulation tool that uses direct single shooting and direct multiple shooting state control systems. 

- $u\in\Reals^q,x\in\Reals^n$

- dynamical equation constraint：
  $$
  \dot x=f(x,u,t)
  $$

- path constraint in inequality
  $$
  g(x,u,t)\le 0
  $$

- boundary constraint as equality
  $$
  c(x(t_0),x(t_f),u(t_0),u(t_f))=0
  $$

- mixed performance index OR Bolza's objective
  $$
  min_{u(t)}J:=\phi(x,t)+\int_{t_0}^{t_f}L(x,u,t)\text d t
  $$

- terminal time $t_f$ free or pre-defined, the latter is transformed into terminal constraint

## algorithm

DSS algorithm solves a problem in this form:

1. discrete $u(t)$, into column vector $(p*N)$ by 1, control is assumed constant in each segment
2. in each segment, $x(t)=\Phi(x,u,t)$ are solved by RK4 integration
3. no gradient info needed
4. 



DMS algorithm solves a problem in this form:

1. $N$ segment to discrete control-state variables
2. discrete $u(t)$, into column vector $(p*N)$ by 1, control is assumed constant in each segment
3.  $x(t)$ are also discretized in each segment, and dynamical constraint is emphasized by constraint and solving $x(k+1)-\Phi(x(k),u_k,t)=0$ this equation by RK4 integration
4. gradient info needed or not, but if provided, better convergence will be when used to highly nonlinear system
5. 



## Examples to test

### example 1: inversed pendulum

initial state known$x_0=[\theta,\omega]^{\mathrm T}=[\pi,0]$，terminal state constraint$[\theta_f,\omega_f]^{\mathrm T}=[0,0]$，terminal time $t_f$free ，and a quadratic criterion
$$
\min_{u(t)}=\frac1 2\mathbf x_f^{\mathrm T}\text Q\mathbf x_f+\frac1 2\int_0^{t_f}\text R\mathbf u^2\\
\text{s.t.}
\begin{bmatrix}\dot\theta\\\dot\omega\end{bmatrix}=\begin{bmatrix}\omega\\1/I(-mgl\sin\theta+b\omega+u)\end{bmatrix}
$$


### examples 2:  moon lander in 1D

a moon lander pinpoint landing problem, terminal time $t_f$ unknown. To find thrust profile $T(t)$ subject to:
$$
\left\{\begin{aligned}
\dot{h}(t) &=v(t) \\
\dot{v}(t) &=-g+\frac{\alpha(t)}{m(t)} \\
\dot{m}(t) &=-\frac{\alpha(t)}{g_e I_{sp}}
\end{aligned}\right.
$$

$$
\max_{T(t)}{ -m_f}\\
s.t. \Psi(\mathbf{x}(t_f),t_f)=[h,v]|_{t_f}=0\\
a\in[1500,7500]N
$$

This is a simple problem and you can test your algorithm. DSS and DMS are all well suited.

### example 3: fastest turning of ship

To find the control angle $\alpha(t)$ of ship and to satisfy terminal state constraints.
$$
\left\{\begin{aligned}
\dot{x}(t) &=v_x(t) \\
\dot{y}(t) &=v_y(t) \\
\dot{v}_x(t) &=3\cos\alpha(t) \\
\dot{v}_y(t) &=3\sin\alpha(t) 
\end{aligned}\right.
$$

$$
\min_{\alpha(t)}{ t_f}\\
s.t. x(0)=y(0)=v_x(0)=0,v_y(0)=5\\
x(t_f)=y(t_f)=v_y(t_f)=5,v_x(0)=0\\
\alpha\in[-\pi,\pi]rad
$$

As the convergence problem is a little poor, DSS algorithm has to be set with initial value very close to optimal, while DMS is suited to solve it.