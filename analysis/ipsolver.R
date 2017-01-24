# Simple implementation of a primal-dual interior-point solver for
# convex programs with convex inequality constraints (it does not
# handle equality constraints). Precisely speaking, it will compute
# the solution to the following optimization problem:
#
#     minimize    f(x)
#     subject to  c(x) < 0,
#
# where f(x) is a convex objective and c(x) is a vector-valued
# function with outputs that are convex in x. This code is mostly
# based on the descriptions provided in this reference:
#
#    Paul Armand, Jean C. Gilbert, Sophie Jan-Jegou. A Feasible BFGS
#    Interior Point Algorithm for Solving Convex Minimization Problems. 
#    SIAM Journal on Optimization, Vol. 11, No. 1, pp. 199-222.
#
# Input x0 is the initial point for the solver. It must be an n x 1
# matrix, where n is the number of (primal) optimization variables.
# DESCENTDIR must be either: 'newton' for the Newton search direction,
# 'bfgs' for the quasi-Newton search direction with the
# Broyden-Fletcher-Goldfarb-Shanno (BFGS) update, or 'steepest' for the
# steepest descent direction. The steepest descent direction is often quite
# bad, and the solver may fail to converge to the solution if you take this
# option. For the Newton direction, you must be able to compute the the
# Hessian of the objective. Also note that we form a quasi-Newton
# approximation to the objective, not to the Lagrangian (as is usually
# done). This means that you will always have to provide second-order
# information about the inequality constraint functions.
# 
# TOL is the tolerance of the convergence criterion; it determines when the
# solver should stop. MAXITER is the maximum number of iterations. And the
# final input, VERBOSE, must be set to true or false depending on whether
# you would like to see the progress of the solver.
#
# The inputs OBJ, GRAD, CONSTR and JACOBIAN must all be function handles. If
# you don't know what function handles are, type HELP FUNCTION_HANDLE in
# MATLAB.
# 
#    * OBJ must be a handle to a function that takes 1 input, the vector
#      of optimization variables, and returns the value of the function at
#      the given point. The function definition should look something like F
#      = OBJECTIVE(X).
#
#    * GRAD is a pointer to a function of the form G = GRADIENT(X), where
#      G is the n x 1 gradient vector of the objective, or [G H] =
#      GRADIENT(X) if the Newton step is used, in which case H is the n x n
#      Hessian of the objective.
#
#    * CONSTR is a handle to a function of the form C = CONSTRAINTS(X),
#      where C is the m x 1 vector of constraint responses at X.
#
#    * JACOBIAN is a handle to a function of the form [J W] =
#      JACOBIAN(X,Z). The inputs are the primal variables X and the m x 1
#      vector of dual variables Z. The return values are the m x n
#      Jacobian matrix (containing the first-order partial derivatives
#      of the inequality constraint functions), and W is the n x n
#      Hessian of the Lagrangian (minus the Hessian of the objective), 
#      which is basically equal to
#
#          W = z(1)*W1 + z(2)*W2 + ... + z(m)*Wm,
#
#      where Wi is the Hessian of the ith constraint.
#
# If you set VERBOSE to true, then at each iteration the solver will output
# the following information (from left to right): 1. the iteration number,
# 2. the value of the objective, 3. the barrier parameter mu, 4. the
# centering parameter sigma, 4. the residuals of the perturbed
# Karush-Kuhn-Tucker system (rx, rc), 5. the step size, and the number of
# iterations in the line search before we found a suitable descent step.
#
# If your optimization problem is large (i.e. it involves a lot of
# optimization variables or inequality constraints) it might speed up the
# solver to output sparse matrices. (Type HELP SPARSE in the MATLAB console
# for more information on sparse matrices.)
#
# As a final note, the interior-point solver may not work very well if your
# problem is very poorly scaled (i.e. the Hessian of the objective or the
# Hessian of one of the constraint functions is poorly conditioned). It
# is up to you to make sure you look at the conditioning of your problem.
ipsolver <- function (x, objective, gradient, constraints, jacobian, 
                      descentdir, tolerance, maxiter, verbose)

  # Some algorithm parameters.
  eps      <- 1e-8   # A number close to zero.
  sigmamax <- 0.5    # The maximum centering parameter.
  etamax   <- 0.25   # The maximum forcing number.
  mumin    <- 1e-9   # Minimum barrier parameter.
  alphamax <- 0.995  # Maximum step size.
  alphamin <- 1e-6   # Minimum step size.
  beta     <- 0.75   # Granularity of backtracking search.
  tau      <- 0.01   # Amount of actual decrease we will accept in 
                      # line search.

  # INITIALIZATION
  # Get the number of primal variables (n), the number of constraints
  # (m), and the total number of primal-dual optimization variables
  # (nv). Initialize the Lagrange multipliers. Initialize the
  # second-order information.
  c  <- constraints(x)
  n  <- length(x)
  m  <- length(c)
  nv <- n + m
  z  <- ones(m,1)
  
  
  # Repeat while the convergence criterion has not been satisfied, and
  # we haven't reached the maximum number of iterations.
  alpha = 0;
  ls    = 0;
  # TO DO: Explain what these columns mean.
  cat("  i f(x)       lg(mu) sigma   ||rx||  ||rc||  alpha   #ls\n")
  for (iter in seq(1,maxiter)) {

    # COMPUTE OBJECTIVE, GRADIENT, CONSTRAINTS, ETC.  
    # Compute the response of the objective function, the gradient of the
    # objective, the response of the inequality constraints, the Jacobian of
    # the inequality constraints, the Hessian of the Lagrangian (minus the
    # Hessian of the objective) and, optionally, the Hessian of the
    # objective.
    f     = objective(x);
    c     = constraints(x);
    [J W] = jacobian(x,z);
    if strcmp(descentdir,'newton')
      [g B] = gradient(x);
    else
      g = gradient(x);
    end
    
    % Compute the responses of the unperturbed Karush-Kuhn-Tucker
    % optimality conditions.
    rx = g + J'*z;  % Dual residual.
    rc = c.*z;      % Complementarity.
    r0 = [rx; rc]; 
    
    % Set some parameters that affect convergence of the primal-dual
    % interior-point method.
    eta        = min(etamax,norm(r0)/nv);
    sigma      = min(sigmamax,sqrt(norm(r0)/nv));
    dualitygap = -c'*z;
    mu         = max(mumin,sigma*dualitygap/m);
    
    % Print the status of the algorithm.
    if verbose
      fprintf('%3d %+0.3e  %+5.2f %0.1e %0.1e %0.1e %0.1e %3d\n',...
	      iter,f,log10(mu),sigma,norm(rx),norm(rc),alpha,ls);
    end

    % CONVERGENCE CHECK.
    % If the norm of the responses is less than the specified tolerance,
    % we are done. 
    if norm(r0)/nv < tolerance
      break
    end
    
    % Update the BFGS approximation to the Hessian of the objective.
    if strcmp(descentdir,'bfgs') & iter > 1
      B = bfgsupdate(B,alpha*px,g-gprev);
    end

    % SOLUTION TO PERTURBED KKT SYSTEM.
    % Compute the search direction of x and z.
    S  = diag(sparse(z./(c-eps)));
    gb = g - mu*J'*(1./(c-eps));
    px = (B + W - J'*S*J) \ (-gb);
    pz = -(z + mu./(c-eps) + S*J*px);
    
    % BACKTRACKING LINE SEARCH.
    % To ensure global convergence, execute backtracking line search to
    % determine the step length. First, we have to find the largest step
    % size which ensures that z remains feasible. Next, we perform
    % backtracking line search.
    alpha = alphamax;
    is    = find(z + pz < 0);
    if length(is)
      alpha = alphamax * min(1,min(z(is) ./ -pz(is)));
    end
    
    % Compute the response of the merit function and the directional
    % gradient at the current point and search direction.
    psi  = merit(x,z,f,c,mu,eps);
    dpsi = gradmerit(x,z,px,pz,g,c,J,mu,eps);
    ls   = 0;
    while true

      % Compute the candidate point, the constraints, and the response of
      % the objective function and merit function at the candidate point.
      ls     = ls + 1;
      xnew   = x + alpha * px;
      znew   = z + alpha * pz;
      f      = objective(xnew);
      c      = constraints(xnew);
      psinew = merit(xnew,znew,f,c,mu,eps);
      
      % Stop backtracking search if we've found a candidate point that
      % sufficiently decreases the merit function and satisfies all the
      % constraints.
      if sum(c > 0) == 0 & psinew < psi + tau*eta*alpha*dpsi
	x     = xnew;
	z     = znew;
	gprev = g;
	break
      end
      
      % The candidate point does not meet our criteria, so decrease the step
      % size for 0 < beta < 1.
      alpha = alpha * beta;
      if alpha < alphamin
	error('Step size too small');
      end
    end
  end
  
% ------------------------------------------------------------------
% Compute the response of the merit function at (x,z).
function psi = merit (x, z, f, c, mu, eps)
  psi = f - c'*z - mu*sum(log(c.^2.*z+eps));
  
% ------------------------------------------------------------------
% Compute the directional derivative of the merit function at (x,z).
function dpsi = gradmerit (x, z, px, pz, g, c, J, mu, eps)
  dpsi = px'*(g - J'*z - 2*mu*J'*(1./(c-eps))) - pz'*(c + mu./(z+eps));

% ------------------------------------------------------------------
% Update the quasi-Newton approximation using Broyden-Fletcher-
% Goldfarb-Shanno (BFGS) formula.
function B = bfgsupdate (B, s, y)  
  if y'*s < 0
    error('dot(y,s) > 0 is not satisfied');
  end
  x = B*s;
  B = B - x*x'/(x'*s) + y*y'/(y'*s);
  
