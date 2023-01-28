template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
	    VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func)
/* Given an n-dimensional point xold[0..n-1], the value of the function and
 gradient there, fold and g[0..n-1], and a direction p[0..n-1], finds a new 
 point x[0..n-1] along the direction p from xold where the function or functor 
 func has decreased "sufficiently".The new function value is returnd in f. 
 stpmax is an input quantity that limits the length of the steps so that you do
 not try to evaluate the function in regions where it is undefined or subject 
 to overflow. p is usually the Newton direction. The output quantity check is 
 false on a normal exit. It is truw when x is too close to xold. In a 
 minimisation algorithm, this usually signals convergence and can be ignored. 
 However, in a zero-finding algorithm the calling program should check whether 
 the convergence is spurious. */
{
  // ALF ensures sufficient decrease in function value.
  const Doub ALF=1.0e-4;

  // TOLX is the convergence criterion on Deltax
  Doub TOLX;
  TOLX = numeric_limits<Doub>::epsilon();

  Doub a, alam, alam2=0.0, alamin, b, disc, f2=0.0;
  Doub rhs1, rhs2, slope=0.0, sum=0.0, temp, test, tmplam;
  Int i, n=xold.size();
  check=false;
  
  for (i=0;i<n;i++) sum += p[i]*p[i];
  sum = sqrt(sum);

  // Scale if attempted step is too big
  if (sum > stpmax) 
    for (i=0;i<n;i++)
      p[i] *= stpmax/sum;

  for (i=0;i<n;i++)
    slope += g[i]*p[i];

  if (slope >= 0.0) throw("Roundoff problem in lnsrch");

  // Compute lambda_min
  test = 0.0;
  for (i=0;i<n;i++) {
    temp = abs(p[i])/MAX(abs(xold[i]),1.0);
    if (temp > test) test = temp;
  }
  alamin = TOLX/test;

  // Always try a full Newton step first
  alam = 1.0;

  // Start of iteration loop.
  for (;;) {
    for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
    f = func(x);
    // Convergence on Delta x.
    // For zero finding, the calling program should verify the convergence.
    if (alam < alamin) {
      for (i=0;i<n;i++) x[i]=xold[i];
      check = true;
      return;
    } else if (f <= fold+ALF*alam*slope) return;
    // Sufficient function decrease.
    else {
      // Backtrack.
      if (alam == 1.0)
	// First time
	tmplam = -slope/(2.0*(f-fold-slope));
      else {
	// Subsequent backtracks.
	rhs1 = f-fold-alam*slope;
	rhs2 = f2-fold-alam2*slope;
	a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc = b*b-3.0*a*slope;
	  if (disc < 0.0) tmplam =0.5*alam;
	  else if (b <= 0.0) tmplam =(-b+sqrt(disc))/(3.0*a);
	  else tmplam=-slope/(b+sqrt(disc));
	}
	// lambda >= 0.5 lambda_1
	if (tmplam > 0.5*alam)
	  tmplam = 0.5*alam;
      }
    }
    
    alam2 = alam;
    f2 = f;
    alam = MAX(tmplam,0.1*alam);
    // lambda >= 0.1 lambda_1
    // Try again.
  }
}


template <class T>
struct NRfdjac {
  /* Computes forward difference approximation to Jacobian. */
  // EPS set to approximate square root of machine precision.                   
  const Doub EPS;
  T &func;

  /* Initialise with user supplied function or functor that returns the vector 
     of functions to be zeroed. */
NRfdjac(T &funcc) : EPS(1.0e-8),func(funcc) {}

  MatDoub operator() (VecDoub_I &x, VecDoub_I &fvec) {
    /* Returns the Jacobian array df[0..n-1][0..n-1]. On input, x[0..n-1] is    
       the point at which the Jacobian is to be evaluated and fvec[0..n-1] is   
       the vector of function values at the point. */
    Int n=x.size();
    MatDoub df(n,n);
    VecDoub xh=x;
    for (Int j=0;j<n;j++) {
      Doub temp = xh[j];
      Doub h = EPS*abs(temp);
      if (h == 0.0) h=EPS;
      // Trick to reduce finite precision error.                                
      xh[j] = temp+h;
      h = xh[j]-temp;
      VecDoub f = func(xh);
      xh[j] = temp;
      // Forward difference formula.                                            
      for (Int i=0;i<n;i++)
        df[i][j] = (f[i]-fvec[i])/h;
    }
    return df;
  }
};

template <class T>
struct NRfmin {
  /* Returns f = 1/2 F \dot F. Also stores values of F in fvec. */
  VecDoub fvec;
  T &func;
  Int n;

  /* Intialise with user supplied function or functor that returns the vector   
     of functions to be zeroed. */
NRfmin(T &funcc) : func(funcc){}

  Doub operator() (VecDoub_I &x) {
    /* Returns f at x and stores F(x) in fvec. */
    n=x.size();
    Doub sum=0.0;
    fvec = func(x);
    for (Int i=0;i<n;i++) sum += SQR(fvec[i]);
    return 0.5*sum;
  }
};

template <class T>
void newt(VecDoub_IO &x, Bool &check, T &vecfunc) 
/* Given an initial guess x[0..n-1] for a root in n dimensions, find the roots
   by a globally convergent Newton's method. The vector of functions to be 
   zeroed, called fvec[0..n-1] in the routine below, is returned by the user-
   supplied function or functor vecfunc (see text). The output quantity check 
   is false on a normal return and true if the routine has converged to a local
   minimum of the function fmin defined below. In this case try restarting from
   a different initial guess.

   newt assumes that the typical values of all components of x and F are of 
   order unity, and it will fail if this assumption is badly violated. You 
   should rescale variables before invoking newt if this occurs.

   Here MAXITS is the maximum number of iterations; TOLF sets the convergence 
   criterion on function values; TOLMIN sets the criterion for deciding whether
   spurious conergence to a minimum of fmin has occured; STPMX is the scaled
   maximum step length allowed in the searches; and TOLX is the convergence 
   criterion on delta x. */
{
  const Int MAXITS=200;
  const Doub TOLF=1.0e-8,TOLMIN=1.02-12,STPMX=100.0;
  const Doub TOLX=numeric_limits<Doub>::epsilon();

  Int i,j,its,n=x.size();
  Doub den,f,fold,stpmax,sum,temp,test;
  VecDoub g(n),p(n),xold(n);
  MatDoub fjac(n,n);

  // Declare NRfmin object
  NRfmin<T> fmin(vecfunc);
 
  // Declare NRfdjac object
  NRfdjac<T> fdjac(vecfunc);

  // Make an alias to simplify coding
  VecDoub &fvec=fmin.fvec;
  
  // fvec is also computed by this call
  f = fmin(x);
  
  // Test for initial guess being a root, use more stringent test than TOLF.
  test = 0.0;
  for (i=0;i<n;i++) 
    if (abs(fvec[i]) > test) test=abs(fvec[i]);
  if (test < 0.01*TOLF) {
    check = false;
    return;
  }

  // Calculate stpmax for line searches.
  sum = 0.0;
  for (i=0;i<n;i++) sum += SQR(x[i]);
  stpmax = STPMX*MAX(sqrt(sum),Doub(n));

  // Start of iteration loop.
  for (its=0;its<MAXITS;its++) {
    fjac = fdjac(x,fvec);
    
    /* If analytic Jacobian is available, you can replace the struct NRfdjac 
       below with your own struct. */
    // Computer Delta f for the line search.
    for (i=0;i<n;i++) {
      sum = 0.0;
      for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    
    // Store x, and f
    for (i=0;i<n;i++) xold[i]=x[i];
    fold = f;
    
    // RHS for linear equations
    for (i=0;i<n;i++) p[i]=-fvec[i];
    
    // Solve linear equations by LU decomposition.
    LUdcmp alu(fjac);
    alu.solve(p,p);
    
    /* lnsrch returns new x and f. It also calculates fvec at new the x when it
       calls fin. */
    lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);

    // Test for convergence on function values.
    test=0.0;
    for (i=0;i<n;i++)
      if (abs(fvec[i]) > test) test=abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      return;
    }

    /* Check for gradient of f zero, i.e. spurious convergence. The relavant 
       gradient is estimated as (deltaf/f)/(deltax/x) 
       where deltaf ~= Nablaf deltax. */
    if (check) {
      test = 0.0;
      den = MAX(f,0.5*n);
      for (i=0;i<n;i++) {
	temp = abs(g[i]) * MAX(abs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      check = (test < TOLMIN);
      return;
    }

    // Test for convergence on deltax
    test=0.0;
    for (i=0;i<n;i++) {
      temp = (abs(x[i]-xold[i]))/MAX(abs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX)
      return;
  }

  throw("MAXITS exceeded in newt");
}

template <class T>
void broydn(VecDoub_IO &x, Bool &check, T &vecfunc)
/* Given an inital guess x[0..n-1] for a root in n dimensions, find the root by
   Broyden's method embedded in a globally convergent strategy. The vector of 
   functions to be zeroed, called fvec[0..n-1] in the routine below, is 
   returned by the user-supplied function or functor vecfunc. The routines 
   NRfdjac and NRfmin from newt are used. The output quantity check is false on
   a normal return and true if the routine had converged to a local minimum of
   the function fmin or if Broyden's method can make no further progress. In 
   this case try restarting from a different initial guess. 

   Here MAXITS is the maximum number of iterations; EPS is the machine 
   precision
*/
{
  const Int MAXITS=200;
  const Doub EPS=numeric_limits<Doub>::epsilon();
  const Doub TOLF=1.0e-8, TOLX=EPS, STPMX=100.0, TOLMIN=1.0e-12;

  Bool restrt, skip;
  Int i, its, j, n=x.size();
  Doub den, f, fold, stpmax, sum, temp, test;
  VecDoub fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n);
  QRdcmp *qr;

  // Set up NRfmin and NRfdjac objects.
  NRfmin<T> fmin(vecfunc);
  NRfdjac<T> fdjac(vecfunc);

  // Make alias to simplify coding.
  VecDoub &fvec=fmin.fvec;

  // The vector fvec is also computed by this call.
  f = fmin(x);
  
  // Test for initial guess being a root.
  // Use more stringent test than simply TOLF.
  test=0.0;
  for (i=0;i<n;i++)
    if (abs(fvec[i]) > test) test=abs(fvec[i]);
  if (test < 0.01*TOLF) {
    check = false;
    return;
  }
	
  // Calculate stpmax for line searches.
  for (sum=0.0,i=0;i<n;i++) sum += SQR(x[i]);
  stpmax = STPMX*MAX(sqrt(sum),Doub(n));

  // Ensure initial Jacobian gets computed
  restrt=true;

  // Start of iteration loop.
  for (its=1;its<=MAXITS;its++) {
    // Initialise, or reinitialise, Jacobian and QR decompose it.
    if (restrt) {
      qr = new QRdcmp(fdjac(x,fvec));
      // Reinitialise the Jacobian with the unit matrix
      if (qr->sing) {
	MatDoub one(n,n,0.0);
	for (i=0;i<n;i++) one[i][i]=1.0;
	delete qr;
	qr = new QRdcmp(one);
      }
    } else {
      // Carry out Broyden update.
      // s = deltax
      for (i=0;i<n;i++) s[i]=x[i]-xold[i];
      // t = R \dot s
      for (i=0;i<n;i++) {
	for (sum=0.0,j=i;j<n;j++) sum += qr->r[i][j]*s[j];
	t[i]=sum;
      }
      skip=true;
      // w = delta F - B \dot s
      for (i=0;i<n;i++) {
	for (sum=0.0,j=0;j<n;j++) sum += qr->qt[j][i]*t[j];
	w[i] = fvec[i]-fvcold[i]-sum;
	// Don't update with noisy components of w
	if (abs(w[i]) >= EPS*(abs(fvec[i])+abs(fvcold[i]))) skip=false;
	else w[i]=0.0;
      }
      if (!skip) {
	// t = Q^T \dot w
	qr-> qtmult(w,t);
	for (den=0.0,i=0;i<n;i++) den += SQR(s[i]);
	// Store s/(s \dot s) in s
	for (i=0;i<n;i++) s[i] /= den;
	// Update R and Q^T
	qr->update(t,s);
	if (qr->sing) throw("singular update in broydn");
      }
    }
    qr->qtmult(fvec,p);
    // RHS for linear equations is -Q^T \dot F.
    for (i=0;i<n;i++) p[i] = -p[i];
    // Compute \Nabla f ~= (Q \dot R)^T \dot F for the line search.
    for (i=n-1;i>=0;i--) {
      for (sum=0.0,j=0;j<=i;j++) sum -= qr->r[j][i]*p[j];
      g[i]=sum;
    }
    // Store x and F.
    for (i=0;i<n;i++) {
      xold[i]=x[i];
      fvcold[i]=fvec[i];
    }
    // Store f.
    fold = f;
    // Solve linear equations.
    qr->rsolve(p,p);
    // Check that the slope is negative.
    Doub slope=0.0;
    for (i=0;i<n;i++) slope += g[i]*p[i];
    if (slope >= 0.0) {
      restrt=true;
      continue;
    }
    // lnsrch returns new x and f. 
    // It also calculates fvec at the new x when it calls fmin.
    lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);

    // Test for convergence on function values.
    test=0.0;
    for (i=0;i<n;i++)
      if (abs(fvec[i]) > test) test=abs(fvec[i]);
    if (test < TOLF) {
      check=false;
      delete qr;
      return;
    }
    if (check) {
      // True if search failed to find new x.
      if (restrt) {
	// Failure; already tried reinitializing the Jacobian.
	delete qr;
	return;
      } else {
	// Check for gradient of f zero, i.e. spurious convergence.
	test=0.0;
	den=MAX(f,0.5*n);
	// The relative gradient is estimated as (deltaf/f)/(deltax/x),
	// where deltaf ~= \Nabla f deltax.
	for (i=0;i<n;i++) {
	  temp=abs(g[i])*MAX(abs(x[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
	if (test < TOLMIN) {
	  delete qr;
	  return;
	}
	// Try reinitialising the Jacobian.
	else restrt=true;
      }
    } else {
      // Successful step; will use Boyden update for next step.
      restrt=false;
      // Test for convergence on \delta x.
      test=0.0;
      for (i=0;i<n;i++) {
	temp=(abs(x[i]-xold[i]))/MAX(abs(x[i]),1.0);
	if (temp > test) test=temp;
      }
      if (test < TOLX) {
	delete qr;
	return;
      }
    }
  }
  throw("MAXITS exceeded in broydn");
}
  
