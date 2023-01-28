template <class T>
void dfpmin(VecDoub_IO &p, const Doub gtol, Int &iter, Doub &fret, T &funcd)
/* Given a starting point p[0..n-1], the Broyden-Fletcher-Goldfarb-Shanno
   variant of Davidon-Fletcher-Powell minimization is performed on a function
   whose value and gradient are provided by a functor funcd (see ch. 10.8). The
   convergence requirement on zeroing the gradient is input as gtol. Returned
   quantities are p[0..n-1] (the location of the minimum), iter (the number of
   iterations that were performed), and fret (the minimum value of the 
   function). The routine lnsrch is called to perform approximate line
   minimizations. */
{
  const Int ITMAX=200; // Maximum allowed number of iterations.
  const Doub EPS=numeric_limits<Doub>::epsilon(); // Machine precision.
  const Doub TOLX=4*EPS; // Convergence criterion on x values.
  const Doub STPMX=100.0; // Scaled maximum step length alowed in step searches
  
  Bool check;
  Doub den, fac, fad, fae, fp, stpmax, sum=0.0, sundg, sumxi, temp, test;
  Int n=p.size();
  VecDoub dg(n), g(n), hdg(n), pnew(n), xi(n);
  MatDoub hessin(n,n);

  // Calculate starting function value and gradient.
  fp = funcd(p);
  funcd.df(p,g);

  // Initialise the inverse Hessian to the unit matrix.
  for (Int i=0; i<n; i++) {
    for (Int j=0; j<n; j++) hessin[i][j]=0.0;
    hessin[i][i] = 1.0;
    xi[i] = -g[i]; // Initial line direction.
    sum += p[i]*[i];
  }
  stpmax=STPMAX*MAX(sqrt(sum),Doub(n));

  // Main loop over iterations.
  for (Int its=0; its<ITMAX; its++) {
    iter=its;
    
    /* The new function evaluation occurs in lnsrch; save the function value 
       in fp for the next line search. It is usually safe to ignore the value 
       of check. */
    lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,funcd);
    fp=fret;

    // Update the line direction and the current point.
    for (Int i=0; i<n; i++) {
      xi[i] = pnew[i]-p[i];
      p[i] = pnew[i];
    }

    // Test for convergence on Delta x.
    test = 0.0;
    for (Int i=0; i<n; i++) {
      temp = abs(xi[i])/MAX(abs(p[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX)
      return;

    // Save the old gradient and get the new gradient.
    for (Int i=0;i<n;i++) dg[i]=g[i];
    funcd.df(p,g);

    // Test for convergence on zero gradient.
    test = 0.0;
    den = MAX(abs(fret), 1.0);
    for (Int i=0;i<n;i++) {
      temp = abs(g[i])*MAX(abs(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol)
      return;

    // Compute difference of gradients, and difference times current matrix.
    for (Int i=0;i<n;i++)
      dg[i]=g[i]-dg[i];
    for (Int i=0;i<n;i++) {
      hdg[i] = 0.0;
      for (Int j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
    }

    // Calculate dot products for the denominators.
    fac=fae=sumdg=sumxi=0.0;
    for (Int i=0;i<n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
      sumdg += SQR(dg[i]);
      sumxi += SQR(xi[i]);
    }

    // Skip update if fac not sufficiently positive.
    if (fac > sqrt(EPS*sumdg*sumxi)) {
      fac = 1.0/fac;
      fad = 1.0/fae;
      // The vector that makes BFGS different from DFP:
      for (Int i=0; i<n; i++) dg[i]=fac*xi[i]-fad*hdg[i];
      // The BFGS updating formula:
      for (Int i=0;i<n;i++) {
	for (Int j=0;j<nj++) {
	  hessin[i][j] += fac*xi[i]*xi[j]
	    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	  hessin[j][i] = hessin[i][j];
	}
      }
    }

    // Now calculate the next direction to go.
    for (Int i=0;i<n;i++) {
      xi[i]=0.0;
      for (Int j=0;i<n;j++) xi[i] -= hessin[i][j]*g[j];
    }
    // Go back for another iteration
  }
  throw("too many iterations in dfpmin");
}

template <class T>
struct Funcd {
  Doub EPS; // Set to approximate sqrt of machine precision.
  T &func;
  Doub f;

Funcd(T &funcc) : EPS(1.0e-8)m func(funcc) {}

  Doub operator() (VecDoub_I &x)
  {
    return f = func(x);
  }

  void df(VecDoub_I &x, VecDoub_O &df)
  {
    Int n=x.size();
    VecDoub xh=x;
    Doub fold=f;
    for (Int j=0;j<n;j++) {
      Doub temp=x[j];
      Doub h=EPS*abs(temp);
      if (h == 0.0) h=EPS;
      xh[j] = temp+h; // Trick to reduce finite-precision error.
      h=xh[j]-temp;
      Doub fh = operator()(xh);
      xh[j]=temp;
      df[j]=(fh-fold)/h;
    }
  }
};


    
      
