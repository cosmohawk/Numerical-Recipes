template <class T>
struct F1dim {
  /* Must accompany linmin in Linemethod. */

  const VecDoub &p;
  const VecDoub &xi;
  Int n;
  T &func;
  VecDoub xt;

F1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcc) : p(pp),
    xi(xii), n(pp.size()), func(funcc), xt(n) {}
  /* Constructor takes as inputs as n-dimensional point p[0..n-1] and an n-
     dimensional direction xi[0..n-1] from linmin, as well as the function or 
     functor that takes a vector argument. */

  Doub operator() (const Doub x)
  /* Functor returning value of given function along a one-dimensional line. */
  {
    for (Int j=0;j<n;j++)
      xt[j]=p[j]+x*xi[j];
    return func(xt);
  }
};

template <class T>
struct Linemethod {
  /* Base class for line minimization algorithms. Provides the 
     line-minimization routine linmin. */
  VecDoub p;
  VecDoub xi;
  T &func;
  Int n;

Linemethod(T &funcc) : func(funcc) {}
  /* Constructor argument is the user supplied function or functor to be 
     minimized. */

  Doub lmin()
  /* Line minimization routine. Given an n-dimensional point p[0..n-1] and an 
     n-dimensional direction xi[0..n-1], moves and resets p to where the 
     function or functor func(p) takes on a minimum along the direction xi from
     p, and replaces xi by the actual vector displacement that p was moved. 
     Also returns the value of func at the returned location p. This is 
     actually all accomplished by calling the routines bracket and minimize of
     Brent. */
  { 
    Doub ax, xx, xmin;
    n = p.size();
    F1dim<T> f1dim(p,xi,func);

    // Initial guess for brackets
    ax=0.0; 
    xx=1.0;

    Brent brent;
    brent.bracket(ax,xx,f1dim);
    xmin = brent.minimize(f1dim);
    
    // Construct the vector results to return.
    for (Int j=0;j<n;j++) {
      xi[j] *= xmin;
      p[j] += xi[j];
    }
    return brent.fmin;
  }
};

template <class T>
struct Powell : Linemethod<T> {
  /* Multidimensional minimization by Powell's method. */

  Int iter;
  Doub fret; // Value of the function at the minimum.

  // Variables from a templated base class are not automatically inherited.
  using Linemethod<T>::func;
  using Linemethod<T>::linmin;
  using Linemethod<T>::p;
  using Linemethod<T>::xi;

  const Doub ftol;
  
 Powell(T &func, const Doub ftoll=3.0e-8) : Linemethod<T>(func), ftol(ftoll) {}
  /* Constructor arguments are func, the function or functor to be minimized, 
     and an optional argument ftoll, the fractional tolerance in the function
     value such that failure to decrease by more than this amount on one 
     iteration signals doneness. */

  VecDoub minimize(VecDoub_I &pp)
    /* Minimization of a functor or function of n variables. Input consists of
       an initial starting point pp[0..n-1]. The initial matrix ximat[0..n-1]
       [0..n-1], whose columns contain the inital pp[0..n-1]. The initial 
       matrix ximat whose columns contain the initial set of directions, is
       set to the n unit vectors. Returned is the best point found, at which 
       point fret is the minimum function value and iter is the number of 
       iterations taken. */
  {
    Int n=pp.size();
    MatDoub ximat(n,n,0.0);
    for (Int i=0;i<n;i++) ximat[i][i] = 1.0;
    return minimize(pp,ximat);
  }

  VecDoub minimize(VecDoub_I &pp, MatDoub_IO &ximat)
    /* Alternative interface: Input consists of the initial starting point 
       pp[0..n-1] and an initial matrix ximat[0..n-1][0..n-1], whose columns 
       conatin the initial set of directions. On output ximat is the then-
       current direction set. */
  {
    const Int ITMAX=200; // Maximum allowed iterations.
    const Doub TINY=1.0e-25; // A small number.

    Doub fptt;
    Int n=pp.size();
    p = pp;
    VecDoub pt(n),ptt(n);
    xi.resize(n);
    fret=func(p);

    // Save the initial point.
    for (Int j=0; j<n;j++) pt[j]=p[j];

    for (iter=0;;++iter) {
      Doub fp=fret;
      Int ibig=0;
      Doub del = 0.0; // Will be the biggest function decrease.
      // In each iteration, loop over all directions in the set.
      for (Int i=0;i<n;i++) {
	// Copy the direction,
	for (Int j=0;j<n;j++) xi[j]=ximat[j][i];
	fptt = fret;
	// minimize along it,
	fret = linmin();
	// and record it if it is the largest decrease so far.
	if (fptt-fret > del) {
	  del = fptt-fret;
	  ibig = i+1;
	}
      }
      // The termination criterion:
      if (2.0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) {
	return p;
      }
      if (iter == ITMAX) throw("powell exceeding maximum iterations.");
      // Construct the extrapolated point and the average direction moved.
      // Save the old starting point.
      for (Int j=0;j<n;j++) {
	ptt[j]=2.0*p[j]-pt[j];
	xi[j]=p[j]-pt[j];
	pt[j]=p[j];
      }
      fptt=func(ptt); // Function value at extrapolated point.
      if (fptt < fp) {
	Doub t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
	if (t < 0.0) {
	  /* Move to the minimum of the new direction, and save the new 
	     direction. */
	  fret = linmin();
	  for (Int j=0;j<n;j++) {
	    ximat[j][ibig-1]=ximat[j][n-1];
	    ximat[j][n-1]=xi[j];
	  }
	}
      }
    }
  }
};

template <class T>
struct Frprmn : Linemethod<T> {
  // Multidimensional minimization by the Flethcer-Reeves-Polak-Ribiere method.

  Int iter;
  Doub fret; // Value of the function at the minimum.

  // Variables from a templated base class are not automatically inherited.
  using Linemethod<T>::func;
  using Linemethod<T>::linmin;
  using Linemethod<T>::p;
  using Linemethod<T>::xi;

  const Doub ftol;
  
 Frprmn(T &funcd, const Doub ftol=3.0e-8) : Linemethod<T>(funcd), 
    ftol(ftoll) {}
  /* Constructor arguments are funcd, the function or functor to be minimized, 
     and an optional argument ftoll, the fractional tolerance in the function 
     value such that failure to decrease by more than this amount on one 
     iteration signals doneness. */

  VecDoub minimize(VecDoub_I &pp)
    /* Given a starting point pp[0..n-1], performs the minimization on a 
       function whose value and gradient are provided by a functor funcd. */
  {
    const Int ITMAX=200;
    const Doub EPS=1.03-18;
    const Doub GTOL=1.0e-8;
    /* Here ITMAX is the maximum allowed number of iterations; EPS is a small
       number to rectify the special case of converging to exactly zero 
       function value; and GTOL is the convergence criterion for the zero
       gradient test. */
    Doub gg,dgg;
    Int n=pp.size(); // Initialisations
    p=pp;
    VecDoub g(n),h(n);
    xi.resize();
    Doub fp=func(p);
    func.df(p,xi);
    for (Int j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j] = h[j]=g[j];
    }

    // Loop over iterations.
    for (Int its=0;its<ITMAX;its++) {
      iter=its;
      fret=linmin();
      // One possible return:
      if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS))
	return p;
      fp = fret;
      func.df(p,xi);
      // Test for convergence on zero gradient.
      Doub test=0.0;
      Doub den=MAX(abs(fp),1.0);
      for (Int j=0;j<n;j++) {
	Doub temp=abs(xi[j])*MAX(abs(p[j]),1.0)/den;
	if (temp > test) test=temp;
      }
      // The other possible return:
      if (test < GTOL) return p;
      dgg = gg = 0.0;
      for (Int j=0; j<n; j++) {
	gg += g[j]*g[j];
	// dgg += xi[j]*xi[j]; // This statement for Fletcher-Reeves.
	dgg += (xi[j]+g[j])*xi[j]; // This statement for Polak-Ribiere.
      }
      // If gradient is exactly zero then we are already done, unlikely.
      if (gg == 0.0)
	return p;
      Doub gam=dgg/gg;
      for (Int j=0;j<n;j++) {
	g[j] = -xi[j];
	xi[j] = h[j]=g[j]+gam*h[j];
      }
    }
    throw("Too many iterations in frprmn");
  }
};

template <class T>
struct Df1dim {
  // Must accompany linmin in Dlinemethod.
  const VecDoub &p;
  const VecDoub &xi;
  Int n;
  T &funcd;
  VecDoub xt;
  VecDoub dft;

Df1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcdd) : p(pp),
    xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
  /* Constructo takes as inputs as n-dimensional point p[0..n-1] and an 
     n-dimensional direction xi[0..n-1] from linmin, as well as the functor 
     funcd. */

  Doub operator()(const Doub x)
  // Functor returning value of the given function along a one-dimensional line
  {
    for (Int j=0;j<n;j++)
      xt[j]=p[j]+x*xi[j];
    return funcd(xt);
  }

  Doub df(const Doub x)
  // Returns derivative along the line.
  {
    /* Dbrent always evaluates the derivative at the same value as the function
       so xt is unchanged. */
    Doub df1=0.0;
    funcd.df(xt,dft);
    for (Int j=0;j<n;j++)
      df1 += dft[j]*xi[j];
    return df1;
  }
};

template <class T>
struct Dlinemethod {
  /* Base class for line-minimization using derivative information. Provides 
     the line minimization routine linmin. */
  VecDoub p;
  VecDoub xi;
  T &func;
  Int n;
  
Dlinemethod(T &funcc) : func(funcc) {}
  /* Constructor argument is the user-supplied function or functor to be 
     minimized. */

  Doub linmin()
  /* Lin minimization routine. Given an n-dimensional point p[0..n-1] and an 
     n-dimensional direction xi[0..n-1], moves and resets p to where the 
     function or functor func(p) takes on a minimum along the direction xi from
     p, and replaces xi by the actual vector displacement that p was moved.
     Also returns the value of func at the returned location p. All of this is
     actually accomplished by calling the routines bracket and minimize of
     Dbrent. */
  {
    Doub ax,xx,xmin;
    n=p.size();
    Df1dim<T> df1dim(p,xi,func);
    
    // Initial guess for brackets.
    ax=0.0;
    xx=1.0;

    Dbrent dbrent;
    dbrent.bracket(ax,xx,df1dim);
    xmin = dbrent.minimize(df1dim);

    // Construct the vector results to return.
    for (Int j=0; j<n; j++) {
      xi[j] *= xmin;
      p[j] += xi[j];
    }

    return dbrent.fmin;
  }
};


    
