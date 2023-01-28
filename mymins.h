struct Bracketmethod
/* Base class for one-dimensional minimisation routines. Provides a routine
   to bracket a minimum and several utility functions. */
{
  Doub ax, bx, cx, fa, fb, fc;
  
  template <class T>
  void bracket(const Doub a, const Doub b, T &func)
  /* Given a function or functor func, and given distinct inital points ax and
     bx, this routine searches in the downhill direction (defined by the 
     function evaluated at the initial points) and returns new points ax, bx, 
     cx that bracket the minimum of the function. Also returned are the 
     function values at the three points, fa, fb ,fc. */
  { 
    const Doub GOLD=1.618034, GLIMIT=100.0, TINY=1.0e-20;
    /* Here GOLD is the default ratio by which successive intervals are 
       magnified and GLIMIT is the maximum magnification allowed for a
       parabolic fit step. */
    ax=a; bx=b;
    Doub fu;
    fa = func(ax);
    fb = func(bx);

    // Switch roles of a and b so that we can go downhill in the direction a b
    if (fb > fa) {
      SWAP(ax,bx);
      SWAP(fb,fa);
    }

    // First guess for c
    cx =bx + GOLD*(bx-ax);
    fc = func(cx);

    // Keep returning here until we bracket.
    while (fb > fc) {
      // Compute u by parabolic extrapolation from a,b,c.
      // TINY is used to prevent any possible division  by zero.
      Doub r=(bx-ax)*(fb-fc);
      Doub q=(bx-cx)*(fb-fa);
      Doub u=bx-((bx-cx)*q-(bx-ax)*r)/
	(2.0*SIGN(MAX(abs(q-r),TINY),q-r));
      // We wont go further than ulim.
      Doub ulim=bx+GLIMIT*(cx-bx);
      // Test various possibilities.
      // Is parabolic u between b and c?
      if ((bx-u)*(u-cx) > 0.0) {
	fu = func(u);
	if (fu < fc) {
	  // Got a minimum between b and c.
	  ax=bx;
	  bx=u;
	  fa=fb;
	  fb=fu;
	  return;
	} else if (fu > fb) {
	  // Got a minimum between a and u.
	  cx = u;
	  fc = fu;
	  return;
	}
	// Parabolic fit was no use, use default magnification.
	u = cx+GOLD*(cx-bx);
	fu = func(u);
      } else if ((cx-u)*(u-ulim) > 0.0) {
	// Parabolic fit is between c and its allowed limit.
	fu = func(u);
	if (fu < fc) {
	  shft3(bx,cx,u,u+GOLD*(u-cx));
	  shft3(fb,fc,fu,func(u));
	}
      } else if ((u-ulim)*(ulim-cx) >= 0.0) {
	// Limit parabolic u to maximum allowed value.
	u = ulim;
	fu = func(u);
      } else {
	// Reject parabolic u. Use default magnification.
	u = cx+GOLD*(cx-bx);
	fu = func(u);
      }
      // Eliminate oldest point and continue.
      shft3(ax,bx,cx,u);
      shft3(fa,fb,fc,fu);
    }
  }

  // Utility unctions used in this struture or others derived from it.
  inline void shft2(Doub &a, Doub &b, const Doub c) {
    a=b; b=c;
  }

  inline void shft3(Doub &a, Doub &b, Doub&c, const Doub d) {
    a=b; b=c; c=d;
  }

  inline void mov3(Doub &a, Doub &b, Doub&c, const Doub d, 
		   const Doub e, const Doub f) {
    a=d; b=e; c=f;
  }
};
	
struct Golden : Bracketmethod {
  /* Golden section search for minimum. */
  Doub xmin, fmin;
  const Doub tol;
 
 Golden(const Doub toll=3.0e-8) : tol(toll) {}
  template <class T>
    Doub minimize(T &func)
    /* Given a function or functor f, and given a bracketing triplet of 
       abscissas ax, bx, cx (such that bx is between ax and cx and that f(bx)
       is less than both f(ax) and f(cx)), this routine performs a golden 
       section search for the minimum, isolating it to a fractional precision
       of about tol. The abscissa of the minimum is returned as xmin, and the 
       function value at the minimum is returned as min, the returned function
       value. */
    {
      const Doub R=0.61803399, C=1.0-R; // The golden ratios.
      
      // At any given time we will keep track of four points, x0,x1,x2,x3.
      Doub x1,x2;
      Doub x0=ax;
      Doub x3=cx;
      
      // Make x0 to x1 the smaller segement, and fill the new point to be tried
      if (abs(cx-bx) > abs(bx-ax)) {
	x1=bx;
	x2=bx+C*(cx-bx);
      } else {
	x2=bx;
	x1=bx-C*(bx-ax);
      }

      /* The initial function evaluations. Note that we never need to evaluate
	 the function at the original end points. */
      Doub f1=func(x1);
      Doub f2=func(x2);
      while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
	// One possible outcome, its housekeeping and a function evaluation.
	if (f2<f1) {
	  shft3(x0,x1,x2,R*x2+C*x3);
	  shft2(f1,f2,func(x2));
	} else {
	  // The other outcome and its function evaluation.
	  shft3(x3,x2,x1,R*x1+C*x0);
	  shft2(f2,f1,func(x1));
	}
      }

      // Back to see if we are done.
      if (f1 < f2) {
	// We are done, output the best of the two current values.
	xmin=x1;
	fmin=f1;
      } else {
	xmin = x2;
	fmin = f2;
      }
      return xmin;
    }
};

struct Brent : Bracketmethod {
  // Brent's method to find a minimum.
  Doub xmin, fmin;
  const Doub tol;
 
 Brent(const Doub toll=3.0e-8) : tol(toll) {}
  template <class T>
    Doub minimize(T &func)
  /* Given a function or functor f, and given a bracketing triplet of       
     abscissas ax, bx, cx (such that bx is between ax and cx and that f(bx) is 
     less than both f(ax) and f(cx)), this routine isolates the minimum to a 
     fractional precision of about tol, using Brent's method. The abscissa of 
     the minimum is returned as xmin, and thefunction value at the minimum is 
     returned as min, the returned function value. */
  {
    const Int ITMAX=100;
    const Doub CGOLD=0.3819660;
    const Doub ZEPS=numeric_limits<Doub>::epsilon()*1.0e-3;
    /* Here ITMAX is the maximum allowed number of iterations; CGOLD is the 
       golden ratio; and ZEPS is a small number that protects against trying to
       acheive fractional accuracy for a minimum that happens to be exactly 
       zero. */
    Doub a,b,d=0.0,etemp,fu,fv,fw,fx;
    Doub p,q,r,tol1,tol2,u,v,w,x,xm;

    // The distance moved on the penultimate step.
    Doub e=0.0;

    // a and b must be in ascending order, but input ascissas need not be.
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    
    // Initialisations.
    x=w=v=bx;
    fw=fv=fx=func(x);
    
    // Main program loop.
    for (Int iter=0;iter<ITMAX;iter++) {
      xm = 0.5*(a+b);
      tol2=2.0*(tol1=tol*abs(x)+ZEPS);
      // Test for done here.
      if (abs(x-xm) <= (tol2-0.5*(b-a))) {
	fmin=fx;
	return xmin=x;
      }
      // Construct a trial parabolic fit.
      if (abs(e) > tol1) {
	r=(x-w)*(fx-fv);
	q=(x-v)*(fx-fw);
	p=(x-v)*q-(x-w)*r;
	q=2.0*(q-r);
	if (q > 0.0) p = -p;
	q=abs(q);
	etemp=e;
	e=d;

	// Conditions that determine the acceptability of the parabolic fit. 
	if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	  // Take the golden section step into the larger of the two segments.
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	else {
	  // Take the parabolic step.
	  d = p/q;
	  u = x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=SIGN(tol1,xm-x);
	}
      } else {
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u = (abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      // This is the only function evaluation per iteration.
      fu = func(u);
      // Now decide what to do with our function evaluation.
      if (fu <= fx) {
	if (u >= x) a=x; else b=x;
	// Housekeeping follows
	shft3(v,w,x,u);
	shft3(fv,fw,fx,fu);
      } else {
	if (u < x) a=u; else b=u;
	if (fu <= fw || w ==x) {
	  v=w;
	  w=u;
	  fv=fw;
	  fw=fu;
	} else if (fu <= fv || v == x || v == w) {
	  v = u;
	  fv = fu;
	}
      }
      // Done with housekeeping, back for another iteration.
    }
    throw("Too many iterations in brent");
  }
};

struct Dbrent : Bracketmethod {
  /* Brent's method to find a minimum, modified to use derivative. */
  Doub xmin,fmin;
  const Doub tol;
 Dbrent(const Doub toll=3.0e-8) : tol(toll) {}
  template <class T>
    Doub minimize(T & funcd)
    /* Given a function or functor f, and given a bracketing triplet of       
    abscissas ax, bx, cx (such that bx is between ax and cx and that f(bx) is 
    less than both f(ax) and f(cx)), this routine isolates the minimum to a 
    fractional precision of about tol, using a modification of Brent's method 
    that include derivatives. The abscissa of the minimum is returned as xmin,
    and thefunction value at the minimum is returned as min, the returned 
    function value. */

    // Comments point out differences with Brent. Read Brent first.
    {
      const Int ITMAX=100;
      const Doub ZEPS=numeric_limits<Doub>::epsilon()*1.0e-3;
      
      // Flags for whether proposed step is acceptable or not.
      Bool ok1, ok2;

      Doub a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
      Doub fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
      
      a=(ax < cx ? ax : cx);
      b=(ax > cx ? ax : cx);
      x=w=v=bx;
      fw=fv=fx=funcd(x);
      dw=dv=dx=funcd.df(x);
      /* All housekeeping chores are doubles due to the necessity of moving 
	 around derivatives as well as function values. */
      for (Int iter=0;iter<ITMAX;iter++) {
	xm = 0.5*(a+b);
	tol1=tol*abs(x)+ZEPS;
	tol2=2.0*tol1;
	if (abs(x-xm) <= (tol2-0.5*(b-a))) {
	  fmin=fx;
	  return xmin=x;
	}
	if (abs(e) > tol1) {
	  // Initialise these ds to an out of bracket value.
	  d1=2.0*(b-a);
	  d2=d1;
	  // Secant methos with one point.
	  if (dw != dx) d1=(w-x)*dx/(dw-dw);
	  if (dv != dx) d2=(v-x)*dx/(dx-dv);
	  // Which of these two estimates of d shall we take?
	  // We must insist that they be within the bracket,
	  // and on the side pointed to by the derivative at x.
	  u1 = x+d1;
	  u2=x+d2;
	  ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	  ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
	  // Movement on the penultimate step.
	  olde=e;
	  e=d;
	  // Take only an acceptable d, if both are acceptable take smallest.
	  if (ok1 || ok2) {
	    if (ok1 && ok2)
	      d = (abs(d1) < abs(d2) ? d1 : d2);
	    else if (ok1)
	      d=d1;
	    else
	      d = d2;
	    if (abs(d) <= abs(0.5*olde)) {
	      u = x+ d;
	      if (u-a < tol2 || b-u < tol2)
		d = SIGN(tol1,xm-x);
	    } else {
	      // Bisect, not golden section.
	      d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
	      // Decide which segment by the sign of the derivative.
	    }
	  } else {
	    d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
	  }
	} else { 
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
	if (abs(d) >= tol1) {
	  u = x+d;
	  fu = funcd(u);
	} else {
	  u = x+SIGN(tol1,d);
	  fu = funcd(u);
	  // If the minimum step in the downhill direction takes us uphill,
	  // then we are done.
	  if (fu > fx) {
	    fmin=fx;
	    return xmin=x;
	  }
	}
	// Housekwwping.
	du=funcd.df(u);
	if (fu <= fx) {
	  if (u >= x) a=x; else b=x;
	  mov3(v,fv,dv,w,fw,dw);
	  mov3(w,fw,dw,x,fx,dx);
	  mov3(x,fx,dx,u,fu,du);
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w==x) {
	    mov3(v,fv,dv,w,fw,dw);
	    mov3(w,fw,dw,u,fu,du);
	  } else if (fu < fv || v == x || v == w) {
	    mov3(v,fv,dv,u,fu,du);
	  }
	}
      }
      throw("Too many iterations in routine dbrent");
    }
};    


  
