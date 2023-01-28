template <class T>
Bool zbrac(T &func, Doub &x1, Doub &x2)
/* Given a function or functor func and an initial guessed range x1 to x2, the
   routine expands the range geometrically until a root is bracketed by the 
   returned values x1 and x2 (in which case zbrac returns true) or until the
   range becomes unacceptably large (in which case zbrac returns false). */
{
  const Int NTRY=50;
  const Doub FACTOR=1.6;

  if (x1 == x2) throw("Bad initial range in zbrac");
  Doub f1 = func(x1);
  Doub f2 = func(x2);
  for (Int j=0; j<NTRY; j++) {
    if (f1*f2 < 0.0) return true;
    if (abs(f1) < abs(f2))
      f1 = func(x1 += FACTOR*(x1-x2));
    else
      f2 = func(x2 += FACTOR*(x2-x1));
  }
  return false;
}

template <class T>
void zbrak(T &fx, const Doub x1, const Doub x2, const Int n, VecDoub_O &xb1,
	   VecDoub_O &xb2, Int &nroot)
/* Given a function or functor fx defined on the interval [x1,x2], subdivide
   the interval into n equally spaced segments, and search for zero crossings 
   of the function. nroot will be set to the number of bracketing pairs found.
   If it is positive, the arrays xb1[0..nroot-1] and xb2[0..nroot-1] will be 
   filled sequentially with any bracketing pairs that are found. On input, 
   these vectors may have any size, including zero; they will be resized to >= 
   nroot. */
{
  Int nb = 20;
  
  xb1.resize(nb);
  xb2.resize(nb);
  nroot = 0;
  
  // Determine the spacing appropriate to the mesh.
  Doub dx = (x2-x1)/n;
  Doub x = x1;
  Doub fp = fx(x1);
  
  // Loop over all intervals.
  for (Int i=0; i<n; i++) {
    Doub fc = fx(x += dx);
    // If a sign change occurs then record values for the bounds.
    if (fc*fp <= 0.0) {
      xb1[nroot] = x-dx;
      xb2[nroot++] = x;
      if (nroot == nb) {
	VecDoub tempvec1(xb1), tempvec2(xb2);
	xb1.resize(2*nb);
	xb2.resize(2*nb);
	for (Int j=0; j<nb; j++) {
	  xb1[j] = tempvec1[j];
	  xb2[j] = tempvec2[j];
	}
	nb *= 2;
      }
    }
    fp = fc;
  }
}

template <class T>
Doub rtbis(T &func, const Doub x1, const Doub x2, const Doub xacc)
  /* Using bisection, return the root of a function or functor func known to 
     lie between x1 and x2. The root will be refined until its accuracy is 
     +/- xacc. */
{
  // Maximum number of bisections.
  const Int JMAX=50;

  Doub dx, xmid, rtb;
  Doub f = func(x1);
  Doub fmid = func(x2);
  
  if (f*fmid >= 0.0) throw("Root must be bracketed for bisection in rtbis");

  // Orient search so that f > 0 lies at x + dx.
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  
  // Bisection loop.
  for (Int j=0; j<JMAX; j++) {
    fmid = func(xmid=rtb+(dx *= 0.5));
    if (fmid <= 0.00) rtb = xmid;
    if (abs(dx) < xacc || fmid == 0.0) return rtb;
  }

  throw("Too many bisections in rtbis");
}
  
template <class T>
Doub rtflsp(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  /* Using the false-position method, return the root of a function or functor
     func known to lie between x1 and x2. The root is refined until its
     accuracy is +/- xacc. */

  // Set the maximum allowed number of iterations
  const Int MAXIT=50;

  Doub xl, xh, del;
  
  // Be sure that the interval brackets a root
  Doub fl=func(x1);
  Doub fh=func(x2);
  if (fl*fh > 0.0) throw("Root must be bracketed in rtflsp");
  
  // Identify the limits so that xl corresponds to the low side.
  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    xl = x2;
    xh = x1;
    SWAP(fl,fh);
  }
  Doub dx = xh-xl;

  // False position loop.
  for (Int j=0; j<MAXIT; j++) {
    
    // Increment w.r.t. latest value
    Doub rtf = xl+dx*fl/(fl-fh);
    Doub f = func(rtf);

    // Replace appropriate limit.
    if (f < 0.0) {
      del = xl-rtf;
      xl = rtf;
      fl = f;
    } else {
      del = xh -rtf;
      xh = rtf;
      fh=f;
    }
    dx = xh -xl;

    // Convergence.
    if (abs(del) < xacc || f == 0.0) return rtf;
  }
  throw("Maximum number of iterations exceeded in rtfflsp");
}

template <class T>
Doub rtsec(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  /* Using the secant method, return the root of a function or functor
     func known to lie between x1 and x2. The root is refined until its
     accuracy is +/- xacc. */

  // Set the maximum allowed number of iterations
  const Int MAXIT=50;

  Doub xl, rts;

  Doub fl=func(x1);
  Doub f=func(x2);

  // Pick the bound with the smaller function as the most recent guess.
  if (abs(fl) < abs(f)) {
    rts=x1;
    xl=x2;
    SWAP(fl,f);
  } else {
    xl=x1;
    rts=x2;
  }

  // Secant loop
  for (Int j=0;j<MAXIT;j++) {
    
    // Increment w.r.t. latest value.
    Doub dx = (xl-rts)*f/(f-fl);
    xl=rts;
    fl=f;
    rts += dx;
    f = func(rts);

    // Convergence.
    if (abs(dx) < xacc || f == 0.0) return rts;
  }
  throw("Maximum number of iterations exceeded in rtsec");
}

template <class T>
Doub zriddr(T &func,const Doub x1, const Doub x2, const Doub xacc) {
  /* Using Ridder's method, return the root of a function or functor func known
     to lie between x1 and x2. The root will be refined to an approximate
     accuracy xacc. */

  // Set the maximum allowed number of iterations
  const Int MAXIT=60;

  Doub fl=func(x1);
  Doub fh=func(x2);

  if ((fl > 0.0 && fh < 0.0) || (fl <0.0 && fh > 0.0)) {
    Doub xl = x1;
    Doub xh = x2;

    // A highly unlikely value to simplify the logic below.
    Doub ans=-9.99e99;

    for (Int j=0;j<MAXIT;j++) {
      Doub xm = 0.5*(xl+xh);
      // First of two function evaluations per iteration.
      Doub fm = func(xm);
      Doub s = sqrt(fm*fm-fl*fh);
      if (s == 0.0) return ans;
      // Updating formula.
      Doub xnew = xm +(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (abs(xnew-ans) <= xacc) return ans;
      ans = xnew;
      // Second of two function evaluations per iteration.
      Doub fnew = func(ans);
      if (fnew == 0.0) return ans;
      // Book keeping to keep the root bracketed on next iteration.
      if (SIGN(fm,fnew) != fm) {
	xl = xm;
	fl = fm;
	xh = ans;
	fh = fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh = ans;
	fh = fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl = ans;
	fl = fnew;
      } else throw("never get here.");
      if (abs(xh-xl) <= xacc) return ans;
    }
    throw("zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    throw("root must be bracketed in zriddr.");
  }
}

template <class T>
Doub zbrent(T &func,const Doub x1, const Doub x2, const Doub tol) {
  /* Using Brent's method, return the root of a function or functor func known
     to lie between x1 and x2. The root will be refined to an approximate
     accuracy tol. */
  // Maximum number of allowed iterations
  const Int ITMAX=100;
  
  // Machine floating point precision
  const Doub EPS=numeric_limits<Doub>::epsilon();

  Doub a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb >0.0) || (fb < 0.0 && fc < 0.0))
    throw("root must be bracketed in zbrent");
  fc = fb;

  for (Int iter=0; iter<ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      // Rename a,b,c and adjust bounding interval d.
      c = a;
      fc = fa;
      e = d = b-1;
    }
    if (abs(fc) < abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    
    // Convergence check.
    tol1 = 2.0*EPS*abs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (abs(xm) <= tol1 || fb == 0.0) return b;
    
    // Attempt to inverse quadratic interpolation.
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
      s = fb/fa;
      if (a == c){
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      // Check whether in bounds.
      if (p >0.0) q=-q;
      p=abs(p);
      Doub min1=3.0*xm*q-abs(tol1*q);
      Doub min2=abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	// Accept interpolation
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      // Bounds decreasing too slowly, use bisection.
      d=xm;
      e=d;
    }
    // Move last best guess to a.
    a=b;
    fa=fb;
    // Evaluate new trial root.
    if (abs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=func(b);
  }
  throw("Maximum number of iterations exceeded in zbrent");
}

template <class T>
Doub rtnewt(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
  /* Using the Newton-Raphson method, return the root of a function or functor
     func known to lie between x1 and x2. The root will be refined to an 
     approximate accuracy +/- xacc. funcd is a user supplied struct that 
     returns the function value as a functor and the first derivative of the
     function at the point x as the function df. */

  // Set meximum number of iterations.
  const Int JMAX=20;

  // Initial guess.
  Doub rtn = 0.5*(x1+x2);

  for (Int j=0; j<JMAX; j++) {
    Doub f = funcd(rtn);
    Doub df = funcd.df(rtn);
    Doub dx=f/df;
    rtn -= dx;
    if ((x1-rtn)*(rtn-x2) < 0.0)
      throw("Jumped out of brackets in rtnewt");
    // Check convergence.
    if (abs(dx) < xacc) return rtn;
  }
  throw("Maximum number of iterations exceeded in rtnewt.");
}

template <class T>
Doub rtsafe(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
  /* Using a combination of Newton-Raphson and bisection, return the root of a
     function or functor func known to lie between x1 and x2. The root will be
     refined to an approximate accuracy +/- xacc. funcd is a user supplied 
     struct that returns the function value as a functor and the first 
     derivative of the function at the point x as the function df. */

  // Maximum allowed number of iterations.
  const Int MAXIT=100;

  Doub xh, xl;
  Doub fl = funcd(x1);
  Doub fh = funcd(x2);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    throw("Root must be bracketed in rtsafe");
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  // Orient search so that f(xl) < 0.
  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    xh = x1;
    xl = x2;
  }
  
  // Initialise the guess for root, the penultimate step size and the last step
  Doub rts = 0.5*(x1+x2);
  Doub dxold = abs(x2-x1);
  Doub dx = dxold;
  Doub f = funcd(rts);
  Doub df = funcd.df(rts);
  
  // Loop over allowed iterations
  for (Int j=0; j<MAXIT; j++) {
    
    // Bisect if Newton out of range or not decreasing fast enough. 
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
	|| (abs(2.0*f) > abs(dxold*df))) {
      dxold = dx;
      dx = 0.5*(xh-xl);
      rts = xl+dx;
      // Convergence check
      if (xl == rts) return rts;
    } else {
      // Newton step acceptable, take it.
      dxold = dx;
      dx=f/df;
      Doub temp=rts;
      rts -= dx;
      // Convergence check.
      if (temp == rts) return rts;
    }

    // Convergence criterion.
    if (abs(dx) < xacc) return rts;
    
    // The one new function evaluation per iteration.
    f = funcd(rts);
    df = funcd.df(rts);

    // Maintain the bracket on the root.
    if (f < 0.0)
      xl = rts;
    else 
      xh = rts;
  }
  throw("Maximum number of iterations exceeded in rtsafe");
}


    
    
