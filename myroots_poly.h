void laguer(VecComplex_I &a, Complex &x, Int &its) {
  /* Given the m+1 complex coefficients a[0..m] of the polynomial
     \sum_{i=0}^m a[i]x^i, and given a complex value x, this routine improves
     x by Laguer's methd until it converges, within the achievable roundoff
     limit, to a root of the given polynomial. The number of iterations taken
     is returned as its. */

  const Int MR=8,MT=10,MAXIT=MT*MR;
  const Doub EPS=numeric_limits<Doub>::epsilon();
  
  /* Here, EPS is the estimated fractional roundoff error. We try to break 
     (rare) limit cycles with MR different fractional values, once every MT 
     steps, for MAXIT total allowed iterations. */
  
  // Fractions used to break a limit cycle.
  static const Doub frac[MR+1] = 
    {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
  
  Complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
  Int m=a.size()-1;
  
  // Loop over iterations up to allowed maximum.
  for (Int iter=1;iter<=MAXIT;iter++) {
    its = iter;
    b=a[m];
    Doub err=abs(b);
    d=f=0.0;
    Doub abx=abs(x);
    
    // Efficient computation of the polynomial and its first two derivatives.
    for (Int j=m-1;j>=0;j--){
      // f stores P''/2
      f = x*f+d;
      d = x*d+b;
      b = x*b+a[j];
      err = abs(b)+abx*err;
    }

    // Estimate the roundoff error in evaluating polynomial
    err *= EPS;

    // Check for convergence on root
    if (abs(b) <= err) return;

    // The generic case: Use Laguerre's formula.
    g = d/b;
    g2 = g*g;
    h = g2 - 2.0*f/b;
    sq = sqrt(Doub(m-1)*(Doub(m)*h-g2));
    gp = g+sq;
    gm = g-sq;
    Doub abp = abs(gp);
    Doub abm = abs(gm);
    if (abp < abm) gp=gm;
    dx = MAX(abp,abm) > 0.0 ? Doub(m)/gp : polar(1+abx,Doub(iter));
    x1 = x-dx;
    // Convergence check
    if (x == x1) return;
    if (iter % MT != 0) x = x1;
    // Every so often we take a fractional step to break any limit cycle (rare)
    else x -= frac[iter/MT]*dx;
  }
  throw("too many iterations in laguer");
  // Very unusual; can occur only for complex roots. Try a different start.
}

void zroots(VecComplex_I &a, VecComplex_O &roots, const Bool &polish) {
  /* Given the m+1 complex coefficients a[0..m] of the polynomial 
     \sum_{i=0}^m a[i]x^i, this routine successively calls laguer and finds all
     m complex roots in roots[0..m-1]. The boolean variable polish should be
     input as true if polishing (also by Laguerre's method) is desired, false 
     if the roots wil be subsequently polished by other means. */

  // A small number
  const Doub EPS=1.0e-14;

  Int i,its;
  Complex x,b,c;
  Int m=a.size()-1;
  VecComplex ad(m+1);

  // Copy of coefficients for successive deflation
  for (Int j=0; j<=m; j++) ad[j]=a[j];

  // Loop over each root to be found.
  for (Int j=m-1;j>=0;j--) {

    // Start at zero to favour convergence to smallest remaining root.
    x = 0.0;
    VecComplex ad_v(j+2);
    for (Int jj=0;jj<j+2;jj++) ad_v[jj]=ad[jj];
    laguer(ad_v,x,its);
    if (abs(imag(x)) <= 2.0*EPS*abs(real(x)))
      x = Complex(real(x),0.0);
    roots[j]=x;
    b=ad[j+1];
    for (Int jj=j;jj>=0;jj--) {
      c = ad[jj];
      ad[jj] = b;
      b = x*b+c;
    }
  }

  // Polish the roots using the undeflated coefficients.
  if (polish) 
    for (Int j=0; j<m; j++) 
      laguer(a,roots[j],its);

  // Sort roots by their real parts by straight insertion.
  for (Int j=1;j<m;j++) {
    x=roots[j];
    for (i=j-1;i>=0;i--) {
      if (real(roots[i]) <= real(x)) break;
      roots[i+1]=roots[i];
    }
    roots[i+1]=x;
  }
}

