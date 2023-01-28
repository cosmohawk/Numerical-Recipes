void qroot(VecDoub_I &p, Doub &b, Doub &c, const Doub eps)
/* Given n+1 coefficients p[0..n] of a polynomial of degree n, and trial values
   for the coefficients of a quadratic factor bx^2 + cx, improve the solution
   until the coefficients b, c change by less than eps. The routine poldiv from
   Sec. 5.1 is used. */
{
  // Maximum number of iterations.
  const Int ITMAX=20;
  const Doub TINY=1.0e-14;
  Doub sc, sb, s, rc, rb, r, dv, delc, delb;
  Int n=p.size() -1;
  VecDoub d(3),q(n+1),qq(n+1),rem(n+1);
  d[2]=1.0;
  
  for (Int iter=0;iter<ITMAX;iter++) {
    d[1]=b;
    d[0]=c;
  
    // First division, r, s
    poldiv(p,d,q,rem);
    s=rem[0];
    r=rem[1];
    
    // Second division, parital r, s with respect to c
    poldiv(q,d,qq,rem);
    sb = -c*(rc = -rem[1]);
    rb = -b*rc+(sc = -rem[0]);

    // Solve 2x2 equation
    dv = 1.0/(sb*rc - sc*rb);
    delb = (r*sc - s*rc)*dv;
    delc = (-r*sb+s*rb)*dv;

    b += (delb=(r*sc-s*rc)*dv);
    c += (delc=(-r*sb+s*rb)*dv);

    if ((abs(delb) <= eps*abs(b) || abs(b) < TINY)
	&& (abs(delc) <= eps*abs(c) || abs(c) < TINY)) {
      // Coefficients converged.
      return;
    }
  }
  throw("Too many iterations in routine qroot");
}


