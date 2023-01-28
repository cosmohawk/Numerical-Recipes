struct MCintegrate {
  /* Object for Monte Carlo integration of one or more functions in an 
     ndim-dimensional region. */

  // Number of dimensions, functions, and points sampled.
  Int ndim, nfun, n;

  // Answers: The integrals and their standard error.
  VecDoub ff, fferr;
  VecDoub xlo, xhi, x, xx, fn, sf, sferr;

  // Volume of the box V.
  Doub vol;

  // Pointers to the user supplied functions.
  VecDoub (*funcsp)(const VecDoub &);
  VecDoub (*xmapp)(const VecDoub &);
  Bool (*inregionp)(const VecDoub &);
  
  // Random number generator.
  Ran ran;

  /* Constructor. The arguments are in the order described in the itemized
     list above. */
  MCintegrate(const VecDoub &xlow, const VecDoub &xhigh,
	      VecDoub funcs(const VecDoub &), Bool inregion(const VecDoub &),
	      VecDoub xmap(const VecDoub &), Int randseed);

  // Sample an additional nstep points, accumulating the various sums.
  void step(Int nstep);

  // Calculate answers ff and fferr using the current sums.
  void calcanswers();
};

/* To use the MCintegrate object, we first write functions that describe the
   integrands and the regions of integration W inside the box V. */

VecDoub torusfuncs(const VecDoub &x) {
  // Return the integrands in eq 7.7.5 with rho=1.

  Doub den = 1.;
  VecDoub f(4);
  f[0] = den;
  for (Int i=1; i<4; i++) f[i] = x[i-1]*den;
  return f;
}

Bool torusregion(const VecDoub &x) {
  // Return the inequality 7.7.3.
  return SQR(x[2])+SQR(sqrt(SQR(x[0])+SQR(x[1]))-3.) <= 1;
}

VecDoub torusmap(const VecDoub &s) {
  /* Return the mapping from s to z defined by the last equation in (7.7.7)
     mapping the other coordinates by the identity map. */
  VecDoub xx(s);
  xx[2] = 0.2*log(5.*s[2]);
  return xx;
}

MCintegrate::MCintegrate(const VecDoub &xlow, const VecDoub &xhigh,
			 VecDoub funcs(const VecDoub &), 
			 Bool inregion(const VecDoub &),
			 VecDoub xmap(const VecDoub &), Int ranseed) : 
ndim(xlow.size()), n(0), xlo(xlow), xhi(xhigh), x(ndim), xx(ndim), 
  funcsp(funcs), xmapp(xmap), inregionp(inregion), vol(1.), ran(ranseed) {
  /* Constructor stores pointers to the user functions. */ 

  // Determine size of returned vector.
  if (xmapp) nfun = funcs(xmapp(xlo)).size();
  else nfun = funcs(xlo).size();

  // Size the sum and answer vectors accordingly.
  ff.resize(nfun);
  fferr.resize(nfun);
  fn.resize(nfun);
  sf.assign(nfun,0.);
  sferr.assign(nfun,0.);
  for (Int j=0;j<ndim;j++) vol *= abs(xhi[j]-xlo[j]);
}

// Eq 7.7.1
void MCintegrate::step(Int nstep) {
  Int i,j;
  for (i=0;i<nstep;i++) {
    for (j=0;j<ndim;j++)
      x[j] = xlo[j]+(xhi[j]-xlo[j])*ran.doub();
    if (xmapp) xx = (*xmapp)(x);
    else xx =x;
    if ((*inregionp)(xx)) {
      fn = (*funcsp)(xx);
      for (j=0;j<nfun;j++) {
	sf[j] += fn[j];
	sferr[j] += SQR(fn[j]);
      }
    }
  } 
  n+= nstep;
}

// Eq 7.7.2
void MCintegrate::calcanswers(){
  for (Int j=0;j<nfun;j++) {
    ff[j] = vol*sf[j]/n;
    fferr[j] = vol*sqrt((sferr[j]/n-SQR(sf[j]/n))/n);
  }
}
