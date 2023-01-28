/* ch 7.3. Deviates from Other Distributions */

struct Expondev : Ran {
  // Structure for exponential deviates.
  
  Doub beta;
  
 Expondev(Doub bbeta, Ullong i) : Ran(i), beta(bbeta) {}
  // Constructor arguments are beta and a random sequence seed.
  
  Doub dev() {
    // Return an exponential deviate.
    Doub u;
    do u = doub(); while (u == 0.);
    return -log(u)/beta;
  }
};

struct Logisticdev : Ran {
  // Structure for logistic deviates.
  
  Doub mu, sig;

 Logisticdev(Doub mmu, Doub ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}
  //Constructor arguments are mu, sig, and a random sequence seed.

  Doub dev() {
    // Return a logistic deviate.
    Doub u;
    do u = doub(); while (u*(1.-u) == 0.);
    return mu + 0.551328895421792050*sig*log(u/(1.-u));
  }
};

struct Normaldev_BM : Ran {
  /* Structure for normal deviates. */

  Doub mu, sig;
  Doub storedval;

  Normaldev_BM(Doub mmu, Doub ssig, Ullong i)
    : Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
  // Constructor arguments are mu, sigma, and a random sequence seed.

  Doub dev() {
    // Return a normal deviate.
    Doub v1, v2, rsq, fac;
    /* We don't have an extra deviate handy, so pick two uniform
       random numbers in the square extending from -1 to +1 in
       each direction, see if they are in the unit circle,
       or try again. */
    if (storedval == 0.) {
      /* We don't have an extra deviate handy, so pick two uniform
	 random numbers in the square extending from -1 to +1 in 
	 each direction, see if they are in the unit circle, 
	 or try again. */
      do {
	v1=2.0*doub()-1.0; 
	v2=2.0*doub()-1.0;
	rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      /* Now we make the Box-Muller transformation to get two normal
	 deviates. Return one and save the other for next time. */
      fac = sqrt(-2.0*log(rsq)/rsq);
      storedval = v1*fac;
      return mu + sig*v2*fac;
    } else {
      /* We have an extra deviate handy, so return it */
      fac = storedval;
      storedval = 0.;
      return mu + sig*fac;
    }
  }
};

struct Cauchydev : Ran {
  // Structure for Cauchy deviates.
  Doub mu, sig;
  
 Cauchydev(Doub mmu, Doub ssig, Ullong i) 
    : Ran(i), mu(mmu), sig(ssig) {}
  // Constructor arguments are mu, sigma, and a random sequence seed
  
  Doub dev() {
    // Return a Cauchy deviate.
    Doub v1, v2;
    // Find a random point in the unit semicircle
    do {
      v1=2.0*doub()-1.0;
      v2=doub();
    } while (SQR(v1)+SQR(v2) >= 1. || v2 == 0.);
    // Ratio of its coordinates is the tangent of a random angle
    return mu + sig*v1/v2;
  }
};

struct Normaldev : Ran {
  // Structure for normal deviates.

  Doub mu, sig;
  
 Normaldev(Doub mmu, Doub ssig, Ullong i)
   : Ran(i), mu(mmu), sig(ssig) {}
  // Constructor arguments are mu, sigma, and a random sequence seed.

  Doub dev() {
    // Return a normal deviate.
    Doub u,v,x,y,q;
    do {
      u = doub();
      v = 1.7156*(doub()-0.5);
      x = u - 0.449871;
      y = abs(v) + 0.386595;
      q = SQR(x) + y*(0.19600*y-0.25472*x);
    } while (q > 0.27597
	     && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
    return mu + sig*v/u;
  }
};

struct Gammadev : Normaldev {
  // Structure for gamma deviates.

  Doub alph, oalph, bet;
  Doub a1, a2;

 Gammadev(Doub aalph, Doub bbet, Ullong i)
   : Normaldev(0.,1.,i), alph(aalph), oalph(aalph), bet(bbet) {
    // Constructor arguments are alpha, beta, and a random sequence seed.
    if (alph <= 0.) throw(" bad alph in Gammadev");
    if (alph < 1.) alph += 1.;
    a1 = alph -1./3.;
    a2 = 1./sqrt(9.*a1);
  }
  
  Doub dev() {
    // Return a gamma deviate by the method of Marsaglia and Tang.
    Doub u,v,x;
    do {
      do {
	x = Normaldev::dev();
	v = 1. + a2*x;
      } while (v <= 0.);
      v = v*v*v;
      u = doub();
    } while (u > 1. - 0.331*SQR(SQR(x)) &&
	     log(u) > 0.5*x*x + a1*(1.-v+log(v))); // Rarely evaluated
    if (alph == oalph) return a1*v/bet;
    else {
      // Case where alpha < 1, per Ripley.
      do u =doub(); while (u==0.);
      return pow(u,1./oalph)*a1*v/bet;
    }
  }
};

struct Poissondev : Ran {
  /* Structure for Poisson deviates */
  
  Doub lambda, sqlam, loglam, lamexp, lambold;
  VecDoub logfact;
  
 Poissondev(Doub llambda, Ullong i) : Ran(i), lambda(llambda),
    logfact(1024,-1.), lambold(-1.) {}
  // Constructor arguments are lambda and a random sequence seed.

  Int dev() {
    // Return a Poisson deviate using the most recently set value of lambda.
    Doub u,u2,v,v2,p,t,lfac;
    Int k;

    if (lambda < 5.) {
      // Will use product of uniforms method.
      if (lambda != lambold) lamexp=exp(-lambda);
      k = -1;
      t = 1.;
      do {
	++k;
	t *= doub();
      } while (t > lamexp);
    } else {
      // Will use ratio-of-uniforms method.
      if (lambda != lambold) {
	sqlam = sqrt(lambda);
	loglam = log(lambda);
      }
      for (;;) {
	u = 0.64*doub();
	v = -0.68 + 1.28*doub();
	if (lambda > 13.5) {
	  // Outer squeeze for fast rejection.
	  v2 = v*v;
	  if (v >= 0.) {if (v2 > 6.5*u*(0.64-u)*(u+0.2)) continue;}
	  else {if (v2 > 9.6*u*(0.66-u)*(u+0.07)) continue;}
	}
	k = Int(floor(sqlam*(v/u)+lambda+0.5));
	if (k < 0) continue;
	u2 = SQR(u);
	if (lambda > 13.5) {
	  // Inner squeeze for fast acceptance.
	  if (v >= 0.) {if (v2 < 15.2*u2*(0.61-u)*(0.8-u)) break;}
	  else {if (v2 < 6.76*u2*(0.62-u)*(1.4-u)) break;}
	}
	if (k < 1024) {
	  if (logfact[k] < 0.) logfact[k] = gammln(k+1.);
	  lfac = logfact[k];
	} else lfac = gammln(k+1.);
	// Only when we must
	p = sqlam*exp(-lambda +k*loglam -lfac);
	if (u2 < p) break;
      }
    }
    lambold = lambda;
    return k;
  }
  Int dev(Doub llambda) {
    // Reset lambda and then return a Poisson deviate.
    lambda = llambda;
    return dev();
  }
};
	
struct Binomialdev : Ran {
  // Structure for binomial deviates.
  
  Doub pp,p,pb,expnp,np,glnp,plog,pclog,sq;
  Int n,swch;
  Ullong uz,uo,unfin,diff,rltp;
  Int pbits[5];
  Doub cdf[64];
  Doub logfact[1024];

 Binomialdev(Int nn, Doub ppp, Ullong i) : Ran(i), pp(ppp), n(nn) {
    // Constructor arguments are n, p, and a random sequence seed.
    
    Int j;
    pb = p = (pp <= 0.5 ? pp : 1.0-pp);
    if (n <= 64) { 
      // Will use bit-parallel direct method.
      uz = 0;
      uo = 0xffffffffffffffffLL;
      rltp = 0;
      for (j=0;j<5;j++) pbits[j]=1 & ((Int)(pb *= 2.));
      // Leading bits of p (above) and remaining fraction.
      p -= floor(pb);
      swch = 0;
    } else if (n*p < 30.) {
      // Will use precomputed cdf table.
      cdf[0] = exp(n*log(1-p));
      for (j=1;j<64;j++) { 
	cdf[j] = cdf[j-1] + exp(gammln(n+1.)
				- gammln(j+1.)-gammln(n-j+1.)
				+ j*log(p)+(n-j)*log(1.-p));
      }
      swch = 1;
    } else {
      // Will use ratio of uniforms method.
      np = n*p;
      glnp=gammln(n+1.);
      plog=log(p);
      pclog=log(1.-p);
      sq=sqrt(np*(1.-p));
      if (n<1024) for (j=0;j<=n;j++) logfact[j] = gammln(j+1.);
      swch=2;
    }
  }

  Int dev() {
    // Return a binomial deviate.
    Int j,k,kl,km;
    Doub y,u,v,u2,v2,b;
    if (swch == 0) {
      // Mark all bits as "unfinished".
      unfin = uo;
      // Compare with first 5 bits of p.
      for (j=0;j<5;j++) {
	// Mask of diff.
	diff = unfin & (int64()^(pbits[j]? uo : uz));
	// Set bits to 1, meaning ran < p.
	if (pbits[j]) rltp |= diff;
	// Set bits to 0, meaning ran > p.
	else rltp=rltp & ~diff;
	// Update unfinished status.
	unfin = unfin & ~diff;
      }
      // Now we just count the events.
      k=0;
      for (j=0;j<n;j++) {
	// Clean up unresolved cases,
	if (unfin & 1) {if (doub() < pb) ++k;}
	// or use bit answer.
	else {if (rltp & 1) ++k;}
	unfin >>= 1;
      }
    } else if (swch == 1) {
      // Use stored cdf.
      y = doub();
      kl = -1;
      k = 64;
      while (k-kl>1) {
	km = (kl+k)/2;
	if (y < cdf[km]) k = km;
	else kl = km;
      }
    } else {
      // Use ratio of uniforms method.
      for (;;) {
	u = 0.645*doub();
	v = -0.63 + 1.25*doub();
	v2 = SQR(v);
	// Try squeeze for fast rejection.
	if (v >= 0.) {if (v2 > 6.5*u*(0.645-u)*(u+0.2)) continue;}
	else {if (v2 > 8.4*u*(0.645-u)*(u+0.1)) continue;}
	k = Int(floor(sq*(v/u)+np+0.5));
	if (k < 0) continue;
	u2 = SQR(u);
	// Try squeeze for fast acceptance.
	if (v >= 0.) {if (v2 < 12.25*u2*(0.615-u)*(0.92-u)) break;}
	else {if (v2 < 7.84*u2*(0.615-u)*(1.2-u)) break;}
	// Only when we must
	b = sq*exp(glnp+k*plog+(n-k)*pclog 
		   - (n < 1024 ? logfact[k]+logfact[n-k]
		      : gammln(k+1.)+gammln(n-k+1.)));
	if (u2 < b) break;
      }
    }
    if (p != pp) k = n - k;
    return k;
  }
};

      
