struct Multinormaldev : Ran {
  /* Structure for multivariate normal deviates */

  // Declare variables.
  Int mm;
  VecDoub mean;
  MatDoub var;
  Cholesky chol;
  VecDoub spt, pt;

 Multinormaldev(Ullong j, VecDoub &mmean, MatDoub &vvar) :
  Ran(j), mm(mmean.size()), mean(mmean), var(vvar), chol(var),
    spt(mm), pt(mm) {
    /* Constructor. Arguments are the random generator seed, the (vector) mean,
       and the (matrix) covariance. Cholesky decomposition of the covariance
       is done here. */
    if (var.ncols() != mm || var.nrows() != mm) throw("bad sizes");
  }

  VecDoub &dev() {
    // Return a multivariate normal deviate.
    Int i;
    Doub u,v,x,y,q;
    
    // Fill a vector of independent normal deviates.
    for (i=0;i<mm;i++) {
      do { 
	u = doub();
	v = 1.7156*(doub()-0.5);
	x = u - 0.449871;
	y = abs(v) + 0.386595;
	q = SQR(x) + y*(0.19600*y-0.25472*x);
      } while (q > 0.27597
	       && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
      spt[i] = v/u;
    }

    // Apply equation (7.4.3).
    chol.elmult(spt,pt);
    for (i=0;i<mm;i++) {pt[i] += mean[i];}
    return pt;
  }

};
  
