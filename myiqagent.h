struct IQagent {
  /* Object for estimating arbitrary quantile values from a continuous stream
     of data values. */
  
  // Batch size. x10 if > 10^6 data values are expected.
  static const Int nbuf = 1000;
  
  Int nq, nt, nd;
  VecDoub pval, dbuf, qile;
  Doub q0, qm;
  
IQagent() : nq(251), nt(0), nd(0), pval(nq), dbuf(nbuf), 
    qile(nq,0.), q0(1.e99), qm(-1.e99) {
  // Constructor. No arguments.
  for (Int j=85; j<=165; j++) pval[j] = (j-75.)/100.;
  // Set general purpose array of p-values ranging from 10^-6 to 1-10^-6.
  // You can change this if you want.
  for (Int j=84; j>=0; j--) {
    pval[j] = 0.87191909*pval[j+1];
    pval[250-j] = 1.-pval[j];
  }
}

  void add(Doub datum) { 
    // Assimilate a new value from the stream.
    dbuf[nd++] = datum;
    if (datum < q0) {q0 = datum;}
    if (datum > qm) {qm = datum;}
    if (nd == nbuf) update();
    // Time for a batch update.
  }

  void update() {
    /* Batch update, as shown in Figure 8.5.1. This function is called by add
       or report and should not be called directly by the user. */
    Int jd=0, jq=1, iq;
    Doub target, told=0., tnew=0., qold, qnew;
    // Will be new quantiles after update.
    VecDoub newqile(nq);
    sort(dbuf,nd);
    // Set lowest and highest to min and max values seen so far.
    // Set compatible p-values.
    qold = qnew = qile[0] = newqile[0] = q0;
    qile[nq-1] = newqile[nq-1] = qm;
    pval[0] = min(0.5/(nt+nd),0.5*pval[1]);
    pval[nq-1] = max(1.-0.5/(nt+nd),0.5*(1.+pval[nq-2]));
    // Main loop over target p-values for interpolation.
    for (iq=1;iq<nq-1;iq++) {
      target = (nt+nd)*pval[iq];
      if (tnew < target) for (;;) {
	  /* Here's the guts: We locate a succession of abscissa-ordinate
	     pairs (qnew,tnew) that are the discontinuities of value or slope
	     in Figure 8.5.1(c), breaking to perform an interpolation as we 
	     cross each target. */
	  if (jq < nq && (jd >= nd || qile[jq] < dbuf[jd])) {
	    // Found slope discontinuity from old CDF.
	    qnew = qile[jq];
	    tnew = jd + nt*pval[jq++];
	    if (tnew >= target) break;
	  } else {
	    // found value discontinuity from batch data CDF.
	    qnew = dbuf[jd];
	    tnew = told;
	    if (qile[jq]>qile[jq-1]) tnew += nt*(pval[jq]-pval[jq-1])
				       *(qnew-qold)/(qile[jq]-qile[jq-1]);
	    jd++;
	    if (tnew >= target) break;
	    told = tnew++;
	    qold = qnew;
	    if (tnew >= target) break;
	  }
	  told = tnew;
	  qold = qnew;
	}
      // Break to here and perform the new interpolation.
      if (tnew == told) newqile[iq] = 0.5*(qold+qnew);
      else newqile[iq] = qold + (qnew-qold)*(target-told)/(tnew-told);
      told = tnew;
      qold = qnew;
    }
    qile = newqile;
    nt += nd;
    nd = 0;
  }

  Doub report(Doub p) {
    /* Return estimated p-quantile for the data seen so far. 
       E.g. p = 0.5 for the median. */
    Doub q;
    if (nd > 0) update(); // You may want to remove this line, see text.
    Int jl=0, jh=nq-1,j;
    // Locate place in table by bisection.
    while (jh-jl>1) {
      j = (jh+jl)>>1;
      if (p > pval[j]) jl=j;
      else jh=j;
    }
    // Interpolate.
    j = jl;
    q = qile[j] + (qile[j+1]-qile[j])*(p-pval[j])/(pval[j+1]-pval[j]);
    return MAX(qile[0],MIN(qile[nq-1],q));
  }
};

