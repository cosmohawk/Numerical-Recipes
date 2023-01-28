void usrfun(VecDoub_I &x, VecDoub_O &fvec, MatDoub_O &fjac);

void mnewt(const Int ntrial, VecDoub_IO &x, const Doub tolx, const Doub tolf)
/* Given an initial guess x[0..n-1] for a root in n dimensions, take ntrial 
   Newton-Raphson steps to improve the root. Stop if the root converges in 
   either summed absolute variable increments tolx or summed absolute function
   values tolf. */
{
  Int i,n=x.size();
  VecDoub p(n),fvec(n);
  MatDoub fjac(n,n);
  
  for (Int k=0;k<ntrial;k++) {
    
    /* User function supplies function values at x in fvec and Jacobian matrix
       in fjac. */
    usrfun(x,fvec,fjac);
    
    // Check function convergence;
    Doub errf = 0.0;
    for (i=0;i<n;i++) errf += abs(fvec[i]);
    if (errf <= tolf) return;

    // RHS of linear equations.
    for (i=0;i<n;i++) p[i] = -fvec[i];

    // Solve linear equations using LU decomposition.
    LUdcmp alu(fjac);
    alu.solve(p,p);

    // Update solution.
    Doub errx = 0.0;
    for (i=0;i<n;i++) {
      errx += abs(p[i]);
      x[i] += p[i];
    }

    // Check root convergence
    if (errx <= tolx) return;
  }
  return;
}

