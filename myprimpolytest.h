struct Primpolytest {
  /* Test polynomials over the integers mod 2 for primitiveness. */

  Int N, nfactors;
  VecUllong factors;
  VecInt t,a,p;

Primpolytest() : N(32), nfactors(5), factors(nfactors), t(N*N), a(N*N), p(N*N) {
  // Constructor. The constants are specific to 32-bit LFSRs.
  Ullong factordata[5] = {3,5,17,257,65537};
  for (Int i=0;i<nfactors;i++) factors[i] = factordata[i];
}

  // Utility to test if p is the identity matrix.
  Int ispident() {
    Int i,j;
    for (i=0; i<N; i++) for (j=0; j<N; j++) {
	if (i==j) {if (p[i*N+j] !=1) return 0;}
	else {if (p[i*N+j] != 0) return 0;}
      }
    return 1;
  }

  // Utility for a *= b on matrices a and b.
  void mattimeseq(VecInt &a, VecInt &b) {
    Int i,j,k,sum;
    VecInt tmp(N*N);
    for (i=0;i<N;i++) for (j=0;j<N;j++) {
	sum=0;
	for (k=0;k<N;k++) sum += a[i*N+k] * b[k*N+j];
	tmp[i*N+j] = sum & 1;
      }
    for (k=0;k<N*N; k++) a[k] = tmp[k];
  }

  // Utility for matrix p = a^n by successive squares.
  void matpow(Ullong n) {
    Int k;
    for (k=0; k<N*N; k++) p[k] = 0;
    for (k=0; k<N; k++) p[k*N+k] = 1;
    while (1) {
      if (n & 1) mattimeseq(p,a);
      n >>= 1;
      if (n == 0) break;
      mattimeseq(a,a);
    }
  }

  /* Main test routine. 
     Returns 1 if the polynomial with serial number n is primitive,
     0 otherwise */
  Int test(Ullong n) {
    Int i,j,k;
    Ullong pow, tnm1, nn = n;
    tnm1 = ((Ullong)1 << N) -1;
    
    if (n > (tnm1>>1)) throw("not a polynomial of degree N");
    
    // Construct the update matrix in t.
    for (k=0;k<N*N; k++) t[k]=0;
    for (i=1;i<N;i++) t[i*N+(i-1)]=1;
    j=0;
    while (nn) {
      if (nn & 1) t[j] = 1;
      nn >>= 1;
      j++;
    }
    t[N-1] = 1;
    
    // Test that t^tnm1 is the indentity matrix.
    for (k=0; k<N*N; k++) a[k] = t[k];
    matpow(tnm1);
    if (ispident() != 1) return 0;
    
    // Test that the t to the required submultiple 
    // powers is not the identity matrix.
    for (i=0; i<nfactors; i++) {
      pow = tnm1/factors[i];
      for (k=0; k<N*N; k++) a[k] = t[k];
      matpow(pow);
      if (ispident() == 1) return 0;
    }
    return 1;
  }
};

