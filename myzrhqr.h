void zrhqr(VecDoub_I &a, VecComplex_O &rt)
/* Find all roots of a polynomial with real coeficients, 
   \sum_{i=0}^m a[i] x^i, given the coefficients a[0..m]. The method is to 
   construct an upper Hessenberg matrix whose eigenvalues are the desired roots
   and then use the routine Unsymeig. The roots are returned in the complex 
   vector rt[0..m-1], sorted in descending order by their real parts. */
{
  Int m=a.size()-1;
  MatDoub hess(m,m);

  // Construct the matrix
  for (Int k=0;k<m;k++) {
    hess[0][k] = -a[m-k-1]/a[m];
    for (Int j=1;j<m;j++) hess[j][k]=0.0;
    if (k!=m-1) hess[k+1][k]=1.0;
  }

  // Find its eigenvalues.
  Unsymmeig h(hess, false, true);
  for (Int j=0;j<m;j++)
    rt[j]=h.wri[j];
}

