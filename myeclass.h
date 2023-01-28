void eclass(VecInt_O &nf, VecInt_I &lista, VecInt_I &listb)
/* Given m equivalences between pairs of n individual elements in the form of 
   the input arrays lista[0..m-1] and listb[0..m-1], this routine returns in 
   nf[0..n-1] the number of the equivalence class of each of the n elements, 
   integers between 0 and n-1 (not all such integers used). */
{ 

  Int l,k,j,n=nf.size(),m=lista.size();

  // Initialise each element its own class.
  for (k=0;k<n;k++) nf[k]=k;

  // For each piece of input information...
  for (l=0; l<m; l++) {
    j = lista[l];
    // Track first element up to its ancestor.
    while (nf[j] != j) j = nf[j];
    k = listb[l];
    // Track second element up to its ancestor.
    while (nf[k] != k) k = nf[k];
    // If they are not already related, make them so.
    if (j != k) nf[j] = k;
  }

  // Final sweep up to highest ancestors.
  for (j=0; j<n; j++) 
    while (nf[j] != nf[nf[j]]) nf[j] = nf[nf[j]];
}

void eclazz(VecInt_O &nf, Bool equiv(const Int, const Int))
/* Given a user-supplied boolean function equiv that tells whether a pair of 
   elements, each in the range 0..n-1, are related, return in nf[0..n-1] 
   equivalence class numbers for each element. */
{ 
  Int kk,jj,n=nf.size();
  nf[0] = 0;
  
  // Loop over first element of all pairs.
  for (jj=1; jj<n; jj++) {
    nf[jj] = jj;
    
    // Loop over second element of all pairs.
    for (kk=0;kk<jj;kk++) {
      // Sweep it up this much.
      nf[kk] = nf[nf[kk]];
      if (equiv(jj+1,kk+1)) nf[nf[nf[kk]]]=jj;
      /* Good exercise for the reader to figure out why this much ancestry is 
	 necessary. */
    }
  }

  // Only this much sweeping is needed finally.
  for (jj=0;jj<n;jj++) nf[jj] = nf[nf[jj]];
}


    
