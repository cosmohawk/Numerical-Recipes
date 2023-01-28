template<class T>
void piksrt(NRvector<T> &arr)
/* Sort an array arr[0..n-1] into ascending numerical order, by straight 
   insertion. arr is replaced on output by its sorted rearrangement. Should 
   only really be used by arrays of small N ~< 20. */
{
  Int i, j, n=arr.size();
  T a;
  
  // Pick out each element in turn.
  for (j=1;j<n;j++) {
    a = arr[j];
    i = j;
    
    // Look for a place to insert it.
    while (i > 0 && arr[i-1] > a) {
      arr[i] = arr[i-1];
      i--;
    }
    
    // Insert it.
    arr[i] = a;
  }
}

template<class T, class U>
  void piksr2(NRvector<T> &arr, NRvector<U> &brr)
  /* Sort an array arr[0..n-1] into ascending numerical order, by straight 
     insertion, whilst making the corresponding rearrangement of the array 
     brr[0..n-1] */
{
  Int i, j, n=arr.size();
  T a;
  U b;

  // Pick out each element in turn.
  for (j=1; j<n; j++) {
    a = arr[j];
    b = brr[j];
    i = j;
    // Look for the place to insert it.
    while (i > 0 && arr[i-1] > a) {
      arr[i] = arr[i-1];
      brr[i] = brr[i-1];
      i--;
    }

    // Insert it.
    arr[i] = a;
    brr[i] = b;
  }
}

template<class T>
void shell(NRvector<T> &a, Int m=-1)
/* Sort an array a[0..n-1] into ascending numerical order by Shell's method
   (diminishing increment sort). a is replaced on output by its sorted 
   rearrangement. Normally, the optional argument m should be omitted, but if
   it is set to a +ve value then only the first m elements of a are sorted. */
{
  Int i, j, inc, n=a.size();
  T v;
  
  // Use optional argument.
  if (m>0) n=MIN(m,n);

  // Determine the starting increment.
  inc=1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  
  // Loop over partial sorts.
  do {
    inc /= 3;
    // Outer loop of straight insertion.
    for (i=inc;i<n;i++) {
      v = a[i];
      j = i;
      // Inner loop of straight insertion.
      while (a[j-inc] > v) {
	a[j] = a[j-inc];
	j -= inc;
	if (j < inc) break;
      }
      a[j] = v;
    }
  } while (inc > 1);
}

template<class T>
void sort(NRvector<T> &arr, Int m=-1)
/* Sort an array arr[0..n-1] into ascending numerical order using the
   Quicksort algorithm. arr is replaced on output by its sorted rearrangement. 
   Normally, the optional argument m should be omitted, but if it is set to a 
   positive value, then only the first m elements of arr are sorted. */
{
  static const Int M=7, NSTACK=64;
  /* Here M is the size of subarrays sorted by straight insertion and NSTACK is
     the required auxiliary storage. */
  Int i,ir,j,k,jstack=-1,l=0,n=arr.size();
  T a;
  VecInt istack(NSTACK);

  // Use optional argument.
  if (m>0) n = MIN(m,n);
  ir=n-1;
  for (;;) {
    if (ir-l < M) {
      // Insertion sort when subarray small enough.
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack < 0) break;
      // Pop stack and begin new round of partitioning.
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      /* Choose median of left, centre, and right elements as partitioning
	 element a. Also rearrange so that a[l] <= a[l+1] <= a[ir] */
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      // Initialise pointers for positioning.
      i=l+1;
      j=ir;
      // Partitioning element.
      a=arr[l+1];
      // Beginning of innermost loop.
      for (;;) {
	// Scan up to find element > a.
	do i++; while (arr[i] < a);
	// Scan down to find element < a.
	do j--; while (arr[j] > a);
	if (j < i) break;
	// Pointers crossed. Partitioning complete.
	// Exchange elements.
	SWAP(arr[i],arr[j]);
      }
      // Insert partitioning element.
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack += 2;
      /* Push pointers to larger subarray on stack; 
	 process smaller subarray immediately. */
      if (jstack >= NSTACK) throw("NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}
  
template<class T, class U>
void sort2(NRvector<T> &arr, NRvector<U> &brr)
/* Sort an array arr[0..n-1] into ascending order using Quicksort, whilst 
   making the corresponding rearrangement of the array brr[0..n-1] */
{
  static const Int M=7, NSTACK=64;
  /* Here M is the size of subarrays sorted by straight insertion and NSTACK is
     the required auxiliary storage. */
  Int i, ir, j, k, jstack=-1,l=0,n=arr.size();
  T a;
  U b;
  VecInt istack(NSTACK);
  ir=n-1;
  for (;;) {
    if (ir-l < M) {
      // Insertion sort when subarray small enough.
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
	b=brr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
        }
        arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack < 0) break;
      // Pop stack and begin new round of partitioning.
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      /* Choose median of left, centre, and right elements as partitioning
         element a. Also rearrange so that a[l] <= a[l+1] <= a[ir] */
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      SWAP(brr[k],brr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
	SWAP(brr[l],brr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
	SWAP(brr[l+1],brr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
	SWAP(brr[l],brr[l+1]);
      }
      // Initialise pointers for positioning.
      i=l+1;
      j=ir;
      // Partitioning element.
      a=arr[l+1];
      b=brr[l+1];
      // Beginning of innermost loop. 
      for (;;) {
        // Scan up to find element > a.
	do i++; while (arr[i] < a);
        // Scan down to find element < a.
	do j--; while (arr[j] > a);
        if (j < i) break;
        // Pointers crossed. Partitioning complete.
	// Exchange elements of both arrays.
	SWAP(arr[i],arr[j]);
	SWAP(brr[i],brr[j]);
      }
      // Insert partitioning element in both arrays.
      arr[l+1] = arr[j];
      brr[l+1] = brr[j];
      arr[j] = a;
      brr[j] = b;
      jstack += 2;
      /* Push pointers to larger subarray on stack; 
	 process smaller subarray immediately. */
      if (jstack >= NSTACK) throw("NSTACK too small in sort2.");
      if (ir-i+1 >= j-1) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template<class T, class U, class V>
  void sort3(NRvector<T> &arr, NRvector<U> &brr, NRvector<V> &crr)
  /* Sort an array arr[0..n-1] into ascending order using Quicksort, whilst 
     making the corresponding rearrangement of the two arrays brr[0..n-1]
     and crr[0..n-1]*/
{
  static const Int M=7, NSTACK=64;
  /* Here M is the size of subarrays sorted by straight insertion and NSTACK is
     the required auxiliary storage. */
  Int i, ir, j, k, jstack=-1,l=0,n=arr.size();
  T a;
  U b;
  V c;
  VecInt istack(NSTACK);
  ir=n-1;
  for (;;) {
    if (ir-l < M) {
      // Insertion sort when subarray small enough.
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
	b=brr[j];
	c=crr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
          brr[i+1]=brr[i];
	  crr[i+1]=crr[i];
        }
        arr[i+1]=a;
        brr[i+1]=b;
	crr[i+1]=c;
      }
      if (jstack < 0) break;
      // Pop stack and begin new round of partitioning.
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      /* Choose median of left, centre, and right elements as partitioning
         element a. Also rearrange so that a[l] <= a[l+1] <= a[ir] */
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      SWAP(brr[k],brr[l+1]);
      SWAP(crr[k],crr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        SWAP(brr[l],brr[ir]);
	SWAP(crr[l],crr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        SWAP(brr[l+1],brr[ir]);
	SWAP(crr[l+1],crr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        SWAP(brr[l],brr[l+1]);
	SWAP(crr[l],crr[l+1]);
      }
      // Initialise pointers for positioning.
      i=l+1;
      j=ir;
      // Partitioning element.
      a=arr[l+1];
      b=brr[l+1];
      c=crr[l+1];
      // Beginning of innermost loop. 
      for (;;) {
        // Scan up to find element > a.
        do i++; while (arr[i] < a);
        // Scan down to find element < a.
        do j--; while (arr[j] > a);
        if (j < i) break;
        // Pointers crossed. Partitioning complete.
        // Exchange elements of both arrays.
        SWAP(arr[i],arr[j]);
        SWAP(brr[i],brr[j]);
	SWAP(crr[i],crr[j]);
      }
      // Insert partitioning element in both arrays.
      arr[l+1] = arr[j];
      brr[l+1] = brr[j];
      crr[l+1] = crr[j];
      arr[j] = a;
      brr[j] = b;
      crr[j] = c;
      jstack += 2;
      /* Push pointers to larger subarray on stack; 
         process smaller subarray immediately. */
      if (jstack >= NSTACK) throw("NSTACK too small in sort3.");
      if (ir-i+1 >= j-1) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}


namespace hpsort_util {
  template<class T>
    void sift_down(NRvector<T> &ra, const Int l, const Int r)
    /* Carry out the sift-down on element ra(l) to maintain the heap structure.
       l and r determine the "left" and "right" range of the sift-down. */
    {
      Int j,jold;
      T a;
      a = ra[l];
      jold=l;
      j=2*l+1;
      while (j <= r) {
	// Compare to the better underling.
	if (j < r && ra[j] < ra[j+1]) j++;
	// Find a's level. Terminate the sift-down.
	if (a >= ra[j]) break;
	// Otherwise, demote a and continue.
	ra[jold] = ra[j];
	jold = j;
	j = 2*j+1;
      }
      // Put a into its slot.
      ra[jold]=a;
    }
}


template<class T>
void hpsort(NRvector<T> &ra)
/* Sort in array ra[0..n-1] into ascending numerical order using the Heapsort 
   algorithm. ra is replaced on output by its sorted rearrangement. */
{
  Int i, n=ra.size();
  for (i=n/2-1; i>=0; i--)
    /* The index i, which here determines the "left" range of the sift-down, 
       i.e., the element to be sifted down, is decremented from n/2-1 down to
       0 during the "hiring" (heap creation) phase. */
    hpsort_util::sift_down(ra,i,n-1);
  for (i=n-1; i>0; i--) {
    /* Here the "right" range of the sift-down is decremented from n-2 down to
       0 during the "retirement-and-promotion" (heap selection) phase. */
    // Clear a space at the end of the array.
    SWAP(ra[0],ra[i]);
    // Retire the top of the heap into it.
    hpsort_util::sift_down(ra,0,i-1);
  }
}

struct Indexx {

  Int n;
  VecInt indx;

  template<class T> Indexx(const NRvector<T> &arr) {
    // Constructor. Calls index and stores an index to the array arr[0..n-1].
    index(&arr[0],arr.size());
  }
  Indexx() {} // Empty constructor. See text.

  template<class T> void sort(NRvector<T> &brr) {
    /* Sort an array brr[0..n-1] into the order defined by the stored index. 
       brr is replace on output by its sorted rearrangement. */
    if (brr.size() != n) throw("bad size in Index sort");
    NRvector<T> tmp(brr);
    for (Int j=0;j<n;j++) brr[j] = tmp[indx[j]];
  }

  template<class T> inline const T & el(NRvector<T> &brr, Int j) const {
    /* This function, and the next, return the element brr that woulb be in
       sorted position j according to the stored index. The vector brr is not
       changed. */
    return brr[indx[j]];
  }
  template<class T> inline T & el(NRvector<T> &brr, Int j) {
    // Same but return an l-value.
    return brr[indx[j]];
  }

  template<class T> void index(const T *arr, Int nn);
  /* This does the actual work of indexing. Normally not called by the user,
     but see text for exceptions. */

  void rank(VecInt_O &irank) {
    /* Returns a rank table, whose jth element is the rank of arr[j], where 
       arr is the vector originally indexed. The smallest arr[j] has rank 0. */
    irank.resize(n);
    for (Int j=0;j<n;j++) irank[indx[j]] = j;
  }

};

template<class T>
void Indexx::index(const T *arr, Int nn)
/* Indexes an array arr[0..nn-1], i.e. resized and sets indx[0..nn-1] such that
   arr[indx[j]] is in ascending order for j=0,1,...,nn-1. Also sets member 
   value n. The input array arr is not changed. */
{
  const Int M=7, NSTACK=64;
  Int i, indxt,ir,j,k,jstack=-1,l=0;
  T a;
  VecInt istack(NSTACK);
  n = nn;
  indx.resize(n);
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt = indx[j];
	a = arr[indxt];
	for (i=j-1;i>=l;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1] = indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack <0) break;
      ir = istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	  SWAP(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
	SWAP(indx[l],indx[l+1]);
      }
      i = l+1;
      j = ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) throw("NSTACK too small index.");
      if (ir-i+1 >= j-1) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}

template<class T>
T select(const Int k, NRvector<T> &arr)
/* Given k in [0..n-1] returns an array value from arr[0..n-1] such that k 
   array values are less than or requal to the one returned. The input array 
   will be rearranged to have this value in location arr[k], with all smaller
   elements moved to arr[0..k-1] (in arbitrary order) and all larger elements 
   in arr[k+1..n-1] (also in arbitrary order). */
{
  Int i, ir, j, l, mid, n=arr.size();
  T a;
  l=0;
  ir=n-1;
  for (;;) {
    if (ir <= l+1) {
      // Active partition contains 1 or 2 elements.
      if (ir == l+1 && arr[ir] < arr[l])
	// Case of 2 elements.
	SWAP(arr[l],arr[ir]);
      return arr[k];
    } else {
      // Choose median of left, centre and right elements as partitioning one a
      mid=(l+ir) >> 1;
      // Also rearrange so that arr[l] <= arr[l+1], arr[ir] >= arr[l+1].
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir])
	SWAP(arr[l],arr[ir]);
      if (arr[l+1] > arr[ir])
	SWAP(arr[l+1],arr[ir]);
      if (arr[l] > arr[l+1])
	SWAP(arr[l],arr[l+1]);
      // Initialise pointers for partitioning.
      i=l+1;
      j=ir;
      // Partitioning element
      a=arr[l+1];
      // Beginning of innermost loop.
      for (;;) {
	// Scan up to find element > a.
	do i++; while (arr[i] <a);
	// Scan down to find element < a.
	do j--; while (arr[j] > a);
	// Pointers crossed? Partitioning complete.
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      // End of innermost loop.
      // Insert partitioning element.
      arr[l+1] = arr[j];
      arr[j] = a;
      // Keep the partition that conatins the kth element active.
      if (j >= k) ir = j-1;
      if (j <= k) l=i;
    }
  }
}

struct Heapselect {
  /* Object for tracking the m largest values seen thus far in a stream of
     values. */
  Int m,n,srtd;
  VecDoub heap;

Heapselect(Int mm) : m(mm), n(0), srtd(0), heap(mm,1.e99) {}
  // Constructor. The argument is the number of largest values to track.

  void add(Doub val) {
    // Assimilate a new value from the stream.
    Int j,k;
    if (n<m) {
      // Heap not yet filled
      heap[n++] = val;
      // Create initial heap by overkill!
      if (n==m) sort(heap);
    } else {
      if (val > heap[0]) {
	// Put it on heap?
	heap[0] = val;
	// Sift down.
	for (j=0;;) {
	  k=(j << 1) + 1;
	  if (k > m-1) break;
	  if (k != (m-1) && heap[k] > heap[k+1]) k++;
	  if (heap[j] <= heap[k]) break;
	  SWAP(heap[k],heap[j]);
	  j=k;
	}
      }
      n++;
    }
    // Mark heap as "unsorted".
    srtd = 0;
  }

  Doub report(Int k) {
    /* Return the kth largest values seen so far, k=0 returns the largest value
       seen, k=1 the second largest, ..., k=m-1 the last position tracked.
       Also, k must be less than the number of previous values assimilated. */
    Int mm = MIN(m,n);
    if (k > mm-1) throw("Heapselect k too big");
    if (k == m-1) return heap[0];
    // Always free, since top of heap.
    // Otherwise need to sort the heap.
    if (! srtd) { sort(heap); srtd =1;}
    return heap[mm-1-k];
  }
};

