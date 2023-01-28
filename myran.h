struct Ran {
  
  /* Implementation of the highest quality recommended generator. The constructor is called with an integer seed and creates an instance of the generator. The member functions int64, doub, and int32 return the next values in the random sequence, as a variable type indicated by their names. The period of the generator is ~3.138e57 */
  
  unsigned long long u, v, w;
  
Ran(unsigned long long j) : v(4101842887655102017LL), w(1) {
  /* Constructor. Call with any integer seed (except value of v above). */
  u = j ^ v; int64();
  v = u; int64();
  w = v; int64();
}
  
  inline unsigned long long int64() {
    /* Return 64 bit random integer. */
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }
  
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  /* Return a random double-precision floating value in the range 0. to 1. */
  
  inline unsigned int int32() { return(unsigned int)int64(); }
  /* Return 32-bit random integer */
  
};

struct Ranq1 {
  /* Recommended generator for everyday use. The period is ~ 1.8e19. Calling conventions same as Ran, above. */

  Ullong v;

Ranq1(Ullong j) : v(4101842887655102017LL) {
  v ^= j;
  v = int64();
}

  inline Ullong int64() {
    v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
    return v * 2685821657736338717LL;
  }
  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
};

struct Ranq2 {
  /* Backup generator if Ranq1 has too short a period and Ran is too slow. 
     The period is ~8.5e37. Calling conventions same as Ran, above. */

  Ullong v,w;
  
Ranq2(Ullong j) : v(4101842887655102017LL), w(1) {
  v ^= j;
  w = int64();
  v = int64();
}

  inline Ullong int64() {
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    return v ^ w;
  }

  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint) int64(); }
};

/* The following code, added to any of the above generators, augments them with an int8() method. 
   (Be sure to initilize bc to zero in the constructor.) */

// Ullong breg;
// Int bc;
// inline unsigned char int8() {
//     if (bc--) return (unsigned char)(breg >>= 8);
//     breg = int64();
//     bc = 7;
//     return (unsigned char) breg;
// }

/* Random Hashes and Random Bytes */

struct Ranhash {
  /* High-quality random hash of an integer into several numeric types. */
  
  inline Ullong int64(Ullong u) {
    // Returns hash of u as a 64-bit integer.
    
    Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
    v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
    v *= 4768777513237032717LL;
    v ^= v << 20; v ^= v >> 41; v ^= v << 5;
    return v;
  }

  inline Uint int32(Ullong u)
  // Returns hash of u as a 32-bit integer.
  { return (Uint)(int64(u) & 0xffffffff) ; }

  inline Doub doub(Ullong u)
  // Returns hash of u as a double-precision floating value between 0. and 1.
  { return 5.4210108624275221E-20 * int64(u); }

};

struct Ranbyte {
  /* Generator for random bytes using the algorithm generally known as RC4 */

  Int s[256],i,j,ss;
  Uint v;

  Ranbyte(Int u) {
    // Constructor. Call with any integer seed.
    v = 2244614371U ^ u;
    for (i=0; i<256; i++) {s[i] = i;}
    for (j=0, i=0; i<256; i++) {
      ss = s[i];
      j = (j + ss + (v >> 24)) & 0xff;
      s[i] = s[j]; s[j] = ss;
      v = (v << 24) | (v >> 8);
    }
    i = j = 0;
  }

  inline unsigned char int8() {
    // Returns the next random byte in the sequence.
    i = (i+1) & 0xff;
    ss = s[i];
    j = (j+ss) & 0xff;
    s[i] = s[j]; s[j] = ss;
    return (unsigned char)(s[(s[i]+s[j]) & 0xff]);
  }

  Uint int32() {
    // Returns a random 32-bit integer constructed from 4 random bytes. Slow!
    v = 0;
    for (int k=0; k<4; k++) {
      i = (i+1) & 0xff;
      ss = s[i];
      j = (j+ss) & 0xff;
      s[i] = s[j]; s[j] = ss;
      v = (v << 8) | s[(s[i]+s[j]) & 0xff];
    }
    return v;
  }

  Doub doub() { 
    // Returns a random double precision floating value between 0. and 1. Slow!!
    return 2.328306443653869629E-10 * ( int32() +
					2.328306443653869629E-10 * int32() );
  }
};

struct Ranfib {
  /* Implements Knuth's subtractive generator using only floating operations. */

  Doub dtab[55], dd;
  Int inext, inextp;
  
Ranfib(Ullong j) : inext(0), inextp(31) {
  // Constructor. Call with any integer seed. Uses Ranq1 to initialize.
  Ranq1 init(j);
  for (int k=0; k<55; k++) dtab[k] = init.doub();
}

  Doub doub() {
    // Returns random double-precision floating value between 0. and 1.
    if (++inext == 55) inext = 0;
    if (++inextp == 55) inextp = 0;
    dd = dtab[inext] - dtab[inextp];
    if (dd < 0) dd += 1.0;
    return (dtab[inext] = dd);
  }

  inline unsigned long int32()
  // Returns random 32-bit integer. Recommended only for testing purposes.
  { return (unsigned long)(doub() * 4294967295.0);}
};
