#include <stdio.h>
#include <stdlib.h>

int euclidean(int a, int b) {
  // Compute the gcd of a and b fast.
    int temp = 0;
  while (b) {
    temp = b;
    b = a%b;
    a = temp;
  }
  return abs(a);
}

int positive_mod(int a, int p) {
  // helper code to always return positive moduli
  return (a%p+p)%p;
}

int inverse_modp(int a, int p) {
  // Return the inverse of a in the Galois field Fp.
  // Return zero if the inverse does not exist, i.e.
  // if p is not prime.
  int u = 1;
  int r = a;
  int x = 0;
  int y = p;
  int q, t, v, s;

  while (y!=0) {
    q = r/y;
    t = r%y;
    s = u - q*x;
    u = x;
    r = y;
    x = s;
    y = t;
  }
  
  v = (r - a*u)/p;
  // r = au + bv and r = gcd(a, b)
  if (r==1)
    return positive_mod(u, p);
  else
    return 0;
}

void superinc_ss(int a[], int x[], int n, int w) {
  // a is the superincreasing set, x is the soln vector
  // n is length(a), w is the sum. Solve fast, linear time.
  // Importantly, this code assumes the solution exists!
  // If not, the contents of x are undefined, probs useless
  int s = w;
  for (int i=n-1; i>=0; i--) {
    if (s >= a[i]) {
      x[i] = 1;
      s -= a[i];
    }
    else
      x[i] = 0;
  }
}
  
void sis_to_pub(int x[], int n, int p, int b) {
  // converts Bob's secret superincreasing set x into
  // the public set. p is a prime greater than the
  // sum of the terms in x, b is a chosen number
  // between 2 and p.
  for (int i=0; i<n; i++) {
    x[i] = positive_mod(x[i]*b, p);
  }
}

int encrypt(int pub_set[], int msg[], int n) {
  // pub_set is the public set published by Bob
  // msg is a binary vector representation of the
  // message, n is the length
  int sum = 0;
  for (int i=0; i<n; i++) {
    sum += msg[i]*pub_set[i];
  }
  return sum;
}

void decrypt(int a[], int x[], int sum, int n, int p, int b) {
  // a is the superincreasing set, x is the output vector,
  // sum is the subset sum (encrypted message), n is length(a),
  // p is the chosen prime, b is the chosen private multiplier
  int truesum = positive_mod(sum*inverse_modp(b, p), p);
  superinc_ss(a, x, n, truesum);
}

void print_array(int x[], int n) {
  // helper debug code
  for (int i=0; i<n; i++) {
    printf("%d\n", x[i]);
  }
}

int main() {
  // some demo code to show off
  int msg[8] = {1, 0, 1, 1, 0, 1, 0, 0};
  int sup_set[8] = {3, 5, 9, 23, 51, 132, 263, 498};
  int pub_set[8] = {3, 5, 9, 23, 51, 132, 263, 498};
  int n = 8;
  int p = 1097;
  int b = 249;
  sis_to_pub(pub_set, n, p, b);
  // now pub_set is the public set to be used for encryption
  print_array(pub_set, n);
  int enc = encrypt(pub_set, msg, n);
  printf("sum: %d\n", enc);
  // to decrypt:
  int dec[8];
  decrypt(sup_set, dec, enc, n, p, b);
  print_array(dec, n);
  return 0;
}
