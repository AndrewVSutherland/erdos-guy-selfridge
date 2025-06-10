#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <memory.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <primecount.h>
#include <primesieve.h>
#include <primesieve/iterator.h>

#define MAXN    ((1L<<48)-1)

/*
     MCUTOFF is a tunable parameter to choose cutoff between enumerating/counting primes
     We will seive p until p > t^(1-
*/

static double MCUTOFF = 0.225;      // tunable parameter to control prime enumeration/counting cutff
                                    // should be around 0.25, but primeisieve is very fast so a bit lower is better

/*
    The tables P,PI,F,M are independent of N < 2^40 and computed at startup (about 0.5s)
    and are read-only after that.  They use about 300MB memory (this is most of what we will use)
*/

static int32_t *P;                  // P[n] is the nth prime for n up to MAXPI, and we put P[0] = 1
static int32_t *PI;                 // PI[n] = pi(n) for n <= MAXP (in particular, P[PI[p]]=p for primes p)
static int32_t MAXP = 310248233;    // we require pi(MAXP) <= MAXPI
static int32_t MAXPI = ((1<<24)-1);

// Factorizations of cofactors m are zero-terminated lists of pp's
// We only consider m that are MAXP-smooth, so pi fits in 16-bits (could extend to 24 bits and make e 8bits)
static struct pp {
    unsigned pi : 24;               // index into P
    unsigned e : 8;
} *F,*Fend;                         // concatenation of zero-terminated factorizations in descending order by pi
static uint32_t *M;                 // F[M[m]] holds the factorization of m <= MAXM for MAXP-smooth m, M[m]=0 ow
static uint32_t MAXM = 0x7FFFFFFF;  // Largest m for which M[m] is valid (so M has length (MAXM+1)), set at startup time.


static inline double get_time (void) // accurate to at least 10ms
    { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (double)t.tv_sec + (double)t.tv_nsec / 1000000000.0L; }

static inline int64_t cdiv (int64_t a, int64_t b)     // ceil(a/b)
    { return (a+b-1) / b; }

static inline int v2 (uint32_t x) { return __builtin_ctz(x); }

static inline int64_t pi (int64_t n)
    { return n <= MAXP ? PI[n] : primecount_pi(n); }

static inline int64_t max (int64_t a, int64_t b)
    { return a<b?b:a; }

static inline int64_t min (int64_t a, int64_t b)
    { return a>b?b:a; }

static inline primesieve_iterator primesieve_start(int64_t minp, int64_t maxp)
{
    primesieve_iterator ctx;
    primesieve_init (&ctx);
    primesieve_jump_to (&ctx, minp, maxp);
    return ctx;
}

static inline void primesieve_stop (primesieve_iterator *ctx)
    { primesieve_free_iterator (ctx); }

void setup (int64_t maxp, int64_t maxm)
{
    assert (maxp <= MAXP);
    MAXP = maxp;
    int64_t maxpi = primecount_pi(MAXP);
    assert (maxpi <= MAXPI);
    MAXPI = maxpi;
    if ( !(maxm&1) ) maxm += 1; // make maxm even for convenience
    assert (maxm <= MAXM);
    assert (maxm < (uint64_t)MAXP*MAXP);
    MAXM = maxm;
    // compute the list of primes p <= MAXP, storing the nth prime in P[n] (we set P[0]=1 for convenience)
    P = malloc((MAXPI+1)*sizeof(*P)); assert(P);  P[0] = 1;
    PI = calloc((MAXP+1),sizeof(*PI)); assert (PI);
    primesieve_iterator ctx = primesieve_start (0,MAXP);
    for ( int32_t p = 0, n = 1 ; (p = primesieve_next_prime(&ctx)) <= MAXP ; ) { P[n] = p; PI[p] = n++; }
    primesieve_stop(&ctx);

    // set PI[n] to pi(n) for all n (using values already set for n prime)
    for ( int32_t n = 1 ; n <= MAXP ; n++ ) if ( !PI[n] ) PI[n] = PI[n-1];

    // set M[m]=n where p_n is the largest odd prime divisor of M up to MAXP (and zero if no such p_n exists)
    // we don't use a windowed sieve here as we cap maxm at a modest value for the standard greedy algorithm
    // and for the fast variant we are only sieving out to N^(5/8)
    M = calloc(maxm+1,sizeof(*M));
    for ( int32_t pi = 1 ; pi <= MAXPI ; pi++ ) for ( int32_t p=P[pi], q = p ; q <= maxm ; q+=p ) M[q] = pi;

    // Compute F and update M so that F[M[m]] holds the factorization of m for all MAXP-smooth m <= MAXM
    // Factorizations are zero terminated lists of pp in reverse order by prime

    // This step could be sped up signficantly (e.g by avoiding divisions), but we expect to reuse this
    // data for many values of N and t so that the amortized cost is negligible and for the fast variant
    // it is negligible even for a single (N,t), since maxm = O(N^(2/3))
    int64_t Fsize = max(4*maxm,1<<10);
    F = malloc(Fsize*sizeof(*F));
    struct pp *f = F; (f++)->pi = 0;                // skip offset 0
    for ( int64_t m = maxm ; m > 1 ; m-= 2 ) {      // handle odd m first
        struct pp *g = f;
        g->pi = M[m]; g->e = 0;
        int64_t q = m;
        for ( ; M[q] ; q /= P[M[q]] ) {
            if ( M[q] == g->pi ) { g->e++; }
            else { (++g)->pi = M[q]; g->e = 1; }
        }
        if (q!=1) { M[m] = 0; continue; }
        g++; (g++)->pi = 0;
        M[m] = f-F; f = g;
        if ( f-F > Fsize-16 ) { int64_t d = f-F; Fsize = (5*Fsize)/4; F = realloc(F,Fsize*sizeof(*F)); f = F+d; }
    }

    M[1] = f-F; (f++)->pi = 0;
    for ( int64_t m = maxm-1 ; m > 1 ; m-= 2 ) {    // now handle even m
        int64_t e = v2(m), q = m>>e;
        if ( !M[q] ) { M[m] = 0; continue; }
        struct pp *g = F+M[q];
        M[m] = f-F;
        while (g->pi) *f++ = *g++;
        f->pi = 1; (f++)->e = e; (f++)->pi = 0;
        if ( f-F > Fsize-16 ) { int64_t d = f-F; Fsize = (5*Fsize)/4; F = realloc(F,Fsize*sizeof(*F)); f = F+d; }
    }

    Fsize = f-F;
    assert (Fsize < (1L<<32));
    F = realloc (F, Fsize*sizeof(*F));
    Fend = F+Fsize;
}

static void unsetup (void)
{
    free(F); free(M); free(P); free(PI);
}

// computes min(e,v_m(P^E)) where m has factorization f and v_m(n) is largest r for which m^r|n
static inline int64_t fcnt (int64_t *E, int64_t e, struct pp *f) 
{
    for ( ; f->pi ; f++ ) e = min(e,E[f->pi]/f->e);
    return e;
}

// computes v_{pm}(P^E) where p = P[i], m has factorization f, and v_{pm}(n) is largest r for which (pm)^r|n
// note that p|m is allowed here
static inline int64_t fcnti (int64_t *E, int32_t i, struct pp *f) 
{
    assert (f >= F && f < Fend);
    assert (f->pi <= i);
    int64_t e = i > f->pi ? E[i] : E[i]/(f->e+1);
    for ( ; f->pi ; f++ ) e = min(e,E[f->pi]/f->e);
    assert (f <= Fend);
    return e;
}

static int32_t fac_s (int64_t t)
{
    int64_t s = sqrt(t); assert(s*(s-1) < t); while ( s*(s-1) < t ) s++;
    assert (s <= MAXP);
    return s;
}

struct fac_item { int64_t n; int64_t m; struct pp f[12]; int64_t p; int64_t q; int64_t c; };

struct fac {
    struct fac_item *L;
    int32_t size;
    int32_t cnt;
    int64_t N;
    int64_t t;
};

static struct fac *fac_alloc (int64_t N, int64_t t)
{
    struct fac *v = malloc(sizeof(*v));
    v->size = 1024; v->cnt = 0; v->N = N; v->t = t;
    v->L = malloc(v->size*sizeof(*v->L));
    return v;
}

static void fac_free (struct fac *v)
{
    free (v->L);
    free (v);
}

static inline void fac_extend (struct fac *v, int64_t n, int64_t m, struct pp *f, int64_t p, int64_t q, int64_t c, int verbosity)
{
    assert(m*(p+1) >= v->t);
    if ( verbosity > 3 ) {
        if ( p+1==q) fprintf(stderr,"factor: (%ld*%ld)^%ld\n",m,q,n);
        else fprintf(stderr,"factor: (%ld*p)^%ld for p in (%ld,%ld]\n",m,n,p,q);
    }
    if (v->cnt >= v->size ) { v->size *= 2; v->L = realloc(v->L,v->size*sizeof(*v->L)); }
    struct fac_item *r = v->L + v->cnt;
    int i = 0; while (f[i++].pi);
    r->n = n;
    r->m = m;
    memcpy(r->f,f,i*sizeof(*f));
    r->p = p;
    r->q = q;
    r->c = c;
    v->cnt++;
}

static inline void fac_extend_mp (struct fac *v, int64_t n, int64_t m, struct pp *f, int64_t p, int verbosity)
    { fac_extend(v,n,m,f,p-1,p,1,verbosity); }

static inline void fac_extend_mp2 (struct fac *v, int64_t n, int64_t m, struct pp *f, int i, int verbosity)
{
    assert (i>f->pi && i <= MAXPI);
    struct pp h[16] = {{i, 1}}, *x = h+1;
    for ( ; f->pi ; *x++ = *f++ );
    x->pi = 0;
    fac_extend_mp(v,n,m*P[i],h,P[i],verbosity);
}

static inline void fac_extend_m (struct fac *v, int64_t n, int64_t m, struct pp *f, int verbosity)
{
    int64_t p = P[f->pi];
    f->e--;
    fac_extend_mp(v,n,m/p,f->e?f:f+1,p,verbosity);
    f->e++;
}

static int fac_verify (struct fac *v, int verbose)
{
    int64_t N = v->N, sqrtN = (int64_t)sqrt(N), t = v->t;
    int32_t s = fac_s(t);
    int32_t maxpi = PI[s-1], maxp = P[maxpi];
    int64_t *E = calloc(maxpi+1,sizeof(*E));
    for ( int32_t i = 1 ; i <= maxpi ; i++ ) for ( int64_t q = P[i] ; q <= N ; q *= P[i] ) E[i] += N/q;
    int64_t cnt = 0, lastp = 0, lastpi = 0, nextpi = 0;
    for ( struct fac_item *r = v->L ; r < v->L + v->cnt ; r++ ) {
        assert (r->n && r->p < r->q && r->q <= N);
        assert (r->m*(r->p+1) >= t);
        if ( r->q <= maxp ) {
            int32_t x = 0;
            for ( int32_t pi = PI[r->p]+1 ; pi <= PI[r->q] ; pi++ ) { E[pi] -= r->n; x += r->n; }
            for ( struct pp *f = r->f ; f->pi ; f++ ) E[f->pi] -= x*f->e;
            cnt += x;
        } else {
            lastpi = ( r->p == lastp && nextpi ) ? nextpi : pi(r->p);
            nextpi = pi(r->q);
            assert (r->p+1 > maxp);
            if ( r->q <= sqrtN ) {
                assert (N/r->q + N/(r->q*r->q) == r->n);
                assert (N/(r->p+1) + N/((r->p+1)*(r->p+1)) == r->n);
            } else {
                assert (N/(r->p+1) == r->n && N/r->q == r->n);
            }
            int64_t x = r->n*(nextpi-lastpi);
            for ( struct pp *f = F+M[r->m] ; f->pi ; f++ ) E[f->pi] -= x*f->e;
            cnt += x; 
          
        }
        lastp = r->q;
    }
    for ( int32_t pi = 1 ; pi <= maxpi ; pi++ ) assert (E[pi] >= 0);
    if ( verbose >= 0 ) fprintf(stderr,"Verified factorization of %ld! into %ld factors >= %ld\n", N, cnt, t);
    free (E);
    return 1;
}

static int fac_dump (struct fac *v, char *filename)
{
    FILE *fp = fopen(filename, "w"); assert(fp);
    for ( struct fac_item *r = v->L ; r < v->L + v->cnt ; r++ ) fprintf(fp,"%ld,%ld,%ld,%ld\n",r->n,r->m,r->p,r->q);
    fclose(fp);
    return 1;
}

int64_t tfac (int64_t N, int64_t t, int fast, int feasible, int verbosity, int verify, char *dumpfile)
{
    if ( verbosity > 1 ) fprintf(stderr,"tfac(%ld,%ld) %s%s%s\n",N,t,fast?"fast":"greedy",feasible?" feasibility test":"",verify?" verification on":"");

    double start = get_time();
    assert (N >= 10 && N < MAXN && 4*t > N && 2*t < N);
    int32_t sqrtN = (int32_t)sqrt(N);
    int32_t s = fac_s(t);
    int32_t maxpi = PI[s-1];
    
    // Compute p-adic valuations of N! for p < s by setting E[pi(p)] = v_p(N!)
    // We could cache and resue this but it takes negligible time to recompute it (under 1ms for N <= 10^12)
    int64_t *E = calloc(maxpi+1,sizeof(*E));
    for ( int32_t i = 1 ; i <= maxpi ; i++ ) for ( int64_t q = P[i] ; q <= N ; q *= P[i] ) E[i] += N/q;

    // Compute upper bound maxm on the smooth cofactors m we will ever use.
    // For the standard greedy algorithm these could be as large as t-1
    // For the fast variant we set maxm = t^(5/8), which works better in practice than t^(2/3)
    if (!fast) assert(t<=MAXM+1);
    int64_t maxm = fast ? pow(t,0.625) : t-1;
    assert (maxm <= MAXM);
    uint32_t *Ms = malloc((maxm+1)*sizeof(*Ms)); Ms[0] = 0;
    for ( uint32_t m = 1, *p = Ms+1 ; m < s ; m++ ) *p++ = m;
    uint32_t numm = s;
    // For the standard greedy algorithm essentially all of the post-setup time is spent here (if fast is set, this step takes negligible time)
    for ( uint32_t m = s ; m <= maxm ; m++ ) if ( F[M[m]].pi && F[M[m]].pi <= (fast ? PI[t/m] : maxpi) ) Ms[numm++] = m;
    Ms = realloc(Ms, numm*sizeof(*Ms));
    numm--; maxm = Ms[numm];

    struct fac *v = verify ? fac_alloc(N,t) : 0;

    /*
      We begin by constructing factors m*p >= t with m minimal for all p >= s >= sqrt(t).
      In this range we can always make m = cdiv(t,p) optimal for p.
      For each prime p in [s,N] we will have n=v_p(N!) identical factors m*p.
      We don't really care about p, we just need the values of m and n.
      We will treat m = ceil(t/s),...,1 in descending order.  When m is large, p will
      be small, and the value of n will be changing rapidly, so we just enumerate primes p
      and compute update m and n as we go.  But when m is small it is much more
      efficient to determine the exact range of p applicable to each (m,n) pair and
      then count the primes in this range, rather than enumerate them (we can us a
      precomputed table of pi(x) values to do this very quickly).

      The crossover point between "large" and "small" m is a performance parameter depending
      on the relative speed of primesieve vs primecount; the choice won't impact the results.
      t^(1/5) seems like a reasonable heuristic value for large t
    */

    int64_t m = cdiv(t,s);                                  // largest possible m for p >= s
    assert (m <= maxm && Ms[m] == m);

    // mid determines the crossover point between sieving primes and counting them
    // the asymptotically correct choice for the implementation of pi(x) we are using is t^(1/6)
    // but we use a slightly larger t^(1/5) here to better balance the time spent sieving/counting primes
    int64_t mid = min((int64_t)pow(t,MCUTOFF),(t-1)/sqrtN); // we will enumerate p up to (t-1)/mid, then switch to counting primes
    if ( (int64_t)sqrtN*mid >= t ) mid = (t-1)/sqrtN;       // force p > sqrtN for m < mid (for convenience only)

    if ( verbosity > 2 ) fprintf(stderr,"N=%ld, t=%ld, sqrt(N)=%d, s=%d, maxpi=%d, maxm=%ld, numm=%d, mid=%ld (%.6fs)\n", N, t, sqrtN, s, maxpi, maxm, numm, mid, get_time()-start);

    primesieve_iterator ctx = primesieve_start(s,(t-1)/mid);
    int64_t p; // p will fit in 32 bits but we want to mults at 64-bits
    int64_t cnt = 0; // running count of the number of factors (goal is to get this >= N)

    // handle primes in [s,sqrt(N)] (here we are happy to recompute n for each p, this phase takes no time)
    while ( (p = primesieve_next_prime(&ctx)) <= sqrtN) {
        while ( (m-1)*p >= t ) m--;     // minimize m with m*p >= t
        int64_t n = N / p + N / (p*p);  // compute n = v_p(N!) (expensive but negligible)
        // update valuations in E to reflect our use of n cofactors m
        for ( struct pp *f = F+M[m] ; f->pi ; f++ ) E[f->pi] -= n*f->e;
        cnt += n;
        if ( v ) fac_extend_mp (v,n,m,F+M[m],p,verbosity);
    }

    if ( verbosity > 2 ) fprintf(stderr,"cnt=%ld for p in [s,sqrt(N)], m=%ld (%.6fs)\n", cnt, m, get_time()-start);

    int64_t pmmax = (t-1) / (m-1);  // largest p for this m
    assert (p>pmmax || m==cdiv(t,p));
    int64_t n = N / (sqrtN+1);      // smallest n for p > sqrtN
    int64_t pnmax = N / n;          // largest p for this n
    int64_t plmmax = (t-1) / mid;   // largest p for m > mid

    // handle primes in (sqrtN,plmmax] with large m using primesieve (this should take about half the time if mid is optimal)
    int64_t pmin = p-1;
    for ( ; p <= plmmax ;) {
        while ( p > pmmax ) { m--; pmmax = (t-1)/(m-1); }   // update m
        while ( p > pnmax ) { n--; pnmax = N/n; }           // update n
        int64_t pmax = min(pmmax,pnmax); assert (p <= pmax);
        int64_t c = 1;
        while ( (p = primesieve_next_prime(&ctx)) <= pmax ) c++;
        for ( struct pp *q = F+M[m] ; q->pi ; q++ ) E[q->pi] -= c*n*q->e;
        cnt += c*n;
        if ( v ) { fac_extend (v,n,m,F+M[m],pmin,pmax,c,verbosity); pmin = p-1; }
    }
    primesieve_stop(&ctx);
    int64_t lastpi = pi(plmmax), nextpi;
    pmin = plmmax;
    if ( verbosity > 2 ) fprintf(stderr,"cnt=%ld for %ld p >= s with m < mid (%.6fs)\n", cnt, lastpi-maxpi, get_time()-start);

    // handle primes in (plmax,t] with small m in [mid,2] using primecount (this should take about half the time if mid is optimal)
    // here we iterate over m rather than p
    for ( m = mid ; m > 1 ; m-- ) {
        int64_t p = cdiv(t,m), pmax = (t-1)/(m-1);       // p in (pmin,pmax] are the p for this m
        n = N/p; pnmax = min(N/n,pmax);                  // p in (pmin,pnmax] are the p for this n
        while ( pmin < pmax ) {
            nextpi = pi(pnmax);
            int64_t c = (nextpi-lastpi); cnt += c*n;     // number of factors for this m and n
            for ( struct pp *q = F+M[m] ; q->pi ; q++ ) E[q->pi] -= c*n*q->e;
            if ( v ) fac_extend (v,n,m,F+M[m],pmin,pnmax,c,verbosity);
            pmin = pnmax; n--; pnmax = min(N/n,pmax);    // proceed to the next n
            lastpi = nextpi;
        }
    }
    if ( verbosity > 2 ) fprintf(stderr,"cnt=%ld for %ld p in [s,t) (%.6fs)\n", cnt, lastpi-maxpi, get_time()-start);
    assert (lastpi == pi(t-1));

    // Finally, handle primes p  in [t,N]
    if ( 3*t <= N ) {
        nextpi = pi(N/3);
        cnt += 3*(nextpi-lastpi);
        if ( v ) fac_extend (v,3,1,F+M[1],t-1,N/3,nextpi-lastpi,verbosity);
        lastpi = nextpi;
    }
    nextpi = pi(N/2);
    cnt += 2*(nextpi-lastpi);
    if ( v ) fac_extend (v,2,1,F+M[1],max(t-1,N/3),N/2,nextpi-lastpi,verbosity);
    lastpi = nextpi; nextpi = pi(N);
    cnt += nextpi-lastpi;
    if ( v ) fac_extend (v,1,1,F+M[1],N/2,N,nextpi-lastpi,verbosity);
    if ( verbosity > 2 ) fprintf(stderr,"cnt=%ld for %ld p in [s,N] (%.6fs)\n", cnt, nextpi-maxpi, get_time()-start);

    // Verify our assumption that we can use the optimal m for all p >= s (this will be verified again at the end but
    // the check is cheap so we do it now before a possible feasibility check).
    for ( int32_t i = 1 ; i <= maxpi ; i++ ) assert(E[i] >= 0);

    if ( feasible ) {
        long double ebits = 0, epsilon = 0.0000000000000001L; // 10^{-15} < 2^52
        for ( int32_t i = 1 ; i <= maxpi ; i++ ) ebits += E[i]*log(P[i]+epsilon); // make sure we get an upper bound
        return cnt + floorl(ebits/log(t-epsilon));
    }

    /*
        At this point we have dealt with all prime factors p >= s of N! using the optimal m = cdiv(t,p)
        for each p and updated E so that

            P^E := prod_{1<=i<=pimax} P[i]^E[i]

        is the divisor of N! that we still need to factor.

        We now process primes p=P[i] in descending order with cofactors m=M[j] in ascending order.
        We need p*m >= t, and m must be p-smooth (because there are no primes > p left).
    
        Up to this point the standard and fast greedy algorithms are identical but now they diverge
    */
    if ( !fast ) { // original greedy (simpler but slower)
        int64_t pcnt = 0;
        for ( int i = 0 ; i <= maxpi ; i++ ) pcnt += E[i];
        for ( int32_t i = maxpi, j = cdiv(t,s) ; i ; ) {
            // Update j so that all m' >= m=Ms[j] are valid for use with p=P[i]
            while ( j <= numm && ((int64_t)P[i]*Ms[j] < t || F[M[Ms[j]]].pi > i) ) j++;
            struct pp *f = F+M[Ms[j]];
            if ( j > numm ) break;
            int64_t e = fcnti(E,i,f);
            if ( !e ) {
                if ( pcnt < 40 ) { // if there are less than 40 primes left (with multiplicity), check if we are done
                    int64_t q = 1;
                    // P[ii] <= MAXP and t < N <= MAXN guarantees q will not overflow sicne MAXP*MAXN < 2^60
                    for ( int ii = i ; ii && q < t ; ii-- ) for ( int x = 0 ; x < E[ii] && q < t ; x++ ) q *= P[ii];
                    if ( q < t ) break;
                }
                j++; continue;
            }
            cnt += e; for ( E[i] -= e, pcnt -= e ; f->pi ; f++ ) { E[f->pi] -= e*f->e; pcnt -= e*f->e; }
            if ( v ) fac_extend_mp (v,e,Ms[j],f,P[i],verbosity);
            while ( i && !E[i] ) i--;
        }
    } else {    // fast greedy
        /*
            In the fast greedy algorithm we require m to be (p-1)-smooth.
            This ensures that the exponents in E associated to p and m are disjoint, and if
            m >= cdiv(t,p) is (p-1)-smooth, so is every candidate m' > m (by construction).

            Our cofactors are m bounded by t^(5/8), which means that as p gets smaller we
            cannot achieve p*m >= t and will need to start using products of primes.

            We proceed in two phases, in the first we handle p for which maxm*p > t where
            we might be able to use p*m as a factor.
        */
        int32_t pimin = pi(cdiv(t,maxm))+1;
        for ( int32_t i = maxpi, j = cdiv(t,s) ; i >= pimin ; i-- ) {
            // Update j so that all m' >= m=Ms[j] are valid for use with p=P[i]
            while ( (int64_t)P[i]*Ms[j] < t || F[M[Ms[j]]].pi >= i ) j++;
            struct pp *f = F+M[Ms[j]];
            int64_t e = fcnt(E,E[i],f);    
            if ( e < E[i] ) { // if we cannot completely remove p from P^e using m, try p^2 and m=cdiv(t,p^2)
                int64_t m = cdiv(t,(int64_t)P[i]*P[i]);
                struct pp *g = F+M[m];
                e = fcnt(E,E[i]/2,g);
                if ( e ) {
                    if ( v ) fac_extend_mp2 (v,e,m,g,i,verbosity);
                    cnt += e; for ( E[i] -= 2*e ; g->pi ; g++) E[g->pi] -= e*g->e;
                }
                e = fcnt(E,E[i],f); // recompute e (may have changed)
            }
            if ( e ) {
                if ( v ) fac_extend_mp (v,e,Ms[j],f,P[i],verbosity);
                cnt += e; for ( E[i] -= e ; f->pi ; f++ ) E[f->pi] -= e*f->e;
            }
            if ( E[i] ) { // p still divides P^E?
                // try a larger m
                e = 0;
                for ( int32_t k = j+1 ; k <= numm ; k++ ) {
                    struct pp *g = F+M[Ms[k]];
                    int64_t x = fcnt(E,E[i],g);
                    if ( x > e ) {  e = x; f = g; m = Ms[k]; if ( e == E[i] ) break; }
                }
                if ( e ) {
                    if ( v ) fac_extend_mp (v,e,m,f,P[i],verbosity);
                    cnt += e; for ( E[i] -= e ; f->pi ; f++ ) E[f->pi] -= e*f->e;
                }
                if ( E[i] ) { // p still divides P^E?
                    // try a larger m for p^2, here we assume m=cdiv(t,p^2) < s so that Ms[m] = m
                    e = 0; f = 0; int64_t m = cdiv(t,(int64_t)P[i]*P[i])+1; assert (Ms[m]==m);
                    for ( int32_t k = m ; k <= numm ; k++ ) {
                        struct pp *g = F+M[Ms[k]];
                        int64_t x = fcnt(E,E[i]/2,g);
                        if ( x > e ) {  e = x; f=g; m = Ms[k]; if ( e == E[i] ) break; }
                    }
                    if ( e ) {
                        if ( v ) fac_extend_mp2 (v,e,m,f,i,verbosity);
                        cnt += e; for ( E[i] -= 2*e ; f->pi ; f++ ) E[f->pi] -= e*f->e;
                    }
                    // we may still have E[i] > 0 but we will usually have E[i] <= 1
                }
            }
        }

        if ( verbosity > 2 ) fprintf(stderr,"cnt=%ld after initial pass of p in (cdiv(t,maxm),s) (%.6fs)\n", cnt, get_time()-start);
        while ( maxpi && !E[maxpi] ) maxpi--;

        // Now we just want to use up whatever is left the best we can
        // This should consist almost entirely of primes < t^3/8 and takes very little time
        int64_t good = 5*cdiv(t,4);  // we will settle for factors in [t,good]
        struct pp c[16]; c->pi = 0; c->e = 0; // factorization of product of primes we are assembling (which may exceed maxm)
        while (maxpi) {
            while ( maxpi && !E[maxpi] ) maxpi--;
            if ( !maxpi ) break;
            int32_t i = maxpi;
            int64_t q = P[i];
            struct pp *f=c; f->pi = i; f->e = 1; (++f)->pi = 0;
            E[i]--; // update E as we go, we will undo the update if we get stuck
            while ( i && !E[i] ) i--;
            if ( !i ) break;
            while ( i && q*P[i] < good ) {
                q *= P[i]; E[i]--;
                if ( (f-1)->pi == i ) { (f-1)->e++; } else { f->pi = i; f->e = 1; (++f)->pi = 0; }
                while ( i && !E[i] ) i--;
            }
            if ( !i && q < t ) break;
            int64_t e = 1 + fcnt(E,E[c->pi]/c->e,c+1);    // 1+v_q(P^E) (we aleady removed one factor of q)
            if ( q < t ) {
                assert ( q > s);
                int64_t b = 0; struct pp *g = 0;
                // first look for a cofactor smaller than the smallest prime divisor of q
                for ( m = cdiv(t,q) ; m < P[(f-1)->pi] ; m++ ) {
                    int64_t x = fcnt(E,e,F+M[m]);
                    if ( x > b ) { b = x; g = F+M[m]; }
                    if ( x == e ) break;
                }
                if ( b ) {
                    // incorporate factorization of m into c
                    while ( g->pi ) { E[g->pi] -= g->e; *f++ = *g++; }
                    f->pi = 0;
                    q *= m;
                } else {
                    if (!i) break;
                    q *= P[i]; E[i]--;
                    if ( (f-1)->pi == i ) { (f-1)->e++; } else { f->pi = i; f->e = 1; (++f)->pi = 0; }
                    b = 1 + fcnt(E,E[c->pi]/c->e,c+1);
                    assert(b);
                }
                assert (q >= t);
                e = b;
            }
            if ( v ) fac_extend_m (v,e,q,c,verbosity);
            cnt += e--;
            for ( f = c ; f->pi ; f++ ) E[f->pi] -= e*f->e;
            c->pi = 0;
            maxpi = i;
        }
        for ( struct pp *f = c ; f->pi ; f++ ) E[f->pi] += f->e; // restore any patial factor we did not remove so we can report/check remainder
    }

    while ( maxpi && !E[maxpi] ) maxpi--;
    free(Ms);

    // MAXP * MAXN < 2^60 guarantees that q will not overflow in the line below
    int64_t q = 1; for ( int32_t i = 1 ; i <= maxpi ; i++ ) { assert(E[i] >= 0); for ( int32_t e = 0 ; e < E[i] ; e++ ) { q *= P[i]; assert (q<t); } }
    if ( verify ) fac_verify(v,verbosity);
    if ( verbosity > 1 ) fprintf(stderr,"%ld factors >= %ld with remainder %ld (%.6fs)\n", cnt, t, q, get_time()-start);
    if ( v && dumpfile ) {
        fprintf(stderr,"Dumping factorization certificate to %s ...\n", dumpfile);
        fac_dump(v,dumpfile);
    }
    if ( v ) fac_free(v);

    free (E);
    return cnt;
}


// returns a value of t >= aN/b that yields a good lower bound on t(N), or 0 if no such t can be found
int64_t tbound (int64_t N, int a, int b, int fast, int exhaustive, int verbosity, int verify)
{
    assert (a*5 <= 2*b && a*4 >= b);
    int64_t t = cdiv(a*N,b);
    int64_t cnt = tfac(N,t,fast,0,verbosity,verify,0);
    while ( cnt < N ) cnt = tfac(N,--t,fast,0,verbosity,verify,0);
    int64_t tmin = t, tmax = (2*N)/5;

    /*
        We use a modified bisection search of [tmin,tmax) for t with tfac(N,t) >= N but tfac(N,t+1) < N that
        uses the excess/deficit (tfac(N,t)-N) to choose the next bisection point.
    */
    while ( tmin < tmax-1 ) {
        if ( cnt >= N ) tmin = max(t,tmin); else tmax = min(t,tmax);
        if ( verbosity > 1 ) fprintf (stderr,"t=%ld gave %ld extra factors, new t-range is [%ld,%ld)\n", t, cnt-N, tmin, tmax);
        t = round(exp(log(t)+(cnt-N)*log(t)/N));
        if ( t <= tmin ) t = max((3*tmin+tmax)/4,tmin+1);
        if ( t >= tmax ) t = min((tmin+3*tmax)/4,tmax-1);
        cnt = tfac(N,t,fast,0,verbosity,verify,0);
    }
    assert (tmax < (2*N)/5);
    if ( ! exhaustive ) return tmin;
    if ( verbosity > 0 ) fprintf(stderr,"t(%ld) >= %ld proved\n", N,tmin);

    /* Now use a binary search to get an upper bound on the best possiible t that tfac(N,t) could return */
    int64_t low = tmin;
    int64_t high = (2*N)/5;
    cnt = tfac(N,high,fast,1,verbosity,0,0);
    assert (cnt < N);
    while ( low < high-1 ) {
        int64_t mid = (low+high)/2;
        cnt = tfac(N,mid,fast,1,verbosity,0,0);
        if ( cnt < N ) { high = mid; tmax = mid; } else { low = mid; }
    }
    assert (tmax > tmin);
    if ( verbosity > 0 ) fprintf(stderr,"t(%ld) >= %ld cannot be proved\n",N,tmax);
    int threads = omp_get_max_threads();
    if ( verbosity > 0 ) fprintf(stderr,"checking %ld values of t in (%ld,%ld) using %d threads\n",tmax-tmin-1,tmin,tmax,threads);
    #pragma omp parallel num_threads(threads)
    {
    int32_t tid = omp_get_thread_num();
    for ( int64_t t = tmin+1 ; t < tmax ; t++ ) {
        if ( (t%threads) != tid ) continue;
        if ( tfac(N,t,fast,0,verbosity,0,0) >= N ) {
            #pragma omp critical(foo)
            {
                tmin = max(tmin,t);
                if ( verbosity >= 0 ) fprintf(stderr,"\rt(%ld) >= %ld proved\r", N,tmin);
            }
        }
    }
    }
    if ( verify ) tfac(N,tmin,fast,0,verbosity,verify,0);
    return tmin;
}

static void usage (void)
{
    fprintf(stderr,
        "Usage: egs [-v level] [-h filename] [-d filename] [-r] [-c] [-e] [-f] N-range [t or t/N ratio]\n"
        "       -v level      integer verbosity level -1 to 4 (optional, default is 0)\n"
        "       -h filename   hint-file with records N:t (required if range of N is specified)\n"
        "       -d filename   output-file to dump factorization to (one factor per line, only valid if t is specified)\n"
        "       -r            verify factorization (set automatically if dump is specified)\n"
        "       -c            create hint-file rather than reading it (must be specified in combination with -h)\n"
        "       -e            use the best t for which the algorithm can prove t(N) >= t (optional)\n"
        "       -f            use fast version of greedy algorithm\n"
        "       -m            exponent for primecount/primesieve cutuff, must lie in [1/6,1/3]\n"
        "       N-range       integer N or range of integers minN-maxN (required, scientific notation supported)\n"
        "       t             integer t to use for single N (optional, a good t will be determined if unspecified)\n"
        "       t/N ratio     a/b with integers a,b>0 specifying t = ceil(aN/b), set to 1/3 if unspecified\n");

}

int main (int argc, char *argv[])
{
    if ( argc < 2 ) { usage(); return 0; }

    int verbosity=0, exhaustive=0, create=0, fast=0, verify=0;
    char *hintfile = 0, *dumpfile = 0;
    int64_t minN=0, maxN=0, t=0;
    int a=1, b=3;

    for ( int i = 1 ; i < argc ; i++ ) {
        char *s = argv[i];
        if ( maxN > minN && t ) { fprintf(stderr, "ignoring extraneous argument %s\n",s); continue; }
        if ( *s == '-' ) {
            switch(*(s+1)) {
            case 'v': { if ( i+1 >= argc) { usage(); return -1; } verbosity = atoi(argv[i+1]); i++; break; }
            case 'h': { if ( i+1 >= argc || hintfile ) { usage(); return -1; } hintfile = argv[i+1]; assert(hintfile[0] != '-'); i++; break; }
            case 'd': { if ( i+1 >= argc || dumpfile ) { usage(); return -1; } dumpfile = argv[i+1]; assert(dumpfile[0] != '-'); i++; verify=1; break; }
            case 'r': { verify = 1; break; }
            case 'c': { create = 1; break; }
            case 'm': { double x = atof(argv[i+1]); assert (x >= 0.2 && x <= 0.3); MCUTOFF = x; i++; break; }
            case 'e': { exhaustive = 1; break; }
            case 'f': { fast = 1; break; }
            default: { printf("unrecognized option %s\n", s); usage(); return -1; }
            }
        } else {
            if ( !minN ) {
                if (*s=='[')s++;
                double x = strtod(s,&s);
                if ( (x-round(x)) > 0.0001 ) { fprintf(stderr, "N=%f must be an integer.\n",x); usage(); return -1; }
                minN = round(x);
                if ( s && *s ) {
                    if ( *s == '.' ) { while (*(++s)=='.'); }
                    else if ( *s == '-' || *s == ',' ) { s++; }
                    else { puts(s); usage(); return -1; }
                    maxN = strtold(s,0);
                    assert(maxN >= minN);
                } else {
                    maxN = minN;
                }
            } else {
                char *x;
                if ( (x=strchr(s,'/')) ) {
                    a = atoi(s); b = atoi(x+1);
                    assert (a > 0 && b > 0 && 4*a >= b && 5*a <= 2*b);
                } else {
                    if ( maxN > minN && t ) { fprintf(stderr, "For a range of N you need to specify the t/N ratio (e.g. 1/3) not a fixed value of t\n"); return -1; }
                    t = strtold(s,0);
                }
            }
        }
    }
    if ( minN < 14 || maxN > MAXN ) { fprintf(stderr,"N-range [%ld,%ld] must be contained in [14,%ld)\n", minN, maxN, MAXN); return -1; }
    if ( t && 4*t <= maxN ) { fprintf(stderr,"t=%ld must be greater than N/4\n", t); return -1; }

    double start = get_time();
    int64_t maxt = 2*maxN/5;
    int64_t maxp = fac_s(maxt);
    // Set the global upper bound on maxm for this run, sufficient to handle any N <= maxN and t <= maxt
    // For the fast variant we use t^(5/8) rather than t^(2/3) because this is faster in practice (possibly due to caching) and doesn't really change the quality of the bounds
    int64_t maxm = fast ? pow(maxt,5.0/8) : maxt-1;
    if ( maxp > MAXP || maxm > MAXM ) {
        assert (!fast);
        fprintf(stderr,"N=%ld is too large for for this implementation of the standard greedy algorithm.  Use the -f option to switch to fast variant.\n", maxN); return -1;
    }

    setup(maxp,maxm);
    if ( verbosity > 0 ) fprintf(stderr,"Computed %d-smooth factorizations of m <= %d using %.3fMB of memory (%.3fs)\n", MAXP,MAXM,4.0*(MAXM+(Fend-F))/(1<<20),get_time()-start);

    char rbuf[32];
    if ( a == 1 ) sprintf(rbuf,"ceil(N/%d)",b); else sprintf(rbuf,"%dN/%d",a,b);

    start = get_time();
    if ( maxN > minN ) {
        if ( create || !hintfile) {
            if ( create && !hintfile ) { fprintf(stderr, "You must use the -h parameter to specify the hint-file to be created.\n"); return -1; }
            if ( !hintfile) fprintf(stderr, "No hint-file specified so t-choices will not be save (use -c -h hintfile.txt options to save them).\n");
            FILE *fp = hintfile ? fopen(hintfile, "w") : 0;
            if ( hintfile && !fp ) { fprintf(stderr, "Error creating hint-file %s\n", hintfile); return -1; }
            int64_t N = minN;
            while ( N <= maxN ) {
                t = tbound(N,a,b,fast,exhaustive,verbosity,verify);
                if ( b*t < a*N ) break;
                if ( verbosity >= 0 ) fprintf (stderr,"t(%ld) >= %ld (t-%s >= %ld) (%.3fs)\n", N, t, rbuf, t-cdiv(a*N,b), get_time()-start);
                if ( fp ) fprintf(fp,"%ld:%ld\n",N,t);
                N = b*t/a+1;
            }
            if ( N > maxN ) {
                fprintf (stderr,"Verified the %s Erdős-Guy-Selfridge conjecture for all N in [%ld,%ld] (%.3fs)\n",rbuf,minN,maxN,get_time()-start);
            }
            else if ( N == minN ) fprintf (stderr,"Unable to verify the %s Erdős-Guy-Selfridge conjecture for N=%ld in (%.3fs)\n",rbuf,minN,get_time()-start);
            else fprintf (stderr,"Only able to verify the %s Erdős-Guy-Selfridge conjecture for N in [%ld,%ld] (%.3fs)\n",rbuf,minN,N-1,get_time()-start);
        } else {
            FILE *fp = fopen(hintfile,"r"); if (!fp) { fprintf(stderr, "Error opening hint-file %s\n", hintfile); return -1; }
            char buf[256];
            int64_t minV = 0, maxV = 0;
            while ( fgets(buf,sizeof(buf),fp) ) {
                int64_t N = atol(buf);
                char *s = strchr(buf,':'); if (!s) { fprintf(stderr, "Error parsing line %s\n", buf); return -1; }
                t = atol(s+1);
                if ( b*t < a*N ) { fprintf(stderr, "Invalid N:t in hint file: %d*%ld < %d*%ld\n", b, t, a, N); return -1; }
                double timer = get_time();
                if ( tfac(N,t,fast,0,verbosity,verify,0) < N ) { fprintf (stderr, "Failed to verify t(%ld) >= %ld !\n", N,t); return -1; }
                if (!minV) {
                    if ( N > minN ) { fprintf(stderr, "Hint file starting N=%ld above range minimum %ld\n", N, minN); return -1; }
                    minV = N;
                    maxV = b*t/a;
                } else {
                    if ( N > maxV+1 ) { fprintf(stderr, "Hint file starting N=%ld leaves a gap!\n", N); return -1; }
                    if ( b*t <= a*maxV ) { fprintf (stderr, "Hint at N=%ld did not extend verified range!\n", N); return -1; }
                    maxV = b*t/a;
                }
                if ( verbosity >= 0 ) printf ("t(%ld) >= %ld (%.3fs)\n", N, t, get_time()-timer);
                if ( maxV >= maxN ) break;
            }
            if ( minV > minN || maxV < maxN ) { fprintf (stderr, "Hint file only allowed verification [%ld,%ld]\n", minV,maxV); return -1; }
            fprintf (stderr,"Verified the %s Erdős-Guy-Selfridg conjecture for N in [%ld,%ld] (%.3fs)\n",rbuf,minN,maxN,get_time()-start);
        }
    } else {
        int64_t N = minN;
        if ( t && exhaustive ) { t=0; fprintf(stderr,"Ignoring specified value of t and searching for optimal value\n"); }
        if ( !t ) {
            t = tbound(N,a,b,fast,exhaustive,verbosity,verify);
            if ( t ) printf("t(%ld) >= %ld (%s %s) with t-%s = %ld (%.3fs)\n", N, t, exhaustive ? "exhaustive" : "heuristic", fast ? "fast" : "greedy",rbuf,t-cdiv(a*N,b),get_time()-start);
            else fprintf(stderr,"failed to prove t(%ld) >= %ld (%.3fs)\n", N, cdiv(N,3), get_time()-start);
        } else {
            int64_t cnt = tfac(N,t,fast,0,verbosity,verify,dumpfile);
            if ( cnt >= N ) printf("t(%ld) >= %ld with %ld extra factors (%.3fs)\n", N, t, cnt - N, get_time()-start);
            else fprintf (stderr,"failed to prove t(%ld) >= %ld with %ld missing factors (%.3fs)\n", N, t, N - cnt, get_time()-start);
        }
    }
    unsetup();
    return 0;
}
