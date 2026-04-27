This directory contains a fast, memory efficient implementation of a variant of the greedy algorithm that yields slightly worse lower bounds on $$t(N)$$ in general, but is able to verify the Erdős-Guy-Selfridge conjecture for $$N \in [10^6,10^{14}]$$ in less than ten minutes using just a few GB of memory (the space complexity of our implementation is $$O(N^{\frac{5}{8}+\epsilon})$$).

An implementation of the standard greedy algorithm is included in the same program, but its range is restricted to $$N\le 10^9$$ to keep the memory requirement low.

This implementation depends on Kim Walisch's [primesieve](https://github.com/kimwalisch/primesieve) and [primecount](https://github.com/kimwalisch/primecount) libraries for enumerating and counting primes, which must be installed before you compile the code in this directory.

These two libraries are used to efficiently handle factors divisible by primes in the range $$[\sqrt{t},N]$$. In this range the greedy algorithm can always use the optimal cofactor $$m = \lceil\frac{t}{p}\rceil$$ to construct $$n = v_p(N!)$$ factors $$mp \ge t$$ that are as small as possible (given that they are divisible by $p$).  To improve efficiency, the algorithm subdivides the interval $$[\sqrt{t},N]$$ into regions within which the values of $m$ and $n$ are constant and simply counts the primes in each region (for small regions it uses `primesieve` to enumerate primes, for large regions $$[a,b]$$ it uses `primecount` to compute $$\pi(b)-\pi(a-1))$$.  See the paper for more details.

To verify the Erdos-Guy-Selfridge conjecture for $$N$$ in $$[67425,10^{11}]$$ (which is sufficient for the portion of our proof of Theorem 1(iii) that depends on the greedy algorithm), after installing [primesieve](https://github.com/kimwalisch/primesieve) (version 12.6 or later, note that many libprimesieve-dev apt packages are older than this) and [primecount](https://github.com/kimwalisch/primecount) (version 7.14 or later), and compiling `egs.c` using `build.sh`, you can type
```
./egs -h hint_67425_1e6_exhaustive_greedy.txt 67425-1e6
./egs -f -h hint_1e6_1e11_heuristic_fast.txt 1e6-1e11
```
This should take 10-20 seconds (depending on your CPU). Note that you need to include the `-f` option to use the fast variant of the greedy algorithm in the second line, as the standard greedy algorithm will be used by default, and our C implementation of the standard greedy algorithm only supports $$N\le 10^9$$.

To verify the range $$[10^{11},10^{14}]$$ use
```
./egs -f -h hint_1e11_1e14_heuristic_fast.txt 1e11-1e14
```
which should take 5-10 minutes (again, depending on your CPU).

The *hint-files* contain lists of pairs $$(N,t)$$ such that the algorithm is known to produce at least $$N$$ factors of size at least $$t \ge N/3$$ on these inputs, such that the range of $$N$$ is completely covered by the intervals $$[N,3t]$$. Here we are exploiting the fact that $$M\ge N$$ implies $$t(M) \ge t(N)$$, so $$t(N)\ge t$$ implies $$t(M)\ge t \ge M/3$$ for $$N\le M\le 3t$$.

To recreate the hint-file `hint_1e6_1e11_heuristic_fast.txt` (or create a new one for a different range) you can type
```
./egs -f -c -h hint_1e6_1e11_heuristic_fast.txt 1e6-1e11
```
which should take a few minutes.

You can compute lower bounds on $$t(N)$$ for particular $$N$$ by typing `egs N` or `egs -f N` (this will use a heuristic search for the best bound that can be easily found)  or use `egs -e N` or `egs -f -e N` to do an exhaustive search for the best bound that the algorithm can prove for the given $$N$$ (this computation will take longer and use all available cores).

```
$ ./egs 1e6
t(1000000) >= 340790 (heuristic greedy) with (t-ceil(N/3)) = 7456 (0.012s)
$ ./egs -f 1e6
t(1000000) >= 337642 (heuristic fast) with (t-ceil(N/3)) = 4308 (0.003s)
$ ./egs -e 1e6
t(1000000) >= 342303 (exhaustive greedy) with (t-ceil(N/3)) = 8969 (0.272s)
$ ./egs -e -f 1e6
t(1000000) >= 338782 (exhaustive fast) with (t-ceil(N/3)) = 5448 (0.041s)

```

The following files contain bounds on $$t(N)$$ produced as above:
- `tbounds_heuristic_greedy_1e4_1e9.txt` contains bounds on $$t(c\cdot 10^n)$$ in the range $$10^4\le 10^9$$ obtained via heuristic search with the greedy algorithm (`egs N`)
- `tbounds_heuristic_fast_1e4_1e14.txt` contains bounds on $$t(c\cdot 10^n)$$ in the range $$10^4\le 10^{14}$$ obtained via heuristic search method with the fast variant of the greedy algorithm (`egs -f N`)
- `tbounds_exhaustive_greedy_1e4_1e9.txt` contains bounds on $$t(c\cdot 10^n)$$ in the range $$10^4\le 10^9$$ obtained via exhaustive search with the greedy algorithm (`egs -e N`)
- `tbounds_exhaustive_fast_1e4_1e9.txt` contains bounds on $$t(c\cdot 10^n)$$ in the range $$10^4\le 10^9$$ obtained via exhaustive search method with the fast variant of the greedy algorithm (`egs -e -f N`)

To see all the command line options that are supported, run `egs` without specifying any arguments:
```
$ ./egs
Usage: egs [-v level] [-h filename] [-d filename] [-r] [-c] [-e] [-f] N-range [t or t/N ratio]
       -v level      integer verbosity level -1 to 4 (optional, default is 0)
       -h filename   hint-file with records N:t (required if range of N is specified)
       -d filename   output-file to dump factorization to (one factor per line, only valid if t is specified)
       -r            verify factorization (set automatically if dump is specified)
       -c            create hint-file rather than reading it (must be specified in combination with -h)
       -e            use the best t for which the algorithm can prove t(N) >= t (optional)
       -f            use fast version of greedy algorithm
       -m            exponent for primecount/primesieve cutoff, must lie in [1/6,1/3]
       N-range       integer N or range of integers minN-maxN (required, scientific notation supported)
       t             integer t to use for single N (optional, a good t will be determined if unspecified)
       t/N ratio     a/b with integers a,b>0 specifying t = ceil(aN/b) (optional, set to 1/3 if unspecified)
```
