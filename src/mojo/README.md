# smoothfac

This directory contains an implementation of the `smoothfac` algorithm for
calculating factorizations of $`N!`$ with large factors.


## Getting Started

The program is written in [Mojo](https://docs.modular.com/mojo/manual/). To
install Mojo, follow the install instructions in the [getting started
guide](https://docs.modular.com/mojo/manual/get-started).

To run the program, change into this directory and start a Magic shell
```sh
cd src/mojo/
magic shell
```

The default mode of the program is to optimize lower and upper bounds for
$`t(N)`$. To find bounds for $`N=300000`$, run the following command
```sh
mojo run smoothfac.mojo 300000
```
You should see the output  
`OEIS result: N = 300000  t(N) = 101903  T_ub = 101906`  
as the final output line.

When the program is invoked with two arguments, the second argument is a filename
and the program writes the lower bound factorization to disk.

To run the program for a single fixed threshold $`t`$, invoke it with the
environment variable `SMTHFC_FIX_THRES` set to the desired value $`t`$:
```sh
SMTHFC_FIX_THRES=100000 mojo run smoothfac.mojo 300000
```
You should see the output  
`Total factors: 300445`  
as the final output line.


## Algorithm

There is a [PDF document](./smoothfac.pdf) containing some information on the
algorithm and implementation details.
