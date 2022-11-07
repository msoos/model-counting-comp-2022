# GPMC

GPMC is an exact model counter for CNF formulas. The current version of GPMC supports model counting, weighted model counting, projected model counting, and weighted projected model counting.

The source codes of this software are based on those of MiniSat-based Sat solver [Glucose 3.0](https://www.labri.fr/perso/lsimon/glucose/) and [SharpSAT](https://github.com/marcthurley/sharpSAT) 12.08.1.
And also, inspired by [sharpSAT-TD](https://github.com/Laakeri/sharpsat-td), the current version incorporates 
a preprocessor for simplifying input formulas, and uses tree decomposition in decision heuristics. 
In the implementation, we use a part of source files of the preprocessor in sharpSAT-TD, and use [FlowCutter](https://github.com/kit-algo/flow-cutter-pace17) for tree decomposition.

## Installation
See `INSTALL.md`.

## Usage
```
$ ./gpmc [options] [file]
```
Options:

-  -mode=<0..3>
    - 0 : Model Counting
    - 1 : Weighted Model Counting
    - 2 : Projected Model Counting
    - 3 : Weighted Projected Model Counting 

-  -cs=\<1..intmax>  
    set maximum component cache size (MB) (default: 4000)  
    NOTE: This bound is only for component cache. GPMC may use much more memory.

-  -bj / -no-bj  
    turn on startup limited backjumping (default: on)

## Input file format
Input boolean formulas should be in DIMACS CNF format, together with additional information, weights and projection variables.
GPMC supports the [MCC2021 format](https://mccompetition.org/assets/files/2021/competition2021.pdf). 

Example:

```
p cnf 3 4
c p show 2 3 0
c p weight 2 0.6 0
c p weight -2 0.4 0
-1 2 0
-2 3 0
2 -3 0
1 -2 3 0
```
The set of projection variables is { 2, 3 }.
The weight of the positive literal of var 2 is 0.6, and that of the negative one is 0.4.

## Author
[Kenji Hashimoto](https://www.trs.cm.is.nagoya-u.ac.jp/~k-hasimt/index-e.html),
Nagoya University, Japan.
