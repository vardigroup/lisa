Overview
=======

Lisa is an (a) LTLf to DFA conversion tool, and (b) an LTLf synthesis tool. 

It is publicly available under the license GNU GPL v3.


Requirements
-----------------------------------

Lisa requires a C++14-compliant compiler.  G++ 5.x or later should work.

Third-party dependencies
-----------------------------------

* [Spot model checker version>=2.9](https://spot.lrde.epita.fr/)

* [MONA](https://github.com/liyong31/MONA)

* [CUDD library](https://github.com/KavrakiLab/cudd.git)

Lisa relies on Spot and MONA to construct a DFA from a small LTLf formula.
When constructing a DFA from an LTLf formula with MONA, Lisa translates an LTLf formula to a formula in first order logic, which is then fed into MONA.

Complilation steps
=======

In the following we assume that we will compile Lisa on a Ubuntu system.

1. Install Spot

    Lisa needs Spot to convert an LTLf to a DFA and to perform intersection of DFAs with explicit state representation.

    * Get the latest version of Spot from https://spot.lrde.epita.fr/install.html.

    * Uncompress Spot and follow install intructions in README to install Spot. Note that the compilation of Spot may take a while.
    
            ./configure && make && sudo make install

    * Type ltl2tgba -f "F a" in command line, you expect to see some output starting with "HOA: v1".

2. Install CUDD

    Lisa employs CUDD for symbolic DFA minimization.

    * Uncompress cudd.zip or get CUDD from https://github.com/KavrakiLab/cudd.git.

    * Install CUDD:

            ./configure --enable-silent-rules --enable-obj --enable-dddmp --prefix=[install location]

            sudo make install

        If you encounter an error about aclocal, this might be fixed as follows.
        * Not having automake:
        
            sudo apt-get install automake
        * Needing to reconfigure, do this before configuring:

            autoreconf -i
        * Using a version of aclocal other than 1.14:

            modify the version 1.14 in configure file accordingly.

3. Install MONA

    Lisa needs MONA to convert a formula in first order logic to a DFA.

    * Go to mona-1.4-17 directory and follow the install instructions in INSTALL.
    
            ./configure && make && sudo make install-strip

    **NOTE** 
    In BDD/bdd.h, the original table size is defined as #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000 (=2^24), which is too small for large DFA generation.
    We modify it to #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000000 (=2^36), so to allow MONA have larger table size during DFA construction.
    Note that MONA has explicit state representation but encodes the labels on transition symbolically.
    For more details on the representation of DFA in MONA, we refer to https://www.brics.dk/mona/mona14.pdf.
    
6. Compile Lisa

    * Compile Lisa with Make:
    
            make T1

        or compile Lisa in command line:

            g++ lisa.cc dfwavar.cc dfwa.cc spotutil.cc ltlf2fol.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc -o lisa -lspot -lbddx -lcudd  -O3

Input format
=======

Lisa accepts LTLf formulas given as a .ltlf file written in SPOT format. 

For synthesis, it also requires a .part file. The .part file indicates the input and output propostitions for the synthesis task. 

Example .ltltf file

```
((COUNTER0 <-> INITCOUNTER0)) && (G (CARRY0 <-> INC)) && (G ((X COUNTER0 -> !(COUNTER0 <-> CARRY0)) && (X !COUNTER0 -> (COUNTER0 <-> CARRY0)))) && ((G ((!INC -> X INC)) -> F (!COUNTER0)))
```

Example .part file

```
.inputs INITCOUNTER0 INC
.outputs COUNTER0 CARRY0
```

Command line usage
=======

If you type ./lisa -h in command line, you should see the following command line usage:

```
Usage: lisa [OPTION...] [FILENAME[/COL]...]
Read a formula file and output the number of states of the constructed DFA

 Input options:
 -h                    show this help page
 -exp                  use only explicit method (default false)
 -min                  minimize the last symbolic DFA (default false)
 -syn                  synthesize after DFA construction (default false)
 -bdd                  use buddy for DFA minimization
 -cdd                  use cudd for DFA minimization
 -nap  <int>           number of atomic propositions for calling mona (default 7)
 -npr  <int>           number of products for calling minimization (default 6)
 -nia  <int>           number of states of individual DFA for calling symbolic approach (default 800)
 -npa  <int>           number of states of product DFA for calling symbolic approach (default 2500)
 -lst  <int>           number of last automata for calling symbolic approach (default -1)
 -out                  print out the wining strategy if realizable
 -part <file>          the file specifying the input and output propositions
 -ltlf <file>          the file specifying the input LTLf formula
 -env                  environment plays first
```

For LTLf to DFA construction
==
1. To use the default setting to construct a DFA from an LTLf formula, type
    
        ./lisa -ltlf ./examples/ltlf3377.ltlf
    
    You are expected to see the output ending with "Number of states (or nodes) is: 78626".
    In the DFA generation, the symbolic representation of DFA has been triggered.
    As we usually do not count the number of states in a symbolic DFA, we only output the number of nodes in the BDD representation of the transition relation of the output DFA.
    
    Note that the two parameters t1 and t2 mentioned in the submission correspond to respectively the options -nia <int> and -npa <int>.
    Recall that the switch from explicit-state form to symbolic-state form is triggered if either the smallest minimal DFA has more than t1 states, or if the product of the number of states in the two smallest minimal DFAs is more than t2.
    For example, the following command
    
        ./lisa -ltlf ./examples/ltlf3377.ltlf -nia 0 -npa 0
    
    corresponds to pure compositional symbolic DFA generation. You are expected to see the output ending with "Number of states (or nodes) is: 40745".

<!--  
2. In order to know how many states in the minimal DFA of the last symbolic DFA, type

        ./lisa -ltlf ./examples/ltlf3377.ltlf -lst
    
    to minimize the output symbolic DFA. The minimization of the DFA may take a while to terminate.
    You are expected to see the output ending with "Number of states in the minimal DFA is: 3376".
    That is, the number of states in the minimal DFA is 3376.
-->

2. In order to use only explicit state representation in DFA generation, type
    
        ./lisa -ltlf ./examples/ltlf3377.ltlf -exp
    
    The DFA generation with explicit states always terminates and returns a minimal DFA corresponding to the input formula.
    You are expected to see the output ending with "Number of states (or nodes) is: 3377".
    That is, the number of states in the output explicit DFA is 3377.
    Note that the DFA with explicit state representation has one more sink state than the DFA with symbolic representation.
    This is because that we use Spot data structure to store DFAs in explicit-state form but Spot only supports automata accepting infinite words.
    Therefore, we need to use an extra sink state for the Spot automata as an indicator that it is the end of a finite trace for DFAs.

For LTLf synthesis
==

1. To perform the synthesis step after DFA generation, type
    
        ./lisa -ltlf ./examples/ltlf3377.ltlf -part ./examples/ltlf3377.part -syn
    
    where the input and output variables in the LTLf formula are specified in the file ltlf3377.part.
    You are expected to see that Lisa outputs "UNREALIZABLE", as this specification is unrealizable.

2. To let the environment move first in the DFA game for LTLf synthesis, type

        ./lisa -ltlf ./examples/ltlf3377.ltlf -part ./examples/ltlf3377.part -syn -env

Syntax
==

The Linear Temporal Logic over finite traces (LTLf) has the same syntax as LTL.
Given a set P of propositions, the syntax of LTLf formulas supported by Spot is as follows:
```
φ ::= 1 | 0 | p | !φ | φ1 && φ2 | φ1 || φ2 | φ1 -> φ2 
      | φ1 <-> φ2 | X φ | X[!] φ | F φ | G φ | φ1 U φ2 | φ1 R φ2 | φ1 W φ2

```
where p ∈ P. Here 1 and 0 represent *true* and *false* respectively.
X (weak Next), X[!] (strong Next), F (Finally), G (Globally), U (Until), R (Release) and W (weak Until) are temporal operators.
We have that X[!] φ ≡ ! (X !φ), F φ = !(G !φ), φ1 U φ2 ≡ !(!φ1 R !φ2) and φ1 W φ2 ≡ G φ1 || (φ1 U φ2).
As usuall, we also have that F φ ≡ 1 U φ and G φ ≡ 0 R φ.

For the semantics of LTLf formula, we refer to [IJCAI13 paper](https://www.cs.rice.edu/~vardi/papers/ijcai13.pdf).
Specially, Spot supports a weak next and a strong next.
    
Weak next: *X a* is true if *a* holds at next step or if there is no next step.
In particular, *X(0)* is true iff there is no successor.
    
Strong next: *X[!] a* is true if *a* holds at next step and there must be a next step.
In particular *X[!]1* is true iff there is a successor.

## Acknowledgment
- Alexandre Duret-Lutz : [Spot](https://spot.lrde.epita.fr/)
- Jørn Lind-Nielsen: [BuDDy](http://vlsicad.eecs.umich.edu/BK/Slots/cache/www.itu.dk/research/buddy/)
- Fabio Somenzi: [CUDD](https://github.com/ivmai/cudd)
- Tom van Dijk: [Sylvan](https://github.com/trolando/sylvan)
- Shufang Zhu: [Syft](https://github.com/saffiepig/Syft)
    

