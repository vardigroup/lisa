Overview
=======

Lisa is an LTLf synthesis tool publicly available under the license GNU GPL v3.

Compilation instructions
=======

Requirements
------------

Lisa requires a C++14-compliant compiler.  G++ 5.x or later should work.

Third-party dependencies
----------------------------------

* Spot (Available under the license GNU GPL v3 at https://spot.lrde.epita.fr/)

* Syft (Available under the MIT license at https://github.com/liyong31/Syft)

* MONA (Available under the license GNU GPL v2 at https://www.brics.dk/mona/)

* CUDD (Available from https://github.com/ivmai/cudd)

Lisa relies on Spot and MONA to construct a DFA from a small LTLf formula.
When constructing a DFA from an LTLf formula with MONA, Lisa requires Syft to
translate an LTLf formula to a formula in first order logic, which is then fed into MONA.

Complilation steps
----------------------------------

In the following we assume that we will compile Lisa on a Ubuntu system.

==== 1) Install Spot ====

Lisa needs Spot to convert an LTLf to a DFA and to perform intersection of DFAs with explicit state representation.

1. Get the latest version of Spot from https://spot.lrde.epita.fr/install.html.

2. Uncompress Spot and follow install intructions in README to install Spot. Note that the compilation of Spot may take a while.
   For instance ./configure && make && sudo make install

3. Type ltl2tgba -f "F a" in command line, you expect to see some output starting with "HOA: v1".

==== 2) Install CUDD ====

Syft needs CUDD to perform BDD operations and Lisa employs CUDD for symbolic DFA minimization.

1. Uncompress cudd.zip or get CUDD from https://github.com/ivmai/cudd.

2. Install CUDD:
    ./configure --enable-silent-rules --enable-obj --enable-dddmp --prefix=[install location]
    sudo make install

    If you encounter an error about aclocal, this might can be fixed as follows.
    a. Not having automake:
       sudo apt-get install automake
    b. Needing to reconfigure, do this before configuring:
       autoreconf -i
    c. Using a version of aclocal other than 1.14:
       modify the version 1.14 in configure file accordingly.

==== 3) Install MONA ====

Lisa needs MONA to convert a formula in first order logic to a DFA.

1. Go to mona-1.4-17 directory and follow the install instructions in INSTALL.
    For instance ./configure && make && sudo make install-strip

**NOTE**
In BDD/bdd.h, the original table size is defined as #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000 (=2^24), which is too small for the DFA generation comparison.
For a fair comparison with MONA, we modify it to #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000000 (=2^36), so to allow MONA have larger table size during DFA construction.
Note that MONA has explicit state representation but encodes the labels on transition symbolically.
For more details on the representation of DFA in MONA, we refer to https://www.brics.dk/mona/mona14.pdf.
    
==== 4) Compile Syft ====

Lisa needs the ltlf2fol tool in Syft to rewrite an LTLf formula as a formula in first order logic.
Please make sure that you have network connection during the installation of flex, bison and boost.

1. Install flex and bison:
    sudo apt-get install flex bison

2. Install BOOST:
    sudo apt-get install libboost-dev

3. Go to Syft directory and make build folder under Syft folder:

   mkdir build && cd build

2. Run CMake to generate the makefile:

   cmake ..

3. Compile using the generated makefile:

   make
   
Then you should be able to find two executables ltlf2fol and Syft inside the folder Syft/build/bin.

**NOTE**
For a fair comparison with Syft, we have enabled dynamic variable ordering in Syft, as Lisa and Part tool also use this dynamic option.
As mentioned in the submission, the improved version of Syft can exhibit as much as 75% reduction in runtime compared to Syft.

==== 5) Compile Lisa ====

1. Copy the executable file ltlf2fol to lisa folder.

2. Compile Lisa with Make:
    
    make

    or compile Lisa in command line:

    g++ lisa.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc  -o lisa -lspot -lbddx -lcudd -O3

*** Command line usage ***

If you type ./lisa -h in command line, you should see the following command line usage:

Usage: lisa [OPTION...] [FILENAME[/COL]...]
Read a formula file and output the number of states (nodes) of the constructed DFA

 Input options:
 -h                    show this help page
 -exp                  use only explicit method (default false)
 -syn                  synthesize after DFA construction (default false)
 -lst                  minimize the output DFA (default false)
 -nia  <int>           number of states of individual DFA for calling symbolic approach (default 800)
 -npa  <int>           number of states of product DFA for calling symbolic approach (default 2500)
 -part <file>          the file specifying the input and output propositions
 -ltlf <file>          the file specifying the input LTLf formula
 -env                  environment plays first in synthesis

1. To use the default setting to construct a DFA from an LTLf formula, type
    
    ./lisa -ltlf ../examples/ltlf3377.ltlf
    
    You are expected to see the output ending with "Number of states (or nodes) is: 78626".
    In the DFA generation, the symbolic representation of DFA has been triggered.
    As we usually do not count the number of states in a symbolic DFA, we only output the number of nodes in the BDD representation of the transition relation of the output DFA.
    
    Note that the two parameters t1 and t2 mentioned in the submission correspond to respectively the options -nia <int> and -npa <int>.
    Recall that the switch from explicit-state form to symbolic-state form is triggered if either the smallest minimal DFA has more than t1 states, or
    if the product of the number of states in the two smallest minimal DFAs is more than t2.
    For example, the following command
    
    ./lisa -ltlf ../examples/ltlf3377.ltlf -nia 0 -npa 0
    
    corresponds to pure compositional symbolic DFA generation. You are expected to see the output ending with "Number of states (or nodes) is: 40745".
    
2. In order to know how many states in the minimal DFA of the last symbolic DFA, type

    ./lisa -ltlf ../examples/ltlf3377.ltlf -lst
    
    to minimize the output symbolic DFA. The minimization of the DFA may take a while to terminate.
    You are expected to see the output ending with "Number of states in the minimal DFA is: 3376".
    That is, the number of states in the minimal DFA is 3376.
    
3. In order to use only explicit state representation in DFA generation, type
    
    ./lisa -ltlf ../examples/ltlf3377.ltlf -exp
    
    The DFA generation with explicit states always terminates and returns a minimal DFA corresponding to the input formula.
    You are expected to see the output ending with "Number of states (or nodes) is: 3377".
    That is, the number of states in the output explicit DFA is 3377.
    Note that the DFA with explicit state representation has one more sink state than the DFA with symbolic representation.
    This is because that we use Spot data structure to store DFAs in explicit-state form but Spot only supports automata accepting infinite words.
    Therefore, we need to use an extra sink state for the Spot automata as an indicator that it is the end of a finite trace for DFAs.

4. To perform the synthesis step after DFA generation, type
    
    ./lisa -ltlf ../examples/ltlf3377.ltlf -part ../examples/ltlf3377.part -syn
    
    where the input and output variables in the LTLf formula are specified in the file ltlf3377.part.
    You are expected to see that Lisa outputs "UNREALIZABLE", as this specification is unrealizable.

5. To let the environment move first in the DFA game for LTLf synthesis, type

    ./lisa -ltlf ../examples/ltlf3377.ltlf -part ../examples/ltlf3377.part -syn -env
    
    Note that in the experiments conducted for the submission, we always let environment move first in the DFA game.
    Still, Lisa outputs "UNREALIZABLE", as this specification is unrealizable.

*** Execute the experiments ***

Due to the space limit, we are unable to provide the raw data of the experiments as they amount to 1.7GB.
We can provide them upon request via an anonymous link durin gauthor response period.

For reviewers who wish to run the experiments on their own machine, we provide two shell scripts in benchmarks folder for running the experiments.
We assume that all executables are available in their respective folders as follows:

lisa      --- lisa/
ltlffol   --- lisa/ and Syft/build/bin
Syft      --- Syft/build/bin

==== DFA generation ====

In this experiment, we compare Lisa to Lisa-Explicit, Lisa-Symbolic and MONA-based method.

1. Go to benchmarks directory.

    type chmod +x ./run_dfa_conversion.sh && ./run_dfa_conversion.sh
    
    Currently the time limit is set to 3600 through a line of code time=3600 (1 hour).

==== LTLf synthesis ====

In this experiment, we compare LisaSynt to Syft+.

1. Go to benchmarks directory.

    type chmod +x ./run_synthesis.sh && ./run_synthesis.sh
    
    Currently the time limit is set to 28800 (8 hour).

The log files for each benchmark set can be found in their respective folders in the log folder inside the benchmarks folder.

