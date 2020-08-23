/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <stdlib.h>
#include <error.h>
#include <argp.h>

// This file is adapted from Argp example #4 at https://www.gnu.org/software/libc/manual/html_node/Argp-Example-4.html

const char *argp_program_version =
  "lisa v0.1";

const char *argp_program_bug_address =
  "<liyong@ios.ac.cn>";

/* Program documentation. */
static char doc[] =
  "Translate an LTLf formula to a DFA and perform synthesis if necessary";
/*
  \
options\
\vThis part of the documentation comes *after* the options;\
 note that the text is automatically filled, but it's possible\
 to force a line-break, e.g.\n<-- here.";
*/

/* A description of the arguments we accept. */
static char args_doc[] = "[SPECIFICATION...]";

/* Keys for options without short-options. */
#define OPT_ABORT  1            /* --abort */

/* The options we understand. */
static struct argp_option options[] = {
    
  // second argument is the key for parsing arguments
  /*-------------------------------------------*/
  { 0,0,0,0, "Input:", 1 },
  {"formula",     'f', "STRING" , 0, "Process the specification STRING"},
  // use json format
  {"file",     'F', "FILE" , 0, "Process the specification in FILE"},
  {"verbose",  'v', 0,       0, "Produce verbose output" },
  {"quiet",    'q', 0,       0, "Don't produce any output" },
  //{"silent",   's', 0,       OPTION_ALIAS },
  {"output",   'o', "FILE",  OPTION_ARG_OPTIONAL, "Output to FILE instead of standard output" },

  /*-------------------------------------------*/
  {0,0,0,0, "DFA construction:", 2},
  {"explicit",   'e', 0, OPTION_ARG_OPTIONAL, "Only explicit approach for DFA construction\n(Default: false)"},
  {"individual",  'i', "INT", 0, "Switch to symbolic approach when the number of states of an individual DFA exceeds INT\n(Default: 800)"},
  {"product",    'p', "INT", 0, "Switch to symbolic approach when the number of states of a product DFA exceeds INT\n(Default: 2500)"},
  
  /*-------------------------------------------*/
  {0,0,0,0, "Synthesis:", 3 },
  {"synthesize",   's', 0, OPTION_ARG_OPTIONAL,
   "Synthesize a strategy from the specification"},
   
  /*-------------------------------------------*/
  {0,0,0,0, "DFA minimization:", 4 },
  {"minimize",   'm', 0, 0,
   "Minimize the DFA for the specification"},
   
  /*-------------------------------------------*/
  {0,0,0,0, "BDD choice:", 5},
   {"cudd",   'c', 0, 0,
   "Apply CUDD for DFA minimization"},
   {"buddy",   'b', 0, 0,
   "Apply BuDDy for DFA minimization"},
   {"sylvan",   'y', 0, 0,
   "Apply Sylvan for DFA minimization"},
  
  /*-------------------------------------------*/
  {0,0,0,0, "Miscellaneous options:", -1},
  { nullptr, 0, nullptr, 0, nullptr, 0 }
  //{"help",  'h', 0,       0, "print this help page" },
  //{0}
};

// command line parser
static struct opt_t
{
	const char* _ltlfile_name = nullptr;
	const char* _parfile_name = nullptr;

	bool _symbolic = true;
	bool _minimization = false;

	unsigned _num_ap_for_mona = 7;
	unsigned _num_product = 6;
	unsigned _num_st_for_single = 800;
	unsigned _num_st_for_product = 2500;
	int _num_last_automata = -1;

	bool _synthesis = false;
	bool _out_start = false;
	bool _env_first = false;

	uint8_t _bdd = 0;

}* opt;

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *arg1;                   /* arg1 */
  char **strings;               /* [string...] */
  int silent, verbose, abort;   /* `-s', `-v', `--abort' */
  char *output_file;            /* file arg to `--output' */
  int repeat_count;             /* count arg to `--repeat' */
};


