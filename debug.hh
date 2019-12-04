#pragma once

#ifndef DEBUG 
//#define DEBUG 0 // set debug mode
#endif

// debug functions
#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)

#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)

#endif
