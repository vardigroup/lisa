/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <nlohmann/json.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>      // std::setw

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>

// for convenience
using json = nlohmann::json;
using namespace std;
using namespace spot;

// description for the specification
const char* DESCR = "description";
// assumption for the specification
const char* ASSUMPTIONS = "assumptions";
// guarantees for the specification
const char* GUARANTEES = "guarantees";
// inputs for the specification
const char* INPUTS = "inputs";
// outputs for the specification
const char* OUTPUTS = "outputs";
// types for the synthesized strategies
const char* TYPE = "type";
// semantics for the synthesized strategies
const char* SEMANTICS = "semantics";
// unobservable inputs 
const char* UNOBSERS = "unobservable";
// 
//const char* WINNING = "winning";
//
const char* DOMINANT = "dominant";

const char* MEALY = "mealy";


enum class strategy_type { WINNING, DOMINANT };
enum class strategy_semantics { MOORE, MEALY };

struct spec
{
    string descr;

    strategy_type start_type;
    strategy_semantics start_semantics;

    vector<string> input_aps;
    vector<string> output_aps;
    vector<string> unobservable_aps;

    vector<formula> assumptions;
    vector<formula> guarantees;


    friend std::ostream& operator<<(std::ostream& os, const spec& obj);

};

typedef spec* spec_ptr;

spec_ptr
parse_spec(const char* file_name);





