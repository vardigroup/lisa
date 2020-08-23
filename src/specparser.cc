/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "specparser.hh"

spec_ptr
parse_spec(const char* file_name)
{
    spec_ptr p = new spec;

    json j;
    ifstream ispecfile(file_name);
    ispecfile >> j;
    cout << j << endl;

    p->descr = j[DESCR];

    vector<string> assumps = j[ASSUMPTIONS];
    cout << "Assumptions: " << endl;
    for(string assump: assumps)
    {
        auto pf1 = spot::parse_infix_psl(assump.c_str());
        formula conjunct = pf1.f;
        p->assumptions.push_back(conjunct);
        cout << pf1.f << endl;
    }
    cout << "Guarantees: " << endl;
    vector<string> guarantees = j[GUARANTEES];
    for(string guarantee: guarantees)
    {
        auto pf1 = spot::parse_infix_psl(guarantee.c_str());
        formula conjunct = pf1.f;
        p->guarantees.push_back(conjunct);
        cout << pf1.f << endl;
    }

    for(string input_ap: j[INPUTS])
    {
        p->input_aps.push_back(input_ap);
        cout << "in: " << input_ap << endl;
    }

    for(string output_ap: j[OUTPUTS])
    {
        p->output_aps.push_back(output_ap);
        cout << "out: " << output_ap << endl;
    }

    for(string ubobs_ap: j[UNOBSERS])
    {
        p->unobservable_aps.push_back(ubobs_ap);
        cout << "unobservable: " << ubobs_ap << endl;
    }

    if(j[TYPE] == DOMINANT)
    {
        p->start_type = strategy_type::DOMINANT;
    }else
    {
        p->start_type = strategy_type::WINNING;
        cout << "Winning" << endl;
    }

    if(j[SEMANTICS] == MEALY)
    {
        p->start_semantics = strategy_semantics::MEALY;
    }else
    {
        p->start_semantics = strategy_semantics::MOORE;
        cout << "Moore" << endl;
    }
    
    return p;
}


int
main(int argc, char** args)
{
    char* name = args[1];
    spec_ptr sp = parse_spec(name);

    delete sp;
    return 0;
}