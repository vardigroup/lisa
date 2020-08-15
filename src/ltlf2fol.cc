/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#include "ltlf2fol.hh"

/*-------------------------------------------------------------------*/
// The code for this translation is adapted from Syft by Shufang Zhu
// Translate an LTLf formula to a first order logic
/*-------------------------------------------------------------------*/
// The input formula should be a formula in normal form
void ltlf_to_fol(ostream &os, formula &f)
{
  int c = 1;
  set<formula> aps;
  get_formula_aps(f, aps);
  // output the atomic propositions
  if (!aps.empty())
  {
    os << "m2l-str;" << endl;
    os << "var2 ";
    int count = 0;
    for (formula ap : aps)
    {
      //cout << "ap: " << ap << endl;
      if (count == 0)
      {
        os << ap.ap_name();
      }
      else
      {
        os << ", " << ap.ap_name();
      }
      count++;
    }
    os << ";" << endl;
  }
  // translate ltlf formulas to FOL formulas
  string res = translate2fol(f, 0, c);
  os << res << ";" << endl;
}

/*-------------------------------------------------------------------*/
// translate ltlf formula to its equivalent first order logic
// make sure input is in backus/negation normal form
// @f the input formula
// @t current position t
// @c unknown parameter for the code from Syft; it seems to be unused
/*-------------------------------------------------------------------*/
inline string
get_step_str(int t)
{
  if (t == 0)
  {
    return to_string(t);
  }
  else
  {
    return "x" + to_string(t);
  }
}
string translate2fol(formula &f, int t, int &c)
{
  string curs, ts;
  string exs, alls;
  string res;
  int cur;
  int count;
  formula r, lft, rgt;
  switch (f.kind())
  {
  // Not operation
  // (¬fol(φ, x))
  case op::Not:
    res = "~(";
    r = f[0];
    res += translate2fol(r, t, c);
    res += ")";
    break;
  // strong Next
  case op::strong_X:
    // X[!](0) = 0
    // X[!] (1) means that there must be a successor
    // Note: with finite semantics X[!](1)≠1.
    // next step is y = x(t+1)
    exs = "x" + to_string(t + 1);
    // current step t
    ts = get_step_str(t);
    // ((∃y)((y = x + 1) ∧ fol(φ, y)))
    res = "(ex1 " + exs + ": (" + exs + "=" + ts + "+1 & (";
    r = f[0];
    res += translate2fol(r, t + 1, c);
    res += ")))";
    break;
  // weak Next
  case op::X:
    //// X(1) = 1
    // X (0) means that there is no successor
    // We do not have X(0)=0 because that
    // is not true with finite semantics.
    // next step y = x(t+1)
    exs = "x" + to_string(t + 1);
    // current step t
    ts = get_step_str(t);
    //((x = last) ∨ ((∃y)((y = x + 1) ∧ fol(φ, y))))
    res = "((ex1 " + exs + ": (" + exs + "=" + ts + "+1 & (";
    r = f[0];
    res += translate2fol(r, t + 1, c);
    res += "))) | (" + ts + " = max $))";
    break;
  case op::F:
    exs = "x" + to_string(t + 1);
    ts = get_step_str(t);
    // (∃y)((x ≤ y ≤ last) ∧ fol(φ, y)
    res = "(ex1 " + exs + ": (" + ts + " <= " + exs + " & (";
    r = f[0];
    res += translate2fol(r, t + 1, c);
    res += ")))";
    break;
  case op::G:
    alls = "x" + to_string(t + 1);
    ts = get_step_str(t);
    //(∀z)((x ≤ z ≤ last) → fol(φ, z)))
    res = "(all1 " + alls + ": ((" + ts + " <= " + alls + ") => (";
    r = f[0];
    res += translate2fol(r, t + 1, c);
    res += ")))";
    break;
  case op::U:
    exs = "x" + to_string(t + 1);
    alls = "x" + to_string(t + 2);
    ts = get_step_str(t);
    lft = f[0];
    rgt = f[1];
    // ((∃y)((x ≤ y ≤ last) ∧ fol(φ2 , y) ∧ (∀z)((x ≤ z < y) → fol(φ1 , z))))
    res = "(ex1 " + exs + ": (" + ts + " <= " + exs + " & (";
    res += translate2fol(rgt, t + 1, c);
    res += ") & (all1 " + alls + ": (" + ts + " <= " + alls + " & " + alls;
    res += " < " + exs + " => (";
    res += translate2fol(lft, t + 2, c);
    res += ")))))";
    break;
  case op::R: //New
    exs = "x" + to_string(t + 1);
    alls = "x" + to_string(t + 2);
    ts = get_step_str(t);
    lft = f[0];
    rgt = f[1];
    // (((∃y)((x ≤ y ≤ last) ∧ fol(φ1 , y) ∧ (∀z)((x ≤ z ≤ y) → fol(φ 2 , z)))) ∨ ((∀z)((x ≤ z ≤ last) → fol(φ 2 , z))))
    res = "((ex1 " + exs + ": (" + ts + " <= " + exs + " & (";
    res += translate2fol(lft, t + 1, c);
    res += ") & (all1 " + alls + ": (" + ts + " <= " + alls + " & " + alls;
    res += " <= " + exs + " => (";
    curs = translate2fol(rgt, t + 2, c);
    res += curs;
    res += ")))))";
    res += "| (all1 " + alls + ": ((" + ts + " <= " + alls + " & " + alls + " <= max $) => (";
    res += curs;
    res += "))))";
    break;
  case op::Or:
    // list of Or operands
    //(fol(φ 1 , x) ∨ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + translate2fol(child, t, c) + ")";
      }
      else
      {
        res += " | (" + translate2fol(child, t, c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::And:
    // list of And operands
    // (fol(φ 1 , x) ∧ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + translate2fol(child, t, c) + ")";
      }
      else
      {
        res += " & (" + translate2fol(child, t, c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::tt:
    res += "(true)";
    break;
  case op::ff:
    res += "(false)";
    break;
  // atomic propositions
  case op::ap:
    res = "(";
    ts = get_step_str(t);
    res += ts + " in ";
    //const std::string& str = f.ap_name();
    res += "" + f.ap_name();
    res += ")";
    break;
  default:
    cerr << "Formula: " << f << ". ";
    throw runtime_error("Error formula in translate2fol()");
    exit(-1);
  }
  // cout<<res<<endl;
  return res;
}

/*-------------------------------------------------------------------*/
// The code for this translation is adapted from Syft by Shufang Zhu
// Translate an LTLf formula to a past first order logic
/*-------------------------------------------------------------------*/
// The input formula should be a formula in normal form
void ltlf_to_pfol(ostream &os, formula &f)
{
  int c = 1;
  set<formula> aps;
  get_formula_aps(f, aps);
  // output the atomic propositions
  if (!aps.empty())
  {
    os << "m2l-str;" << endl;
    os << "var2 ";
    int count = 0;
    for (formula ap : aps)
    {
      //cout << "ap: " << ap << endl;
      if (count == 0)
      {
        os << ap.ap_name();
      }
      else
      {
        os << ", " << ap.ap_name();
      }
      count++;
    }
    os << ";" << endl;
  }
  // translate ltlf formulas to FOL formulas
  string res = translate2pfol(f, "max $", c);
  os << res << ";" << endl;
}
// obtain a formula in backus normal form
inline formula
negate_formula(formula &r)
{
  if (r.kind() == op::Not)
  {
    return r[0];
  }
  else
  {
    return formula::unop(op::Not, r);
  }
}
/*-------------------------------------------------------------------*/
// translate ltlf formula to a past first order logic accepting reverse language
// make sure input is in backus/negation normal form
// @f the input formula interpreted as a past LTLf formula
// @t current position t
// @c unknown parameter for the code from Syft; it seems to be unused
/*-------------------------------------------------------------------*/
string translate2pfol(formula &f, string t, int &c)
{
  string curs, ts;
  string exs, alls;
  string res;
  int cur;
  int count;
  int intt;
  formula r, lft, rgt;

  // last position
  if (t == "max $")
  {
    ts = t;
    // then current index is 0
    intt = 0;
  }
  else
  {
    ts = "x" + t;
    // to integer number
    intt = stoi(t);
  }
  if (intt == 0)
  {
    ts = "max $";
  }

  // ts is current position = x

  switch (f.kind())
  {
  // Not operation
  // (¬folp(φ, x))
  case op::Not:
    res = "~(";
    r = f[0];
    res += translate2pfol(r, t, c);
    res += ")";
    break;
  // strong Next
  case op::strong_X:
    // X[!](0) = 0
    // X[!] (1) means that there must be a successor
    // Note: with finite semantics X[!](1)≠1.
    // next step is y = x(t+1)
    // X[!] φ => x: . φ ....
    // Y φ => .... φ .: x
    exs = "x" + to_string(intt + 1); // y ts = current position
    r = f[0];
    // current step t
    // ((∃exs)((exs = ts - 1) ∧ (exs >= 0) ∧ folp(φ, exs)))
    res = "(ex1 " + exs + ": ((" + exs + "=" + ts + "-1) & (" + ts + " > 0) & (";
    res += translate2pfol(r, to_string(intt + 1), c);
    res += ")))";
    break;
  // weak Next
  case op::X:
    //// X(1) = 1
    // X (0) means that there is no successor
    // We do not have X(0)=0 because that
    // is not true with finite semantics.
    // next step y = x(t+1)
    // current step t
    //((ts = 0) ∨ ((∃exs)((exs = ts - 1) ∧ (exs >= 0) ∧ folp(φ, exs))))
    r = f[0];
    exs = "x" + to_string(intt + 1);
    res = "((ex1 " + exs + ": (" + exs + "=" + ts + "-1 & (" + ts + " > 0) & (";
    res += translate2pfol(r, to_string(intt + 1), c);
    res += "))) | (" + ts + " = 0))";
    break;
  case op::F:
    // F φ => x: .... φ ....
    // P φ => .... φ ....: x
    // (∃y)((0 <= y <= x ) ∧ folp(φ, y)
    r = f[0];
    exs = "x" + to_string(intt + 1);
    res = "(ex1 " + exs + ": (" + exs + " <= " + ts + " & 0 <= " + exs + " & (";
    res += translate2pfol(r, to_string(intt + 1), c);
    res += ")))";
    break;
  //TODO
  case op::G:
    // G φ => x: φ...φ
    // G φ => φ...φ: x
    // G φ = ! (F !φ)
    r = f[0];
    // negate formula
    lft = negate_formula(r);
    rgt = formula::unop(op::F, lft);
    rgt = negate_formula(rgt);
    res = translate2pfol(rgt, t, c);
    break;
  // U corresponds to S
  case op::U:
    lft = f[0];
    rgt = f[1];
    // folp (ψ1 S ψ2 , x) = ((∃y)((0 ≤ y ≤ x) ∧ folp (ψ2 , y) ∧ (∀z)((y < z ≤ x) → folp (ψ1 , z))))
    // ψ2 holds before x and after that ψ1 holds until current position
    // ψ1 U ψ2 => x: ψ1 ψ1 ψ1 ... ψ1 ψ2 ....
    // ψ1 S ψ2 => .... ψ2 ψ1 .... ψ1 ψ1 ψ1: x
    exs = "x" + to_string(intt + 1);  // y
    alls = "x" + to_string(intt + 2); // z
    // (0 ≤ y ≤ x) ∧ folp (ψ2 , y)  ts = x current position
    res = "(ex1 " + exs + ": (" + exs + " <= " + ts + " & 0 <= " + exs + " & (";
    res += translate2pfol(rgt, to_string(intt + 1), c);
    // (y < z ≤ x) → folp (ψ1 , z)
    res += ") & (all1 " + alls + ": (" + alls + " <= " + ts + " & " + exs;
    res += " < " + alls + " => (";
    res += translate2pfol(lft, to_string(intt + 2), c);
    res += ")))))";
    break;
  //TODO a R b = !(!a U !b)
  case op::R: //New
    lft = f[0];
    rgt = f[1];
    lft = negate_formula(lft); // ! a
    rgt = negate_formula(rgt); // ! b
    rgt = formula::U(lft, rgt); // !a U !b
    rgt = negate_formula(rgt);  // !(!a U !b)
    res = translate2pfol(rgt, t, c);
    break;
  case op::Or:
    // list of Or operands
    //(fol(φ 1 , x) ∨ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + translate2pfol(child, to_string(intt), c) + ")";
      }
      else
      {
        res += " | (" + translate2pfol(child, to_string(intt), c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::And:
    // list of And operands
    // (fol(φ 1 , x) ∧ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + translate2pfol(child, to_string(intt), c) + ")";
      }
      else
      {
        res += " & (" + translate2pfol(child, to_string(intt), c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::tt:
    res += "(true)";
    break;
  case op::ff:
    res += "(false)";
    break;
  // atomic propositions
  case op::ap:
    res = "(";
    //ts = get_step_str(t);
    res += ts + " in ";
    //const std::string& str = f.ap_name();
    res += "" + f.ap_name();
    res += ")";
    break;
  default:
    cerr << "Formula: " << f << ". ";
    throw runtime_error("Error formula in translate2pfol()");
    exit(-1);
  }
  // cout<<res<<endl;
  return res;
}

void trans_prefixltlf2fol(ostream &os, formula &f)
{
  int c = 1;
  set<formula> aps;
  get_formula_aps(f, aps);
  // output the atomic propositions
  if (!aps.empty())
  {
    os << "m2l-str;" << endl;
    os << "var2 ";
    int count = 0;
    for (formula ap : aps)
    {
      //cout << "ap: " << ap << endl;
      if (count == 0)
      {
        os << ap.ap_name();
      }
      else
      {
        os << ", " << ap.ap_name();
      }
      count++;
    }
    os << ";" << endl;
  }
  // translate ltlf formulas to FOL formulas
  string res = get_prefix2fol(f, 0, c);
  os << "(ex1 yend: (yend <= max $ & (ex1 yy: (yend <= yy & yy<= max $) => (~(false))) & (" << res << ")))"
     << ";" << endl;
}

/*-------------------------------------------------------------------*/
// translate ltlf formula to a first order logic whose prefix satisfies the ltlf formula
// USE 'yend' to indicate the end of a word for the ltlf formula
// make sure input is in backus/negation normal form
// @f the input formula
// @t current position t
// @c unknown parameter for the code from Syft; it seems to be unused
/*-------------------------------------------------------------------*/
string
get_prefix2fol(formula &f, int t, int &c)
{
  string curs, ts;
  string exs, alls;
  string res;
  int cur;
  int count;
  formula r, lft, rgt;
  switch (f.kind())
  {
  // Not operation
  // (¬fol(φ, x, y))
  case op::Not:
    res = "~(";
    r = f[0];
    res += get_prefix2fol(r, t, c);
    res += ")";
    break;
  // strong Next
  case op::strong_X:
    // X[!](0) = 0
    // X[!] (1) means that there must be a successor
    // Note: with finite semantics X[!](1)≠1.
    // next step is y = x(t+1)
    exs = "x" + to_string(t + 1);
    // current step t
    ts = get_step_str(t);
    // ((∃y)((y = x + 1) ∧ (y <= yend) ∧fol(φ, y)))
    res = "(ex1 " + exs + ": (" + exs + "=" + ts + "+1 & " + exs + " <= yend" + " & (";
    r = f[0];
    res += get_prefix2fol(r, t + 1, c);
    res += ")))";
    break;
  // weak Next
  case op::X:
    //// X(1) = 1
    // X (0) means that there is no successor
    // We do not have X(0)=0 because that
    // is not true with finite semantics.
    // next step y = x(t+1)
    exs = "x" + to_string(t + 1);
    // current step t
    ts = get_step_str(t);
    //((x = last) ∨ ((∃y)((y = x + 1) ∧ (y <= yend) ∧ fol(φ, y))))
    res = "((ex1 " + exs + ": (" + exs + "=" + ts + "+1 & " + exs + " <= yend" + " & (";
    r = f[0];
    res += get_prefix2fol(r, t + 1, c);
    res += "))) | (" + ts + " = max $))";
    break;
  case op::F:
    exs = "x" + to_string(t + 1);
    ts = get_step_str(t);
    // (∃y)((x ≤ y ≤ last) ∧ fol(φ, y)
    res = "(ex1 " + exs + ": (" + ts + " <= " + exs + " & " + exs + " <= yend" + " & (";
    r = f[0];
    res += get_prefix2fol(r, t + 1, c);
    res += ")))";
    break;
  case op::G:
    alls = "x" + to_string(t + 1);
    ts = get_step_str(t);
    //(∀z)((x ≤ z ≤ last) → fol(φ, z)))
    res = "(all1 " + alls + ": ((" + ts + " <= " + alls + " & " + alls + " <= yend) => (";
    r = f[0];
    res += get_prefix2fol(r, t + 1, c);
    res += ")))";
    break;
  case op::U:
    exs = "x" + to_string(t + 1);
    alls = "x" + to_string(t + 2);
    ts = get_step_str(t);
    lft = f[0];
    rgt = f[1];
    // ((∃y)((x ≤ y ≤ last) ∧ fol(φ2 , y) ∧ (∀z)((x ≤ z < y) → fol(φ1 , z))))
    res = "(ex1 " + exs + ": (" + ts + " <= " + exs + " & " + exs + " <= yend & (";
    res += get_prefix2fol(rgt, t + 1, c);
    res += ") & (all1 " + alls + ": (" + ts + " <= " + alls + " & " + alls;
    res += " < " + exs + " => (";
    res += get_prefix2fol(lft, t + 2, c);
    res += ")))))";
    break;
  case op::R: //New
    exs = "x" + to_string(t + 1);
    alls = "x" + to_string(t + 2);
    ts = get_step_str(t);
    lft = f[0];
    rgt = f[1];
    // (((∃y)((x ≤ y ≤ last) ∧ fol(φ1 , y) ∧ (∀z)((x ≤ z ≤ y) → fol(φ 2 , z)))) ∨ ((∀z)((x ≤ z ≤ last) → fol(φ 2 , z))))
    res = "((ex1 " + exs + ": (" + ts + " <= " + exs + " & " + exs + " <= yend & (";
    res += get_prefix2fol(lft, t + 1, c);
    res += ") & (all1 " + alls + ": (" + ts + " <= " + alls + " & " + alls;
    res += " <= " + exs + " => (";
    res += get_prefix2fol(rgt, t + 2, c);
    res += ")))))";
    res += "| (all1 " + alls + ": ((" + ts + " <= " + alls + " & " + alls + " <= yend) => (";
    res += get_prefix2fol(rgt, t + 2, c);
    res += "))))";
    break;
  case op::Or:
    // list of Or operands
    //(fol(φ 1 , x) ∨ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + get_prefix2fol(child, t, c) + ")";
      }
      else
      {
        res += " | (" + get_prefix2fol(child, t, c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::And:
    // list of And operands
    // (fol(φ 1 , x) ∧ fol(φ 2 , x))
    count = 0;
    // size must larger than 2
    res += "(";
    for (formula child : f)
    {
      if (count == 0)
      {
        res += "(" + get_prefix2fol(child, t, c) + ")";
      }
      else
      {
        res += " & (" + get_prefix2fol(child, t, c) + ")";
      }
      count++;
      //cout << "subformula: " << child << endl;
    }
    res += ")";
    break;
  case op::tt:
    res += "(true)";
    break;
  case op::ff:
    res += "(false)";
    break;
  // atomic propositions
  case op::ap:
    res = "(";
    ts = get_step_str(t);
    res += ts + " in ";
    //const std::string& str = f.ap_name();
    res += "" + f.ap_name();
    res += ")";
    break;
  default:
    cerr << "Formula: " << f << ". ";
    throw runtime_error("Error formula in translate2fol()");
    exit(-1);
  }
  // cout<<res<<endl;
  return res;
}
// push negation into inner formulas
formula
push_not_in(formula &f)
{
  formula res;
  if (f.kind() == op::ap)
  {
    res = formula::unop(op::Not, f);
  }
  else if (f.is_ff())
  {
    res = f.tt();
  }
  else if (f.is_tt())
  {
    res = f.ff();
  }
  else
  {
    formula lft, rgt, l, r, res2;
    std::vector<formula> lst;
    switch (f.kind())
    {
    case op::Not:
      lft = f[0];
      res = get_nnf(lft);
      break;
    case op::strong_X:
      // !(X[!] a) = X !a
      lft = f[0];
      r = push_not_in(lft);
      res = formula::unop(op::X, r);
      break;
    case op::X:
      // !(X a) = X[!] !a
      lft = f[0];
      r = push_not_in(lft);
      res = formula::unop(op::strong_X, r);
      break;
    case op::G:
      // !(G a) = F !a
      lft = f[0];
      r = push_not_in(lft);
      res = formula::unop(op::F, r);
      break;
    case op::F:
      // !(F b) = G !b
      lft = f[0];
      r = push_not_in(lft);
      res = formula::unop(op::G, r);
      break;
    case op::U:
      // !(a U b) = !a R !b
      lft = f[0];
      rgt = f[1];
      l = push_not_in(lft);
      r = push_not_in(rgt);
      res = formula::R(l, r);
      break;
    // weak Until
    // φ W ψ ≡ (φ U ψ) ∨ G φ ≡ φ U (ψ ∨ G φ) ≡ ψ R (ψ ∨ φ)
    // !(p W q) = !(G p | (p U q)) = (F !p ) & (!p R !q)
    case op::W:
      l = f[0];
      r = f[1];
      l = push_not_in(l);
      r = push_not_in(r);
      // res2 = ! p;
      res2 = formula::unop(op::F, l); // F !l
      res = formula::R(l, r);
      res = formula::multop(op::And, {res2, res});
      break;
    case op::R:
      // !(a R b) = !a U !b
      lft = f[0];
      rgt = f[1];
      l = push_not_in(lft);
      r = push_not_in(rgt);
      res = formula::U(l, r);
      break;
    case op::And:
      // !(a & b & ... & c) = !a | !b | ... | !c
      for (formula child : f)
      {
        res = push_not_in(child);
        lst.push_back(res);
      }
      res = formula::multop(op::Or, lst);
      break;
    case op::Or:
      // !(a | b | ... | c) = !a & !b & ... & !c
      for (formula child : f)
      {
        res = push_not_in(child);
        lst.push_back(res);
      }
      res = formula::multop(op::And, lst);
      break;
    case op::Implies:
      // !(a -> b) = a & !b
      lft = f[0];
      rgt = f[1];
      l = get_nnf(lft);
      r = push_not_in(rgt);
      res = formula::multop(op::And, {l, r});
      break;
    case op::Equiv:
      // !(a <-> b) = !(a->b) | !(b->a) = (a | b) & (! a | !b)
      lft = f[0];
      rgt = f[1];
      l = formula::binop(op::Implies, lft, rgt);
      r = formula::binop(op::Implies, rgt, lft);
      res2 = formula::multop(op::And, {l, r});
      res = push_not_in(res2);
      break;
    default:
      throw runtime_error("Error formula in push_not_in()");
      exit(-1);
    }
  }
  return res;
}

// get negation normal form of a formula
formula
get_nnf(formula &f)
{
  formula res;
  if (f.kind() == op::ap || f.kind() == op::tt || f.kind() == op::ff)
  {
    res = f;
  }
  else
  {
    formula l = f[0], r, lft, rgt, res2;
    formula lt, rt;
    std::vector<formula> lst;
    switch (f.kind())
    {
    case op::Not:
      res = push_not_in(l);
      break;
    // weak Next
    // X p = !(X[!] !p)
    case op::X:
      r = get_nnf(l);
      res = formula::unop(op::X, r);
      break;
    // strong Next
    case op::strong_X:
      r = get_nnf(l);
      res = formula::unop(op::strong_X, r);
      break;
    case op::G:
      r = get_nnf(l);
      res = formula::unop(op::G, r);
      break;
    case op::F:
      r = get_nnf(l);
      res = formula::unop(op::F, r);
      break;
    case op::U:
      r = f[1];
      lft = get_nnf(l);
      rgt = get_nnf(r);
      res = formula::U(lft, rgt);
      break;
    // weak Until
    // φ W ψ ≡ (φ U ψ) ∨ G φ ≡ φ U (ψ ∨ G φ) ≡ ψ R (ψ ∨ φ)
    // p W q = G p | (p U q) = (! F !p ) | (p U q)
    case op::W:
      r = f[1];
      l = get_nnf(l);
      r = get_nnf(r);
      // res2 = ! p;
      res2 = formula::unop(op::G, l);
      res = formula::U(l, r);
      res = formula::multop(op::Or, {res2, res});
      break;
    case op::R:
      r = f[1];
      lft = get_nnf(l);
      rgt = get_nnf(r);
      res = formula::R(lft, rgt);
      break;
    case op::And:
      for (formula child : f)
      {
        l = get_nnf(child);
        lst.push_back(l);
      }
      res = formula::multop(op::And, lst);
      break;
    case op::Or:
      for (formula child : f)
      {
        l = get_nnf(child);
        lst.push_back(l);
      }
      res = formula::multop(op::Or, lst);
      break;
    case op::Implies:
      r = f[1];
      lft = formula::Not(l);
      rgt = formula::multop(op::Or, {lft, r});
      res = get_nnf(rgt);
      break;
    case op::Equiv:
      // a <-> b = (a->b) & (b->a)
      r = f[1];
      lt = formula::binop(op::Implies, l, r);
      lft = get_nnf(lt);
      rt = formula::binop(op::Implies, r, l);
      rgt = get_nnf(rt);
      res = formula::multop(op::And, {lft, rgt});
      break;
    default:
      cerr << "Formula: " << f << ". ";
      throw runtime_error("Error formula in get_nnf()");
      exit(-1);
    }
  }
  return res;
}

// get a formula in Boolean normal form (BNF)
formula
get_bnf(formula &f)
{
  formula res;
  // propositions or constant
  if (f.kind() == op::ap || f.kind() == op::tt || f.kind() == op::ff)
  {
    res = f;
  }
  else
  {
    formula l = f[0], r, res2, lft, rgt;
    vector<formula> lst;
    switch (f.kind())
    {
    case op::Not:
      r = get_bnf(l);
      res = negate_formula(r);
      break;
    case op::strong_X:
      r = get_bnf(l);
      res = formula::unop(op::strong_X, r);
      break;
    // only use strong X, negation is cheaper than or operation
    // X p = !(X[!] !p)
    case op::X:
      r = get_bnf(l);
      // get !l
      r = negate_formula(r);
      res = formula::unop(op::strong_X, r);
      res = formula::unop(op::Not, res);
      break;
    // Mona only allows for extential quantification
    // G p = !(F !p)
    case op::G:
      r = get_bnf(l);
      // get !l
      r = negate_formula(r);
      res = formula::unop(op::F, r);
      res = formula::unop(op::Not, res);
      break;
    // Finally
    case op::F:
      r = get_bnf(l);
      res = formula::unop(op::F, r);
      break;
    // Until
    case op::U:
      r = f[1];
      l = get_bnf(l);
      r = get_bnf(r);
      res = formula::U(l, r);
      break;
    // weak Until
    // φ W ψ ≡ (φ U ψ) ∨ G φ ≡ φ U (ψ ∨ G φ) ≡ ψ R (ψ ∨ φ)
    // p W q = G p | (p U q) = (! F !p ) | (p U q)
    case op::W:
      r = f[1];
      l = get_bnf(l);
      r = get_bnf(r);
      // res2 = ! p;
      res2 = res = negate_formula(l);
      res = formula::unop(op::F, res2);  // F ! p
      res = formula::unop(op::Not, res); // ! F (!p)
      r = formula::U(l, r);              // p U q
      res = formula::multop(op::Or, {res, r});
      break;
    // Release, to Until
    // l R r = !(!l U !r)
    case op::R:
      r = f[1];
      l = get_bnf(l);
      r = get_bnf(r);
      // !l
      l = negate_formula(l);
      // !r
      r = negate_formula(r);
      res = formula::U(l, r);
      res = formula::unop(op::Not, res);
      break;
    // l1 & l2 & ... = ! (!l1 | !l2 | ... )
    case op::And:
      for (formula child : f)
      {
        l = get_bnf(child);
        l = negate_formula(l);
        // ! l
        lst.push_back(l);
      }
      res = formula::multop(op::Or, lst);
      res = formula::unop(op::Not, res);
      break;
    case op::Or:
      for (formula child : f)
      {
        l = get_bnf(child);
        lst.push_back(l);
      }
      res = formula::multop(op::Or, lst);
      break;
    case op::Implies:
      l = formula::unop(op::Not, l);
      r = f[1];
      r = formula::multop(op::Or, {l, r});
      res = get_bnf(r);
      break;
    case op::Equiv:
      r = f[1];
      lft = formula::binop(op::Implies, l, r);
      rgt = formula::binop(op::Implies, r, l);
      res2 = formula::multop(op::And, {lft, rgt});
      res = get_bnf(res2);
      break;
    default:
      cerr << "Formula: " << f << ". ";
      throw runtime_error("Error formula in get_bnf()");
      exit(0);
    }
  }
  return res;
}
