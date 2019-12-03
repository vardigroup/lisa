
#ifdef DEBUG

void
trans_ltlf2fol(stringstream &ss, formula& f)
{
    // make sure f is negation normal form
    
}

/*-------------------------------------------------------------------*/
// translate ltlf formula to first order logic
/*-------------------------------------------------------------------*/
string trans_fol(formula& f, int t, int& c){
  string curs, ts;
  string exs, alls;
  string res;
  int cur;
  switch(f.kind())
  {
        case op::Not:
          res = "~(";
          formula r = f[0]
          res += trans_fol(r, t, c);
          res += ")";
          break;
        case op::X:
          exs = "x" + to_string(t+1);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x" + to_string(t);
          formula r = f[0]
          res = "(ex1 " + exs + ": ("+ exs + "=" + ts + "+1 & (";
          res += trans_fol(r, t+1, c);
          res += ")))";
          break;
        // weak Next is always replace with negation of next
        /*
        case eWNEXT:
          exs = "x"+to_string(t+1);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x"+to_string(t);
          res = "((ex1 "+exs+": ("+exs+"="+ts+"+1 & (";
          res += trans_fol(root->_right, t+1, c);
          res += "))) | ("+ts+" = max $))";
          break;
        */
        case op::F:
          exs = "x"+to_string(t+1);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x"+to_string(t);
          res = "(ex1 "+exs+": ("+ts+" <= "+exs+" & (";
          formula r = f[0]
          res += trans_fol(r, t+1, c);
          res += ")))";
          break;
        case op::G:
          alls = "x"+to_string(t+1);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x"+to_string(t);
          res = "(all1 "+alls+": (("+ts+" <= "+alls+") => (";
          formula r = f[0]
          res += trans_fol(r, t+1, c);
          res += ")))";
          break;
        case op::U:
          exs = "x"+to_string(t+1);
          alls = "x"+to_string(t+2);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x"+to_string(t);
          formula lft = f[0], rgt = f[1];
          res = "(ex1 "+exs+": ("+ts+" <= "+exs+" & (";
          res += trans_fol(rgt, t+1, c);
          res += ") & (all1 "+alls+": ("+ts+" <= "+alls+" & "+alls;
          res += " < "+exs+" => (";
          res += trans_fol(lft, t+2, c);
          res += ")))))";
          break;
        case op::R: //New
          exs = "x"+to_string(t+1);
          alls = "x"+to_string(t+2);
          if (t == 0)
            ts = to_string(t);
          else
            ts = "x"+to_string(t);
          formula lft = f[0], rgt = f[1];
          res = "((ex1 "+exs+": ("+ts+" <= "+exs+" & (";
          res += trans_fol(lft, t+1, c);
          res += ") & (all1 "+alls+": ("+ts+" <= "+alls+" & "+alls;
          res += " <= "+exs+" => (";
          res += trans_fol(rgt, t+2, c);
          res += ")))))";
          res += "| (all1 "+alls+": (("+ts+" <= "+alls+" & "+alls+" <= max $) => (";
          res += trans_fol(rgt, t+2, c);
          res += "))))";
          break;
        case op::Or:
        // list of Or operands
          res += "(("+trans_fol(root->_right, t, c);
          res += ") | (";
          res += trans_fol(root->_left, t, c)+"))";
          break;
        case op::And:
        // list of And operands
          res += "(("+trans_fol(root->_right, t, c);
          res += ") & (";
          res += trans_fol(root->_left, t, c)+"))";
          break;
        case op::tt:
          res += "(true)";
          break;
        case op::ff:
          res += "(false)";
          break;
        case 3:
          if (t == 0)
            ts = "("+to_string(t);
          else
            ts = "(x"+to_string(t);
          res += ts+" in ";
          res += alphabet_no_comma(root);
          res +=")";
          break;
        default:
          break;
  }
  // cout<<res<<endl;
  return res;
}

#endif
