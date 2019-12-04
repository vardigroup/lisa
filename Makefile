all: T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12

T1: 
	g++ explicit.cc mona.cc spotutil.cc dfwavar.cc dfwa.cc -o explicit -lspot -lbddx

T2: minimize.cc test_min.cc
	g++ test_min.cc minimize.cc -o test_min -lspot -lbddx
	
T3: test_eq.cc mona.cc spotutil.cc
	g++ test_eq.cc mona.cc spotutil.cc -o test_eq -lspot -lbddx
	
T4: test_dfwa.cc 
	g++ test_dfwa.cc dfwavar.cc spotutil.cc dfwa.cc mona.cc -o test_dfwa -lspot -lbddx

T5: test_synt.cc synt.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc
	g++ test_synt.cc synt.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc -o test_synt -lspot -lbddx

T6: ltlfsynt.cc synt.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc
	g++ ltlfsynt.cc strategy.cc synt.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc -o ltlfsynt -lspot -lbddx
	
T7: test_dfwanew.cc 
	g++ test_dfwanew.cc dfwavar.cc dfwanew.cc mona.cc dfwa.cc spotutil.cc -o test_dfwanew -lspot -lbddx

T8: 
	g++ test_dfwamin.cc dfwavar.cc minimize.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc -o test_dfwamin -lspot -lbddx -lcudd

T9: 
	g++ symboltrans.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc -o symbolictrans -lspot -lbddx -lcudd
	
T10: 
	g++ symboltrans2.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin2.cc dfwamin3.cc  -o symbolictrans2 -lspot -lbddx -lcudd -lsylvan

T11: 
	g++ symboltrans3.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin2.cc dfwamin3.cc  -o symbolictrans3 -lspot -lbddx -lcudd -lsylvan
	
T12: 
	g++ lisa.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc dfwamin3.cc  -o lisa -lspot -lbddx -lcudd -lsylvan
	
T13: 
	g++ lisanfa.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc dfwamin3.cc  -o lisanfa -lspot -lbddx -lcudd -lsylvan
#------------------------------------------------------
clean: #clean
	rm -f *.o main *.cc~ *.h~ Makefile~
#------------------------------------------------------

