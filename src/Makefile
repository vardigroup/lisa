all: T1 T2

T1: 
	g++ lisa.cc  dfwavar.cc dfwa.cc spotutil.cc ltlf2fol.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc   -o lisa -lspot -lbddx -lcudd -O3
	
T2:
	g++ testltlf2fa.cc aututil.cc ltlf2fol.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc  -o testltlf2fol -lspot -lbddx -lcudd

#------------------------------------------------------
clean: #clean
	rm -f *.o main *.cc~ *.h~ Makefile~
#------------------------------------------------------

