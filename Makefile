all: T1

T1: 
	g++ lisa.cc minimize.cc dfwavar.cc dfwa.cc spotutil.cc mona.cc dfwamin.cc synt.cc strategy.cc dfwamin2.cc   -o lisa -lspot -lbddx -lcudd -lsylvan

#------------------------------------------------------
clean: #clean
	rm -f *.o main *.cc~ *.h~ Makefile~
#------------------------------------------------------

