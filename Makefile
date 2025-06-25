CFLAGS= -Wall -m64 -g -w
CXX=g++
ILOG=/opt/ibm/ILOG/CPLEX_Studio2211
CPPFLAGS= -DIL_STD -I$(ILOG)/cplex/include -I$(ILOG)/concert/include
CPLEXLIB=-L$(ILOG)/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread

comp-cplex:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o cs NewCutandSolve.cpp CuttingPlane.cpp SolveSparseProblem.cpp LocalSearch.cpp minknap.cpp Functions.cpp  $(CPLEXLIB) -ldl
run:comp-pc
	./cs
clean:
	rm -f  *.out *.aux *.log *.nav *.snm *.out *.toc 

