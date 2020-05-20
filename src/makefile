TARG = ../SC

VPATH= ../

CPP = g++

CPPFLAGS = -O3 -W -Wall -I/opt/local/include/
#CPPFLAGS = -O3 -W -Wall --static
#CPPFLAGS = -g3 -W -Wall -pg -I/opt/local/include/

OFILES = datas.o additional_functions.o
 
OMAIN=  main_kg.cpp
LDFLAGS = -L/opt/local/lib/ -lgsl -lgslcblas
OMAIN1=  main_test.cpp


default : all

all : $(TARG) $(TARG1)

$(TARG) : $(OFILES) $(OMAIN)
	$(CPP) $(CPPFLAGS) -o $(TARG) $(OFILES) $(OMAIN) $(LDFLAGS)

$(TARG1) : $(OFILES) $(OMAIN1)
	$(CPP) $(CPPFLAGS) -o $(TARG1) $(OFILES) $(OMAIN1) $(LDFLAGS)

clean :
	\rm -f ./*.o $(TARG)  *.Rout

