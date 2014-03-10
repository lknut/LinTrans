include $(LAWA_HOME)/Makefile.common


SOURCEFILES  = $(wildcard *.cpp)
#TARGETS  = $(patsubst %.cpp,%,$(wildcard *.cpp))
TARGETS      = $(patsubst %.cpp,%,$(SOURCEFILES))

#before: LDFLAGS additions as below, but did not work
#LDFLAGS += -L$(LAWA_HOME) -llawa -lextensionsflens -lextensionssparsegrid 
EX_LDFLAGS += -bind_at_load
#ADDITONAL_FLAGS += -llawa -lextensionsflens -llapack -lf77blas -lcblas -latlas -lgsl -lgslcblas
ADDITONAL_FLAGS += -llawa -lextensionsflens -llapack -lf77blas -lcblas -latlas -lgsl -lgslcblas 

ifdef GMP_HOME
  CXXFLAGS += -DMPFR -I$(MPFR_HOME)/include
  LDFLAGS += -L$(MPFR_HOME)/lib -lgmpxx -lgmp -llawa -lextensionsflens 
endif

.Phony: print_vars 


$(TARGETS): 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)  $(SOURCEFILES) -o $@  $(ADDITONAL_FLAGS)
# das würde eine .o erstellen aus allen .cpp files
#$(TARGETS): %.o : %.cpp 
#$(TARGETS): 
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(SOURCEFILES) $< -o $@  $(ADDITONAL_FLAGS)
clean:
	$(RM) *.o
	$(RM) $(TARGETS)
	$(RMDIR) *.dSYM

ifneq ($(MAKECMDGOALS),clean)
$(MAKECMDGOALS): $(MAKECMDGOALS).cpp 
#print_vars 
#waveletgalerkinoptionpricer1d.cpp
#$(MAKECMDGOALS).cpp 
#print_vars
#%: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)  $(MAKECMDGOALS).cpp -o $@  $(ADDITONAL_FLAGS)
endif

#waveletgalerkinoptionpricer1d: print_vars  
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) waveletgalerkinoptionpricer1d.cpp -o waveletgalerkinoptionpricer1d $(ADDITONAL_FLAGS)

print_vars:
	@echo '(TARGETS) = $(TARGETS)'
	@echo '< = $<'
	@echo '@ = $@'
	@echo 'dollar 1= $(1)'
	@echo 'MAKECMDGOALS= $(MAKECMDGOALS)'
	@echo 'SOURCEFILES= $(SOURCEFILES)'
	@echo 'CXXFLAGS= $(CXXFLAGS)'
	@echo 'LDFLAGS= $(LDFLAGS)'
	@echo 'CPPFLAGS= $(CPPFLAGS)'
	@echo 'LIBS= $(LIBS)'
	@echo 'LDFLAGS= $(LDFLAGS)'
	@echo 'MPFR_HOME= $(MPFR_HOME)'
	@echo 'LAWA_HOME= $(LAWA_HOME)'






