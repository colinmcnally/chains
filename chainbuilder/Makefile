export OPENGL=0
export OPENMP=0
export PROFILING=0
include Makefile.defs

all: 
	@echo ""
	@echo "Compiling Rebound problem file ..."
	$(CC) $(INC) -I/opt/local/include -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -L/opt/local/lib -lhdf5 -lhdf5_hl -o rebound
	$(CC) -c -I/opt/local/include -I../../src/  $(OPT) $(PREDEF) problem.c -fPIC -o chainbuilder.o
	$(CC) $(INC) -L. -lrebound $(LIB) -L/opt/local/lib -lhdf5 -lhdf5_hl -shared -o chainbuilder.so chainbuilder.o
	@echo ""
	@echo "REBOUND compiled successfully."

# replace with a link to librebound from local dir
#librebound: 
#	@echo "Compiling shared library librebound.so ..."
#	$(MAKE) -C ../../src/
#	@-rm -f librebound.so
#	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
