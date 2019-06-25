export OPENGL=0
export OPENMP=0
export PROFILING=0
include Makefile.defs

all: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I/opt/local/include -I../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -L/opt/local/lib -lhdf5 -lhdf5_hl -o rebound
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
