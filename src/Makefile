all: MI00.so MI01.so
	# Comment here preserves the prior tab
MI00.so: MI00.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('MI00.cpp','-O0 -g')"
MI01.so: MI01.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('MI01.cpp','-O0 -g')"

clean:
	rm -rf *o
