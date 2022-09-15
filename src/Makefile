all: MI00.so MI01.so vls.so anc_hiv.so hiv_anc.so
	# Comment here preserves the prior tab
MI00.so: MI00.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('MI00.cpp','-O0 -g')"
MI01.so: MI01.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('MI01.cpp','-O0 -g')"
vls.so: vls.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('vls.cpp','-O0 -g')"

anc_hiv.so: anc_hiv.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('anc_hiv.cpp','-O0 -g')"

hiv_anc.so: hiv_anc.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('hiv_anc.cpp','-O0 -g')"

clean:
	rm -rf *o
