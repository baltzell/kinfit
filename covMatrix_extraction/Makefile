CFILES := covMatrix_single_particle covMatrix_pip_pim_test

ROOTFLAGS := `root-config --cflags --glibs`

HIPOFLAGS := -I$(HIPO)/hipo4 -L$(HIPO)/slib -lhipo4 

CXXFLAGS := $(ROOTFLAGS) $(HIPOFLAGS)

CXX := g++

test:
	@test -d bin || mkdir -p bin
	@for t in $(CFILES); do\
		$(CXX) $$t.C -o bin/$$t.exe $< $(CXXFLAGS); \
	done

clean:
	@rm -rf bin
	@rm -rf data
	@rm -f *.so
	@rm -f *.d
	@rm -f *.pcm
	@rm -f *.pdf
