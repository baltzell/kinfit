
TESTS := test_4C test_MissingParticle test_InvMass test_MissAndInv

TESTURL := https://clasweb.jlab.org/clas12offline/distribution/coatjava/validation_files/eb/5.1-fid-r11/

ROOTFLAGS := `root-config --cflags --glibs`

CXX := g++

test:
	@test -d bin || mkdir -p bin
	@for t in $(TESTS); do\
		$(CXX) $$t.C -o bin/$$t.exe $< $(ROOTFLAGS); \
	done

clean:
	@rm -rf bin
	@rm -rf data
	@rm -f *.pdf

download:
	@test -d data || mkdir -p data
	@cd data && wget --progress=dot:mega -N --no-check-certificate $(TESTURL)/electronFTgamma.hipo

