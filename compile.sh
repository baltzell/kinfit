g++ test_MissingParticle.C -o test_MissingParticle.exe `root-config --cflags --glibs`
g++ test_InvMass.C -o test_InvMass.exe `root-config --cflags --glibs`
g++ test_4C.C -o test_4C.exe `root-config --cflags --glibs`
g++ test_MissAndInv.C -o test_MissAndInv.exe `root-config --cflags --glibs`

#g++ -c -o KinFitter.o KinFitter.h `root-config --cflags --glibs`
