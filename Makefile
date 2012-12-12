CC=g++
CCFLAGS=-D__STDC_LIMIT_MACROS -D__cplusplus -O3 -march=native
LDFLAGS=-lboost_program_options -lgmp -lgmpxx
SRCDIR=src
INC=-I./include

all: qme-ng
Exception.o: $(SRCDIR)/Exception.cpp
	$(CC) $(CCFLAGS) $< $(INC) -c
greenSizeHash.o: $(SRCDIR)/greenSizeHash.cpp
	$(CC) $(CCFLAGS) $< $(INC) -c
carquois.o: $(SRCDIR)/carquois.cpp Exception.o
	$(CC) $(CCFLAGS) $< $(INC) -c
principalExtension.o: $(SRCDIR)/principalExtension.cpp Exception.o
	$(CC) $(CCFLAGS) $< $(INC) -c
mutexplorator.o: $(SRCDIR)/mutexplorator.cpp carquois.o Exception.o
	$(CC) $(CCFLAGS) $< $(INC) -c
greenexplorator.o: $(SRCDIR)/greenexplorator.cpp greenSizeHash.o principalExtension.o Exception.o
	$(CC) $(CCFLAGS) $< $(INC) -c
greenfinder.o: $(SRCDIR)/greenfinder.cpp principalExtension.o Exception.o
	$(CC) $(CCFLAGS) $< $(INC) -c
mutexploratorSeq.o: $(SRCDIR)/mutexploratorSeq.cpp carquois.o Exception.o mutexplorator.o
	$(CC) $(CCFLAGS) $< $(INC) -c
rng.o: $(SRCDIR)/rng.c
	$(CC) $(CCFLAGS) $< $(INC) -c
naututil.o: $(SRCDIR)/naututil.c
	$(CC) $(CCFLAGS) $< $(INC) -c
nauty.o: $(SRCDIR)/nauty.c
	$(CC) $(CCFLAGS) $< $(INC) -c
nautil.o: $(SRCDIR)/nautil.c
	$(CC) $(CCFLAGS) $< $(INC) -c
nausparse.o : $(SRCDIR)/nausparse.c
	$(CC) $(CCFLAGS) $< $(INC) -c
naugraph.o : $(SRCDIR)/naugraph.c
	$(CC) $(CCFLAGS) $< $(INC) -c
nautinv.o : $(SRCDIR)/nautinv.c
	$(CC) $(CCFLAGS) $< $(INC) -c

qme-ng: qme-ng.cpp greenexplorator.o greenfinder.o mutexploratorSeq.o mutexplorator.o greenSizeHash.o principalExtension.o carquois.o Exception.o nautil.o rng.o nauty.o naututil.o nausparse.o naugraph.o nautinv.o
	$(CC) $(CCFLAGS) $(INC) $^ $(LDFLAGS) -o $@

clean:
	-rm -f *.o qme-ng *.cpp~ src/*.cpp~ include/*.hpp~ src/*.c~ include/*.h~
