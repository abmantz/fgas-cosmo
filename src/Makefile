CPP = g++
FLAGS = -O3

FLAGS += -fPIC
#FLAGS += -arch i386 # mac-specific

OBJECTS =  clusters.o fgas.o lensing.o wrapper.o xray.o ConfigFile.o
SEPOBJ = fwrapper.o

ARCHIVES = libfgas.a

default: objects archives

listobjects:
	@echo $(OBJECTS)

objects: $(OBJECTS) $(SEPOBJ)

archives: $(ARCHIVES)

$(ARCHIVES): $(OBJECTS)

clusters.o: clusters.hpp util.hpp
fgas.o: clusters.hpp fgas.hpp lensing.hpp util.hpp xray.hpp
lensing.o: clusters.hpp lensing.hpp  util.hpp ConfigFile.h
wrapper.o: clusters.hpp fgas.hpp lensing.hpp wrapper.hpp xray.hpp
xray.o: clusters.hpp util.hpp xray.hpp



%.o: %.cpp
	$(CPP) -o $@ -c $(FLAGS) $(INCLUDE) $*.cpp

%.a:
	rm -f $@; ar -r $@ $(OBJECTS)

clean:
	rm -f $(EXE) $(OBJECTS) $(SEPOBJ) $(ARCHIVES)
