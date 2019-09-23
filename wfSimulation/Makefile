ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
#EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
#LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .cc,.C
.PREFIXES: ./obj/

#----------------------------------------------------#

runApp: runTest runUsbwclib 

runTest: runTest.cpp obj/wfSim.o obj/waveform.o
	g++ $(CXXFLAGS) -o runTest runTest.cpp $(OUTLIB)/*.o `root-config --cflags --glibs`

runUsbwclib: runUsbwclib.cpp obj/usbwclib.o
	g++ -o runUsbwclib runUsbwclib.cpp $(OUTLIB)/*.o

obj/wfSim.o: src/wfSim.cpp src/wfSim.hh src/waveform.cpp src/waveform.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)wfSim.o $<

obj/waveform.o: src/waveform.cpp src/waveform.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)waveform.o $<

obj/waveformSimple.o: src/waveformSimple.cpp src/waveformSimple.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)waveform.o $<

#obj/usbwclib.o: src/usbwclib.cpp src/usbwclib.hh
#	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)usbwclib.o $<

obj/usbwclib.o: src/usbwclib.cpp src/usbwclib.hh
	g++ -c -o ./obj/usbwclib.o src/usbwclib.cpp


clean:
	rm -f *~
	rm -f src/*~
	rm -f $(OUTLIB)*.o
	rm -f runTest
	rm -f runUsbwclib
