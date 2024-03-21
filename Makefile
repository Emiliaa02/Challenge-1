CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations -I/home/emilia_farina/pacs-examples/Examples/include -I/home/emilia_farina/pacs-examples/Examples/include/muparser

EXEC     = main
LDFLAGS ?= -L/home/emilia_farina/pacs-examples/Examples/lib -Wl,-rpath=/home/emilia_farina/pacs-examples/Examples/lib/libmuparser.so.2
LIBS  ?= -lmuparser

all: $(EXEC)

%.o: %.cpp Gradientmethod.hpp Point.hpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

$(EXEC): main.o Gradientmethod.o Point.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

clean:
	$(RM) *.o

distclean: clean
	$(RM) *~
