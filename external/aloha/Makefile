
OBJS := $(patsubst %.c,%.o,$(wildcard *.c))
OBJS2 := $(patsubst %.cc,%.o,$(wildcard *.cc))

CXX=g++
CXXFLAGS=-Wall -pedantic -ffast-math -O3

libaloha.a: $(OBJS) $(OBJS2)
	ar rcs $@ $(OBJS) $(OBJS2)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f libaloha.a
