# для сборки под Windows с mingw32-make
CXXFLAGS = -Wall -g -std=c++11 -lm
CXXLIBS = 
SOURCE = main.cpp physics.cpp calculations.cpp parse.cpp
DFLAGS = -D__WIN__

OBJS = $(SOURCE:.cpp=.o)
DEPS = $(SOURCE:.cpp=.d)

release: $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o main.exe $(CXXLIBS)

%.o: %.cpp %.d
	$(CXX) $(DFLAGS) -c -MD $< -o $@ $(CXXFLAGS)

%.d: %.cpp
	$(CXX) $(DFLAGS) -c -MD $<

clean:
	$(RM) $(OBJS) $(DEPS)

-include $(DEPS)
