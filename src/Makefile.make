# для сборки под Windows с mingw32-make
CXXFLAGS = -Wall -g -std=c++11
CXXLIBS = 
SOURCE = main.cpp physics.cpp calculations.cpp parse.cpp

OBJS = $(SOURCE:.cpp=.o)
DEPS = $(SOURCE:.cpp=.d)

release: $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o main.exe $(CXXLIBS)

%.o: %.cpp %.d
	$(CXX) -c -MD $< -o $@ $(CXXFLAGS)

%.d: %.cpp
	$(CXX) -c -MD $<

clean:
	$(RM) $(OBJS) $(DEPS)

-include $(DEPS)
