CXXFLAGS := -Wall -g -std=c++11
CXXLIBS := -lm
SOURCE := spectrum.cpp plots.cpp physics.cpp calculations.cpp parse.cpp monte-carlo.cpp INIReader.cpp

OBJS := $(patsubst %.cpp, %.o, $(SOURCE))
DEPS := $(patsubst %.o, %.d, $(OBJS))

release: $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o spectrum $(CXXLIBS)

%.o: %.cpp %.d
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.d: %.cpp
	@set -e; $(CXX) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) $(PINEFM) $(OBJS) $(DEPS) *~

-include $(DEPS)
