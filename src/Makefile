CXXFLAGS := -Wall -g -std=c++11
CXXLIBS := -lm
SOURCE := main.cpp physics.cpp calculations.cpp

OBJS := $(patsubst %.cpp, %.o, $(SOURCE))
DEPS := $(patsubst %.o, %.d, $(OBJS))

release: $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o main $(CXXLIBS)

%.o: %.cpp %.d
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.d: %.cpp
	@set -e; $(CXX) $(GTKMM_FLAGS) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) $(PINEFM) $(OBJS) $(DEPS) *~

-include $(DEPS)
