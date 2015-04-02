CXXFLAGS := -Wall -g -std=c++11
CXXLIBS := -lm

SOURCE := spectrum.cpp plots.cpp physics.cpp calculations.cpp parse.cpp INIReader.cpp
GUI_A_S := analytics.cpp $(SOURCE)
GUI_M_S := mc.cpp $(SOURCE)

GUI_A_O := $(patsubst %.cpp, %.o, $(GUI_A_S))
GUI_A_D := $(patsubst %.o, %.d, $(GUI_A_O))

GUI_M_O := $(patsubst %.cpp, %.o, $(GUI_M_S))
GUI_M_D := $(patsubst %.o, %.d, $(GUI_M_O))

analytics: $(GUI_A_O) $(GUI_A_D)
	$(CXX) $(GUI_A_O) -o analytics $(CXXLIBS)

mc: $(GUI_M_O) $(GUI_M_D)
	$(CXX) $(GUI_M_O) -o monte-carlo $(CXXLIBS)

%.o: %.cpp %.d
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.d: %.cpp
	@set -e; $(CXX) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) *.o *.d *~

-include $(DEPS)
