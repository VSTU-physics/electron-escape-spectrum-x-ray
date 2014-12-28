CXXFLAGS := -Wall -g -std=c++11
CXXLIBS := -lm
GTKFLAGS := $(shell pkg-config --cflags gtk+-3.0)
GTKLIBS := $(shell pkg-config --libs gtk+-3.0)
SOURCE := spectrum.cpp plots.cpp physics.cpp calculations.cpp parse.cpp monte-carlo.cpp INIReader.cpp
TEST_S := spectrum_test.cpp $(SOURCE)
GUI_A_S := analytics.cpp $(SOURCE)

TEST_O := $(patsubst %.cpp, %.o, $(TEST_S))
TEST_D := $(patsubst %.o, %.d, $(TEST_O))
GUI_A_O := $(patsubst %.cpp, %.o, $(GUI_A_S))
GUI_A_D := $(patsubst %.o, %.d, $(GUI_A_O))

test: $(TEST_O) $(TEST_D)
	$(CXX) $(TEST_O) -o test $(CXXLIBS)

analytics: $(GUI_A_O) $(GUI_A_D)
	$(CXX) $(GUI_A_O) -o analytics $(CXXLIBS) $(GTKLIBS)

%.o: %.cpp %.d
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(GTKFLAGS)

%.d: %.cpp
	@set -e; $(CXX) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) *.o *.d *~

-include $(DEPS)
