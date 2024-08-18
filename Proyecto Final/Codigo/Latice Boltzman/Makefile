GCC = g++
CFLAGS = -O3 --std=c++17

TARGET = LB_AdvectionDifusion_D2Q9copy
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: ejecutar gif clean

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

ejecutar: $(TARGET).x
	./$<	

gif: 
	convert -delay 25 -loop 0 frames/*.png AdvectionDifusion.gif
.PHONY: clean

clean:
	rm -f *.x frames/*.png data/*.dat