GCC = g++
CFLAGS = -O3 --std=c++17


TARGET = LB_AdvectionDiffusion_D2Q9
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: ejecutar gif 

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

ejecutar: $(TARGET).x
	./$< 0.01 0.016	

gif:
	convert -delay 25 -loop 0 $(shell ls frames/*.png | sort -V) AdvectionDifusion.gif

.PHONY: clean

clean:
	rm *.x  *.gp
	rm -r frames data