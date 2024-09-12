GCC = g++
CFLAGS = -O3 --std=c++17
#S_glob = 0.001
D = 0.016

TARGET = LB_AdvectionDiffusion_D2Q9
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: gif clean

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $< 

ejecutar: $(TARGET).x
	./$< $(D)
	

gif: ejecutar
	convert -delay 25 -loop 0 $(shell ls frames/*.png | sort -V) AdvectionDifusion.gif

.PHONY: clean

clean:
	rm *.x  *.gp 
	rm -r frames data