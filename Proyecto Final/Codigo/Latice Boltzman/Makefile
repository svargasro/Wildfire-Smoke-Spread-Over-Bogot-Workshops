GCC = g++
CFLAGS = -O3 --std=c++17


TARGET = LB_AdvectionDifusion_D2Q9
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: ejecutar gif clean

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

ejecutar: $(TARGET).x
	./$<	

gif:
	convert -delay 25 -loop 0 $(ls frames/*.png | sort -V) AdvectionDifusion.gif 
	convert -delay 25 -loop 0 $(ls frames/*.png | sort -V) AdvectionDifusion.gif

.PHONY: clean

clean:
	rm -f *.x  *.gp
	rm -f -r frames data