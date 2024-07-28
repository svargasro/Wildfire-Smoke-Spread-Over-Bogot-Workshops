GCC = g++
CFLAGS = -O3 --std=c++17

TARGET = $(patsubst %.cpp,%,$(wildcard *.cpp))
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: ejecutar graficar gif

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

ejecutar: $(TARGET).x
	./$<	

graficar: 
	gnuplot grafica.gp

gif: 
	mkdir -p frames
	gnuplot frames.gp
	convert -delay 25 -loop 0 frames/*.png AdvectionDifusion.gif
	rm -rf frames
.PHONY: clean

clean:
	rm -f *.x *.dat 