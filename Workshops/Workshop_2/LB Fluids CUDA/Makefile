GCC = g++
CFLAGS = -O3 --std=c++17

TARGET = $(patsubst %.cpp,%,$(wildcard *.cpp))
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: graficar

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

graficar: $(TARGET).x
	./$<
	gnuplot grafica.gp
.PHONY: clean

clean:
	rm -f *.x *.dat 