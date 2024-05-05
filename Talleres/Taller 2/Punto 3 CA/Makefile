GCC = g++
CFLAGS = -O3 --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = $(patsubst %.cpp,%,$(wildcard *.cpp))

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(FSANITIZE) $(TARGET).cpp -o $(TARGET).out

$(TARGET).txt: $(TARGET).out
	./$(TARGET).out > $(TARGET).txt

$(TARGET).png: $(TARGET).txt $(TARGET).py
	python $(TARGET).py


clean: 
	rm -f *.out *.txt *.x