CC = g++
CFLAGS = -Wall -O3 --std=c++11

Target = main
Objects = main.o functions.o

all: animation.gif

%.o: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

$(Target).x: $(Objects)
	$(CC) $(Objects) -o $(Target).x


animation.gif: $(Target).x
	./$(Target).x > $(Target).gp
	gnuplot $(Target).gp

.PHONY: clean

clean:
	rm -f $(Objects) $(Target)

