all: plot.png

%.x: %.o
	g++ $^ -o $@

data.txt: punto4.x
	./punto4.x

plot.png: data.txt plot.py
	python3 plot.py


# plot.png: data.txt plot.gp
# 	gnuplot plot.gp
# 	open plot.png

clean:
	rm -f *.o *.x *.txt
