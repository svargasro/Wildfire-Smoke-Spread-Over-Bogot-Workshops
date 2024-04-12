all: plot.pdf

%.x: %.o
	g++ $^ -o $@

data.txt: punto3.x
	./punto3.x

plot.pdf: data.txt plot.py
	python3 plot.py

# plot.pdf: data.txt plot.gp
# 	gnuplot plot.gp
# 	open plot.pdf

clean:
	rm -f *.o *.x *.txt
