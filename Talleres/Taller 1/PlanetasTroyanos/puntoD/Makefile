all: plot.pdf

%.x: %.o
	g++ $^ -o $@

data.txt: punto4.x
	./punto4.x

plot.pdf: data.txt plot.py
	python3 plot.py
	open plot.pdf

# plot.pdf: data.txt plot.gp
# 	gnuplot plot.gp
# 	open plot.pdf

clean:
	rm -f *.o *.x
