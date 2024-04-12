all: grafico4.1.png 4.1.txt gif4.2.gif 4.3.txt grafico4.3.png

4.1.txt: 4.1.out
	./4.1.out > 4.1.txt

grafico4.1.png: 4.1.gp 4.1.txt
	gnuplot 4.1.gp

4.1.out: 4.1_GasDeLennard-Jones.cpp
	g++ 4.1_GasDeLennard-Jones.cpp -o 4.1.out

gif4.2.gif: 4.2.gp
	gnuplot 4.2.gp

4.2.gp: 4.2.out
	./4.2.out > 4.2.gp

4.2.out: 4.2_GasDeLennard-Jones.cpp
	g++ 4.2_GasDeLennard-Jones.cpp -o 4.2.out

4.3.txt: 4.3.out
	./4.3.out > 4.3.txt

grafico4.3.png: 4.3.gp 4.3.txt
	gnuplot 4.3.gp

4.3.out: 4.3_GasDeLennard-Jones.cpp
	g++ 4.3_GasDeLennard-Jones.cpp -o 4.3.out

clean: 
	rm -f *.out *.txt *.x