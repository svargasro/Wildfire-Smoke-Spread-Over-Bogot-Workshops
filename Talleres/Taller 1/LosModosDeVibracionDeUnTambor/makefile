all: final

final: datos_a Grafica_2a datos_b Grafica_2b Grafica_2c Grafica_2d clean

Grafica_2d: Ejercicio_2d.py
	python3 Ejercicio_2d.py

Grafica_2c: datos_b.dat Ejercicio_2c.py
	python3 Ejercicio_2c.py

Grafica_2b: datos_b.dat plot_2b.py
	python3 plot_2b.py

datos_b: Ejercicio_2b.cpp
	g++ Ejercicio_2b.cpp -o Ejercicio_2b.exe
	./Ejercicio_2b.exe 
	
Grafica_2a: datos_a.dat plot_2a.py
	python3 plot_2a.py

datos_a: Ejercicio_2a.cpp
	g++ Ejercicio_2a.cpp -o Ejercicio_2a.exe
	./Ejercicio_2a.exe 
clean:  
	rm -f *.exe *.dat
