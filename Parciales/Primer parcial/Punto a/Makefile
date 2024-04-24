CC = g++
CXXFLAGS = -Wall -O3 -std=c++17

# Encuentra todos los archivos .cpp en el directorio
CPP_FILES := $(wildcard *.cpp)

# Elimina el sufijo .cpp de cada archivo para obtener el nombre base
TARGETS := $(patsubst %.cpp,%,$(CPP_FILES))

# Primer target predeterminado
all: $(TARGETS:=.gif) pendulo.png clean1

# Reglas para convertir un .cpp en un .x (ejecutable)
%.x: %.cpp
	$(CC) $(CXXFLAGS) $< -o $@

# Reglas para ejecutar un .x y guardar su salida en un .txt
%.txt: %.x
	parallel -j 1 --line-buffer './$<' ::: {1..10} > $@

# Reglas para convertir un .txt en un .gif usando gnuplot
%.gif: %.txt
	gnuplot $<
	
#Regla para convertir un .txt en un .png usando gnuplot
pendulo.png: graph.gp
	gnuplot graph.gp

# Regla para limpiar todos excepto Benchmark_Pendulos.txt, los .gif y el .png
clean1:
	find . -type f \( -name "*.x" -o -name "*.txt" -o -name "*.gp" \) -not -name "Benchmark_Pendulos.txt" -exec rm -f {} +


# Marca 'clean' como un target ficticio
.PHONY: clean

# Comando para limpiar todos los archivos generados
clean:
	rm -f *.x *.txt *.gp *.gif *.png
#find . -type f \( -name 
	