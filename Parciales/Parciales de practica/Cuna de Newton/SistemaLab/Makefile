CC = g++
CXXFLAGS = -Wall -O3 -std=c++17

# Encuentra todos los archivos .cpp en el directorio
CPP_FILES := $(wildcard *.cpp)

# Elimina el sufijo .cpp de cada archivo para obtener el nombre base
TARGETS := $(patsubst %.cpp,%,$(CPP_FILES))

# Primer target predeterminado
all: $(TARGETS:=.gif)

# Reglas para convertir un .cpp en un .x (ejecutable)
%.x: %.cpp
	$(CC) $(CXXFLAGS) $< -o $@

# Reglas para ejecutar un .x y guardar su salida en un .txt
%.txt: %.x
	./$< > $@

# Reglas para convertir un .txt en un .gif usando gnuplot
%.gif: %.txt
	gnuplot $<

# Marca 'clean' como un target ficticio
.PHONY: clean

# Comando para limpiar todos los archivos generados
clean:
	rm -f *.x *.gp *.pdf *.txt *.gif
