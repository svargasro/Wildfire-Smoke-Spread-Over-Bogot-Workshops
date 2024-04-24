CC = g++
CXXFLAGS = -Wall -O3 -std=c++17

# Encuentra todos los archivos .cpp en el directorio
CPP_FILES := CunaNewton.cpp

# Elimina el sufijo .cpp de cada archivo para obtener el nombre base
TARGETS := $(patsubst %.cpp,%,$(CPP_FILES))
TARGETS_POS := $(wildcard PenduloXvsY_K=*.txt)
COMMAND = ./$< {} > $(TARGETS)_K={}.txt
TIME = /usr/bin/time --format=$(FORMAT)
FORMAT = 'Real Time: %e\nUser Time: %U\nSystem Time: %S\nCPU usage: %P\nMax memory: %MkB\n'

# PARAMS := 0.1 
# PARAMS := 0.1 1 10
PARAMS := 0.1 0.2 0.5 1 2 5 10

TXT_FILES :=$(patsubst %,$(TARGETS)_K=%.txt,$(PARAMS))
# GP_FILES_POS := $(patsubst %.txt,%_graph.gp,$(TARGETS_POS))


# Primer target predeterminado
all: gifs pngs punto_c.png 
#clean1
# all: pngs 

# Reglas para convertir un .cpp en un .x (ejecutable)
$(TARGETS).x: $(TARGETS).cpp
	$(CC) $(CXXFLAGS) $< -o $@

# Regla para ejecutar un .x con parallel y guardar su salida en un .txt
# Genera archivos individuales para cada valor dt automáticamente
$(TXT_FILES): $(TARGETS).x
	$(TIME)	parallel -j 8 --line-buffer '$(COMMAND)' ::: $(PARAMS) >> Benchmark_Total.txt


# Regla para convertir todos los .txt en un .gif usando gnuplot
# Dependencias deben ser ajustadas para cada archivo específico generado
gifs: $(TXT_FILES)
	gnuplot $(TXT_FILES)

#Regla para convertir un .txt en un .png usando gnuplot
pngs: 
	gnuplot $(patsubst %.txt,%_graph.gp,$(TARGETS_POS))

#Regla para convertir los .txt generados en el numeral c
#en los respectivos .png usando gnuplot
punto_c.png: 
	gnuplot graph_c.gp

# Regla para limpiar todos excepto Benchmark_Pendulos.txt, los .gif y el .png
clean1:
	find . -type f \( -name "*.x" -o -name "*.txt" -o -name "*.gp"  \) -not -name "Benchmark_CunaNewton.txt" -not -name "graph_c.gp" -exec rm -f {} +


# Marca 'clean' como un target ficticio
.PHONY: clean

# Comando para limpiar todos los archivos generados
clean:
	find . -type f \( -name "*.x" -o -name "*.txt"  -o -name "*.gif" -o -name "*.png" \)  -not -name "graph_c.gp" -exec rm -f {} +

	