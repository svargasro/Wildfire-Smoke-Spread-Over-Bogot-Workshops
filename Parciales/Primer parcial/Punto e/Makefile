CC = g++
CXXFLAGS = -Wall -O3 -std=c++17

# Encuentra todos los archivos .cpp en el directorio
CPP_FILES := CunaNewton.cpp

# Elimina el sufijo .cpp de cada archivo para obtener el nombre base
TARGETS := $(patsubst %.cpp,%,$(CPP_FILES))

# Define los valores de dt para usar en simulaciones
#PARAMS := 0.1 
PARAMS := 0.1 0.2 0.5 1 2 5 10

TXT_FILES :=$(patsubst %,$(TARGETS)_K=%.txt,$(PARAMS))
TARGETS_MODS := $(patsubst %.txt,%,$(TXT_FILES))


# Primer target predeterminado
#all: $(TARGETS_MODS).gif pendulo.png punto_c.png punto_e.png 
all: punto_e.png  clean1

# Reglas para convertir un .cpp en un .x (ejecutable)
$(TARGETS).x: $(TARGETS).cpp
	$(CC) $(CXXFLAGS) $< -o $@

# Regla para ejecutar un .x con parallel y guardar su salida en un .txt
# Genera archivos individuales para cada valor dt automáticamente
$(TXT_FILES): $(TARGETS).x
	parallel -j 4 --line-buffer './$< {} > $(TARGETS)_K={}.txt' ::: $(PARAMS)


# Regla para convertir todos los .txt en un .gif usando gnuplot
# Dependencias deben ser ajustadas para cada archivo específico generado
%.gif: $(TXT_FILES)
	gnuplot $(TXT_FILES)

#Regla para convertir un .txt en un .png usando gnuplot
pendulo.png: graph.gp $(TXT_FILES)
	gnuplot graph.gp


#Regla para convertir los .txt generados en el numeral e
#en los respectivos .png usando gnuplot
punto_e.png: graph_e.gp $(TXT_FILES)
	gnuplot graph_e.gp


# Regla para limpiar todos excepto Benchmark_Pendulos.txt, los .gif y el .png
clean1:
	find . -type f \( -name "*.x" -o -name "*.txt" \) -not -name Benchmark_CunaNewton.txt -exec rm -f {} +


# Marca 'clean' como un target ficticio
.PHONY: clean

# Comando para limpiar todos los archivos generados
clean:
	rm -f *.x *.txt *.gif *.png

	