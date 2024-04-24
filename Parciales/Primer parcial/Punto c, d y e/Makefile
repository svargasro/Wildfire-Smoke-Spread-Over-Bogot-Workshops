CC = g++
CXXFLAGS = -Wall -O3 -std=c++17

# Encuentra todos los archivos .cpp en el directorio
CPP_FILES := $(wildcard *.cpp)

# Elimina el sufijo .cpp de cada archivo para obtener el nombre base
TARGETS := $(patsubst %.cpp,%,$(CPP_FILES))

# Define los valores de dt para usar en simulaciones
PARAMS := 1.0

TXT_FILES :=$(patsubst %,$(TARGETS)_K=%.txt,$(PARAMS))
TARGETS_MODS := $(patsubst %.txt,%,$(TXT_FILES))
BENCHMARK_FILES = $(patsubst %,Benchmark_CunaNewton_K=%e+10.txt,$(PARAMS))

$(info $(BENCHMARK_FILES))

# Primer target predeterminado
all: $(TARGETS_MODS).gif pendulo.png clean1

# Reglas para convertir un .cpp en un .x (ejecutable)
$(TARGETS).x: $(TARGETS).cpp
	$(CC) $(CXXFLAGS) $< -o $@

# Regla para ejecutar un .x con parallel y guardar su salida en un .txt
# Genera archivos individuales para cada valor dt automáticamente
$(TARGETS_MODS).txt: $(TARGETS).x
	parallel -j 8 --line-buffer './$< {} > $(TARGETS)_K={}.txt' ::: $(PARAMS)


# Regla para convertir todos los .txt en un .gif usando gnuplot
# Dependencias deben ser ajustadas para cada archivo específico generado
%.gif: $(TXT_FILES)
	gnuplot $(TXT_FILES)

#Regla para convertir los .txt generados en el numeral c
#en los respectivos .png usando gnuplot
punto_c.png: graph_c.gp
	gnuplot graph_c.gp

#Regla para convertir los .txt generados en el numeral d
#en los respectivos .png usando gnuplot
punto_d.png: graph_d.gp
	gnuplot graph_d.gp

#Regla para convertir los .txt generados en el numeral e
#en los respectivos .png usando gnuplot
punto_e.png: graph_e.gp
	gnuplot graph_e.gp

# Regla para limpiar todos excepto Benchmark_Pendulos.txt, los .gif y el .png
clean1:
	find . -type f \( -name "*." -o -name "*.x" -o -name "*.txt" -o -name "*.gp" \) -not -name $(BENCHMARK_FILES) -exec rm -f {} +


# Marca 'clean' como un target ficticio
.PHONY: clean

# Comando para limpiar todos los archivos generados
clean:
	rm -f *.x *.txt *.gp *.gif *.png
#find . -type f \( -name "*.x" -o -name "*.txt" -o -name "*.gp" -o -name "*.png" -o -name "*.gif" \) -not -name "Benchmark_Pendulos.txt" -exec rm -f {} +	
	