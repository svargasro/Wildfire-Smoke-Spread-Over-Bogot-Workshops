all: run_graph

data_to_table: punto_d.cpp
	g++ punto_d.cpp -o data_to_table.x

run_table: data_to_table
	./data_to_table.x

run_graph: run_table
	gnuplot plot_script.gp

clean:
	rm -f data_to_table.x 
