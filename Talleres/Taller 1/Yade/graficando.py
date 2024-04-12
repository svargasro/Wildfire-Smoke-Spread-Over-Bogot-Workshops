import pandas as pd
import matplotlib.pyplot as plt
import os  # Importa la biblioteca os para trabajar con archivos del sistema
import numpy as np 
# Initialize a figure and axis
fig, ax = plt.subplots()

# Lista de frecuencias corregida
W = [0.0025,0.003,0.0035,0.004,0.0045, 0.005, 0.0055]

# Iterate over frequency values (frec from 1 to 10)
for frec in W:
    # Verifica si el archivo CSV existe antes de intentar leerlo
    filename = f'Data/reloj20radio-{frec}.csv'
    if os.path.exists(filename):
        df = pd.read_csv(filename)
        
        # Convert DataFrame columns to numpy arrays
        t_values = df['t'].values
        z_values = df['z'].values
        
        # Plot the data
        ax.plot(t_values, z_values, label=f'radio {frec}')
        filtered_flows = z_values[z_values > 0.008]
        if len(filtered_flows) > 0:
            average_flow = np.mean(filtered_flows)
            print(f"Flujo promedio para el radio {frec}: {average_flow}")
        else:
            print(f"No hay valores en z_values mayores que 0.01 para el radio {frec}")
        
    else:
        print(f"File '{filename}' not found.")

# Set labels and title
ax.set_xlabel('t(s)')
ax.set_ylabel('W(kg/s)')
ax.set_title('Plot de W vs t para diferentes radios')

# Add legend
ax.legend()

# Show the plot
plt.show()

