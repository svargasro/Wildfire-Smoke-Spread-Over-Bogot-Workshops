import os
from PIL import Image

# Ruta de la carpeta que contiene los GIFs
carpeta_gifs = "."

# Función para extraer el número al final del nombre del archivo
def obtener_numero(archivo):
    return int(archivo.split('_')[-1].split('.')[0])

# Obtener lista de archivos GIF en la carpeta y ordenarlos por el número al final del nombre
gifs = sorted([f for f in os.listdir(carpeta_gifs) if f.endswith(".gif")], key=obtener_numero)

# Cargar los GIFs como imágenes
imagenes = [Image.open(os.path.join(carpeta_gifs, gif)) for gif in gifs]

# Guardar el primer GIF con el resto de las imágenes como frames adicionales
imagenes[0].save("AdvectionDifusion.gif", save_all=True, append_images=imagenes[1:], duration=100, loop=0)

print("GIF combinado guardado como AdvectionDifusion.gif")
