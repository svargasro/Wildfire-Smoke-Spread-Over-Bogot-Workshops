import os
from PIL import Image

# Ruta de la carpeta que contiene los GIFs
carpeta_gifs = "."

# Obtener lista de archivos GIF en la carpeta y ordenarlos
gifs = sorted([f for f in os.listdir(carpeta_gifs) if f.endswith(".gif")])

# Cargar los GIFs como imágenes
imagenes = [Image.open(os.path.join(carpeta_gifs, gif)) for gif in gifs]

# Guardar el primer GIF con el resto de las imágenes como frames adicionales
imagenes[0].save("AdvectionDifusion.gif", save_all=True, append_images=imagenes[1:], duration=100, loop=0)

print("GIF combinado guardado como output.gif")
