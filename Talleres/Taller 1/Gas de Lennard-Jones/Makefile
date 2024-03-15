gcc -O3 flags.c -o flags.out

# Leer el número de potencias de 2 a calcular
read -p "Ingrese el numero de potencias de 2 a imprimir: " n

# Inicializar la potencia de 2 en 1 (2^0)
potencia=1

for ((i = 1; i <= n; i++))
do
    echo "2^$i = $potencia"
    # Calcular la siguiente potencia de 2 multiplicando por 2
    potencia=$((potencia * 2))

    # echo "n = 512  n_size = $potencia"
    # echo 512  $potencia | ./flags.out

    # echo "n = $((potencia * 2))  n_size = $potencia"
    # echo $potencia+1  $potencia | ./flags.out

    echo "n = $((potencia * 16))   n_size = $potencia"
    echo $((potencia * 16)) | ./flags.out
    # make run nsize
done
