import ODE

import numpy as np

import matplotlib.animation as animation

import argparse




if __name__ == "__main__":


    
    # Argumentos de consola
    
    parser = argparse.ArgumentParser(description='Resolución de la ecuación de onda en un dominio finito con condición radiativa de Sommerfeld')

    parser.add_argument('--dt', help='Paso de tiempo', type = float, default = 0.01)

    parser.add_argument('-n', help='Cantidad de puntos', type = int, default = 400)

    parser.add_argument('--upwind', help='Upwind de primer orden en los bordes', action = 'store_true', dest='borde')

    args = parser.parse_args()



    # Inicializacion de vectores y matrices para la solucion

    V = np.zeros( (2*args.n-2, 1) )

    M = np.zeros( (2*args.n-2, 2*args.n-2) )


    
    
    # Condicion de borde en M

    if args.borde:
    
        M[0][0] = -C/dx;
        M[0][1] =  C/dx;
        M[N-1][N-2] =  C/dx;
        M[N-1][N-1]   = -C/dx;


# #   Condiciones de borde para u, en x=0. Segundo Orden
#   M(1,1) = -3*C/(2*dx);
#   M(1,2) =  4*C/(2*dx);
#   M(1,3) = -1*C/(2*dx);

#  ###  Condiciones de borde para u, en x=4. Segundo Orden
#   M(N,N-2) = -1*C/(2*dx);
#   M(N,N-1) =  4*C/(2*dx);
#   M(N,N)   = -3*C/(2*dx);    


    
    pass

