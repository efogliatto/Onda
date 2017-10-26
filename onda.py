import numpy as np

import matplotlib.animation as animation

import matplotlib.pyplot as plt

import argparse




def vmedio( x, modo = 'uniforme'):

    """
    Velocidad c para la posicion x
    """
    
    c = 1.


    if modo == 'uniforme':

        c = 1.

    elif modo == 'transicion':

        c = 1.25 - 0.25 * np.tanh(40.*(0.75-x))

    elif modo == 'alien':

        c = 1.5 - np.exp(-300.*(x-1.75)**2)

        

    return c
    




if __name__ == "__main__":


    
    # Argumentos de consola
    
    parser = argparse.ArgumentParser(description='Resolución de la ecuación de onda en un dominio finito con condición radiativa de Sommerfeld')

    parser.add_argument('--dt', help='Paso de tiempo', type = float, default = 0.01)

    parser.add_argument('-n', help='Cantidad de puntos', type = int, default = 400)

    parser.add_argument('--tf', help='Tiempo final', type = float, default = 7.)

    parser.add_argument('--upwind', help='Upwind de primer orden en los bordes', action = 'store_true', dest='borde')

    parser.add_argument('-c','--velocidad', help='Perfil de velocidad del medio', dest = 'vel', default = 'uniforme', choices = ['uniforme','transicion','alien'])

    args = parser.parse_args()



    # Inicializacion de vectores y matrices para la solucion

    nt = int( args.tf / args.dt )

    V = np.zeros( (2*args.n-2, nt + 1) )

    M = np.zeros( (2*args.n-2, 2*args.n-2) )

    dx = 4. / (args.n - 1)

    
    
    # Condicion de borde en M
    
    if args.borde:
    
        M[0][0]     = -vmedio(0., modo = args.vel) / dx;

        M[0][1]     =  vmedio(0., modo = args.vel) / dx;

        M[args.n-1][args.n-2] =  vmedio(4., modo = args.vel) / dx;

        M[args.n-1][args.n-1] = -vmedio(4., modo = args.vel) / dx;


    else:

        M[0][0] = -3. * vmedio(0., modo = args.vel) / (2.*dx);
        
        M[0][1] =  4. * vmedio(0., modo = args.vel) / (2.*dx);

        M[0][2] = -1. * vmedio(0., modo = args.vel) / (2.*dx);

        M[args.n-1][args.n-3] = -1.* vmedio(4., modo = args.vel) / (2*dx);
        
        M[args.n-1][args.n-2] =  4.* vmedio(4., modo = args.vel) / (2*dx);

        M[args.n-1][args.n-1] = -3.* vmedio(4., modo = args.vel) / (2*dx);    




    # Coeficientes para u en el interior (du/dt = v)
    
    for i in range( 1, (args.n-1) ):

        M[i][i+args.n-1] = 1.



    # Coeficientes para v en el interior (dv/dt = c * d2u/d2x)
    
    for i in range(args.n-2):
        
        M[i+args.n][i]   = vmedio(i*dx, modo = args.vel) / (dx*dx)
        
        M[i+args.n][i+1] = -2 * vmedio(i*dx, modo = args.vel) / (dx*dx)
        
        M[i+args.n][i+2] = vmedio(i*dx, modo = args.vel) / (dx*dx)





    # Condicion inicial para u

    xdata = np.array( [i*dx for i in range(args.n)] )

    V[0:args.n,0] = np.exp( -200.*(xdata-0.25)**2 )


        
        
    # Resolucion con RK4
    
    k1 = np.zeros( (2*args.n-2, 1) )
    
    k2 = np.zeros( (2*args.n-2, 1) )

    k3 = np.zeros( (2*args.n-2, 1) )

    k4 = np.zeros( (2*args.n-2, 1) )
    
    
    for i in range(nt):

        k1 = args.dt * np.matmul(M, V[:,i])

        k2 = args.dt * np.matmul(M, ( V[:,i] + 0.5 * k1 ))

        k3 = args.dt * np.matmul(M, ( V[:,i] + 0.5 * k2 ))

        k4 = args.dt * np.matmul(M, ( V[:,i] + k3 ) )
        
        V[:,i+1] = V[:,i] + (1./6.) * k1 + (1./3.)*(k2+k3)  + (1./6.) * k4



    


    # fig, ax = plt.subplots()

    # # time_template = 'time = %.1f'
        
    # ax.set_xlim(0, 4)

    # ax.set_ylim(-0.5, 1.1)
        
    # # time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    # l1, = plt.plot([], [], '-', animated=True)


    # # plt.xticks([])

    # # plt.yticks([])
            

    # def update(frame, line1):            
    
    #     line1.set_data(xdata, V[0:args.n,frame])
            
    #     # ttext.set_text(time_template%(frame*args.dt))
                        
    #     return line1


        
    # ani = animation.FuncAnimation(fig, update, np.arange(1,nt), fargs = (l1), interval=1,  blit=True)



        
    # plt.show()

    
        
    
    pass

