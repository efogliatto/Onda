import numpy as np


def theta( myfunc, y0, h, nt, **kwargs ):

    """
    Resuelve la ode definida por myfunc (f') usando el metodo theta.

    Argumentos
    h: paso de tiempo
    nt: numero de intervalos temporales
    th: theta

    Devuelve un array con la solucion con todos los pasos de tiempo
    """

    y = np.zeros( (int(nt)+1, len(y0)))

    y[0] = y0


    th = 0.

    for key in kwargs:

        if key == 'th':

            th = kwargs[key]


    

    
    
    # Explicit Euler

    if th == 0:

        for i in range(nt):

            y[i+1] = y[i] + h * myfunc(y[i], **kwargs)



    # Other implicit
            
    else:
        
        for i in range(nt):

            # dif = 1.

            # y_aux = y[i]


            # while (dif > 1.e-15):

            #     y[i+1] = y[i]   +   h * (1-th) * myfunc(y[i],g,l)   +   h * th * myfunc(y_aux,g,l)

            #     dif = np.linalg.norm( (y[i+1]-y_aux) )

            #     y_aux = y[i+1]
                
            
            for k in range(50):
                            
                y[i+1] = y[i]   +   h * (1-th) * myfunc(y[i], **kwargs)   +   h * th * myfunc(y[i+1], **kwargs)
            
        





    return y


    


def RK4( myfunc, y0, h, nt, **kwargs ):


    """
    Resuelve la ode definida por myfunc (f') usando el metodo theta.

    Argumentos
    h: paso de tiempo
    nt: numero de intervalos temporales
    th: theta

    Devuelve un array con la solucion con todos los pasos de tiempo
    """

    y = np.zeros( (int(nt)+1, len(y0)))

    y[0] = y0


    # Arrays auxiliares
    
    k1 = np.zeros( len(y0) )
    
    k2 = np.zeros( len(y0) )

    k3 = np.zeros( len(y0) )

    k4 = np.zeros( len(y0) )


    
    
    for i in range(nt):

        k1 = h * myfunc( y[i], **kwargs )

        k2 = h * myfunc( y[i] + 0.5 * k1, **kwargs )

        k3 = h * myfunc( y[i] + 0.5 * k2, **kwargs )

        k4 = h * myfunc( y[i] + k3, **kwargs )

        y[i+1] = y[i]   +   ( k1 + 2.*k2 + 2.*k3 + k4 ) / 6.

        

    return y



