# -*- coding: cp1252 -*-

import numpy as np
import matplotlib.pyplot as plt
import methods as mt


def modelisation_1_maillon(g, l):
    """Modélisation d'un pendiule  1 maillon
    retourne une fonction qui traduit l'equation d'un pendule a 1 maillon"""
    f = lambda t,Y: np.array([Y[1],- (g/l) *np.sin(Y[0])])
    return f
    
def theta_1_maillon(theta0, g, l, t0, tf, eps, meth):
    """renvoie le tableau de couples (theta, thetap)"""
    f = modelisation_1_maillon(g, l)
    y0 = np.array([theta0,0.])
    y = mt.meth_epsilon(y0,t0,tf,eps,f,meth)
    return y
    
def frequence_1_maillon(theta0, g, l, t0, tf, eps, meth):
    """renvoie la fréquence d'oscillation du systeme dun pendule � 1 maillon"""
    y = theta_1_maillon(theta0, g, l, t0, tf, eps, meth) #tableau des valeurs de theta et thetaprime
    i = 0
    n = float(len(y))
    pas = (tf - t0)/n
    
    deb_var = abs(y[i+1][0] - y[i][0])/(y[i+1][0] - y[i][0]) #variation initiale
    var = deb_var
    while ( var == deb_var and ( i < len(y) - 1)): #tant que la variation est dans le meme sens
        i = i + 1
        a = y[i+1][0] - y[i][0]
        var = abs(a)/a #mise a jour de la variation
    T = i * pas
    if T == 0:
        return 0
    return 1/float(2*T) #retourne la frequence
        

def affiche_frequence(g, l, t0, tf, eps, meth):
    """affiche la fr�quence d'oscillation d'un pendule � 1 maillon en fonction de theta initial """
    THETA = np.arange(-np.pi/2.,np.pi/2.,1e-2) #tableau des valeurs de theta
    FREQUENCE = []
    FREQUENCE = [ frequence_1_maillon(THETA[k], g, l, t0, tf, eps, meth) for k in range(len(THETA))] #tableau des fr�quences en fonction de theta initial
    plt.plot(THETA, FREQUENCE)
    F = [ (1./(2*np.pi))*np.sqrt(g/l) for k in THETA] #frequence au voisinage de theta
    plt.plot(THETA, F)
    plt.xlabel('Frequence')
    plt.ylabel('Theta initial')
    plt.savefig('rapport/methode.png', transparent = True)
    plt.show()
    

def modelisation_2_maillons(g, l1, l2, m1, m2):
    """Mod�lisation d'un pendule � deux maillon
    retourne une fonction qui traduit l'equation d'un pendule a deux maillons"""
    f = lambda t, Y: np.array([Y[1],
                                   (-g*(2*m1+m2)*np.sin(Y[0])-m2*g*np.sin(Y[0]-2*Y[2])-2*np.sin(Y[0]-Y[2])*m2*(l2*Y[3]**2+l1*np.cos(Y[0]-Y[2])*Y[1]**2))/(l1*(2*m1+m2-m2*np.cos(2*Y[0]-2*Y[2]))),
                                   Y[3],
                                   (2*np.sin(Y[0]-Y[2])*((m1+m2)*l1*Y[1]**2+g*(m1+m2)*np.cos(Y[0])+l2*m2*np.cos(Y[0]-Y[2])*Y[3]**2))/(l2*(2*m1+m2-m2*np.cos(2*Y[0]-2*Y[2])))])
    return f
    
def theta_2_maillons(theta01, theta02, g, l1, l2, m1, m2, t0, tf, eps, meth):
    """renvoie le tableau des quadruplets (theta1, theta1p, theta2, theta2p)"""
    f = modelisation_2_maillons(g, l1, l2, m1, m2)
    y0 = np.array([theta01,0.,theta02,0.])
    y = mt.meth_epsilon(y0,t0,tf,eps,f,meth)
    return y
    
def affiche_traj_pos_double(theta01, theta02, g, l1, l2, m1, m2, t0, tf, eps, meth):
    """affiche la position du pendule au bout en fonction du temps"""
    plt.clf()
    theta = theta_2_maillons(theta01, theta02, g, l1, l2, m1, m2, t0, tf, eps, meth)
    n = len(theta)
    
    theta1 = [theta[i][0] for i in range(n)] #tableau de theta1 en fonction du temps
    theta2 = [theta[i][2] for i in range(n)] #tableau de theta2 en fonction du temps
    
    x = [l1*np.sin(theta1[i])-l2*np.sin(theta2[i]) for i in range(n)] #projection sur x
    y = [-l1*np.cos(theta1[i])-l2*np.cos(theta2[i]) for i in range(n)] #projection sur y

    plt.plot(x, y, linewidth=1.0)
    plt.xlabel('Axe des x')
    plt.ylabel('Axs des y')
    plt.show()
    
def temps_retournement(tab):
    """Retourne le temps du premier retournement a partir d'une courbe donnee"""
    i = 1
    l = len(tab)
    while ( ( i < l - 1 ) and np.pi - abs(tab[i][0]) > 0. and np.pi - abs(tab[i][2]) > 0. ): #condition de retournement
        i = i + 1
    if ( i == l - 1 ):
        return 0
    return i
    
def affiche_temps_retournement(nmax, g, l1, l2, m1, m2, t0, tmax, meth):
    """Graphe du temps mis par le pendule double pour se retourner"""
    t = np.arange(0.0,tmax+tmax/nmax,tmax/nmax) #discr�tisation du temps
    T = np.zeros([nmax,nmax])
    theta = (2*np.pi)/T.shape[0] #initialisation des angles
    theta0 = -np.pi #initialisation des angles
    p,q = T.shape[0], T.shape[1]
    for i in range(p):
        for j in range(q):
            y0 = np.array([i*theta+theta0, 0.0, j*theta+theta0, 0.0])
            f = modelisation_2_maillons(g, l1, l2, m1, m2)
            res = mt.meth_n_step(y0,t0,nmax,10.0/nmax,f,meth) #tableau de theta1, theta1p, theta2, theta2p en fonction du temps
            tmps = temps_retournement(res) #temps de retournement
            T[i][j] = t[tmps] #valeur la la matrice
    
    plt.clf()      
    plt.imshow(T)

    plt.gca().invert_yaxis()
    plt.xlabel('Axe des theta1')
    plt.ylabel('Axe des theta2')
    plt.gcf()
    plt.clim()
    cb = plt.colorbar()
    cb.set_label('Temps en seconde')
    plt.savefig('rapport/temps.png', transparent = True)
    plt.show()
    
    


if __name__=="__main__":
    affiche_frequence(9.81, 1., 0, 20, 1e-1, mt.step_rk4)
    affiche_traj_pos_double(50, 70, 9.81, 2, 3, 2, 4, 0, 20, 1e-1, mt.step_rk4)
    affiche_traj_pos_double(51, 70, 9.81, 2, 3, 2, 4, 0, 20, 1e-1, mt.step_rk4)    
    affiche_traj_pos_double(51, 71, 9.81, 2, 3, 2, 4, 0, 20, 1e-1, mt.step_rk4)
    affiche_temps_retournement(100, 9.81, 2, 3, 2, 4, 0.0, 10.0, mt.step_rk4)