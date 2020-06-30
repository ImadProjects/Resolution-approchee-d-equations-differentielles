from methods import *
import numpy as np
import matplotlib.pyplot as plt
import math

def values(x0, xf, N):
    """
        Renvoie un tableau de valeurs uniformément réparties entre x0 et xf
    """
    h = ( xf - x0 ) / ( N-1 )
    res = np.zeros(N)
    for i in range(N):
        res[i] = x0 + h * i
    return res

def modele_malthusien(t0,tf,eps,nb_individus,gamma,meth):
    """
         Resolution avec le modele Malthusien
    """  
    return meth_epsilon(nb_individus, t0, tf, eps,lambda y, t: gamma*t, step_rk4);

def modele_de_verhulst(t0,tf,k,eps,nb_individus,gamma,meth):
    """
         Resolution avec le modele de Verhulst
    """
    return meth_epsilon(nb_individus, t0, tf, eps,lambda y, t: gamma * t * (1 - t/k),step_rk4);

def malthusien_draw():
    """
         Affiche la courbe de la variation de la population selon le modèle Malthusien
    """  
    t0 = 0.
    tf = 3.
    eps = 1e-2
    nb_individus = 1000.

    pop_growth = modele_malthusien(t0,tf,eps,nb_individus,2,meth_epsilon)
    plt.plot(values(t0, tf, len(pop_growth)), pop_growth,"grey", label="naissances > morts")

    cst_pop = modele_malthusien(t0,tf,eps,nb_individus,0,meth_epsilon)
    plt.plot(values(t0, tf, len(cst_pop)), cst_pop,"blue", label="naissances = morts")

    pop_decrease = modele_malthusien(t0,tf,eps,nb_individus,-2,meth_epsilon)
    plt.plot(values(t0, tf, len(pop_decrease)), pop_decrease,"black", label="naissance < morts")

    plt.xlabel("Temps")
    plt.ylabel("Nombre d'individus")
    plt.axis([t0,tf,0,5000])
    plt.legend()
    plt.title("Évolution de la population avec le modele malthusien")
    plt.savefig('rapport/malthusien.png')
    plt.show()

def verhulst_draw():
    """
         Affiche la courbe de la variation de la population selon le modèle de Verhulst
    """    
    t0 = 0.
    tf = 10.
    eps = 1e-2
    nb_individus = 1000.
    gamma = 1

    k = 2000
    pop_growth = modele_de_verhulst(t0,tf,k,eps,nb_individus,gamma,meth_epsilon)
    plt.plot(values(t0, tf, len(pop_growth)), pop_growth,"grey", label="nombre d'individu < nombre limite d'individu")

    k = 600
    pop_decrease = modele_de_verhulst(t0,tf,k,eps,nb_individus,gamma,meth_epsilon)
    plt.plot(values(t0, tf, len(pop_decrease)), pop_decrease,"black", label="nombre d'individu > nombre limite d'individu")

    k = 1000
    cst_pop = modele_de_verhulst(t0,tf,k,eps,nb_individus,gamma,meth_epsilon)
    plt.plot(values(t0, tf, len(cst_pop)), cst_pop,"blue", label="nombre d'individu = nombre limite d'individu")

    plt.xlabel("Temps")
    plt.ylabel("Nombre d'individus")
    plt.axis([t0,tf,0,2000])
    plt.legend()
    plt.title("Évolution de la population avec le modèle de Verhulst")
    plt.savefig('rapport/verhulst.png')
    plt.show()

def derivee_Nt_Pt(a, b, c, d):
    return lambda y, t : np.array([t[0]*(a-b*t[1]),t[1]*(c*t[0]-d)])

def distinction_Nt_Pt(y):
    y_len = len(y)
    Nt = np.zeros(y_len)
    Pt = np.zeros(y_len)
    for i in range(y_len):
        Nt[i] = y[i][0]
        Pt[i] = y[i][1]
    return (Nt, Pt)    

def Pt_f_Nt(t0,tf,eps,a,b,c,d,nb_individus):
    modele = meth_epsilon(nb_individus, t0, tf, eps, derivee_Nt_Pt(a, b, c, d), step_rk4);
    (modele_proie, modele_pred) = distinction_Nt_Pt(modele)
    
    plt.plot(modele_proie, modele_pred)
    plt.xlabel("Nombre proies N(t)")
    plt.ylabel("Nombre predateurs P(t)")
    plt.title("Affichage de P(t) en fonction de N(t) en utilisant le modele de Lotka-Volterra." )
    plt.savefig('rapport/P_N.png')
    plt.show()


def modele_proie_predateur(t0,tf,eps,a,b,c,d,nb_individus):
    """
         Affiche les courbes modélisants un écosystème proie/prédateurs
    """  
    modele = meth_epsilon(nb_individus, t0, tf, eps, derivee_Nt_Pt(a, b, c, d), step_rk4);
    (modele_proie, modele_pred) = distinction_Nt_Pt(modele)

    plt.plot(values(t0, tf, len(modele_proie)), modele_proie,"green", label="N(t), nombre proies")
    plt.plot(values(t0, tf, len(modele_pred)), modele_pred,"purple", label="P(t), nombre predateurs")
    plt.legend()
    plt.xlabel("Temps")
    plt.ylabel("Nombre d'individus")
    plt.title("Modèle Lotka-Volterra")
    plt.savefig('rapport/lotka_volterra.png')
    plt.show()

def diagramme_solutions(nb_individus_init, var, N, a, b, c, d, t0, tf, eps):
    h = var/(2*N)
    for i in range (-N, N+1):
        for j in range (-N, N+1):
            nv_nb_individus = nb_individus_init + np.array([h*i,h*j])
            modele = meth_epsilon(nv_nb_individus, t0, tf, eps, derivee_Nt_Pt(a, b, c, d), step_rk4);
            (modele_proie, modele_predateur) = distinction_Nt_Pt(modele)
            plt.plot(modele_proie, modele_predateur,"grey")


    plt.xlabel("N(t): nombre de proies")
    plt.ylabel("P(t): nombre de predateurs")
    plt.title("Variations entre des solutions partant de conditions initiales\n très proches :%s" % nb_individus_init)
    plt.axis([35,45,5,7])
    plt.savefig('rapport/variations.png')
    plt.show()

def draw_figures():
    t0 = 0
    tf = 365
    a = 0.25
    b = 0.05
    c = 0.01 
    d = 0.1
    nb_individus = np.array([20., 15.])
    eps = 0.01
    Pt_f_Nt(t0,tf,eps,a,b,c,d,nb_individus)
    malthusien_draw()
    verhulst_draw()
    modele_proie_predateur(t0,tf,eps,a,b,c,d,nb_individus)
    diagramme_solutions(nb_individus, 6, 3,a, b, c, d, t0, tf, eps)

if __name__ == "__main__":
    draw_figures()
