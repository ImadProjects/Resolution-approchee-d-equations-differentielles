import numpy as np
import matplotlib.pyplot as plt
import methods as mt

#Comparaison à un pas donné des valeurs y calculées et des valeurs de la solution à ses points communs avec y. Retourne le maximum des écarts.

def step_error(y,a,b,sol):
	h=(b-a)/len(y)
	if type(y[0]) is np.ndarray:
		sol_values=np.zeros((len(y),y[0].shape[0]))
	else:
		sol_values=np.zeros(len(y))
	for i in range(len(y)):
		sol_values[i]=sol(a+i*h)
	if type(y[0]) is np.ndarray:
		diff_values=[np.linalg.norm(y[i]-sol_values[i]) for i in range(len(y))]
	else:
		diff_values=[abs(y[i]-sol_values[i]) for i in range(len(y))]
	diff=max(diff_values)
	return diff 
	
#Tabulation des écarts maximums pour différents pas (en croissance logarithmique), sur un intervalle [a,b], construisant le tableau de valeurs avec la méthode meth. f et y0 sont nécessaires pour l'exécution de meth_n_step.
	
def errors_over_n(sol,f,meth,y0,a,b):
	space=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
	diffs=np.zeros(len(space))
	for i in range(len(diffs)):
		h=(b-a)/space[i]
		y=mt.meth_n_step(y0,a,space[i],h,f,meth)
		diffs[i]=step_error(y,a,b,sol)
	return diffs
	
#Création de graphes d'erreurs pour les quatre méthodes sur le cas défini par sol, f, y0, et sur l'intervalle [a,b].

def plot_errors(sol,f,y0,a,b):
	diffs_euler=errors_over_n(sol,f,mt.step_euler,y0,a,b)
	diffs_middle=errors_over_n(sol,f,mt.step_middle,y0,a,b)
	diffs_heun=errors_over_n(sol,f,mt.step_heun,y0,a,b)
	diffs_rk4=errors_over_n(sol,f,mt.step_rk4,y0,a,b)
	
	plt.clf()
	X=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
	plt.plot(X,diffs_euler,label="Méthode d'Euler")
	plt.plot(X,diffs_middle,label="Méthode du point milieu")
	plt.plot(X,diffs_heun,label="Méthode de Heun")
	plt.plot(X,diffs_rk4,label="Méthode de Runge-Kutta d'ordre 4")
	plt.legend()
	
	plt.xlabel("Pas")
	plt.ylabel("Erreur relative maximale")
	plt.title("Graphe de comparaison des erreurs entre les quatre méthodes")
	plt.xscale("log")
	plt.grid(True,which="both",linestyle='-.')
	plt.show()
	

if __name__=='__main__':
	f=lambda t,y: y/(1+t**2)
	sol=lambda t: np.exp(np.arctan(t))
	y0=1
	a=0
	b=5
	plot_errors(sol,f,y0,a,b)
