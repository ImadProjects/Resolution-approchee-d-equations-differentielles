import numpy as np
import matplotlib.pyplot as plt
import sys
from locale import atof

#Problème de Cauchy : (y0,f)
#Ordre des paramètres de f : t,y
#Les y seront des scalaires, ou en dimension supérieure, des arrays (et pas des listes)

#Calculs d'une valeur yp d'une fonction à partir de la valeur y précédente en fonction f de celle-ci, du point d'évaluation t, et du pas h par différentes méthodes

def step_euler(y,t,h,f):
	yp=y+h*f(t,y)
	return yp
	
def step_middle(y,t,h,f):
	y_half=y+(h/2.)*f(t,y)
	p=f(t+h/2.,y_half)
	yp=y+h*p
	return yp
	
def step_heun(y,t,h,f):
	p1=f(t,y)
	y2=y+h*p1
	p2=f(t+h,y2)
	yp=y+h*(p1+p2)/2.
	return yp
	
def step_rk4(y,t,h,f):
	p1=f(t,y)
	y2=y+(h/2.)*p1
	p2=f(t+h/2.,y2)
	y3=y+(h/2.)*p2
	p3=f(t+h/2.,y3)
	y4=y+h*p3
	p4=f(t+h,y4)
	yp=y+(h/6.)*(p1+2*p2+2*p3+p4)
	return yp
	
#Construit un tableau de n valeurs par appels successifs d'une méthode à un pas

def meth_n_step(y0,t0,N,h,f,meth):
	if (type(y0)) is np.ndarray:
		values=np.zeros((N,y0.shape[0]))
	else:
		values=np.zeros(N)
	t=t0
	values[0]=y0
	for i in range(1,N):
		t+=h
		values[i]=meth(values[i-1],t,h,f)
	return values
	
#Compare le maximum des différences entre les points communs de deux arrays y1 et y2 à un epsilon donné.

def step_compare(y1,y2,eps):
	y=np.zeros(len(y1))
	for i in range(len(y1)):
		if type(y1[0]) is np.ndarray:
			y[i]=np.linalg.norm(y1[i]-y2[2*i])
		else:
			y[i]=abs(y1[i]-y2[2*i])
	return max(y)>eps
	
#Calcule avec un N doublé à chaque itération N valeurs y sur l'intervalle |t0,tf] jusqu'à ce que le maximum des différences de deux itérations soit inférieur ou égal à epsilon.

def meth_epsilon(y0,t0,tf,eps,f,meth):
	N=10
	h=(tf-t0)/N
	y1=meth_n_step(y0,t0,N,h,f,meth)
	N*=2
	h=(tf-t0)/N
	y2=meth_n_step(y0,t0,N,h,f,meth)
	while(step_compare(y1,y2,eps)):
		N*=2
		h/=2.
		y1=y2
		y2=meth_n_step(y0,t0,N,h,f,meth)
	return y2
	
#Récupération des arguments (epsilon et choix de l'affichage d'erreur)
	
def opts():
	if len(sys.argv)==1:
		eps=1e-2
		print("Epsilon par défaut : ",eps)
		print("Pour une valeur autre : python 3", sys.argv[0]," eps avec epsilon le critère de convergence souhaité")
	else:
		eps=atof(sys.argv[1])
		print("Epsilon choisi : ",eps)
	error=input("Afficher le graphe d'erreur des méthodes ? Y/N ")
	if error=="Y" or error=="y":
		error=True
	elif error=="N" or error=="n":
		error=False
	else:
		print("Réponse invalide")
		exit()
	return eps,error
	
#Tracé des valeurs en fonction du temps (réelles et approchées), éventuellement des graphes d'erreur
	
def compare_plot(values,sol,start,end,title,error=False):
	X=np.linspace(start,end,len(values))
	c_values=[sol(X[i]) for i in range(len(X))]
	plt.clf()
	plt.plot(X,c_values,label="Solution réelle")
	plt.plot(X,values,label="Solution approchée")
	plt.legend()
	plt.xlabel("Temps")
	plt.ylabel("Valeurs des coordonnées")
	plt.title(title)
	plt.show()
	if error==True:
		plt.clf()
		errors=[abs(c_values[i]-values[i]) for i in range(len(X))]
		plt.plot(X,errors,label="Erreur sur l'intervalle")
		plt.legend()
		plt.xlabel("Temps")
		plt.ylabel("Erreur")
		plt.title("Graphe d'erreur")
		plt.show()
		
#Tracé d'un portrait de phase en deux dimensions

def phase_graph(values,sol,start,end,title):
	X=np.linspace(0,end,len(values))
	c_values=[sol(X[i]) for i in range(len(X))]
	values_x=[values[i][0] for i in range(len(values))]
	values_y=[values[i][1] for i in range(len(values))]
	c_values_x=[c_values[i][0] for i in range(len(values))]
	c_values_y=[c_values[i][1] for i in range(len(values))]
	plt.clf()
	plt.plot(values_x,values_y,label="Portrait de phase approché")
	plt.plot(c_values_x,c_values_y,label="Portrait de phase réel")
	plt.xlabel("Coordonnée X")
	plt.ylabel("Coordonnée Y")
	plt.legend()
	plt.title(title)
	plt.show()
	
#Tracé du champ des tangentes en dimension (dim) 1 ou 2, avec une abscisse X (coïncidant en dimension 1 avec les temps t), entre a et b sur les deux axes, et selon la fonction f. title est le titre du graphe produit.

def champ_tangentes(dim,a,b,K,f,title):
	N=K*1j

	Y, X = np.mgrid[a:b:N, a:b:N]
	if(dim==1):
		U=X
		V=f(X,Y)
	else:
		U,V=f(0,[X,Y])

	normes = np.sqrt(U**2 + V**2)
	UU = U/normes
	VV = V/normes
	plt.axis([a, b, a, b])

	plt.quiver(X, Y, UU, VV,headlength=5, width=0.002)
	plt.title(title)
	plt.show()
	
if __name__=="__main__":
	(eps,error)=opts()
	print("\nCas de test en dimension 1 : (1,y'=y/(1+t**2))")
	f=lambda t,y: y/(1+t**2) #Prudence : écrire des f compatibles avec des arrays
	title="Champ des tangentes avec f(t,y)=y/(1+t^2)"
	champ_tangentes(1,1e-15,5,20,f,title) #Pour ne pas donner zéro à la division...
	sol=lambda t: np.exp(np.arctan(t))
	title="Résultat avec méthode d'Euler"
	print(title)
	values=meth_epsilon(1,0,5,eps,f,step_euler)
	compare_plot(values,sol,0,5,title,error)
	title="Résultat avec méthode du point milieu"
	print(title)
	values=meth_epsilon(1,0,5,eps,f,step_middle)
	compare_plot(values,sol,0,5,title,error)
	title="Résultat avec méthode de Heun"
	print(title)
	values=meth_epsilon(1,0,5,eps,f,step_heun)
	compare_plot(values,sol,0,5,title,error)
	title="Résultat avec méthode de Runge-Kutta d'ordre 4"
	print(title)
	values=meth_epsilon(1,0,5,eps,f,step_rk4)
	compare_plot(values,sol,0,5,title,error)

	print("\nCas de test en dimension 2 : ((1,0),y'=(-y2,y1))")
	f=lambda t,y: np.array([-y[1],y[0]])
	title="Champ des tangentes avec f(t,Y)=(-Y[1],Y[0])"
	print(title)
	champ_tangentes(2,-2*np.pi,2*np.pi,20,f,title)
	sol=lambda t: [np.cos(t),np.sin(t)]
	title="Résultat avec méthode d'Euler"
	print(title)
	values=meth_epsilon(np.array([1,0]),0,2*np.pi,eps,f,step_euler)
	compare_plot(values,sol,0,2*np.pi,title,error)
	phase_graph(values,sol,0,2*np.pi,title)
	title="Résultat avec méthode du point milieu"
	print(title)
	values=meth_epsilon(np.array([1,0]),0,2*np.pi,eps,f,step_middle)
	compare_plot(values,sol,0,2*np.pi,title,error)
	phase_graph(values,sol,0,2*np.pi,title)
	title="Résultat avec méthode de Heun"
	print(title)
	values=meth_epsilon(np.array([1,0]),0,2*np.pi,eps,f,step_heun)
	compare_plot(values,sol,0,2*np.pi,title,error)
	phase_graph(values,sol,0,2*np.pi,title)
	title="Résultat avec méthode de Runge-Kutta d'ordre 4"
	print(title)
	values=meth_epsilon(np.array([1,0]),0,2*np.pi,eps,f,step_rk4)
	compare_plot(values,sol,0,2*np.pi,title,error)
	phase_graph(values,sol,0,2*np.pi,title)
