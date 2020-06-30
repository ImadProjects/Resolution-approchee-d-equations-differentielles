# Resolution-approchee-d-equations-differentielles-Modelisation-de-systemes-dynamiques

                    Projet algorithmique numérique numéro 6
                    

Partie 1: methods.py et graphe_erreur.py

methods.py est un fichier contenant les méthodes de base permettant la résolution d'équations différentielles.
Elles sont présentes sous quatre fonctions itératives de base (step_euler, step_middle, step_heun et step_rk4), itérées par meth_n_step qui, grâce au critère de convergence fourni par step_compare, permet de satisfaire une proximité epsilon avec la solution, par appel répétés dans meth_epsilon.

L'exécution du fichier lance une série de tests avec un epsilon par défaut de 0.01, qui peut être modifié avec un argument de ligne de commande.
Il est demandé à l'utilisateur s'il souhaite, pour le cas en dimension 1, observer les graphes d'erreur associés à chaque méthode.

Les méthodes sont alors testées sur les deux équations différentielles proposées en première partie du sujet. Sont affichés les champs des tangentes, les résultats approchés et exacts, et éventuellement les erreurs.

-----

graphe_erreur.py, une fois exécuté, produit un graphe comparatif d'erreur sur le cas de test en dimension 1, pour des pas allant de 10 à 1000 en progression logarithmique. 

Partie 2: Nmaillons.py

  Cette partie traite de l'application à un système de pendule à N maillons.

Pour obtenir les figures présentes dans le rapport :
python3 Nmaillons.py

Partie 3: Lotka-Volterra.py

  Cette partie met en oeuvre et compare 3 méthodes: le modèle malthusien,
le modèle de Verhulst et le modèle de Lotka-Volterra

Pour compiler:
python3 Lotka-Volterra.py

![alt text](lotka_volterra.png)
![alt text](malthusien.png)