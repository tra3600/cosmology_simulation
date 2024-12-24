# cosmology_simulation
programme C++ qui simule l expansion de l'univers en prenant en compte les plus récentes avancées en matière de modèle cosmologique

Pour simuler l'expansion de l'univers en C++, nous devons prendre en compte les modèles cosmologiques actuels, tels que le modèle Lambda-CDM (Lambda Cold Dark Matter). Ce modèle est le modèle standard de la cosmologie, qui inclut la constante cosmologique (Λ) associée à l'énergie sombre et la matière noire froide.

Concepts Clés
Paramètres Cosmologiques : La constante de Hubble (H0), la densité de matière (Ωm), la densité d'énergie sombre (ΩΛ), etc.
Équation de Friedmann : Décrit l'évolution de l'univers en fonction de ces paramètres.
Redshift (z) : Utilisé pour mesurer la distance et l'expansion de l'univers.
Étapes
Définir les paramètres cosmologiques.
Utiliser l'équation de Friedmann pour calculer l'évolution de l'univers.
Simuler l'expansion de l'univers sur une période de temps donnée.

Explications
Cosmology Class:

Constructeur : Initialise les paramètres cosmologiques H0, Ωm, et ΩΛ.
scaleFactor : Calcule le facteur d'échelle à un moment donné.
hubbleParameter : Calcule le paramètre de Hubble à un moment donné.
simulateExpansion : Simule l'expansion de l'univers sur une période donnée.
Main Function:

Initialise les paramètres du modèle Lambda-CDM.
Simule l'expansion de l'univers de 0 à 13,8 milliards d'années avec un pas de temps de 0,1 milliard d'années.
Affiche les résultats.
Compilation et Exécution
Assurez-vous que tous les fichiers sont dans le même répertoire, puis compilez et exécutez le programme avec :

g++ -o cosmology_simulation main.cpp Cosmology.cpp
./cosmology_simulation

Ce programme est une simplification du modèle cosmologique Lambda-CDM et ne prend pas en compte tous les facteurs complexes de l'expansion de l'univers. Pour une simulation plus précise, des modèles plus avancés et des données réelles seraient nécessaires.
