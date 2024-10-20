# projet_NGS : ELABORATION D'UN PROGRAMME D'EVALUATION DE FICHIER SAM.
########################################################################################################################################################################################################
#                                                                         COQUERELLE MICKAEL - HOMERO SANCHEZ 18.10.2024                                                                                #      
#			        		                                     --  SCORE : SAM Characterization and Observational Report for Evaluations  --                                                              
########################################################################################################################################################################################################

# 1) Fichier main.sh 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Dans le fichier main.sh on élabore une stratégie pour faire quelques contrôles préliminaires sur la validité du fichier passé en paramètre au programme :
- On construit des expressions booléennes pour les vérifications courantes du type fichier vide, format de fichier, etc ...
- Avec une instruction conditionnelle on test les opérandes booléennes et renvois l'éxécutions des différents scripts pythons si tout est OK. Le cas échéant on retourne un message avec les erreurs .

# 2) Fichier NGS.py
- Dans le fichier NGS.py on a crée une fonction liste_flags qui prend en entrée une chaîne de caractères correspondant au chemin d'un fichier SAM et retourne une liste contenant l'ensemble des flags du fichier
- On a crée une fonction Dico_flags prend en entrée une liste de flags et retourne un dictionnaire avec comme clés les flags et comme valeurs le nombre de fois que chaque flag est présent dans la liste.
-En utilisant la biblioteque Panda, on crée un tableau avec les flags et le nombre de fois que chaque flag est présent dans la liste.

# 3) Fichier Translate_FLAG.py
- L'utilité de ce script python est de décoder les flags SAM (Sequence Alignment/Map) en interprétant la valeur binaire pour chaque flag du fichier, on évalue bit à bit l'expression binaire à partir de la valeur décimale données dans le fichier. On propose ensuite en sortie une liste de commentaire associé à chaque flag et donc à chaque Reads. L'intérêt sous jacent sera d'exploiter cette fonction pour l'alloué à chaque Flag F données pour N occurence du fichier dans un tableau de sortie vers une colonne dédiée.

- 
