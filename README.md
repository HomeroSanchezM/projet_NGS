# projet_NGS

Dans le fichier main.sh on élabore une stratégie pour faire quelques contrôles préliminaires sur la validité du fichier passé en paramètre au programme :
- On construit des expressions booléennes pour les vérifications courantes du type fichier vide, format de fichier, etc ...
- Avec une instruction conditionnelle on test les opérandes booléennes et renvois l'éxécutions des différents scripts pythons si tout est OK. Le cas échéant on retourne un message avec les erreurs .

Dans le fichier NGS.py on a crée une fonction liste_flags qui prend en entrée une chaîne de caractères correspondant au chemin d'un fichier SAM et retourne une liste contenant l'ensemble des flags du fichier

Dans le fichier NGS.py on a crée une fonction Dico_flags prend en entrée une liste de flags et retourne un dictionnaire avec comme clés les flags et comme valeurs le nombre de fois que chaque flag est présent dans la liste.

Dans le fichier NGS.py en utilisant la biblioteque Panda, on crée un tableau avec les flags et le nombre de fois que chaque flag est présent dans la liste. 