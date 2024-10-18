import sys

#La fonction liste_flags prend en entrée une chaîne de caractères 
#correspondant au chemin d'un fichier SAM et retourne une liste contenant
#l'ensemble des flags du fichier

def liste_flags(ficher_sam):
    file = open(ficher_sam, 'r') #ouverture en mode lecture
    flags = [] #creation de liste vide pour contenir les flags
    for ligne in file :
        if ligne[0]!="@*": #verifie que la ligne ne commence pas par @
            colonnes = ligne.split()  #colonne correspont a une liste des element de chaque ligne qui etait separée par des tabulation
            flags.append(colonnes[1]) #on ajoute a la liste flag la valeur de la deuxieme colonne de chaque ligne
    file.close()
    return flags

 
print(liste_flags(sys.argv[1]))#Prend le chemin du chemin du fichier sam donne en parametre