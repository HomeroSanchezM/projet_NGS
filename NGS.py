import sys #pour donner des parametres lors de l'appel de la fonction sur le terminal ou dans un .sh
import pandas as pd #pour faire des tableau

####################
#HOMERO 18/10/2024#
####################
#La fonction liste_flags prend en entrée une chaîne de caractères 
#correspondant au chemin d'un fichier SAM et retourne une liste contenant
#l'ensemble des flags du fichier

def liste_flags(ficher_sam):
    file = open(ficher_sam, 'r') #ouverture en mode lecture
    l_flags = [] #creation de liste vide pour contenir les flags
    for i_ligne in file :
        if i_ligne[0]!="@": #verifie que la ligne ne commence pas par @
            l_colonnes = i_ligne.split()  #colonne correspond a une liste des element de chaque ligne qui etait separée par des tabulation
            l_flags.append(l_colonnes[1]) #on ajoute a la liste flag la valeur de la deuxieme colonne de chaque ligne
    file.close()#on referme le fichier SAM
    return l_flags

 
#print(liste_flags(sys.argv[1]))#Prend le chemin du chemin du fichier sam donne en parametre


#La fonction Dico_flags prend en entrée une liste de flags et
#retourne un dictionnaire avec comme clés les flags et comme 
#valeurs le nombre de fois que chaque flag est présent dans la liste.

def Dico_flags(l_flags):
    d_flags = {} #creation d'un dictionnaire vide
    for i_flag in l_flags:
        if i_flag in d_flags: #verifie si une clé est déja dans le dictionnaire
            d_flags[i_flag] += 1 #si elle exite deja, on l'incremente de 1
        else:
            d_flags[i_flag] = 1 #si la clé n'existe pas on la crée et on lui donne la valeur 1
    return d_flags

print(Dico_flags(liste_flags(sys.argv[1]))) # afficher le dictionnaire 


###################
#MICKAEL 19/10/2024
###################

#On va élaborer un dictionnaire de correspondances des flags établir une correspondance pour pouvoir 

#création d'un tableau avec le dictionnaire retourner par Di	co_flags
d_flags =Dico_flags(liste_flags(sys.argv[1]))
#Créer un DataFrame à partir du dictionnaire
t_flags = pd.DataFrame(list(d_flags.items()), columns=['Flag', 'nb de fois présent'])
#Afficher le tableau
print(t_flags)

