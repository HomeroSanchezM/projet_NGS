import sys #pour donner des parametres lors de l'appel de la fonction sur le terminal ou dans un .sh
import pandas as pd #pour faire des tableau

####################
#HOMERO 5/11/2024
####################

#format du dictionnaire d'extraction:

#DICOEXTRACTION1
#{
#    "QNAME1": {
#        "FLAG": FLAG1,
#        "RNAME": RNAME1,
#        "POS": POS1,
#        "MAPQ": MAPQ1,
#        "CIGAR": CIGAR1,
#        "RNEXT": RNEXT1,
#        "PNEXT": PNEXT1,
#        "TLEN": TLEN1,
#        "SEQ": SEQ1,
#        "QUAL": QUAL1
#    },
#    "QNAME2": {
#        "FLAG": FLAG2,
#        "RNAME": RNAME2,
#        # etc.
#    },
#    # etc.
#}

fichier_sam = sys.argv[1]

#La fonction dico_extraction1 prend en entrée une chaîne de caractères 
#correspondant au chemin d'un fichier SAM et retourne une dico de la forme du DICOEXTRACTION1
def dico_extraction1(fichier_sam):
    file = open(fichier_sam, 'r') #ouverture en mode lecture
    d_sam = {} #creation du dico vide pour contenir les info des reads
    id_ligne = 1
    for i_ligne in file :
        if i_ligne[0]!="@": #verifie que la ligne ne commence pas par @
            l_colonnes = i_ligne.split()  #colonne correspond a une liste des elements de chaque ligne qui etait separée par des tabulations (découpe la ligne en colonnes)
            # Extraire les champs du fichier SAM
            QNAME = l_colonnes[0]  
            FLAG = l_colonnes[1]
            RNAME = l_colonnes[2]
            POS = l_colonnes[3]
            MAPQ = l_colonnes[4]
            CIGAR = l_colonnes[5]
            RNEXT = l_colonnes[6]
            PNEXT = l_colonnes[7]
            TLEN = l_colonnes[8]
            SEQ = l_colonnes[9]

            # Ajouter les informations dans le dictionnaire
            d_sam[id_ligne] = {
                "QNAME":QNAME,
                "FLAG": FLAG,
                "RNAME": RNAME,
                "POS": POS,
                "MAPQ": MAPQ,
                "CIGAR": CIGAR,
                "RNEXT": RNEXT,
                "PNEXT": PNEXT,
                "TLEN": TLEN,
                "SEQ": SEQ
            }
            id_ligne+=1
    file.close()#on referme le fichier SAM
    return d_sam
    

#print(dico_extraction1(sys.argv[1])) #decommenter pour tester

#################
# Mickael 07/11/24
#################
def analyse_SEQ(d_sam):  
    comptes_base = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    
    for read in d_sam.values():
        Seq_d_sam = read["SEQ"]
        
        for base in Seq_d_sam:
            if base in comptes_base:
                comptes_base[base] += 1
    
    total_BASE = sum(comptes_base.values())
    
    pourcentages_BASE = {
        'A': (comptes_base['A'] / total_BASE * 100) if total_BASE > 0 else 0,
        'T': (comptes_base['T'] / total_BASE * 100) if total_BASE > 0 else 0,
        'G': (comptes_base['G'] / total_BASE * 100) if total_BASE > 0 else 0,
        'C': (comptes_base['C'] / total_BASE * 100) if total_BASE > 0 else 0
    }
    
    return comptes_base, pourcentages_BASE, total_BASE


# Appel de la fonction dico_extraction1 pour obtenir le dictionnaire d_sam
d_sam = dico_extraction1(fichier_sam)

# Appel de la fonction analyse_SEQ avec d_sam comme argument
comptes, pourcentages, total = analyse_SEQ(d_sam)

# Affichage des résultats
print("Comptes des bases :", comptes)
print("")
print("Total des bases : ", total)

print("")  # Ligne vide pour séparer les deux parties

# Affichage des pourcentages avec 2 chiffres après la virgule
print("Pourcentages des bases :")
print("--------------------------------------------------------")

for base, pourcent in pourcentages.items():
    print(f"  {base}: {pourcent:.2f}%")
print("--------------------------------------------------------")

#################
#Mickael 06/11/24
#################

def analyse_CIGAR(d_sam):
    # Initialiser un dictionnaire pour compter les opérations CIGAR
    comptes_CIGAR = {
        'M': 0,  # Match (alignement de base)
        'I': 0,  # Insertion (délétion dans la séquence de référence)
        'D': 0,  # Deletion (insertion dans la séquence cible)
        'N': 0,  # Skipping (indication de l'existence d'un intron)
        'S': 0,  # Soft clipping (bases coupées mais conservées dans le fichier SAM)
        'H': 0,  # Hard clipping (bases coupées et non enregistrées)
        'P': 0,  # Padding (espaces réservés dans l'alignement)
        '=': 0,  # Identique (opérations qui correspondent aux bases de la séquence)
        'X': 0,  # Mismatch (mismatch entre les bases d'alignement)
        '*': 0   # Opérations non spécifiées ou indéfinies
    }

    # Parcourir chaque CIGAR dans le dictionnaire d_sam
    for read in d_sam.values():
        CIGAR_d_sam = read["CIGAR"]  # Extraire le CIGAR du dico sam

        iCig = 0  # Initialisation de notre indice pour la chaîne CIGAR
        
        while iCig < len(CIGAR_d_sam):
            # Trouver la taille de l'opération CIGAR
            jCig = iCig
            
            while jCig < len(CIGAR_d_sam) and CIGAR_d_sam[jCig].isdigit(): #Vérifier qu'on sort pas du CIGAR et que l'on traite le nombre d'occurence de l'événement
                jCig += 1
       
            taille_str = CIGAR_d_sam[iCig:jCig]  # Extraction de l'événement
            
            if not taille_str:  # Gère nos cas où l'événement est unique dans le CIGAR
                iCig = jCig + 1  # Passer à l'élément suivant pour éviter une boucle infinie
                continue  # Passer à la prochaine itération du while principal

            # Vérifier si la taille est un nombre entier valide
            if taille_str.isdigit():
                taille = int(taille_str)  # Convertir la taille en entier
            else:
                iCig = jCig + 1  # Passer à l'élément suivant si ce n'est pas un nombre valide
                continue

            operation = CIGAR_d_sam[jCig]  # Extraire l'opération

            comptes_CIGAR[operation] += taille  # Somme des comptes pour l'élément courant cigar

            # Passer à l'élément suivant du CIGAR
            iCig = jCig + 1

    # Calculer le total des opérations et les pourcentages en une seule étape
    total_operations = sum(comptes_CIGAR.values())

    # Initialisation du dictionnaire des pourcentages CIGAR
    pourcentages_CIGAR = {
        'Alignés match ou mismatch(M)': (comptes_CIGAR['M'] / total_operations * 100) if total_operations > 0 else 0,
        'Insertions (I)': (comptes_CIGAR['I'] / total_operations * 100) if total_operations > 0 else 0,
        'Délétions (D)': (comptes_CIGAR['D'] / total_operations * 100) if total_operations > 0 else 0,
        'Sauts de bases (N)': (comptes_CIGAR['N'] / total_operations * 100) if total_operations > 0 else 0,
        'Soft Clipping (S)': (comptes_CIGAR['S'] / total_operations * 100) if total_operations > 0 else 0,
        'Hard Clipping (H)': (comptes_CIGAR['H'] / total_operations * 100) if total_operations > 0 else 0,
        'Complétion (P)': (comptes_CIGAR['P'] / total_operations * 100) if total_operations > 0 else 0,
        'base match (=)': (comptes_CIGAR['='] / total_operations * 100) if total_operations > 0 else 0,
        'bases mismatch (X)': (comptes_CIGAR['X'] / total_operations * 100) if total_operations > 0 else 0,
        'Événements non spécifiés (*)': (comptes_CIGAR['*'] / total_operations * 100) if total_operations > 0 else 0
    }

    print("Distribution des données CIGAR pour évaluation de l'alignement")
    # Afficher les résultats
    for operation, pourcentage in pourcentages_CIGAR.items():
        print(f"{operation}: {pourcentage:.3f}%")

    print(f"\nTotal des opérations : {total_operations}")
    
#Affichage des informations de d'appariements et FLAG couplés aux informations de distributions des évènements stockés dans les CIGARs ; 

print(analyse_CIGAR(dico_extraction1(fichier_sam)))


#----------------------------------------------------------------------------------------#

# Pour HOMERO 07/11 : Ici il faut que l'on fasse un camembert avec maplot lib ou intégrer du code R (ggplot) que ce soit plus rigoureux et adaptés.
# Ce que je te propose sur la base de ce que j'ai fais pour les CIGARS cette nuit : 
#   - Je reproduis une structure itérative équivalente pour criblé les lectures et sortir les occurences de chaque base, evènements indel ou autre.   
#  - Je réfléchis à la manière de représenter les distributions pour CIGAR et les %ATGC dans R avec ggplot
#   - Je veux bien réfléchir à Jupyternotebook en dernier lieux

# Ouvert a tes autres suggestions si tu veux faire un travail avec un jeu de boucles sur les autres items du dictionnaires qui te parait pertinent pour représenter les données?
# Pour la partie Pourcentage CIGAR je me suis pas cassé la tete, au début j'avais fait ca avec op pour simplifier la synthaxe du dictionnaire mais c'était pas cohérent donc j'ai modifier les produits en croix et j'ai fais copier coller pour chaque évènement.


<<<<<<< HEAD
#DICOEXTRACTION2
#{
#    "FLAG": {
#        "QNAME1": FLAG1,
#        "QNAME2": FLAG2,
#        "QNAME3": FLAG3,
#        "QNAME4": FLAG4,
#        "QNAME5": FLAG5,
#        "QNAME6": FLAG6,
#         ect.
#    },
#    "RNAME": {
#        "QNAME1": RNAME1,
#        "QNAME2": RNAME2,
#         etc.
#    },
#      etc.
#}


#La fonction dico_extraction1 prend en entrée une chaîne de caractères 
#correspondant au chemin d'un fichier SAM et retourne une dico de la forme du DICOEXTRACTION2
def dico_extraction2(ficher_sam):
    file = open(ficher_sam, 'r') #ouverture en mode lecture
    # Initialiser le dictionnaire principal
    d_sam = {
        "FLAG": {},
        "RNAME": {},
        "POS": {},
        "MAPQ": {},
        "CIGAR": {},
        "RNEXT": {},
        "PNEXT": {},
        "TLEN": {},
    }

    for i_ligne in file :
        if i_ligne[0]!="@": #verifie que la ligne ne commence pas par @
            l_colonnes = i_ligne.split()  #colonne correspond a une liste des element de chaque ligne qui etait separée par des tabulation
            # Extraire les champs du fichier SAM
            QNAME = l_colonnes[0]
            FLAG = l_colonnes[1]
            RNAME = l_colonnes[2]
            POS = l_colonnes[3]
            MAPQ = l_colonnes[4]
            CIGAR = l_colonnes[5]
            RNEXT = l_colonnes[6]
            PNEXT = l_colonnes[7]
            TLEN = l_colonnes[8]

            # Ajouter les informations dans le dictionnaire principal
            d_sam["FLAG"][QNAME] = FLAG
            d_sam["RNAME"][QNAME] = RNAME
            d_sam["POS"][QNAME] = POS
            d_sam["MAPQ"][QNAME] = MAPQ
            d_sam["CIGAR"][QNAME] = CIGAR
            d_sam["RNEXT"][QNAME] = RNEXT
            d_sam["PNEXT"][QNAME] = PNEXT
            d_sam["TLEN"][QNAME] = TLEN
    file.close()
    return d_sam

#print(dico_extraction2(sys.argv[1])) #decomenter pour tester



####################
#HOMERO 18/10/2024
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

#print(Dico_flags(liste_flags(sys.argv[1]))) # afficher le dictionnaire 



###################
#MICKAEL 19/10/2024
###################

#On va élaborer un dictionnaire de correspondances des flags établir une correspondance pour pouvoir 

#création d'un tableau avec le dictionnaire retourner par Dico_flags
d_flags =Dico_flags(liste_flags(sys.argv[1]))
#Créer un DataFrame à partir du dictionnaire
t_flags = pd.DataFrame(list(d_flags.items()), columns=['Flag', 'nb de fois présent'])
#Afficher le tableau
#print(t_flags)

=======
>>>>>>> 02bd13f1aaface7058be7db16c95f9b7cadbfd11

###################
# MICKAEL 19/10/2024
###################

# Dictionnaire qui stocke des bits comme clés et des commentaires comme valeurs
d_Binary_sam = {
    1: "- Read apparié.",
    2: "- Segment apparié correctement selon les critères de l'aligneur.",
    4: "- Segment particulier non aligné.",
    8: "- Segment complémentaire aligné sur le brin négatif.",
    16: "- Segment aligné sur le brin négatif.",
    32: "- L'autre read est aligné en réverse sur le brin positif.",
    64: "- Il s'agit du premier read d'une paire sur le brin positif (5'->3').",
    128: "- Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').",
    256: "- Alignement secondaire (non spécifique, alignement multiple).",
    512: "- Read qui n'a pas passé les filtres de qualité.",
    1024: "- Duplication due à la PCR ou au processus optique.",
    2048: "- Alignement supplémentaire (non spécifique, alignement multiple).",
}

# Fonction pour décoder la valeur d'un flag en affichant les commentaires correspondants
def decodage_flags(valeur_du_flag):
    l_synthese = []

    # Chaque bit représente un flag spécifique et son commentaire associé :
    for i_bit, s_commentaire in d_Binary_sam.items():
        if valeur_du_flag & i_bit : # En gros ici on vérifie l'activation ou non de chaque bit pour la valeur du flag
            l_synthese.append(f"{s_commentaire}")
            
    return l_synthese
    

####################
#HOMERO 28/10/2024#
####################
#non utiliser
#création d'un  dictionnaire retourner par Dico_flags
#d_flags =Dico_flags(liste_flags(sys.argv[1]))
#
#ajout a ce dictionnaire les commmentaire de decodate_flags
#
#for i_flag in d_flags:
#    l_decodage = decodage_flags(int(i_flag)) #création d'une liste avec les commentaires correspondant a chaque flag 
#    i_nb_fois_present = d_flags[i_flag] #recupere la valeur du nombre de flag present associe a chaque flag 
#    d_flags[i_flag] = [i_nb_fois_present, l_decodage ] #associe a chaque flag une liste avec le nombre de fois que le flag est present et une liste avec les commentaires correspondant au flag
#
# Conversion du dictionnaire en une liste de tuples [(clé, valeur1, valeur2), ...]
#data = [(cle, valeurs[0], valeurs[1]) for cle, valeurs in d_flags.items()]
#
# Création du DataFrame avec les colonnes spécifiées
#t_flags = pd.DataFrame(data, columns=['Flag', 'nb de fois présent', 'decodage'])
#
# Affichage du DataFrame
#print(t_flags)
#
####################
#HOMERO 11/11/2024
###################
#Analyse des flag en utilisant les info extraite par dico_extraction1 contenu par d_sam (cf ligne 99)

def analyse_flag(d_sam):  
    d_flags = {}
    for read in d_sam.values():
        flag_d_sam = read["FLAG"]
        
        if flag_d_sam in d_flags: #verifie si une clé est déja dans le dictionnaire
            d_flags[flag_d_sam] += 1 #si elle exite deja, on l'incremente de 1
        else:
            d_flags[flag_d_sam] = 1 #si la clé n'existe pas on la crée et on lui donne la valeur 1
    return d_flags

#print(analyse_flag(d_sam)) # afficher le dictionnaire 

#création d'un  dictionnaire retourner par analyse_flag
d_flags =analyse_flag(d_sam)

#ajout a ce dictionnaire les commmentaire de decodate_flags

for i_flag in d_flags:
    l_decodage = decodage_flags(int(i_flag)) #création d'une liste avec les commentaires correspondant a chaque flag 
    i_nb_fois_present = d_flags[i_flag] #recupere la valeur du nombre de flag present associe a chaque flag 
    d_flags[i_flag] = [i_nb_fois_present, l_decodage ] #associe a chaque flag une liste avec le nombre de fois que le flag est present et une liste avec les commentaires correspondant au flag

# Conversion du dictionnaire en une liste de tuples [(clé, valeur1, valeur2), ...]
data = [(cle, valeurs[0], valeurs[1]) for cle, valeurs in d_flags.items()]

# Création du DataFrame avec les colonnes spécifiées
t_flags = pd.DataFrame(data, columns=['Flag', 'nb de fois présent', 'decodage'])

# Affichage du DataFrame
print(t_flags)


<<<<<<< HEAD

=======
#couverture: nb de read mappé/taille de chromosome
>>>>>>> 02bd13f1aaface7058be7db16c95f9b7cadbfd11
