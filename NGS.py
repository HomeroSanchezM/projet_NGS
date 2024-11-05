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

            # Ajouter les informations dans le dictionnaire
            d_sam[QNAME] = {
                "FLAG": FLAG,
                "RNAME": RNAME,
                "POS": POS,
                "MAPQ": MAPQ,
                "CIGAR": CIGAR,
                "RNEXT": RNEXT,
                "PNEXT": PNEXT,
                "TLEN": TLEN,
            }
            
    file.close()#on referme le fichier SAM
    return d_sam


#print(dico_extraction1(sys.argv[1])) #decommenter pour tester

#################
#Mickael 06/11/24
#################

def analyse_CIGAR(d_sam):
    # Initialiser un dictionnaire pour compter les opérations CIGAR
    comptes_CIGAR = {op: 0 for op in 'MIDNSHP=X*'}

    # Parcourir chaque lecture dans le dictionnaire d_sam
    for read in d_sam.values():
        CIGAR_d_sam = read["CIGAR"]  # Extraire le CIGAR de la lecture

        iCig = 0  # Initialiser l'indice de la chaîne CIGAR
        while iCig < len(CIGAR_d_sam):
            # Trouver la taille de l'opération
            jCig = iCig
            while jCig < len(CIGAR_d_sam) and CIGAR_d_sam[jCig].isdigit():
                jCig += 1
            
            taille_str = CIGAR_d_sam[iCig:jCig]  # Extraire la taille

            if not taille_str:  # Si la taille est vide
                iCig = jCig + 1  # Passer à l'élément suivant pour éviter une boucle infinie
                continue  # Passer à la prochaine itération du while principal

            # Vérifier si la taille est un nombre entier valide
            if taille_str.isdigit():
                taille = int(taille_str)  # Convertir la taille en entier
            else:
                iCig = jCig + 1  # Passer à l'élément suivant si ce n'est pas un nombre valide
                continue

            operation = CIGAR_d_sam[jCig]  # Extraire l'opération

            comptes_CIGAR[operation] += taille  # Incrémenter le compte de l'opération

            # Passer à l'élément suivant du CIGAR
            iCig = jCig + 1

    # Calculer le total des opérations et les pourcentages en une seule étape
    total_operations = sum(comptes_CIGAR.values())
    pourcentages_CIGAR = {op: (compte / total_operations * 100) if total_operations > 0 else 0
                          for op, compte in comptes_CIGAR.items()}

    # Afficher les résultats
    for operation, pourcentage in pourcentages_CIGAR.items():
        print(f"{operation}: {pourcentage:.3f}%")
    print(f"\nTotal des opérations : {total_operations}")

    #return comptes_CIGAR, pourcentages_CIGAR, total_operations

# Exemple d'appel
print(analyse_CIGAR(dico_extraction1(fichier_sam)))


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


###################
# MICKAEL 19/10/2024
###################

# Dictionnaire qui stocke des bits comme clés et des commentaires comme valeurs
d_Binary_sam = {
    1: "Read apparié.",
    2: "Segment apparié correctement selon les critères de l'aligneur.",
    4: "Segment particulier non aligné.",
    8: "Segment complémentaire aligné sur le brin négatif.",
    16: "Segment aligné sur le brin négatif.",
    32: "L'autre read est aligné en réverse sur le brin positif.",
    64: "Il s'agit du premier read d'une paire sur le brin positif (5'->3').",
    128: "Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').",
    256: "Alignement secondaire (non spécifique, alignement multiple).",
    512: "Read qui n'a pas passé les filtres de qualité.",
    1024: "Duplication due à la PCR ou au processus optique.",
    2048: "Alignement supplémentaire (non spécifique, alignement multiple).",
}

# Fonction pour décoder la valeur d'un flag en affichant les commentaires correspondants
def decodage_flags(valeur_du_flag):
    l_synthese = []

    # Chaque bit représente un flag spécifique et son commentaire associé
    for i_bit, s_commentaire in d_Binary_sam.items():
        if valeur_du_flag & i_bit : # En gros ici on vérifie l'activation ou non de chaque bit pour la valeur du flag
            l_synthese.append(f"- {s_commentaire}")

    return "\n".join(l_synthese)
    

#i_test = 99
#print(decodage_flags(i_test))

####################
#HOMERO 28/10/2024#
####################

#création d'un  dictionnaire retourner par Dico_flags
d_flags =Dico_flags(liste_flags(sys.argv[1]))

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

