import sys          # Pour donner des parametres lors de l'appel de la fonction sur le terminal ou dans un .sh
import pandas as pd # Pour faire des tableau
import re           # Exploitation des regex pour extraire les motifs du CIGAR



# DANS LE README :
# DIRE CE QUE LE PROGRAMME FAIT
# FAIRE UN REQUIREMENT TXT AVEC TOUS LES PREREQUIS ET SURTOUT LENVIRONNEMENT DETRAVAIL
# QUELLE TYPE DE FICHIER EN ENTREE ET OU IL DOIT SE TROUVER ET SOUS QUEL FORMAT
# COMMENT LE SCRIPT FONCTIONNE.
# BIEN EXPLIQUE L'OUTPUT.
# LICENCE POUR LES GENS QUI UTILISE NOTRE CODE ET LA VERSION, 

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
import re
import sys

fichier_sam = sys.argv[1]
Sep = ("-" * 70)+ "\n"
#__________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                1. CREATION DU DICTIONNAIRE PRINCIPAL PAR STRUCTURE ITERATIVE FOR                                                                         #
#__________________________________________________________________________________________________________________________________________________________________________________________________________#

# La fonction dico_extraction1 prend en entrée une chaîne de caractères 
# correspondant au chemin d'un fichier SAM et retourne un dictionnaire 
# sous la forme du DICOEXTRACTION1.

def dico_extraction1(fichier_sam):
    file = open(fichier_sam, 'r')  # Ouverture en mode lecture
    d_sam = {}  # Création du dictionnaire vide pour contenir les infos des reads
    id_ligne = 1
    for i_ligne in file:
        if i_ligne[0] != "@":  # Vérifie que la ligne ne commence pas par @
            l_colonnes = i_ligne.split()  # Découpe la ligne en colonnes séparées par des tabulations
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
                "QNAME": QNAME,
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
            id_ligne += 1
    file.close()  # On referme le fichier SAM
    return d_sam

# Appeler la fonction dico_extraction1 pour obtenir le dictionnaire d_sam
d_sam = dico_extraction1(fichier_sam)

#__________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                      2. TRAITEMENT DES CIGARS : EXTRACTION DES VALEURS ABSOLUES, CALCUL DES SOMMES ET VALEURS RELATIVES (REGEX)                                                          #
#__________________________________________________________________________________________________________________________________________________________________________________________________________#

# Analyse des CIGAR
def analyse_CIGAR(d_sam):
    comptes_CIGAR = {
        'M': ["Alignés", 0],        'I': ["Insertions", 0],    'D': ["Délétions", 0],
        'N': ["Sauts de bases", 0], 'S': ["Soft Clipping", 0], 'H': ["Hard Clipping", 0],  
        'P': ["Complétion", 0],     '=': ["Match exact", 0],   'X': ["Mismatch", 0] }
        
    REGEX_CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
    TOTAL_OPE_CIG = 0

    for read in d_sam.values():
        matches = REGEX_CIGAR.findall(read["CIGAR"])
        for SUM_OPE_CIG_str, OPE_CIG in matches:
            SUM_OPE_CIG = int(SUM_OPE_CIG_str)
            if OPE_CIG in comptes_CIGAR:
                comptes_CIGAR[OPE_CIG][1] += SUM_OPE_CIG
                TOTAL_OPE_CIG += SUM_OPE_CIG

    # Calcul des pourcentages
    for OPE_CIG, (commentaire, count) in comptes_CIGAR.items():
        pourcentage = (count / TOTAL_OPE_CIG * 100) if TOTAL_OPE_CIG > 0 else 0
        print(f"{commentaire}: {count} ({pourcentage:.2f}%)")

    return comptes_CIGAR

analyse_CIGAR(d_sam)


#__________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                  2. ANALYSE NUCLEOTIDIQUE : COMPTAGES ET DISTRIBUTIONS RELATIVES                                                                         #
#__________________________________________________________________________________________________________________________________________________________________________________________________________#

def analyse_SEQ(d_sam):  
    comptes_base = {'A': [0, "Adénine"], 'T': [0, "Thymine"], 'G': [0, "Guanine"], 'C': [0, "Cytosine"]}
    
    for read in d_sam.values():  # Pour chaque séquence
        Seq_d_sam = read["SEQ"]
        for base in Seq_d_sam:  # Pour chaque base de chaque séquence
            if base in comptes_base:
                comptes_base[base][0] += 1  # Mise à jour du comptage de la base trouvée

    total_BASE = sum(count[0] for count in comptes_base.values())  # Total des bases comptées
    
    PRCT_BASES = {  # Distribution relative des bases dans les séquences
        base: (count[0] / total_BASE * 100) if total_BASE > 0 else 0
        for base, count in comptes_base.items()
    }
    
    return comptes_base, PRCT_BASES, total_BASE


# Appel de la fonction analyse_SEQ
comptes, pourcentages, total = analyse_SEQ(d_sam)

# Affichage des résultats des bases
data_bases = [(cle, valeur[1], valeur[0], f"{pourcentages[cle]:.2f} %") for cle, valeur in comptes.items()]


t_data_bases = pd.DataFrame(data_bases, columns=["Motif   ", " Nom ", "Valeur absolue    ", "   Valeur relative"])

print(" ",Sep,"3. ANALYSE NUCLEOTIDIQUE : COMPTAGES ET DISTRIBUTIONS RELATIVES  \n",Sep,t_data_bases)
print(Sep)

#__________________________________________________________________________________________________________________________________________________________________________________________________________# #                                      4 TRAITEMENT DES FLAGS : TRADUCTION ET DISTRIBUTIONS                                                                                                                #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#




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
#HOMERO 11/11/2024
###################
#Analyse des flag en utilisant les info extraite par dico_extraction1 contenu par d_sam

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
t_flags = pd.DataFrame(data, columns=['Flag', 'Occurences', 'Decodage'])

# Affichage du DataFrame
print(t_flags)


    # Écrire les pourcentages dans un fichier CSV
    #with open("Output_DATA.csv", mode="w", newline="") as file:
     #   write_CIGAR = csv.writer(file)
      #  writer.writerow(["Operation", "Pourcentage"])
       # for operation, pourcentage in pourcentages_CIGAR.items():
        #    writer.writerow([operation, pourcentage])

    #print("Les pourcentages CIGAR ont été exportés vers pourcentages_CIGAR.csv")


