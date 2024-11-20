import sys  # Pour donner des parametres lors de l'appel de la fonction sur le terminal ou dans un .sh
import pandas as pd  # Pour faire des tableau
import re  # Exploitation des regex pour extraire les motifs du CIGAR en autre.
import matplotlib.pyplot as plt  # Pour faire des représentation graphiques des Outputs
import argparse  # Pour gérer le système d'options en argument dans bash.


# DANS LE README :
# DIRE CE QUE LE PROGRAMME FAIT
# FAIRE UN REQUIREMENT TXT AVEC TOUS LES PREREQUIS ET SURTOUT LENVIRONNEMENT DETRAVAIL
# QUELLE TYPE DE FICHIER EN ENTREE ET OU IL DOIT SE TROUVER ET SOUS QUEL FORMAT
# COMMENT LE SCRIPT FONCTIONNE.
# BIEN EXPLIQUE L'OUTPUT.
# LICENCE POUR LES GENS QUI UTILISE NOTRE CODE ET LA VERSION,

# DICOEXTRACTION1
# {
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
# }

fichier_sam = sys.argv[1]
Sep = ("-" * 70) + "\n"


# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                1. CREATION DU DICTIONNAIRE PRINCIPAL PAR STRUCTURE ITERATIVE FOR                                                                         #
# __________________________________________________________________________________________________________________________________________________________________________________________________________#

# La fonction dico_extraction1 prend en entrée une chaîne de caractères  correspondant au chemin d'un fichier SAM et retourne un dictionnaire  sous la forme du DICOEXTRACTION1.

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


# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                    2. TRAITEMENT DES CIGARS : EXTRACTION DES VALEURS ABSOLUES, CALCUL DES SOMMES ET VALEURS RELATIVES (REGEX)                                            #
# __________________________________________________________________________________________________________________________________________________________________________________________________________#
def analyse_CIGAR(d_sam):
    comptes_CIGAR = {
        'M': ["Alignés", 0], 'I': ["Insertions", 0], 'D': ["Délétions", 0],
        'N': ["Sauts de bases", 0], 'S': ["Soft Clipping", 0], 'H': ["Hard Clipping", 0],
        'P': ["Complétion", 0], '=': ["Match exact", 0], 'X': ["Mismatch", 0]}

    REGEX_CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
    TOTAL_OPE_CIG = 0

    # Calcul des sommes pour chaque motif CIGAR
    for read in d_sam.values():
        matches = REGEX_CIGAR.findall(read["CIGAR"])
        for SUM_OPE_CIG_str, OPE_CIG in matches:
            SUM_OPE_CIG = int(SUM_OPE_CIG_str)
            if OPE_CIG in comptes_CIGAR:
                comptes_CIGAR[OPE_CIG][1] += SUM_OPE_CIG
                TOTAL_OPE_CIG += SUM_OPE_CIG

    # Calcul des pourcentages et mise à jour du dictionnaire
    for OPE_CIG, (commentaire, count) in comptes_CIGAR.items():
        pourcentage = (count / TOTAL_OPE_CIG * 100) if TOTAL_OPE_CIG > 0 else 0
        comptes_CIGAR[OPE_CIG].append(pourcentage)  # Ajouter le pourcentage au dictionnaire

    return comptes_CIGAR


# Analyser les CIGAR
comptes_CIGAR = analyse_CIGAR(d_sam)

# Construction du Dataframe pour l'affichage dans le terminal.
Data_cigar = [(f"{cle_CIG}|", f"{Val_CIG[0]} |", f"{Val_CIG[1]} |", f"{Val_CIG[2]:.3f} % |") for cle_CIG, Val_CIG in
              comptes_CIGAR.items()]
t_Data_cigar = pd.DataFrame(Data_cigar,
                            columns=["Motif       ", "Nom            ", "Occurences          ", "Valeur relative"])


# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                  3. ANALYSE NUCLEOTIDIQUE : COMPTAGES ET DISTRIBUTIONS RELATIVES                                                                         #
# __________________________________________________________________________________________________________________________________________________________________________________________________________#

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
data_bases = [(f"{cle} |", f"{Val_SEQ[1]} |", f"{Val_SEQ[0]} |", f"{pourcentages[cle]:.2f} %") for cle, Val_SEQ in
              comptes.items()]
t_data_bases = pd.DataFrame(data_bases, columns=["Motif       ", "Nom          ", "Occurences    ", "Valeur relative"])


# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                           4 TRAITEMENT DES FLAGS : TRADUCTION ET DISTRIBUTIONS                                                                           #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#
# Fonction pour décoder la valeur d'un flag en affichant les commentaires correspondants

def decodage_flags(valeur_du_flag):
    # Dictionnaire qui stocke des bits comme clés et des commentaires comme valeurs
    d_Binary_sam = {
        1: "A",  # Read apparié.
        2: "B",  # Segment apparié correctement selon les critères de l'aligneur.
        4: "C",  # Segment particulier non aligné
        8: "D",  # Segment complémentaire aligné sur le brin négatif.
        16: "E",  # Segment aligné sur le brin négatif.
        32: "F",  # L'autre read est aligné en réverse sur le brin positif.
        64: "G",  # Il s'agit du premier read d'une paire sur le brin positif (5'->3').
        128: "H",  # Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').
        256: "I",  # Alignement secondaire (non spécifique, alignement multiple).\n"
        512: "J",  # Read qui n'a pas passé les filtres de qualité.
        1024: "K",  # Duplication due à la PCR ou au processus optique.
        2048: "L"  # Alignement supplémentaire (non spécifique, alignement multiple).
    }
    l_synthese = []

    # Chaque bit représente un flag spécifique et son commentaire associé :
    for i_bit, s_commentaire in d_Binary_sam.items():
        if valeur_du_flag & i_bit:  # En gros ici on vérifie l'activation ou non de chaque bit pour la valeur du flag
            l_synthese.append(s_commentaire)

    return l_synthese

# Analyse des flag en utilisant les info extraite par dico_extraction1 contenu par d_sam

def analyse_flag(d_sam):
    d_flags = {}
    for read in d_sam.values():
        flag_d_sam = read["FLAG"]

        if flag_d_sam in d_flags:  # verifie si une clé est déja dans le dictionnaire
            d_flags[flag_d_sam] += 1  # si elle exite deja, on l'incremente de 1
        else:
            d_flags[flag_d_sam] = 1  # si la clé n'existe pas on la crée et on lui donne la valeur 1
    return d_flags


# création d'un  dictionnaire retourner par analyse_flag
d_flags = analyse_flag(d_sam)

# ajout a ce dictionnaire les commmentaire de decodate_flags

for i_flag in d_flags:
    l_decodage = decodage_flags(int(i_flag))  # création d'une liste avec les commentaires correspondant a chaque flag
    i_nb_fois_present = d_flags[i_flag]  # recupere la valeur du nombre de flag present associe a chaque flag
    d_flags[i_flag] = [i_nb_fois_present,
                       l_decodage]  # associe a chaque flag une liste avec le nombre de fois que le flag est present et une liste avec les commentaires correspondant au flag

# Conversion du dictionnaire en une liste de tuples [(clé, valeur1, valeur2), ...]

data = [(f"{cle}  |", f"{valeurs[0]} |", f"{valeurs[1]}  |") for cle, valeurs in d_flags.items()]
t_flags = pd.DataFrame(data, columns=["Flag    ", "Occurences       ", "Decodage"])


# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                             5 TRAITEMENT DES POSITIONS CHROMOSOMIQUES :                                                                             #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#
def analyse_Dpos(d_pos):
    d_posD = {}
    for read in d_pos.values():
        POS = read["POS"]
        if POS in d_posD:
            d_posD[POS] += 1
        else:
            d_posD[POS] = 1

    return d_posD


print(f"nombre de positions : {len(analyse_Dpos(d_sam))}")

d_posD = analyse_Dpos(d_sam)

Data_pos = [(f"{cle}   |", f"{val}   |") for cle, val in d_posD.items()]
t_Data_pos = pd.DataFrame(Data_pos, columns=["Position de départ", "Nombre de reads"])

# plt.bar(d_posD.keys(), d_posD.values(), color='g')
# plt.xlabel('Positions de départ')
# plt.ylabel('Nombre de reads')
# plt.title('Distribution des reads par position de départ')
# plt.show()


# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                                       6  A DEFINIR PAR HOMERO                                                                                            #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#

# on donne le pourcentage de read correctement apparié
read_aligné = 0
read_non_aligné = 0
read_aligné_paire_non = 0
for i_flag in d_flags:
    if "B" in d_flags[i_flag][1]:  # si read aligné
        read_aligné += d_flags[i_flag][0]
    if ("C" in d_flags[i_flag][1]):  # si read non aligné
        read_non_aligné += d_flags[i_flag][0]
    if ("B" in d_flags[i_flag][1]) and ("D" in d_flags[i_flag][1]):  # nombre de read aligné avec la paire non aligné
        read_aligné_paire_non += d_flags[i_flag][0]

    # __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                          7. TRAITEMENT DE LA QUALITE DE MAPPING                                                                                           #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#


def analyse_qualité(d_sam):
    d_qual = {}
    for read in d_sam.values():
        qual_d_sam = read["MAPQ"]

        if qual_d_sam in d_qual:  # verifie si une clé est déja dans le dictionnaire
            d_qual[qual_d_sam] += 1  # si elle exite deja, on l'incremente de 1
        else:
            d_qual[qual_d_sam] = 1  # si la clé n'existe pas on la crée et on lui donne la valeur 1
    return d_qual


d_qual = analyse_qualité(d_sam)

data = [(f"{cle}  |", f"{valeurs} |") for cle, valeurs in d_qual.items()]
t_qual = pd.DataFrame(data, columns=["Qualité    ", "Occurences       "])

# ____________________________________________________________________________________________________________________________________________________________________________________________________ #
#                                                    DEFINITIONS DES OPTIONS A PASSER EN PARAMETRES DU MAIN SH en argument > $1  ET ORDONNANCEMENT                                                   #
# __________________________________________________________________________________________________________________________________________________________________________________________________ #
parser = argparse.ArgumentParser(description="Analyse du fichier SAM.")
parser.add_argument("sam_file", help="Chemin vers le fichier SAM.")
parser.add_argument("--all", action="store_true", help="Exécuter toutes les analyses (par défaut).")
parser.add_argument("--cigar", action="store_true", help="Analyser des motifs CIGAR.")
parser.add_argument("--base", action="store_true", help="Analyse de la distribution nucléotidiques du fichier.")
parser.add_argument("--flag", action="store_true", help="Traduction des flags ")
parser.add_argument("--pos", action="store_true", help="Distribution des reads en fonction de leurs positions")
parser.add_argument("--ali", action="store_true", help="Analyse de l'appariement")
parser.add_argument("--qual", action="store_true", help="Analyse de la qualité de l'alignement")
args = parser.parse_args()


# On définit la structure conditionnelle pour réaliser l'intégralité des analyses si pas d'arguments renseignés. Sinon, on exécute les analyses demandées :

if args.all or not any([args.cigar, args.base, args.flag, args.pos, args.qual]):
    print("Rapport intégrale de l'analyse du fichier SAM : \n")
    print(" ", Sep, "1. CREATION DU DICTIONNAIRE PRINCIPAL : EXTRACTIONS DES ELEMENTS TABULES DU FICHIER SAM ")
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(" ", Sep, "2. ANALYSE DES CIGARS : COMPTAGES DES MOTIFS ET DISTRIBUTION RELATIVE \n", Sep, "\n", t_Data_cigar,
          "\n", )
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(" ", Sep, "3. ANALYSE NUCLEOTIDIQUE : COMPTAGES ET DISTRIBUTION RELATIVE  \n", Sep, t_data_bases, "\n")
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(Sep, "4. ANALYSE DES FLAGS : OCCURENCES ET TRADUCTION \n", Sep, "différent commentaires possibles: \n",
          "\u2022 A: Read apparié.\n",
          "\u2022 B: Segments appariés correctement selon les critères de l'aligneur. \n",
          "\u2022C: Segment particulier non aligné. \n",
          "\u2022 D: Segment complémentaire non aligné sur le brin négatif. \n",
          "\u2022 E: Segment aligné sur le brin négatif.\n",
          "\u2022 F: L'autre read est aligné en réverse sur le brin négatif.\n",
          "\u2022 G: Il s'agit du premier read d'une paire sur le brin positif (5'->3').\n",
          "\u2022 H: Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').\n",
          "\u2022 I: Alignement secondaire (non spécifique, alignement multiple).\n",
          "\u2022 J: Read qui n'a pas passé les filtres de qualité.\n",
          "\u2022 K: Duplication due à la PCR ou au processus optique.\n",
          "\u2022 L: Alignement supplémentaire (non spécifique, alignement multiple).\n", "\n", t_flags)
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(" ", Sep, "5. DISTRIBUTIONS DES READS PAR POSITIONS DE DEPART  \n", Sep, t_Data_pos, "\n")
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(Sep, "6. ANALYSE DE L'ALIGNEMENT  \n", Sep, "\n")

    print("nombre de reads : ", len(d_sam), "\n")
    print("single read:", "\n")

    print("\u2022 read mappés", "\n")
    print("nombre de reads apparié correctement selon les critères de l'aligneur : ", read_aligné)
    print("pourcentage de read correctement apparié : ", format((read_aligné / len(d_sam)) * 100, '.3f'), " %", "\n")

    print("\u2022 read non mappés", "\n")
    print("nombre de reads non aligné : ", read_non_aligné)
    print("pourcentage de read correctement apparié : ", format((read_non_aligné / len(d_sam)) * 100, '.3f'), " %",
          "\n")

    print("\u2022 les paires de reads où un seul read de la paire est entierement mappé et l’autre non mappé", "\n")
    print("nombre de read aligné avec la paire non aligné", read_aligné_paire_non)
    print("pourcentage de read correctement apparié : ", format((read_aligné_paire_non / len(d_sam)) * 100, '.3f'),
          " %",
          "\n")

    print("pair read", "\n")
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    print(Sep, "7. ANALYSE DE LA QUALITÉ DE MAPPING\n", Sep, "\n",t_qual)

elif args.cigar :
    print(" ", Sep, "2. ANALYSE DES CIGARS : COMPTAGES DES MOTIFS ET DISTRIBUTION RELATIVE \n", Sep, "\n", t_Data_cigar,
          "\n", )
elif args.base :
    print(" ", Sep, "3. ANALYSE NUCLEOTIDIQUE : COMPTAGES ET DISTRIBUTION RELATIVE  \n", Sep, t_data_bases, "\n")
elif args.flag :
    print(Sep, "4. ANALYSE DES FLAGS : OCCURENCES ET TRADUCTION \n", Sep, "différent commentaires possibles: \n",
          "\u2022 A: Read apparié.\n",
          "\u2022 B: Segments appariés correctement selon les critères de l'aligneur. \n",
          "\u2022C: Segment particulier non aligné. \n",
          "\u2022 D: Segment complémentaire non aligné sur le brin négatif. \n",
          "\u2022 E: Segment aligné sur le brin négatif.\n",
          "\u2022 F: L'autre read est aligné en réverse sur le brin négatif.\n",
          "\u2022 G: Il s'agit du premier read d'une paire sur le brin positif (5'->3').\n",
          "\u2022 H: Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').\n",
          "\u2022 I: Alignement secondaire (non spécifique, alignement multiple).\n",
          "\u2022 J: Read qui n'a pas passé les filtres de qualité.\n",
          "\u2022 K: Duplication due à la PCR ou au processus optique.\n",
          "\u2022 L: Alignement supplémentaire (non spécifique, alignement multiple).\n", "\n", t_flags)
elif args.pos :
    print(" ", Sep, "5. DISTRIBUTIONS DES READS PAR POSITIONS DE DEPART  \n", Sep, t_Data_pos, "\n")
elif args.ali :
    print(Sep, "6. ANALYSE DE L'ALIGNEMENT  \n", Sep, "\n")

    print("nombre de reads : ", len(d_sam), "\n")
    print("single read:", "\n")

    print("\u2022 read mappés", "\n")
    print("nombre de reads apparié correctement selon les critères de l'aligneur : ", read_aligné)
    print("pourcentage de read correctement apparié : ", format((read_aligné / len(d_sam)) * 100, '.3f'), " %", "\n")

    print("\u2022 read non mappés", "\n")
    print("nombre de reads non aligné : ", read_non_aligné)
    print("pourcentage de read correctement apparié : ", format((read_non_aligné / len(d_sam)) * 100, '.3f'), " %",
          "\n")

    print("\u2022 les paires de reads où un seul read de la paire est entierement mappé et l’autre non mappé", "\n")
    print("nombre de read aligné avec la paire non aligné", read_aligné_paire_non)
    print("pourcentage de read correctement apparié : ", format((read_aligné_paire_non / len(d_sam)) * 100, '.3f'),
          " %",
          "\n")

    print("pair read", "\n")
elif args.qual :
    print(Sep, "7. ANALYSE DE LA QUALITÉ DE MAPPING\n", Sep, "\n", t_qual)





