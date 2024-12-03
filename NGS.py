import sys 
import pandas as pd 
import re 
import argparse 

fichier_sam = sys.argv[1]
Sep = ("-" * 70) + "\n"

# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                1. CREATION OF THE MAIN DICTIONARY BY ITERATIVE STRUCTURE FOR                                                                              #
# __________________________________________________________________________________________________________________________________________________________________________________________________________#

# The dico_extraction1 function takes as input a character string corresponding to the path of a SAM file and returns a dictionary in the form of DICOEXTRACTION1.

def dico_extraction1(fichier_sam):
    file = open(fichier_sam, 'r')         # Opening in play mode
    d_sam = {}                            # Creation of an empty dictionary to contain the information from the reads
    id_ligne = 1
    for i_ligne in file:
        if i_ligne[0] != "@":             # Check that the line does not begin with @.
            l_colonnes = i_ligne.split()  # Divides the line into columns separated by tabs and extracts the columns
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

            d_sam[id_ligne] = {           # Add information to the dictionary
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
    file.close()                        # We reclose the file sam
    return d_sam 

d_sam = dico_extraction1(fichier_sam)   # To call function dico_extraction1 for get the dictionnary d_sam

# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                    2. CIGARs TREATMENT : extract the absolute values and calculate sum and relative's values                                                              #
# __________________________________________________________________________________________________________________________________________________________________________________________________________#
def analyse_CIGAR(d_sam):
    comptes_CIGAR = {
        'M': ["Alignes", 0], 'I': ["Insertions", 0], 'D': ["Deletions", 0],
        'N': ["Sauts de bases", 0], 'S': ["Soft Clipping", 0], 'H': ["Hard Clipping", 0],
        'P': ["Completion", 0], '=': ["Match exact", 0], 'X': ["Mismatch", 0]}

    REGEX_CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
    TOTAL_OPE_CIG = 0

    # Calculate the sums for each motif CIGAR
    for read in d_sam.values():
        matches = REGEX_CIGAR.findall(read["CIGAR"])
        for SUM_OPE_CIG_str, OPE_CIG in matches:
            SUM_OPE_CIG = int(SUM_OPE_CIG_str)
            if OPE_CIG in comptes_CIGAR:
                comptes_CIGAR[OPE_CIG][1] += SUM_OPE_CIG
                TOTAL_OPE_CIG += SUM_OPE_CIG

    # Calculate the percentage and update the dictionnary
    for OPE_CIG, (commentaire, count) in comptes_CIGAR.items():
        pourcentage = (count / TOTAL_OPE_CIG * 100) if TOTAL_OPE_CIG > 0 else 0
        comptes_CIGAR[OPE_CIG].append(pourcentage)  # add percentage at the dictionnary

    return comptes_CIGAR
    
# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                        3. Relative distribution of the nucleotides in the sequences                                         
# __________________________________________________________________________________________________________________________________________________________________________________________________________#

def analyse_SEQ(d_sam):
    comptes_base = {'A': [0, "Adenine"], 'T': [0, "Thymine"], 'G': [0, "Guanine"], 'C': [0, "Cytosine"]}

    for read in d_sam.values():             # For each sequence
        Seq_d_sam = read["SEQ"]
        for base in Seq_d_sam:              # For each base of each sequence
            if base in comptes_base:
                comptes_base[base][0] += 1  # Update count of the find nucleotide

    total_BASE = sum(count[0] for count in comptes_base.values())  # Total of nucleotides sum

    PRCT_BASES = {  # Relative Distribution of Bases in Sequences
        base: (count[0] / total_BASE * 100) if total_BASE > 0 else 0
        for base, count in comptes_base.items()
    }

    return comptes_base, PRCT_BASES, total_BASE

# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                          4. Treatments FLAGs : translate and relative distribution                                                                   #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#

def decodage_flags(flag_value):
    # Dictionary storing bits as keys and comments as values
    d_Binary_sam = {
        1: "A",    # Paired read
        2: "B",    # Segment paired correctly according to alignment criteria
        4: "C",    # Segment not aligned
        8: "D",    # Mate is not aligned
        16: "E",   # Read is mapped to the reverse strand
        32: "F",   # Mate is mapped to the reverse strand
        64: "G",   # First read in a pair (5'->3')
        128: "H",  # Second read in a pair on the reverse strand
        256: "I",  # Secondary alignment (non-specific, multiple alignment)
        512: "J",  # Read failed quality control filters
        1024: "K", # PCR or optical duplicate
        2048: "L"  # Supplemental alignment (non-specific, multiple alignment)
    }

    l_synthese = []
    # Each bit represents a specific flag and its associated comment:
    for i_bit, s_commentaire in d_Binary_sam.items():
        if flag_value & i_bit:  # Check whether each bit is activated or not for the flag value
            l_synthese.append(s_commentaire)

    return l_synthese

# Analyze the flags using the information extracted by dico_extraction1 contained in d_sam

def analyse_flag(d_sam):
    d_flags = {}
    for read in d_sam.values():
        flag_d_sam = read["FLAG"]

        if flag_d_sam in d_flags: #Check if a key already exists in the dictionary
            d_flags[flag_d_sam] += 1 # If it already exists, increment its value by 1
        else:
            d_flags[flag_d_sam] = 1  # If the key does not exist, create it and assign a value of 1
    return d_flags


# Create a dictionary returned by analyse_flag
d_flags = analyse_flag(d_sam)

# Add the comments from decodage_flags to this dictionary

for i_flag in d_flags:
    l_decodage = decodage_flags(int(i_flag))  # Create a list with the comments corresponding to each flag
    i_nb_fois_present = d_flags[i_flag] # Retrieve the value representing the number of flags present, associated with each flag
    d_flags[i_flag] = [i_nb_fois_present,l_decodage] # Associate each flag with a list containing the number of times the flag is present and a list with the corresponding comments for that flag

# Convert the dictionnary into a list of tuples [(Key, value1, value2), ...]
data = [(f"{cle}  |", f"{valeurs[0]} |", f"{valeurs[1]}  |") for cle, valeurs in d_flags.items()]
t_flags = pd.DataFrame(data, columns=["Flag    ", "Occurences       ", "Decodage"])

# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                                     5 CHROMOSOMAL POSITION PROCESSING:                                                                                   #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#
def analyse_Dpos(d_pos):
    d_posD = {}
    for read in d_pos.values():
        POS = read["POS"]
        if POS in d_posD:
            d_posD[POS] += 1
        else:
            d_posD[POS] = 1

    return d_posD
# __________________________________________________________________________________________________________________________________________________________________________________________________________# #                                                                                       6 ALIGNEMENT ANALISYS                                                                                              #	#__________________________________________________________________________________________________________________________________________________________________________________________________________#

# Provide the percentage of correctly paired reads
read_aligne = 0
read_non_aligne = 0
# Pair read
read_aligne_paire_non = 0

for i_flag in d_flags:
    if "B" in d_flags[i_flag][1]:  
        read_aligne += d_flags[i_flag][0]
    if ("C" in d_flags[i_flag][1]): 
        read_non_aligne += d_flags[i_flag][0]
    if ("B" in d_flags[i_flag][1]) and ("D" in d_flags[i_flag][1]):  # Number of reads aligned with the unaligned pair
        read_aligne_paire_non += d_flags[i_flag][0]

# __________________________________________________________________________________________________________________________________________________________________________________________________________#
#                                                                                7.  MAPPING QUALITY ANALYSIS                                                                                    #
# __________________________________________________________________________________________________________________________________________________________________________________________________________# 
def analyse_qualite(d_sam):
    d_qual = {}
    for read in d_sam.values():
        qual_d_sam = read["MAPQ"]

        if qual_d_sam in d_qual: # Check if a key already exists in the dictionary
            d_qual[qual_d_sam] += 1  # If it already exists, increment its value by 1
        else:
            d_qual[qual_d_sam] = 1   # If the key does not exist, create it and assign a value of 1 
    return d_qual

# __________________________________________________________________________________________________________________________________________________________________________________________________ #
#                                                    8. DEFINITIONS OF OPTIONS TO BE PASSED AS PARAMETERS TO THE MAIN SHELL SCRIPT > $1 AND SCHEDULING                                                #
# __________________________________________________________________________________________________________________________________________________________________________________________________ #

parser = argparse.ArgumentParser(description="Analyse du fichier SAM.")
parser.add_argument("sam_file", help="Chemin vers le fichier SAM.")
parser.add_argument("--all", action="store_true", help="Executer toutes les analyses (par defaut).")
parser.add_argument("--cigar", action="store_true", help="Analyser des motifs CIGAR.")
parser.add_argument("--base", action="store_true", help="Analyse de la distribution nucleotidiques du fichier.")
parser.add_argument("--flag", action="store_true", help="Traduction des flags ")
parser.add_argument("--pos", action="store_true", help="Distribution des reads en fonction de leurs positions")
parser.add_argument("--ali", action="store_true", help="Analyse de l'appariement")
parser.add_argument("--qual", action="store_true", help="Analyse de la qualite de l'alignement")
args = parser.parse_args()

# Define the conditional structure to perform all analyses if no arguments are provided. Otherwise, execute the requested analyses:

if args.all or not any([args.cigar, args.base, args.flag, args.pos, args.qual, args.ali]):
    print("Rapport integrale de l'analyse du fichier SAM : \n")
   # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    # preliminary analysis
    print(Sep, "PRELIMINARY ANALYSIS  \n", Sep, "\n")

    print("number of segments : ", len(d_sam), "\n")
    print("single read:", "\n")

    print("\u2022 segment aligned", "\n")
    print("segments properly aligned according to the aligner : ", read_aligne)
    print("percentage of segments properly aligned : ", format((read_aligne / len(d_sam)) * 100, '.3f'), " %", "\n")

    print("\u2022 segments unmapped", "\n")
    print("number of segment unmapped : ", read_non_aligne)
    print("percentage of segments unmapped : ", format((read_non_aligne / len(d_sam)) * 100, '.3f'), " %",
          "\n")
    #print("pair read", "\n")

    #print("\u2022 les paires de reads où un seul read de la paire est entierement mappe et l’autre non mappe", "\n")
    #print("nombre de read aligne avec la paire non aligne", read_aligne_paire_non)
    #print("pourcentage de read correctement apparie : ", format((read_aligne_paire_non / len(d_sam)) * 100, '.3f'),
    #      " %",
    #      "\n")

    
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    # Call the analyse_SEQ function
    comptes, pourcentages, total = analyse_SEQ(d_sam)
    
    # Display the results of the bases
    data_bases = [(f"{cle} |", f"{Val_SEQ[1]} |", f"{Val_SEQ[0]} |", f"{pourcentages[cle]:.2f} %") for cle, Val_SEQ in comptes.items()]
    t_data_bases = pd.DataFrame(data_bases, columns=["Motif       ", "Nom          ", "Occurences    ", "Valeur relative"])

    print(" ", Sep, "NUCLEOTIDES ANALYSIS : Relative distribution of the nucleotides in the segments  \n", Sep, t_data_bases, "\n")
    print("percentage of GC of all the segments : ",f"{pourcentages['C']+pourcentages['G']:.2f} %" )
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    # Analysis of the CIGAR
    comptes_CIGAR = analyse_CIGAR(d_sam)

    # Construct the DataFrame for display in the terminal.
    Data_cigar = [(f"{cle_CIG}|", f"{Val_CIG[0]} |", f"{Val_CIG[1]} |", f"{Val_CIG[2]:.3f} % |") for cle_CIG, Val_CIG in
              comptes_CIGAR.items()]
              
    t_Data_cigar = pd.DataFrame(Data_cigar,
                            columns=["Motif       ", "Nom            ", "Occurences          ", "Valeur relative"])
                            
    print(" ", Sep, "1. ANALYSE DES CIGARS : COMPTAGES DES MOTIFS ET DISTRIBUTION RELATIVE \n", Sep, "\n", t_Data_cigar, "\n")
	
   
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    # Analysis of the flags
    data = [(f"{cle}  |", f"{valeurs[0]} |", f"{valeurs[1]}  |") for cle, valeurs in d_flags.items()]
    t_flags = pd.DataFrame(data, columns=["Flag    ", "Occurences       ", "Decodage"])

    print(Sep, "ANALYSE DES FLAGS : OCCURENCES ET TRADUCTION \n", Sep, "different commentaires possibles: \n",
          "\u2022 A: Read apparie.\n",
          "\u2022 B: Segments apparies correctement selon les critères de l'aligneur. \n",
          "\u2022 C: Segment particulier non aligne. \n",
          "\u2022 D: Segment complementaire non aligne sur le brin negatif. \n",
          "\u2022 E: Segment est reverse complement.\n",
          "\u2022 F: Segment complementaire est reverse complement.\n",
          "\u2022 G: Il s'agit du premier read d'une paire sur le brin positif (5'->3').\n",
          "\u2022 H: Il s'agit du second read d'une paire sur le brin negatif (5' -> 3').\n",
          "\u2022 I: Alignement secondaire (non specifique, alignement multiple).\n",
          "\u2022 J: Read qui n'a pas passe les filtres de qualite.\n",
          "\u2022 K: Duplication due à la PCR ou au processus optique.\n",
          "\u2022 L: Alignement supplementaire (non specifique, alignement multiple).\n", "\n", t_flags)
    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    #starts positions
    print(f"nombre de positions : {len(analyse_Dpos(d_sam))}")

    d_posD = analyse_Dpos(d_sam)

    Data_pos = [(f"{cle}   |", f"{val}   |") for cle, val in d_posD.items()]
    t_Data_pos = pd.DataFrame(Data_pos, columns=["Position de depart", "Nombre de reads"])

    print(" ", Sep, "DISTRIBUTIONS DES READS PAR POSITIONS DE DEPART  \n", Sep, t_Data_pos, "\n")
    

   

    # __________________________________________________________________________________________________________________________________________________________________________________________________ #
    #Analysis of the quality
    d_qual = analyse_qualite(d_sam)

    data = [(f"{cle}  |", f"{valeurs} |") for cle, valeurs in d_qual.items()]
    t_qual = pd.DataFrame(data, columns=["Qualite    ", "Occurences       "])
    print(Sep, "ANALYSE DE LA QUALITe DE MAPPING\n", Sep, "\n",t_qual)

elif args.cigar :
    # Analyze the CIGAR string :
    comptes_CIGAR = analyse_CIGAR(d_sam)

    # Construct the DataFrame for display in the terminal
    Data_cigar = [(f"{cle_CIG}|", f"{Val_CIG[0]} |", f"{Val_CIG[1]} |", f"{Val_CIG[2]:.3f} % |") for cle_CIG, Val_CIG in
              comptes_CIGAR.items()]
    t_Data_cigar = pd.DataFrame(Data_cigar,
                            columns=["Motif       ", "Nom            ", "Occurences          ", "Valeur relative"])
                            
    print(" ", Sep, "ANALYSE DES CIGARS : COMPTAGES DES MOTIFS ET DISTRIBUTION RELATIVE \n", Sep, "\n", t_Data_cigar,
          "\n" )
elif args.base :
    # Call the analyse_SEQ function
    comptes, pourcentages, total = analyse_SEQ(d_sam)
    
    # Display the results of the bases
    data_bases = [(f"{cle} |", f"{Val_SEQ[1]} |", f"{Val_SEQ[0]} |", f"{pourcentages[cle]:.2f} %") for cle, Val_SEQ in comptes.items()]
    t_data_bases = pd.DataFrame(data_bases, columns=["Motif       ", "Nom          ", "Occurences    ", "Valeur relative"])

    print(" ", Sep, "NUCLEOTIDES ANALYSIS : Relative distribution of the nucleotides in the segments  \n", Sep, t_data_bases, "\n")
    print("percentage of GC of all the segments : ",f"{pourcentages['C']+pourcentages['G']:.2f} %" )


elif args.flag :
    data = [(f"{cle}  |", f"{valeurs[0]} |", f"{valeurs[1]}  |") for cle, valeurs in d_flags.items()]
    t_flags = pd.DataFrame(data, columns=["Flag    ", "Occurences       ", "Decodage"])

    print(Sep, "ANALYSE DES FLAGS : OCCURENCES ET TRADUCTION \n", Sep, "different commentaires possibles: \n",
          "\u2022 A: Read apparie.\n",
          "\u2022 B: Segments apparies correctement selon les critères de l'aligneur. \n",
          "\u2022 C: Segment particulier non aligne. \n",
          "\u2022 D: Segment complementaire non aligne sur le brin negatif. \n",
          "\u2022 E: Segment est reverse complement.\n",
          "\u2022 F: Segment complementaire est reverse complement.\n",
          "\u2022 G: Il s'agit du premier read d'une paire sur le brin positif (5'->3').\n",
          "\u2022 H: Il s'agit du second read d'une paire sur le brin negatif (5' -> 3').\n",
          "\u2022 I: Alignement secondaire (non specifique, alignement multiple).\n",
          "\u2022 J: Read qui n'a pas passe les filtres de qualite.\n",
          "\u2022 K: Duplication due à la PCR ou au processus optique.\n",
          "\u2022 L: Alignement supplementaire (non specifique, alignement multiple).\n", "\n", t_flags)
elif args.pos :
    print(f"nombre de positions : {len(analyse_Dpos(d_sam))}")

    d_posD = analyse_Dpos(d_sam)

    Data_pos = [(f"{cle}   |", f"{val}   |") for cle, val in d_posD.items()]
    t_Data_pos = pd.DataFrame(Data_pos, columns=["Position de depart", "Nombre de reads"])

    print(" ", Sep, "DISTRIBUTIONS DES READS PAR POSITIONS DE DEPART  \n", Sep, t_Data_pos, "\n")
elif args.ali :
    print(Sep, "PRELIMINARY ANALYSIS  \n", Sep, "\n")

    print("number of segments : ", len(d_sam), "\n")
    print("single read:", "\n")

    print("\u2022 segment aligned", "\n")
    print("segments properly aligned according to the aligner : ", read_aligne)
    print("percentage of segments properly aligned : ", format((read_aligne / len(d_sam)) * 100, '.3f'), " %", "\n")

    print("\u2022 segment unmapped", "\n")
    print("number of segment unmapped : ", read_non_aligne)
    print("percentagz of segments unmapped : ", format((read_non_aligne / len(d_sam)) * 100, '.3f'), " %",
          "\n")
    
    #print("pair read", "\n")
    
    #print("\u2022 les paires de reads où un seul read de la paire est entierement mappe et l’autre non mappe", "\n")
    #print("nombre de read aligne avec la paire non aligne", read_aligne_paire_non)
    #print("pourcentage de read correctement apparie : ", format((read_aligne_paire_non / len(d_sam)) * 100, '.3f'),
    #      " %",
    #      "\n")

elif args.qual :
    d_qual = analyse_qualite(d_sam)

    data = [(f"{cle}  |", f"{valeurs} |") for cle, valeurs in d_qual.items()]
    t_qual = pd.DataFrame(data, columns=["Qualite    ", "Occurences       "])

    print(Sep, "ANALYSE DE LA QUALITe DE MAPPING\n", Sep, "\n", t_qual)
