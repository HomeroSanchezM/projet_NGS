import os, sys

sam = sys.argv[1]
file = open(sam, 'r')
for ligne in file :
    print(ligne, end="")
file.close()

