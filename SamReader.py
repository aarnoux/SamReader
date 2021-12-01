#!/usr/bin/env python3
#-*- coding utf-8 -*-
import os, sys, csv, re, colorama
from colorama import Fore, Back

csv.field_size_limit(sys.maxsize)   # augmentation de la taille max des champs dans le fichier csv afin de pouvoir analyser des alignements de reads > 131kb

## 1. Vérification que l'argument est un fichier SAM non vide

# on vérifie que l'argument est un fichier :
if os.path.isfile(sys.argv[1]) ==  False:
    print("File error : L'argument n'est pas un fichier.")
    exit()
# on teste si le fichier est vide :
elif os.path.getsize(sys.argv[1]) == 0:
    print("File error : Le fichier est vide.")
    exit()
# on check l'extension du fichier :
elif sys.argv[1].endswith(".sam") == False: 
    print("Format error : seul le format SAM est accepté.")
    exit()

## 2. Préparation de l'analyse du fichier

## 2.1 Ouverture et lecture du fichier

# on lit le fichier dans une variable "file", on indique que les champs sont séparés par des tabulations :
fi = open(sys.argv[1])
file = csv.reader(fi, delimiter = '\t', quoting = csv.QUOTE_NONE)
fi.close

## 2.2 Initialisations

# compilation d'expression régulière :
regex = re.compile('^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$')

# variables :
i = 0
j = 0
analyseComplete = 0
unmapped_count = 0
partially_mapped_count = 0
totally_mapped_count = 0

# dictionnaires :
dicoCigar = {'M':0,'I':0,'D':0,'S':0,'H':0,'N':0,'P':0,'X':0,'=':0}
percentCigar = {'M':0,'I':0,'D':0,'S':0,'H':0,'N':0,'P':0,'X':0,'=':0}

# tuple :
headers = ('@HD','@SQ','@RG','@PG','@CO')

## 3. Définition d'une fonction pour mettre le FLAG en binaire

def flagBinary(flag):
        flagB = bin(int(flag))  # Transform the integer into a binary.
        flagB = flagB[2:]   # Remove '0b' Example: '0b1001101' > '1001101'
        flagB = list(flagB) 
        if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
            add = 12 - len(flagB)   # We compute the difference between the maximal flag size (12) and the length of the binary flag.
            for m in range(add):
                flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
        return(flagB)

## 4. Analyse du fichier

# on prépare les fichiers de sortie :
with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file_um, open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file_pm:
# analyse ligne par ligne du SAM :
    for ligne in file:
        i = i + 1   # incrémentation du numéro de ligne à chaque tour de boucle (début à 1 et non 0, plus logique)
        k = 0
        if ligne[0].startswith(headers) == False:   # si la ligne analysée ne fait pas partie du header on analyse les champs

# la recherche d'erreurs dans les champs débute par la recherche de l'expression régulière (RE) du champ, puis :
# - si la RE est correcte on incrémente le compteur 'k' de 1
# - si la RE est incorrecte, on affiche un message d'erreur et on finit d'analyser la ligne

# on pourrait imbriquer les conditions 'if...else...' afin de gagner du temps lors de l'analyse du fichier, mais on perd de la précision lors de l'affichage des erreurs (on peut le changer en fonction du profil d'utilisateur visé)

            if (re.search('^[!-?A-~]{1,254}$', ligne[0])):
                k = k + 1
            else:
                print(Fore.RED+"QNAME field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ QNAME a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[0])+"\n"+str(ligne))
            if (re.search('^[0-9]{1,4}$', ligne[1])):
                k = k + 1
            else:
                print(Fore.RED+"FLAG field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ FLAG a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[1])+"\n"+str(ligne))
            if (regex.search(ligne[2]) or re.search('\*', ligne[2])):
                k = k + 1
            else:
                print(Fore.RED+"RNAME field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ RNAME a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[2])+"\n"+str(ligne))
            if (re.search('^[0-9]*$', ligne[3])):
                k = k + 1
            else:
                print(Fore.RED+"POS field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ POS a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[3])+"\n"+str(ligne))
            if (re.search('^[0-5][0-9]$''|^[0-9]$''|^(60)$''|^(255)$', ligne[4])):
                k = k + 1
            else:
                print(Fore.RED+"MAPQ field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ MAPQ a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[4])+"\n"+str(ligne))
            if (re.search('^\*$''|^([0-9]+[MIDNSHPX=])+$', ligne[5])):
                k = k + 1
            else:
                print(Fore.RED+"CIGAR field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ CIGAR a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[5])+"\n"+str(ligne))                             
            if (regex.search(ligne[6]) or re.search('^\*$', ligne[6]) or re.search('^=$', ligne[6])):
                k = k + 1
            else:
                print(Fore.RED+"RNEXT field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ RNEXT a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[6])+"\n"+str(ligne))
            if (re.search('^[0-9]*$', ligne[7])):
                k = k + 1
            else:
                print(Fore.RED+"PNEXT field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ PNEXT a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[7])+"\n"+str(ligne))
            if (re.search('^^-?[0-9]{1,30}$''|^-?20{31}$', ligne[8])):
                k = k + 1
            else:
                print(Fore.RED+"TLEN field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ TLEN a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[8])+"\n"+str(ligne))
            if (re.search('^\*|[A-Za-z=.]+$', ligne[9])):
                k = k + 1
            else:
                print(Fore.RED+"SEQ field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ SEQ a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[9])+"\n"+str(ligne))
            if (re.search('^[!-~]+$', ligne[10])):
                k = k + 1
            else:
                print(Fore.RED+"QUAL field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ QUAL a une expression régulière non conforme, ou un des champ de la ligne n'est pas délimité par une tabulation.\nChamp non conforme : "+str(ligne[10])+"\n"+str(ligne))

# si la ligne ne présente pas d'erreurs alors 'k' vaut 11

            if k == 11:
                flag = flagBinary(ligne[1]) # on appelle la fonction pour mettre le FLAG en binaire
                if int(flag[-3]) == 1:  # si '0100' est présent dans le binaire ça correspond à un FLAG de 4 (= read non mappé)
                    unmapped_count += 1
                    unmapped_fasta.write(str(flag))
                if int(flag[-2]) == 1:  # si '0010' est présent dans le binaire ça correspond à un FLAG de 2 (= read bien aligné)
                    if re.match('(?![0-9]+M$)', ligne[5]):  # dans le CIGAR, vérifier que la RE "XXXM" (= 100% de match) n'est pas présente
                        partially_mapped_count += 1
                        partillay_mapped_fasta.write(str(flag))
                    else:
                        totally_mapped_count +=1
                    ext = re.findall('[0-9]+\D', ligne[5])
                    for info in range(len(ext)):
                        temp = ext[info]
                        if temp[-1] == 'M' and info != 0:    # si le 'M' correspond à un mismatch, stocker la valeur dans la clé 'X' et pas 'M' (car 'M' c'est pour les matchs)
                            dicoCigar['X'] += int(temp[:(len(temp)-1)])
                        else:
                            dicoCigar[temp[-1]] += int(temp[:(len(temp)-1)])
            elif j == 0:
                j = j + 1
                k = 11 - k
                analyseComplete = input(Back.RED+"Document non analysable"+Back.RESET+" : "+str(k)+" erreur(s) d'expressions régulières trouvée(s) à la ligne "+str(i)+", souhaitez-vous rechercher les erreurs dans le reste du fichier ?\ny/n\n")
                if analyseComplete == "n":
                    exit()
                elif analyseComplete != "y":
                    print(Fore.RED+"Input error : arrêt du script"+Fore.RESET)
                    exit()            

    if analyseComplete == "y":
        print(Fore.YELLOW+"Fin de la recherche d'erreur."+Fore.RESET)
        exit()

    vtotal = 0
    for v in dicoCigar:
        vtotal += dicoCigar[v]
    for w in dicoCigar.keys():
        percentCigar[w] = round((dicoCigar[w]*100)/vtotal, 2)

    print("\ntotal reads count : "+str(i))
    print("\t-> properly aligned reads count : "+Fore.GREEN+str(totally_mapped_count+partially_mapped_count)+Fore.RESET+" ("+str(round((totally_mapped_count+partially_mapped_count)*100/i))+"% of total reads)")
    print("\t\t-> totally mapped reads count : "+str(totally_mapped_count))
    summary_file_pm.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
    print("\t\t-> partially mapped reads count : "+str(partially_mapped_count))
    summary_file_um.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
    print("\t-> unmapped reads count : "+str(unmapped_count)+"\n")
    print(dicoCigar)
    print(percentCigar)

with open("Final_Cigar_table.txt", "w") as FinalCigar:
    FinalCigar.write("Global cigar mutation observed :"+"\n"
                    +"Alignlent Match : "+str(round(percentCigar['M'],2))+"%\n"
                    +"Insertion : "+str(round(percentCigar['I'],2))+"%\n"
                    +"Deletion : "+str(round(percentCigar['D'],2))+"%\n"
                    +"Skipped region : "+str(round(percentCigar['S'],2))+"%\n"
                    +"Soft Clipping : "+str(round(percentCigar['H'],2))+"%\n"
                    +"Hard Clipping : "+str(round(percentCigar['N'],2))+"%\n"
                    +"Padding : "+str(round(percentCigar['P'],2))+"%\n"
                    +"Sequence Match : "+str(round(percentCigar['X'],2))+"%\n"
                    +"Sequence Mismatch : "+str(round(percentCigar['='],2))+"%\n")