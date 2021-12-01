#!/usr/bin/env python3
#-*- coding utf-8 -*-

__author__ = ("Alizée ARNOUX")
__contact__ = ("alizee.arnoux@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/14/2021"
__licence__ ="This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."
        
import os, sys, csv, re, colorama, numpy
from colorama import Fore, Back

# ********************
# ** HELP FUNCTIONS **
# ********************

def helpFLAG():
    print("Bit\t\tDescription")
    print("----------------------------------------------------------------------")
    print("1\t0x1\ttemplate having multiple segments in sequencing")
    print("2\t0x2\teach segment properly aligned according to the aligner")
    print("4\t0x4\tsegment unmapped")
    print("8\t0x8\tnext segment in the template unmapped")
    print("16\t0x10\tSEQ being reverse complemented")
    print("32\t0x20\tSEQ of the next segment in the template being reverse complemented")
    print("64\t0x40\tthe first segment in the template")
    print("128\t0x80\tthe last segment in the template")
    print("256\t0x100\tsecondary alignment")
    print("512\t0x200\tnot passing filters, such as platform/vendor quality controls")
    print("1024\t0x400\tPCR or optical duplicate")
    print("2048\t0x800\tsupplementary alignment")
    exit()

def helpRequirements():
    print("This program runs with "+Fore.RED+"python 3"+Fore.RESET+", and needs the following packages:")
    print(Fore.YELLOW+"os\tsys\tcsv\tre\tnumpy\tcolorama"+Fore.RESET)
    print("\nYou can install them using the following commands:")
    print("conda install *package name*\nor\npip install *package name*")
    exit()

def helpCigar():
    print("Op\tBAM\tDescription")
    print("----------------------------------------------------------------------")
    print("M\t0\talignment match (can be a sequence match or mismatch)")
    print("I\t1\tinsertion to the reference")
    print("D\t2\tdeletion from the reference)")
    print("N\t3\tskipped region from the reference")
    print("S\t4\tsoft clipping (clipped sequences present in SEQ)")
    print("H\t5\thard clipping (clipped sequences NOT present in SEQ)")
    print("P\t6\tpadding (silent deletion from padded reference)")
    print("=\t7\tsequence match")
    print("X\t8\tsequence mismatch")
    exit()

def helpSAM():
    print("Col\tField\tType\tRegex/Range\t\t\tDescription")
    print("----------------------------------------------------------------------")
    print("1\tQNAME\tString\t[!-?A-~]{1,254}\t\t\tQuery template NAME")
    print("2\tFLAG\tInt\t[0,2e16 − 1]\t\t\tbitwise FLAG")
    print("3\tRNAME\tString\t\*|[:rname:∧ *=][:rname:]*\tReference sequence NAME")
    print("4\tPOS\tInt\t[0,2e31 − 1]\t\t\t1-based leftmost mapping POSition")
    print("5\tMAPQ\tInt\t[0,2e8 − 1]\t\t\tMAPping Quality")
    print("6\tCIGAR\tString\t\*|([0-9]+[MIDNSHPX=])+\t\tCIGAR string")
    print("7\tRNEXT\tString\t\*|=|[:rname: *=][:rname:]*\tReference name of the mate/next read")
    print("8\tPNEXT\tInt\t[0,2e31 − 1]\t\t\tPosition of the mate/next read")
    print("9\tTLEN\tInt\t[−2e31 + 1,2e31 − 1]\t\tobserved Template LENgth")
    print("10\tSEQ\tString\t\*|[A-Za-z=.]+\t\t\tsegment SEQuence")
    print("11\tQUAL\tString\t[!-~]+\t\t\t\tASCII of Phred-scaled base QUALity+33")
    exit()

def helpProgramm():
    print("\nSamReader is a small programm that analyses SAM files. It can analyse multiple files in one go.")
    print("Choosing from one of the following options is MANDATORY :\n-t	show the results in the terminal, without saving them in a file\n-o	save the results in a file which name is given by the user, or in the default file 'summary_data_file.txt'\n-t -o	show the results int he terminal AND save them in a file")
    exit()

# on vérifie si l'utilisateur a appelé la section Help, et le cas échéant quelle section il désire voir :
helpResponse = 0
if re.search("^-h[a-x]?$", sys.argv[1]):
        if sys.argv[1] == "-h":
            print("What do you need help for ?")
            print("FLAG field ? -hf \nCIGAR field ? -hc\nSystem requirements ? -hr\nSAM file ? -hs")
            helpResponse = input()
        if sys.argv[1] == "-hf" or helpResponse == "-hf":helpFLAG()
        if sys.argv[1] == "-hc" or helpResponse == "-hc":helpCigar()
        if sys.argv[1] == "-hr" or helpResponse == "-hr":helpRequirements()
        if sys.argv[1] == "-hs" or helpResponse == "-hs":helpSAM()
        if sys.argv[1] == "-hp" or helpResponse == "-hp":helpProgramm()

# *******************
# ** FILE ANALYSIS **
# *******************

csv.field_size_limit(sys.maxsize)   # augmentation de la taille max des champs dans le fichier csv afin de pouvoir analyser des alignements de reads > 131kb

## Définition d'une fonction pour mettre le FLAG en binaire

def flagBinary(flag):
        flagB = bin(int(flag))  # Transform the integer into a binary.
        flagB = flagB[2:]   # Remove '0b' Example: '0b1001101' > '1001101'
        flagB = list(flagB) 
        if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
            add = 12 - len(flagB)   # We compute the difference between the maximal flag size (12) and the length of the binary flag.
            for m in range(add):
                flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
        return(flagB)

## 1. Initialisation des éléments non modifiés par la suite
# matrices :
Header_Line = numpy.array([['VN','SO','GO','SS'],['Format version','Sorting order of alignments','Grouping of alignments','Sub-sorting order of alignments']])
Reference_sequence_dictionary = numpy.array([['SN','LN','AH','AN','AS','DS','M5','SP','TP','UR'],['Reference sequence name','Reference sequence length','Alternate locus','Alternative reference sequence names','Genome assembly identifier','Description','MD5 checksum of the sequence','Species','Molecule topology','URI of the sequence']])
Read_group = numpy.array([['ID','BC','CN','DS','DT','FO','KS','LB','PG','PI','PL','PM','PU','SM'],['Read group identifier','Barcode sequence','Name of sequencing center','Description','Date the run was produced','Flow order','Array of nucleotide bases','Library','Processing programs','Predicted median insert size','Platform/technology','Platform model','Platform unit','Sample']])
Program = numpy.array([['ID','PN','CL','PP','DS','VN'],['Program record identifier','Program name','Command line','Previous @PG-ID','Description','Program version']])
Comments = numpy.array([['CO'],['Commentaire(s)']])
# dictionnaires :
infoHeaders = {0:Header_Line,1:Reference_sequence_dictionary,2:Read_group,3:Program,4:Comments,5:'Header_Line',6:'Reference_sequence_dictionary',7:'Read_group',8:'Program',9:'Comments'}
# compilation d'expression régulière :
regex = re.compile('^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$')
# tuple :
headers = ('@HD','@SQ','@RG','@PG','@CO')

## 2. Gestion des options utilisateur
# on cherche le nombre de fichiers donnés en input :
fileNumber = 0
for o in range(len(sys.argv)):
    if len(sys.argv[o]) > 2 and re.search('.sam$', sys.argv[o]):
        fileNumber += 1
# on laisse le choix à l'utilisateur de demander les données dans le terminal, ou écrites dans un fichier output (dont le nom est donné ou non), ou les 2 :
options = len(sys.argv)-fileNumber-1
if len(sys.argv)-(fileNumber+1) <= 0:
    answer = input("Error : no options were given, '-hp' for help about how to run the programm.\n")
    if answer == '-hp': helpProgramm()
    else: exit()
if sys.argv[fileNumber+1] == "-t":
    if len(sys.argv)-(fileNumber+2) > 0 and sys.argv[fileNumber+2] == "-o":
        options = 1
    else: options = 0

for samFile in range(1,fileNumber+1):
    if options != 0:
        summary_data_file = open("summary_data_file"+str(samFile)+".txt", "w")
        summary_data_file.write('Informations :\n\n')
    
    ## 1. Vérification que l'argument est un fichier SAM non vide

    # on vérifie que l'argument est un fichier :
    if os.path.isfile(sys.argv[samFile]) ==  False:
        print("File error : L'argument n'est pas un fichier.")
        exit()
    # on teste si le fichier est vide :
    elif os.path.getsize(sys.argv[samFile]) == 0:
        print("File error : Le fichier est vide.")
        exit()
    # on check l'extension du fichier :
    elif sys.argv[samFile].endswith(".sam") == False: 
        print("Format error : seul le format SAM est accepté.")
        exit()

    ## 2. Préparation de l'analyse du fichier

    ## 2.1 Ouverture et lecture du fichier

    # on lit le fichier dans une variable "file", on indique que les champs sont séparés par des tabulations :
    fi = open(sys.argv[samFile])
    file = csv.reader(fi, delimiter = '\t', quoting = csv.QUOTE_NONE)
    fi.close

    ## 2.2 (ré)initialisations

    # variables :
    i = 0
    j = 0
    analyseComplete = 0
    unmapped_count = 0
    partially_mapped_count = 0
    totally_mapped_count = 0

    # dictionnaires :
    dicoCigar = {'M':0,'I':0,'D':0,'S':0,'H':0,'N':0,'P':0,'X':0,'=':0}

    ## 4. Analyse du fichier

    # analyse ligne par ligne du SAM :
    for ligne in file:
        i += 1   # incrémentation du numéro de ligne à chaque tour de boucle
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

                if int(flag[-2]) == 1:  # si '0010' est présent dans le binaire ça correspond à un FLAG de 2 (= read bien aligné), donc à analyser
                    if re.match('(?![0-9]+M$)', ligne[5]):  # dans le CIGAR, vérifier que la RE "XXXM" (= 100% de match) n'est pas présente
                        partially_mapped_count += 1
                    else:
                        totally_mapped_count +=1

                    expCigar = re.findall('[0-9]+\D', ligne[5])  # on isole tous les champs du CIGAR
                    for champ in range(len(expCigar)):
                        temp = expCigar[champ]
                        dicoCigar[temp[-1]] += int(temp[:(len(temp)-1)])

            elif j == 0:    # si k ne vaut pas 11 (= présence d'une erreur dans le fichier)
                j = j + 1   # on ne veut passer qu'une fois dans cette condition
                k = 11 - k  # calcul du nombre d'erreurs trouvées sur la ligne

                # l'utilisateur a le choix de chercher toutes les erreurs du fichier ou de quitter le programme :
                analyseComplete = input(Back.RED+"Document non analysable"+Back.RESET+" : "+str(k)+" erreur(s) d'expressions régulières trouvée(s) à la ligne "+str(i)+", souhaitez-vous rechercher les erreurs dans le reste du fichier ?\ny/n\n")
                if analyseComplete == "n":
                    exit()
                elif analyseComplete != "y":
                    print(Fore.RED+"Input error : arrêt du script"+Fore.RESET)
                    exit()

        # si la ligne anlysée fait partie du header on extrait les informations qu'on stocke dans le fichier "summary_data_file" :
        else:
            for h in range(4):
                if ligne[0] == headers[h]:
                    if options == 0 or options == 1:
                        print("\n"+headers[h]+" - "+str(infoHeaders[h+5]))
                    if options != 0:
                        summary_data_file.write(headers[h]+" - "+str(infoHeaders[h+5])+"\n")
                    for z in range(1,len(ligne)):
                        for y in range(len(infoHeaders[h])):
                            if infoHeaders[h][0,y] == ligne[z][:2]:
                                if options == 0 or options == 1:
                                    print(str(infoHeaders[h][1,y])+" : "+str(ligne[z][3:]))
                                if options != 0:
                                    summary_data_file.write(str(infoHeaders[h][1,y])+" : "+str(ligne[z][3:])+"\n")  
                    if options != 0:
                        summary_data_file.write("\n")

    # si l'utilisateur avait choisi de chercher toutes les erreurs du document : on arrête la recherche à la fin du fichier et on quitte
    if analyseComplete == "y":
        print(Fore.YELLOW+"Fin de la recherche d'erreur."+Fore.RESET)
        exit()

    if options != 0:
        summary_data_file.write("**********************************\ntotal reads count : "+str(i))
        summary_data_file.write("\n\t-> properly aligned reads count : "+str(totally_mapped_count+partially_mapped_count)+" ("+str(round((totally_mapped_count+partially_mapped_count)*100/i))+"% of total reads)")
        summary_data_file.write("\n\t\t-> totally mapped reads count : "+str(totally_mapped_count))
        summary_data_file.write("\n\t\t-> partially mapped reads count : "+str(partially_mapped_count))
        summary_data_file.write("\n\t-> unmapped reads count : "+str(unmapped_count)+"\n")

    vtotal = 0
    for v in dicoCigar:
        vtotal += dicoCigar[v]
    for w in dicoCigar.keys():
        dicoCigar[w] = round((dicoCigar[w]*100)/vtotal, 2)

    if options != 0:
        summary_data_file.write("\n**********************************\nGlobal CIGAR mutations observed on well-aligned sequences:\n\n"
                        +"Alignment Match : "+str(round(dicoCigar['M'],2))+"%\n"
                        +"Insertion : "+str(round(dicoCigar['I'],2))+"%\n"
                        +"Deletion : "+str(round(dicoCigar['D'],2))+"%\n"
                        +"Skipped region : "+str(round(dicoCigar['S'],2))+"%\n"
                        +"Soft Clipping : "+str(round(dicoCigar['H'],2))+"%\n"
                        +"Hard Clipping : "+str(round(dicoCigar['N'],2))+"%\n"
                        +"Padding : "+str(round(dicoCigar['P'],2))+"%\n"
                        +"Sequence Match : "+str(round(dicoCigar['X'],2))+"%\n"
                        +"Sequence Mismatch : "+str(round(dicoCigar['='],2))+"%\n")
        print(Fore.YELLOW+"\nSummary file for "+str(sys.argv[samFile])+" created."+Fore.RESET)
        summary_data_file.close()

    if options != 0 and len(sys.argv[-fileNumber+samFile-1]) > 2:
        os.rename("summary_data_file"+str(samFile)+".txt", sys.argv[-fileNumber+samFile-1])

    if options == 0 or options == 1:
        print("\ntotal reads count : "+str(i))
        print("\t-> properly aligned reads count : "+Fore.GREEN+str(totally_mapped_count+partially_mapped_count)+Fore.RESET+" ("+str(round((totally_mapped_count+partially_mapped_count)*100/i))+"% of total reads)")
        print("\t\t-> totally mapped reads count : "+str(totally_mapped_count))
        print("\t\t-> partially mapped reads count : "+str(partially_mapped_count))
        print("\t-> unmapped reads count : "+str(unmapped_count))
        print("\nCIGAR analysis (percentages of total properly aligned reads):\n"
                    +"Alignment Match : "+str(round(dicoCigar['M'],2))+"%\n"
                    +"Insertion : "+str(round(dicoCigar['I'],2))+"%\n"
                    +"Deletion : "+str(round(dicoCigar['D'],2))+"%\n"
                    +"Skipped region : "+str(round(dicoCigar['S'],2))+"%\n"
                    +"Soft Clipping : "+str(round(dicoCigar['H'],2))+"%\n"
                    +"Hard Clipping : "+str(round(dicoCigar['N'],2))+"%\n"
                    +"Padding : "+str(round(dicoCigar['P'],2))+"%\n"
                    +"Sequence Match : "+str(round(dicoCigar['X'],2))+"%\n"
                    +"Sequence Mismatch : "+str(round(dicoCigar['='],2))+"%\n")
