#!/usr/bin/env python3
#-*- coding utf-8 -*-

__author__ = ("Alizée ARNOUX")
__contact__ = ("alizee.arnoux@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/14/2021"
__licence__ ="This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."
        
import os, sys, csv, re, colorama, numpy
from colorama import Fore, Back

# Matrices
headerLine = numpy.array([['VN','SO','GO','SS'],['Format version','Sorting order of alignments','Grouping of alignments','Sub-sorting order of alignments']])
ReferenceSequenceDictionary = numpy.array([['SN','LN','AH','AN','AS','DS','M5','SP','TP','UR'],['Reference sequence name','Reference sequence length','Alternate locus','Alternative reference sequence names','Genome assembly identifier','Description','MD5 checksum of the sequence','Species','Molecule topology','URI of the sequence']])
readGroup = numpy.array([['ID','BC','CN','DS','DT','FO','KS','LB','PG','PI','PL','PM','PU','SM'],['Read group identifier','Barcode sequence','Name of sequencing center','Description','Date the run was produced','Flow order','Array of nucleotide bases','Library','Processing programs','Predicted median insert size','Platform/technology','Platform model','Platform unit','Sample']])
program = numpy.array([['ID','PN','CL','PP','DS','VN'],['Program record identifier','Program name','Command line','Previous @PG-ID','Description','program version']])
comments = numpy.array([['CO'],['Commentaire(s)']])

# Dictionnary
infoHeader = {0:headerLine,1:ReferenceSequenceDictionary,2:readGroup,3:program,4:comments,5:'Header line',6:'Reference sequence dictionary',7:'Read group',8:'Program',9:'Comments'}

# Regex
regex = re.compile('^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$')

# Tuple
headers = ('@HD','@SQ','@RG','@PG','@CO')

# Constants
ARGUMENTS_LIST = sys.argv
SCRIPT_CALL = 1
MIN_LINE_LENGHT = 11

# Increasing the maximal size of csv fields in order to be able to analyze reads > 131kb
csv.field_size_limit(sys.maxsize)

#******************
#* HELP FUNCTIONS *
#******************

def helpFlag():
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

def helpProgram():
    print("\nSamReader is a small program that analyses SAM files. It can analyse multiple files in one go.")
    print("Choosing from one of the following options is MANDATORY :\n-t	show the results in the terminal, without saving them in a file\n-o	save the results in a file which name is given by the user, or in the default file 'summary_data_file.txt'\n-t -o	show the results int he terminal AND save them in a file")
    exit()

# Help function analyses the input of the user and call the corresponding help function
def help():
    helpQuery = "NULL"
    if re.search("^-h[a-x]?$", sys.argv[1]): # Check if the user is calling the help function
        if sys.argv[1] == "-h":
            print("What do you need help with ?")
            print("FLAG field ? -hf \nCIGAR field ? -hc\nSystem requirements ? -hr\nSAM file ? -hs\nSamReader ? -hp\n")
            helpQuery = input()
        if sys.argv[1] == "-hf" or helpQuery == "-hf":helpFlag()
        if sys.argv[1] == "-hc" or helpQuery == "-hc":helpCigar()
        if sys.argv[1] == "-hr" or helpQuery == "-hr":helpRequirements()
        if sys.argv[1] == "-hs" or helpQuery == "-hs":helpSAM()
        if sys.argv[1] == "-hp" or helpQuery == "-hp":helpProgram()

#*****************
#* FILE ANALYSIS *
#*****************

# Function that returns the number of arguments given as input (sam files + script calling)
def inputFileNumber(ARGUMENTS_LIST):
    fileNumber = SCRIPT_CALL
    for argument in range(len(ARGUMENTS_LIST)):
        if len(ARGUMENTS_LIST[argument]) > 2 and re.search('.sam$', ARGUMENTS_LIST[argument]): fileNumber += 1
    return fileNumber

# Function that lets the user choose the output mode
def userOptionChoice(fileNumber):
    option = "NULL"
    if len(ARGUMENTS_LIST)-(fileNumber) <= 0:   # Not choosing any options ends the script
        if input("Error : no options were given, '-hp' for help about how to run the program.\n") == '-hp': helpProgram()
        else: exit()
    if ARGUMENTS_LIST[fileNumber] == "-t": option = "t" # -t means that the results will only be shown in the terminal
    if len(ARGUMENTS_LIST)-(fileNumber+1) > 0 and ARGUMENTS_LIST[fileNumber+1] == "-o": option = "t+o"
    elif ARGUMENTS_LIST[fileNumber] == "-o": option = "o" # -o means that the results will only be saved in a file
    if option == "NULL":
        if input("Option error : a non-authorized option was given, '-hp' for help about how to run the program.\n") == '-hp-': helpProgram()   # An unauthorized input ends the script
        else: exit()
    return option

def checkFileIntegrity(samFileNumber):
    if os.path.isfile(ARGUMENTS_LIST[samFileNumber]) ==  False:
        print("File error : L'argument n'est pas un fichier.")
        exit()
    elif os.path.getsize(ARGUMENTS_LIST[samFileNumber]) == 0:
        print("File error : Le fichier est vide.")
        exit()
    elif ARGUMENTS_LIST[samFileNumber].endswith(".sam") == False: 
        print("Format error : seul le format SAM est accepté.")
        exit()

def fileHandler(samFileNumber):
    fi = open(ARGUMENTS_LIST[samFileNumber])
    file = csv.reader(fi, delimiter = '\t', quoting = csv.QUOTE_NONE)
    fi.close
    return file

def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB

def fileAnalysis(file, option, samFileNumber, fileNumber):
    i = 0
    researchQuery = False
    errorSearch = "NULL"
    unmappedCount = 0
    partiallyMappedCount = 0
    totallyMappedCount = 0
    dicoCigar = {'M':0,'I':0,'D':0,'S':0,'H':0,'N':0,'P':0,'X':0,'=':0}
    nucleotideNb = []
    mutation = []

    if option != "t":
        summary_data_file = open("summary_data_file"+str(samFileNumber)+".txt", "w")
        summary_data_file.write(str(sys.argv[samFileNumber])+'\n\nInformations :\n')

    if option == "t" or "t+o":
        print("\nAnalyzing :\n"+Back.WHITE+Fore.BLACK+str(sys.argv[samFileNumber])+Fore.RESET+Back.RESET)

    for line in file:
        i += 1
        ERROR_COUNT = 11    # il y a 11 champs obligatoires dans un fichier SAM

        if line[0].startswith(headers) == True:
            for field in range(len(headers)):
                if line[0] == headers[field]:
                    if option != "o":
                        print("\n"+headers[field]+" - "+infoHeader[field+5])
                    elif option != "t":
                        summary_data_file.write("\n"+headers[field]+" - "+infoHeader[field+5]+"\n")
                    for headerSubField in range(1,len(line)):
                        for m in range(len(infoHeader[field])):
                            if infoHeader[field][0,m] == line[headerSubField][:2]:
                                if option != "o":
                                    print(str(infoHeader[field][1,m])+" : "+str(line[headerSubField][3:]))
                                if option != "t":
                                    summary_data_file.write(str(infoHeader[field][1,m])+" : "+str(line[headerSubField][3:]+"\n"))
        else:
            # la recherche d'erreurs dans les champs débute par la recherche de l'expression régulière (RE) du champ, puis :
            # - si la RE est correcte on diminue le compteur 'ERROR_COUNT' de 1
            # - si la RE est incorrecte, on affiche un message d'erreur et on finit d'analyser la line
            if (re.search('^[!-?A-~]{1,254}$', line[0])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"QNAME field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ QNAME a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[0])+"\n"+str(line))
            if (re.search('^[0-9]{1,4}$', line[1])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"FLAG field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ FLAG a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[1])+"\n"+str(line))
            if (regex.search(line[2]) or re.search('\*', line[2])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"RNAME field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ RNAME a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[2])+"\n"+str(line))
            if (re.search('^[0-9]*$', line[3])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"POS field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ POS a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[3])+"\n"+str(line))
            if (re.search('^[0-5][0-9]$''|^[0-9]$''|^(60)$''|^(255)$', line[4])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"MAPQ field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ MAPQ a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[4])+"\n"+str(line))
            if (re.search('^\*$''|^([0-9]+[MIDNSHPX=])+$', line[5])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"CIGAR field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ CIGAR a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[5])+"\n"+str(line))                             
            if (regex.search(line[6]) or re.search('^\*$', line[6]) or re.search('^=$', line[6])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"RNEXT field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ RNEXT a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[6])+"\n"+str(line))
            if (re.search('^[0-9]*$', line[7])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"PNEXT field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ PNEXT a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[7])+"\n"+str(line))
            if (re.search('^^-?[0-9]{1,30}$''|^-?20{31}$', line[8])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"TLEN field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ TLEN a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[8])+"\n"+str(line))
            if (re.search('^\*|[A-Za-z=.]+$', line[9])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"SEQ field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ SEQ a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[9])+"\n"+str(line))
            if (re.search('^[!-~]+$', line[10])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"QUAL field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ QUAL a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[10])+"\n"+str(line))

            # si la line ne présente pas d'erreurs alors 'ERROR_COUNT' vaut 11
            if ERROR_COUNT == 0:
                flag = flagBinary(line[1])

                # si '0100' est présent dans le binaire ça correspond à un FLAG de 4 (= read non mappé)
                if int(flag[-3]) == 1:
                    unmappedCount += 1

                # si '0010' est présent dans le binaire ça correspond à un FLAG de 2 (= read bien aligné)
                if int(flag[-2]) == 1 or int(flag[-1]) == 0:

                    # dans le CIGAR, vérifier qu'il n'y a pas 100% de match (= read partiellement mappé)
                    if re.match('(?![0-9]+M$)', line[5]): partiallyMappedCount += 1 
                    else:
                        totallyMappedCount +=1
                        if len(line)-MIN_LINE_LENGHT > 0:
                            for c in range(MIN_LINE_LENGHT,len(line)):
                                pointerSeq = 0
                                if line[c][:5] == "MD:Z:" and re.match('(?![0-9]+$)', line[c][5:]):
                                    mutations = re.findall('[0-9]+\D', line[c])
                                    for m in range(len(mutations)):
                                        matchNb = int(mutations[m][:-1])
                                        if int(line[8]) < 0:
                                            nucleotideNb.append(int(line[3])-pointerSeq-matchNb)
                                        else:
                                            nucleotideNb.append(int(line[3])+pointerSeq+matchNb)
                                        mutation.append(str(line[9][pointerSeq+matchNb])+" -> "+str(mutations[m][-1]))
                                        pointerSeq += matchNb

                    champsCigar = re.findall('[0-9]+\D', line[5])
                    for champ in range(len(champsCigar)):
                        temp = champsCigar[champ]
                        dicoCigar[temp[-1]] += int(temp[:(len(temp)-1)])
            
            elif researchQuery == False:
                researchQuery = True
                errorSearch = input(Back.RED+"Document non analysable"+Back.RESET+" : "+str(ERROR_COUNT)+" erreur(s) d'expressions régulières trouvée(s) à la line "+str(i)+", souhaitez-vous rechercher les erreurs dans le reste du fichier ?\ny/n\n")
                if errorSearch == "n":
                    exit()
                elif errorSearch != "y":
                    print(Fore.RED+"Input error : arrêt du script"+Fore.RESET)
                    exit()

    if errorSearch == "y":
        print(Fore.YELLOW+"Fin de la recherche d'erreur."+Fore.RESET)
        exit()

    if option != "t":
        summary_data_file.write("\n**********************************\ntotal reads count : "+str(i))
        summary_data_file.write("\n\t-> aligned reads count : "+str(totallyMappedCount+partiallyMappedCount)+" ("+str(round((totallyMappedCount+partiallyMappedCount)*100/i))+"% of total reads)")
        summary_data_file.write("\n\t\t-> totally mapped reads count : "+str(totallyMappedCount))
        summary_data_file.write("\n\t\t-> partially mapped reads count : "+str(partiallyMappedCount))
        summary_data_file.write("\n\t-> unmapped reads count : "+str(unmappedCount)+"\n")

    if option != "o":
            print("\ntotal reads count : "+str(i))
            print("\t-> aligned reads count : "+Fore.GREEN+str(totallyMappedCount+partiallyMappedCount)+Fore.RESET+" ("+str(round((totallyMappedCount+partiallyMappedCount)*100/i))+"% of total reads)")
            print("\t\t-> totally mapped reads count : "+str(totallyMappedCount))
            print("\t\t-> partially mapped reads count : "+str(partiallyMappedCount))
            print("\t-> unmapped reads count : "+str(unmappedCount))

    if totallyMappedCount + partiallyMappedCount > 0:
        totalValue = 0
        for value in dicoCigar: totalValue += dicoCigar[value]
        for key in dicoCigar.keys(): dicoCigar[key] = round((dicoCigar[key]*100/totalValue), 2)

        if option != "t":
            summary_data_file.write("\n**********************************\nGlobal CIGAR mutations observed on aligned sequences:\n\n"
                        +"Alignment Match : "+str(dicoCigar['M'])+"%\n"
                        +"Insertion : "+str(dicoCigar['I'])+"%\n"
                        +"Deletion : "+str(dicoCigar['D'])+"%\n"
                        +"Skipped region : "+str(dicoCigar['N'])+"%\n"
                        +"Soft Clipping : "+str(dicoCigar['S'])+"%\n"
                        +"Hard Clipping : "+str(dicoCigar['H'])+"%\n"
                        +"Padding : "+str(dicoCigar['P'])+"%\n"
                        +"Sequence Match : "+str(dicoCigar['='])+"%\n"
                        +"Sequence Mismatch : "+str(dicoCigar['X'])+"%\n")
                                
            summary_data_file.write("\n**********************************\nList of substitutions in totally mapped reads :\n\nNucleotide N°\t\tMutation\n----------------------------------\n")
            for x in range(len(mutation)):
                summary_data_file.write(str(nucleotideNb[x])+"\t\t\t"+str(mutation[x])+"\n")

            print(Fore.YELLOW+"\nSummary file for "+str(ARGUMENTS_LIST[samFileNumber])+" created."+Fore.RESET)

            summary_data_file.close()
        
        if option != "o":
            print("\nCIGAR analysis (percentages of total aligned reads):\n"
                        +"Alignment Match : "+str(dicoCigar['M'])+"%\n"
                        +"Insertion : "+str(dicoCigar['I'])+"%\n"
                        +"Deletion : "+str(dicoCigar['D'])+"%\n"
                        +"Skipped region : "+str(dicoCigar['N'])+"%\n"
                        +"Soft Clipping : "+str(dicoCigar['S'])+"%\n"
                        +"Hard Clipping : "+str(dicoCigar['H'])+"%\n"
                        +"Padding : "+str(dicoCigar['P'])+"%\n"
                        +"Sequence Match : "+str(dicoCigar['='])+"%\n"
                        +"Sequence Mismatch : "+str(dicoCigar['X'])+"%\n")

    else:
        if option != "t":
            summary_data_file.write("No reads could be analyzed.")
        if option != "o": print(Fore.RED+"No reads could be analyzed"+Fore.RESET)

    if option == "o" and len(ARGUMENTS_LIST) > (fileNumber + samFileNumber) and len(ARGUMENTS_LIST[fileNumber + samFileNumber]) > 2:
        os.rename("summary_data_file"+str(samFileNumber)+".txt", ARGUMENTS_LIST[samFileNumber + fileNumber]+".txt")
    elif option == "t+o" and len(ARGUMENTS_LIST) > (fileNumber + samFileNumber + 1) and len(ARGUMENTS_LIST[fileNumber + samFileNumber + 1]) > 2:
        os.rename("summary_data_file"+str(samFileNumber)+".txt", ARGUMENTS_LIST[samFileNumber + fileNumber + 1]+".txt")

def main(ARGUMENTS_LIST):
    help()
    fileNumber = inputFileNumber(ARGUMENTS_LIST)
    option = userOptionChoice(fileNumber)
    for samFileNumber in range(1,fileNumber):
        checkFileIntegrity(samFileNumber)
        file = fileHandler(samFileNumber)
        fileAnalysis(file, option, samFileNumber, fileNumber)

if __name__ == "__main__":
    main(sys.argv)
