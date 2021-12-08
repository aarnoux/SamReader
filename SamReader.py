#!/usr/bin/env python3
#-*- coding utf-8 -*-

__author__ = ("Alizée ARNOUX")
__contact__ = ("alizee.arnoux@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/14/2021"
__licence__ ="This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."

import os, sys, csv, re, colorama, numpy
from colorama import Fore, Back

# Matrix
mHEADER_LINE = numpy.array([['VN','SO','GO','SS'],['Format version','Sorting order of alignments','Grouping of alignments','Sub-sorting order of alignments']])
mREF_SEQ_DICTIONARY = numpy.array([['SN','LN','AH','AN','AS','DS','M5','SP','TP','UR'],['Reference sequence name','Reference sequence length','Alternate locus','Alternative reference sequence names','Genome assembly identifier','Description','MD5 checksum of the sequence','Species','Molecule topology','URI of the sequence']])
mREAD_GROUP = numpy.array([['ID','BC','CN','DS','DT','FO','KS','LB','PG','PI','PL','PM','PU','SM'],['Read group identifier','Barcode sequence','Name of sequencing center','Description','Date the run was produced','Flow order','Array of nucleotide bases','Library','Processing programs','Predicted median insert size','Platform/technology','Platform model','Platform unit','Sample']])
mPROGRAM = numpy.array([['ID','PN','CL','PP','DS','VN'],['Program record identifier','Program name','Command line','Previous @PG-ID','Description','program version']])
mCOMMENTS = numpy.array([['CO'],['Commentaire(s)']])
mCIGAR_MATRIX = numpy.array([['M','I','D','N','S','H','P','=','X'],['Alignement Match','Insertion','Deletion','Skipped region','Soft Clipping','Hard Clipping','Padding','Sequence Match','Sequence Mismatch']])
mQUAL_INTERPRET = numpy.array([['!','"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I'],['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40']])

# Dictionary
dHEADER_FIELD_CALL = {0:mHEADER_LINE,1:mREF_SEQ_DICTIONARY,2:mREAD_GROUP,3:mPROGRAM,4:mCOMMENTS}
dHEADER_TITLE = {0:'Header line',1:'Reference sequence dictionary',2:'Read group',3:'Program',4:'Comments'}

# Regular expression compile
rEXPRESSION = re.compile('^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$')

# Tuple
tHEADERS = ('@HD','@SQ','@RG','@PG','@CO')

# Constant
cARGUMENTS_LIST = sys.argv
cSCRIPT_CALL = 1
cMIN_LINE_LENGHT = 11

# Increasing the maximal size of csv fields in order to be able to analyze reads > 131kb
csv.field_size_limit(sys.maxsize)

# Help functions
def helpFlag():
    print("\nCombination of bitwise FLAGs. Each bit is explained in the following table:\n\n"
            +"Bit\t\tDescription\n"
            +"----------------------------------------------------------------------\n"
            +"1\t0x1\ttemplate having multiple segments in sequencing\n"
            +"2\t0x2\teach segment properly aligned according to the aligner\n"
            +"4\t0x4\tsegment unmapped\n"
            +"8\t0x8\tnext segment in the template unmapped\n"
            +"16\t0x10\tSEQ being reverse complemented\n"
            +"32\t0x20\tSEQ of the next segment in the template being reverse complemented\n"
            +"64\t0x40\tthe first segment in the template\n"
            +"128\t0x80\tthe last segment in the template\n"
            +"256\t0x100\tsecondary alignment\n"
            +"512\t0x200\tnot passing filters, such as platform/vendor quality controls\n"
            +"1024\t0x400\tPCR or optical duplicate\n"
            +"2048\t0x800\tsupplementary alignment")
    exit()

def helpRequirements():
    print("\nThis program runs with "+Fore.RED+"python 3"+Fore.RESET+", and needs the following packages:\n"
            +Fore.YELLOW+"os\tsys\tcsv\tre\tnumpy\tcolorama"+Fore.RESET+"\n"
            +"You can install them using the following commands:\n\n"
            +"$ conda install <package name>\nor\n$ pip install <package name>")
    exit()

def helpCigar():
    print("\nThe CIGAR operations are given in the following table (set ‘*’ if unavailable):\n\n"
            +"Op\tBAM\tDescription\n"
            +"----------------------------------------------------------------------\n"
            +"M\t0\talignment match (can be a sequence match or mismatch)\n"
            +"I\t1\tinsertion to the reference\n"
            +"D\t2\tdeletion from the reference)\n"
            +"N\t3\tskipped region from the reference\n"
            +"S\t4\tsoft clipping (clipped sequences present in SEQ)\n"
            +"H\t5\thard clipping (clipped sequences NOT present in SEQ)\n"
            +"P\t6\tpadding (silent deletion from padded reference)\n"
            +"=\t7\tsequence match\n"
            +"X\t8\tsequence mismatch")
    exit()

def helpSAM():
    print("\nSAM stands for Sequence Alignment/Map format. It is a TAB-delimited text format consisting of a header"
            +"section, which is optional, and an alignment section. If present, the header must be prior to the alignments. "
            +"Header lines start with ‘@’, while alignment lines do not. Each alignment line has 11 mandatory fields for "
            +"essential alignment information such as mapping position, and variable number of optional fields for flexible "
            +"or aligner specific information. Mandatory fields are presented in the following table :\n\n"
            +"Col\tField\tType\tRegex/Range\t\t\tDescription\n"
            +"----------------------------------------------------------------------\n"
            +"1\tQNAME\tString\t[!-?A-~]{1,254}\t\t\tQuery template NAME\n"
            +"2\tFLAG\tInt\t[0,2e16 − 1]\t\t\tbitwise FLAG\n"
            +"3\tRNAME\tString\t\*|[:rname:∧ *=][:rname:]*\tReference sequence NAME\n"
            +"4\tPOS\tInt\t[0,2e31 − 1]\t\t\t1-based leftmost mapping POSition\n"
            +"5\tMAPQ\tInt\t[0,2e8 − 1]\t\t\tMAPping Quality\n"
            +"6\tCIGAR\tString\t\*|([0-9]+[MIDNSHPX=])+\t\tCIGAR string\n"
            +"7\tRNEXT\tString\t\*|=|[:rname: *=][:rname:]*\tReference name of the mate/next read\n"
            +"8\tPNEXT\tInt\t[0,2e31 − 1]\t\t\tPosition of the mate/next read\n"
            +"9\tTLEN\tInt\t[−2e31 + 1,2e31 − 1]\t\tobserved Template LENgth\n"
            +"10\tSEQ\tString\t\*|[A-Za-z=.]+\t\t\tsegment SEQuence\n"
            +"11\tQUAL\tString\t[!-~]+\t\t\t\tASCII of Phred-scaled base QUALity+33")
    exit()

def helpProgram():
    print("\nSamReader is a small program that analyses SAM files. It takes one or more SAM file in input and gives an output either in a file or directly in the terminal. The syntaxe is as follow :\n"
            +"$ /path/to/SamReader.py <input file 1> <input file 2> option <output file 1> <output file 2>\n\n"
            +"Choosing from one of the following options is MANDATORY :\n"
            +"-s	show the results in the terminal, without saving them in a file\n"
            +"-o	save the results in a file which name is given by the user, or in the default file 'summary_data_file.txt'\n"
            +"-s -o	show the results int he terminal AND save them in a file\n"
            +"\nEnter "+Fore.YELLOW+"$ /path/to/SamReader.py -h"+Fore.RESET+" for help")
    exit()

# The help() function analyses the input of the user and call the corresponding help function
def help():
    helpQuery = "NULL"
    if re.search("^-h[a-x]?$", sys.argv[1]): # Check if the user is calling the help function
        if sys.argv[1] == "-h":
            print("What do you need help with ?\nFLAG field ? -hf \nCIGAR field ? -hc\nSystem requirements ? -hr\nSAM file ? -hs\nSamReader ? -hp\n")
            helpQuery = input("Please enter an option below:\n")
        if sys.argv[1] == "-hf" or helpQuery == "-hf":helpFlag()
        if sys.argv[1] == "-hc" or helpQuery == "-hc":helpCigar()
        if sys.argv[1] == "-hr" or helpQuery == "-hr":helpRequirements()
        if sys.argv[1] == "-hs" or helpQuery == "-hs":helpSAM()
        if sys.argv[1] == "-hp" or helpQuery == "-hp":helpProgram()

# Function that verifies that the files given in input are non-empty SAM files, and returns the number of said files
def checkInput(cARGUMENTS_LIST):
    fileNumber = cSCRIPT_CALL
    for arguments in range(len(cARGUMENTS_LIST)):
        if cARGUMENTS_LIST[arguments] == ("-s" or "-o"):    # extraction of the arguments index corresponding to the first option chosen
            for file in range(1,arguments): # for files that are given as input (= put before the first option)
                if os.path.isfile(cARGUMENTS_LIST[file]) ==  False:
                    print("File error : the input is not a file.")
                    exit()
                elif cARGUMENTS_LIST[file].endswith(".sam") == False: 
                    print("Format error : only SAM file format is accepted as input.")
                    exit()
                elif os.path.getsize(cARGUMENTS_LIST[file]) == 0:
                    print("File error : the file is empty.")
                    exit()
                else: fileNumber += 1
    return fileNumber

# Function that lets the user choose the output mode
def userOptionChoice(fileNumber):
    option = "NULL"
    if len(cARGUMENTS_LIST)-(fileNumber) <= 0:   # not choosing any option returns an error
        if input("Error : no options were given, '-hp' for help about how to run the program.\n") == '-hp': helpProgram()
        else: exit()
    if cARGUMENTS_LIST[fileNumber] == "-s": option = "s" # -s means that the results will only be shown in the terminal
    if len(cARGUMENTS_LIST)-(fileNumber+1) > 0 and cARGUMENTS_LIST[fileNumber+1] == "-o": option = "s+o"
    elif cARGUMENTS_LIST[fileNumber] == "-o": option = "o" # -o means that the results will only be saved in a file
    if option == "NULL":
        if input("Option error : a non-authorized option was given, '-hp' for help about how to run the program.\n") == '-hp': helpProgram()   # An unauthorized input ends the script
        else: exit()
    return option

# Function that reads the SAM file as a csv file, with tabulations delimiters, and returns it to be analyzed 
def fileHandler(samFile):
    fi = open(cARGUMENTS_LIST[samFile])
    file = csv.reader(fi, delimiter = '\t', quoting = csv.QUOTE_NONE)
    fi.close
    return file

# Function that returns the FLAG once converted in binary
def flagBinary(flag) :
    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB

# Function that analyzes the file
def fileAnalysis(file, option, samFile, fileNumber):
    # Variables that need to be reinitialized for each file
    i = 0
    headerCount = 0
    unmappedCount = 0
    badlyMappedCount = 0
    totallyMappedCount = 0
    cigarTotalCount = 0
    researchQuery = False
    errorSearch = "NULL"
    mutatedBaseList = []
    mutationsList = []
    qScoreList = []
    dicoSubstitutions = {}
    dicoCigar = {}

    if option != "s":
        outputFile = open("outputFile"+str(samFile)+".txt", "w")
        outputFile.write(str(sys.argv[samFile])+'\n\nInformations :\n')
    if option != "o":
        print("\nAnalyzing :\n"+Back.WHITE+Fore.BLACK+str(sys.argv[samFile])+Fore.RESET+Back.RESET)

    for line in file:
        i += 1
        ERROR_COUNT = 11    # 11 car il y a 11 champs obligatoires dans un fichier SAM

        # A first sorting is made on whether or not the line is from the header section or not
        if line[0].startswith(tHEADERS) == True:
            headerCount += 1
            for field in range(len(tHEADERS)):
                if line[0] == tHEADERS[field]:   # Identifying the main section title
                    if option != "o":
                        print("\n"+tHEADERS[field]+" - "+dHEADER_TITLE[field])
                    elif option != "s":
                        outputFile.write("\n"+tHEADERS[field]+" - "+dHEADER_TITLE[field]+"\n")
                    for headerSubField in range(1,len(line)):
                        for m in range(len(dHEADER_FIELD_CALL[field][0])):
                            if dHEADER_FIELD_CALL[field][0,m] == line[headerSubField][:2]:  # Identifying subsections titles by browsing the first row of the relevant matrix
                                
                                # When the subsection is found, print or write the extended description from the second row of the relevant matrix, and the corresponding informations of the analyzed SAM file
                                if option != "o":
                                    print(str(dHEADER_FIELD_CALL[field][1,m])+" : "+str(line[headerSubField][3:]))
                                if option != "s":
                                    outputFile.write(str(dHEADER_FIELD_CALL[field][1,m])+" : "+str(line[headerSubField][3:]+"\n"))
        else:
            # Error searching begins by researching the regular expression (RE) from the fields, and then:
            # - if the RE is correct, the counter 'ERROR_COUNT' is decreased by 1
            # - if the RE is incorrect, an error message is printed but the line continues to be analyzed
            if (re.search('^[!-?A-~]{1,254}$', line[0])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"QNAME field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ QNAME a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[0])+"\n"+str(line))
            if (re.search('^[0-9]{1,4}$', line[1])):
                ERROR_COUNT -= 1
            else:
                print(Fore.RED+"FLAG field error"+Fore.RESET+" : Erreur de contenu à la ligne "+str(i)+" , le champ FLAG a une expression régulière non conforme, ou un des champ de la line n'est pas délimité par une tabulation.\nChamp non conforme : "+str(line[1])+"\n"+str(line))
            if (rEXPRESSION.search(line[2]) or re.search('\*', line[2])):
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
            if (rEXPRESSION.search(line[6]) or re.search('^\*$', line[6]) or re.search('^=$', line[6])):
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

            # if there is no errors on the line, ERROR_COUNT equals 0
            if ERROR_COUNT == 0:
                flag = flagBinary(line[1])  # the flag is translated into binary

                # if '0100' is present in the flag field, it means that a '4' is in the flag (= unmapped read)
                if int(flag[-3]) == 1: unmappedCount += 1

                # if '0010' is present in the flag field, it means that a '2' is in the flag (= properly aligned read), a flag of 0 doesn't mean that the flag is unmapped so it has to be taken into account
                else:

                    # in the CIGAR field, extract the badly mapped reads
                    if re.match('(?![0-9]+M$)', line[5]): badlyMappedCount += 1 
                    else:
                        totallyMappedCount +=1
                        
                        # from the totally mapped reads, search for the optional "MD" field
                        if len(line)-cMIN_LINE_LENGHT > 0:
                            for field in range(cMIN_LINE_LENGHT,len(line)):
                                baseCounter = 0
                                if line[field][:5] == "MD:Z:" and re.match('(?![0-9]+$)', line[field][5:]): # check that the "MD" field is present and the sequence match is not 100%
                                    mutations = re.findall('[0-9]+\D', line[field])
                                    for substitution in range(len(mutations)):
                                        matchedBase = int(mutations[substitution][:-1])    # if so, extract all mutations
                                        if int(line[8]) < 0:    # test whether the mutation is present on the complementary read
                                            mutatedBaseList.append(int(line[3])-baseCounter-matchedBase)    # if the mutation is on the complementary read, substract the mutation position from the beginning of the read alignement
                                        else:
                                            mutatedBaseList.append(int(line[3])+baseCounter+matchedBase)    # if the mutation is on the forward read, add the mutation position to the beginning of the read alignement
                                        mutationsList.append(str(line[9][baseCounter+matchedBase])+" -> "+str(mutations[substitution][-1]))  # store the mutation in the form 'N -> N'
                                        qScoreList.append(str(line[10][baseCounter+matchedBase]))   #   store the QUAL value of substitued base
                                        baseCounter += matchedBase + 1

                    champsCigar = re.findall('[0-9]+\D', line[5])   # find all CIGAR information
                    for champ in range(len(champsCigar)):
                        if champsCigar[champ][-1] in dicoCigar.keys():    # check if there is already an entry for the mutation in the dictionary
                            dicoCigar[champsCigar[champ][-1]] += int(champsCigar[champ][:-1])
                        else:
                            dicoCigar[champsCigar[champ][-1]] = int(champsCigar[champ][:-1])
            
            # if there are errors on the line, pass 1 time in the following condition
            elif researchQuery == False:
                researchQuery = True
                errorSearch = input(Back.RED+"Document non analysable"+Back.RESET+" : "+str(ERROR_COUNT)+" erreur(s) d'expressions régulières trouvée(s) à la line "+str(i)+", souhaitez-vous rechercher les erreurs dans le reste du fichier ?\ny/n\n")
                if errorSearch == "n":
                    exit()
                elif errorSearch != "y":
                    print(Fore.RED+"Input error : exit."+Fore.RESET)
                    exit()

    # if the user wanted to search the entire file, after the error search is complete exit the script
    if errorSearch == "y":
        print(Fore.YELLOW+"End of error research."+Fore.RESET)
        exit()

    if totallyMappedCount + badlyMappedCount == 0:
        print(Fore.RED+"No reads could be analyzed"+Fore.RESET)
        exit()

    for value in dicoCigar: cigarTotalCount += dicoCigar[value]

    for mutation in mutationsList:
        if mutation in dicoSubstitutions.keys():
            dicoSubstitutions[mutation] += 1
        else:
            dicoSubstitutions[mutation] = 1
    tDicoSub = sorted(dicoSubstitutions.items(),key=lambda x:x[1], reverse=True)    # sort the dictionary in descending order

    if option != "s":
        outputFile.write("\n**********************************\n\ntotal reads count : "+str(i-headerCount)
            +"\n\t-> aligned reads count : "+str(totallyMappedCount+badlyMappedCount)+" ("+str(round((totallyMappedCount+badlyMappedCount)*100/(i-headerCount),2))+"% of total reads)"
            +"\n\t\t-> totally mapped reads count : "+str(totallyMappedCount)
            +"\n\t\t-> partially mapped reads count : "+str(badlyMappedCount)
            +"\n\t-> unmapped reads count : "+str(unmappedCount)
            +"\n\n**********************************\n\nGlobal CIGAR mutations observed on aligned sequences:\n\n")
        for key in dicoCigar.keys():
            for x in range(len(mCIGAR_MATRIX[0,])):
                if key == mCIGAR_MATRIX[0,x]:
                    outputFile.write(str(mCIGAR_MATRIX[1,x])+" : "+str(round((dicoCigar[key]*100/cigarTotalCount), 4))+"%     ("+str(dicoCigar[key])+" out of "+str(cigarTotalCount)+" mutations)\n")
        outputFile.write("\n**********************************\n\nSummary of nucleotide substitutions :\n\nSubstitution\t\tIteration\n------------------------------------\n")
        for x in range(len(tDicoSub)):
            outputFile.write(str(tDicoSub[x][0])+"\t\t\t"+str(tDicoSub[x][1])+"\n")
        outputFile.write("\n**********************************\n\nList of substitutions in totally mapped reads :\n\nNucleotide N°\t\tMutation\t\tBase call accuracy (%)\n----------------------------------------------------------------------\n")
        for x in range(len(mutationsList)):
            for qScoreListValue in range(len(mQUAL_INTERPRET[0,])):
                if qScoreList[x] == mQUAL_INTERPRET[0,qScoreListValue]:
                    quality = int(mQUAL_INTERPRET[1,qScoreListValue])
                    outputFile.write(str(mutatedBaseList[x])+"\t\t\t"+str(mutationsList[x])+"\t\t\t"+str(round((1-(10**(-quality/10)))*100,2))+"\n")

        print(Fore.YELLOW+"\nOutput file for "+str(cARGUMENTS_LIST[samFile])+" created."+Fore.RESET)
        outputFile.close()

    if option != "o":
        print("\ntotal reads count : "+str(i-headerCount)
            +"\n\t-> aligned reads count : "+str(totallyMappedCount+badlyMappedCount)+" ("+str(round((totallyMappedCount+badlyMappedCount)*100/(i-headerCount),2))+"% of total reads)"
            +"\n\t\t-> totally mapped reads count : "+str(totallyMappedCount)
            +"\n\t\t-> partially mapped reads count : "+str(badlyMappedCount)
            +"\n\t-> unmapped reads count : "+str(unmappedCount)
            +"\n\nMutations analysis:")
        for key in dicoCigar.keys():
            for x in range(len(mCIGAR_MATRIX[0,])):
                if key == mCIGAR_MATRIX[0,x]:
                    print(str(mCIGAR_MATRIX[1,x])+" : "+str(round((dicoCigar[key]*100/cigarTotalCount), 4))+"%     ("+str(dicoCigar[key])+" out of "+str(cigarTotalCount)+" mutations)")
        for mut in mutationsList:
            if mut in dicoSubstitutions.keys():
                dicoSubstitutions[mut] += 1
            else:
                dicoSubstitutions[mut] = 1
        print("\nSummary of nucleotide substitutions :\nSubstitution\t\tIteration\n------------------------------------")
        for x in dicoSubstitutions:
            print(str(x)+"\t\t\t"+str(dicoSubstitutions[x]))

    # rename the output files if the user gave specific names
    if option == "o" and len(cARGUMENTS_LIST) > (fileNumber + samFile) and len(cARGUMENTS_LIST[fileNumber + samFile]) > 2:
        os.rename("outputFile"+str(samFile)+".txt", cARGUMENTS_LIST[samFile + fileNumber]+".txt")
    elif option == "s+o" and len(cARGUMENTS_LIST) > (fileNumber + samFile + 1) and len(cARGUMENTS_LIST[fileNumber + samFile + 1]) > 2:
        os.rename("outputFile"+str(samFile)+".txt", cARGUMENTS_LIST[samFile + fileNumber + 1]+".txt")

# main function
def main(cARGUMENTS_LIST):
    if len(cARGUMENTS_LIST) == 1: helpProgram()
    help()
    fileNumber = checkInput(cARGUMENTS_LIST)
    option = userOptionChoice(fileNumber)
    for samFile in range(1,fileNumber):
        file = fileHandler(samFile)
        fileAnalysis(file, option, samFile, fileNumber)

if __name__ == "__main__":
    main(sys.argv)
