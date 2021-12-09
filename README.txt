**************************************************************************
      _____   ____    __     __
     /____/| /___/\  /_/\   /_/|
    / ___|/ / __ \ \| \  \ / | |
   / /\___ | ||_| | |  \  /  | |
   \ \/__/\| |/_| | |   \/   | | ___    ___   ___   ___    ___  ___
    \___ \ |  __  | | |\  /| | ||  _ \ |  _| / _ \ |  _ \ |  _||  _ \
     ___\ \| | || | | ||\/ | | || |_| || |_ | |_| || | | || |_ | |_| |
    /___/ /| | || | | ||   | | || __ / | |_ |  _  || |_| || |_ | __ /
   |_____/_|_|/ |_|/|_|/   |_|/ |_| \_\|___||_| |_||____/ |___||_| \_\

**************************************************************************

SamReader v0.0.1

author : Alizée ARNOUX

contact : alizee.arnoux@etu.umontpellier.fr

date : 12/14/2021

licence : This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. This program is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.

This README file is organized as follow :
- GENERAL PURPOSE
- EXECUTION
- OPTIONS
- HELP
- OUTPUT


**** GENERAL PURPOSE ****

SamReader is a programm that analyses SAM files. It was developped with Python 3. It can analyse multiple files
in one go. Results are either shown in the terminal or written in an output file. Features include:
- Verifying the integrity of the files given as arguments.
- Mutation analysis includes:
	-> quantification of each type of mutations possible from CIGAR analysis
	For substitution mutations:
	-> calculation of the number of each substitution possibilities (i.e. A->T, A->C, A->G, etc)
	-> list of all the substitution found in the query sequence relative to the reference, with the quality of
	base calling, and whether or not the mutation is synonymous for each possible open reading frames (ORF)
- Output in the terminal (-s), in a file (-o) or both (-s -o). A CSV file containing the list described above is
created independently of the option choosen.

ATTENTION: The three ORF are arbitrarily defined for each mutations, as such you CAN NOT compare the results between different lines from the .csv file.

**** EXECUTION ****

This extra step is sometime necessary on MacOS, does not concern Linux or Windows users:
Add execution permission to the script:
$ chmod +x ~/path/to/file/SamReader.py

Run the following command in the terminal to execute SamReader (output file field is optionnal):
$ ~/path/to/file/SamReader.py <input-file.sam> -options <output-file.txt>

SamReader can analyse mutliple files at once (output files fields are optionnals):
$ ~/path/to/file/SamReader.py <input-file1.sam> <input-file2.sam> -options <output-file1.txt> <outpute-file2.txt>


**** OPTIONS ****

Choosing from one of the following options is MANDATORY :
-s	show the results in the terminal, without saving them in a file
-o	save the results in a file which name is given by the user, or in the default file "outputFile.txt"
-s -o	show the results in the terminal AND save them in a file


**** HELP ****

$ ~/path/to/file/SamReader.py -h	access help to show the different sections available and choose one of them,
					or enter directly in the terminal one of the following terms :

Sections:
-hf	show details about the FLAG field
$ ~/path/to/file/SamReader.py -hf

-hc	show details about the CIGAR field
$ ~/path/to/file/SamReader.py -hc

-hr	show details about the requirements to run the programm
$ ~/path/to/file/SamReader.py -hr

-hs	show details about the SAM file format
$ ~/path/to/file/SamReader.py -hs

-hp	show details about the programm
$ ~/path/to/file/SamReader.py -hp


**** OUTPUT ****

A CSV file that contains the list of all the substitution found in the query sequence relative to the reference, with the quality of
base calling, read ID, and whether or not the mutation is synonymous for each possible reading frames.
The output in the terminal or in file is presented as follow:


Analyzing :
mapping.sam	# file name, and path if given

#extraction of header sections and of the corresponding informations
Informations :

@SQ - Reference_sequence_dictionary
Reference sequence name : Reference

@PG - Program
Program record identifier : bwa

CSV file for mapping.sam created.
Output file for mapping.sam created.	# if option -o is chosen

#reads informations
total reads count : 351330
	-> aligned reads count : 350015 (99.63% of total reads)
		-> totally mapped reads count : 349893
		-> partially mapped reads count : 122
	-> unmapped reads count : 1315

# mutations informations as given by the CIGAR field
Mutations analysis:
Alignement Match : 99.9951%     (35000611 out of 35002311 nucleotides)
Deletion : 0.0002%     (61 out of 35002311 nucleotides)
Insertion : 0.0001%     (18 out of 35002311 nucleotides)
Soft Clipping : 0.0046%     (1621 out of 35002311 nucleotides)

# number of each substitution possibilities sorted in descending order
Summary of nucleotide substitutions :
Substitution		Iteration
------------------------------------
G -> A			5391
T -> A			5363
