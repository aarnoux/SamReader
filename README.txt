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
- Quantification and classification of reads, read pairs, and gap or overlap within read pairs
- Mutation analysis includes:
	-> quantification of each type of mutations possible from CIGAR analysis
	For substitution mutations:
	-> calculation of the number of each substitution possibilities (i.e. A->T, A->C, A->G, etc)
	-> list of all the substitution found in the query sequence relative to the reference, with the quality of
	base calling, and whether or not the mutation is synonymous for each possible open reading frames (ORF)
- Output  in a file (-o) and a CSV file containing the list described above is created independently of the option choosen.

ATTENTION: The three ORF are arbitrarily defined for each mutations, as such you CAN NOT compare the results between different lines from the .csv file.


**** EXECUTION ****

Make sure that you run the version 3 of Python.
SamReader requires the following modules to function correctly:
-re
-numpy
-os
-sys
-csv

You can install them with the following command:
- Linux/MacOS:
$ pip install <package name>
- Windows:
> python -m pip install <package name>

This extra step is sometimes necessary on MacOS, does not concern Linux or Windows users:
Add execution permission to the script:
$ chmod +x ~/path/to/file/SamReader.py

Go to the directory containing the script, or add the path to your .bashrc or equivalent, then run the following command in the terminal to execute SamReader (output file field is optionnal):
- Linus/MacOS:
$ SamReader.py <input-file.sam> -options <output-file.txt>
- Windows:
> py SamReader.py <input-file.sam> -options <output-file.txt>

SamReader can analyse multiple files at once (output files fields are optionnals):
- Linus/MacOS:
$ SamReader.py <input-file1.sam> <input-file2.sam> -options <output-file1.txt> <outpute-file2.txt>
- Windows:
> py SamReader.py <input-file1.sam> <input-file2.sam> -options <output-file1.txt> <outpute-file2.txt>



**** OPTIONS ****

Choosing from one of the following options is MANDATORY :
-o		followed by desired names for output files: save the results in said file.
--check		followed by a number or 'all': indicates how much of the file to screen for SAM field errors. Not putting this option results in no checking
		while putting 'all' results in the verification of the totality of the file. The analysis is approx. 30% slower when checking all than when not 		performing any verifications.


**** HELP ****

The commands are shown for Linux/MacOS users, please adapt for Windows user as presented in section "EXECUTION".

$ SamReader.py -h	access help to show the different sections available and choose one of them

Or enter directly in the terminal one of the following terms :
-hf	show details about the FLAG field
$ SamReader.py -hf

-hc	show details about the CIGAR field
$ SamReader.py -hc

-hr	show details about the requirements to run the programm
$ SamReader.py -hr

-hs	show details about the SAM file format
$ SamReader.py -hs

-hp	show details about the programm
$ SamReader.py -hp


**** OUTPUT ****

A CSV file that contains the list of all the substitution found in the query sequence relative to the reference, with the quality of base calling, read ID, and whether or not the mutation is synonymous for each possible reading frames.

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

# informations about the paired reads
- paired read analysis:
unmapped + unmapped: 0.37%   (650 out of 175665 pairs)
unmapped + badly mapped: 0.0085%   (15 out of 175665 pairs)
badly mapped + totally mapped: 0.0609%   (107 out of 175665 pairs)
totally mapped + totally mapped: 99.5605%   (174893 out of 175665 pairs)

- pairs alignement analysis:	# is printed only if the reads are paired
6924 read pairs are well aligned on both ends of the fragment

52732 read pairs present a gap of [1,10[ nucleotide(s) between the forward and the reverse read

55400 read pairs present an overlap of [1,10[ nucleotide(s) between the forward and the reverse read

# mutations informations as given by the CIGAR field (relative to reference sequences)
Mutations analysis:
Alignement Match : 99.9951%     (35000611 out of 35002311 nucleotides)
Deletion : 0.0002%     (61 out of 35002311 nucleotides)
Insertion : 0.0001%     (18 out of 35002311 nucleotides)
Soft Clipping : 0.0046%     (1621 out of 35002311 nucleotides)

# number of each substitution possibilities sorted in descending order (relative to reference sequences)
Summary of nucleotide substitutions :
Substitution		Iteration
------------------------------------
G -> A			5391
T -> A			5363
