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
- Reads that are unmapped are not taken into account for the calculations.
- Calculations for the percentages of mutations are based on the CIGAR field.


**** EXECUTION ****

Run the following command in the terminal to execute SamReader (output file field is optionnal):
$ ~/path/to/file/SamReader.py <input-file.sam> -options <output-file.txt>

SamReader can analyse mutliple files at once (output files fields are optionnal):
$ ~/path/to/file/SamReader.py <input-file1.sam> <input-file2.sam> -options <output-file1.txt> <outpute-file2.txt>


**** OPTIONS ****

Choosing from one of the following options is MANDATORY :
-t	show the results in the terminal, without saving them in a file
-o	save the results in a file which name is given by the user, or in the default file "summary_data_file.txt"
-t -o	show the results in the terminal AND save them in a file


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

The output is the same weither you choose to store it or not and is presented as follow :


Informations :							       
									 |
@SQ - Reference_sequence_dictionary					 |
Reference sequence name : Reference					 |
Reference sequence length : 1000000					 | 	extraction of header sections 
									 |	and of the corresponding informations
@PG - Program								 |
Program record identifier : bwa						 |
Program name : bwa						         |
									
**********************************				       
total reads count : 351331						 |
	-> properly aligned reads count : 349998 (100% of total reads)	 |
		-> totally mapped reads count : 349891			 |	reads informations
		-> partially mapped reads count : 107			 |
	-> unmapped reads count : 1315					 |
								         |
**********************************			      	
Global cigar mutation observed on well-aligned sequences:	  |
								  |
Raw number of mutations:					  |
Alignment Match : 34999647					  |
Insertion : 18							  |
Deletion : 61							  |
Skipped region : 135						  |
Soft Clipping : 0						  |
Hard Clipping : 0						  |
Padding : 0							  |
Sequence Match : 0						  |	mutations informations
Sequence Mismatch : 0						  | 	deducted from the CIGAR field
								  | 
Percentages:							  |
Alignment Match : 100.0%					  |
Insertion : 0.0%						  |
Deletion : 0.0%							  |
Skipped region : 0.0%						  |
Soft Clipping : 0.0%						  |
Hard Clipping : 0.0%						  |
Padding : 0.0%							  |
Sequence Match : 0.0%						  |
Sequence Mismatch : 0.0%				          |
