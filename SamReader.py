#!/usr/bin/env python3
# -*- coding utf-8 -*-

"""Analyzes SAM files and calls substitution mutations.

Usage:
======
    SamReader.py inputSAM

    inputSAM: SAM file(s) to analyze
"""

import os
import re
import sys
import csv
import numpy as np

__authors__ = ("Alizée ARNOUX")
__contact__ = ("alizee.arnoux@etu.umontpellier.fr")
__version__ = "0.0.2"
__date__ = "12/14/2021"
__licence__ = ("This program is free software: you can redistribute it and/or"
               " modify it under the terms of the GNU General Public License"
               " as published by the Free Software Foundation, either version"
               " 3 of the License, or (at your option) any later version."
               " This program is distributed in the hope that it will be"
               " useful, but WITHOUT ANY WARRANTY; without even the implied"
               " warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR"
               " PURPOSE. See the GNU General Public License for more"
               " details. You should have received a copy of the GNU General"
               " Public License along with this program. If not, see"
               " <https://www.gnu.org/licenses/>.")

# Matrix
HEADER_LINE = np.array([['VN', 'SO', 'GO', 'SS'],
                        ['Format version', 'Sorting order of alignments',
                         'Grouping of alignments',
                         'Sub-sorting order of alignments']])
REF_SEQ_DICTIONARY = np.array([['SN', 'LN', 'AH', 'AN', 'AS', 'DS', 'M5',
                                'SP', 'TP', 'UR'],
                               ['Reference sequence name',
                                'Reference sequence length',
                                'Alternate locus',
                                'Alternative reference sequence names',
                                'Genome assembly identifier', 'Description',
                                'MD5 checksum of the sequence', 'Species',
                                'Molecule topology',
                                'URI of the sequence']])
READ_GROUP = np.array([['ID', 'BC', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB',
                        'PG', 'PI', 'PL', 'PM', 'PU', 'SM'],
                       ['Read group identifier', 'Barcode sequence',
                        'Name of sequencing center', 'Description',
                        'Date the run was produced', 'Flow order',
                        'Array of nucleotide bases', 'Library',
                        'Processing programs', 'Predicted median insert size',
                        'Platform/technology', 'Platform model',
                        'Platform unit', 'Sample']])
PROGRAM = np.array([['ID', 'PN', 'CL', 'PP', 'DS', 'VN'],
                    ['Program record identifier', 'Program name',
                     'Command line', 'Previous @PG-ID', 'Description',
                     'program version']])
COMMENTS = np.array([['CO'], ['Commentaire(s)']])
CIGAR_MATRIX = np.array([['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'],
                         ['Alignement Match', 'Insertion', 'Deletion',
                          'Skipped region', 'Soft Clipping', 'Hard Clipping',
                          'Padding', 'Sequence Match', 'Sequence Mismatch']])
QUAL_INTERPRET = np.array([['!', '"', '#', '$', '%', '&', '\'', '(', ')',
                            '*', '+', ', ', '-', '.', '/', '0', '1', '2', '3',
                            '4', '5', '6', '7', '8', '9', ':', ';', '<', '=',
                            '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
                            'H', 'I'],
                           ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                            '10', '11', '12', '13', '14', '15', '16', '17',
                            '18', '19', '20', '21', '22', '23', '24', '25',
                            '26', '27', '28', '29', '30', '31', '32', '33',
                            '34', '35', '36', '37', '38', '39', '40']])

# Dictionary
HEADER_FIELD = {0: HEADER_LINE, 1: REF_SEQ_DICTIONARY,
                2: READ_GROUP, 3: PROGRAM, 4: COMMENTS}
HEADER_TITLE = {0: 'Header line', 1: 'Reference sequence dictionary',
                2: 'Read group', 3: 'Program', 4: 'Comments'}

# list
AMBIGUITY_CODES = ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']

# dictionary
CODONS = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu',
          'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile',
          'ATA': 'Ile', 'ATG': 'Met', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val',
          'GTG': 'Val', 'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
          'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'ACT': 'Thr',
          'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'GCT': 'Ala', 'GCC': 'Ala',
          'GCA': 'Ala', 'GCG': 'Ala', 'TAT': 'Tyr', 'TAC': 'Tyr', 'CAT': 'His',
          'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln', 'AAT': 'Asn', 'AAC': 'Asn',
          'AAA': 'Lys', 'AAG': 'Lys', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu',
          'GAG': 'Glu', 'TGT': 'Cys', 'TGC': 'Cys', 'TGG': 'Trp', 'CGT': 'Arg',
          'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGT': 'Ser', 'AGC': 'Ser',
          'AGA': 'Arg', 'AGG': 'Arg', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly',
          'GGG': 'Gly', 'TAA': 'STOP', 'TAG': 'STOP', 'TGA': 'STOP'}

# Tuple
HEADERS = ('@HD', '@SQ', '@RG', '@PG', '@CO')
OPTIONS = ('-o', '--check')

# Constant
ARGUMENTS_LIST = sys.argv[1:]
MIN_LINE_LENGHT = 11


def help_flag():
    print("\nCombination of bitwise FLAGs. Each bit is explained in the"
          " following table:\n\n"
          "Bit\t\tDescription\n"
          "---------------------------------------------------------"
          "-------------\n"
          "1\t0x1\ttemplate having multiple segments in sequencing\n"
          "2\t0x2\teach segment properly aligned according to the aligner\n"
          "4\t0x4\tsegment unmapped\n"
          "8\t0x8\tnext segment in the template unmapped\n"
          "16\t0x10\tSEQ being reverse complemented\n"
          "32\t0x20\tSEQ of the next segment in the template being"
          " reverse complemented\n"
          "64\t0x40\tthe first segment in the template\n"
          "128\t0x80\tthe last segment in the template\n"
          "256\t0x100\tsecondary alignment\n"
          "512\t0x200\tnot passing filters,  such as platform/vendor"
          " quality controls\n"
          "1024\t0x400\tPCR or optical duplicate\n"
          "2048\t0x800\tsupplementary alignment")
    exit()


def help_req():
    print("\nThis program runs with python 3 and needs the following"
          " packages: \n"
          "os\tsys\tcsv\tre\tnumpy\n\n"
          "You can install them using the following commands: \n"
          "$ conda install <package name>\nor\n$ pip install <package name>")
    exit()


def help_cigar():
    print("\nThe CIGAR operations are given in the following table (set ‘*’"
          " if unavailable): \n\n"
          "Op\tBAM\tDescription\n"
          "-------------------------------------------------------------------"
          "---\n"
          "M\t0\talignment match (can be a sequence match or mismatch)\n"
          "I\t1\tinsertion to the reference\n"
          "D\t2\tdeletion from the reference)\n"
          "N\t3\tskipped region from the reference\n"
          "S\t4\tsoft clipping (clipped sequences present in SEQ)\n"
          "H\t5\thard clipping (clipped sequences NOT present in SEQ)\n"
          "P\t6\tpadding (silent deletion from padded reference)\n"
          "=\t7\tsequence match\n"
          "X\t8\tsequence mismatch")
    exit()


def help_sam():
    print("\nSAM stands for Sequence Alignment/Map format."
          " It is a TAB-delimited text format consisting of a header"
          "section,  which is optional,  and an alignment section."
          " If present,  the header must be prior to the alignments."
          " Header lines start with ‘@’,  while alignment lines do not."
          " Each alignment line has 11 mandatory fields for"
          " essential alignment information such as mapping position,"
          " and variable number of optional fields for flexible"
          " or aligner specific information. Mandatory fields are presented"
          " in the following table :\n\n"
          "Col\tField\tType\tRegex/Range\t\t\tDescription\n"
          "-------------------------------------------------------------------"
          r"---\n"
          r"1\tQNAME\tString\t[!-?A-~]{1, 254}\t\tQuery template NAME\n"
          r"2\tFLAG\tInt\t[0, 2e16 − 1]\t\t\tbitwise FLAG\n"
          r"3\tRNAME\tString\t\*|[: rname: ∧ *=][: rname: ]*\t"
          r"Reference sequence NAME\n"
          r"4\tPOS\tInt\t[0, 2e31 − 1]\t\t\t1-based leftmost mapping"
          " POSition\n"
          r"5\tMAPQ\tInt\t[0, 2e8 − 1]\t\t\tMAPping Quality\n"
          r"6\tCIGAR\tString\t\*|([0-9]+[MIDNSHPX=])+\t\tCIGAR string\n"
          r"7\tRNEXT\tString\t\*|=|[: rname:  *=][: rname: ]*\t"
          r"Reference name of the mate/next read\n"
          r"8\tPNEXT\tInt\t[0, 2e31 − 1]\t\t\tPosition of the mate/next read\n"
          r"9\tTLEN\tInt\t[−2e31 + 1, 2e31 − 1]\t\tobserved Template LENgth\n"
          r"10\tSEQ\tString\t\*|[A-Za-z=.]+\t\t\tsegment SEQuence\n"
          r"11\tQUAL\tString\t[!-~]+\t\t\t\tASCII of Phred-scaled base"
          " QUALity+33")
    exit()


def help_program():
    print("\nSamReader is a small program that analyses SAM files. It takes"
          " one or more SAM file in input and gives an output either in a file"
          " for"
          " each found reference. The syntaxe is as follow :\n"
          "$ /path/to/SamReader.py <input file 1> <input file 2> -o"
          " <output file 1> <output file 2>\n\n"
          "-o followed by the custom name for the output file, if not provided"
          " the default file name is 'outputFile_{n° of file}_{reference}"
          ".txt'\n\nEnter $ /path/to/SamReader.py -h for help")
    exit()


def help():
    """Analyses user input and call the corresponding help function."""
    helpQuery = "NULL"
    # Check if the user is calling the help function.
    if re.search("^-h[a-x]?$", ARGUMENTS_LIST[0]):
        if ARGUMENTS_LIST[0] == "-h":
            print("What do you need help with ?\n"
                  "FLAG field ? Enter -hf \n"
                  "CIGAR field ? Enter -hc\n"
                  "System requirements ? Enter -hr\n"
                  "SAM file ? Enter -hs\n"
                  "SamReader ? Enter -hp\n")
            helpQuery = input("Please enter an option below: \n")
        if ARGUMENTS_LIST[0] == "-hf" or helpQuery == "-hf":
            help_flag()
        if ARGUMENTS_LIST[0] == "-hc" or helpQuery == "-hc":
            help_cigar()
        if ARGUMENTS_LIST[0] == "-hr" or helpQuery == "-hr":
            help_req()
        if ARGUMENTS_LIST[0] == "-hs" or helpQuery == "-hs":
            help_sam()
        if ARGUMENTS_LIST[0] == "-hp" or helpQuery == "-hp":
            help_program()
        # If input is not an option, exit the program.
        print("INPUT ERROR: exit")
        exit()


def input_error_check(argument):
    """Verifies that the inputs are non-empty SAM files."""
    if os.path.isfile(argument) is False:
        print(f"PATH ERROR for '{argument}': "
              "the input is not a file")
        exit()
    elif argument.endswith('.sam') is False:
        print(f"FORMAT ERROR for '{argument}': "
              "only the SAM file format is accepted as input")
        exit()
    elif os.path.getsize(argument) == 0:
        print(f"FILE ERROR for '{argument}': "
              "the file is empty")
        exit()


def input_file_number():
    """Return the number of files given as input."""
    fileNumber = 0
    for argument in ARGUMENTS_LIST:
        if argument in OPTIONS:
            break
        input_error_check(argument)
        fileNumber += 1

    return fileNumber


def integrity_line_number():
    """Returns the number of line to perform integrity-check on."""
    to_check = 0
    if '--check' in ARGUMENTS_LIST:
        i = ARGUMENTS_LIST.index('--check')
        # Finds the instruction given in parameters.
        if ARGUMENTS_LIST[i+1] == 'all':
            to_check = 'ALL'
        else:
            to_check = int(ARGUMENTS_LIST[i+1])

    return to_check


def analysis_progress(line, file_size, octet_analyzed):
    """Shows the progression of file analysis."""
    for field in line:
        # 1 character = 1 octet, + 1 octet for the tabulation.
        octet_analyzed += len(field) + 1
    analysis_progression = round(octet_analyzed * 100 / file_size)
    print(f"progression of file analysis: {analysis_progression}%", end='\r')

    return octet_analyzed


def file_handler(input_file):
    """Reads SAM file as csv file delimited by tabulations."""
    fi = open(ARGUMENTS_LIST[input_file])
    file = csv.reader(fi, delimiter='\t', quoting=csv.QUOTE_NONE)
    fi.close

    return file


def integrity_check(line, re, line_number, to_check):
    """Verifies the integrity of each field for a given number of lines."""
    if type(to_check) == str or line_number <= to_check:
        # There are 11 mandatory fields.
        ERROR_COUNT = 11
        # Verifying regular expression matching to official documentation.
        if re.fullmatch('[!-?A-~]{1, 254}', line[0]):
            ERROR_COUNT -= 1
        # if the field is incorrect,  an error message is displayed
        else:
            print("QNAME ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"QNAME: {line[0]}\nLine: {line}\n")

        if re.fullmatch('[0-9]*', line[1]) and int(line[1]) <= ((2 ** 16) - 1):
            ERROR_COUNT -= 1
        else:
            print("FLAG ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"FLAG: {line[1]}\nLine: {line}\n")

        if re.fullmatch(r'\*|[0-9A-Za-z!#$%&+./: ;?@^_|~-]'
                        '[0-9A-Za-z!#$%&*+./: ;=?@^_|~-]*', line[2]):
            ERROR_COUNT -= 1
        else:
            print("RNAME ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"RNAME: {line[2]}\nLine: {line}\n")

        if re.fullmatch('[0-9]*', line[3]) and int(line[3]) <= ((2 ** 31) - 1):
            ERROR_COUNT -= 1
        else:
            print("POS ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"POS: {line[3]}\nLine: {line}\n")

        if re.fullmatch('[0-9]*', line[4]) \
        and int(line[4]) >= 0 and int(line[4]) <= ((2**8)-1):
            ERROR_COUNT -= 1
        else:
            print("MAPQ ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"MAPQ: {line[4]}\nLine: {line}\n")

        if re.fullmatch(r'\*|([0-9]+[MIDNSHPX=])+', line[5]):
            ERROR_COUNT -= 1
        else:
            print("CIGAR ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"CIGAR: {line[5]}\nLine: {line}\n")

        if re.fullmatch(r'\*|=|[0-9A-Za-z!#$%&+./: ;?@^_|~-]'
                        '[0-9A-Za-z!#$%&*+./: ;=?@^_|~-]*', line[6]):
            ERROR_COUNT -= 1
        else:
            print("RNEXT ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"RNEXT: {line[6]}\nLine: {line}\n")

        if (re.fullmatch('[0-9]*', line[7])
         and int(line[7]) <= ((2 ** 31) - 1)):
            ERROR_COUNT -= 1
        else:
            print("PNEXT ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"PNEXT: {line[7]}\nLine: {line}\n")

        if (re.fullmatch('-[0-9]*|[0-9]*', line[8]) and
        int(line[8]) >= ((-2 ** 31) + 1) and int(line[8]) <= ((2 ** 31) - 1)):
            ERROR_COUNT -= 1
        else:
            print("TLEN ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"TLEN: {line[8]}\nLine: {line}\n")

        if re.fullmatch(r'\*|[A-Za-z=.]+', line[9]):
            ERROR_COUNT -= 1
        else:
            print("SEQ ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"SEQ: {line[9]}\nLine: {line}\n")

        if re.fullmatch('[!-~]+', line[10]):
            ERROR_COUNT -= 1
        else:
            print("QUAL ERROR: "
                  "an unauthorized regular expression was found.\n"
                  f"QUAL: {line[10]}\nLine: {line}\n")

        return ERROR_COUNT

    ERROR_COUNT = 0

    return ERROR_COUNT


def binary_flag(flag):
    """Returns the flag in binary."""
    MAX_FLAG_SIZE = 12
    flag_bin = bin(int(flag))
    # Removes the 'ob' in front of the binary.
    flag_bin = flag_bin[2:]
    flag_bin = list(flag_bin)
    if len(flag_bin) < MAX_FLAG_SIZE:
        difference = MAX_FLAG_SIZE - len(flag_bin)
        for i in range(difference):
            # Insert 0 to complete until the maximal flag size.
            flag_bin.insert(0, '0')

    return flag_bin


def header_analysis(header_count, line, output_head_list):
    """Analyzes the header section."""
    head = 'no'
    if line[0] in HEADERS:
        head = 'yes'
        header_count += 1
        i = HEADERS.index(line[0])
        # Writes the header section title (abbreviated and full).
        output_head_list.append(f"{HEADERS[i]} - {HEADER_TITLE[i]}")
        for subfield in line[1:]:
            # Searches the matrix for the abbreviated subfield title.
            [n, m] = np.where(HEADER_FIELD[i] == subfield[:2])
            # Writes the corresponding full title and file description.
            output_head_list.append(f"{HEADER_FIELD[i][n + 1, m]}: "
                                    f"{subfield[3:]}")
    return head, header_count, output_head_list


def output_header_info(outputFile, output_head_list):
    """Writes header informations."""
    for i in range(len(output_head_list)):
        outputFile.write(f"{output_head_list[i]}\n")
        if output_head_list[i][3:5] in HEADERS:
            outputFile.write("\n")


def paired_reads(flag, unmap_count, badmap_count, totalmap_count, line,
                 dico_pair, not_paired_count, map_status_first,
                 map_status_sec, read_length):
    """Returns True if the read is mapped and stores the mapping status."""
    mapped = True
    totally_mapped = False
    # bin(4) = 100, '4' in flag means unmapped read.
    if int(flag[-3]) == 1:
        unmap_count += 1
        # If the reads are paired, tests if this is the first or second read.
        if int(flag[-7]) == 1:
            map_status_first = "unmapped"
        else:
            map_status_sec = "unmapped"
        mapped = False
    # For the not unmapped reads, counts the totally and badly mapped.
    if mapped is True:
        if re.match('(?![0-9]+M$)', line[5]):
            badmap_count += 1
            if int(flag[-7]) == 1:
                map_status_first = "badly mapped"
            else:
                map_status_sec = "badly mapped"
        else:
            totally_mapped = True
            totalmap_count += 1
            read_length = line[5][:3]
            if int(flag[-7]) == 1:
                map_status_first = "totally mapped"
            else:
                map_status_sec = "totally mapped"
    # If the read is the second in the pair, increments the corresponding dict.
    if int(flag[-8]) == 1:
        if (f'{map_status_sec} + {map_status_first}') in dico_pair.keys():
            dico_pair[f'{map_status_sec} + {map_status_first}'] += 1
        # As the dict.keys are predefined, reverse order of reads if not found.
        else:
            dico_pair[f'{map_status_first} + {map_status_sec}'] += 1
    # If the read is neither the first nor second, then it is not paired.
    if int(flag[-7]) != 1 and int(flag[-8]) != 1:
        not_paired_count += 1

    return (unmap_count, badmap_count, totalmap_count, read_length,
            not_paired_count, dico_pair, map_status_first, map_status_sec,
            totally_mapped)


def dico_init():
    """Creates dictionaries / lists for the current reference."""
    dico_pair = {"unmapped + unmapped": 0,
                  "unmapped + totally mapped": 0, "unmapped + badly mapped": 0,
                  "badly mapped + totally mapped": 0,
                  "badly mapped + badly mapped": 0,
                  "totally mapped + totally mapped": 0}
    dico_cigar = {}
    dico_align = {-1: 0, 0: 0, 1: 0}

    substitutions_list = []
    dico_substitutions = {}
    call_quality_list = []
    mut_position_list = []
    qname_list = []
    aa_orf1_list = []
    aa_orf2_list = []
    aa_orf3_list = []

    return (dico_pair, dico_cigar, dico_align, substitutions_list,
            dico_substitutions, call_quality_list, mut_position_list,
            qname_list, aa_orf1_list, aa_orf2_list, aa_orf3_list)


def fetch_dico(ref):
    """Returns dictionaries for the reference, if already initiated."""
    dico_pair = globals()[f'dico_pair_{ref}']
    dico_cigar = globals()[f'dico_cigar_{ref}']
    dico_align = globals()[f'dico_align_{ref}']
    substitutions_list = globals()[f'sub_list_{ref}']
    dico_substitutions = globals()[f'dico_sub_{ref}']
    call_quality_list = globals()[f'dico_qual_{ref}']
    mut_position_list = globals()[f'dico_pos_{ref}']
    qname_list = globals()[f'dico_qname_{ref}']
    aa_orf1_list = globals()[f'dico_orf1_{ref}']
    aa_orf2_list = globals()[f'dico_orf2_{ref}']
    aa_orf3_list = globals()[f'dico_orf3_{ref}']

    return (dico_align, dico_cigar, dico_pair, substitutions_list,
           dico_substitutions, call_quality_list, mut_position_list,
            qname_list, aa_orf1_list, aa_orf2_list, aa_orf3_list)


def cigar_analysis(line, dico_cigar):
    """Stores cigar info in a dictionary."""
    info_cigar = re.findall(r'[0-9]+\D', line[5])
    for info in info_cigar:
        if info[-1] in dico_cigar.keys():
            dico_cigar[info[-1]] += int(info[:-1])
        else:
            dico_cigar[info[-1]] = int(info[:-1])

    return dico_cigar


def dynamic_dico(dico_pair, dico_cigar, dico_align, substitutions_list,
                 dico_substitutions, ref, call_quality_list, mut_position_list,
                 qname_list, aa_orf1_list, aa_orf2_list, aa_orf3_list):
    """Stores dict data in correponding dynamic dict."""
    globals()[f'dico_pair_{ref}'] = dico_pair
    globals()[f'dico_cigar_{ref}'] = dico_cigar
    globals()[f'dico_align_{ref}'] = dico_align
    globals()[f'sub_list_{ref}'] = substitutions_list
    globals()[f'dico_sub_{ref}'] = dico_substitutions
    globals()[f'dico_qual_{ref}'] = call_quality_list
    globals()[f'dico_pos_{ref}'] = mut_position_list
    globals()[f'dico_qname_{ref}'] = qname_list
    globals()[f'dico_orf1_{ref}'] = aa_orf1_list
    globals()[f'dico_orf2_{ref}'] = aa_orf2_list
    globals()[f'dico_orf3_{ref}'] = aa_orf3_list

    return (globals()[f'dico_pair_{ref}'],
            globals()[f'dico_cigar_{ref}'], globals()[f'dico_align_{ref}'],
            globals()[f'sub_list_{ref}'], globals()[f'dico_sub_{ref}'],
            globals()[f'dico_qual_{ref}'], globals()[f'dico_pos_{ref}'],
            globals()[f'dico_qname_{ref}'], globals()[f'dico_orf1_{ref}'],
            globals()[f'dico_orf2_{ref}'], globals()[f'dico_orf3_{ref}'])


def sub_analysis(line, read_length, qname_list, mut_position_list,
                 substitutions_list, call_quality_list,
                 totally_mapped, aa_orf1_list, aa_orf2_list, aa_orf3_list):
    """Calls substitutions from 'MD:' field and stores related info."""
    x = 0
    for i in range(len(line)):
        if re.search('^MD:Z:', line[i][:5]) and line[i][5:] != read_length:
            x = i
            break
    if totally_mapped is True and x != 0:
        mut_position = 0
        mutations = re.findall(r'[0-9]+\D', line[x][5:])
        for mut in mutations:
            qname_list.append(line[0])
            mut_position += int(mut[:-1])
            # Tests if the mutation is on the complementary read.
            if int(line[8]) < 0:
                mut_position_list.append(int(line[3]) - mut_position)
            else:
                mut_position_list.append(int(line[3]) + mut_position)
            # Stores the mutation in the form 'X -> X'.
            substitutions_list.append(f'{line[9][mut_position]} -> {mut[-1]}')
            # Stores the QUAL value of the substituted base.
            call_quality_list.append(f'{line[10][mut_position]}')

            aa_orf1_list = compares_orf1(mut_position, line, mut, aa_orf1_list)
            aa_orf2_list = compares_orf2(mut_position, line, mut, aa_orf2_list)
            aa_orf3_list = compares_orf3(mut_position, line, mut, aa_orf3_list)

            # Adds 1 to count the analyzed mutation.
            mut_position += 1

    return (qname_list, mut_position_list, substitutions_list,
            aa_orf1_list, aa_orf2_list, aa_orf3_list)


def compares_orf1(mut_position, line, mut, aa_orf1_list):
    """Compares the substituted amino acid on arbitrary ORF1 (as 'XXN')."""
    # Mutated nucleotide must be at least at position 3.
    if mut_position - 2 >= 0:
        # Extract the codon from both sequences.
        ref_orf1 = (line[9][mut_position - 2], line[9][mut_position - 1],
                    line[9][mut_position])
        query_orf1 = (line[9][mut_position - 2], line[9][mut_position - 1],
                      mut[-1])
        # Associates the codon to the corresponding amino acid.
        ref_orf1 = "".join(ref_orf1)
        query_orf1 = "".join(query_orf1)
        if ref_orf1 in CODONS:
            ref_orf1 = CODONS[ref_orf1]
        else:
            ref_orf1 = 'ND'
        if query_orf1 in CODONS:
            query_orf1 = CODONS[query_orf1]
        else:
            query_orf1 = 'ND'
        if ref_orf1 == query_orf1:
            if ref_orf1 != 'ND':
                aa_orf1_list.append('synonymous')
            else:
                aa_orf1_list.append('ND')
        else:
            aa_orf1_list.append(f'{ref_orf1} to {query_orf1}')
    # If nucleotide is at position 2 or less, amino acid is "not determined".
    else:
        aa_orf1_list.append('ND')

    return(aa_orf1_list)


def compares_orf2(mut_position, line, mut, aa_orf2_list):
    """Compares the substituted amino acid on arbitrary ORF2 (as 'XNX')."""
    # Mutated nucleotide must be at least at position 2 and not the last one.
    if mut_position - 1 >= 0 and mut_position + 1 < len(line[9]):
        ref_orf2 = (line[9][mut_position - 1], line[9][mut_position],
                    line[9][mut_position + 1])
        query_orf2 = (line[9][mut_position - 1], mut[-1],
                      line[9][mut_position + 1])
        # Associates the codon to the corresponding amino acid.
        ref_orf2 = "".join(ref_orf2)
        query_orf2 = "".join(query_orf2)
        if ref_orf2 in CODONS:
            ref_orf2 = CODONS[ref_orf2]
        else:
            ref_orf2 = 'ND'
        if query_orf2 in CODONS:
            query_orf2 = CODONS[query_orf2]
        else:
            query_orf2 = 'ND'
        if ref_orf2 == query_orf2:
            if ref_orf2 != 'ND':
                aa_orf2_list.append('synonymous')
            else:
                aa_orf2_list.append('ND')
        else:
            aa_orf2_list.append(f'{ref_orf2} to {query_orf2}')
    # If nucleotide is at position 2 or less, amino acid is "not determined".
    else:
        aa_orf2_list.append('ND')

    return(aa_orf2_list)


def compares_orf3(mut_position, line, mut, aa_orf3_list):
    """Compares the substituted amino acid on arbitrary ORF3 (as 'NXX')."""
    # Mutated nucleotide must not be the last or second last one.
    if (mut_position + 2) < len(line[9]):
        ref_orf3 = (line[9][mut_position], line[9][mut_position + 1],
                    line[9][mut_position + 2])
        query_orf3 = (mut[-1], line[9][mut_position + 1],
                      line[9][mut_position + 2])
        # Associates the codon to the corresponding amino acid.
        ref_orf3 = "".join(ref_orf3)
        query_orf3 = "".join(query_orf3)
        if ref_orf3 in CODONS:
            ref_orf3 = CODONS[ref_orf3]
        else:
            ref_orf3 = 'ND'
        if query_orf3 in CODONS:
            query_orf3 = CODONS[query_orf3]
        else:
            query_orf3 = 'ND'
        if ref_orf3 == query_orf3:
            if ref_orf3 != 'ND':
                aa_orf3_list.append('synonymous')
            else:
                aa_orf3_list.append('ND')
        else:
            aa_orf3_list.append(f'{ref_orf3} to {query_orf3}')
    # If nucleotide is at position 2 or less, amino acid is "not determined".
    else:
        aa_orf3_list.append('ND')

    return(aa_orf3_list)


def alignement_pairs(line, dico_align, read_length):
    """Returns the size of difference or overlap between the paired reads."""
    # Analyzes only one of the paired reads.
    if int(line[8]) > 0:
        # Determines if there is a difference or overlap.
        difference = int(line[8]) - (2*int(read_length))
        if difference == 0:
            dico_align[0] += 1
        # If there is a gap.
        elif difference < 0:
            gap = difference
            # Rounds the gap to the nearest lower ten.
            for i in range(-int(read_length), 0, 10):
                if i - gap >= 0:
                    if i in dico_align.keys():
                        dico_align[i] += 1
                    else:
                        dico_align[i] = 1
                    break
                # If the gap is 10 nucleotide or less add to a predefined key.
                elif gap > -10:
                    dico_align[-1] += 1
                    break
        # If there is an overlap.
        else:
            overlap = difference
            # Rounds the overlap to the nearest lower ten.
            for i in range(0, int(read_length), 10):
                if i - overlap >= 0:
                    if i in dico_align.keys():
                        dico_align[i] += 1
                    else:
                        dico_align[i] = 1
                    break
                # If overlap is 10 nucleotide or less add to a predefined key.
                elif overlap < 10:
                    dico_align[1] += 1
                    break

    return dico_align


def error_input(ERROR_COUNT, line_number):
    """Determines what to do when errors are found, based on user input."""
    error_search = input(f"Document non analysable: {ERROR_COUNT} erreur(s)"
                         " d'expressions régulières trouvée(s) à la line"
                         f" {line_number}, souhaitez-vous rechercher les"
                         " erreurs dans le reste du fichier ?\ny/n\n")
    if error_search == "n":
        exit()
    elif error_search == 'y':
        to_check = 'NULL'
    else:
        print("INPUT ERROR: exit.")
        exit()

    return error_search, to_check


def cigar_total_count(dico_cigar):
    """"Calculates the total of all cigar value for each reference."""
    cigar_total = 0
    for value in dico_cigar:
        cigar_total += dico_cigar[value]

    return cigar_total


def paired_total_count(dico_pair):
    """Calculates the total pairs number for each reference."""
    paired_total = 0
    for value in dico_pair:
        paired_total += dico_pair[value]

    return paired_total


def substitution_count(substitutions_list, dico_substitutions):
    """Counts each possible substitutions found (A->T,  A->C,  A->G,  etc)."""
    for mutation in substitutions_list:
        if mutation in dico_substitutions.keys():
            dico_substitutions[mutation] += 1
        else:
            dico_substitutions[mutation] = 1
    # Sorts the dictionary in descending order.
    sorted_dico_sub = sorted(dico_substitutions.items(), key=lambda x: x[1],
                             reverse=True)

    return sorted_dico_sub


def output_read_info(outputFile, header_count, totalmap_count, badmap_count,
                     unmap_count, line_number):
    """Writes info regarding reads."""
    percent = round((totalmap_count + badmap_count) * 100 /
                    (line_number - header_count), 2)
    outputFile.write(f"\n\ntotal read count: {line_number - header_count}"
                     "\n\t-> aligned read count:"
                     f" {totalmap_count + badmap_count} ({percent}%"
                     " of total read)"
                     f"\n\t\t-> totally mapped read count: {totalmap_count}"
                     f"\n\t\t-> badly mapped read count: {badmap_count}"
                     f"\n\t-> unmapped read count: {unmap_count}")


def output_pairs_info(outputFile, dico_pair, not_paired_count,
                      paired_total, ref):
    """Writes info about read pairs."""
    outputFile.write("\n\n-> Analysis of paired-reads:")
    for x in dico_pair:
        if dico_pair[x] != 0:
            outputFile.write(f"\n{x}: "
                             f"{round(dico_pair[x] * 100 / (paired_total), 4)}"
                             "%"
                             f" ({dico_pair[x]} out of {paired_total} pairs)")
    if not_paired_count != 0:
        outputFile.write(f"-> {not_paired_count} reads are not paired.")


def output_cigar(outputFile, CIGAR_MATRIX, dico_cigar, cigar_total):
    """Writes info about cigar mutation."""
    outputFile.write("\n\n-> Global CIGAR mutations observed on aligned"
                     " sequences:\n\n")
    for key in dico_cigar.keys():
        [n, m] = np.where(CIGAR_MATRIX == key)
        outputFile.write(f"{CIGAR_MATRIX[n + 1, m]}:"
                         f" {round((dico_cigar[key] * 100 / cigar_total), 4)}%"
                         f" ({dico_cigar[key]} out of {cigar_total}"
                         " nucleotides)\n")


def output_sub(outputFile, sorted_dico_sub):
    """Writes substitutions information."""
    outputFile.write("\n\n")
    # Tests whether the dictionary is not empty.
    if sorted_dico_sub:
        outputFile.write("-> Summary of nucleotide substitutions"
                         " :\n\nSubstitution\t\tIteration"
                         "\n------------------------------------\n")
        for x in range(len(sorted_dico_sub)):
            outputFile.write(f"{sorted_dico_sub[x][0]}\t\t\t"
                             f"{sorted_dico_sub[x][1]}\n")
    else:
        outputFile.write("-> No substitutions were found.")


def rename_file(input_file, ref):
    """Renames the output file if a name is given by the user."""
    if '-o' in ARGUMENTS_LIST:
        outputIndex = ARGUMENTS_LIST.index('-o')
        os.rename(f'outputFile_{input_file}_{ref}.txt',
                  f'{ARGUMENTS_LIST[outputIndex + input_file]}_{ref}.txt')
        print(f"OUTPUT FILE: {ARGUMENTS_LIST[outputIndex + input_file]}_{ref}"
              ".txt created.")
    else:
        print(f"OUTPUT FILE: outputFile_{input_file}_{ref}.txt created.")


def csv_sub_writes(input_file, substitutions_list, call_quality_list,
                   qname_list, mut_position_list, aa_orf1_list, aa_orf2_list,
                   aa_orf3_list, ref):
    """Compile in a csv file info about substitutions."""
    outCsv = open(f'mutationFile_{input_file}_{ref}.csv', 'w')
    outCsv.write(f"{ARGUMENTS_LIST[input_file]}\nn°read, Position,"
                 "Mutation, Base call accuracy (%), ORF1, ORF2, ORF3\n")
    for x in range(len(substitutions_list)):
        for call in range(len(QUAL_INTERPRET[0, ])):
            if call_quality_list[x] == QUAL_INTERPRET[0, call]:
                quality = int(QUAL_INTERPRET[1, call])
                outCsv.write(f"{qname_list[x]}, {mut_position_list[x]},"
                             f" {substitutions_list[x]}, "
                             f"{round((1-(10**(-quality/10)))*100, 2)}"
                             f", {aa_orf1_list[x]}, {aa_orf2_list[x]},"
                             f" {aa_orf3_list[x]}\n")
    outCsv.close()
    print(f"\nCSV FILE:  mutationFile_{input_file}_{ref}.csv created.")


def output_align_reads(dico_align, outputFile):
    """Writes the number of well aligned read pairs."""
    if dico_align[0] != 0:
        outputFile.write(f"{dico_align[0]} read pairs are well aligned"
                         " on both ends of the fragment.\n")


def output_gap_reads(dico_align, outputFile):
    """Writes the number of paired-reads with a gap between them."""
    if dico_align[1] != 0:
        outputFile.write(f"{dico_align[1]} read pair(s) present a"
                         " gap of [1,10[ nucleotides between the forward"
                         " and the reverse read.\n")
    for gap in dico_align.keys():
        if gap > 0 and gap != 1:
            outputFile.write(f"{dico_align[gap]} read pair(s) present"
                             f" a gap of [{gap},{gap + 10}[ nucleotides"
                             " between the forward and the reverse read.\n")


def output_overlap_reads(dico_align, outputFile):
    """Writes the number of paired-read with an overlap between them."""
    if dico_align[-1] != 0:
        outputFile.write(f"{dico_align[-1]} read pair(s) present an"
                         " overlap of [1,10[ nucleotides between the"
                         " forward and the reverse read.\n")
    for overlap in dico_align.keys():
        if overlap < 0 and overlap != -1:
            outputFile.write(f"{dico_align[overlap]} read pair(s) present"
                             f" an overlap of [{-overlap},{-overlap + 10}["
                             " nucleotides between the forward and the reverse"
                             " read.\n")


def main():
    if len(sys.argv) == 1:
        help_program()
    help()
    fileNumber = input_file_number()
    map_status_first = 'NULL'
    map_status_sec = 'NULL'
    for input_file in range(fileNumber):
        current_file = ARGUMENTS_LIST[input_file]
        file_size = os.path.getsize(current_file)
        file = file_handler(input_file)
        to_check = integrity_line_number()
        # variables to reset between each file
        line_number = 0
        octet_analyzed = 0
        header_count = 0
        read_length = 0
        unmap_count = 0
        badmap_count = 0
        totalmap_count = 0
        paired_total = 0
        not_paired_count = 0
        research_query = False
        ref = 'NULL'
        error_search = 'NULL'
        references = []
        output_head_list = []

        print(f"\nAnalyzing:\n{current_file}\n")

        for line in file:
            line_number += 1

            (head, header_count,
             output_head_list) = header_analysis(header_count, line,
                                                 output_head_list)

            if head == 'no':
                ERROR_COUNT = integrity_check(line, re, line_number,
                                              to_check)

            # If there is no errors on the line, ERROR_COUNT equals 0.
            if head == 'no' and ERROR_COUNT == 0 and research_query is False:
                if line[2] != ref:
                    ref = line[2]
                    if ref not in references:
                        references.append(ref)
                        (dico_pair, dico_cigar, dico_align, substitutions_list,
                         dico_substitutions, call_quality_list,
                         mut_position_list, qname_list, aa_orf1_list,
                         aa_orf2_list, aa_orf3_list) = dico_init()
                    else:
                        (dico_align, dico_cigar, dico_pair, substitutions_list,
                         dico_substitutions, call_quality_list,
                         mut_position_list, qname_list, aa_orf1_list,
                         aa_orf2_list, aa_orf3_list) = fetch_dico(ref)

                octet_analyzed = analysis_progress(line, file_size,
                                                   octet_analyzed)

                flag = binary_flag(line[1])

                (unmap_count, badmap_count, totalmap_count, read_length,
                 not_paired_count, dico_pair, map_status_first,
                 map_status_sec,
                 totally_mapped) = paired_reads(flag,
                                                unmap_count, badmap_count,
                                                totalmap_count, line,
                                                dico_pair, not_paired_count,
                                                map_status_first,
                                                map_status_sec, read_length)

                dico_cigar = cigar_analysis(line, dico_cigar)

                (qname_list, mut_position_list,
                 substitutions_list, aa_orf1_list, aa_orf2_list,
                 aa_orf3_list) = sub_analysis(line, read_length, qname_list,
                                              mut_position_list,
                                              substitutions_list,
                                              call_quality_list,
                                              totally_mapped,
                                              aa_orf1_list, aa_orf2_list,
                                              aa_orf3_list)

                dico_align = alignement_pairs(line, dico_align,
                                              read_length)

                (globals()[f'dico_pair_{ref}'], globals()[f'dico_cigar_{ref}'],
                 globals()[f'dico_align_{ref}'], globals()[f'sub_list_{ref}'],
                 globals()[f'dico_sub_{ref}'],globals()[f'dico_qual_{ref}'],
                 globals()[f'dico_pos_{ref}'], globals()[f'dico_qname_{ref}'],
                 globals()[f'dico_orf1_{ref}'], globals()[f'dico_orf2_{ref}'],
                 globals()[f'dico_orf3_{ref}']) = dynamic_dico(dico_pair,
                                                    dico_cigar, dico_align,
                                                    substitutions_list,
                                                    dico_substitutions, ref,
                                                    call_quality_list,
                                                    mut_position_list,
                                                    qname_list, aa_orf1_list,
                                                    aa_orf2_list, aa_orf3_list)

            # If there is error on the line, passe once in the condition.
            if head == 'no' and ERROR_COUNT != 0 and research_query is False:
                research_query = True
                error_search, to_check = error_input(ERROR_COUNT,
                                                     line_number)

        if error_search == "y":
            print("End of error research.")
            exit()

        if totalmap_count + badmap_count == 0:
            print("No reads could be analyzed")
            exit()

        sys.stdout.write("\033[F")  # Cursor up one line.
        sys.stdout.write("\033[K")  # Clear the entire line.

        for ref in references:
            (dico_align, dico_cigar, dico_pair, substitutions_list,
             dico_substitutions, call_quality_list, mut_position_list,
             qname_list, aa_orf1_list, aa_orf2_list,
             aa_orf3_list) = fetch_dico(ref)

            sorted_dico_sub = substitution_count(substitutions_list,
                                                 dico_substitutions)

            cigar_total = cigar_total_count(dico_cigar)
            paired_total = paired_total_count(dico_pair)

            csv_sub_writes(input_file, substitutions_list, call_quality_list,
                           qname_list, mut_position_list, aa_orf1_list,
                           aa_orf2_list, aa_orf3_list, ref)

            with open(f'outputFile_{input_file}_{ref}.txt', 'w') as outputFile:
                outputFile.write(f"{current_file}\n\nFile informations:\n\n\n")
                output_header_info(outputFile, output_head_list)
                output_read_info(outputFile, header_count, totalmap_count,
                                 badmap_count, unmap_count, line_number)

                outputFile.write("\n\n**********************************\n"
                                 f"\nInformations relative to reference: {ref}"
                                 "\n")
                output_pairs_info(outputFile, dico_pair, not_paired_count,
                                  paired_total, ref)
                if paired_total != 0:
                    outputFile.write("\n\n-> Pairs alignement analysis:\n")
                    output_align_reads(dico_align, outputFile)
                    output_gap_reads(dico_align, outputFile)
                    output_overlap_reads(dico_align, outputFile)
                output_cigar(outputFile, CIGAR_MATRIX, dico_cigar,
                             cigar_total)
                output_sub(outputFile, sorted_dico_sub)

            rename_file(input_file, ref)


if __name__ == "__main__":
    main()
