#!/usr/bin/env python3

import re
from datetime import datetime


def main():
    genome = input("Please type in a DNA sequence: ").lstrip().rstrip().upper()
    motive = input("Please type in a DNA string: ").lstrip().rstrip().upper()

    filename1 = "dna-sequencer.fasta"
    filename2 = "dna-sequencer.txt"
    timestamp = datetime.now()

    if dna_sequence_validator(genome, motive) == True:
        valid = "Valid input of DNA sequence and DNA string!"
        print(valid)

        function1 = dna_sequence_length(genome)
        function2 = calculate_gc_content(genome)
        function3 = reverse_complement(genome)
        function4 = translate_dna_to_protein(genome)
        function5 = find_motive(genome, motive)

        print(function1)
        print(function2)
        print(function3)
        print(function4)
        print(function5)

        with open(filename2, "a") as file:
            file.write(f"{timestamp} - Status: {valid}\n")
            file.write(f"Your DNA sequence: {genome}\n")
            file.write(f"Your DNA string: {motive}\n")
            file.write(f"{function1}\n")
            file.write(f"{function2}\n")
            file.write(f"{function3}\n")
            file.write(f"{function4}\n")
            file.write(f"{function5}\n")
            file.write("\n")

        print("Textfile with your specific output saved - see program folder!")

        exit()

    else:
        invalid = "Invalid input of DNA sequence or DNA string!"
        print(invalid)

        with open(filename2, "a") as file:
            file.write(f"{timestamp} - Status: {invalid}\n")
            file.write(f"Your DNA sequence: {genome}\n")
            file.write(f"Your DNA string: {motive}\n")
            file.write("\n")

        exit()

    # Sample input 1: ATGCGTACGTAGCTAGCTAGCTGATCGTAGCTAGCGTACGTAGCTAGCTGACTGACTGATCGTAGCTGACTGATCGTAGCGTACGTGCTAGCGTCTAGCTGATCGTAGCTAGCTGATCGTGCTAGCTAGCTAGCTGACTGATCGTACGTGATCGTAGCTGATCGTAGCTAGCTGACTGACGTACGTAGA
    # Sample input 2: GAT / GTA / TAGC / CTGAC / GTACG / etc. (choose one of inbetween "/")


def dna_sequence_validator(genome: str, motive: str):
    pattern = r"^[ATGC]+$"

    if re.match(pattern, genome) and re.match(pattern, motive):
        return True

    elif genome == "" or motive == "":
        return False

    else:
        return False


def dna_sequence_length(genome: str):
    return f"Length of your DNA sequence: {len(genome)}"


def calculate_gc_content(genome: str):
    gccount = genome.count("G") + genome.count("C")
    gccontent = (gccount / len(genome)) * 100

    return f"GC content of your DNA sequence: {gccontent:.2f}%"


def reverse_complement(genome: str):
    comp = ""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    reversecomplement = comp.join([complement[base] for base in genome[::-1]])

    return f"Reverse complement of your DNA sequence: {reversecomplement}"


def translate_dna_to_protein(genome: str):
    codontable = {
        "TCA": "S",  # Serina
        "TCC": "S",  # Serina
        "TCG": "S",  # Serina
        "TCT": "S",  # Serina
        "TTC": "F",  # Fenilalanina
        "TTT": "F",  # Fenilalanina
        "TTA": "L",  # Leucina
        "TTG": "L",  # Leucina
        "TAC": "Y",  # Tirosina
        "TAT": "Y",  # Tirosina
        "TAA": "*",  # Stop
        "TAG": "*",  # Stop
        "TGC": "C",  # Cisteina
        "TGT": "C",  # Cisteina
        "TGA": "*",  # Stop
        "TGG": "W",  # Triptofano
        "CTA": "L",  # Leucina
        "CTC": "L",  # Leucina
        "CTG": "L",  # Leucina
        "CTT": "L",  # Leucina
        "CCA": "P",  # Prolina
        "CCC": "P",  # Prolina
        "CCG": "P",  # Prolina
        "CCT": "P",  # Prolina
        "CAC": "H",  # Histidina
        "CAT": "H",  # Histidina
        "CAA": "Q",  # Glutamina
        "CAG": "Q",  # Glutamina
        "CGA": "R",  # Arginina
        "CGC": "R",  # Arginina
        "CGG": "R",  # Arginina
        "CGT": "R",  # Arginina
        "ATA": "I",  # Isoleucina
        "ATC": "I",  # Isoleucina
        "ATT": "I",  # Isoleucina
        "ATG": "M",  # Methionina
        "ACA": "T",  # Treonina
        "ACC": "T",  # Treonina
        "ACG": "T",  # Treonina
        "ACT": "T",  # Treonina
        "AAC": "N",  # Asparagina
        "AAT": "N",  # Asparagina
        "AAA": "K",  # Lisina
        "AAG": "K",  # Lisina
        "AGC": "S",  # Serina
        "AGT": "S",  # Serina
        "AGA": "R",  # Arginina
        "AGG": "R",  # Arginina
        "GTA": "V",  # Valina
        "GTC": "V",  # Valina
        "GTG": "V",  # Valina
        "GTT": "V",  # Valina
        "GCA": "A",  # Alanina
        "GCC": "A",  # Alanina
        "GCG": "A",  # Alanina
        "GCT": "A",  # Alanina
        "GAC": "D",  # Acido Aspartico
        "GAT": "D",  # Acido Aspartico
        "GAA": "E",  # Acido Glutamico
        "GAG": "E",  # Acido Glutamico
        "GGA": "G",  # Glicina
        "GGC": "G",  # Glicina
        "GGG": "G",  # Glicina
        "GGT": "G",  # Glicina
    }

    # Source of the Codon Table dictionary: https://gist.github.com/juanfal/09d7fb53bd367742127e17284b9c47bf

    protein = ""

    for x in range(0, len(genome), 3):
        codon = genome[x : x + 3]
        protein += codontable.get(codon, "?")

    return f"Protein sequence of your DNA sequence: {protein}"


def find_motive(genome: str, motive: str):
    index = []

    for sub in re.finditer(motive, genome):
        indices = sub.start()
        index.append(indices)

    if index == []:
        return f"The indices of your DNA string inside the DNA sequence: {None}"

    elif index != []:
        return f"The indices of your DNA string inside the DNA sequence: {index}"


if __name__ == "__main__":
    main()
