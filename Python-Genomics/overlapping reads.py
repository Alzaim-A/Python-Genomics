"""This script is meant to take in reads from a text file (one read per line), create contigs by finding overlapping read, outputting those contigs in a FASTA file and finallay,
create a line graph showing coverage along each contig.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from difflib import SequenceMatcher
from collections import defaultdict
import matplotlib.pyplot as plt


def read_reads(input_file):
    """Reads in the input file and creats a list of reads
    param: input_file - one line per read txt
    return: reads_list - list of reads"""
    reads_list = []
    with open(input_file, 'r') as inputfile:
        for line in inputfile:
            reads_list.append(line.strip())
    return reads_list


def write_contigs(contigs, output_file):
    """Takes contigs and writes them to a FASTA file
    param: contigs - list of whole contigs
    param: output_file - name of FASTA file where contigs will be placed"""
    contig_list = []
    for i, contig in enumerate(contigs, 1):
        record = SeqRecord(
            Seq(contig),
            id=f"contig_{i}",
            description=""
        )
        contig_list.append(record)
    
    with open(output_file, 'w') as fastafile:
        SeqIO.write(contig_list, fastafile, "fasta")


def find_overlap(seq_a, seq_b, k):
    """Finds overlapping sequences and determines match_type
    param: seq_a - current contig(or first read)
    param: seq_b - current read
    return: amount of overlap, match type"""

    #checking for middle match first
    if seq_b in seq_a:
        return (len(seq_b), "middle")    
    
    #set up SequenceMatcher object
    matcher = SequenceMatcher(None, seq_a, seq_b)
    #find_longest_match() returns a named tuple containing (index_a, index_b, overlap)
    match = matcher.find_longest_match(0, len(seq_a), 0, len(seq_b))
    
    #check that overlap >= k and determine forward or reverse match_type
    if match.size >= k:
        if match.a + match.size == len(seq_a) and match.b == 0:
            return (match.size, "forward")
        elif match.b + match.size == len(seq_b) and match.a == 0:
            return (match.size, "reverse")
    return (0, None)


def assemble_reads(reads_list, k):
    """Assembles contigs from reads. Calls find_overlap() to get match_type and amount of overlap. 
    Middle match means the read is fully represented in the contig, duplicate reads would also be handled in this way.
    Forward and reverse matches represent where the read is being added to the contif (front or back).
    param: reads_list - list of reads
    param: k - overlap parameter
    return: contigs - list of assembled whole contigs"""

    reads_list = reads_list.copy()   
    contigs = []
    
    #starts loop with first read as new contig and removes
    while reads_list:
        current_contig = reads_list.pop(0)
        merged = True
        #builds contig with reads based on match_type and removes read
        #will keep looping until no more matches are found
        while merged and reads_list:
            merged = False
            i = 0
            while i < len(reads_list):
                overlap, match_type = find_overlap(current_contig, reads_list[i], k)
                if match_type == "middle":
                    reads_list.pop(i)
                    merged = True   
                elif match_type == "forward":
                    current_contig = current_contig + reads_list[i][overlap:]
                    reads_list.pop(i)
                    merged = True                    
                elif match_type == "reverse":
                    current_contig = reads_list[i] + current_contig[overlap:]
                    reads_list.pop(i)
                    merged = True                    
                else:
                    i += 1
                    
        contigs.append(current_contig)
    
    return contigs


def coverage_plot(reads_list, contigs):
    """Takes reads and maps them back to assembled contigs to create a coverage plot
    param: reads_list -list of reads
    param: contigs - list of assembled whole contigs"""
   
    coverage = defaultdict(list)

    #matches reads, finding start and end position in contig
    for i, contig in enumerate(contigs):
        for read in reads_list:
            if read in contig:
                start = contig.find(read)
                end = start + len(read)
                coverage[f"contig_{i+1}"].append((start, end))
    
    #uses coverage dict to loop through contig and positions
    #adds to line count where read is present
    plt.figure(figsize=(13, 8))
    for contig_name, positions in coverage.items():
        #initlaize coverage_line with max position of mapped read
        coverage_line = [0] * (max([x for _, x in positions]) + 1)
        for start, end in positions:
            for i in range(start, end):
                coverage_line[i] += 1
        plt.plot(coverage_line, label=contig_name)

    plt.xlabel("Position in Contig")
    plt.ylabel("Amount of Coverage")
    plt.title("Contig Coverage")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    input_file = "seqReadFile2023.txt"
    output_file = "Adam_Alzaim.fasta"

    reads_list = read_reads(input_file)
    assembled_contigs = assemble_reads(reads_list, 10)

    write_contigs(assembled_contigs, output_file)
    coverage_plot(reads_list, assembled_contigs)