"""This script simulates sequence generation using a hidden markov model.
It generates sequences based on a set of states, transition probabilities, and emission probabilities.
The transition probabilities are adjusted dynamically based on the simulation step to increase the 
likelihood of reaching a stop state as the sequence length increases.
"""

import random
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

random.seed(123)

#hidden states
states = ["Start", "Internal", "Stop"]

#initial state probs(k0)
k_0 = {
    "Start": 0.05,
    "Internal": 0.90,
    "Stop": 0.05
}

#base transition matrix(A)
A = {
    "Start":    {"Start": 0.01,  "Internal": 0.98, "Stop": 0.01},
    "Internal": {"Start": 0.1, "Internal": 0.89, "Stop": 0.01},  #will be adjusted dynamically
    "Stop":     {"Start": 0.1, "Internal": 0.89, "Stop": 0.01}
}

#emission probs (E_k)
E_k = {
    "Start": {"ATG": 1.0},
    "Stop": {"TAA": 0.33, "TAG": 0.34, "TGA": 0.33}
}

#internal codons (excluding start and stop codons)
internal_codons = [
    'TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
    'ATT', 'ATC', 'ATA', 'GTT', 'GTC', 'GTA', 'GTG', 'TCT',
    'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'CCT', 'CCC', 'CCA',
    'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA',
    'GCG', 'TAT', 'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT',
    'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
    'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'GGT',
    'GGC', 'GGA', 'GGG', 'TGG'
]

#calculate uniform prob for internal codons
prob_internal = 1.0 / len(internal_codons)
E_k["Internal"] = {codon: prob_internal for codon in internal_codons}

#global variables
num_seqs = 250  #number of sequences to generate   
x = 100 #number of emmisions/codons in sequence
stop_mult = 10 #scaling factor for stop probability increase

def stop_multiplier(t, x):
    """Calculate a multiplier that increases linearly as the sequence lengthens.
    The multiplier is computed as 1 plus the product of a global scaling factor
    (stop_mult) and the ratio of the current step to the max codon length. This
    gradually boosts the probability of transitioning to the stop state over time.
    param: t - current step
    param: x - ref number of emmisions
    return: multiplier value for stop probability
    """
    return 1 + stop_mult * (t / x)

def get_transition_probs(current_state, step, orf_started):
    """Calculate transition probs for the next state given the current state and step.
    For the internal state the stop probability is increased dynamically using the
    stop multiplier and the internal probability is reduced accordingly to keep the overall
    probability normalized. Because the initial prob of stop is never 0, and never 1, all ORF legths up to x-1 
    are possible.
    param: current_state - string current state
    param: step - current step
    return: dict of normalized transition probs for next state
    """
    if current_state == "Internal" and orf_started:
        #when state = internal adjust transition prob for stop based on step
        base = A["Internal"]
        mult = stop_multiplier(step, x)  
        #prob stop adjusted by multiplier
        p_stop = base["Stop"] * mult      
        #reduce internal prob to offset
        p_internal = base["Internal"] - (p_stop - base["Stop"])
        p_start = base["Start"]
        #new transition probs
        return {
            "Start": p_start,
            "Internal": p_internal,
            "Stop": p_stop 
        }
    else:
        #for non-internal states, use fixed transition probs
        return A[current_state]

def write_output(sequences, output_file):
    """Write sequences to a txt in fasta format with orf length.
    param: sequences - tuple sequence string, orf length, sequence index)
    param: output_file - output name
    """
    records = []
    for seq_str, orf_len, i in sequences:
        record = SeqRecord(
            Seq(seq_str),
            id=f"SimSeq_{i+1}",
            description=f"ORF length: {orf_len}"
        )
        records.append(record)
    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta")

def get_first_orf_length(seq):
    """Determine the length of the first ORF in the sequence.
    param: seq - string of nucleotide sequence
    return: orf length in codons 
    """
    start_found = False
    orf_start_idx = 0
    stop_codons = {"TAA", "TAG", "TGA"}
    #iterate over sequence in codon steps 
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if not start_found:
            if codon == "ATG":
                start_found = True 
                orf_start_idx = i 
        else:
            if codon in stop_codons:
                #calculate orf length in codons
                return (i - orf_start_idx) // 3 + 1
    return 0

def choose_initial_state(k_0):
    """Choose an initial state based on the initial probs.
    param: k_0 - dict of initial state probs
    return: string of state
    """
    state_list = list(k_0.keys())
    probs = list(k_0.values())
    return random.choices(state_list, weights=probs, k=1)[0]

def choose_next_state(current_state, step, orf_started):
    """Choose the next state based on current state and simulation step.
    param: current_state - string of current state
    param: step - current step
    return: string next state
    """
    probs_dict = get_transition_probs(current_state, step, orf_started)
    next_states = list(probs_dict.keys())
    probs = list(probs_dict.values())
    return random.choices(next_states, weights=probs, k=1)[0]

def emit_codon(state, E_k):
    """Emit a codon based on the emission probs for the state
    param: state - string current state
    param: E_k - dict of emission probs
    return: codon
    """
    codons = list(E_k[state].keys())
    probs = list(E_k[state].values())
    return random.choices(codons, weights=probs, k=1)[0]

def ORF_generator(k_0, E_k, seq_index):
    """Create ORF based on state, step and emission probs.
    param: k_0 - dict of initial state probs
    param: E_k - dict of emission probs
    param: seq_index - index of sequence
    return: tuple sequence string, orf length, sequence index
    """
    sequence = []
    state = choose_initial_state(k_0)
    #initialize ORF started flag
    orf_started = (state == "Start")
    #simulate sequence codon by codon for a fixed range x
    for step in range(x):
        codon = emit_codon(state, E_k)
        sequence.append(codon)
        #check if an ORF has started
        if not orf_started and codon == "ATG":
            orf_started = True
        state = choose_next_state(state, step, orf_started)
    
    sequence_str = "".join(sequence)
    orf_length = get_first_orf_length(sequence_str)
    return (sequence_str, orf_length, seq_index)

def plot_orf_length_distribution(sequences, plot_file):
    """plot a histogram of ORF lengths.
    param: sequences - tuples of sequence, orf length, and sequence index
    param: plot_file - plot file naem
    return: none
    """
    orf_lengths = [orf_len for _, orf_len, _ in sequences]
    plt.figure(figsize=(12, 7))
    plt.hist(orf_lengths, bins=range(0, max(orf_lengths) + 5, 3), edgecolor='black')
    plt.title("Distribution of First ORF Lengths")
    plt.xlabel("ORF Length (codons)")
    plt.ylabel("Frequency")
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.show()

def main():
    output_file = "aalzaim_ORFoutput.txt"
    plot_file = "orf_length_distribution.png"
    sequences = []

    for i in range(num_seqs):
        seq_tuple = ORF_generator(k_0, E_k, i)
        sequences.append(seq_tuple)

    write_output(sequences, output_file) 
    plot_orf_length_distribution(sequences, plot_file)  

if __name__ == "__main__":
    main()
