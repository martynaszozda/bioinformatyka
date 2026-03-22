import random
import os
from fasta_parser import save_fasta

def generate_base_sequence(length, seq_type):
    if seq_type == 'DNA':
        alphabet = ['A', 'C', 'G', 'T']
    else:
        # Standard 20 amino acids
        alphabet = list("ARNDCQEGHILKMFPSTWYV")
    return ''.join(random.choices(alphabet, k=length))

def mutate_sequence(seq, seq_type, mutation_rate=0.1, indel_rate=0.05):
    if seq_type == 'DNA':
        alphabet = ['A', 'C', 'G', 'T']
    else:
        alphabet = list("ARNDCQEGHILKMFPSTWYV")
        
    new_seq = []
    for char in seq:
        r = random.random()
        if r < indel_rate / 2:
            # Deletion (skip character)
            pass
        elif r < indel_rate:
            # Insertion
            new_seq.append(char)
            new_seq.append(random.choice(alphabet))
        elif r < indel_rate + mutation_rate:
            # Substitution
            choices = [c for c in alphabet if c != char]
            new_seq.append(random.choice(choices))
        else:
            new_seq.append(char)
            
    return ''.join(new_seq)

def create_dataset(num_seqs, base_length, seq_type, similarity, output_file):
    if similarity == 'high':
        mut_rate, indel_rate = 0.05, 0.02
    elif similarity == 'medium':
        mut_rate, indel_rate = 0.15, 0.05
    else: # low
        mut_rate, indel_rate = 0.30, 0.10
        
    base_seq = generate_base_sequence(base_length, seq_type)
    sequences = {}
    for i in range(num_seqs):
        mutated = mutate_sequence(base_seq, seq_type, mut_rate, indel_rate)
        sequences[f"seq_{i+1}"] = mutated
        
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    save_fasta(sequences, output_file)

def generate_all_datasets(base_dir="../datasets"):
    # DNA Datasets
    configs = [
        ('small', 'DNA', 5, 30),
        ('medium', 'DNA', 7, 60),
        ('large', 'DNA', 9, 100),
        ('small', 'PROTEIN', 5, 40),
        ('medium', 'PROTEIN', 7, 70),
    ]
    
    similarities = ['high', 'medium', 'low']
    
    for size, stype, n_seqs, length in configs:
        for sim in similarities:
            for i in range(3): # Generate 3 instances per config
                filename = f"{base_dir}/{stype.lower()}_{size}_{sim}_{i+1}.fasta"
                create_dataset(n_seqs, length, stype, sim, filename)
                print(f"Generated {filename}")

if __name__ == '__main__':
    generate_all_datasets()
