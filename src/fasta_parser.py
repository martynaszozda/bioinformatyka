def read_fasta(file_path):
    '''Reads FASTA files into a dictionary.'''
    sequences = {}
    current_name = None
    current_seq = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
            
    return sequences

def save_fasta(sequences, file_path):
    '''Saves a dictionary of sequences to a FASTA file.'''
    with open(file_path, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            f.write(f"{seq}\n")
