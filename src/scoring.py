from blosum62 import get_blosum62_score

def sp_score_dna(aligned_sequences, match=1, mismatch=-1, gap_open=-2, gap_extend=-1):
    '''Calculates SP score for DNA sequences.'''
    return _sp_score(aligned_sequences, seq_type='DNA', match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)

def sp_score_protein(aligned_sequences, gap_open=-12, gap_extend=-1):
    '''Calculates SP score for Protein sequences.'''
    return _sp_score(aligned_sequences, seq_type='PROTEIN', match=None, mismatch=None, gap_open=gap_open, gap_extend=gap_extend)

def _sp_score(aligned_sequences, seq_type, match, mismatch, gap_open, gap_extend):
    num_seqs = len(aligned_sequences)
    if num_seqs < 2:
        return 0
    total_score = 0
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            total_score += _pair_score(aligned_sequences[i], aligned_sequences[j], seq_type, match, mismatch, gap_open, gap_extend)
    return total_score

def _pair_score(seq1, seq2, seq_type, match, mismatch, gap_open, gap_extend):
    score = 0
    gap_state1 = False
    gap_state2 = False
    
    for c1, c2 in zip(seq1, seq2):
        if c1 == '-' and c2 == '-':
            gap_state1 = False
            gap_state2 = False
            continue
            
        is_gap1 = (c1 == '-')
        is_gap2 = (c2 == '-')
        
        if is_gap1:
            if not gap_state1:
                score += gap_open
            else:
                score += gap_extend
            gap_state1 = True
        else:
            gap_state1 = False
            
        if is_gap2:
            if not gap_state2:
                score += gap_open
            else:
                score += gap_extend
            gap_state2 = True
        else:
            gap_state2 = False
            
        if not is_gap1 and not is_gap2:
            if seq_type == 'PROTEIN':
                score += get_blosum62_score(c1, c2)
            else:
                if c1 == c2:
                    score += match
                else:
                    score += mismatch
    return score

def print_alignment(sequences):
    for name, seq in sequences.items():
        print(f"{name: <10} {seq}")
