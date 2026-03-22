from aco import ACO_MSA
from scoring import print_alignment

def test_aco():
    seqs = {
        's1': 'AGTCGTAG',
        's2': 'ATCGTCG',
        's3': 'GTAG',
        's4': 'GTAGAG'
    }
    aco = ACO_MSA(seqs, seq_type='DNA', num_ants=5, iterations=5, 
                  match=1, mismatch=-1, gap_open=-2, gap_extend=-1)
    score, alignment = aco.run()
    
    print(f"Final best score: {score}")
    print_alignment({f's{i+1}': alignment[i] for i in range(len(alignment))})

if __name__ == '__main__':
    test_aco()
