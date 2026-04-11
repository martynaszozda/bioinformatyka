import os
import glob
import time
import csv
from fasta_parser import read_fasta
from aco import ACO_MSA

def run_experiment(fasta_file, seq_type, out_csv, num_ants, iterations, alpha, beta, rho):
    seqs = read_fasta(fasta_file)
    if not seqs:
        return
        
    start_time = time.time()
    
    if seq_type == 'DNA':
        gap_open, gap_extend = -2, -1
        match, mismatch = 1, -1
    else: # PROTEIN
        gap_open, gap_extend = -12, -1
        match, mismatch = None, None
        
    aco = ACO_MSA(seqs, seq_type=seq_type, num_ants=num_ants, iterations=iterations,
                  alpha=alpha, beta=beta, rho=rho, 
                  match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)
                  
    score, alignment = aco.run()
    elapsed = time.time() - start_time
    
    file_exists = os.path.isfile(out_csv)
    with open(out_csv, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['Dataset', 'SeqType', 'NumAnts', 'Iterations', 'Alpha', 'Beta', 'Rho', 'Score', 'Time(s)'])
        dataset_name = os.path.basename(fasta_file)
        writer.writerow([dataset_name, seq_type, num_ants, iterations, alpha, beta, rho, score, f"{elapsed:.2f}"])
    
    print(f"Finished {dataset_name}: Score {score} in {elapsed:.2f}s")
    return score

def run_all_sweeps():
    base_path = os.path.dirname(os.path.abspath(__file__))
    datasets_dir = os.path.join(base_path, "..", "datasets")
    results_file = os.path.join(base_path, "..", "results.csv")
    
    with open(results_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Dataset', 'SeqType', 'NumAnts', 'Iterations', 'Alpha', 'Beta', 'Rho', 'Score', 'Time(s)'])

    files = glob.glob(f"{datasets_dir}/*.fasta")
    if not files:
        print(f"No datasets found in {datasets_dir}")
        return
        
    # We will test a parameter sweep on one small dataset to demonstrate
    test_files = [f for f in files if 'small' in f and 'medium' in f]
    if not test_files:
        test_files = files[:2]

    num_repeats = 10

    for f_path in test_files:
        f_name = os.path.basename(f_path)
        seq_type = 'DNA' if 'dna' in f_name.lower() else 'PROTEIN'
        
        for ants in [5, 10]:
            for r in range(num_repeats):
                run_experiment(f_path, seq_type, results_file, 
                               num_ants=ants, iterations=20, 
                               alpha=1.0, beta=1.0, rho=0.1)
                print(end="", flush=True)

        for r_val in [0.05, 0.2]:
            for r in range(num_repeats):
                run_experiment(f_path, seq_type, results_file, 
                               num_ants=10, iterations=20, 
                               alpha=1.0, beta=1.0, rho=r_val)
                print(end="", flush=True)

if __name__ == '__main__':
    run_all_sweeps()
