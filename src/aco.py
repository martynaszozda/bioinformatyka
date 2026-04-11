import math
import random
import itertools
from blosum62 import get_blosum62_score
from scoring import _sp_score

class ACO_MSA:
    def __init__(self, sequences, seq_type='PROTEIN', 
                 num_ants=10, iterations=50, 
                 alpha=1.0, beta=1.0, rho=0.1, 
                 match=1, mismatch=-1, gap_open=-12, gap_extend=-1):
        self.sequences = sequences
        self.seq_names = list(sequences.keys())
        self.seqs = [sequences[name] for name in self.seq_names]
        self.k = len(self.seqs)
        self.seq_lengths = [len(seq) for seq in self.seqs]
        
        self.seq_type = seq_type
        
        # ACO params
        self.num_ants = num_ants
        self.iterations = iterations
        self.alpha = alpha
        self.beta = beta
        self.rho = rho # evaporation rate
        
        # Scoring params
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
        self.pheromones = {} # key: (state, move), value: tau
        self.tau_0 = 1.0
        
        # All 2^k - 1 possible moves
        self.all_moves = list(itertools.product([0, 1], repeat=self.k))
        self.all_moves.remove(tuple([0]*self.k)) # Remove all-zero move
        
        self.best_alignment = None
        self.best_score = -float('inf')

    def run(self):
        for it in range(self.iterations):
            ant_alignments = []
            ant_paths = []
            
            for ant in range(self.num_ants):
                path, alignment = self._construct_solution()
                score = _sp_score(alignment, self.seq_type, self.match, self.mismatch, self.gap_open, self.gap_extend)
                
                ant_alignments.append((score, alignment))
                ant_paths.append((score, path))
                
                if score > self.best_score:
                    self.best_score = score
                    self.best_alignment = alignment
            
            self._update_pheromones(ant_paths)
        
        return self.best_score, self.best_alignment

    def _construct_solution(self):
        current_state = tuple([0] * self.k)
        target_state = tuple(self.seq_lengths)
        
        alignment_cols = []
        path_edges = []
        
        prev_move = None
        
        while current_state != target_state:
            valid_moves = []
            for m in self.all_moves:
                # Check if move is valid 
                invalid = False
                for i in range(self.k):
                    if m[i] == 1 and current_state[i] == target_state[i]:
                        invalid = True
                        break
                if not invalid:
                    valid_moves.append(m)
            
            if not valid_moves:
                break
                
            # Calculate probabilities
            move_scores = []
            
            for m in valid_moves:
                tau = self.pheromones.get((current_state, m), self.tau_0)
                eta = self._calculate_heuristic(current_state, m, prev_move)
                move_scores.append((m, tau, eta))
            
            # Max-min bound eta for exp formulation if needed, or simply linear
            max_eta = max([x[2] for x in move_scores])
            min_eta = min([x[2] for x in move_scores])
            
            probs = []
            for m, tau, eta in move_scores:
                if max_eta == min_eta:
                    norm_eta = 1.0
                else:
                    norm_eta = (eta - min_eta) / (max_eta - min_eta)
                
                # add a tiny offset so log probability is well-behaved
                weight = (tau ** self.alpha) * ((norm_eta + 0.01) ** self.beta)
                probs.append(weight)
                
            sum_probs = sum(probs)
            if sum_probs == 0:
                probs = [1.0/len(probs)] * len(probs)
            else:
                probs = [p / sum_probs for p in probs]
                
            # Select move
            r = random.random()
            cumulative = 0.0
            selected_move = valid_moves[-1]
            for i, p in enumerate(probs):
                cumulative += p
                if r <= cumulative:
                    selected_move = valid_moves[i]
                    break
                    
            # Apply move
            path_edges.append((current_state, selected_move))
            
            # Construct column
            col = []
            new_state = list(current_state)
            for i in range(self.k):
                if selected_move[i] == 1:
                    col.append(self.seqs[i][current_state[i]])
                    new_state[i] += 1
                else:
                    col.append('-')
            
            alignment_cols.append(col)
            current_state = tuple(new_state)
            prev_move = selected_move
            
        # Reconstruct full alignments
        final_alignment = []
        for i in range(self.k):
            seq_aligned = ''.join([col[i] for col in alignment_cols])
            final_alignment.append(seq_aligned)
            
        return path_edges, final_alignment

    def _calculate_heuristic(self, state, move, prev_move):
        score = 0
        chars = []
        for i in range(self.k):
            if move[i] == 1:
                chars.append(self.seqs[i][state[i]])
            else:
                chars.append('-')
                
        # Pairwise column score
        for i in range(self.k):
            for j in range(i+1, self.k):
                c1 = chars[i]
                c2 = chars[j]
                if c1 == '-' and c2 == '-':
                    continue
                
                is_gap1 = (c1 == '-')
                is_gap2 = (c2 == '-')
                
                if is_gap1 and not is_gap2:
                    if prev_move is not None and prev_move[i] == 0:
                        score += self.gap_extend
                    else:
                        score += self.gap_open
                elif is_gap2 and not is_gap1:
                    if prev_move is not None and prev_move[j] == 0:
                        score += self.gap_extend
                    else:
                        score += self.gap_open
                elif not is_gap1 and not is_gap2:
                    if self.seq_type == 'PROTEIN':
                        score += get_blosum62_score(c1, c2)
                    else:
                        score += self.match if c1 == c2 else self.mismatch
        return score

    def _update_pheromones(self, ant_paths):
        # Evaporation
        keys_to_delete = []
        for key in self.pheromones:
            self.pheromones[key] *= (1.0 - self.rho)
            if self.pheromones[key] < 0.001:
                keys_to_delete.append(key)
        for key in keys_to_delete:
            del self.pheromones[key]
            
        # Pheromone deposit (Elitist updating to speed up, using top 10%)
        # Normalization of scores to properly weigh Delta Tau
        ant_paths.sort(key=lambda x: x[0], reverse=True)
        # Use top N ants
        top_n = max(1, self.num_ants // 4)
        best_ants = ant_paths[:top_n]
        
        # Max score for scaling
        max_s = best_ants[0][0]
        min_s = min([p[0] for p in ant_paths])
        diff = max_s - min_s if max_s != min_s else 1.0
        
        for score, path in best_ants:
            # We map score to a positive quantity
            delta_tau = (score - min_s) / diff * 2.0 + 0.1 # bounded [0.1, 2.1]
            for state, move in path:
                if (state, move) not in self.pheromones:
                    self.pheromones[(state, move)] = self.tau_0
                self.pheromones[(state, move)] += delta_tau
