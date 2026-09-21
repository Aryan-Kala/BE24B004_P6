from pathlib import Path

import numpy as np

def calculate_entropy_robust(filename):
    sequences = {}
    is_fasta = False
    current_id = None
    
    with open(filename, 'r') as f:
        for line in f:
            line_str = line.rstrip()
            if not line_str:
                continue
            
            # Check if it's a FASTA file
            if line_str.startswith('>'):
                is_fasta = True
                current_id = line_str[1:].split()[0]
                sequences[current_id] = ""
                continue
                
            if is_fasta:
                if current_id:
                    sequences[current_id] += line_str.strip()
            else:
                # Clustal format block parsing
                if (line_str.startswith('CLUSTAL') or 
                    line_str.startswith('*') or 
                    line_str.startswith(':') or 
                    line_str.startswith('.') or
                    set(line_str.strip()) <= set(' *:.')):
                    continue
                
                parts = line_str.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    seq_chunk = parts[1]
                    # Ensure it's sequence data and not a stray number or header
                    if any(c.isalpha() or c == '-' for c in seq_chunk):
                        if seq_id not in sequences:
                            sequences[seq_id] = ""
                        sequences[seq_id] += seq_chunk

    msa = list(sequences.values())
    if not msa:
        raise ValueError(f"Could not extract sequences from {filename}")
        
    n = len(msa)
    length = len(msa[0])
    residues = ["-", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    
    entropy_scores = []
    for i in range(length):
        freq = {residue: 0 for residue in residues}
        count = 0
        for j in range(n):
            if i < len(msa[j]) and msa[j][i] != "-":
                freq[msa[j][i]] += 1
                count += 1
        
        if count > 0:
            freq = {residue: freq[residue] / count for residue in residues}
        
        entropy_score = 0
        for residue in residues:
            if freq[residue] != 0:
                entropy_score -= freq[residue] * np.log(freq[residue])
                
        entropy_scores.append(round(entropy_score, 3))
        
    return entropy_scores

script_dir = Path(__file__).resolve().parent
clustal_file = script_dir / "clustal_set1.fasta.aln-clustal_num"
mafft_file = script_dir / "mafft_set1.fasta.fasta"
muscle_file = script_dir / "muscle_set1.fasta.aln-clustalw"

c_scores = calculate_entropy_robust(clustal_file)
m_scores = calculate_entropy_robust(mafft_file)
mu_scores = calculate_entropy_robust(muscle_file)

print("--- POSITION-BY-POSITION COMPARISON ---")
similar_count = 0
diff_count = 0

for pos, (c, m, mu) in enumerate(zip(c_scores, m_scores, mu_scores)):
    max_diff = max(abs(c - m), abs(c - mu), abs(m - mu))
    
    if max_diff <= 0.05:
        print(f"Position {pos + 1}: SIMILAR (Clustal: {c}, MAFFT: {m}, MUSCLE: {mu})")
        similar_count += 1
    else:
        print(f"Position {pos + 1}: DIFFERENT (Clustal: {c}, MAFFT: {m}, MUSCLE: {mu})")
        diff_count += 1

print(f"\nTotal Similar Positions: {similar_count}")
print(f"Total Different Positions: {diff_count}")