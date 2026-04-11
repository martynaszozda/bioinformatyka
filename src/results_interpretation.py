import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('results.csv')

def parse_dataset_info(name):
    parts = name.replace('.fasta', '').split('_')
    size = parts[1] if len(parts) > 1 else 'unknown'
    similarity = parts[2] if len(parts) > 2 else 'medium' 
    return pd.Series([size, similarity])

df[['Size', 'Similarity']] = df['Dataset'].apply(parse_dataset_info)

df['Size'] = pd.Categorical(df['Size'], categories=['small', 'medium', 'large'], ordered=True)
df['Similarity'] = pd.Categorical(df['Similarity'], categories=['high', 'medium', 'low'], ordered=True)

sns.set_theme(style="whitegrid")

# --- WYKRES 1: Skalowalność (Small vs Medium vs Large) ---
plt.figure(figsize=(10, 6))
sns.boxplot(data=df, x='Size', y='Score', hue='SeqType')
plt.title('Skalowalność algorytmu: Jakość dopasowania względem rozmiaru sekwencji')
plt.ylabel('SP-Score')
plt.savefig('size_comparison.png')

# --- WYKRES 2: Trudność biologiczna (High vs Medium vs Low Similarity) ---
plt.figure(figsize=(10, 6))
sns.violinplot(data=df, x='Similarity', y='Score', hue='SeqType', split=True)
plt.title('Wpływ ewolucyjnego podobieństwa sekwencji na wyniki ACO')
plt.ylabel('SP-Score')
plt.savefig('similarity_impact.png')

# --- WYKRES 3: Efektywność czasowa (Ile sekund na jaki rozmiar) ---
plt.figure(figsize=(10, 6))
sns.lineplot(data=df, x='Size', y='Time(s)', hue='SeqType', marker='o')
plt.title('Koszt obliczeniowy: Czas wykonania vs Rozmiar instancji')
plt.yscale('log')

plt.savefig('execution_time.png')

print("Nowe analizy zapisane: size_comparison.png, similarity_impact.png, execution_time.png")