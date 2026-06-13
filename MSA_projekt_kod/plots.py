#!/usr/bin/env python3
# Agregacja wynikow (CSV z programu C++) -> wykresy PNG + tabele tekstowe.
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = "results"; FIG = "figures"
os.makedirs(FIG, exist_ok=True)

# --- styl ---
ACCENT = "#1f8a70"; ACCENT2 = "#d98324"; ACCENT3 = "#3b6ea5"; GREY="#444"
plt.rcParams.update({
    "figure.dpi": 130, "font.size": 11, "axes.grid": True,
    "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans"
})

def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

def agg(rows, key, field):
    """grupuj po key, zwroc {val: np.array(field)} (posortowane po kluczu numerycznym jesli mozna)"""
    d = {}
    for r in rows:
        d.setdefault(r[key], []).append(float(r[field]))
    out = {k: np.array(v) for k, v in d.items()}
    return out

def numkey(k):
    try: return float(k)
    except: return k

tables = []  # zbieramy teksty tabel

# ===========================================================================
# 1-3. SWEEPY PARAMETROW (DNA)
# ===========================================================================
def sweep_plot(csvname, xlabel, title, outpng, color):
    rows = read_csv(f"{RES}/{csvname}")
    finals = agg(rows, "value", "final_best")
    improv = agg(rows, "value", "improvement")
    xs = sorted(finals.keys(), key=numkey)
    xn = [numkey(x) for x in xs]
    means = [finals[x].mean() for x in xs]
    stds  = [finals[x].std()  for x in xs]
    bests = [finals[x].max()  for x in xs]
    nruns = sum(len(finals[x]) for x in xs)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.errorbar(xn, means, yerr=stds, marker="o", color=color, capsize=4,
                lw=2, label="srednia ± odch. std.")
    ax.plot(xn, bests, marker="s", ls="--", color=GREY, lw=1.2, alpha=0.8,
            label="najlepszy")
    ax.set_xlabel(xlabel); ax.set_ylabel("koncowy SP-score")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/{outpng}"); plt.close(fig)

    # tabela
    lines = [f"\n### {title}",
             "wartosc | srednia SP | odch.std | najlepszy SP | srednia poprawa"]
    for x in xs:
        lines.append(f"{x} | {finals[x].mean():.1f} | {finals[x].std():.1f} | "
                     f"{finals[x].max():.0f} | {improv[x].mean():.1f}")
    lines.append(f"(zestaw strojeniowy: 8 instancji x {len(finals[xs[0]])//8 if len(xs)>0 else 0} "
                 f"uruchomien = {len(finals[xs[0]])} wynikow na wartosc)")
    tables.append("\n".join(lines))
    return xs, means, stds

sweep_plot("sweep_pop.csv", "liczność populacji",
           "Wpływ liczności populacji na końcowy SP\n(DNA, zestaw strojeniowy, 8 inst. x 10 uruch., budżet 40000 ocen)",
           "fig_sweep_pop.png", ACCENT)
sweep_plot("sweep_pm.csv", "prawdopodobieństwo mutacji  p_mut",
           "Wpływ p_mut na końcowy SP\n(DNA, zestaw strojeniowy, 8 inst. x 10 uruch., budżet 40000 ocen)",
           "fig_sweep_pm.png", ACCENT2)
sweep_plot("sweep_pc.csv", "prawdopodobieństwo krzyżowania  p_cross",
           "Wpływ p_cross na końcowy SP\n(DNA, zestaw strojeniowy, 8 inst. x 10 uruch., budżet 40000 ocen)",
           "fig_sweep_pc.png", ACCENT3)

# ===========================================================================
# 4. KLASY PODOBIENSTWA
# ===========================================================================
def class_plot(csvname, order, labels, title, outpng, ylab="koncowy SP-score",
               field="final_best"):
    rows = read_csv(f"{RES}/{csvname}")
    d = agg(rows, "class", field)
    means = [d[c].mean() for c in order]
    stds  = [d[c].std()  for c in order]
    bests = [d[c].max()  for c in order]
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    bars = ax.bar(x, means, yerr=stds, capsize=5, color=[ACCENT,ACCENT2,ACCENT3][:len(order)],
                  alpha=0.85, edgecolor="black", lw=0.6)
    ax.scatter(x, bests, color=GREY, marker="D", zorder=5, label="najlepszy")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(ylab); ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/{outpng}"); plt.close(fig)
    # tabela
    impr = agg(rows, "class", "improvement")
    lines = [f"\n### {title}",
             "klasa | srednia SP | odch.std | najlepszy SP | srednia poprawa | n"]
    for c,l in zip(order,labels):
        lines.append(f"{l} | {d[c].mean():.1f} | {d[c].std():.1f} | {d[c].max():.0f} | "
                     f"{impr[c].mean():.1f} | {len(d[c])}")
    tables.append("\n".join(lines))

class_plot("cls_similarity.csv", ["simEasy","simMed","simHard"],
           ["łatwe","średnie","trudne"],
           "Jakość końcowa wg podobieństwa sekwencji\n(DNA k=6, L=40, 6 inst./klasę x 10 uruch.)",
           "fig_cls_similarity.png")

# ===========================================================================
# 5. KLASY STRUKTURY LUK
# ===========================================================================
class_plot("cls_gap.csv", ["gapShort","gapLong"],
           ["wiele krótkich przerw","nieliczne długie bloki"],
           "Jakość końcowa wg struktury luk\n(DNA k=6, L=50, 6 inst./klasę x 10 uruch.)",
           "fig_cls_gap.png")

# ===========================================================================
# 6. KLASY ROZMIARU  (SP nieporownywalny -> pokazujemy poprawe)
# ===========================================================================
class_plot("cls_size.csv", ["szSmall","szMed","szLarge"],
           ["małe\n(k6,L30)","średnie\n(k8,L80)","duże\n(k10,L130)"],
           "Średnia poprawa względem rozwiązania startowego wg rozmiaru\n(DNA, 6/6/4 inst. x 10 uruch., budżet 60000 ocen)",
           "fig_cls_size.png", ylab="poprawa SP (final − init)", field="improvement")

# ===========================================================================
# 7-8. KRZYWE ZBIEZNOSCI
# ===========================================================================
def conv_plot(csvname, title, outpng, color):
    rows = read_csv(f"{RES}/{csvname}")
    ev = np.array([float(r["evals"]) for r in rows])
    m  = np.array([float(r["best_mean"]) for r in rows])
    s  = np.array([float(r["best_std"]) for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(ev, m, color=color, lw=2, label="średni best-so-far")
    ax.fill_between(ev, m-s, m+s, color=color, alpha=0.18, label="± odch. std.")
    ax.set_xlabel("liczba ocen funkcji celu"); ax.set_ylabel("najlepszy SP-score")
    ax.set_title(title, fontsize=11); ax.legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/{outpng}"); plt.close(fig)

conv_plot("conv_dna.csv",
          "Zbieżność GA (DNA)\n(zestaw strojeniowy, 8 inst. x 10 uruch., uśrednione)",
          "fig_conv_dna.png", ACCENT)
conv_plot("conv_prot.csv",
          "Zbieżność GA (białka, BLOSUM62)\n(protSet, 6 inst. x 10 uruch., uśrednione)",
          "fig_conv_prot.png", ACCENT2)

# ===========================================================================
# 9. BIALKA: sweep mutacji
# ===========================================================================
sweep_plot("sweep_pm_prot.csv", "prawdopodobieństwo mutacji  p_mut",
           "Wpływ p_mut na końcowy SP (białka, BLOSUM62)\n(protSet, 6 inst. x 10 uruch., budżet 50000 ocen)",
           "fig_sweep_pm_prot.png", ACCENT2)

# ===========================================================================
# Tabela poprawy startowej (DNA sim classes) + timing
# ===========================================================================
tim = read_csv(f"{RES}/timing.csv")
lines = ["\n### Reprezentatywne czasy działania (1 rdzeń)",
         "klasa | budżet ocen | średni czas [s]"]
for r in tim:
    lines.append(f"{r['klasa']} | {r['budzet_evals']} | {r['sredni_czas_s']}")
tables.append("\n".join(lines))

with open(f"{RES}/tables.txt", "w") as f:
    f.write("\n".join(tables))

print("Wygenerowano wykresy:")
for p in sorted(os.listdir(FIG)): print("  figures/"+p)
print("\n=== TABELE ZBIORCZE ===")
print("\n".join(tables))
