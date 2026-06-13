# Projekt MSA — algorytm genetyczny

Metaheurystyka (algorytm genetyczny) dla problemu dopasowania wielu sekwencji
(Multiple Sequence Alignment). Rdzeń w C++17; Python tylko do wykresów.

## Kompilacja
    g++ -O2 -std=c++17 -o msa src/main.cpp

## Uruchomienie
    ./msa validate            # walidacja funkcji celu (przyklady z instrukcji)
    ./msa gen <args>          # generowanie instancji
    ./msa run instancja.txt   # pojedyncze uruchomienie GA
    ./msa sweep <args>        # testy parametrow
    ./msa classes <args>      # testy klas instancji
    ./msa converge <args>     # historia zbieznosci
    python3 plots.py          # wykresy PNG z plikow CSV (results/ -> figures/)

## Struktura
    src/scoring.hpp   funkcja celu SP, kary afiniczne, DNA + bialka
    src/matrices.hpp  BLOSUM62 + uproszczona tabela do walidacji
    src/align.hpp     reprezentacja, init, krzyzowanie, mutacje, naprawa
    src/ga.hpp        petla GA: selekcja, elityzm, budzet ocen
    src/main.cpp      interfejs CLI
    plots.py          generowanie wykresow
    gen_all.sh        generowanie pelnego zestawu instancji
    instances/        64 wygenerowane instancje
    results/          wyniki CSV + tabele + przyklady
    figures/          wykresy PNG

## Format pliku instancji
Pierwszy wiersz: DNA lub PROT. Kolejne wiersze: sekwencje (po jednej w wierszu).
