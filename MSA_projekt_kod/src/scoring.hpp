// ============================================================================
//  scoring.hpp
//  Funkcja celu dla problemu MSA: SP-score (sum-of-pairs) z karami afinicznymi.
//
//  Implementacja jest zgodna z definicjami z instrukcji prowadzacego:
//    - match / mismatch (DNA) lub macierz podstawien BLOSUM62 (bialka),
//    - gapopen / gapextend (kary afiniczne za luki),
//    - gap vs gap = 0,
//    - ocena calego dopasowania = suma ocen wszystkich par wierszy (k po 2).
//
//  KLUCZOWA ZASADA KAR AFINICZNYCH (dla pary wierszy a,b w kolumnie c):
//    Jesli dokladnie jeden z wierszy ma luke (drugi ma litere), to jest to
//    gapextend WTEDY I TYLKO WTEDY, gdy w kolumnie c-1 TEN SAM wiersz mial luke,
//    a drugi wiersz mial litere. W przeciwnym razie jest to gapopen.
//    (Tym samym "kreska z kreska" w kolumnie poprzedniej zeruje serie luk –
//     pierwsza litera po takiej kresce to nowe otwarcie, nie przedluzenie.)
// ============================================================================
#pragma once
#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace msa {

constexpr char GAP = '-';

// ---------------------------------------------------------------------------
//  Konfiguracja punktacji
// ---------------------------------------------------------------------------
enum class ScoreMode { DNA, PROTEIN, PROTEIN_SIMPLE };

struct Scoring {
    ScoreMode mode = ScoreMode::DNA;
    // Parametry dla DNA:
    int match    =  1;
    int mismatch = -1;
    // Kary za luki (wspolne dla obu trybow):
    int gopen    = -2;   // otwarcie luki
    int gext     = -1;   // przedluzenie luki
    // Tablica podstawien (uzywana w trybach PROTEIN / PROTEIN_SIMPLE):
    // indeks po kodzie znaku (0..127), wartosc to indeks w macierzy lub -1.
    std::array<int8_t,128> aaIndex;
    std::vector<std::vector<int>> sub; // macierz podstawien NxN

    Scoring() { aaIndex.fill(-1); }
};

// Substytucja litera-litera (zaklada, ze oba znaki to NIE luka).
inline int substScore(char x, char y, const Scoring& s) {
    if (s.mode == ScoreMode::DNA) {
        return (x == y) ? s.match : s.mismatch;
    }
    int ix = s.aaIndex[(unsigned char)x];
    int iy = s.aaIndex[(unsigned char)y];
    return s.sub[ix][iy];
}

// ---------------------------------------------------------------------------
//  Ocena pojedynczej pary wierszy (dwa napisy rownej dlugosci, z lukami '-').
//  To jest dokladnie pairwise SP z karami afinicznymi.
// ---------------------------------------------------------------------------
inline long scorePair(const std::string& a, const std::string& b, const Scoring& s) {
    long sc = 0;
    const size_t L = a.size();
    for (size_t c = 0; c < L; ++c) {
        const char x = a[c], y = b[c];
        const bool gx = (x == GAP), gy = (y == GAP);
        if (!gx && !gy) {                 // litera vs litera
            sc += substScore(x, y, s);
        } else if (gx && gy) {            // luka vs luka
            sc += 0;
        } else {                          // dokladnie jedna luka
            // ustal, ktory wiersz ma luke (gapped), a ktory litere
            bool extend;
            if (c == 0) {
                extend = false;           // pierwsza kolumna => zawsze otwarcie
            } else {
                if (gx) {
                    // x ma luke teraz; przedluzenie tylko gdy x mial luke w c-1
                    // ORAZ y mial litere w c-1
                    extend = (a[c-1] == GAP) && (b[c-1] != GAP);
                } else {
                    // y ma luke teraz
                    extend = (b[c-1] == GAP) && (a[c-1] != GAP);
                }
            }
            sc += extend ? s.gext : s.gopen;
        }
    }
    return sc;
}

// ---------------------------------------------------------------------------
//  SP-score calego dopasowania: suma po wszystkich parach wierszy.
// ---------------------------------------------------------------------------
inline long spScore(const std::vector<std::string>& rows, const Scoring& s) {
    long total = 0;
    const size_t k = rows.size();
    for (size_t i = 0; i < k; ++i)
        for (size_t j = i + 1; j < k; ++j)
            total += scorePair(rows[i], rows[j], s);
    return total;
}

} // namespace msa
