#pragma once
#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace msa {

constexpr char GAP = '-';

//  Konfiguracja punktacji
enum class ScoreMode { DNA, PROTEIN, PROTEIN_SIMPLE };

struct Scoring {
    ScoreMode mode = ScoreMode::DNA;
    // Parametry dla DNA:
    int match    =  1;
    int mismatch = -1;
    // Kary za luki
    int gopen    = -2;   // otwarcie luki
    int gext     = -1;   // przedluzenie luki
    // Tablica podstawien 
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

//  Ocena pojedynczej pary wierszy (dwa napisy rownej dlugosci, z lukami '-').

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
                extend = false;     
            } else {
                if (gx) {
                    extend = (a[c-1] == GAP) && (b[c-1] != GAP);
                } else {
                    extend = (b[c-1] == GAP) && (a[c-1] != GAP);
                }
            }
            sc += extend ? s.gext : s.gopen;
        }
    }
    return sc;
}


//  SP-score calego dopasowania

inline long spScore(const std::vector<std::string>& rows, const Scoring& s) {
    long total = 0;
    const size_t k = rows.size();
    for (size_t i = 0; i < k; ++i)
        for (size_t j = i + 1; j < k; ++j)
            total += scorePair(rows[i], rows[j], s);
    return total;
}

}
