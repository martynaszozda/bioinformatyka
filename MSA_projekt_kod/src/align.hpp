#pragma once
#include "scoring.hpp"
#include <string>
#include <vector>
#include <random>
#include <algorithm>

namespace msa {

using Align = std::vector<std::string>;
using RNG   = std::mt19937;

// ---- pomocnicze ----
inline int countLetters(const std::string& r) {
    int c = 0; for (char ch : r) if (ch != GAP) ++c; return c;
}
inline size_t alnLen(const Align& A) { return A.empty() ? 0 : A[0].size(); }

// Wyrownaj dlugosci wierszy przez dopelnienie lukami z prawej.
inline void padToEqual(Align& A) {
    size_t L = 0; for (auto& r : A) L = std::max(L, r.size());
    for (auto& r : A) if (r.size() < L) r.append(L - r.size(), GAP);
}

// Usun kolumny zlozone wylacznie z luk
inline void removeAllGapColumns(Align& A) {
    if (A.empty()) return;
    size_t L = A[0].size();
    std::string keep; keep.reserve(L);
    for (size_t c = 0; c < L; ++c) {
        bool allGap = true;
        for (auto& r : A) if (r[c] != GAP) { allGap = false; break; }
        keep.push_back(allGap ? '0' : '1');
    }
    for (auto& r : A) {
        std::string nr; nr.reserve(L);
        for (size_t c = 0; c < L; ++c) if (keep[c] == '1') nr.push_back(r[c]);
        r.swap(nr);
    }
}

// Pelna naprawa rozwiazania po operatorach.
inline void repair(Align& A) {
    padToEqual(A); //wyrównanie wierszy do najdłuższego, dodając przerwy z prawej
    removeAllGapColumns(A);  //usunięcie pustych kolumn (same przerwy)
    padToEqual(A);
}

// Sprawdzenie poprawnosci wzgledem sekwencji wejsciowych (do testow)
inline bool isValid(const Align& A, const std::vector<std::string>& seqs) {
    if (A.size() != seqs.size()) return false;
    size_t L = alnLen(A);
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i].size() != L) return false;
        std::string stripped;
        for (char ch : A[i]) if (ch != GAP) stripped.push_back(ch);
        if (stripped != seqs[i]) return false;
    }
    // brak kolumn samych luk
    for (size_t c = 0; c < L; ++c) {
        bool allGap = true;
        for (auto& r : A) if (r[c] != GAP) { allGap = false; break; }
        if (allGap) return false;
    }
    return true;
}

//inicjalizacja
inline Align randomAlign(const std::vector<std::string>& seqs, RNG& rng, double slack = 1.30) {
    size_t maxLen = 0; for (auto& s : seqs) maxLen = std::max(maxLen, s.size());
    size_t L = std::max<size_t>(maxLen + 1, (size_t)(maxLen * slack));
    Align A(seqs.size());
    for (size_t i = 0; i < seqs.size(); ++i) {
        const std::string& s = seqs[i];
        int n = (int)s.size();
        // wybierz n rosnacych pozycji z [0, L)
        std::vector<int> pos(L); for (size_t j = 0; j < L; ++j) pos[j] = (int)j;
        std::shuffle(pos.begin(), pos.end(), rng);
        pos.resize(n);
        std::sort(pos.begin(), pos.end());
        std::string row(L, GAP);
        for (int j = 0; j < n; ++j) row[pos[j]] = s[j];
        A[i] = std::move(row);
    }
    repair(A);
    return A;
}

// krzyżowanie
inline Align crossoverOnePoint(const Align& P1, const Align& P2, RNG& rng) {
    size_t L1 = alnLen(P1);
    if (L1 < 2) return P1;
    std::uniform_int_distribution<size_t> dc(1, L1 - 1);
    size_t c = dc(rng);
    Align child(P1.size());
    for (size_t i = 0; i < P1.size(); ++i) {
        std::string prefix = P1[i].substr(0, c);
        int ni = countLetters(prefix);
        int cnt = 0; size_t p = 0;
        const std::string& r2 = P2[i];
        while (p < r2.size() && cnt < ni) { if (r2[p] != GAP) ++cnt; ++p; }
        std::string suffix = r2.substr(p);
        child[i] = prefix + suffix;
    }
    repair(child);
    return child;
}

// mutacje
inline void mutShift(Align& A, RNG& rng) {
    size_t k = A.size(), L = alnLen(A);
    if (L < 2) return;
    std::uniform_int_distribution<size_t> dr(0, k - 1);
    size_t i = dr(rng);
    std::vector<std::pair<size_t,size_t>> sw;
    for (size_t c = 0; c + 1 < L; ++c) {
        if (A[i][c] == GAP && A[i][c+1] != GAP) sw.push_back({c, c+1});
        if (A[i][c] != GAP && A[i][c+1] == GAP) sw.push_back({c+1, c});
    }
    if (sw.empty()) return;
    auto pr = sw[std::uniform_int_distribution<size_t>(0, sw.size()-1)(rng)];
    std::swap(A[i][pr.first], A[i][pr.second]);
}

// M2: relokacja bloku luk
inline void mutGapBlockMove(Align& A, RNG& rng) {
    size_t k = A.size(), L = alnLen(A);
    if (L < 2) return;
    std::uniform_int_distribution<size_t> dr(0, k - 1);
    size_t i = dr(rng);
    std::vector<std::pair<size_t,size_t>> blocks;
    size_t c = 0;
    while (c < L) {
        if (A[i][c] == GAP) {
            size_t st = c; while (c < L && A[i][c] == GAP) ++c;
            blocks.push_back({st, c - st});
        } else ++c;
    }
    if (blocks.empty()) return;
    auto blk = blocks[std::uniform_int_distribution<size_t>(0, blocks.size()-1)(rng)];
    std::string letters; for (char ch : A[i]) if (ch != GAP) letters.push_back(ch);
    int n = (int)letters.size();
    int glen = (int)blk.second;
    int insPos = std::uniform_int_distribution<int>(0, n)(rng);
    std::string row;
    row.reserve(L);
    for (int j = 0; j < insPos; ++j) row.push_back(letters[j]);
    row.append(glen, GAP);
    for (int j = insPos; j < n; ++j) row.push_back(letters[j]);
    if (row.size() < L) row.append(L - row.size(), GAP);
    A[i] = row;
}

// M3: gap insertion – nowa luka w losowym podzbiorze wierszy, a w pozostalych wstaw luke w innym losowym miejscu (wyrownanie dlugosci)
inline void mutGapInsertion(Align& A, RNG& rng) {
    size_t k = A.size(), L = alnLen(A);
    if (k < 2) return;
    std::uniform_int_distribution<size_t> dpos(0, L);
    size_t p1 = dpos(rng), p2 = dpos(rng);
    std::uniform_int_distribution<int> coin(0, 1);
    std::vector<char> inS(k);
    int cntS = 0;
    for (size_t i = 0; i < k; ++i) { inS[i] = (char)coin(rng); cntS += inS[i]; }
    if (cntS == 0) inS[0] = 1;
    if (cntS == (int)k) inS[0] = 0;
    for (size_t i = 0; i < k; ++i) {
        size_t p = inS[i] ? p1 : p2;
        if (p > A[i].size()) p = A[i].size();
        A[i].insert(A[i].begin() + p, GAP);
    }
}


inline void mutate(Align& A, RNG& rng) {
    int r = std::uniform_int_distribution<int>(0, 99)(rng);
    if      (r < 45) mutShift(A, rng);        // 45% drobne przesuniecia
    else if (r < 75) mutGapBlockMove(A, rng); // 30% relokacja bloku luk
    else             mutGapInsertion(A, rng); // 25% wstawienie nowej luki
    repair(A);
}

}

