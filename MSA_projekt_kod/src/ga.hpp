// ============================================================================
//  ga.hpp
//  Algorytm genetyczny dla MSA.
//
//  Budzet obliczeniowy mierzymy LICZBA OCEN FUNKCJI CELU (evaluations), a nie
//  liczba generacji. Dzieki temu porownania miedzy roznymi wielkosciami populacji
//  sa uczciwe (kazda konfiguracja dostaje tyle samo wywolan funkcji celu).
// ============================================================================
#pragma once
#include "align.hpp"
#include "scoring.hpp"
#include <vector>
#include <random>
#include <algorithm>
#include <limits>

namespace msa {

struct GAParams {
    int    popSize    = 60;
    double pCross     = 0.85;
    double pMut       = 0.30;
    int    elitism    = 2;
    int    tournament = 3;
    long   maxEvals   = 40000;   // budzet: liczba ocen funkcji celu
    int    maxGen     = 100000;  // twardy limit generacji (zwykle nieaktywny)
};

struct GAResult {
    long initBest  = std::numeric_limits<long>::min();
    long finalBest = std::numeric_limits<long>::min();
    Align bestAlign;
    int  gens  = 0;
    long evals = 0;
    std::vector<long> histBest;   // best-so-far po kazdej generacji
    std::vector<long> histEvals;  // skumulowane oceny po kazdej generacji
};

// selekcja turniejowa: zwraca indeks zwyciezcy
inline int tournamentPick(const std::vector<long>& fit, int tsize, RNG& rng) {
    std::uniform_int_distribution<int> di(0, (int)fit.size() - 1);
    int best = di(rng);
    for (int t = 1; t < tsize; ++t) {
        int c = di(rng);
        if (fit[c] > fit[best]) best = c;
    }
    return best;
}

inline GAResult runGA(const std::vector<std::string>& seqs,
                      const Scoring& sc, const GAParams& P, RNG& rng) {
    GAResult R;
    const int N = P.popSize;

    std::vector<Align> pop(N);
    std::vector<long>  fit(N);
    for (int i = 0; i < N; ++i) {
        pop[i] = randomAlign(seqs, rng);
        fit[i] = spScore(pop[i], sc);
    }
    R.evals = N;

    auto bestIdx = [&]() {
        int b = 0; for (int i = 1; i < N; ++i) if (fit[i] > fit[b]) b = i; return b;
    };
    int bi = bestIdx();
    R.initBest  = fit[bi];
    R.finalBest = fit[bi];
    R.bestAlign = pop[bi];

    std::uniform_real_distribution<double> u01(0.0, 1.0);

    while (R.evals < P.maxEvals && R.gens < P.maxGen) {
        // posortuj indeksy malejaco po fitness (dla elityzmu)
        std::vector<int> ord(N);
        for (int i = 0; i < N; ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int a, int b){ return fit[a] > fit[b]; });

        std::vector<Align> next; next.reserve(N);
        std::vector<long>  nfit; nfit.reserve(N);

        // elityzm: przepisz najlepszych bez zmian (bez ponownej oceny)
        for (int e = 0; e < P.elitism && e < N; ++e) {
            next.push_back(pop[ord[e]]);
            nfit.push_back(fit[ord[e]]);
        }
        // reszta przez selekcje + krzyzowanie + mutacje
        while ((int)next.size() < N) {
            int a = tournamentPick(fit, P.tournament, rng);
            Align child;
            if (u01(rng) < P.pCross) {
                int b = tournamentPick(fit, P.tournament, rng);
                child = crossoverOnePoint(pop[a], pop[b], rng);
            } else {
                child = pop[a];
            }
            if (u01(rng) < P.pMut) mutate(child, rng);
            else repair(child); // kopie tez normalizujemy
            long f = spScore(child, sc);
            R.evals += 1;
            next.push_back(std::move(child));
            nfit.push_back(f);
        }

        pop.swap(next);
        fit.swap(nfit);
        R.gens += 1;

        int nb = bestIdx();
        if (fit[nb] > R.finalBest) { R.finalBest = fit[nb]; R.bestAlign = pop[nb]; }
        R.histBest.push_back(R.finalBest);
        R.histEvals.push_back(R.evals);
    }
    return R;
}

} // namespace msa
