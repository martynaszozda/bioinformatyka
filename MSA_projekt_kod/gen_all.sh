#!/bin/sh
M=./msa
# --- Zestaw STROJENIA parametrow (maly, ale zroznicowany; DNA k=6 len=40) ---
$M gen --prefix tuneA --count 3 --k 6 --len 40 --psub 0.06 --pind 0.02 --seed 101
$M gen --prefix tuneB --count 3 --k 6 --len 40 --psub 0.15 --pind 0.05 --seed 202
$M gen --prefix tuneC --count 2 --k 6 --len 40 --psub 0.28 --pind 0.09 --seed 303
# --- Klasy PODOBIENSTWA (DNA k=6 len=40) ---
$M gen --prefix simEasy --count 6 --k 6 --len 40 --psub 0.05 --pind 0.02 --gaplen 1 --seed 111
$M gen --prefix simMed  --count 6 --k 6 --len 40 --psub 0.15 --pind 0.05 --gaplen 1 --seed 222
$M gen --prefix simHard --count 6 --k 6 --len 40 --psub 0.30 --pind 0.10 --gaplen 1 --seed 333
# --- Klasy STRUKTURY LUK (DNA k=6 len=50, podobne podobienstwo) ---
$M gen --prefix gapShort --count 6 --k 6 --len 50 --psub 0.12 --pind 0.12 --gaplen 1 --seed 444
$M gen --prefix gapLong  --count 6 --k 6 --len 50 --psub 0.12 --pind 0.03 --gaplen 5 --seed 555
# --- Klasy ROZMIARU (DNA) ---
$M gen --prefix szSmall --count 6 --k 6  --len 30  --psub 0.15 --pind 0.05 --seed 666
$M gen --prefix szMed   --count 6 --k 8  --len 80  --psub 0.15 --pind 0.05 --seed 777
$M gen --prefix szLarge --count 4 --k 10 --len 130 --psub 0.15 --pind 0.05 --seed 888
# --- Zestaw BIALKOWY (wersja 5.0; PROT k=6 len=50) ---
$M gen --prefix protSet --mode PROT --count 6 --k 6 --len 50 --psub 0.18 --pind 0.05 --seed 909
