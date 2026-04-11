# Multiple Sequence Alignment


wersja 1.0. 1 21. 03 .20 26

## 1. Wprowadzenie do problemu dopasowania wielu sekwencji

## (MSA)

Jednym z podstawowych typów danych w bioinformatyce są sekwencje biologiczne: sekwencje DNA,
RNA oraz białek. Z punktu widzenia informatyka można je traktować jako napisy zbudowane nad
niewielkim alfabetem. Dla DNA alfabetem jest zwykle zbiór {A, C, G, T}, dla RNA {A, C, G, U}, a dla białek
zbiór dwudziestu podstawowych aminokwasów. Choć taka reprezentacja wydaje się bardzo prosta,
analiza i porównywanie sekwencji pozwala wnioskować o funkcji genów i białek, pokrewieństwie
organizmów, konserwowanych fragmentach ewolucyjnych, a także o możliwych mutacjach związanych
z chorobami. Z tego powodu dopasowywanie sekwencji należy do najważniejszych problemów
obliczeniowych w bioinformatyce.

Najprostszą wersją tego zadania jest dopasowanie pary sekwencji. Jego celem jest takie ustawienie
dwóch napisów względem siebie, aby podobne lub identyczne pozycje znalazły się w tych samych
kolumnach. W praktyce wolno przy tym wstawiać znaki przerwy (gaps), najczęściej oznaczane jako „-”,
które reprezentują insercje lub delecje powstałe w procesie ewolucji. Dopasowanie ocenia się przy
użyciu funkcji punktowej: zgodność znaków zwykle zwiększa wynik, niezgodność go obniża, a
wstawienie luki jest karane. Już na tym poziomie widać wyraźnie informatyczny charakter problemu:
trzeba zdefiniować reprezentację rozwiązania, funkcję celu oraz znaleźć algorytm wyznaczający
rozwiązanie najlepsze lub przynajmniej bardzo dobre.

W wielu zastosowaniach porównanie tylko dwóch sekwencji nie wystarcza. W biologii i medycynie
często dysponujemy całymi rodzinami sekwencji pochodzących od różnych organizmów, genów lub
wariantów tego samego białka. Wtedy naturalnym rozszerzeniem jest **MSA (Multiple Sequence
Alignment)** czyli dopasowanie wielu sekwencji jednocześnie. Celem nie jest już tylko stwierdzenie, czy
dwie sekwencje są podobne, ale znalezienie wspólnego układu kolumn dla zbioru danych. Dzięki temu
można wykrywać pozycje silnie konserwatywne, identyfikować motywy funkcjonalne a także
analizować zmienność genetyczną w większych populacjach.

Szczególnie ważne jest spojrzenie na MSA jako na problem optymalizacyjny. Rozwiązaniem jest tutaj
pewne uporządkowanie sekwencji z dopuszczonymi lukami, natomiast jakość rozwiązania mierzy się
za pomocą funkcji oceny. Jedną z najczęściej stosowanych jest metoda sum-of-pairs, w której całe
dopasowanie ocenia się jako sumę ocen wszystkich dopasowań parami w każdej kolumnie. Taki model
jest intuicyjny i dobrze nadaje się do analizy algorytmicznej, ale szybko okazuje się, że przestrzeń
możliwych rozwiązań rośnie bardzo gwałtownie wraz z liczbą sekwencji i ich długością. O ile dla dwóch
sekwencji istnieją klasyczne algorytmy dynamiczne wyznaczające rozwiązanie optymalne, o tyle dla
wielu sekwencji problem staje się znacznie trudniejszy obliczeniowo i w praktyce wymaga stosowania
heurystyk oraz metod przybliżonych.


Z tego powodu w bioinformatyce ogromne znaczenie mają algorytmy, które nie gwarantują optimum
globalnego, ale potrafią znaleźć dobre dopasowania w rozsądnym czasie. Wśród typowych podejść
znajdują się metody progresywne, strategie oparte na profilach, algorytmy iteracyjne oraz różne
techniki inspirowane optymalizacją kombinatoryczną. Z informatycznego punktu widzenia MSA jest
więc bardzo ciekawym przykładem zadania, w którym spotykają się: modelowanie danych
symbolicznych, projektowanie funkcji oceny, analiza złożoności obliczeniowej oraz konstruowanie
praktycznych heurystyk.

Warto też podkreślić, że dopasowanie sekwencji nie jest celem samym w sobie. Jest ono zwykle etapem
pośrednim w większym procesie analitycznym. Wynik MSA może służyć jako wejście do dalszych metod:
wykrywania konserwowanych regionów, przewidywania domen funkcjonalnych, analizy zależności
ewolucyjnych czy budowy modeli statystycznych opisujących rodzinę sekwencji. Innymi słowy,
poprawność i jakość dopasowania wpływa później na jakość kolejnych wniosków biologicznych.
Dlatego problem MSA ma znaczenie zarówno praktyczne, jak i dydaktyczne: dobrze pokazuje, w jaki
sposób narzędzia informatyczne są wykorzystywane do rozwiązywania realnych problemów nauk
przyrodniczych.

W dalszej części materiału problem MSA będzie traktowany formalnie jako zadanie optymalizacyjne:
zostanie zdefiniowana reprezentacja dopasowania, funkcja oceny oraz sposób obliczania wartości
rozwiązania na konkretnym przykładzie. Takie ujęcie pozwala przejść od intuicyjnego porównywania
sekwencji do ścisłego modelu, który można analizować, implementować i rozwiązywać metodami
algorytmicznymi.

## 2. Przykłady

Założmy, że mamy dopasować takie 4 sekwencje ( _w1-w4 : wiersze od 1 do 4)._

AGTCGTAG w1 **A G T C G T A G**
ATCGTCG ➔ dopasowujemy, np. tak: ➔ w2 **A – T C G T C G**
GTAG w3 **– – – – G T A G**
GTAGAG w4 **– G T A G – A G**

Należy teraz jakoś ocenić jakość tego dopasowania. Potrzebujemy częściowych punktów:

- match, +1 punkt, _np. A z A, G z G, itd._
- mismatch, -1 punkt, _np. A z G_
- gapopen, gopen = -2 punkty – „ _otwierający” znak – , np. wiersze 3 i 4 z pierwszej kolumny, ale też_
    **_druga_** _kolumna wiersza w2 oraz_ **_szósta_** _kolumna dla wiersza w4 – to są „otwarcia” przerwy_
    _pomiędzy nukleotydami_
- gapextend, gext = -1 punkt, _w przykładzie to druga, trzecia i czwarta kolumna wiersza w_
- gap vs gap: 0 (zero) punktów

No dobrze, ale jak oceniamy? Ocena za całe dopasowanie to suma ocen częściowych każdy-z-
każdym, określa tą liczbę dwumian Newtona, w naszym przypadku „4 na 2”, gdzie 4 to liczba wierszy.
Musimy ocenić dopasowania w1-w2 potem w1-w3, w1-w4, w2-w3, w2-w4 i finalnie w3-w4.


**Ocena za dopasowania w1 z w2:** ( **A G T C G T A G** oraz **A – T C G T C G** )

match A-A: +
gopen G „–” : -2 punkty
match T-T: +
match C-C: +

```
match G-G: +
match T-T: +
mismatch A-C: - 1
match G-G: +
```
Czyli suma w1-w2 to: +1 – 2 +1 +1 +1 +1 - 1 +1 = 3.

**Ocena za dopasowania w1 z w3:** ( **A G T C G T A G** oraz **– – – – G T A G** )

gopen A „–”: -2 punkty
gext G „–”: -1 punkt
gext T „–”: - 1
gext C „–”: - 1

```
match G-G : +
match T-T : +
match A-A : +
match G-G : + 1
```
Suma w1-w3: - 2 - 1 - 1 - 1 +1 +1 +1 +1 = -5 + 4 = - 1

**Ocena za dopasowania w1 z w4:** ( **A G T C G T A G** oraz **– G T A G – A G** )

gopen A i „–” : -2 punkty
match G-G: +
match T-T: +
mismatch C-A: - 1

```
match G-G: +
gopen T i „–” : -2 punkty
match A-A: +
match G-G: + 1
```
Suma w1-w4: - 2 +1 +1 -1 +1 - 2 +1 +1 = 0

**Ocena za dopasowania w2 z w3:** ( **A – T C G T C G** oraz **– – – – G T A G** )

gopen A i „–” : -2 punkty
porównanie „–” z „–” : 0 punktów
gopen T i „–” : -2 punkty
gext C i „–” : -1 punkt

```
match G-G: +
match T-T: +
mismatch C-A: - 1
match G-G: +
```
Suma w2-w3: - 2 +0 - 2 - 1 +1 +1 - 1 +1 = - 3

**Ocena za dopasowania w2 z w4:** ( **A – T C G T C G** oraz **– G T A G – A G** )

gopen A i „–” : -2 punkty
gopen „–” i G : -2 punkty
match T-T: +
mismatch C-A: - 1

```
match G-G: +
gopen T i „–” : -2 punkty
mismatch C-A: - 1
match G-G: +
```
Suma w2-w4: - 2 - 2 +1 -1 +1 - 2 - 1 +1 = - 5

**Ocena za dopasowania w3 z w4:** ( **– – – – G T A G** oraz **– G T A G – A G** )

porównanie „–” z „–” : 0 punktów
gopen „–” i G : -2 punkty
gext „–” i T : -1 punkt
gext „–” i A : -1 punkt

```
match G-G: +
gopen T i „–” : -2 punkty
match A-A: +
match G-G: +
```
Suma w3-w4: +0 - 2 - 1 - 1 +1 - 2 +1 +1 = - 3

Ocena calkowita dopasowania: **SP** = 3 – 1 + 0 – 3 – 5 – 3 = **- 9**


Czy to dużo czy to mało? Nie wiemy, możemy co najwyżej wymyślić inne dopasowania i porównywać
wyniki. Naszym celem jest oczywiście maksymalizacja, czyli najlepsze dopasowania to te z
największą liczbą, nawet gdyby była ujemna.

### 2.1 MSA dla aminokwasów

W przypadku sekwencji białkowych ocena dopasowania jest zwykle bardziej realistyczna niż dla
prostego schematu „zgodność/niezgodność”, ponieważ nie wszystkie zamiany aminokwasów są
biologicznie jednakowo prawdopodobne. Z tego powodu stosuje się macierze podstawień takie jak
**BLOSUM62** , które przypisują punktację każdej parze aminokwasów. Wysokie wartości dodatnie
odpowiadają zamianom częstym i biologicznie „akceptowalnym”, natomiast wartości ujemne
oznaczają podstawienia rzadsze lub mniej prawdopodobne. Kary za luki nie są częścią samej macierzy
i muszą zostać określone osobno. W poniższym przykładzie pokazano sposób oceny dopasowania
sekwencji białkowych z użyciem fragmentu macierzy BLOSUM62 oraz przyjętego modelu kar za luki, tak
aby zilustrować, jak zmienia się punktacja w porównaniu z prostszymi przykładami dla sekwencji DNA.

**BLOSUM62 –** uproszczona tabela dla przykładu z kilkoma aminokwasami:

#### A^ C^ G^ N^ P^ S^ T^ –^

```
A 4 0 0 - 2 - 1 1 0 - 4
C 0 9 - 3 - 3 - 3 - 1 - 1 - 4
G 0 - 3 6 0 - 2 0 - 2 - 4
N - 2 - 3 0 6 - 2 1 0 - 4
P - 1 - 3 - 2 - 2 7 - 1 - 1 - 4
S 1 - 1 0 1 - 1 4 1 - 4
T 0 - 1 - 2 0 - 1 1 5 - 4
```
**- -**^4 **-**^4 **-**^4 **-**^4 **-**^4 **-**^4 **-**^4 

**Przykład:** 5 sekwencji, przykładowe dopasowanie.

#### AGTCGTAGNPST w1: A G T C G T A G N P S T

#### ASTCGTAGPST w2: A S T C G T A G – P S T

#### GTPGAGNST w3: – G T P G – A G N – S T

#### ANCGTNPT w4: A – N C G T – – N P – T

#### GATAPS w5: – G – – A T A – – P S –


Pierwsze cztery dopasowania z (^52 )= 10 (każdy wiersz z każdym):

**Ocena za dopasowanie w1 z w2:** ( A G T C G T A G N P S T oraz A S T C G T A G – P S T)

A-A: +
G-S: 0
T-T: +
C-C: +

```
G-G: +
T-T: +
A-A: +
G-G: +
```
```
gap N i „–” : -4 punkty
P-P: +
S-S: +
T-T: +
```
Czyli suma w1-w2: +4 +0 +5 +9 +6 +5 +4 +6 -4 +7 +4 +5 = **51**

**Dopasowanie w1 z w3:** ( A G T C G T A G N P S T oraz – G T P G – A G N – S T)

gap A i „–” : -4 punkty
G-G: +
T-T: +
C-P: - 3

```
G-G: +
gap T i „–” : -4 punkty
A-A: +
G-G: +
```
```
N-N: +
gap P i „–” : -4 punkty
S-S: +
T-T: +
```
Czyli suma w1-w3: - 4 +6 +5 -3 +6 -4 +4 +6 +6 -4 +4 +5 = **27**

**Dopasowanie w1 z w4:** ( A G T C G T A G N P S T oraz A – N C G T – – N P – T)

A-A: +
gap G i „–” : -4 punkty
T-N: 0
C-C: +

```
G-G: +
T-T: +
gap A i „–” : -4 punkty
gap G i „–” : -4 punkty
```
```
N-N: +
P-P: +
gap S i „–” : -4 punkty
T-T: +
```
Czyli suma w1-w4: +4 -4 +0 +9 +6 +5 - 4 - 4 +6 +7 -4 +5 = **26**

**Dopasowanie w1 z w5:** ( A G T C G T A G N P S T oraz – G – – A T A – – P S – )

gap A i „–” : -4 punkty
G-G: +
gap T i „–” : -4 punkty
gap C i „–” : -4 punkty

```
G-A: 0
T-T: +
A-A: +
gap G i „–” : -4 punkty
```
```
gap N i „–” : -4 punkty
P-P: +
S-S: +
gap T i „–” : -4 punkty
```
Czyli suma w1-w5: - 4 +6 - 4 - 4 +0 +5 +4 - 4 - 4 +7 +4 -4 = **2**

**Pozostałe pary wierszy:**

- suma **w2-w3** = **11**
- suma **w2-w4** = **16**
    - suma **w2-w5** = **0**
    - suma **w3-w4** = **- 14**
       - suma **w3-w5** = **- 14**
       - suma **w4-w5** = **- 20**

**Wszystkie “części” dopasowania** (dla 5 sekwencji mamy łącznie 10 par):

- w1-w2 = 51
- w1-w3 = 27
- w1-w4 = 26
    - w1-w5 = 2
    - w2-w3 = 11
    - w2-w4 = 16
       - w2-w5 = 0
       - w3-w4 = - 14
       - w3-w5 = - 14
          - w4-w5 = - 20

Zatem wartość całego dopasowania wynosi: **SP = 51 + 27 + 26 + 2 + 11 + 16 + 0 - 14 - 14 - 20 = 85**


#### A teraz przyjmiemy: gopen = -12, gext = - 1

**Ocena za dopasowanie w1 z w2:** ( A G T C G T A G N P S T oraz A S T C G T A G – P S T)

A-A: +
G-S: +
T-T: +
C-C: +

```
G-G: +
T-T: +
A-A: +
G-G: +
```
```
gopen N i „–” : -12 punktów
P-P: +
S-S: +
T-T: +
```
Czyli suma w1-w2: +4 +0 +5 +9 +6 +5 +4 +6 -12 +7 +4 +5 = **43**

**Dopasowanie w1 z w3:** ( A G T C G T A G N P S T oraz – G T P G – A G N – S T)

gopen A i „–” : -12 punktów
G-G: +
T-T: +
C-P: - 3

```
G-G: +
gopen T i „–” : -12 punktów
A-A: +
G-G: +
```
```
N-N: +
gopen P i „–” : -12 punktów
S-S: +
T-T: +
```
Czyli suma w1-w3: - 12 +6 +5 -3 +6 -12 +4 +6 +6 -12 +4 +5 = **3**

**Dopasowanie w1 z w4:** ( A G T C G T A G N P S T oraz A – N C G T – – N P – T)

A-A: +
gopen G i „–” : -12 punktów
T-N: +
C-C: +

```
G-G: +
T-T: +
gopen A i „–” : -12 punktów
gext G i „–” : -1 punkt
```
```
N-N: +
P-P: +
gopen S i „–” : -12 punktów
T-T: +
```
Czyli suma w1-w4: +4 -12 +0 +9 +6 +5 - 12 - 1 +6 +7 -12 +5 = **5**

**Dopasowanie w1 z w5:** ( A G T C G T A G N P S T oraz – G – – A T A – – P S – )

gopen A i „–” : -12 punktów
G-G: +
gopen T i „–” : -12 punktów
gext C i „–” : -1 punkt

```
G-A: +
T-T: +
A-A: +
gopen G i „–” : -12 punktów
```
```
gext N i „–” : -1 punkt
P-P: +
S-S: +
gopen T i „–” : -12 punktów
```
Czyli suma w1-w5: - 12 +6 - 12 - 1 +0 +5 +4 - 12 - 1 +7 +4 -12 = **- 24**

**Pozostałe pary wierszy:**

- suma **w2-w3** = **- 21**
- suma **w2-w4** = **- 13**
- suma **w2-w5** = **- 29**
    - suma **w3-w4** = **- 59**
    - suma **w3-w5** = **- 48**
    - suma **w4-w5** = **- 73**

**Wszystkie elementy częściowe dopasowania:**

- w1-w2 = 43
- w1-w3 = 3
- w1-w4 = 5
    - w1-w5 = - 24
    - w2-w3 = - 21
    - w2-w4 = - 13
       - w2-w5 = - 29
       - w3-w4 = - 59
       - w3-w5 = - 48
          - w4-w5 = - 73

Zatem: **SP** = 43 + 3 + 5 - 24 - 21 - 13 - 29 - 59 - 48 - 73 = **- 216**. To jest formalnie poprawne dla konwencji
BLAST/NCBI.


Wracamy do DNA:

Spójrzmy jeszcze raz na dopasowanie:

w1 **A G T C G T A G**
w2 **A – T C G T C G**
w3 **– – – – G T A G**
w4 **– G T A G – A G**

Macierz binarna:

1 1 1 1 1 1 1 1
1 0 1 1 1 1 1 1
0 0 0 0 1 1 1 1
0 1 1 1 1 0 1 1

Czyli nasze dopasowanie jest opisane macierzą w której jedynka oznacza literę, a zero oznacza kreskę.
Tylko tyle i aż tyle. Ale uwaga, z tego wynika sposób generowania rozwiązań losowych:

- musimy wygenerować tyle wierszy ile jest sekwencji do dopasowania
- wiersz opisujący daną sekwencję wejściową musi mieć dokładnie tyle jedynek, ile dane
    sekwencja ma liter.

Kolumny w tej macierzy możemy sobie nazywać jak chcemy, np. ruchami. W tym znaczenie ruch to coś,
co prowadzi nas dalej, do konsktrukcji rozwiązania dopuszczalnego. Rozwiązanie dopuszczalne jest
własnie zdefiniowane tymi dwoma powyższymi punktami nad tym akapitem – tyle wierszy ile sekwencji
+ liczba jedynek w wierszu = liczba liter w danej sekwencji (reszta to zera). Należy tutaj zwrócić uwagę
na to, że kolumna składająca się z samych zer nie ma wielkiego sensu – można by ją dodać wiele razy,
ale poza obniżaniem wyniku funkcji celu nie służyłaby ona do niczego pozytywnego...

Należy też zauważyć, że każdy ruch można oceniać (prawie) niezależnie. „Prawie” ponieważ dla
określania gapopen / gapextend w aktualnie analizowanej kolumnie, musimy odwołać się do tego co jest w
danych wiersza w kolumnie bezpośrednio wcześniej, czyli tej po lewej stronie od aktualnej.

Np. pierwsza kolumna to „1 1 0 0”. Jest tam +1 punkt za match nukletydów A z pierwszego i drugiego
wiersza. A także cztery raz po -2 punkty za „gapopen” (wiersz w1 z w3 (-2 pkt.) i i w4 (-2 pkt.), oraz wiersz
w2 z w3 (-2 pkt) i w2 z w4 (-2). Razem: +1 – 2 – 2 – 2 – 2 = -7.

W drugiej kolumnie „1 0 0 1” mamy „gapopen” (w1 vs. w2, w 2 vs. w4 oraz w3 vs. w4 (!), każdy karany - 2
pkt.), dwa razy „gapextend” (w1 z w3 oraz w3 z w4, każde karane -1 pkt), 1 match (w1 z w4, +1 pkt).
Porównanie gap’ów (w2 z w3) jest za zaro punktów. Czyli - 2 - 1 +1 +0 - 2 – 2 = - 6 (w1-w2; w1-w3; w1-w4;
w2-w3, w2-w4, w3-w4). Każdy ruch/kolumnę można więc oceniać niezależnie kolumna po kolumnie,
tak jak na poprzedniej stronie było pokazane, jak wiersz po wierszu liczyliśmy funkcję celu/jakość
dopasowania.


Wracając do przykładu z aminokwasami:

**AGTCGTAGNPST w1: A G T C G T A G N P S T**

**ASTCGTAGPST w2: A S T C G T A G – P S T**

**GTPGAGNST w3: – G T P G – A G N – S T**

**ANCGTNPT w4: A – N C G T – – N P – T**

**GATAPS w5: – G – – A T A – – P S –**

Macierz binarna:

```
1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 0 1 1 1
0 1 1 1 1 0 1 1 1 0 1 1
1 0 1 1 1 1 0 0 1 1 0 1
0 1 0 0 1 1 1 0 0 1 1 0
```
**Wersja 1: prosta kara za gap = - 4**

Dla przykładu aminokwasowego rozważmy kolumny siódmą i ósmą. W kolumnie siódmej mamy układ
**„A A A – A”**. Występuje tu sześć porównań typu **A** – **A** (w1 z w2, w1 z w3, w1 z w5, w2 z w3, w2 z w5,
w3 z w5), każde punktowane **+4** , co daje razem **+24 punktów**. Ponadto są cztery porównania
aminokwasu z luką (w1 z w4, w2 z w4, w3 z w4 oraz w4 z w5), każde po **- 4 punkty** , czyli razem **- 16
punktów**.

#### Suma dla całej kolumny wynosi więc +24 - 16 = +8. W kolumnie ósmej mamy układ „G G G – – ”. Są tu

trzy porównania **G** – **G** (w1 z w2, w1 z w3, w2 z w3), każde po **+6 punktów** , czyli razem **+18 punktów**.
Dodatkowo występuje 6 porównań aminokwasu z luką (w1 z w4, w1 z w5, w2 z w4, w2 z w5, w3 z w4,
w3 z w5), każde po **- 4 punkty** , co daje **- 24 punkty** , a porównanie dwóch luk (w4 z w5) ma wartość **0**.

Ostatecznie dla tej kolumny otrzymujemy **+18 - 24 + 0 = - 6**. W tym uproszczonym modelu każdą kolumnę
można więc oceniać niezależnie, ponieważ kara za przerwę nie zależy od tego, co znajdowało się w
kolumnie poprzedniej.

**Wersja 2: gapopen = -12, gapextend = - 1**

Dla tego samego przykładu aminokwasowego, ale przy punktacji z rozróżnieniem na **gapopen = - 12** oraz
**gapextend = - 1** , sytuacja jest bardziej złożona, ponieważ do oceny danej kolumny trzeba znać także
kolumnę poprzednią. Rozważmy ponownie kolumny siódmą i ósmą. Kolumna siódma ma postać
**„A A A – A”** , natomiast kolumna szósta, bezpośrednio ją poprzedzająca, to **„T T – T T”**.

W kolumnie siódmej znów mamy sześć porównań **A** – **A** , każde po **+4 punkty** , czyli razem **+24 punkty**.
Cztery porównania z luką (w1 z w4, w2 z w4, w3 z w4 oraz w4 z w5) są tutaj traktowane jako **gapopen** ,
ponieważ w kolumnie szóstej nie było dla tych par tej samej, już trwającej luki; w szczególności dla pary
w3 z w4 przerwa „przechodzi na drugą stronę”, więc także liczy się jako nowe otwarcie.

Każde z tych czterech porównań daje więc **- 12 punktów** , łącznie **- 48 punktów** , a zatem suma dla
kolumny siódmej wynosi **+24 - 48 = - 24**. W kolumnie ósmej mamy układ **„G G G – – ”**. Trzy porównania


G-G (w1 z w2, w1 z w3, w2 z w3) dają razem **+18 punktów**. Dla par w1 z w4, w2 z w4 oraz w3 z w4 luka
w w4 trwa już od kolumny siódmej, więc są to trzy przypadki **gapextend** , po **- 1 punkcie** każdy. Natomiast
dla par w1 z w5, w2 z w5 oraz w3 z w5 luka pojawia się dopiero teraz, więc są to trzy przypadki **gapopen** ,
po **- 12 punktów** każdy. Porównanie w4 z w5 to znowu przerwa z przerwą, czyli **0 punktów**. Ostatecznie
dla kolumny ósmej otrzymujemy więc **+18 - 1 - 1 - 1 - 12 - 12 - 12 + 0 = - 21**. Widać zatem, że przy modelu
gapopen/gapextend kolumn nie można oceniać całkowicie niezależnie: trzeba pamiętać, czy w danej parze
sekwencji przerw właśnie się zaczyna, czy jest już tylko przedłużana.

## 3. Metaheurystyki

### 3.1. Algorytm mrowkowy (ACO) w problemie MSA

Algorytm mrowkowy, czyli **Ant Colony Optimization (ACO)** , jest metaheurystyką inspirowaną
zachowaniem prawdziwych mrówek poszukujących pożywienia. W wersji obliczeniowej zamiast
prawdziwych owadów mamy zbiór prostych agentów, nazywanych mrówkami, które budują kolejne
rozwiązania problemu. Każda mrówka konstruuje rozwiązanie krok po kroku, poruszając się po
odpowiednio zdefiniowanym grafie (lub innej strukturze danych opisującej sposoby konstrukcji
rozwiązania). Wybór kolejnego ruchu nie jest całkowicie losowy: zależy z jednej strony od śladu
feromonowego **,** czyli informacji zgromadzonej przez wcześniejsze dobre rozwiązania, a z drugiej strony
lokalnej oceny atrakcyjności danego ruchu. Po zakończeniu iteracji feromon na mniej obiecujących
decyzjach częściowo zanika (parowanie), a na elementach należących do lepszych rozwiązań zostaje
wzmocniony (bo rozwiązania w kolejnych iteracjach zawierają ruchy, które wcześniej były użyte i
zostawiły ślad feromonowy). Dzięki temu algorytm stopniowo uczy się, które decyzje warto preferować.

W problemie **MSA (Multiple Sequence Alignment)** ideę tę trzeba przełożyć na język dopasowania
sekwencji. Tutaj rozwiązaniem nie jest trasa komiwojażera ani ścieżka w klasycznym grafie, lecz pełne
dopasowanie wielu sekwencji **,** czyli układ kolumn zawierających symbole biologiczne oraz ewentualne
luki. Z punktu widzenia ACO trzeba więc odpowiedzieć na trzy pytania:

- po pierwsze, **jak reprezentować rozwiązanie** ,
- po drugie, **jaki ruch wykonuje mrówka** ,
- po trzecie, **na czym odkładać feromon**.

Tutaj pojawia się ważna cecha ACO dla MSA **:** nie istnieje jeden jedyny kanoniczny wariant algorytmu**.**
W literaturze spotyka się kilka sposobów modelowania problemu. Jedne podejścia budują
dopasowanie bezpośrednio jako ścieżkę w wielowymiarowym grafie, który uwzględnia możliwe
wstawienia luk i kolejne decyzje wyrównania. Inne traktują ACO jako metodę udoskonalania już
istniejącego dopasowania **,** na przykład wygenerowanego wcześniej przez heurystykę. Są też warianty
hybrydowe, w których ACO działa jako lokalne przeszukiwanie osadzone w większym algorytmie, na
przykład genetycznym.


Najbardziej intuicyjna interpretacja wydaje się następująca. Wyobraźmy sobie, że mrówka buduje
dopasowanie kolumna po kolumnie. W każdej chwili musi zdecydować, które sekwencje „przesunąć”
dalej o jeden znak, a w których wstawić lukę. Taki wybór odpowiada zbudowaniu następnej kolumny
dopasowania. Dla 𝑘sekwencji każda kolumna może być opisana jako wektor binarny długości 𝑘:
jedynka oznacza pobranie kolejnego symbolu (litery oznaczającej nukleotyd lub aminokwas) z danej
sekwencji, a zero oznacza wstawienie luki. Nie wszystkie wektory są oczywiście dopuszczalne, ale
każdy dopuszczalny wektor można potraktować jako lokalny ruch w przestrzeni rozwiązań. Mrówka
powtarza taki wybór wielokrotnie, aż wszystkie sekwencje zostaną w pełni wykorzystane, a my
otrzymujemy kompletne dopasowanie. To ujęcie bardzo dobrze łączy się z wcześniejszą reprezentacją
macierzy binarnej dla MSA: ACO można po prostu rozumieć jako metodę inteligentnego konstruowania
kolejnych kolumn tej macierzy.

W praktyce decyzja mrówki nie jest czysto losowa. Każdy możliwy ruch ma dwa składniki oceny.
Pierwszy to feromon: jeżeli w poprzednich iteracjach pewne decyzje często występowały w dobrych
dopasowaniach, to ich atrakcyjność rośnie. Drugi składnik to informacja heurystyczna: na przykład
lokalne podobieństwo znaków, przewidywany koszt luki, zgodność z macierzą podstawień albo
poprawa funkcji celu. Właśnie połączenie tych dwóch mechanizmów daje równowagę między
eksploracją a eksploatacją. Eksploracja oznacza sprawdzanie nowych układów kolumn i nowych
miejsc wstawiania luk, natomiast eksploatacja oznacza częstsze wybieranie decyzji, które już
wcześniej doprowadzały do dobrych wyników (czyli większe bazowanie na feromonach). Po każdej
iteracji najlepsze albo wszystkie zbudowane dopasowania wpływają na aktualizację feromonu, zaś
proces parowania feromonu zapobiega temu, aby algorytm zbyt wcześnie „uwierzył” w jedną, być może
tylko lokalnie dobrą, strukturę rozwiązania.

W MSA jakość rozwiązania mierzy się zwykle przy użyciu pewnej funkcji oceny, ewentualnie z karami za
luki i macierzą podstawień dla białek. Po skonstruowaniu dopasowania przez mrówkę można obliczyć
jego wartość i na tej podstawie zdecydować, ile feromonu należy dopisać do decyzji, które do tego
dopasowania doprowadziły. W bardziej zaawansowanych wersjach ACO dla MSA nie buduje się całego
dopasowania od zera, lecz bierze istniejące ustawienie sekwencji i poprawia tylko fragmenty ocenione
jako słabe, na przykład bloki z dużą liczbą niezgodności albo źle ustawione luki. Takie podejście
występuje zarówno w starszych pracach hybrydowych, gdzie ACO działa jako lokalne przeszukiwanie
wewnątrz algorytmu genetycznego, jak i w nowszych metodach, w których mrówki dopracowują już
wygenerowane dopasowanie startowe.

Najważniejsze jest to, że ACO dla MSA nie zagwarantuje optimum globalnego, ale dobrze wpisuje się w
charakter tego problemu. MSA ma ogromną przestrzeń rozwiązań, a liczba możliwych układów luk i
kolumn rośnie bardzo szybko wraz z liczbą sekwencji i ich długością. Dokładne metody szybko stają się
zbyt kosztowne obliczeniowo, dlatego w praktyce stosuje się heurystyki i metaheurystyki. ACO jest tu
atrakcyjne, ponieważ łączy uczenie się z doświadczenia populacji **z** losowością kontrolowaną przez
heurystykę. Innymi słowy, algorytm nie przeszukuje przestrzeni ślepo, ale też nie zamyka się od razu w
jednym wariancie dopasowania. Dobrze nadaje się więc do pokazania studentom, jak biologiczny
problem można zamienić w problem optymalizacyjny, a następnie rozwiązywać go przez iteracyjne
konstruowanie i ocenianie rozwiązań.


W najprostszej wersji opis algorytmu ACO dla MSA można streścić następująco: najpierw definiujemy
sposób reprezentacji dopasowania, następnie wiele mrówek konstruuje własne kandydackie
dopasowania, każdą decyzję podejmując na podstawie feromonu i heurystyki. Po obliczeniu wartości
funkcji celu feromon na elementach należących do lepszych rozwiązań zostaje zwiększony, a część
starego feromonu odparowuje. Proces powtarza się wielokrotnie, aż osiągnięty zostanie limit iteracji
albo przestanie pojawiać się zauważalna poprawa. W rezultacie otrzymujemy dobre przybliżone
dopasowanie wielu sekwencji, choć niekoniecznie najlepsze możliwe.

### 3.2 Algorytm genetyczny

Dla MSA algorytm genetyczny (GA) najlepiej opisywać jako GA z operatorami przystosowanymi do
problemu dopasowania. To ważne, bo zwykłe bitowe mutacje i krzyżowanie bardzo łatwo psują
poprawność dopasowania: zmieniają długości wierszy, zmieniają kolejność znaków w sekwencji albo
tworzą pełne kolumny kresek (gaps). Najprostszy wariant można oprzeć na reprezentacji: **osobnik** =
**pełne dopasowanie MSA** , czyli macierz o 𝑘wierszach i wspólnej długości 𝐿, gdzie każdy wiersz zawiera
symbole sekwencji oraz przerwy. Chromosomem nie jest tu pojedyncza sekwencja, tylko całe
dopasowanie. Funkcja fitness może być dokładnie ta sama jak z przykładów, sum-of-pairs z
odpowiednią punktacją, np. match/mismatch/gap albo macierz BLOSUM + kary za przerwy.

Schemat algorytmu. Najpierw generujemy populację startową, na przykład losowo (ale można jakąś
prostszą heurystyką, bo im lepsze wejście, tym lepszy wynik). Następnie w każdej generacji liczymy
dopasowanie, wybieramy rodziców, tworzymy rozwiązania potomne przez krzyżowanie i mutację,
naprawiamy dopasowanie jeżeli potrzeba, usuwamy pełne kolumny przerw, a potem wybieramy
osobniki do następnej generacji. W praktyce bardzo pomaga strategia elitaryzmu, czyli przepisywanie
kilku najlepszych dopasowań bez zmian do kolejnej populacji.

Najważniejsza część to **krzyżowanie**. W MSA sensowny crossover nie powinien mieszać znaków „po
komórkach”, tylko łączyć fragmenty całych dopasowań. Najprostsza wersja to jednopunktowe
krzyżowanie (dla dopasowań): wybieramy kolumnę cięcia w rodzicu nr 1 a potem w drugim, dobieramy
miejsce tak, aby po sklejeniu lewej części jednego i prawej części drugiego zachować poprawną
kolejność znaków w każdej sekwencji.

Druga sensowna wersja to _uniform/block crossover_ , w którym rozwiązanie potomne składa się z
naprzemiennie wybieranych bloków pochodzących od różnych rodziców, ale tylko między pozycjami
„zgodnymi”, czyli takimi, gdzie połączenie nie łamie spójności sekwencji. Jeszcze prościej: krzyżowanie
= wymiana bloków kolumn między dwoma dopasowaniami, ale z warunkiem, że po tej wymianie każdy
wiersz nadal zawiera symbole swojej sekwencji w tej samej kolejności co przed dopasowaniem. Jeżeli
po sklejeniu długości wierszy się różnią, trzeba dołożyć przerwy wyrównujące; jeżeli powstaną kolumny
z samych przerw, trzeba je usunąć.

**Mutacja** daje algorytmowi możliwość „przesuwania” luk i lokalnego poprawiania dopasowań.
Najprostsza mutacja to _gap insertion_ : wybieramy grupę sekwencji i wstawiamy w nich przerwę tego
samego rozmiaru, a w pozostałych sekwencjach wyrównujemy długość odpowiednim wstawieniem w
innym miejscu. Taki operator generuje nowe konfiguracje przerw których wcześniej nie było.

Druga dość “naturalna” mutacja to tzw. _block shuffling_ albo prościej: przesuwanie bloków przerw lub
bloków znaków w lewo albo w prawo. Idea jest taka, że po wstawieniu przerwy optymalny układ często


różni się tylko tym, że blok powinien być przesunięty o 1–3 kolumny. Podsumowując taki sposób:
mutacja nie zmienia biologicznej treści sekwencji, tylko zmienia położenie przerw lub całych lokalnych
fragmentów dopasowania.

Nadają się tutaj takie mutacje jak na przykład:

- **Shift mutation** : wybieramy jedną przerwę albo blok przerw w danym wierszu i przesuwamy go o
    jedną kolumnę w lewo lub prawo, jeśli nie zmienia to kolejności znaków.
- **Merge/split mutation** : dwa sąsiednie bloki przerw łączymy w jeden, albo odwrotnie: duży blok
    dzielimy na dwa mniejsze.
- **Column cleanup mutation** : usuwamy pełne kolumny przerw (to nawet nie mutacja, tylko
    korekta rozwiązania) lub kolumny prawie puste, a następnie naprawiamy długości wierszy. Taki
    zestaw jest dość praktyczny, bo lokalnie eksploruje przestrzeń rozwiązań bez gwałtownego
    niszczenia dobrych fragmentów dopasowania. Ogólna idea „gap-oriented mutations” jest
    zgodna z kierunkiem rozwoju operatorów w GA dla MSA, gdzie mutacje projektuje się wokół
    przerw i bloków, a nie wokół pojedynczych liter.

Bardzo ważny element techniczny to naprawa rozwiązania po krzyżowaniu (lub mutacji). Po każdym
operatorze sprawdzić trzy rzeczy:

1. Upewnić się, że w każdym wierszu znaki występują w poprawnej kolejności.
2. Usunąć pełne kolumny przerw.
3. Ewentualnie przeliczyć i skompresować dopasowanie do kanonicznej postaci. To nie jest
    kosmetyka tylko warunek sensownego działania GA. Już w reprezentacjach z literatury pojawia
    się założenie, że alignment nie powinien zawierać pełnych kolumn przerw albo że takie kolumny
    są tylko technicznym artefaktem reprezentacji i należy je usuwać.

**Selekcja** rodziców może być zwykła: turniejowa albo rankingowa, albo jeszcze inna, np. ruletka.

GA dla MSA działa na całych dopasowaniach jako osobnikach; krzyżowanie miesza dobre bloki kolumn
z dwóch rodziców, a mutacja lokalnie przesuwa i rekonfiguruje przerwy, po czym najlepsze
dopasowania przechodzą do następnych generacji oceniane selekcją według funkcji celu.


## 4. Specyfikacja projektu/sprawozdanie

**Projekt zaliczeniowy: metaheurystyka dla problemu MSA**

Celem projektu jest zaimplementowanie i przetestowanie wybranej metaheurystyki dla problemu MSA
(Multiple Sequence Alignment) czyli dopasowania wielu sekwencji. Problem ten jest trudny
obliczeniowo, dla typowych funkcji celu takich jak sum-of-pairs wyznaczenie rozwiązania optymalnego
jest problemem NP-trudnym/NP-zupełnym, dlatego w praktyce stosuje się heurystyki i metaheurystyki.

Projekt wykonują zespoły dwuosobowe. Można wybrać algorytm ACO, GA albo inną sensowną
metaheurystykę, ale niezależnie od wyboru trzeba pokazać, że zespół rozumie zarówno sam problem
MSA, jak i działanie zastosowanej metody optymalizacyjnej.

Wariant podstawowy projektu dotyczy sekwencji DNA i pozwala uzyskać ocenę do 4 – 4.5. Wariant
rozszerzony, pozwalający ubiegać się o ocenę 5.0, powinien dotyczyć sekwencji białkowych, z użyciem
sensowniejszej funkcji punktacji, np. BLOSUM62 + kary za luki.

Należy zaimplementować program, który:

- przyjmuje zestaw sekwencji wejściowych,
- tworzy dopasowanie wielu sekwencji,
- ocenia jego jakość za pomocą z góry ustalonej funkcji celu,
- poprawia rozwiązanie za pomocą wybranej metaheurystyki,
- pozwala uruchamiać serię eksperymentów dla różnych parametrów.

Nie oczekuję że algorytm będzie działać bardzo dobrze, najważniejsze jest jednak, aby: działał
poprawnie, był dobrze opisany, dawał się sensownie testować, pokazywał jakąkolwiek poprawę jakości
rozwiązania w zadanym czasie. Sprawozdanie powinno mieć dwie główne części.

**1. Opis problemu i algorytmu**

Ta część powinna mieć około 2 – 3 stron i powinna być napisana tak, abym mógł ocenić, czy grupa
naprawdę rozumie, co implementuje :) Nie chodzi o kopiowanie definicji z internetu, tylko o własny,
techniczny opis. Powinny się tam znaleźć:

- krótkie wyjaśnienie, czym jest MSA i co oznacza jakość dopasowania,
- definicja użytej funkcji celu,
- sposób reprezentacji rozwiązania w programie,
- opis wybranej metaheurystyki,
- opis najważniejszych operatorów lub kroków algorytmu,
- warunek stopu,
- lista najważniejszych parametrów.


Dla ACO interesujące będą zwłaszcza: liczba mrówek, sposób przechowywania feromonu, parowanie,
heurystyka lokalna, reguła wyboru ruchu. Dla GA: reprezentacja chromosomu, sposób tworzenia
populacji początkowej, selekcja, krzyżowanie, mutacja, elityzm, itd. Ta część ma pokazać, że zespół
rozumie zarówno problem biologiczny/informatyczny, jak i mechanikę zastosowanej metaheurystyki.

**2. Część eksperymentalna**

To ma być najważniejsza część sprawozdania (z punktu widzenia oceny). Interesuje mnie nie tylko to, że
program „coś liczy”, ale też czy zespół potrafi przeprowadzić sensowne testy i wyciągać wnioski. W tej
części powinny się znaleźć:

- opis zbiorów testowych,
- opis sposobu pomiaru wyników,
- testy parametrów algorytmu,
- testy na różnych klasach instancji,
- wykresy lub tabele,
- komentarz i wnioski.

Główna miara sprawozdania, czyli ocena, leży w eksperymentach i ich prezentacji. Podstawową miarą
powinna być wartość funkcji celu, czyli np. końcowy SP-score po ustalonym limicie obliczeniowym. To
jest chyba najbardziej naturalne, bo algorytm właśnie tę wartość optymalizuje. Sama końcowa wartość
funkcji celu nie zawsze wystarczy, dlatego dobrze, aby w sprawozdaniu pojawiły się także:

- średni wynik końcowy z wielu uruchomień,
- odchylenie standardowe albo przynajmniej informacja o rozrzucie,
- najlepszy wynik,
- średnia poprawa względem rozwiązania początkowego,
- czas działania lub liczba iteracji / ocen funkcji celu.

Dla algorytmów tego typu pojedynczy eksperyment nie jest wiarygodny. Dlatego każdy punkt na
wykresie lub każdy wpis w tabeli powinien być liczony na podstawie wielu uruchomień, np. 10 albo 20.

Jak mają wyglądać testy parametrów: w sprawozdaniu powinny znaleźć się 2–3 testy parametrów. Nie
trzeba stroić wszystkiego naraz, wystarczy wybrać 2–3 najważniejsze parametry i pokazać, jak wpływają
na wynik. Przykładowo dla **ACO** :

- liczba mrówek,
- współczynnik parowania feromonu,
- relacja wpływu feromonu i heurystyki,
- liczba iteracji.

Dla **GA** :


- liczność populacji,
- prawdopodobieństwo mutacji,
- prawdopodobieństwo krzyżowania,
- poziom elityzmu.

Każdy taki test powinien wyglądać mniej więcej tak:

- wybieramy **jeden parametr** ,
- dla niego sprawdzamy kilka wartości, np. 4–6,
- pozostałe parametry zostają stałe w danym teście,
- uruchamiamy algorytm na **zbiorze testowym** , a nie na jednej instancji (!!!),
- każdy eksperyment powtarzamy wielokrotnie,
- raportujemy średni wynik, ewentualnie rozrzut, i wyciągamy wnioski.

Bardzo ważne: nie należy robić wykresu „parametr vs wynik” na jednej instancji i wyciągać z tego daleko
idących wniosków (o ile jakiekolwiek by można...). To jest za mało wiarygodne. Parametry należy
oceniać na małym, ale zróżnicowanym zbiorze instancji.

Tu nie chodzi o bicie rekordów, tylko o sensowną metodykę. Przy projektach realizowanych na zwykłych
laptopach i w językach typu **C++, Java, C#** dużo lepiej sprawdza się zestaw **krótszych i średnich
instancji** , uruchamianych wiele razy, niż jedna ogromna instancja liczona przez wiele godzin. Jako
praktyczną rekomendację proponuję:

**Wariant dla DNA.** Do strojenia parametrów:

- **małe instancje** : 5–6 sekwencji, długość 20–40,
- **średnie instancje** : 6–8 sekwencji, długość 40–80,
- **większe instancje** : 8–10 sekwencji, długość 80–150.

Dobrze byłoby przygotować:

- po 5–10 instancji w każdej klasie,
- czyli łącznie około 15–30 instancji do eksperymentów.

Przy takim zakresie rozsądny jest czas rzędu:

- 10 – 20 s na małą instancję,
- 20 – 30 s na średnią,
- 30 – 60 s na większą.

To nie są wartości „jedynie słuszne”, tylko rekomendacja. Chodzi o to, aby algorytm miał czas wykonać
pewną optymalizację, ale żeby cały pakiet testów był jeszcze realny do przeprowadzenia.

**Wariant białkowy.** Dla wersji rozszerzonej proponuję:


- 5 – 8 sekwencji,
- długości około 30–100 aminokwasów,
- punktacja typu BLOSUM62 + kara za luki.

**Jak sensownie różnicować klasy instancji:** samo zwiększanie liczby sekwencji i ich długości to trochę
za mało. Dobrze, aby w sprawozdaniu pojawiły się też różne klasy trudności. Polecam dwa proste i
(chyba) sensowne pomysły.

**1. Klasy podobieństwa sekwencji**

Można generować dane z jednej „sekwencji przodka”, a potem wprowadzać mutacje i insercje/delecje
z różnym nasileniem. Dzięki temu dostajemy klasy:

- **łatwe** : sekwencje bardzo podobne,
- **średnie** : umiarkowanie podobne,
- **trudne** : mało podobne.

To jest o tyle dobry pomysł, że wtedy można pokazać, jak jakość działania algorytmu spada wraz ze
wzrostem „odległości ewolucyjnej” między sekwencjami.

**2. Klasy o różnej strukturze luk**

Można utrzymać podobną długość i podobny poziom podobieństwa, ale zmieniać charakter przerw:

- wiele krótkich przerw,
- nieliczne, ale długie bloki przerw.

To też daje bardzo sensowne porównanie, bo dla metaheurystyk przesuwanie i ustawianie długich
bloków przerw bywa trudniejsze niż obsługa krótkich, rozproszonych zmian. Opcjonalnie można
dołożyć trzecią klasę:

- sekwencje z **jednym lub dwoma silnie konserwatywnymi motywami** otoczonymi „szumem”.

Taka klasa jest ciekawa, bo pozwala zobaczyć, czy algorytm potrafi dobrze ustawiać wspólne bloki, a
nie tylko poprawiać wynik globalny przypadkowymi przesunięciami.

**Jak powinny wyglądać wykresy i tabele**

Oczekuję prostych, czytelnych wykresów i tabel. Nie musi być ich dużo, ale powinny być sensowne.

Przykłady:

- wartość funkcji celu w zależności od parametru,
- średnia poprawa względem rozwiązania początkowego,
- jakość końcowa dla różnych klas instancji,
- przebieg zbieżności w czasie lub po iteracjach,
- tabela średni wynik / najlepszy wynik / odchylenie / czas.


Dobrze, jeśli dla każdego wykresu będzie jasno napisane:

- jaki był zbiór instancji,
- ile było uruchomień,
- co dokładnie przedstawia oś Y,
- czy pokazywana jest średnia, mediana czy najlepszy wynik.

**Co będzie miało największy wpływ na ocenę** :

1. **poprawność i kompletność implementacji** ,
2. **zrozumienie problemu i algorytmu** ,
3. **jakość części eksperymentalnej** ,
4. **sensowność wniosków**.

Nie będę oceniał projektu po tym, czy uzyskaliście bardzo wysoki wynik funkcji celu. Znacznie
ważniejsze jest, czy potraficie:

- dobrze zdefiniować problem,
- zaprojektować reprezentację rozwiązania,
- poprawnie zaimplementować metaheurystykę,
- przeprowadzić wiarygodne testy,
- uzasadnić swoje decyzje.

W praktyce:

- na ocenę **3 – 3.5** wystarczy działający program i podstawowy raport,
- na **4 – 4.5** potrzebne będą już sensowne testy parametrów i poprawna analiza wyników,
- na **5.0** oczekuję wersji białkowej albo wyraźnie ambitniejszej metodologii eksperymentów,
    najlepiej z bardzo dobrą analizą oraz starannym porównaniem klas instancji.

Lepiej zrobić mniejszy zakres ale porządnie, niż próbować rozwiązywać bardzo duże instancje bez
sensownej analizy.

