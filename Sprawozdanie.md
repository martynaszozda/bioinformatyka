# Raport z Projektu: Dopasowanie Wielu Sekwencji (MSA) z wykorzystaniem Algorytmu Mrówkowego (ACO)

## 1. Opis problemu i zastosowanego algorytmu

### 1.1. Wprowadzenie do problemu Multiple Sequence Alignment (MSA)
Dopasowywanie wielu sekwencji (MSA) stanowi jedno z fundamentalnych zagadnień obliczeniowych we współczesnej bioinformatyce. Proces ten polega na jednoczesnym zestawieniu co najmniej trzech biologicznych sekwencji (w tym układzie peptydowych/aminokwasowych), w sposób maksymalizujący wyłonienie ewolucyjnie konserwowanych obszarów, wspólnych motywów strukturalnych oraz relacji filogenetycznych. Cel ten osiągany jest przez przestrzenne rozmieszczenie liter oraz wprowadzanie znaków przerwy (kodowanych jako luki lub *gaps*), co biologicznie odwzorowuje mutacje o charakterze insercji i delecji, nagromadzonych historycznie w badanych liniach populacyjnych.

Dlaczego w tym projekcie zdecydowano się na operowanie danymi białkowymi zamiast bazowaniem na prostych łańcuchach kwasów nukleinowych (DNA)? Kwasy nukleinowe operują zaledwie na czteroliterowym alfabecie (A, C, G, T), a znaczna część punktowych mutacji genowych na tym poziomie jest w istocie mutacjami synonimicznymi (z racji na degenerację kodu genetycznego, czyli zjawiska determinującego białka przez liczne tripletów-nie-skorelowane). Przejście na badanie złożonych białek, które składają się z 20 zróżnicowanych biochemicznie aminokwasów, stwarza model o wiele dojrzalszy z punktu widzenia wnioskowań. Zastąpienie aminokwasu hydrofobowego polarmościowym ma bowiem niebagatelny i wyraźny wpływ na zróżnicowanie badanych osobników na etapie całego ustroju.

### 1.2. Definicja funkcji oceny dopasowania
Do determinacji globalnej jakości wyrównania wykorzystano klasyczną funkcję sumowania par (ang. *Sum-of-Pairs*, SP). Opiera się ona na niezależnym ewaluowaniu kosztów sumy dopasowań pojedynczych par wyizolowanych z całości zbioru w obrębie wszystkich sekwencji. Ewaluacja każdego dopasowania obłożona jest punktacją podyktowaną przez empiryczną macierz probabilistyczną podstawień **BLOSUM62** oraz model kar afinicznych dla badanych przerw (gdzie koszt inauguracji nowej formacji wyniósł *-12*, zaś obciążalność za ewentualną kontynuację podjętej luki ujęto przy *-1* parametrze rozszerzającym gapextend).

Przyjęcie takiej kombinacji punktowej (tj. SP z BLOSUM62 i karami afinicznymi) jest obiektywnie ugruntowane w standardach nowoczesnej biologii molekularnej oraz pakietu ujęć statystycznych BLAST. Macierz uwzględnia zmienność fenotypową organizmów. Zastosowanie obniżonej kary za długą i litą lukę zamiast szatkowania z perspektywy stochastycznej, lepiej charakteryzuje nienaturalne wyłomy translacyjne. Powstanie ewolucyjne pojedynczego braku w matrycy dziesięciu aminokwasów jest zjawiskiem znacznie częstszym od dziesięciu izolowanych mutacyjnych utrat osadzonych z rzadka.

### 1.3. Reprezentacja i dyskretyzacja przestrzeni rozwiązania
Postawiony zbiór trudności obliczeniowej umodelowano jako iteracyjny proces kreacji i składania matrycy dopasowującej – w sekwencji kolumna po kolumnie. Pojedynczym ruchem algorytmicznym nazwano dyskretny wybór wektora binarnego nałożonego na aktualną rozdzielczość. Odczyt „1” dekoduje zgodę na „inkorporowanie” kolejnego fizycznego znaku z danego łańcucha białkowego do wspólnego pola odczytu wyrównania, rzutując odwrotnie znak „0” na operację spowolnienia tej inicjatywy, poprzez wpisanie znaku luki („-”).

Implementacja binarnych mask z restrykcją dla wyłączonych zer (wektorów samych luk) posiada fundamentalną zaletę dla poprawności konstruktu strukturalnego. Deterministycznie zabezpiecza on architekturę sekwencji przed naruszeniem (brakiem, wycięciem lub podwójnym zaliczeniem symbolu wejściowego z bazy), nierzadko trapiącym operatorów ucięcia stosowanych np. w krzyżowaniu w algorytmach ewolucyjnych lub genetycznych (GA). Traktując to ujęcie ze ściśle matematycznej precyzji: graficzna dyskretyzacja przekształca optymalizację rozkładu z problemu stochastycznego tasowania do odnalezienia racjonalnej odległości i ścieżki w multigrafie tranzycji.

### 1.4. Metaheurystyka: Algorytm mrówkowy (Ant Colony Optimization - ACO)
Ant Colony Optimization jest wyrafinowanym przedstawieniem koncepcji inteligencji roju i algorytmów populacyjnych. Swoje fundamenty czerpie w badaniu etologii kolonii rzeczywistych owadów i pozostawianych śladach infochemicznych (feromonowych). Wymierny zapach szlaku implikuje zagęszczenie poruszających się jednostek, potęgując sprzężenie zwrotne dla obiecującej odnogi systemu.

Wytypowanie tego narzędzia dla wyznaczania kompromisu w roztrzyganiu wyrównania MSA opiera się o twardy szacunek jego złożoności. Ze względu na fakt, że dogłębne wytyczanie najmniejszej ścieżki algorytmem eksaktnym w zagadnieniu o profilu NP-trudnym pożarłoby ogrom procesorów (oraz czasu) badaczy w relacji przestrzeni kombinatorycznej. Zastosowana metaheurystyka proponuje w tej domenie nienarzucający, harmonijny i niezwykle płynny bilans między celową eksploatacją wiedzy zgromadzonej z iteracji powtórzeń uprzednich, a zuchwałą eksploracją i losowym odchyleniem bazującym na lokalnej ocenie stymulowanej odgórnym parametrem ścierania starszych śladów z gałęzi nieposiadających perspektyw ulepszających. Prawdopodobieństwo optymalnego i sprawnego odnalezienia sub-optymalnych tras silnie zyskuje pod ujęciem kolonii względem technik w pełni ślepego losowego dobierania lub naiwnych struktur heurystyk zachłannych.

### 1.5. Architektura głównych interakcji i elita stada
Fundament selekcyjny przy komponowaniu kolejnego wektora kolumnowego poddany jest rachunkom opartym na rozdziale zmiennych: stężeniu lokalnego feromonu (pamięci ogólnej) $\tau^{\alpha}$ pomnożonym przez pre-zdefiniowane oddziaływanie wskaźnika oceniającego bieżący koszt, wyrażonej heurystyką lokalną $\eta^{\beta}$. Przesłanką heurystyczną jest uproszczona, a zarazem dynamiczna kalkulacja częściowa na funkcji SP nowej kolumny – by przed ukończeniem pętli, mrówka już operowała intuicją wykluczania fatalnej pod kątem kar konfiguracji.

Na szczególną obłokację analityczną w systemie narzuconym środowisku algorytmu zwraca jednakże system "elitaryzmu" w dyslokacji i modyfikacji ścieżek ($\Delta \tau$). Dopuszczalną partycypację we wlewaniu świeżego czynnika feromonowego obarczono ekskluzywnym blokiem najlepszej frakcji pracującej ułamka liczebności algorytmicznego cyklu operacyjnego (zwykle operującej w zakresie $10\%-25\%$ iteracji górnego kwartyla w zależności i na miarę populacji). Implementacja owego restrykcyjnego progu diametralnie rzutuje na redukcję jałowych szmerów informacyjnych wprowadzanych po węzłach mrówek fatalnych. Konwergencja do globalnego atraktora minimum osiągana pod elitarną frakcją roju przebiega drastycznie szybciej i gwarancyjnie bez gwałtownej zatraty eksploracyjnej wynikłej z powszechnej obecności $\rho$ odparowywań na zapomnianych węzłach.

### 1.6. Parametryzacyjna konstrukcja oraz warunek stopu
Osią krytyczną i nieprzerywalną w narzuconym algorytmie operującym programem ustanowiono klasyczny dyktat iteracyjny pętli mrówek pojęciowy, narzucający ściśle graniczny pułap powtórzeń działania optymalizacyjnego po populacji roju. Ten zbieg omijający kalkulację zbieżności w systemach asygnowanych względem limitów i okien poprawności gwarantuje przejrzyste estymowanie wydajności i kosztów algorytmu połączonych z realnym wymiarem i wolumenem środowiskowo-pamięciowym testowanej struktury. 
Ostateczny werdykt jakości opiera się m.in na:
- **num_ants**: wolumen mrówek w populacji determinujący rzut zrównolegleń losowych poszukiwań bazy
- **iterations**: graniczne sprzężenie wyżej przedstawionego stopu
- $\alpha$, $\beta$: sprzężone modyfikatory relacji priorytetów sił ufności stada a samodzielnej analizy fizyki środowiska
- $\rho$: wskaźnik wymierania i ścierania feromonowego dający algorytmowi wyzwolenia od skamieniałych minimów regionalnych po starych ścieżkach węzłów.

Odrzucono koncepcje "braku poprawy przez $N$ populacji" zważając, iż w ujęciu losowości natury rzędu mrówkowego i złożoności samego oceniania wyliczalnego funkcji SP przy kilkunastu polach przerw o wymiarze dziesiątek aminokwasów, prewidencja gwarancyjna iteracji pozwala optymalniej dystrybuować eksperymenty w ograniczonym przedziale czasowym badacza dla maszyn obliczeniowych osobistego rzędu. 

---

## 2. Metodyka Badawcza i Eksperymentalna

### 2.1. Klasyfikacja oraz kompozycja zbiorów instancji testowych
Walor sprawdzający ujęto obiektywnym we współczesnej metodologii badawczej tworzeniem pakietowanej grupy wirtualnych homologów. Autorskie wykreowanie łańcucha narzuca system sekwencji ewolucyjnej, od predefiniowanego sztucznie szczepu (wspólnego prapoczątku referencyjnego peptydu) wymuszającego z odpowiednim kwartylem wariancji biologicznej szumy mutacyjne (wywołujące translokacyjne substytucje oraz pożądane Indels).

Rozpiętość ujęto podziałem dwoistym dla rozmiaru populacji (grupa testowa *Small* oparta na wektorze 5 nici do 40 rzędowości oraz obciążalnością grupy *Medium* posiadającej 7 nici objętych po 70 cząsteczek wektora bazy). Najbardziej fundamentalne z punktu widzenia wnioskowań poznawczych założenia naświetla stopień separacji homologicznej. Generacja posiłkuje podział w przedziałach: 
- `High similarity` – ewolucyjne gałęzie bliźniacze, o wciąż znacznej spójności,
- `Medium similarity` – przeciętnie dotknięte dryftem mutacyjnym łańcuchy,
- `Low similarity` – zaświadczające o dalekim filogenetycznym rozejściu struktury ewolucyjnie odrębnych szczepów. Badanie tego typu szczepów dla zjawiska algorytmu optymalizacyjnego i heurystyk uwydania sprawność mechanizmów stochastycznych roju do znajdowania zbieżności między głęboką erozją danych a prawidłowym wykryciem odległych, wspólnych konserwacji biologicznych po rekolonizacji wielkimi modelami kar z tytułu luk i przesunięć.

### 2.2. Protokół uśrednień oraz optymalizacja miar eksperymentów
Trafność badań ewaluowana pomiarami i miarami oparta została stricte o ocenę jakości funkcji finalnego kosztu SP (SP-score) narzuconą wyznacznikami i statusem czasochłonności wyrażonym zrzutami interwału operacyjnego wyrażonego metryką chronologiczną (w sekundach układów roboczych i optymalizacyjnych środowiska interpretera w systemie wykonawczym). Prezentowany reżim nakłada konfrontację parametrów decyzyjnych algorytmu – dyslokacja i siła kolonii `num_ants` a potęga gęstości parującej $\rho$ z jednoczesną retencją wymogów dla środowiska stałych wymiarów algorytmicznych instancji kontrolnych by dowieźć sprawczości, bądź obalać mity uśrednionych popularyzacyjnie pre-definicji.


*(TU PÓŹNIEJ ZOSTANĄ DOŁĄCZONE WŁAŚCIWE ZAŁĄCZNIKI (TABELE/WYKRESY) BADANE I WYDEDUKOWANE NA BAZIE POWTARZANYCH EKSPERYMENTÓW I UŚREDNIEŃ KORELACJI URUCHAMIANYCH SKRYPTEM BADACZY W REŻIMIE UŚREDNIAJĄCYM LUB POMIARACH PROBABILISTYCZNYCH ODCHYLEŃ O KILKUKROTNYM ZRZUCIE GENERATORA DANYCH DO PLIKÓW ARCHIWUM REJESTRACYJNEGO RESULT)*.
