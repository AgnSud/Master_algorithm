# Ścisłe metody numeryczne dla równań różniczkowych zwyczajnych oparte na wielomianach Czebyszewa

Algorytm implementujący iteracyjną metodę rozwiązywania zagadnień brzegowych dla układów dynamicznych z wykorzystaniem rozwinięć Czebyszewa i arytmetyki przedziałowej. Projekt powstał w ramach pracy magisterskiej na Uniwersytecie Jagiellońskim.

## 2. Opis / Streszczenie

Projekt zawiera implementację ścisłej metody numerycznej opartej na rozwinięciach Czebyszewa do rozwiązywania zagadnień brzegowych dla równań różniczkowych zwyczajnych.  
W odróżnieniu od podejścia pracy autorstwa Jan Bouwe van den Berg oraz Ray Sheombarsing 
_Rigorous numerics for ODEs using Chebyshev series and domain decomposition_, algorytm działa iteracyjnie – pojedynczy krok odpowiada przejściu między sekcjami Poincaré i wykorzystuje jedno rozwinięcie Czebyszewa. Dzięki temu metoda pozwala analizować nie tylko orbity okresowe, ale również dowolne trajektorie układów dynamicznych.

Główne etapy algorytmu:
1. Wyznaczenie przybliżonego rozwiązania równania $F(x^*) = 0$ metodą Newtona.
2. Dobór parametrów: stopnia rozwinięcia Czebyszewa $N$ oraz kroku czasowego $\Delta t$.
3. Obliczenie ograniczeń $Y$ i $Z$ w arytmetyce przedziałowej.
4. Konstrukcja wielomianów radii $p_j$
5. Wyznaczenie przedziału promieni $r$, w którym istnieje dokładne rozwiązanie.
6. Powtórzenie procedury dla kolejnych kroków trajektorii, przyjmując przedział końcowy
z poprzedniego kroku jako przedział początkowy, aż do osiągnięcia zadanego czasu całkowitego $r_t$.

Program został przetestowany na klasycznym **układzie Lorenza**.  
Porównano trajektorie uzyskane metodą Czebyszewa z referencyjnymi wynikami otrzymanymi metodą Taylora (CAPD). Eksperymenty wykazały:
- wysoką zgodność trajektorii obu metod,
- wpływ stopnia wielomianu Czebyszewa na dokładność oraz stabilność obliczeń,
- możliwość kontrolowania długości kroków czasowych przez dobór parametrów $N$, $\nu$ oraz $\Delta t$,
- istnienie optymalnych ustawień, które minimalizują promień kuli zawierającej rozwiązanie przy rozsądnym czasie obliczeń.

Dzięki zastosowaniu arytmetyki przedziałowej algorytm nie tylko przybliża rozwiązanie, ale dostarcza także **ściśle zweryfikowanych przedziałów** zawierających rozwiązanie dokładne.

## 3. Założenia i kontekst teoretyczny

Algorytm opiera się na połączeniu **arytmetyki przedziałowej** oraz rozwinięć w **wielomianach Czebyszewa**.

### Arytmetyka przedziałowa
W obliczeniach numerycznych pojawia się problem ograniczonej precyzji reprezentacji liczb zmiennoprzecinkowych.  
Aby zapewnić kontrolę błędów i uzyskać wyniki ściśle zweryfikowane, zastosowano arytmetykę przedziałową, 
w której niepewne wartości reprezentowane są w postaci przedziałów ```Interval``` z CAPD.

Dzięki temu można zagwarantować, że dokładne rozwiązanie znajduje się wewnątrz otrzymanego przedziału.

### Wielomiany Czebyszewa
Podstawą algorytmu jest rozwinięcie funkcji w szereg Czebyszewa:
\[ 
y = a_0 + 2 \sum_{k=1}^{N} a_k T_k
\]

gdzie `T_k` oznacza wielomian Czebyszewa stopnia *k*.  

Do obliczania wartości szeregu Czebyszewa wykorzystano algorytm **Clenshawa**, który jest wydajną i stabilną metodą sumowania.
Używając $(a_k)_{k=0}^N$ skontruowano operator Czebyszewa $F(\omega, a)$, którego zero szukane jest w algorytmie. 
Jego znalezienie zapewni rozwiązanie zadanego problemu brzegowego, zgodnie z _Twierdzeniem 4.1_ z pracy. 
Pełne omówienie teorii zawarte jest w Rozdziale 3 i 4 pracy magisterskiej.

## 4. Wymagania systemowe

### System operacyjny
- Linux (Ubuntu 20.04 lub nowszy testowany)  
  *(możliwe uruchomienie na innych systemach po odpowiedniej konfiguracji CAPD)*

### Kompilator
- g++ z obsługą **C++20**

### Narzędzia budowania
- **CMake** w wersji co najmniej 3.24

### Zewnętrzne biblioteki
Projekt korzysta z biblioteki **CAPD (Computer Assisted Proofs in Dynamics)**:
- `capdDynSys`
- `capdAlg`
- `capdAux`
- `capdExt`

Ścieżki do bibliotek muszą być skonfigurowane w systemie przed kompilacją (wskazanie folderów `include` i `lib` w CMakeLists.txt).

## 5. Instrukcja instalacji

### Krok 1. Pobranie repozytorium
Sklonuj repozytorium z GitHuba:
```bash
git clone https://github.com/AgnSud/Master_algorithm/tree/master
```
### Krok 2. Kompilacja projektu

W katalogu głównym repozytorium utwórz folder `build` i skompiluj projekt:
```bash
mkdir build
cd build
cmake ..
make
```


Po pomyślnej kompilacji plik wykonywalny `master_algorithm` znajdzie się w katalogu `build`.

## 6. Struktura katalogów i plików

```bash
master_algorithm/
├─ source/
│ └─ ChebyshevSeries/
│ ├─ ChebyshevSeries.hpp # Klasa reprezentujaca rozwiniecie Czebyszewa 
│ ├─ ChebyshevSeries.tpp 
│ ├─ Norm.hpp # Normy ważone w $\ella_n^1$ i związane operacje 
│ ├─ norm.tpp 
│ ├─ ChebyshevOperatorFinite.hpp # Operator F_N (skończony wymiar operatora Czebyszewa): f0, f1, DF, A itp.
│ ├─ ChebyshevOperatorFinite.tpp
│ ├─ RadiiPolynomials.hpp # Konstrukcja i znalezienie szukanego promienia dla wielomianów radii
│ ├─ RadiiPolynomials.tpp
├─ helpers.hpp # Funkcje pomocnicze: generacja multiindeksów, defineFunctionG,
│ # konwersje, checkSolution, I/O itp.
├─ testing.hpp # Dodatkowe narzędzia testowe / walidacyjne
├─ main.cpp # Główny program:
│ # - interpolacja Taylorem i Czebyszewem na [a,b]
│ # - Newton w skończonym wymiarze $(\omega, a)$
│ # - przejście do arytmetyki przedziałowej (CAPD Interval)
│ # - liczenie Y/Z, wielomianów radii i wyznaczanie promienia r
│ # - dalsze kroki iteracyjnie do kolejnych sekcji
│ # - zapis wyników do CSV
├─ CMakeLists.txt # Konfiguracja budowania (C++20, CAPD, filib++)
├─ all_plots.ipynb # Plik w którym są wszystkie rysujące wykresy 
└─ README.md # Dokumentacja projektu
```


### Pliki wynikowe (artefakty programu)

Podczas działania program zapisuje wyniki w plikach `.csv` w katalogu uruchomienia:

- **`trajectoryN_<N>_step_<k>_base_time_step_<Δt>_nu_<ν>_attempt_nr_<id>.csv`**  
  Porównanie przebiegów: czas (Chebyshev/Taylor), punkty `(x,y,z)` dla obu metod oraz norma różnicy.

- **`u0_interval_and_midpoint_N<N>_base_time_step<Δt>_global_step_<K>_nu_<ν>_attempt_nr_<id>.csv`**  
  Przedziały (lewy/prawy) i środki kolejnych punktów startowych `u0` po weryfikacji promieniem `r`, wraz z faktycznym `dt` per krok.

- **`radii_polynomial_N<N>_base_time_step<Δt>_global_step_<K>_nu_<ν>_attempt_nr_<id>.csv`**  
  Promienie `r` wyznaczone w kolejnych iteracjach.

> ℹ️ **Uwaga**: katalogi `cmake-build-*` to foldery robocze generowane przez CMake/IDE i zwykle powinny być dodane do `.gitignore`.
