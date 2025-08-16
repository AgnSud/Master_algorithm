# Ścisłe metody numeryczne dla równań różniczkowych zwyczajnych oparte na wielomianach Czebyszewa

Algorytm implementujący iteracyjną metodę rozwiązywania zagadnień brzegowych dla układów dynamicznych z wykorzystaniem rozwinięć Czebyszewa i arytmetyki przedziałowej. Projekt powstał w ramach pracy magisterskiej na Uniwersytecie Jagiellońskim.

## 2. Opis / Streszczenie

Projekt zawiera implementację ścisłej metody numerycznej opartej na rozwinięciach Czebyszewa do rozwiązywania zagadnień brzegowych dla równań różniczkowych zwyczajnych.  
W odróżnieniu od podejścia pracy autorstwa Jan Bouwe van den Berg oraz Ray Sheombarsing 
\textit{Rigorous numerics for ODEs using Chebyshev series and domain decomposition}, algorytm działa iteracyjnie – pojedynczy krok odpowiada przejściu między sekcjami Poincaré i wykorzystuje jedno rozwinięcie Czebyszewa. Dzięki temu metoda pozwala analizować nie tylko orbity okresowe, ale również dowolne trajektorie układów dynamicznych.

Główne etapy algorytmu:
1. Wyznaczenie przybliżonego rozwiązania równania $F(x^*) = 0$ metodą Newtona.
2. Dobór parametrów: stopnia rozwinięcia Czebyszewa $N$ oraz wagi $\nu$.
3. Obliczenie ograniczeń $Y$ i $Z$ w arytmetyce przedziałowej.
4. Konstrukcja wielomianów radii $p_j$ i wyznaczenie przedziału promieni $r$, w którym istnieje dokładne rozwiązanie.
5. Potwierdzenie istnienia rozwiązania w kuli o promieniu $r$ oraz wyznaczenie nowego punktu początkowego.
6. Powtórzenie procedury dla kolejnych kroków trajektorii aż do osiągnięcia zadanego czasu całkowitego $r_t$.

Program został przetestowany na klasycznym **układzie Lorenza**, znanym z chaotycznego zachowania i wrażliwości na warunki początkowe.  
Porównano trajektorie uzyskane metodą Czebyszewa z referencyjnymi wynikami otrzymanymi metodą Taylora (CAPD). Eksperymenty wykazały:
- wysoką zgodność trajektorii obu metod,
- wpływ stopnia wielomianu Czebyszewa na dokładność oraz stabilność obliczeń,
- możliwość kontrolowania długości kroków czasowych przez dobór parametrów $N$, $\nu$ oraz $\Delta t$,
- istnienie optymalnych ustawień, które minimalizują promień kuli zawierającej rozwiązanie przy rozsądnym czasie obliczeń.

Dzięki zastosowaniu arytmetyki przedziałowej algorytm nie tylko przybliża rozwiązanie, ale dostarcza także **ściśle zweryfikowanych przedziałów** zawierających rozwiązanie dokładne.
