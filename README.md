# Dynamic_Modelling
Package voor het modelleren, simuleren en fitten van tumorgroei met verschillende klassieke groeimodellen.

## Installatie en quick start

### 1. Repository clonen

```bash
git clone https://github.com/jouwnaam/Dynamic_Modelling.git
cd Dynamic_Modelling
```

```python
from tumor_models import Solver, Searcher
from tumor_models import GompertzLesModel, LogisticGrowthModel
```
# Biologische achtergrond van tumorgroei

Tumorgroei ontstaat wanneer cellen ongecontroleerd gaan delen door genetische en epigenetische veranderingen. Een tumor is niet simpelweg een massa van delende cellen, maar een dynamisch systeem waarin de balans tussen celdeling, celdood, aanvoer van voedingsstoffen en interacties met de micro-omgeving voortdurend verandert.

In de vroege fase van tumorgroei is meestal voldoende zuurstof en ruimte aanwezig, waardoor de groei relatief snel verloopt. Naarmate de tumor groter wordt, ontstaat er een gebrek aan nutriënten en zuurstof doordat diffusie via het oppervlak onvoldoende wordt. Dit leidt tot afname van de netto groeisnelheid, necrose in de kern en uiteindelijk de behoefte aan angiogenese (vorming van nieuwe bloedvaten). Hierdoor vertonen tumoren typisch drie groeifasen:

Initiële exponentiële fase — snelle celdeling, weinig beperking.

Sub-exponentiële fase — groei wordt vertraagd door ruimte- en nutrientlimitaties.

Plateau fase — netto groei ≈ 0 doordat celdeling en celdood in evenwicht komen.

Biologische processen zoals hypoxie, mechanische stress, immunologische interacties en therapierespons versterken deze dynamiek. Tumorgroeimodellen proberen deze biologische principes wiskundig te beschrijven. Simpele modellen zoals exponentiële groei passen bij vroege tumorontwikkeling, terwijl modellen als Gompertz, logistiek of Von Bertalanffy beter aansluiten bij tumoren waar groei wordt afgeremd door biologische beperkingen.
# Moddelen 
## Lineaire groei
[Linaire groei](https://www.wiskundeacademie.nl/lesmethode/getal-en-ruimte/havo-a/10-groei/10-1-lineaire-en-exponentiele-groei) is een groei die per tijdseenheid dezelfde helling heeft woordoor de groei consant is. 
### Winkunde model
$$
\frac{\text{d}V}{\text{d}t} = c
$$s
waarbij

- V(t) het volume (of populatie, tumorvolume, etc.) op tijd t voorstelt,
- c een constante groeisnelheid is (positief voor groei, negatief voor krimp).

## Exponentieel toenemende groei

Bij [exponentiële groei](https://www.wiskundeacademie.nl/lesmethode/getal-en-ruimte/havo-a/10-groei/10-1-lineaire-en-exponentiele-groei) neemt de hoeveelheid niet met een constant bedrag toe, maar met een vast percentage per tijdseenheid.
Hoe groter het volume
V wordt, hoe sneller de groei plaatsvindt. Dit maakt exponentiële groei veel sneller dan lineaire groei.
### Wiskundig model
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V
$$
waarbij
- V(t) het volume op tijd t is,
- c de groeiconstante is (positief voor groei, negatief voor afname).

## Mendelsohn-groei

[Het Mendelsohn-model](https://aacrjournals.org/cancerres/article/73/8/2407/591536/The-Model-Muddle-In-Search-of-Tumor-Growth-LawsThe) beschrijft een groeiproces waarbij de groeisnelheid afhankelijk is van een machtsfunctie van het volume.
Het is een gegeneraliseerde vorm van sub-exponentiële groei: de tumor groeit sneller dan lineair, maar langzamer dan exponentieel (voor 
0<d<1
0<d<1).

### Wiskundig model

De algemene differentiaalvergelijking is:
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V^d
$$
waarbij

- V(t): het tumorvolume op tijd t,
- c: groeiconstante,
- d: machts-exponent (typisch 
0<d<1
0<d<1).
## Exponentieel afvlakkende groei
[Exponentieel afvlakkende groei](https://nl.wikipedia.org/wiki/Exponenti%C3%ABle_groei) beschrijft een proces waarbij de groeisnelheid afneemt naarmate het volume toeneemt.
Naarmate V(t) dichter bij de maximale waarde 
Vmax komt, wordt de netto groei kleiner. Dit model wordt veel gebruikt wanneer er een natuurlijke limiet is aan de groei (bijv. door beperkte middelen of omgevingsfactoren).
### Wiskundig model
$$
\frac{dV}{dt} = c \cdot \left( V_{\text{max}} - V \right)
$$ 
waarbij
- V(t): het volume op tijd t
- c: groeiconstante,
- Vmax	​: de maximale omvang (carrying capacity)
## Logistische groei 
[Logistische groei](https://www.hhofstede.nl/modules/logistischegroei.htm) beschrijft een proces waarbij een populatie (of tumorvolume) eerst snel groeit, maar waarvan de groei later afremt doordat er een maximale capaciteit bestaat.
### Wiskundig model

$$
\frac{dV}{dt} = c \cdot V \cdot \left( V_{\text{max}} - V \right)
$$

Hierbij:
- V(t): het volume op tijd t
- c: groeiconstante
- Vmax: maximale waarde (carrying capacity / limiet)
Interpretatie van de termen
- V: maakt de groei proportioneel aan de huidige hoeveelheid
(Vmax−V): remt de groei zodraV dichter bij Vmax komt

## Montroll groei 
[Dit model](https://www.ime.unicamp.br/~biomat/Bio26_art10.pdf) is een gegeneraliseerde vorm van logistische groei waarin zowel het volume als de maximale capaciteit worden opgehoogd tot een exponent d. Hierdoor ontstaat een flexibele S-curve waarvan de vorm wordt bepaald door de waarde vand. Het model beschrijft groei die eerst versnelt en later vertraagt naarmate V de verzadigingswaarde 
Vmax nadert, maar met een afvlakking die afhankelijk is van de gekozen exponent.
## Wiakundig model 
$$
\frac{dV}{dt}
 = c \cdot V \cdot \left( V_{\text{max}}^{d} - V^{d} \right)
$$

## Allee-effect (Allee Effects)

Het Allee-effect beschrijft dat individuen in een populatie minder goed overleven, groeien of voortplanten wanneer de populatie te klein is. Dit betekent dat bij lage dichtheid de “fitness” van individuen afneemt.
[Nature Education / Scitable – Allee Effects](https://www.nature.com/scitable/knowledge/library/allee-effects-19699394/)

Normaal verwacht je dat een kleine populatie sneller groeit door weinig competitie, maar het Allee-effect laat zien dat te weinig individuen juist problemen veroorzaakt. Voorbeelden hiervan zijn:

- moeite om partners te vinden 

- verlies van groepsbescherming

- te weinig sociale interacties

- lagere genetische diversiteit

Hierdoor kan de populatie bij lage dichtheid langzaam groeien of zelfs krimpen.

Er bestaan twee vormen:

🔹 ****Zwak Allee-effect****

De groei wordt trager, maar blijft nog steeds positief.
De populatie kan overleven, maar moeilijker.

🔹 ****Sterk Allee-effect****

Er bestaat een kritische drempel (Allee-drempel).
Als de populatie onder die drempel komt, is de groei negatief -> uitsterven wordt heel waarschijnlijk.

### Wiskundig model (symmetrische vorm)

Een veelgebruikte vereenvoudigde vorm van het Allee-effect is:
$$
\frac{dV}{dt} = c \cdot (V - V_{\text{min}}) \cdot (V_{\text{max}} - V)
$$
### Parameters
- Vmin — minimale drempel (Allee-drempel).
Onder deze waarde is de groei negatief.

- Vmax — maximale populatie / draagkracht (carrying capacity).
De groei stopt wanneer V = Vmax.

- c — groeiconstante.
Bepaalt hoe snel de populatie groeit of krimpt.

- V — populatiegrootte of tumorvolume.
De waarde die over de tijd verandert.

## Lineair gelimiteerde groei
Bij [lineair gelimiteerde groei](https://pmc.ncbi.nlm.nih.gov/articles/PMC6813171/) neemt de groeisnelheid toe met het volume, maar wordt deze direct beperkt door een saturatieterm in de noemer. Hierdoor groeit het systeem snel bij kleine waarden van 
V maar vlakt de groei geleidelijk af wanneer V groter wordt. Het model beschrijft dus een vorm van mild begrensde groei, zonder een expliciete maximumcapaciteit zoals bij logistieke groei.
### Wiskundig model
$$
\frac{\text{d}V}{\text{d}t} = c \cdot \frac{V}{V + d}
$$
De parameter 
d bepaalt hoe snel de afremming inzet: een groter d betekent dat de groei sterker wordt gedempt. Het model wordt soms gebruikt als eenvoudige rationele saturatie­functie wanneer een beperkt remmend effect gewenst is
## Oppervlakte-gelimiteerde groei

[Dit model](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12450) gaat ervan uit dat groei afhankelijk is van het oppervlak van de tumor. Naarmate de tumor groter wordt, neemt de toevoer van zuurstof en nutriënten relatief af, waardoor de groei afremt.

### Wiskundig model
$$
\frac{\text{d}V}{\text{d}t} = c \cdot \frac{V}{\left( V + d \right)^\frac{1}{3}}
$$
Kernidee: groei vertraagt doordat het oppervlak (∝ V2/3V2/3) relatief kleiner wordt ten opzichte van het volume.

## Von Bertalanffy groei

[Een klassiek biologisch model](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12450) waarbij groei evenredig is aan oppervlak en verlies evenredig aan volume.

### Wiskundig model
kernidee:

- cV2/3: aanvoer via oppervlak (groei)

- dV: metabolisch verlies (krimp)-resulteert in een geleidelijke afvlakking naar een evenwichtsvolume.
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V^\frac{2}{3} - d \cdot V
$$
## Gompertz groei

Een van de meest gebruikte [modellen(https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12450)] voor tumorgroei. Beschrijft snelle initiële groei, gevolgd door een steeds sterker wordende logaritmische remming.

### Wiskundig model
$$
\frac{dV}{dt} = c \cdot V \cdot \ln\left(\frac{V_{\text{max}}}{V}\right)
$$

Betekenis van de termen

- c – groeiconstante
- V – huidig tumorvolume
- Vmax – maximale waarde / verzadigingsniveau

### Interpretatie
De groeisnelheid is evenredig met het huidige volume V.
De factor
ln(VVmax) wordt kleiner naarmate V groter wordt.

Hierdoor daalt de groei exponentieel en nadert V langzaam 
Vmax⁡.

# Fitting van tumorgroeimodellen aan data

Om te bepalen welk groeimodel het best aansluit bij experimentele of klinische data, wordt elk model numeriek gesimuleerd en vervolgens vergeleken met de gemeten waarden. Dit gebeurt in meerdere stappen:

## Numerieke oplossing van het model

Elk model definieert een differentiaalvergelijking

$$
\frac{dV}{dt} = f(V, t; \theta)
$$


Omdat deze vergelijking niet analytisch oplosbaar is voor de meeste modellen, wordt de oplossing benaderd met numerieke methoden:

- Euler

- Heun (verbeterde Euler)

- Runge–Kutta 4 (standaardmethode)

Deze solver berekent een vloeiende tijdreeks 
Vmodel(t)
## Vergelijken met data (log-MSE)

Omdat tumorvolumes vaak sterk variëren (soms over meerdere ordes van grootte), gebruiken we logaritmische Mean Squared Error (log-MSE):
$$
\text{MSE}_{\log}
= \frac{1}{N}
\sum_{i=1}^{N}
\left[
\ln\big(V_{\text{model}}(t_i)\big)
-
\ln\big(V_{\text{data}}(t_i)\big)
\right]^2
$$

Voordelen:

- voorkomt dat grote waarden te veel gewicht krijgen

- voorkomt negatieve volumes (log clamp)

- beter geschikt voor biologisch exponentiële processen

Om dit te berekenen wordt de modeloplossing eerst geïnterpoleerd naar de exacte tijdstippen van de datapunten.
## Parameteroptimalisatie met hill-climbing random search

Omdat de foutfunctie vaak niet glad is en meerdere minima kan bevatten, gebruiken we een robuust maar eenvoudig optimalisatie-algoritme:

Hill-climbing random search

- begin met een willekeurig parameterpunt binnen grenzen

- voeg kleine willekeurige Gauss-stappen toe aan alle parameters

- accepteer een nieuwe parameter alleen als de fout kleiner wordt

- herhaal dit tot er geen verbetering meer optreedt

Dit algoritme werkt goed voor:

- niet-lineaire groeimodellen

- ruis in data
- modellen zonder analytische oplossing
## Vergelijken van modellen met informatiecriteria

Voor modelselectie worden AIC, BIC en AICc gebruikt:

- AIC: straft extra parameters licht

- BIC: zwaardere straf → kiest simpelere modellen

- AICc: verbeterde AIC voor kleine datasets (zoals bij tumorgroeidata)

Formules:
$$
\text{AIC}
=
n \cdot \ln(\text{MSE})
+
2k
$$

$$
\text{BIC}
=
n \cdot \ln(\text{MSE})
+
k \cdot \ln(n)
$$

$$
\text{AICc}
=
\text{AIC}
+
\frac{2k(k+1)}{\,n - k - 1\,}
$$

waarbij
- n = aantal datapunten

- k = aantal parameters van het model

Het model met de laagste waarde wordt beschouwd als het best passende.
## Samenvatting van het proces

1. ODE-model kiezen

2. Modelparameters optimaliseren via random search

3. Numerieke simulatie uitvoeren

4. MSE berekenen

5. AIC/BIC vergelijken

6. Beste model selecteren en visualiseren