# BLAST & Docking — Cyanopeptides × Protéines membranaires

Pipeline bioinformatique pour
(i) l’identification de protéines membranaires homologues par BLAST
(ii) l’étude de leurs interactions avec des cyanopeptides par docking moléculaire (AutoDock Vina).

## BLAST — Recherche de protéines membranaires homologues
---

### Étape 1 — Décompresser les protéomes

Dézipper les fichiers sont compressés (.zip) dans le dossier `blast/proteomes`.

---

### Étape 2 — Installer l’environnement Conda

Dans le terminal, placez-vous à la racine du projet (là où se trouve `environment.yml`) puis :

```bash
conda env create -f environment.yml
conda activate blast_docking
```
Vérifier que BLAST est bien installé : 

```bash
makeblastdb --version
blastp -version
```
---

### Étape 3 — Lancer le BLAST local

```bash
python blast_local/scripts/blast_final.py
```
Les résultats sont générés automatiquement dans le dossier `blast/resultats/`

Un fichier Excel récapitulatif est ensuite créé, regroupant
17 protéines sélectionnées communes aux 5 organismes étudiés :`./docking/proteins/new/1_proteine_test.xlsx`.

## DOCKING moléculaire

### 1) Génération des structures PDB (AlphaFold)

Téléchargement automatique des structures AlphaFold des protéines sélectionnées :

```bash
python docking/scripts/0_import_pdb.py
```
Les fichiers PDB sont générés dans : `./docking/proteins/new/`

### 2) Docking protéines–ligands connus (jeu de référence)

Docking de complexes protéines–cyanopeptides décrits dans la littérature, utilisés comme référence.

Le script :
	•	nettoie les protéines (suppression eau, ions, ligands, glycans),
	•	convertit protéines et ligands en PDBQT,
	•	génère des structures 3D des ligands,
	•	lance AutoDock Vina.

Paramètres de docking :
	•	Blind docking centré sur la protéine
	•	Taille de boîte :
	•	24 × 24 × 24 Å (enzymes)
	•	30 × 30 × 30 Å (transporteurs)
	•	exhaustiveness = 16
	•	num_modes = 50

```bash
python docking/scripts/run_known.py
```
Les resultats sont dans `docking/results/`


### 3) Docking rapide

Docking rapide des protéines issues du BLAST avec les cyanopeptides étudiés.

Objectif :
	•	criblage large,
	•	identification des interactions potentielles,
	•	sélection des complexes à affiner.

```bash
python docking/scripts/2_run_fast.py
```

### 4) Docking affiné

Les meilleurs complexes issus du docking rapide sont redockés avec des critères plus stricts
(affinité, RMSD, type de protéine).

```bash
python docking/scripts/3_run_refine.py
```
Les résultats finaux sont synthétisés dans des fichiers CSV : 
`results/summary_fast.csv`
`results/summary_refine_all.csv`
`results/summary_refine_filtered.csv`

### 5) Analyse des interactions et des poches

Calcul :
	•	propriétés des ligands (RDKit),
	•	composition des poches de liaison,
	•	liaisons hydrogène (PLIP).

```bash
python docking/scripts/4_ajout_colonne.py
```

### 6) Analyses statistiques et visualisation
Analyses exploratoires, statistiques et PCA :
```bash
docking/scripts/5_analyses.ipynb
```






