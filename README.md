## BLAST : Recherches des proteines membranaires homologues

---

### Étape 1 — Dézipper les protéomes

Dézipper les fichiers sont compressés (.zip) dans le dossier `blast/proteomes`.

---

### Étape 2 — Installer l’environnement Conda

#### 1) Créer et activer l’environnement
Dans le terminal, placez-vous à la racine du projet (là où se trouve `environment.yml`) puis :

```bash
conda env create -f environment.yml
conda activate blast_docking
```
#### 2) Vérifier que BLAST est bien installé

```bash
makeblastdb --version
blastp -version
```
---

### Etape 3 - Lancer le blast

```bash
python blast_local/scripts/blast_final.py
```
Les résultats sont générés automatiquement dans le dossier `resultats/`

---

Nous pouvons ainsi créer un nouveau fichier excel qui regroupe 17 proteines selectionnées communes aux 5 organismes. Ce fichier `1_proteine_test.xlsx` se retrouve dans `./docking/proteins/new/`.

## DOCKING

### 1) Création des fichiers PDB

Nous allons dans un premier temps créer les fichiers PDB de nos proteines selectionnées dans le dossier `./docking/proteins/new`.

```bash
python docking/scripts/0_import_pdb.py
```

### 2) Docking proteines-ligands connu

Nous faisons ensuite un docking sur les complexes proteines-cyanopeptides connus de par la litterature, ceux utilisé pour faire un BLAST.

Ce script prépare les proteines en noettoyant les PDB (suppression de l'eau, ions, ligands, glycans...), les converti en PDBQT, puis prepare les ligands en generant des strcutures 3D, puis convertit en PDBQT. Le docking se fait avec AutoDock Vina avec les paramètres suivants : blind docking centré sur la proteine, taille de la boite 24x24x24 Å pour les enzymes, et 30x30x30 Å pour les transporteurs. Les paramètres sont : exhaustiveness = 16, num_modes = 50.
Les resultats sont dans `docking/results/`.

```bash
python docking/scripts/run_known.py
```

### 3) Docking rapide

Docking rapide des proteines selectionnées issus du blast avec les cyanopeptides partenaires.
Les paramètres sont de : 





