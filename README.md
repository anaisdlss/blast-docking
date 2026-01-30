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





