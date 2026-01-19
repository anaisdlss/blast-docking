# Blast

---

## Étape 1 — Dézipper les protéomes

Dézipper les fichiers sont compressés (.zip) dans le dossier `blast/proteomes`.

---

## Étape 2 — Installer l’environnement Conda

### 1) Créer et activer l’environnement
Dans le terminal, placez-vous à la racine du projet (là où se trouve `environment.yml`) puis :

```bash
conda env create -f environment.yml
conda activate blast_docking
```
### 2) Vérifier que BLAST est bien installé

```bash
makeblastdb --version
blastp -version
```
---

## Etape 3 - Lancer le blast

```bash
python blast_local/scripts/blast_final.py
```
Les résultats sont générés automatiquement dans le dossier `resultats/`

