# 🚀 ADMET Compass - Guide de Démarrage Rapide

## ⚡ Démarrage en 5 minutes

### 1️⃣ Installation des dépendances

Ouvre un terminal dans le dossier `admet-compass` et exécute :

```bash
pip install -r requirements.txt
```

**Note:** Si tu n'as pas encore Python installé, télécharge Python 3.9+ depuis https://www.python.org/downloads/

### 2️⃣ Lancement de l'application

```bash
streamlit run app.py
```

L'application s'ouvrira automatiquement dans ton navigateur à l'adresse :
```
http://localhost:8501
```

### 3️⃣ Test rapide

1. Dans la sidebar, sélectionne **"Demo Gallery"**
2. Choisis **"Aspirin"** dans la liste
3. Observe le profil ADMET généré en temps réel !

---

## 🎯 Premiers pas

### Mode Demo Gallery (recommandé pour commencer)
- Sélectionne parmi 15 blockbusters pharmaceutiques validés
- Chaque molécule affiche son profil ADMET complet
- Compare avec d'autres molécules via le menu déroulant

### Mode Custom SMILES
- Entre ton propre code SMILES
- Exemple : `CC(=O)Oc1ccccc1C(=O)O` (Aspirin)
- L'analyse se met à jour automatiquement

---

## 📊 Comprendre les résultats

### Radar Chart ADMET
- **Zone verte** = Optimal Zone (propriétés idéales)
- **Ligne bleue** = Profil de ta molécule
- Plus la surface bleue est grande et proche du vert, meilleure est la molécule

### Métriques détaillées
Chaque propriété affiche :
- ✅ **Vert** = Dans la plage optimale
- ⚠️ **Orange** = Acceptable mais à surveiller
- ❌ **Rouge** = Hors des seuils recommandés

---

## 🐛 Dépannage

### Erreur "No module named rdkit"
```bash
pip install rdkit
```

Si ça ne fonctionne pas, utilise conda :
```bash
conda install -c conda-forge rdkit
```

### L'application ne se lance pas
Vérifie que tu es dans le bon dossier :
```bash
cd admet-compass
ls  # Tu devrais voir app.py, requirements.txt, etc.
```

### Port 8501 déjà utilisé
Streamlit utilise le port 8501 par défaut. Si occupé :
```bash
streamlit run app.py --server.port 8502
```

---

## 📦 Déploiement sur Streamlit Cloud (GRATUIT)

### Étapes pour rendre ton app publique

1. **Créer un compte GitHub** (si tu n'en as pas)
   - Va sur https://github.com
   - Crée un compte gratuit

2. **Upload ton projet sur GitHub**
   ```bash
   git init
   git add .
   git commit -m "ADMET Compass MVP"
   git remote add origin https://github.com/TON_USERNAME/admet-compass.git
   git push -u origin main
   ```

3. **Déployer sur Streamlit Cloud**
   - Va sur https://streamlit.io/cloud
   - Connecte ton compte GitHub
   - Clique sur "New app"
   - Sélectionne ton repo `admet-compass`
   - Clique sur "Deploy"

4. **Obtenir ton URL publique**
   - Tu obtiendras une URL du type : `https://admet-compass-XXXX.streamlit.app`
   - Partage cette URL sur ton CV et LinkedIn !

**Temps total : 10 minutes** ⏱️

---

## 💡 Utilisation en entretien

### Démo rapide (2-3 minutes)
1. "Je vais vous montrer un outil que j'ai développé pour l'optimisation ADMET"
2. Ouvre l'app → Sélectionne Atorvastatin (Lipitor, 13B$ de ventes)
3. Explique le radar chart : "Voici le profil ADMET d'un blockbuster validé"
4. Compare avec une autre molécule
5. "Cet outil permet d'accélérer les décisions en early discovery"

### Points clés à mentionner
- ✅ Développé en 2 jours avec Streamlit
- ✅ RDKit pour les calculs chimiques
- ✅ Dataset validé de 15 blockbusters (>$50B de ventes cumulées)
- ✅ Déployable en production sur cloud gratuit
- ✅ Extensible avec API ChEMBL, batch processing, ML models

---

## 🔄 Prochaines étapes

### Personnalisation rapide

**Change le logo/branding :**
Édite `app.py`, ligne 63 :
```python
st.markdown('<div class="main-header">🧬 TON_NOM - ADMET Compass</div>')
```

**Ajoute tes propres molécules :**
Édite `data/demo_molecules.csv`, ajoute des lignes :
```csv
Ma_Molecule,CC(C)CCO,Mon_Projet,0,Description de ma molécule
```

**Change les couleurs :**
Édite la section CSS dans `app.py` (lignes 27-52)

---

## 📞 Besoin d'aide ?

Si tu bloques sur un point technique :
1. Vérifie le README.md complet
2. Google l'erreur exacte (souvent une solution existe déjà)
3. Teste sur un autre navigateur (Chrome recommandé)

---

**Bon coding ! 🚀**

**Muriel - Biolevate R&D Innovation**
