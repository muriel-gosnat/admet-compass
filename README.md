# 🧬 ADMET Compass

**Interactive ADMET Prediction Dashboard for Drug Discovery**

Rapid screening tool for pharmaceutical R&D teams to evaluate drug-likeness and ADMET properties of small molecules in early-stage discovery.

---

## 🎯 Overview

ADMET Compass enables medicinal chemists and project managers to:
- **Predict** 4 critical ADMET properties in seconds
- **Visualize** molecular profiles with interactive radar charts
- **Benchmark** compounds against validated pharmaceutical blockbusters
- **Document** results for project reviews and regulatory submissions

**Developed by:** Muriel (PhD Chemistry, MBA AI & Data Innovation)  
**Context:** Biolevate - AI Product Specialist in Pharma R&D  
**Use Case:** Hit-to-lead optimization decision support

---

## 🚀 Quick Start

### Installation

1. **Clone or download this repository**
```bash
cd admet-compass
```

2. **Install dependencies**
```bash
pip install -r requirements.txt
```

3. **Run the application**
```bash
streamlit run app.py
```

4. **Open your browser**
```
Local URL: http://localhost:8501
```

---

## 📊 Key Features

### 4 Critical ADMET Properties

1. **Lipophilicity (LogP)**
   - Optimal range: 0-3 (Lipinski Rule)
   - Impact: Membrane permeability, bioavailability

2. **Aqueous Solubility (LogS)**
   - Optimal: > -4
   - Impact: Formulation feasibility, absorption

3. **hERG Cardiotoxicity Risk**
   - Target: < 30% probability
   - Impact: Safety, regulatory approval

4. **Drug-likeness (QED)**
   - Scale: 0-1 (higher is better)
   - Composite score: MW, LogP, HBD, HBA, PSA, rotatable bonds

### Lipinski Rule of Five
- Automatic compliance check
- Visual alerts for violations
- Guidance for optimization

### Interactive Visualizations
- **Radar charts** for ADMET profiles
- **Side-by-side comparison** with reference molecules
- **Real-time updates** on input changes

### Demo Gallery
- 15 pharmaceutical blockbusters (Aspirin, Lipitor, Gleevec, etc.)
- Total market value: >$50B peak sales
- Validated SMILES from ChEMBL/PubChem

---

## 🛠️ Technical Stack

| Component | Technology | Purpose |
|-----------|-----------|---------|
| Frontend | Streamlit | Interactive web interface |
| Chemistry | RDKit | Molecular manipulation, descriptors |
| ML/Predictions | ESOL, QED | ADMET property estimation |
| Visualization | Plotly | Interactive charts |
| Data | Pandas | Dataset management |

---

## 📁 Project Structure

```
admet-compass/
├── app.py                      # Main Streamlit application
├── requirements.txt            # Python dependencies
├── README.md                   # This file
├── data/
│   └── demo_molecules.csv     # 15 validated blockbusters
└── utils/
    ├── admet_predictions.py   # ADMET calculation engine
    └── visualizations.py      # Plotly charts & rendering
```

---

## 📈 Use Cases

### 1. Hit-to-Lead Optimization
Screen compound libraries for drug-likeness before synthesis investment.

### 2. Project Portfolio Decisions
Rank and prioritize drug candidates based on ADMET profiles.

### 3. Regulatory Preparation
Generate preliminary ADMET reports for IND documentation.

### 4. Team Communication
Share visual ADMET comparisons with cross-functional stakeholders (chemistry, DMPK, toxicology).

---

## 🔬 Scientific Background

### ADMET Prediction Methods

**LogP** - Wildman-Crippen method (RDKit native)
- Empirical model from 2,000+ compounds
- Accounts for atom types and connectivity

**LogS** - ESOL (Delaney) formula
- Linear regression on descriptors (LogP, MW, RB, AP)
- R² = 0.87 on training set (>2,800 molecules)

**hERG Risk** - Heuristic scoring
- Based on: LogP, MW, aromatic rings, basicity
- Conservative threshold for cardiotoxicity

**QED** - Bickerton et al. (2012)
- Weighted geometric mean of 8 properties
- Optimized to distinguish drugs from non-drugs

---

## 🎓 Learning Resources

### RDKit Documentation
https://www.rdkit.org/docs/

### Lipinski's Rule of Five
Lipinski CA et al. (2001) Adv Drug Deliv Rev 46:3-26

### QED Drug-likeness
Bickerton GR et al. (2012) Nat Chem 4:90-98

### ESOL Solubility
Delaney JS (2004) J Chem Inf Comput Sci 44:1000-1005

---

## 🚧 Roadmap (Future Enhancements)

### Phase 2
- [ ] Integration with ChEMBL API for live data
- [ ] PDF report generation with branding
- [ ] Batch processing (upload CSV of SMILES)
- [ ] hERG prediction with ML model (DeepChem)

### Phase 3
- [ ] Additional ADMET properties (CYP450, BBB permeability)
- [ ] Structure-Activity Relationship (SAR) analysis
- [ ] Export to ELN systems (e.g., Benchling)
- [ ] Multi-language support (French/English)

---

## 💼 Professional Context

This tool demonstrates:
- ✅ **Domain expertise** in pharma R&D workflows
- ✅ **Technical skills** in Python, ML, data visualization
- ✅ **Product thinking** for user-centric design
- ✅ **Regulatory awareness** (EU AI Act, GxP compliance)

**Target audience for this portfolio piece:**
- Pharma/biotech hiring managers (Product Owner, Program Manager roles)
- R&D Directors evaluating AI governance candidates
- Recruiters screening for technical + scientific background

---

## 📞 Contact

**Muriel**  
Product Specialist IA | Biolevate  
PhD Chemistry | MBA AI & Data Innovation  

LinkedIn: [Your LinkedIn URL]  
GitHub: [Your GitHub URL]

---

## 📄 License

This project is a professional portfolio demonstration.  
Educational and non-commercial use permitted with attribution.

---

## 🙏 Acknowledgments

- **RDKit Community** - Open-source cheminformatics toolkit
- **Streamlit Team** - Rapid prototyping framework
- **Biolevate R&D** - Real-world pharma context and use cases

---

**Last Updated:** February 2026  
**Version:** 1.0.0 (MVP)
