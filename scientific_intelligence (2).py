"""
Scientific Intelligence Module - ADMET Compass
Integrates PubMed, ClinicalTrials.gov, EPO Espacenet and Claude API
to generate a structured scientific digest for a given molecule.
"""

import requests
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import anthropic
import streamlit as st
import json
import time


# ─────────────────────────────────────────────
# 1. PUBCHEM : SMILES → Nom de la molécule
# ─────────────────────────────────────────────

def smiles_to_name(smiles: str) -> str | None:
    """Convertit un SMILES en nom courant via PubChem."""
    try:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON"
        params = {"smiles": smiles}
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()
        props = data["PC_Compounds"][0].get("props", [])
        for prop in props:
            label = prop.get("urn", {}).get("label", "")
            if label == "IUPAC Name":
                return prop["value"]["sval"]
        # Fallback: preferred name
        cid = data["PC_Compounds"][0]["id"]["id"]["cid"]
        r2 = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON",
            timeout=10
        )
        props2 = r2.json()["PropertyTable"]["Properties"][0]
        return props2.get("Title") or props2.get("IUPACName")
    except Exception:
        return None


def resolve_molecule_name(user_input: str) -> str:
    """
    Retourne le nom exploitable pour les requêtes API.
    Si l'input ressemble à un SMILES, tente la conversion via PubChem.
    """
    # Heuristique simple : un SMILES contient =, (, ), # ou des chiffres intercalés
    looks_like_smiles = any(c in user_input for c in ["=", "(", ")", "#", "@"]) or (
        sum(c.isdigit() for c in user_input) > 2
    )
    if looks_like_smiles:
        name = smiles_to_name(user_input)
        return name if name else user_input
    return user_input


# ─────────────────────────────────────────────
# 2. PUBMED : Publications récentes
# ─────────────────────────────────────────────

def fetch_pubmed(molecule_name: str, max_results: int = 20) -> list[dict]:
    """
    Interroge PubMed via E-utilities et retourne une liste de publications.
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    query = f"{molecule_name}[Title/Abstract] AND (ADMET OR pharmacokinetics OR toxicity OR bioavailability)"

    # Recherche des IDs
    search_url = f"{base}esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "sort": "date",
        "retmode": "json",
        "datetype": "pdat",
        "reldate": 1825,  # 5 dernières années
    }
    try:
        r = requests.get(search_url, params=params, timeout=10)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        return []

    if not ids:
        return []

    # Récupération des résumés
    fetch_url = f"{base}efetch.fcgi"
    fetch_params = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "xml",
        "rettype": "abstract",
    }
    try:
        r2 = requests.get(fetch_url, params=fetch_params, timeout=15)
        r2.raise_for_status()
        root = ET.fromstring(r2.content)
    except Exception:
        return []

    articles = []
    for article in root.findall(".//PubmedArticle"):
        try:
            title = article.findtext(".//ArticleTitle") or "N/A"
            abstract = article.findtext(".//AbstractText") or ""
            year = article.findtext(".//PubDate/Year") or "N/A"
            journal = article.findtext(".//Journal/Title") or "N/A"
            pmid = article.findtext(".//PMID") or ""
            articles.append({
                "title": title,
                "abstract": abstract[:500],  # Tronqué pour le LLM
                "year": year,
                "journal": journal,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })
        except Exception:
            continue

    return articles


# ─────────────────────────────────────────────
# 3. CLINICALTRIALS.GOV : Essais cliniques
# ─────────────────────────────────────────────

def fetch_clinical_trials(molecule_name: str, max_results: int = 20) -> list[dict]:
    """
    Interroge ClinicalTrials.gov API v2 et retourne les essais pertinents.
    """
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": molecule_name,
        "pageSize": max_results,
        "fields": "NCTId,BriefTitle,Phase,OverallStatus,Condition,StartDate,LeadSponsorName",
        "sort": "StartDate:desc",
    }
    try:
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return []

    trials = []
    for study in data.get("studies", []):
        proto = study.get("protocolSection", {})
        id_mod = proto.get("identificationModule", {})
        status_mod = proto.get("statusModule", {})
        design_mod = proto.get("designModule", {})
        cond_mod = proto.get("conditionsModule", {})
        sponsor_mod = proto.get("sponsorCollaboratorsModule", {})

        trials.append({
            "nct_id": id_mod.get("nctId", "N/A"),
            "title": id_mod.get("briefTitle", "N/A"),
            "phase": ", ".join(design_mod.get("phases", [])) or "N/A",
            "status": status_mod.get("overallStatus", "N/A"),
            "conditions": ", ".join(cond_mod.get("conditions", [])),
            "start_date": status_mod.get("startDateStruct", {}).get("date", "N/A"),
            "sponsor": sponsor_mod.get("leadSponsor", {}).get("name", "N/A"),
            "url": f"https://clinicaltrials.gov/study/{id_mod.get('nctId', '')}",
        })

    return trials


# ─────────────────────────────────────────────
# 4. EPO ESPACENET : Brevets
# ─────────────────────────────────────────────

def fetch_patents(molecule_name: str, max_results: int = 10) -> list[dict]:
    """
    Interroge l'API Open Patent Services (OPS) de l'EPO.
    Nécessite un compte EPO gratuit pour les clés API.
    Fallback sur une requête publique limitée si pas de clé.
    """
    # API OPS EPO - endpoint public (limité à 25 req/30s)
    url = "https://ops.epo.org/3.2/rest-services/published-data/search"
    headers = {"Accept": "application/json"}
    query = f'ti="{molecule_name}" OR ab="{molecule_name}"'
    params = {
        "q": query,
        "Range": f"1-{min(max_results, 25)}",
    }

    try:
        r = requests.get(url, headers=headers, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()
        results = data.get("ops:world-patent-data", {}) \
                      .get("ops:biblio-search", {}) \
                      .get("ops:search-result", {}) \
                      .get("ops:publication-reference", [])

        # Normalisation en liste si un seul résultat
        if isinstance(results, dict):
            results = [results]

        patents = []
        for p in results:
            doc_id = p.get("@document-id-type", "")
            country = p.get("country", {}).get("$", "")
            doc_number = p.get("doc-number", {}).get("$", "")
            kind = p.get("kind", {}).get("$", "")
            patents.append({
                "id": f"{country}{doc_number}{kind}",
                "country": country,
                "number": doc_number,
                "url": f"https://worldwide.espacenet.com/patent/search?q={country}{doc_number}",
            })
        return patents

    except Exception:
        # Fallback : recherche Google Patents via SerpAPI ou message d'info
        return []


def fetch_patents_count_google(molecule_name: str) -> int:
    """
    Estimation du nombre de brevets via une requête simple sur Google Patents.
    Alternative légère si l'API EPO pose problème.
    """
    # Note : nécessite SerpAPI ou une clé API Google Custom Search
    # Retourne 0 si non configuré, géré dans l'affichage
    return 0


# ─────────────────────────────────────────────
# 5. SYNTHÈSE LLM avec Claude
# ─────────────────────────────────────────────

def generate_digest_with_claude(
    molecule_name: str,
    articles: list[dict],
    trials: list[dict],
    patents: list[dict],
    api_key: str,
) -> dict:
    """
    Envoie les données brutes à Claude et retourne un digest structuré en JSON.
    """
    client = anthropic.Anthropic(api_key=api_key)

    # Construction du contexte
    articles_text = "\n".join([
        f"- [{a['year']}] {a['title']} ({a['journal']}): {a['abstract'][:200]}..."
        for a in articles[:15]
    ])
    trials_text = "\n".join([
        f"- {t['phase']} | {t['status']} | {t['conditions']} | Sponsor: {t['sponsor']}"
        for t in trials[:15]
    ])
    patents_text = f"{len(patents)} patents found" if patents else "Patent data unavailable."

    prompt = f"""You are an expert in Drug Discovery and Pharmaceutical Strategic Intelligence.
Analyze the following data on the molecule **{molecule_name}** and generate a structured digest.

## Recent PubMed Publications ({len(articles)} found)
{articles_text if articles_text else "No publications found."}

## ClinicalTrials.gov Clinical Trials ({len(trials)} found)
{trials_text if trials_text else "No clinical trials found."}

## Patents (EPO Espacenet)
{patents_text}

Return ONLY a valid JSON object with this exact structure:
{{
  "recent_trends": "4-5 sentence synthesis of dominant research directions, new indications, emerging toxicity signals.",
  "clinical_stage": {{
    "summary": "2-3 sentence narrative of the clinical landscape.",
    "phases": {{
      "Phase I": 0,
      "Phase II": 0,
      "Phase III": 0,
      "Phase IV": 0,
      "Not specified": 0
    }},
    "main_conditions": ["condition1", "condition2", "condition3"],
    "main_sponsors": ["sponsor1", "sponsor2"]
  }},
  "patent_landscape": {{
    "summary": "2-3 sentence narrative on IP status.",
    "recent_activity": "Low / Moderate / High",
    "ip_status": "Active patent / Expired patent / Mixed / Unknown"
  }},
  "gaps_opportunities": "3-4 sentences identifying underexplored angles, unmet needs, and strategic opportunities.",
  "maturity_score": {{
    "score": 3,
    "justification": "1-2 sentence explanation of the score out of 5."
  }}
}}

The maturity_score must be an integer between 1 and 5:
1 = Very preliminary, 2 = Exploratory, 3 = Active development, 4 = Mature, 5 = Very mature/commercialized.
"""

    try:
        message = client.messages.create(
            model="claude-opus-4-5",
            max_tokens=1500,
            messages=[{"role": "user", "content": prompt}]
        )
        raw = message.content[0].text.strip()
        # Nettoyage éventuel des backticks markdown
        if raw.startswith("```"):
            raw = raw.split("```")[1]
            if raw.startswith("json"):
                raw = raw[4:]
        return json.loads(raw)
    except Exception as e:
        return {"error": str(e)}


# ─────────────────────────────────────────────
# 6. AFFICHAGE STREAMLIT
# ─────────────────────────────────────────────

def display_scientific_intelligence(molecule_input: str, api_key: str):
    """
    Fonction principale à appeler depuis app.py dans l'onglet Scientific Intelligence.
    
    Usage dans app.py :
        from scientific_intelligence import display_scientific_intelligence
        
        with tab_sci_intel:
            display_scientific_intelligence(molecule_input, ANTHROPIC_API_KEY)
    """

    st.markdown("## 🔬 Scientific Intelligence")
    st.caption(f"Real-time analysis for: **{molecule_input}**")

    if st.button("🚀 Generate Scientific Digest", type="primary"):

        # Resolve molecule name
        with st.spinner("Identifying molecule..."):
            molecule_name = resolve_molecule_name(molecule_input)
            st.info(f"Molecule analyzed: **{molecule_name}**")

        # Data collection
        col1, col2, col3 = st.columns(3)

        with col1:
            with st.spinner("Searching PubMed..."):
                articles = fetch_pubmed(molecule_name)
            st.metric("Publications found", len(articles))

        with col2:
            with st.spinner("Searching ClinicalTrials.gov..."):
                trials = fetch_clinical_trials(molecule_name)
            st.metric("Clinical trials", len(trials))

        with col3:
            with st.spinner("Searching EPO patents..."):
                patents = fetch_patents(molecule_name)
            st.metric("Patents (EPO)", len(patents) if patents else "N/A")

        # Claude synthesis
        with st.spinner("Generating digest with Claude..."):
            digest = generate_digest_with_claude(molecule_name, articles, trials, patents, api_key)

        if "error" in digest:
            st.error(f"Synthesis error: {digest['error']}")
            return

        st.divider()

        # ── Maturity score ──
        score = digest.get("maturity_score", {})
        score_val = score.get("score", 0)
        score_labels = {1: "Very preliminary", 2: "Exploratory", 3: "Active development", 4: "Mature", 5: "Very mature"}
        score_colors = {1: "🔴", 2: "🟠", 3: "🟡", 4: "🟢", 5: "🟢"}

        st.markdown(f"### {score_colors.get(score_val, '⚪')} Scientific Maturity Score: {score_val}/5 — *{score_labels.get(score_val, '')}*")
        st.caption(score.get("justification", ""))

        st.divider()

        # ── Recent trends ──
        st.markdown("### 📈 Recent Trends")
        st.write(digest.get("recent_trends", "N/A"))

        # ── Clinical stage ──
        st.markdown("### 🏥 Clinical Stage")
        clinical = digest.get("clinical_stage", {})
        st.write(clinical.get("summary", ""))

        phases = clinical.get("phases", {})
        if any(v > 0 for v in phases.values()):
            cols = st.columns(len(phases))
            for col, (phase, count) in zip(cols, phases.items()):
                col.metric(phase, count)

        if clinical.get("main_conditions"):
            st.markdown("**Main conditions:** " + " | ".join(clinical["main_conditions"]))
        if clinical.get("main_sponsors"):
            st.markdown("**Main sponsors:** " + " | ".join(clinical["main_sponsors"]))

        # ── Patent landscape ──
        st.markdown("### 🔒 Patent Landscape")
        patents_digest = digest.get("patent_landscape", {})
        st.write(patents_digest.get("summary", ""))
        col_a, col_b = st.columns(2)
        col_a.metric("Filing activity", patents_digest.get("recent_activity", "N/A"))
        col_b.metric("IP status", patents_digest.get("ip_status", "N/A"))

        # ── Gaps & opportunities ──
        st.markdown("### 💡 Gaps & Opportunities")
        st.info(digest.get("gaps_opportunities", "N/A"))

        # ── Sources ──
        with st.expander("📚 Sources — PubMed Publications"):
            for a in articles[:10]:
                st.markdown(f"- [{a['title']} ({a['year']})]({a['url']})")

        with st.expander("🏥 Sources — Clinical Trials"):
            for t in trials[:10]:
                st.markdown(f"- [{t['title']}]({t['url']}) | {t['phase']} | {t['status']}")

        if patents:
            with st.expander("🔒 Sources — EPO Patents"):
                for p in patents[:10]:
                    st.markdown(f"- [{p['id']}]({p['url']})")

        # ── Export ──
        st.divider()
        export_data = {
            "molecule": molecule_name,
            "date": datetime.now().isoformat(),
            "digest": digest,
            "sources": {
                "pubmed_count": len(articles),
                "trials_count": len(trials),
                "patents_count": len(patents),
            }
        }
        st.download_button(
            label="⬇️ Export digest (JSON)",
            data=json.dumps(export_data, ensure_ascii=False, indent=2),
            file_name=f"digest_{molecule_name.replace(' ', '_')}_{datetime.now().strftime('%Y%m%d')}.json",
            mime="application/json",
        )
