"""
ADMET Compass - Interactive ADMET Prediction Dashboard
Streamlit Application

Author: Muriel (Biolevate)
"""

import streamlit as st
import pandas as pd
from rdkit import Chem
from utils.admet_predictions import (
    calculate_admet_profile, 
    interpret_logp, 
    interpret_logs, 
    interpret_herg, 
    interpret_qed
)
from utils.visualizations import (
    create_radar_chart,
    molecule_to_image,
    create_comparison_chart,
    create_lipinski_violations_chart
)

# Page configuration
st.set_page_config(
    page_title="ADMET Compass",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        padding: 1rem 0;
    }
    .subheader {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        padding-bottom: 2rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .warning-box {
        background-color: #fff3cd;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #ffc107;
    }
    .success-box {
        background-color: #d4edda;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #28a745;
    }
    </style>
    """, unsafe_allow_html=True)

# Load demo data
@st.cache_data
def load_demo_data():
    """Load demo molecules dataset"""
    df = pd.read_csv('data/demo_molecules.csv')
    return df

# Initialize session state
if 'selected_molecules' not in st.session_state:
    st.session_state.selected_molecules = []

def main():
    # Header
    st.markdown('<div class="main-header">🧬 ADMET Compass</div>', unsafe_allow_html=True)
    st.markdown('<div class="subheader">Predict ADMET properties in 30 seconds | Pharma R&D Decision Support Tool</div>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("📥 Input Configuration")
        
        input_mode = st.radio(
            "Select Input Mode:",
            ["Demo Gallery", "Custom SMILES"],
            help="Choose demo molecules or enter your own SMILES"
        )
        
        st.divider()
        
        if input_mode == "Demo Gallery":
            demo_df = load_demo_data()
            
            selected_molecule = st.selectbox(
                "Select Molecule:",
                demo_df['name'].tolist(),
                help="Choose from validated pharmaceutical blockbusters"
            )
            
            molecule_row = demo_df[demo_df['name'] == selected_molecule].iloc[0]
            smiles = molecule_row['smiles']
            molecule_name = molecule_row['name']
            category = molecule_row['category']
            sales = molecule_row['peak_sales_usd']
            description = molecule_row['description']
            
            # Display molecule info
            st.info(f"**Category:** {category}\n\n**Peak Sales:** {sales}\n\n{description}")
            
        else:
            molecule_name = st.text_input(
                "Molecule Name (optional):",
                placeholder="e.g., Compound A"
            )
            
            smiles = st.text_area(
                "Enter SMILES:",
                placeholder="CC(=O)Oc1ccccc1C(=O)O",
                help="Enter the SMILES notation of your molecule"
            )
        
        st.divider()
        
        # About section
        with st.expander("ℹ️ About ADMET Compass"):
            st.markdown("""
            **ADMET Compass** is a rapid screening tool for drug discovery professionals.
            
            **Key Features:**
            - 4 critical ADMET properties
            - Lipinski Rule of Five compliance
            - Interactive visualizations
            - Benchmarking against blockbusters
            
            **Developed by:** Muriel  
            **Context:** Biolevate R&D Portfolio  
            **Use Case:** Early-stage hit-to-lead optimization
            """)
    
    # Main content
    if smiles:
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("❌ Invalid SMILES notation. Please check your input.")
            return
        
        # Calculate ADMET profile
        with st.spinner("Calculating ADMET properties..."):
            admet_profile = calculate_admet_profile(smiles)
        
        if admet_profile is None:
            st.error("❌ Error calculating ADMET properties.")
            return
        
        # Layout: 2 columns
        col1, col2 = st.columns([1, 1.5])
        
        with col1:
            st.subheader("🔬 Molecular Structure")
            
            # Display 2D structure
            img = molecule_to_image(smiles, size=(350, 350))
            if img:
                st.image(img, use_container_width=True)
            
            # Display SMILES
            st.code(smiles, language=None)
            
            # Lipinski violations
            st.subheader("📋 Lipinski Rule of Five")
            violations = admet_profile['lipinski_violations']
            
            if violations == 0:
                st.success("✅ No violations - Excellent oral bioavailability potential")
            elif violations == 1:
                st.warning(f"⚠️ {violations} violation - Acceptable for some drug classes")
            else:
                st.error(f"❌ {violations} violations - Poor drug-like properties")
        
        with col2:
            st.subheader("📊 ADMET Profile Radar")
            
            # Radar chart
            radar_fig = create_radar_chart(admet_profile, molecule_name if molecule_name else "Molecule")
            st.plotly_chart(radar_fig, use_container_width=True)
        
        # Properties section
        st.divider()
        st.subheader("📈 Detailed ADMET Properties")
        
        # 4 columns for metrics
        metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)
        
        with metric_col1:
            st.metric(
                label="Lipophilicity (LogP)",
                value=f"{admet_profile['logp']}",
                help="Optimal: 0-3 (Lipinski)"
            )
            interpretation = interpret_logp(admet_profile['logp'])
            if "✅" in interpretation:
                st.success(interpretation)
            elif "⚠️" in interpretation:
                st.warning(interpretation)
            else:
                st.error(interpretation)
        
        with metric_col2:
            st.metric(
                label="Solubility (LogS)",
                value=f"{admet_profile['logs']}",
                help="Optimal: > -4"
            )
            interpretation = interpret_logs(admet_profile['logs'])
            if "✅" in interpretation:
                st.success(interpretation)
            elif "⚠️" in interpretation:
                st.warning(interpretation)
            else:
                st.error(interpretation)
        
        with metric_col3:
            st.metric(
                label="hERG Risk",
                value=f"{admet_profile['herg_risk']}",
                help="Cardiotoxicity probability (0-1, lower is better)"
            )
            interpretation = interpret_herg(admet_profile['herg_risk'])
            if "✅" in interpretation:
                st.success(interpretation)
            elif "⚠️" in interpretation:
                st.warning(interpretation)
            else:
                st.error(interpretation)
        
        with metric_col4:
            st.metric(
                label="Drug-likeness (QED)",
                value=f"{admet_profile['qed']}",
                help="Quantitative Estimate (0-1, higher is better)"
            )
            interpretation = interpret_qed(admet_profile['qed'])
            if "✅" in interpretation:
                st.success(interpretation)
            elif "⚠️" in interpretation:
                st.warning(interpretation)
            else:
                st.error(interpretation)
        
        # Comparison section
        st.divider()
        st.subheader("🔄 Benchmark Comparison")
        
        comparison_molecule = st.selectbox(
            "Compare with:",
            ["None"] + load_demo_data()['name'].tolist(),
            help="Select a reference molecule for comparison"
        )
        
        if comparison_molecule != "None":
            ref_smiles = load_demo_data()[load_demo_data()['name'] == comparison_molecule].iloc[0]['smiles']
            ref_profile = calculate_admet_profile(ref_smiles)
            
            comp_col1, comp_col2 = st.columns(2)
            
            with comp_col1:
                st.markdown(f"**{molecule_name if molecule_name else 'Your Molecule'}**")
                st.plotly_chart(create_radar_chart(admet_profile, "Your Molecule"), use_container_width=True)
            
            with comp_col2:
                st.markdown(f"**{comparison_molecule} (Reference)**")
                st.plotly_chart(create_radar_chart(ref_profile, comparison_molecule), use_container_width=True)
        
        # Export section
        st.divider()
        
        export_col1, export_col2 = st.columns([3, 1])
        
        with export_col1:
            st.info("💡 **Next Steps:** Export this report for project documentation or integrate into your R&D pipeline.")
        
        with export_col2:
            if st.button("📄 Generate PDF Report", type="primary"):
                st.success("✅ PDF generation coming soon!")
                st.balloons()
    
    else:
        # Landing page when no molecule selected
        st.info("👈 Select a molecule from the sidebar to begin ADMET analysis")
        
        # Display demo gallery preview
        st.subheader("🏆 Demo Gallery - Pharmaceutical Blockbusters")
        
        demo_df = load_demo_data()
        
        # Display table
        display_df = demo_df[['name', 'category', 'peak_sales_usd', 'description']].copy()
        display_df.columns = ['Molecule', 'Category', 'Peak Sales', 'Description']
        
        st.dataframe(display_df, use_container_width=True, hide_index=True)

# Footer
st.divider()
st.markdown("""
    <div style='text-align: center; color: #666; padding: 1rem;'>
        <p><strong>ADMET Compass</strong> | Developed by Muriel | Biolevate R&D Innovation</p>
        <p style='font-size: 0.8rem;'>Powered by RDKit • Streamlit • Python</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
