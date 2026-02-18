"""
Visualizations Module
Creates interactive Plotly charts for ADMET data
"""

import plotly.graph_objects as go
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image


def create_radar_chart(admet_profile, molecule_name=""):
    """
    Create radar chart for ADMET profile
    
    Args:
        admet_profile (dict): ADMET profile with normalized scores
        molecule_name (str): Name of the molecule
    
    Returns:
        plotly.graph_objects.Figure
    """
    if not admet_profile:
        return None
    
    normalized = admet_profile['normalized']
    
    categories = ['Lipophilicity<br>(LogP)', 
                  'Solubility<br>(LogS)', 
                  'hERG Safety', 
                  'Drug-likeness<br>(QED)']
    
    values = [normalized['logp'], 
              normalized['logs'], 
              normalized['herg_safety'], 
              normalized['qed']]
    
    # Add first value at the end to close the polygon
    values_closed = values + [values[0]]
    categories_closed = categories + [categories[0]]
    
    fig = go.Figure()
    
    # Add molecule trace
    fig.add_trace(go.Scatterpolar(
        r=values_closed,
        theta=categories_closed,
        fill='toself',
        name=molecule_name if molecule_name else 'Molecule',
        line=dict(color='rgb(99, 110, 250)', width=2),
        fillcolor='rgba(99, 110, 250, 0.3)'
    ))
    
    # Add optimal zone (reference)
    optimal_values = [0.8, 0.8, 0.8, 0.8, 0.8]
    fig.add_trace(go.Scatterpolar(
        r=optimal_values,
        theta=categories_closed,
        fill='toself',
        name='Optimal Zone',
        line=dict(color='rgb(34, 197, 94)', width=1, dash='dash'),
        fillcolor='rgba(34, 197, 94, 0.1)'
    ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1],
                tickfont=dict(size=10)
            ),
            angularaxis=dict(
                tickfont=dict(size=11)
            )
        ),
        showlegend=True,
        title={
            'text': f"ADMET Profile - {molecule_name}" if molecule_name else "ADMET Profile",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        height=450,
        margin=dict(l=80, r=80, t=80, b=40)
    )
    
    return fig


def create_comparison_chart(molecules_data):
    """
    Create comparison bar chart for multiple molecules
    
    Args:
        molecules_data (list): List of dicts with 'name' and 'admet_profile'
    
    Returns:
        plotly.graph_objects.Figure
    """
    names = [m['name'] for m in molecules_data]
    
    fig = go.Figure()
    
    # QED scores
    qed_scores = [m['admet_profile']['qed'] for m in molecules_data]
    
    fig.add_trace(go.Bar(
        x=names,
        y=qed_scores,
        name='Drug-likeness (QED)',
        marker_color='rgb(99, 110, 250)',
        text=[f"{score:.2f}" for score in qed_scores],
        textposition='outside'
    ))
    
    fig.update_layout(
        title="Drug-likeness Comparison (QED Score)",
        xaxis_title="Molecule",
        yaxis_title="QED Score",
        yaxis=dict(range=[0, 1]),
        height=400,
        showlegend=False
    )
    
    return fig


def create_property_gauge(value, title, min_val, max_val, optimal_range, unit=""):
    """
    Create gauge chart for single property
    
    Args:
        value (float): Property value
        title (str): Property name
        min_val (float): Minimum value for scale
        max_val (float): Maximum value for scale
        optimal_range (tuple): (min_optimal, max_optimal)
        unit (str): Unit of measurement
    
    Returns:
        plotly.graph_objects.Figure
    """
    # Determine color based on optimal range
    if optimal_range[0] <= value <= optimal_range[1]:
        color = "green"
    else:
        color = "orange" if abs(value - np.mean(optimal_range)) < (max_val - min_val) * 0.3 else "red"
    
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=value,
        domain={'x': [0, 1], 'y': [0, 1]},
        title={'text': f"{title} ({unit})" if unit else title},
        gauge={
            'axis': {'range': [min_val, max_val]},
            'bar': {'color': color},
            'steps': [
                {'range': [min_val, optimal_range[0]], 'color': "lightgray"},
                {'range': optimal_range, 'color': "lightgreen"},
                {'range': [optimal_range[1], max_val], 'color': "lightgray"}
            ],
            'threshold': {
                'line': {'color': "black", 'width': 2},
                'thickness': 0.75,
                'value': value
            }
        }
    ))
    
    fig.update_layout(height=250, margin=dict(l=20, r=20, t=50, b=20))
    
    return fig


def molecule_to_image(smiles, size=(300, 300)):
    """
    Convert SMILES to 2D structure image
    
    Args:
        smiles (str): SMILES string
        size (tuple): Image size (width, height)
    
    Returns:
        PIL.Image
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    img = Draw.MolToImage(mol, size=size)
    return img


def create_lipinski_violations_chart(violations_count):
    """
    Create visual for Lipinski violations
    
    Args:
        violations_count (int): Number of violations (0-4)
    
    Returns:
        plotly.graph_objects.Figure
    """
    categories = ['MW ≤ 500', 'LogP ≤ 5', 'HBD ≤ 5', 'HBA ≤ 10']
    
    # This is simplified - you'd need to check each rule individually
    # For now, we'll create a simple indicator
    
    if violations_count == 0:
        color = 'green'
        text = "✅ Compliant with Lipinski's Rule of Five"
    elif violations_count == 1:
        color = 'orange'
        text = f"⚠️ {violations_count} Lipinski violation"
    else:
        color = 'red'
        text = f"❌ {violations_count} Lipinski violations"
    
    fig = go.Figure(go.Indicator(
        mode="number+delta",
        value=violations_count,
        title={'text': "Lipinski Violations"},
        delta={'reference': 0},
        domain={'x': [0, 1], 'y': [0, 1]}
    ))
    
    fig.update_layout(
        height=150,
        annotations=[
            dict(
                text=text,
                x=0.5,
                y=0.3,
                xref='paper',
                yref='paper',
                showarrow=False,
                font=dict(size=12)
            )
        ]
    )
    
    return fig


import numpy as np  # Add this for the gauge function
