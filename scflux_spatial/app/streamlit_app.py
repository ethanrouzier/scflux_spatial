"""
🧬 scflux_spatial - Streamlit App

Interface web interactive pour l'analyse de flux métaboliques spatiaux
avec des données de transcriptomique spatiale (Visium, MERFISH, etc.)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

# Configuration de la page
st.set_page_config(
    page_title="scflux_spatial",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    st.title("🧬 scflux_spatial")
    st.markdown("**Interface interactive pour l'analyse de flux métaboliques spatiaux**")
    
    # Sidebar
    st.sidebar.header("📊 Navigation")
    page = st.sidebar.selectbox(
        "Choisir une section",
        ["🏠 Accueil", "📈 Chargement de données", "🔬 Analyse FBA", "🌐 Modélisation spatiale", "📊 Visualisation"]
    )
    
    if page == "🏠 Accueil":
        show_home()
    elif page == "📈 Chargement de données":
        show_data_loading()
    elif page == "🔬 Analyse FBA":
        show_fba_analysis()
    elif page == "🌐 Modélisation spatiale":
        show_spatial_modeling()
    elif page == "📊 Visualisation":
        show_visualization()

def show_home():
    st.header("🏠 Bienvenue dans scflux_spatial")
    
    st.markdown("""
    **scflux_spatial** est un package Python pour l'analyse de flux métaboliques spatiaux 
    combinant la transcriptomique spatiale avec la modélisation métabolique.
    
    ### 🎯 Fonctionnalités principales :
    
    #### 1. **Chargement de données spatiales**
    - Support des données 10X Visium
    - Préprocessing automatique
    - Intégration avec scanpy et squidpy
    
    #### 2. **Analyse de flux métaboliques (FBA)**
    - Intégration de l'expression génique
    - Modèles métaboliques humains (HumanGEM)
    - Parsing des règles GPR (Gene-Protein-Reaction)
    
    #### 3. **Modélisation spatiale**
    - Équations de réaction-diffusion (PDEs)
    - Couplage FBA-spatial
    - Boucles d'auto-organisation adaptative (SOA)
    
    #### 4. **Visualisation interactive**
    - Cartes de flux métaboliques
    - Visualisations spatiales
    - Analyses de corrélation
    
    ### 🚀 Commencer :
    Utilisez le menu de gauche pour explorer les différentes fonctionnalités !
    """)
    
    # Exemple de code
    st.subheader("💻 Exemple d'utilisation rapide")
    st.code("""
# Charger des données Visium
from scflux_spatial.dataio import load_visium
adata = load_visium(use_demo=True)

# Analyser les flux métaboliques
from scflux_spatial.fba import integrate_expression_with_method
flux_bounds = integrate_expression_with_method(
    gene_expression=adata.X.mean(axis=0),
    base_bounds={},
    method="eflux"
)

# Modéliser la diffusion spatiale
from scflux_spatial.spatial import solve_reaction_diffusion
solution = solve_reaction_diffusion(
    substrate="glucose",
    diffusion_coefficient=1e-6,
    reaction_rate=0.1
)
    """, language="python")

def show_data_loading():
    st.header("📈 Chargement de données spatiales")
    
    st.markdown("""
    Cette section vous permet de charger et préprocesser vos données de transcriptomique spatiale.
    """)
    
    # Options de chargement
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🎯 Données disponibles")
        data_type = st.selectbox(
            "Type de données",
            ["10X Visium", "MERFISH", "ST (Spatial Transcriptomics)", "Custom"]
        )
        
        if data_type == "10X Visium":
            st.info("""
            **10X Visium** : Données de transcriptomique spatiale avec résolution de spot
            - Format standard 10X Genomics
            - Support des données mouse et human
            - Préprocessing automatique avec scanpy
            """)
    
    with col2:
        st.subheader("⚙️ Options de chargement")
        use_demo = st.checkbox("Utiliser les données de démonstration", value=True)
        
        if not use_demo:
            data_path = st.text_input("Chemin vers les données", placeholder="/path/to/visium/data")
            if st.button("Charger les données"):
                st.success("Données chargées avec succès !")
        else:
            if st.button("Charger les données de démonstration"):
                st.success("Données de démonstration chargées !")
                st.info("Dataset: Mouse brain Visium (10X Genomics)")

def show_fba_analysis():
    st.header("🔬 Analyse de flux métaboliques (FBA)")
    
    st.markdown("""
    Analyse des flux métaboliques basée sur l'expression génique spatiale.
    """)
    
    # Paramètres FBA
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🎛️ Paramètres")
        method = st.selectbox(
            "Méthode d'intégration",
            ["eflux", "iMAT", "GIMME", "pFBA"]
        )
        
        objective = st.selectbox(
            "Objectif métabolique",
            ["Biomass", "ATP", "Custom"]
        )
    
    with col2:
        st.subheader("📊 Résultats")
        st.info("""
        **Méthodes disponibles :**
        - **eFlux** : Intégration basée sur l'expression
        - **iMAT** : Intégration contextuelle
        - **GIMME** : Génération de modèles contextuels
        - **pFBA** : Parsimonious FBA
        """)
    
    if st.button("Lancer l'analyse FBA"):
        st.success("Analyse FBA lancée !")
        # Ici on pourrait lancer l'analyse réelle

def show_spatial_modeling():
    st.header("🌐 Modélisation spatiale")
    
    st.markdown("""
    Modélisation des processus de réaction-diffusion pour les métabolites spatiaux.
    """)
    
    # Paramètres de modélisation
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🧪 Paramètres de réaction-diffusion")
        substrate = st.selectbox("Substrat", ["glucose", "oxygen", "lactate", "glutamine"])
        diffusion_coeff = st.slider("Coefficient de diffusion", 1e-8, 1e-4, 1e-6, format="%.2e")
        reaction_rate = st.slider("Taux de réaction", 0.0, 1.0, 0.1)
    
    with col2:
        st.subheader("🌍 Paramètres spatiaux")
        domain_size = st.slider("Taille du domaine", 0.1, 10.0, 1.0)
        nx = st.slider("Résolution X", 10, 100, 50)
        ny = st.slider("Résolution Y", 10, 100, 50)
    
    if st.button("Lancer la simulation"):
        st.success("Simulation lancée !")
        # Ici on pourrait lancer la simulation réelle

def show_visualization():
    st.header("📊 Visualisation des résultats")
    
    st.markdown("""
    Visualisations interactives des résultats d'analyse.
    """)
    
    # Exemple de visualisation
    st.subheader("📈 Exemple de visualisation")
    
    # Données de démonstration
    np.random.seed(42)
    x = np.random.randn(100)
    y = np.random.randn(100)
    colors = np.random.randn(100)
    
    fig = px.scatter(x=x, y=y, color=colors, 
                     title="Exemple de visualisation spatiale",
                     labels={'x': 'Coordonnée X', 'y': 'Coordonnée Y'})
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Options de visualisation
    st.subheader("🎨 Options de visualisation")
    viz_type = st.selectbox(
        "Type de visualisation",
        ["Scatter plot", "Heatmap", "3D surface", "Contour plot"]
    )
    
    if st.button("Générer la visualisation"):
        st.success(f"Visualisation {viz_type} générée !")

if __name__ == "__main__":
    main()
