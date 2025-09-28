"""
Unit conversion utilities for spatial flux analysis.
"""

def flux_to_volumetric_source(v_mmol_gDW_h, rho_gDW_per_m3):
    """
    Convertit un flux FBA (mmol·gDW⁻¹·h⁻¹) en source volumique (mol·m⁻³·s⁻¹).
    Positif = production, négatif = consommation.
    """
    return float(v_mmol_gDW_h) * float(rho_gDW_per_m3) * 1e-3 / 3600.0