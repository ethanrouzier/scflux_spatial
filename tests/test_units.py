"""
Unit tests for unit conversion utilities.
"""

from scflux_spatial.spatial.units import flux_to_volumetric_source, volumetric_source_to_flux


def test_flux_conversion():
    """Test flux to volumetric source conversion."""
    # 1 mmol/gDW/h at 1000 gDW/m^3 -> 1.0 / 3600.0 mol/m^3/s
    expected = 1.0 / 3600.0
    result = flux_to_volumetric_source(1.0, 1000.0)
    assert abs(result - expected) < 1e-12


def test_volumetric_source_to_flux():
    """Test volumetric source to flux conversion."""
    # Round-trip test
    original_flux = 1.0  # mmol/gDW/h
    rho = 1000.0  # gDW/m^3
    
    source = flux_to_volumetric_source(original_flux, rho)
    recovered_flux = volumetric_source_to_flux(source, rho)
    
    assert abs(recovered_flux - original_flux) < 1e-12


def test_negative_flux():
    """Test negative flux (consumption)."""
    # Negative flux should give negative source
    flux = -0.5  # mmol/gDW/h (consumption)
    rho = 2000.0  # gDW/m^3
    
    source = flux_to_volumetric_source(flux, rho)
    assert source < 0  # Consumption should be negative
    assert abs(source - (-0.5 * 2000.0 * 1e-3 / 3600.0)) < 1e-12
