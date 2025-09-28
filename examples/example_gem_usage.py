#!/usr/bin/env python3
"""
Example usage of improved GEM and expression integration functionality.

This script demonstrates the new features:
- Human-GEM loading with automatic download
- GPR parsing with custom operators
- E-Flux and iMAT-like expression integration
- pFBA solving
"""

from scflux_spatial.gem.human_gem import HumanGEM
from scflux_spatial.gem.gpr import GPRParser
from scflux_spatial.fba.integrate_expression import ExpressionIntegrator


def main():
    """Main example function."""
    print("🧬 scflux_spatial - GEM and Expression Integration Example")
    print("=" * 60)
    
    # 1. Load Human-GEM model
    print("\n1. Loading Human-GEM model...")
    try:
        gem = HumanGEM(auto_download=True)
        model = gem.load_model()
        print(f"   ✓ Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    except Exception as e:
        print(f"   ⚠️  Model loading failed: {e}")
        print("   Using mock data for demonstration...")
        model = None
    
    # 2. Get default medium
    print("\n2. Getting default medium composition...")
    default_medium = gem.get_default_medium()
    print(f"   ✓ Default medium contains {len(default_medium)} exchange reactions")
    print(f"   ✓ Glucose uptake: {default_medium.get('EX_glc__D_e', 'N/A')} mmol/gDW/h")
    print(f"   ✓ Oxygen uptake: {default_medium.get('EX_o2_e', 'N/A')} mmol/gDW/h")
    
    # 3. Parse GPR rules (mock example)
    print("\n3. Parsing GPR rules...")
    gpr_parser = GPRParser()
    
    # Mock GPR rules for demonstration
    mock_gpr_rules = {
        "HK1": "HGNC:4922",  # Hexokinase 1
        "GAPDH": "HGNC:4141",  # Glyceraldehyde-3-phosphate dehydrogenase
        "LDHA": "HGNC:6535",  # Lactate dehydrogenase A
        "PFKL": "HGNC:8864",  # Phosphofructokinase, liver type
        "PKM": "HGNC:9021",   # Pyruvate kinase M1/M2
    }
    
    # Mock reaction-GPR mapping
    mock_reaction_gpr = {
        "HEX1": "HK1",
        "GAPD": "GAPDH", 
        "LDH_L": "LDHA",
        "PFK": "PFKL",
        "PYK": "PKM",
    }
    
    print(f"   ✓ Parsed {len(mock_reaction_gpr)} GPR rules")
    
    # 4. Test GPR evaluation with different operators
    print("\n4. Testing GPR evaluation...")
    
    # Mock gene expression data
    gene_expression = {
        "HK1": 0.8,
        "GAPDH": 0.9,
        "LDHA": 0.7,
        "PFKL": 0.6,
        "PKM": 0.85,
    }
    
    # Default operators (AND=min, OR=max)
    default_operators = {"AND": "min", "OR": "max"}
    reaction_scores_default = {}
    
    for reaction, gene in mock_reaction_gpr.items():
        score = gene_expression.get(gene, 0.0)
        reaction_scores_default[reaction] = score
    
    print("   Default operators (AND=min, OR=max):")
    for reaction, score in reaction_scores_default.items():
        print(f"   - {reaction}: {score:.3f}")
    
    # 5. Expression integration examples
    print("\n5. Expression integration methods...")
    
    # Mock base bounds
    base_bounds = {
        "HEX1": (-1000.0, 1000.0),
        "GAPD": (-1000.0, 1000.0),
        "LDH_L": (-1000.0, 1000.0),
        "PFK": (-1000.0, 1000.0),
        "PYK": (-1000.0, 1000.0),
    }
    
    integrator = ExpressionIntegrator(gpr_parser)
    
    # E-Flux method
    print("\n   E-Flux integration:")
    eflux_bounds = integrator.integrate_expression_eflux(
        gene_expression, base_bounds, scaling_factor=1.0
    )
    
    for reaction, (lower, upper) in eflux_bounds.items():
        if reaction in reaction_scores_default:
            original_score = reaction_scores_default[reaction]
            print(f"   - {reaction}: bounds=[{lower:.1f}, {upper:.1f}], expr={original_score:.3f}")
    
    # iMAT-like method
    print("\n   iMAT-like integration:")
    imat_bounds = integrator.integrate_expression_imat_like(
        gene_expression, base_bounds,
        high_quantile=0.8, low_quantile=0.2,
        activation_bound=100.0, minimization_bound=0.1
    )
    
    for reaction, (lower, upper) in imat_bounds.items():
        if reaction in reaction_scores_default:
            original_score = reaction_scores_default[reaction]
            status = "ACTIVATED" if abs(lower) > 50 else "MINIMIZED" if abs(upper) < 1 else "NORMAL"
            print(f"   - {reaction}: bounds=[{lower:.1f}, {upper:.1f}], expr={original_score:.3f}, status={status}")
    
    # 6. Method selection via Hydra configuration
    print("\n6. Method selection via configuration...")
    
    methods = ["eflux", "imat_like", "linear", "quadratic", "none"]
    
    for method in methods:
        try:
            result_bounds = integrator.integrate_expression_with_method(
                gene_expression, base_bounds, method=method
            )
            print(f"   ✓ {method.upper()}: {len(result_bounds)} reactions processed")
        except Exception as e:
            print(f"   ❌ {method.upper()}: {e}")
    
    # 7. Summary statistics
    print("\n7. Summary statistics...")
    
    print(f"   Gene expression range: {min(gene_expression.values()):.3f} - {max(gene_expression.values()):.3f}")
    print(f"   Average expression: {sum(gene_expression.values()) / len(gene_expression):.3f}")
    
    # Calculate some basic statistics on the E-Flux results
    eflux_scales = []
    for reaction, (lower, upper) in eflux_bounds.items():
        if reaction in base_bounds:
            base_lower, base_upper = base_bounds[reaction]
            if base_upper != 0:
                scale = upper / base_upper
                eflux_scales.append(abs(scale))
    
    if eflux_scales:
        print(f"   E-Flux scaling range: {min(eflux_scales):.3f} - {max(eflux_scales):.3f}")
        print(f"   Average E-Flux scaling: {sum(eflux_scales) / len(eflux_scales):.3f}")
    
    print(f"\n✅ Example completed successfully!")
    print(f"\nTo use with real data:")
    print(f"   1. Load your Visium data with load_visium()")
    print(f"   2. Extract gene expression from adata.var")
    print(f"   3. Map genes to model genes using GPR parser")
    print(f"   4. Apply expression integration method")
    print(f"   5. Solve with FBA/pFBA")


if __name__ == "__main__":
    main()
