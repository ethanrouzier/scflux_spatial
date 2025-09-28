"""
Spatial FBA coupling with reaction-diffusion equations.
"""

from typing import Callable, Dict
import numpy as np
from .units import flux_to_volumetric_source


class SpatialFBACoupler:
    """
    Orchestrates operator-splitting between RD and FBA.

    Args
    ----
    rd: RDSolver-like avec .snapshot(), .set_source(name, arr), .step(dt).
    fba_fn: callable conc_dict[str->np.ndarray] -> dict[str->np.ndarray] flux en mmol·gDW⁻¹·h⁻¹
    rho_gDW_per_m3: densité biomasse (gDW·m⁻³)
    alpha: under-relaxation (0..1)
    tol_rel: tolérance convergence (max norme relative)
    max_iters: borne de sécurité (non utilisé si piloté par extérieur)
    """
    def __init__(self, rd, fba_fn: Callable[[Dict[str, np.ndarray]], Dict[str, np.ndarray]],
                 rho_gDW_per_m3: float, alpha: float = 0.4, tol_rel: float = 1e-3, max_iters: int = 20):
        self.rd = rd
        self.fba_fn = fba_fn
        self.rho = float(rho_gDW_per_m3)
        self.alpha = float(alpha)
        self.tol = float(tol_rel)
        self.max_iters = int(max_iters)
        self._prev_fluxes = None
        self._prev_C = None

    def _relax(self, new, old):
        return self.alpha * new + (1.0 - self.alpha) * old

    def iterate(self, dt: float) -> Dict[str, np.ndarray]:
        # 1) Snapshot concentrations (mol·m⁻³)
        C = self.rd.snapshot()

        # 2) FBA -> fluxes (mmol·gDW⁻¹·h⁻¹)
        fluxes = {k: np.asarray(v) for k, v in self.fba_fn(C).items()}

        # 3) Under-relaxation des flux
        if self._prev_fluxes is not None:
            for k in fluxes:
                fluxes[k] = self._relax(fluxes[k], self._prev_fluxes[k])
        self._prev_fluxes = {k: v.copy() for k, v in fluxes.items()}

        # 4) Conversion en sources volumiques (mol·m⁻³·s⁻¹) + clip de sécurité
        for met, v in fluxes.items():
            S = flux_to_volumetric_source(v, self.rho)  # broadcast
            S = np.clip(S, -1e3, 1e3)
            self.rd.set_source(met, S)

        # 5) Pas RD
        self.rd.step(dt)

        # 6) Convergence (max norme relative L2)
        C_new = self.rd.snapshot()
        if self._prev_C is None:
            self._prev_C = {k: v.copy() for k, v in C_new.items()}
            return C_new

        # non-négativité soft
        try:
            for name, var in getattr(self.rd, "vars", {}).items():
                var.setValue(np.maximum(var.value, 0.0))
        except Exception:
            # si pas d'API .vars, ignorer: on assume clamp côté RD
            pass

        rels = []
        for k in C_new:
            num = np.linalg.norm(C_new[k] - self._prev_C[k])
            den = np.linalg.norm(self._prev_C[k]) + 1e-12
            rels.append(num / den)
            self._prev_C[k] = self._relax(C_new[k], self._prev_C[k])

        return C_new  # le pilote externe lit rels via l'historique (voir notebook)