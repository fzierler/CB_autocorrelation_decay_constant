import numpy as np
from scipy import linalg

from bootstrap import BootstrapSampleSet
import fitting

def extract_meson_mass(C_tmp, plateau_start, plateau_end):
    E_fit, A_fit, chisquare = fitting.fit_cosh_bootstrap(
        C_tmp, plateau_start, plateau_end
    )

    return E_fit, A_fit, round(chisquare, 2)


def meson_decay_constant(Css, Csp, plateau_start, plateau_end):
    # load the ensamble info
    lattice_t = np.shape(Css.mean)[1]

    E_fit, A_fit, chisquare = fitting.fit_coshsinh_simultaneous(
        Css, Csp, plateau_start, plateau_end, lattice_t
    )

    return E_fit, A_fit, round(chisquare, 2)
