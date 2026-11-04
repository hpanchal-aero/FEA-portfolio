#!/usr/bin/python3
"""Ti-6Al-4V (Grade 5, annealed) property set for Project 05.

Single source of truth: every script imports MATERIAL from here.
Units: SI (m, kg, s, N, Pa, W, J). Temperatures in degC, stress-free
reference temperature T_REF_C = 20 degC (the CTE below is the 20-100 degC mean).

Sources (accessed 2026-11-01; typical vendor/aggregator values, NOT design
allowables; MMPDS/ASM handbooks were not available):
  rho, k, E, alpha, sigma_y : https://www.ulbrich.com/alloys/ti-6-al-4-v-uns-r56400/
  cp, nu, alpha(T) trend    : https://parts-badger.com/properties-of-grade-5-titanium
  cross-checks              : https://www.makeitfrom.com/material-properties/Grade-5-Ti-6Al-4V-3.7165-R56400-Titanium
                              https://www.voestalpine.com/specialtymetals/eu/app/uploads/sites/251/2026/03/Titanium-Grade-5-Ti6Al4V.pdf

Constant properties (decision P05-D2). Real k and alpha rise with temperature,
E and sigma_y fall; results are valid for a moderate temperature rise only.
"""

T_REF_C = 20.0

MATERIAL = {
    "name": "Ti-6Al-4V (Grade 5, annealed, typical)",
    "E": 113.8e9,        # Pa, Young's modulus
    "nu": 0.33,          # -, Poisson's ratio
    "alpha": 8.6e-6,     # 1/K, mean CTE 20-100 degC
    "rho": 4430.0,       # kg/m3
    "cp": 526.0,         # J/(kg K)
    "k": 6.7,            # W/(m K) at 21 degC
    "sigma_y": 880.0e6,  # Pa, 0.2 % yield, annealed typical, room temperature
}


def derived(m=MATERIAL):
    """Quantities used by the analytical checks and for pulse design."""
    diffusivity = m["k"] / (m["rho"] * m["cp"])                    # m2/s
    sigma_per_K = m["E"] * m["alpha"] / (1.0 - m["nu"])            # Pa/K, biaxially restrained
    return {
        "diffusivity": diffusivity,
        "rho_cp": m["rho"] * m["cp"],
        "sigma_per_K_restrained_biaxial": sigma_per_K,
        "dT_to_yield_restrained_biaxial": m["sigma_y"] / sigma_per_K,
    }


if __name__ == "__main__":
    m = MATERIAL
    d = derived(m)
    print(m["name"])
    print("T_ref            %8.1f degC" % T_REF_C)
    print("E                %8.2f GPa" % (m["E"] / 1e9))
    print("nu               %8.3f" % m["nu"])
    print("alpha            %8.2e 1/K" % m["alpha"])
    print("rho              %8.1f kg/m3" % m["rho"])
    print("cp               %8.1f J/(kg K)" % m["cp"])
    print("k                %8.2f W/(m K)" % m["k"])
    print("sigma_y          %8.1f MPa" % (m["sigma_y"] / 1e6))
    print("rho*cp           %8.3e J/(m3 K)" % d["rho_cp"])
    print("diffusivity      %8.4e m2/s" % d["diffusivity"])
    print("E*alpha/(1-nu)   %8.4f MPa/K (restrained biaxial stress per K)" % (d["sigma_per_K_restrained_biaxial"] / 1e6))
    print("dT to yield      %8.1f K (restrained biaxial, constant properties)" % d["dT_to_yield_restrained_biaxial"])
    for t in (1e-3, 2e-3, 3e-3):
        print("t = %.0f mm: t^2/diffusivity = %.3f s (through-thickness diffusion time scale)" % (t * 1e3, t * t / d["diffusivity"]))
