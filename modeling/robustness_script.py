import numpy as np
import scipy.integrate as integrate
from numba import njit
import os
import time

SEED = int(os.environ.get("SEED", 2024))
@njit
def seed_numba(s):
    np.random.seed(s)   
seed_numba(SEED)

START_TIME = time.time()

@njit
def dynamics(t, logy, etaMatrix, muMatrix, supplyVec, delta, beta, Np, Ns, Nr):
    logphages = logy[:Np]
    logbacteria = logy[Np:Np+Ns]
    logresources = logy[Np+Ns:]
    uptake = muMatrix * np.exp(logresources)
    dlogphage = beta * np.sum(etaMatrix * np.exp(logbacteria), axis=1) - delta
    dlogbacteria = np.sum(uptake, axis=1) - np.sum(etaMatrix.T * np.exp(logphages), axis=1) - delta
    dlogresources = delta * (supplyVec * np.exp(-logresources) - 1.0) - np.sum(muMatrix.T * np.exp(logbacteria), axis=1)
    dlogphage[logphages <= -35] = 0
    dlogbacteria[logbacteria <= -35] = 0
    dlogresources[logresources <= -35] = 0
    return np.concatenate((dlogphage, dlogbacteria, dlogresources))

@njit
def diversity_measures(populations, survthreshold):
    pops = populations.copy()
    pops[pops < 1e-15] = 1e-15
    richness = np.sum(pops > survthreshold)
    freqs = pops / np.sum(pops)
    if np.sum(pops) < 1e-12:
        shannon = -1.0
        hill = -1.0
    else:
        shannon = -np.sum(freqs * np.log(freqs))
        hill = np.exp(shannon)
    return np.array([richness, shannon, hill])

@njit
def generate_eta(Np, Ns, max_infectivity):
    etaMatrix = np.random.gamma(0.5, 2.5, (Np, Ns))
    for i in range(Ns):
        etaMatrix[:, i] = etaMatrix[:, i] / np.sum(etaMatrix[:, i]) * max_infectivity[i]
    return etaMatrix

@njit
def generate_mu(Ns, Nr, max_growth_rate, mu_shape):
    muMatrix = np.random.gamma(mu_shape, 2.5, (Ns, Nr))
    for i in range(Ns):
        muMatrix[i, :] = muMatrix[i, :] / np.sum(muMatrix[i, :]) * max_growth_rate[i]
    return muMatrix

@njit
def generate_traits(Ns, mmg, mmi, corr):
    cov = np.array([[1, corr], [corr, 1]]) * 0.5**2
    L = np.linalg.cholesky(cov)
    Z = np.random.randn(Ns, 2)
    cv = Z @ L.T + np.array([mmg, mmi])
    max_growth_rate = np.exp(cv[:, 0])
    max_infectivity = np.exp(cv[:, 1])
    return max_growth_rate, max_infectivity

def compute_stats(abundances, Np, Ns, Nr, survthreshold):
    n_t = abundances.shape[1]
    cutoff = max(n_t // 10, 2)

    last_10 = abundances[:, -cutoff:]
    mean_abund = np.mean(last_10, axis=1)
    gamma_div = diversity_measures(mean_abund[Np:Np+Ns], survthreshold)
    alpha_div = np.mean(
        [diversity_measures(last_10[Np:Np+Ns, i], survthreshold) for i in range(cutoff)],
        axis=0
    )

    # Late detrended CV (oscillation score)
    endpoint = abundances[Np:Np+Ns, -1]
    survived_late = endpoint > survthreshold
    late_dcv = 0.0
    if np.sum(survived_late) > 0 and cutoff >= 3:
        t_idx = np.arange(cutoff, dtype=np.float64)
        cvs = []
        for j in np.where(survived_late)[0]:
            traj = last_10[Np + j, :]
            trend = np.polyval(np.polyfit(t_idx, traj, 2), t_idx)
            m = np.mean(traj)
            if m > 1e-30:
                cvs.append(np.std(traj - trend) / m)
        if len(cvs) > 0:
            late_dcv = np.mean(cvs)

    return np.array([alpha_div[0], alpha_div[1], alpha_div[2], 0.0, late_dcv])

def integrate_one(Np, Ns, Nr, etaMatrix, muMatrix, supplyVec, delta, beta,
                  t_eval, has_phage):
    if has_phage:
        ic_phage = np.full(Np, 0.1)
    else:
        ic_phage = np.full(Np, 1e-30)
    ic = np.log(np.concatenate((ic_phage, np.full(Ns, 1.0), supplyVec)))

    try:
        sol = integrate.solve_ivp(
            dynamics, (0, t_eval[-1]), ic,
            args=(etaMatrix, muMatrix, supplyVec, delta, beta, Np, Ns, Nr),
            t_eval=t_eval, method='Radau', rtol=1e-4, atol=1e-9
        )

        if sol.success:
            abundance = np.exp(sol.y)
            stats = compute_stats(abundance, Np, Ns, Nr, survthreshold=1e-2)
            mean_abund = np.mean(abundance[Np:Np+Ns, -abundance.shape[1]//10:], axis=1)
            return stats, mean_abund
        else:
            return np.full(5, np.nan), np.full(Ns, np.nan)
    except (ValueError, FloatingPointError):
        return np.full(5, np.nan), np.full(Ns, np.nan)


def run_single_community(Np, Ns, Nr, etaMatrix, muMatrix_sugar, muMatrix_pollen,
                         supplyVec_sugar, supplyVec_pollen, delta, beta, t_eval):
    """Run all 4 conditions for one community with fully specified supply vectors."""
    s_no, _ = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_sugar, supplyVec_sugar,
                            delta, beta, t_eval, False)
    s_ph, _ = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_sugar, supplyVec_sugar,
                            delta, beta, t_eval, True)
    p_no, _ = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_pollen, supplyVec_pollen,
                            delta, beta, t_eval, False)
    p_ph, _ = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_pollen, supplyVec_pollen,
                            delta, beta, t_eval, True)
    return np.array([s_no, s_ph, p_no, p_ph])


def make_supply_vectors(Nr, Nr_sugar, muMatrix_sugar, totalCarbon_sugar, totalCarbon_pollen):
    """Construct supply vectors for sugar and pollen."""
    supplyVec_sugar = np.zeros(Nr)
    top_res = np.argsort(muMatrix_sugar.sum(axis=0))[::-1][:Nr_sugar]
    supplyVec_sugar[top_res] = totalCarbon_sugar * Nr / Nr_sugar
    supplyVec_sugar[supplyVec_sugar == 0] = 1e-30

    supplyVec_pollen = np.full(Nr, totalCarbon_pollen)
    return supplyVec_sugar, supplyVec_pollen


def run_robustness(num_runs, Np, Ns, Nr, Nr_sugar,
                   shape_sugar, shape_pollen,
                   mmg, mmi, corr,
                   delta, beta,
                   totalCarbon_sugar, totalCarbon_pollen,
                   t_eval, label=""):
    """
    General-purpose runner for robustness checks.
    Returns: combined_stats array of shape (num_runs, 4, 5)
             4 conditions: [SW_no, SW_ph, PE_no, PE_ph]
             5 stats: [richness, shannon, hill, early_dcv, late_dcv]
    """
    combined_stats = np.zeros((num_runs, 4, 5))

    for runID in range(num_runs):
        if runID % 10 == 0 and label:
            print(f"  {label}: run {runID}/{num_runs}")

        max_growth_rate, max_infectivity = generate_traits(Ns, mmg, mmi, corr)
        etaMatrix = generate_eta(Np, Ns, max_infectivity)
        muMatrix_sugar = generate_mu(Ns, Nr, max_growth_rate, shape_sugar)
        muMatrix_pollen = generate_mu(Ns, Nr, max_growth_rate, shape_pollen)

        supplyVec_sugar, supplyVec_pollen = make_supply_vectors(
            Nr, Nr_sugar, muMatrix_sugar, totalCarbon_sugar, totalCarbon_pollen
        )

        combined_stats[runID] = run_single_community(
            Np, Ns, Nr, etaMatrix, muMatrix_sugar, muMatrix_pollen,
            supplyVec_sugar, supplyVec_pollen, delta, beta, t_eval
        )

    return combined_stats

## Default parametrs

Np = 24
Ns = 12
Nr = 30
Nr_sugar_default = 6
delta = 0.1
beta_default = 10
totalCarbon_default = 5.0
mmg_default = -1
mmi_default = -3
corr_default = 0
shape_sugar_default = 2
shape_pollen_default = 0.25
t_eval = np.linspace(0, 1000, 100)
NUM_RUNS = 100

os.makedirs("robustness_results", exist_ok=True)


CHECKS_TO_RUN = [1,2,3,4,5] 

## Check 1 - checks if smaller number of resources supplies for sucrose changes or if niche overlap is the real driver

def run_check_1():
    print("CHECK 1: Nr × niche overlap")
    configs = [(6, 0.25, "Nr6_lowOverlap"), (6, 2.0, "Nr6_highOverlap"),
               (30, 0.25, "Nr30_lowOverlap"), (30, 2.0, "Nr30_highOverlap")]
    results = {}
    for Nr_sugar, shape_pollen, label in configs:
        print(f"  {label}")
        results[label] = run_robustness(
            NUM_RUNS, Np, Ns, Nr, Nr_sugar, shape_sugar_default, shape_pollen,
            mmg_default, mmi_default, corr_default, delta, beta_default,
            totalCarbon_default, totalCarbon_default, t_eval, label=label)
    np.savez("robustness_results/check1_Nr_x_overlap.npz", **results)

## Check 2 - what if maximum growth rate was correlated with maximum infectivity (or the reverse)?

def run_check_2():
    print("CHECK 2: Growth-susceptibility correlation")
    results = {}
    for c in [-0.5, 0.0, 0.5, 0.9]:
        label = f"corr_{c:.1f}"
        print(f"  {label}")
        results[label] = run_robustness(
            NUM_RUNS, Np, Ns, Nr, Nr_sugar_default, shape_sugar_default,
            shape_pollen_default, mmg_default, mmi_default, c, delta, beta_default,
            totalCarbon_default, totalCarbon_default, t_eval, label=label)
    results["corr_values"] = np.array([-0.5, 0.0, 0.5, 0.9])
    np.savez("robustness_results/check2_correlation.npz", **results)

## Check 3 - Robustness against phage burst size (beta)

def run_check_3():
    print("CHECK 3: Burst size")
    results = {}
    for b in [5, 10, 50, 100]:
        label = f"beta_{b}"
        print(f"  {label}")
        results[label] = run_robustness(
            NUM_RUNS, Np, Ns, Nr, Nr_sugar_default, shape_sugar_default,
            shape_pollen_default, mmg_default, mmi_default, corr_default, delta, b,
            totalCarbon_default, totalCarbon_default, t_eval, label=label)
    results["beta_values"] = np.array([5, 10, 50, 100])
    np.savez("robustness_results/check3_burst_size.npz", **results)

## Check 4 - Robustness against different maximum growth rates and maximum infectivities

def run_check_4():
    print("CHECK 4: Trait distributions (mmg × mmi)")
    results = {}
    for mg in [-2.0, -1.0, 0.0]:
        for mi in [-4.0, -3.0, -2.0]:
            label = f"mmg_{mg:.0f}_mmi_{mi:.0f}"
            print(f"  {label}")
            results[label] = run_robustness(
                NUM_RUNS, Np, Ns, Nr, Nr_sugar_default, shape_sugar_default,
                shape_pollen_default, mg, mi, corr_default, delta, beta_default,
                totalCarbon_default, totalCarbon_default, t_eval, label=label)
    results["mmg_values"] = np.array([-2.0, -1.0, 0.0])
    results["mmi_values"] = np.array([-4.0, -3.0, -2.0])
    np.savez("robustness_results/check4_trait_distributions.npz", **results)

## Check 5 - Does total supplied carbon against sucrose vs pollen extract matter?

def run_check_5():
    print("CHECK 5: Total carbon ratio")
    results = {}
    for cp in [2.5, 5.0, 10.0]:
        ratio = cp / totalCarbon_default
        label = f"carbon_ratio_{ratio:.1f}"
        print(f"  {label}")
        results[label] = run_robustness(
            NUM_RUNS, Np, Ns, Nr, Nr_sugar_default, shape_sugar_default,
            shape_pollen_default, mmg_default, mmi_default, corr_default, delta,
            beta_default, totalCarbon_default, cp, t_eval, label=label)
    results["carbon_pollen_values"] = np.array([2.5, 5.0, 10.0])
    np.savez("robustness_results/check5_carbon_ratio.npz", **results)

## run dict for the different checks

CHECK_FUNCTIONS = {
    1: run_check_1,   2: run_check_2,   3: run_check_3,
    4: run_check_4,   5: run_check_5
}

for check_num in CHECKS_TO_RUN:
    t0 = time.time()
    print("\n" + "=" * 60)
    CHECK_FUNCTIONS[check_num]()
    print(f"  Time: {time.time() - t0:.0f}s")

print(f"\nTotal time: {time.time() - START_TIME:.0f}s")
