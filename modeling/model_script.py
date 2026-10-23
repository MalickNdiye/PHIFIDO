import numpy as np
import scipy.integrate as integrate
from numba import njit,prange
import time
import os

SEED = int(os.environ.get("SEED", 2024))
@njit
def seed_numba(s):
    np.random.seed(s)   
seed_numba(SEED)

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

    # Diversity from last 10%
    last_10 = abundances[:, -cutoff:]
    mean_abund = np.mean(last_10, axis=1)
    gamma_div = diversity_measures(mean_abund[Np:Np+Ns], survthreshold)
    alpha_div = np.mean([diversity_measures(last_10[Np:Np+Ns, i], survthreshold) for i in range(cutoff)], axis=0)
    assert alpha_div.shape == (3,)

    # Early detrended CV (first 10%)
    first_10 = abundances[:, :cutoff]
    n_tp = first_10.shape[1]
    survived_early = first_10[Np:Np+Ns, :].mean(axis=1) > survthreshold

    early_dcv = 0.0
    if np.sum(survived_early) > 0 and n_tp >= 3:
        t_idx = np.arange(n_tp, dtype=np.float64)
        cvs = []
        for j in np.where(survived_early)[0]:
            traj = first_10[Np + j, :]
            trend = np.polyval(np.polyfit(t_idx, traj, 2), t_idx)
            m = np.mean(traj)
            if m > 1e-30:
                cvs.append(np.std(traj - trend) / m)
        if len(cvs) > 0:
            early_dcv = np.mean(cvs)

    # Late detrended CV (last 10%)
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

    # [richness, shannon, hill, early_dcv, late_dcv]
    return np.array([alpha_div[0], alpha_div[1], alpha_div[2], early_dcv, late_dcv])

def integrate_one(Np, Ns, Nr, etaMatrix, muMatrix, supplyVec, delta, beta,
                  t_eval, has_phage):
    if has_phage:
        ic_phage = np.full(Np, 0.1)
    else:
        ic_phage = np.full(Np, 1e-30)
    ic = np.log(np.concatenate((ic_phage, np.full(Ns, 1.0), supplyVec)))

    sol = integrate.solve_ivp(
        dynamics, (0, t_eval[-1]), ic,
        args=(etaMatrix, muMatrix, supplyVec, delta, beta, Np, Ns, Nr),
        t_eval=t_eval, method='Radau', rtol=1e-4, atol=1e-9
    )
    abundance = np.exp(sol.y)

    if sol.success:
        return compute_stats(abundance, Np, Ns, Nr, survthreshold=1e-2),np.mean(abundance[Np:Np+Ns, -abundance.shape[1]//10:], axis=1)
    else:
        return np.full(5, np.nan),np.nan

Np = 24
Ns = 12
Nr = 30

delta = 0.1
totalCarbon = 5
beta = 10
t = np.linspace(0, 1000, 100)

def run_trial(seed, Nr,t_eval,etaMatrix,muMatrix_sugar, muMatrix_pollen):
    np.random.seed(seed)
    # Sugar water: dense mu, Nr//5 resources, fixed carbon
    n_sugar = max(Nr // 5, 2)
    supplyVec_sugar = np.zeros(Nr)
    top_res_sugar = np.argsort(muMatrix_sugar.sum(axis=0))[::-1][:n_sugar]
    supplyVec_sugar[top_res_sugar] = totalCarbon * Nr / n_sugar
    supplyVec_sugar[supplyVec_sugar == 0] = 1e-30

    # Pollen: variable mu sparsity, all Nr resources, variable carbon
    supplyVec_pollen = np.full(Nr, totalCarbon)

    # Run all 4 conditions
    s_no, s_no_mean = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_sugar, supplyVec_sugar,
                         delta, beta, t_eval, False)
    s_ph, s_ph_mean = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_sugar, supplyVec_sugar,
                         delta, beta, t_eval, True)
    p_no, p_no_mean = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_pollen, supplyVec_pollen,
                         delta, beta, t_eval, False)
    p_ph, p_ph_mean = integrate_one(Np, Ns, Nr, etaMatrix, muMatrix_pollen, supplyVec_pollen,
                         delta, beta, t_eval, True)

    # Return 20 values: [sugar_no(5), sugar_ph(5), pollen_no(5), pollen_ph(5)]
    return (s_no, s_ph, p_no, p_ph), (s_no_mean, s_ph_mean, p_no_mean, p_ph_mean)

def run_communities(num_rums,shape_sugar, shape_pollen, mmg, mmi, corr,metabolic_flag=False):
    combined_stats = np.zeros((num_rums, 4,5))
    combined_traits = np.zeros((num_rums, Ns, 2))  # Store max growth and infectivity for each run
    combined_eta = np.zeros((num_rums, Np, Ns))
    combined_abundances = np.zeros((num_rums, 4, Ns))  # Store full trajectories for each run

    for runID in range(num_rums):
        max_growth_rate, max_infectivity = generate_traits(Ns, mmg, mmi, corr)
        combined_traits[runID, :, 0] = max_growth_rate
        combined_traits[runID, :, 1] = max_infectivity

        if not metabolic_flag:
            muMatrix_sugar = generate_mu(Ns, Nr, max_growth_rate, shape_sugar)
            muMatrix_pollen = generate_mu(Ns, Nr, max_growth_rate, shape_pollen)
        else:
            muMatrix_sugar = generate_mu(Ns, Nr, max_growth_rate, shape_sugar)
            muMatrix_pollen = muMatrix_sugar.copy()
            
        etaMatrix = generate_eta(Np, Ns, max_infectivity)
        combined_eta[runID] = etaMatrix

        stats_community, abundances_community = run_trial(0, Nr, t, etaMatrix, muMatrix_sugar, muMatrix_pollen)
        combined_stats[runID] = stats_community
        combined_abundances[runID] = abundances_community

    return combined_stats, combined_traits, combined_eta, combined_abundances


###########################################################
################### SINGLE SHAPE VALUE ####################
###########################################################

shape_sugar = 2.0
shape_pollen = 0.25
mmg = -1.0
mmi = -3.0
corr = 0.0

num_runs = 100

combined_stats, combined_traits, combined_eta, combined_abundances = run_communities(num_runs, shape_sugar, shape_pollen, mmg, mmi, corr,metabolic_flag=False)

np.save("model_results/combined_stats.npy", combined_stats)
print("Saved combined_stats.npy for a single shape value.")

###########################################################
################### ACROSS SHAPE VALUES ###################
###########################################################

num_runs = 50 ## reduced for speed
shape_pollen_values = np.geomspace(0.1, 2, 20)

across_shape_stats = np.zeros((len(shape_pollen_values), num_runs,4, 5))

for i, shape_pollen in enumerate(shape_pollen_values):
    combined_stats_sparse, _, _,_ = run_communities(num_runs, shape_sugar, shape_pollen, mmg, mmi, corr,metabolic_flag=False)
    across_shape_stats[i] = combined_stats_sparse
    print(f"Completed shape_pollen {shape_pollen:.3f} ({i+1}/{len(shape_pollen_values)})")

np.savez("model_results/across_shape_stats.npz", shape_pollen_values=shape_pollen_values, across_shape_stats=across_shape_stats)

def expected_niche_overlap(num_samples, shape_pollen, mmg, mmi):
    niche_overlaps = np.zeros(num_samples)
    for runID in range(num_samples):
        max_growth_rate, _ = generate_traits(Ns, mmg, mmi, 0)
        muMatrix = generate_mu(Ns, Nr, max_growth_rate, shape_pollen)
        norms = np.sqrt(np.sum(muMatrix**2, axis=1))
        cosine_matrix = (muMatrix @ muMatrix.T) / np.outer(norms, norms)
        mask = ~np.eye(Ns, dtype=bool)
        niche_overlaps[runID] = np.median(cosine_matrix[mask])
    return niche_overlaps

expected_niche_overlap_values = np.zeros(len(shape_pollen_values))
for i, sp in enumerate(shape_pollen_values):
    expected_niche_overlap_values[i] = np.mean(
        expected_niche_overlap(10000, sp, mmg, mmi)
    )

delta_pollen_combined = np.zeros((len(shape_pollen_values), num_runs))
for i in range(len(shape_pollen_values)):
    delta_pollen_combined[i] = across_shape_stats[i, :, 3, 2] - across_shape_stats[i, :, 2, 2]

np.savez("model_results/delta_pollen_values.npz",
         shape_pollen_values=shape_pollen_values,
         delta_pollen_combined=delta_pollen_combined,
         expected_niche_overlap_values=expected_niche_overlap_values)

print("Saved delta_pollen_values.npz with delta pollen and expected niche overlap values.")

