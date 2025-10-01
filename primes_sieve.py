from mpmath import *
import numpy as np
from sympy import primerange, primepi, isprime, mobius
from math import exp, log, pi, sqrt
import time
import argparse
import os

mp.dps = 20

def riemann_R(x, K=50):
    s = mpf(0)
    for n in range(1, K+1):
        mu = mobius(n)
        if mu == 0:
            continue
        s += mu / n * li(x ** (1/n))
    return float(s)

def get_gammas_dynamic(num_zeros):
    return np.array([float(im(zetazero(k))) for k in range(1, num_zeros + 1)])

def segmented_pre_sieve(start_n, end_n, B):
    small_primes = list(primerange(2, B + 1))
    length = end_n - start_n + 1
    is_candidate = [True] * length
    for p in small_primes:
        start_multiple = max(p * p, ((start_n + p - 1) // p) * p)
        if start_multiple > end_n:
            continue
        idx = start_multiple - start_n
        for i in range(idx, length, p):
            is_candidate[i] = False
    return [start_n + i for i in range(length) if is_candidate[i]]

def compute_spectral_scores(candidates, gammas, h):
    log_ns = np.log(candidates)
    osc_plus = np.zeros(len(candidates), dtype=complex)
    osc_minus = np.zeros(len(candidates), dtype=complex)
    for gamma in gammas:
        taper = exp(-0.5 * h**2 * gamma**2)
        shift_plus = (0.5 + 1j * gamma) * h
        shift_minus = (0.5 + 1j * gamma) * (-h)
        phases_plus = np.exp(1j * gamma * log_ns + shift_plus)
        phases_minus = np.exp(1j * gamma * log_ns + shift_minus)
        osc_plus += taper * phases_plus / (0.5 + 1j * gamma)
        osc_minus += taper * phases_minus / (0.5 + 1j * gamma)
    psi_plus = np.array(candidates) * exp(h) - 2 * np.real(osc_plus)
    psi_minus = np.array(candidates) * exp(-h) - 2 * np.real(osc_minus)
    logn_arr = np.log(candidates)
    scores = (psi_plus - psi_minus) / (2 * h * np.array(candidates) * logn_arr)
    return scores

def chudnovsky_like_sieve(N, T=50, K=50, epsilon=1.2):
    gammas = get_gammas_dynamic(T)
    B = int(sqrt(N)) + 1
    mid = N / 2
    h = 0.05 / log(mid)
    approx = riemann_R(N, K)
    M = int(ceil(epsilon * approx))
    candidates = segmented_pre_sieve(2, N, B)
    scores = compute_spectral_scores(candidates, gammas, h)
    mean_s = np.mean(scores)
    sigma_s = max(np.std(scores), 0.01)
    z_scores = (scores - mean_s) / sigma_s
    top_idx = np.argsort(-z_scores)[:M]
    top_candidates = np.array(candidates)[top_idx]
    primes = [int(c) for c in top_candidates if isprime(int(c))]
    # Check
    expected_lower = approx - sqrt(N)
    expected_upper = approx + sqrt(N)
    num_primes = len(primes)
    if num_primes < expected_lower or num_primes > expected_upper:
        print(f"Warning: Found {num_primes}, expected ~{approx}")
        # Refine: e.g., increase K to 100
        # For demo, skip
    return sorted(primes)

def validate_primes(predicted, start_n, end_n):
    """Enhanced validation: Compute metrics and analyze missed primes."""
    true_primes = set(primerange(start_n, end_n + 1))
    predicted_set = set(predicted)
    TP = len(predicted_set & true_primes)
    FP = len(predicted_set - true_primes)
    FN = len(true_primes - predicted_set)
    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall = TP / (TP + FN) if TP + FN > 0 else 0
    missed = sorted(true_primes - predicted_set)[:10]  # Top 10 missed
    gaps = [missed[i+1] - missed[i] for i in range(len(missed)-1)] if len(missed) > 1 else []
    return precision, recall, missed, gaps
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Chudnovsky-like prime sieve demo")
    parser.add_argument("--n", type=int, default=10000, help="Upper bound N (inclusive) for prime search")
    parser.add_argument("--output", type=str, default=None, help="Write generated primes to this file (one per line)")
    args = parser.parse_args()

    N = args.n
    if N < 2:
        raise SystemExit("--n must be an integer >= 2")

    start_time = time.time()
    primes = chudnovsky_like_sieve(N)
    runtime = time.time() - start_time
    print(f"Found {len(primes)} primes up to {N} in {runtime:.2f}s")
    print("Last 5: ", [int(x) for x in primes[-5:]])
    if args.output:
        dirn = os.path.dirname(args.output)
        if dirn:
            os.makedirs(dirn, exist_ok=True)
        with open(args.output, "w") as f:
            f.write("\n".join(str(p) for p in primes))
        print(f"Wrote {len(primes)} primes to {args.output}")
    true_pi = primepi(N)
    print(f"True pi({N}):", true_pi)
    print("Accuracy:", len(primes) == true_pi)
    end_time = time.time()
    print("Runtime: ", end_time - start_time)
    precision, recall, missed, gaps = validate_primes(primes, 2, N)
    print(f"Precision: {precision:.4f}, Recall: {recall:.4f}")
    print(f"Missed primes (first 10): {missed}")
    print(f"Gaps in missed: {gaps}")
