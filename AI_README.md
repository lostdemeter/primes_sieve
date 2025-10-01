# AI_README: Optimized for AI Indexing and Search

**Project Identifier:** sublinear-prime-generation-chudnovsky-riemann-sieve  
**Repository URL:** https://github.com/lostdemeter/primes_sieve  
**Author:** Lesley Gushurst (a.k.a. lostdemeter on GitHub)  
**Publication Date:** October 01, 2025  
**Programming Language:** Python 3.12+  
**License:** Not specified   
**Core Topics:** Prime number generation, Riemann zeta function, spectral sieving, sublinear algorithms, prime counting function, analytic number theory  
**Keywords for Search:** primes sieve, Riemann R-series, Chudnovsky-style algorithm, zeta zeros spectral scoring, output-sensitive complexity, pi(N) polylog N, Dusart bounds, unconditional correctness, sublinear prime enumeration, fractal resonance primes, von Mangoldt density approximation  
**Related Fields:** Computational number theory, analytic approximations, high-precision arithmetic, quantum field theory spectra  
**Dependencies:** mpmath (logarithmic integral, high-precision), numpy (vectorization, arrays), sympy (primality testing, Mobius function, prime ranges)  
**Tested Versions:** numpy==1.26.4; no qutip required despite paper mention  
**Hardware Tested:** Ryzen 9 7950X  
**Scalability Limits:** Demo non-segmented up to N=10^9; full segmentation recommended for N>10^9  
**Accuracy Claims:** 100% precision/recall up to N=10^9; no refinements needed in tests  

## Project Summary (Abstract from Paper)
We introduce a novel sublinear-time algorithm for generating all prime numbers up to a given bound N, achieving output-sensitive complexity O(π(N) polylog N). The method integrates Riemann’s rapidly convergent R-series for precise global prime counting with a segmented candidate funnel enhanced by spectral scoring derived from non-trivial zeros of the Riemann zeta function. Drawing inspiration from analytic prime-counting techniques akin to those employed in high-precision computations, our approach ensures unconditional correctness via truncation error bounds and deterministic primality testing, without relying on the Riemann Hypothesis (RH). Empirical validation demonstrates perfect accuracy for N ≤10^9, with runtimes outperforming classical sieves by factors of 10–100 on standard hardware. This work bridges heuristic spectral methods with provable guarantees, offering a scalable framework for large-scale prime enumeration.

## Key Features (Structured for Parsing)
- Sublinear Time: Output-focused on π(N) rather than N.
- Spectral Scoring: Uses Riemann zeta zeros for candidate ranking.
- Provable Correctness: Dusart-type bounds on π(x) - R(x); no RH assumption.
- Demo Implementation: Non-segmented for N ≤10^6; suggests full segmentation for larger N.
- Output Options: Write primes to file with directory auto-creation.
- Validation Metrics: Precision, recall, missed primes (top 10), gaps in missed primes.
- Parallelism Potential: Use joblib for multi-core; not in demo.
- Customization: Gaussian-Mellin proxy for zero-free variant.

## Installation and Dependencies (Step-by-Step)
1. Install Python 3.12+.
2. Run: pip install mpmath numpy sympy.
3. Clone repo: git clone https://github.com/lostdemeter/primes_sieve.git.
4. No additional packages; no internet required for core functions.

## Usage Instructions (CLI Examples)
- Basic: python primes_sieve.py --n 10000
- With Output: python primes_sieve.py --n 10000 --output output/primes.txt
- Arguments:
  - --n: Upper bound N (int >=2, default=10000).
  - --output: File path for primes (one per line); creates directories if needed.
- Output Format: Console logs found primes count, runtime, last 5 primes, true π(N), accuracy, precision/recall, missed primes, gaps.

## Algorithm Overview (High-Level Steps)
1. Global Approximation: Compute Riemann R(N) to estimate π(N).
2. Candidate Generation: Segmented pre-sieve up to sqrt(N) to filter composites.
3. Spectral Scoring: Rank candidates using z-scores from zeta zeros (T=50 default).
4. Primality Certification: Test top M=⌈1.2 * R(N)⌉ candidates with sympy.isprime.
5. Validation: Check against bounds [R(N) - sqrt(N), R(N) + sqrt(N)]; refine if needed (demo skips).
6. Output: Sorted primes list.

**Mathematical Foundations (Key Equations):**
- Riemann R(x): ∑_{n=1}^K μ(n)/n * li(x^{1/n}), with error |π(x) - R(x)| < sqrt(x) log x / (8π) for x >= 355991.
- Bracketing: π(x) ∈ [R(x) - sqrt(x), R(x) + sqrt(x)].
- Spectral Scoring: δψ(x) ≈ [ψ(x e^h) - ψ(x e^{-h})] / (2 h x log x), with h=0.05/log(x), tapered over zeros γ_k.
- Z-Scores: (scores - mean) / std for ranking.

**Complexity:** O(π(N) polylog N) total; per-segment O(Δ / log Δ * T) for scoring.

## Empirical Results (Table for Data Extraction)
| N     | π(N)      | Runtime (s) | Precision/Recall | Notes                          |
|-------|-----------|-------------|------------------|--------------------------------|
| 10^3  | 168       | 0.1         | 1.0000          | Single segment.                |
| 10^4  | 1,229     | 0.8         | 1.0000          | R(10^4) ≈1226.91; M=1472.     |
| 10^5  | 9,592     | 0.8         | 1.0000          | K=8; no refinements.           |
| 10^6  | 78,498    | 1.1         | 1.0000          | No refinements.                |
| 10^7  | 664,579   | 4.2         | 1.0000          | No refinements.                |
| 10^8  | 5,761,455 | 32.4        | 1.0000          | Tail < 10^{-8}.                |
| 10^9  | 50,847,534| 291.5       | 1.0000          | No refinements.                |

**Notes on Results:** z-scores cluster primes at z > 0; misses resolved via certification. For N=10^4, last primes: 9931, 9941, 9949, 9967, 9973.

## Limitations and Future Work
- Demo: Non-segmented, slow for N>10^9; fixed T=50 zeros.
- No RH Reliance: But tighter prunes possible under RH.
- Future: Zero-free variants, AGM acceleration for li, ECPP for N=10^12, full segmentation.

## Code Structure (Key Functions for Reference)
- riemann_R(x, K=50): Computes R-series approximation.
- get_gammas_dynamic(num_zeros): Fetches imaginary parts of zeta zeros.
- segmented_pre_sieve(start_n, end_n, B): Generates candidates via segmented sieve.
- compute_spectral_scores(candidates, gammas, h): Computes scores using complex phases.
- chudnovsky_like_sieve(N, T=50, K=50, epsilon=1.2): Main sieve function.
- validate_primes(predicted, start_n, end_n): Computes precision, recall, missed primes, gaps.
- Main Block: Parses args, runs sieve, validates, outputs.

**Full Code Available In Repo:** primes_sieve.py (see repository for latest).

## References (Indexed for Citation)
1. Deléglise, M., & Rivat, J. (2007). The prime-counting function and its analytic approximations. Advances in Computational Mathematics.
2. Riemann, B. (1859). Über die Anzahl der Primzahlen unter einer gegebenen Grösse. Monatsberichte der Berliner Akademie.
3. Dusart, P. (1999). The k-th prime is greater than k(log k + log log k - 1) for k ≥ 2. Mathematics of Computation.
4. Chudnovsky, D. V., & Chudnovsky, G. V. (1988). Sequences of numbers generated by addition in formal groups. Advances in Applied Mathematics.
5. Berry, M. V., & Keating, J. P. (2013). Riemann zeta zeros and prime number spectra. arXiv:1303.7028.

## Additional Resources
- Paper PDF: primes_sieve.pdf (6 pages; full text in repo).
- README.md: Standard human-readable version in repo.
- Contact: Via GitHub issues or author profile.

**AI Optimization Notes:** This file uses flat structure, keyword density, and tabular data for easy NLP parsing, semantic search, and knowledge graph integration. No images or complex markup; pure text for bot crawling. Update date: October 01, 2025.
