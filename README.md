# h-Sieve: Numerical Verification

This repository contains the numerical verification code for the paper:

**"The h-Sieve: A Hybrid Weight Approach to Prime Tuples"**
X.A., March 12, 2026
DOI: [10.5281/zenodo.18981267](https://zenodo.org/records/18981267)

---

## What this code verifies

This script numerically verifies two key claims in the paper:

1. **Convergence of rho_2(1)**: The classical Maynard sieve ratio computed via generalized eigenvalue problem with symmetric polynomial basis of degree D.
2. **Threshold table**: The minimum threshold y for which rho_2(h~) = B(y) * rho_2(1) > 4.

---

## Results

### Convergence of rho_2(1) with degree D

| D | rho_2(1)  |
|---|-----------|
| 2 | 1.385653  |
| 4 | 1.385931  |
| 6 | 1.385933  |
| 8 | 1.385933  |

### Threshold table for various y

| Threshold y | Boost B(y) | Sieve Ratio rho_2(h~) | rho > 4? |
|-------------|------------|-----------------------|----------|
| 100         | 2.628665   | 3.643155              | No       |
| 241         | 2.872108   | 3.980550              | No       |
| **263**     | **2.888943**| **4.003882**         | **Yes**  |
| 1000        | 3.203742   | 4.440173              | Yes      |
| 10000       | 3.694422   | 5.120223              | Yes      |
| 30000       | 3.906913   | 5.414721              | Yes      |

y = 263 is the minimum threshold for which rho_2(h~) > 4.
The paper uses y = 30,000 for numerical comfort.

---

## Usage
```
pip install numpy scipy
python verify_h_sieve.py
```

## Dependencies

- Python 3.x
- numpy
- scipy
