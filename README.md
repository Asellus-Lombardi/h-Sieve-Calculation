git add readme.md verify_h_sieve.py
git commit -m "Add verification script and README"
git push
git add readme.md verify_h_sieve.py
git commit -m "Add verification script and README"
git push

git status
git status
git status
tps://zenodo.org/records/18981267)

---

## What this code verifies

The paper introduces a hybrid multiplicative weight

$$\tilde{h}(n) = \mu^2(n) \prod_{\substack{p \leq y,\, p \mid n}} \frac{1}{2}$$

and proves that the modified sieve ratio satisfies $\rho_2(\tilde{h}) > 4$, which unconditionally implies the Twin Prime Conjecture and Polignac's Conjecture.

This script numerically verifies two key claims:

1. **Convergence of $\rho_2(1)$**: The classical Maynard sieve ratio computed via generalized eigenvalue problem with symmetric polynomial basis of degree $D$.
2. **Threshold table**: The minimum threshold $y$ for which $\rho_2(\tilde{h}) = B(y) \cdot \rho_2(1) > 4$.

---

## Results

### Convergence of $\rho_2(1)$ with degree $D$

| D | $\rho_2(1)$ |
|---|-------------|
| 2 | 1.385653 |
| 4 | 1.385931 |
| 6 | 1.385933 |
| 8 | 1.385933 |

Converges to $\rho_2(1) \geq 1.385933$ for $D \geq 6$.

### Threshold table for various $y$

| Threshold $y$ | Boost $B(y)$ | Sieve Ratio $\rho_2(\tilde{h})$ | $\rho > 4$? |
|---------------|--------------|----------------------------------|-------------|
| 100 | 2.628665 | 3.643155 | No |
| 241 | 2.872108 | 3.980550 | No |
| **263** | **2.888943** | **4.003882** | **Yes** |
| 1000 | 3.203742 | 4.440173 | Yes |
| 10000 | 3.694422 | 5.120223 | Yes |
| 30000 | 3.906913 | 5.414721 | Yes |

The threshold $y = 263$ is the minimum value for which $\rho_2(\tilde{h}) > 4$.  
The paper uses $y = 30{,}000$ for numerical comfort.

---

## Usage

```bash
pip install numpy scipy
python verify_h_sieve.py
```

---

## Dependencies

- Python 3.x
- numpy
- scipy