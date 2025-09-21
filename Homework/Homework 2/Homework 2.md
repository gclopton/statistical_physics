

> [!Question]
> The mean value $\mu$ (or the average) of $N$ numbers is:
> 
> $$
> \mu=\frac{1}{N} \sum_{i=1}^N a_i
> $$
> 
> 
> Write a Python function my_mean that takes an array a and returns its mean value.



```Python
import numpy as np

def my_mean(a):
    """
    Calculates the mean value of an array

    Args:
        a (np.array): time series (or any array-like of numbers)
    Returns:
        float: mean value
    """
    arr = np.asarray(a, dtype=float)
    if arr.size == 0:
        raise ValueError("my_mean: empty array has no mean.")
    return float(arr.sum() / arr.size)

```



> [!Question]
> The sample standard deviation $\sigma$ of $N$ numbers is:
> 
> $$
> \sigma=\sqrt{\frac{1}{N-1} \sum_{i=1}^N\left(O_i-\bar{O}\right)^2}=\sqrt{\frac{1}{N-1}} \sqrt{\left(\sum_{i=1}^N O_i^2\right)-\frac{1}{N}\left(\sum_{i=1}^N O_i\right)^2}
> $$
> 
> 
> Write a Python function my_stddev that takes an array a and returns its standard deviation.



```Python
import numpy as np

def my_std(a):
    """
    Computes the sample standard deviation (1/sqrt(N-1) normalization)
    deducting 1 degree of freedom because the mean is estimated.

    Args:
        a (np.array): time series (or any array-like of numbers)
    Returns:
        float: standard deviation of a
    """
    arr = np.asarray(a, dtype=float)
    n = arr.size
    if n < 2:
        raise ValueError("my_std: need at least two values to compute sample standard deviation.")
    mean = arr.sum() / n
    diff = arr - mean
    ssd = np.dot(diff, diff)  # sum of squared deviations
    return float(np.sqrt(ssd / (n - 1)))

```



> [!Question]
> Given a dataset a [] with $N$ elements, mean $\bar{a}$, and standard deviation $\sigma$, we define the autocorrelation time of a as
> 
> $$
> \kappa=1+2 \sum_{t=1}^{t_{\text {cutoff }}-1} C(t),
> $$
> 
> where (using the python convention of arrays going from 0 to $\mathrm{N}-1$ ):
> 
> $$
> C(t)=\frac{1}{\sigma^2} \frac{1}{N-t} \sum_{i=0}^{N-t-1}\left(a_i-\bar{a}\right)\left(a_{i+t}-\bar{a}\right),
> $$
> 
> and $t_{\text {cutoff }}$ is the smallest value of $t$ for which $C(t) \leq 0$. Do not include $C\left(t_{\text {cutoff }}\right)$ when computing $\kappa$.
> Note: This convention of choosing $t_{\text {cutoff }}$ is not universal! Remember this is not a fundamental result, but a practical method to estimate the autocorrelation time from a finite data set. Our approach of choosing $t_{\text {cutoff }}$ at the first zero crossing is straightforward and easy to understand, but there are other, more sophisticated schemes you can come up with.
> 
> Write a Python function my_actime that takes an array a and returns the autocorrelation time.

```Python
import numpy as np
def my_actime(a): 
  """
  Computes the autocorrelation time of a

  Args:
    a (np.array): time series
  Returns:
    float: autocorrelation time of a
  """
  pass
```


```Python
import numpy as np

def my_actime(a):
    """
    Computes the autocorrelation time of a, using:
      kappa = 1 + 2 * sum_{t=1}^{t_cutoff-1} C(t),
    where C(t) = [1/(sigma^2)] * [1/(N-t)] * sum_{i=0}^{N-t-1} (a_i - mean)*(a_{i+t} - mean),
    and t_cutoff is the smallest t with C(t) <= 0 (excluded from the sum).

    Args:
        a (np.array): time series (array-like)
    Returns:
        float: autocorrelation time of a
    """
    x = np.asarray(a, dtype=float)
    N = x.size
    if N < 2:
        raise ValueError("my_actime: need at least two points.")

    # Center and use sample variance (ddof=1) to match the spec
    x = x - x.mean()
    sigma2 = (x @ x) / (N - 1)
    if sigma2 <= 0:  # constant or degenerate series
        return 1.0

    # Unbiased autocovariance for lags t = 0..N-1:
    # cov[t] = (1/(N-t)) * sum_{i=0}^{N-t-1} x[i]*x[i+t]
    acf_full = np.correlate(x, x, mode='full')[N-1:]     # sums for t >= 0
    denom = (N - np.arange(N)).astype(float)
    cov = acf_full / denom

    # Normalized autocorrelation C(t)
    C = cov / sigma2

    # Find first non-positive C(t) for t >= 1
    nz = np.where(C[1:] <= 0)[0]
    t_cutoff = 1 + int(nz[0]) if nz.size else N

    # Integrated autocorrelation time
    kappa = 1.0 + 2.0 * float(C[1:t_cutoff].sum())
    return kappa

```



> [!Question]
> The standard error $\epsilon$ of the mean of $N$ numbers can be calculated from the sample standard deviation $\sigma$ and auto-correlation $\kappa$ :
> 
> $$
> \epsilon=\frac{\sigma}{\sqrt{N / \kappa}}
> $$
> 
> 
> Write a Python function my_stderr that returns the standard error of an array a.



```Python
import numpy as np
def my_stderr(a):
  """
  Computes the standard error of the mean of a

  Args:
    a (np.array):  time series
  Returns:
    float: standard error of a
  """
  pass
```


```Python
import numpy as np

def my_stderr(a):
    """
    Computes the standard error of the mean of a:
      epsilon = sigma * sqrt(kappa / N),
    where sigma is the sample std (ddof=1) and
    kappa = 1 + 2 * sum_{t=1}^{t_cutoff-1} C(t),
    with t_cutoff the first t where C(t) <= 0 (excluded).

    Args:
        a (np.array): time series
    Returns:
        float: standard error of a
    """
    x = np.asarray(a, dtype=float)
    N = x.size
    if N < 2:
        raise ValueError("my_stderr: need at least two values (uses sample std, ddof=1).")

    # Sample std (ddof=1)
    mean = x.mean()
    diff = x - mean
    ssd = float(diff @ diff)
    sigma2 = ssd / (N - 1)
    if sigma2 <= 0.0:
        return 0.0  # constant/degenerate series: zero stderr

    # Autocorrelation function C(t) with unbiased normalization
    acf_sums = np.correlate(diff, diff, mode='full')[N-1:]     # sums for lags t >= 0
    denom = (N - np.arange(N)).astype(float)                  # (N - t)
    cov = acf_sums / denom
    C = cov / sigma2

    # First non-positive C(t) for t >= 1
    nz = np.where(C[1:] <= 0)[0]
    t_cutoff = 1 + int(nz[0]) if nz.size else N

    # Integrated autocorrelation time kappa
    kappa = 1.0 + 2.0 * float(C[1:t_cutoff].sum())
    if not np.isfinite(kappa) or kappa <= 0:
        kappa = 1.0

    # Standard error
    return float(np.sqrt(sigma2) * np.sqrt(kappa / N))

```




> [!Question] 2.5
> ![[Screenshot 2025-09-21 at 12.33.29 AM.png]]






> [!Question] 2.6
> ![[Screenshot 2025-09-21 at 12.34.10 AM.png]]



> [!Question] 2.7
> ![[Screenshot 2025-09-21 at 12.35.18 AM.png]]




> [!Question] 2.8
> ![[Screenshot 2025-09-21 at 12.36.08 AM.png]]




> [!Question] 2.9
> ![[Screenshot 2025-09-21 at 12.36.54 AM.png]]








> [!Question] 2.10
> ![[Screenshot 2025-09-21 at 12.38.38 AM.png]]




> [!Question] 2.11
> ![[Screenshot 2025-09-21 at 12.39.13 AM.png]]





> [!Question] 2.12
> 
> ![[Screenshot 2025-09-21 at 12.39.49 AM.png]]
> 
> 
> 



> [!Question] 2.13
> ![[Screenshot 2025-09-21 at 12.40.27 AM.png]]


