# Fast-Fourier-Transform-using-Cooley-Tukey-algorithm
the most common fast Fourier transform (FFT) algorithm

Cooley–Tukey re-expresses the discrete Fourier transform (DFT) of an arbitrary composite size $N = N_1N_2$ in terms of smaller DFTs of sizes $N_1$ and $N_2$, recursively, in order to reduce the computation time to $O(N log N)$ for highly composite N (smooth numbers).

####The radix-2 DIT case

A radix-2 decimation-in-time (DIT) FFT is the simplest and most common form of the Cooley–Tukey algorithm, although highly optimized Cooley–Tukey implementations typically use other forms of the algorithm as described below. Radix-2 DIT divides a DFT of size N into two interleaved DFTs (hence the name "radix-2") of size $N/2$ with each recursive stage.

The discrete Fourier transform (DFT) is defined by the formula:

$$
X_k = \sum_{n=0}^{N-1} x_n e^{-\frac{2\pi i}{N} nk},
$$

where $k$ is an integer ranging from $0 to N-1$.

Radix-2 DIT first computes the DFTs of the even-indexed inputs $(x_{2m}=x_0, x_2, \ldots, x_{N-2})$ and of the odd-indexed inputs $(x_{2m+1}=x_1, x_3, \ldots, x_{N-1})$, and then combines those two results to produce the DFT of the whole sequence. This idea can then be performed recursively to reduce the overall runtime to $O(N log N)$. This simplified form assumes that N is a power of two; since the number of sample points $N$ can usually be chosen freely by the application, this is often not an important restriction.

The Radix-2 DIT algorithm rearranges the DFT of the function $x_n$ into two parts: a sum over the even-numbered indices $n={2m}$ and a sum over the odd-numbered indices $n={2m+1}$:

$$
  \begin{matrix} X_k & =
& \sum \limits_{m=0}^{N/2-1} x_{2m}e^{-\frac{2\pi i}{N} (2m)k}   +   \sum \limits_{m=0}^{N/2-1} x_{2m+1} e^{-\frac{2\pi i}{N} (2m+1)k}
  \end{matrix}
$$

One can factor a common multiplier $e^{-\frac{2\pi i}{N}k}$ out of the second sum, as shown in the equation below. It is then clear that the two sums are the DFT of the even-indexed part $x_{2m}$ and the DFT of odd-indexed part $x_{2m+1}$ of the function $x_n$. Denote the DFT of the Even-indexed inputs $x_{2m}$ by $E_k$ and the DFT of the Odd-indexed inputs $x_{2m + 1}$ by $O_k$ and we obtain:

$$
\begin{matrix} X_k= \underbrace{\sum \limits_{m=0}^{N/2-1} x_{2m}   e^{-\frac{2\pi i}{N/2} mk}}_{\mathrm{DFT\;of\;even-indexed\;part\;of\;} x_m} {} +  e^{-\frac{2\pi i}{N}k}
 \underbrace{\sum \limits_{m=0}^{N/2-1} x_{2m+1} e^{-\frac{2\pi i}{N/2} mk}}_{\mathrm{DFT\;of\;odd-indexed\;part\;of\;} x_m} =  E_k + e^{-\frac{2\pi i}{N}k} O_k.
\end{matrix}
$$

Thanks to the periodicity of the DFT, we know that

$$
E_{{k} + \frac{N}{2}} = E_k
$$

and

$$
O_{{k} + \frac{N}{2}} = O_k. 
$$

Therefore, we can rewrite the above equation as

$$
\begin{matrix} X_k & = & \left\{
\begin{matrix}
E_k + e^{-\frac{2\pi i}{N}k} O_k & \mbox{for } 0 \leq k < N/2 \\ \\
E_{k-N/2} + e^{-\frac{2\pi i}{N}k} O_{k-N/2} & \mbox{for } N/2 \leq k < N . \\
\end{matrix}
\right. \end{matrix}
$$

We also know that the twiddle factor $e^{-2\pi i k/ N}$ obeys the following relation:

$$
\begin{matrix} e^{\frac{-2\pi i}{N} (k + N/2)} & = & e^{\frac{-2\pi i k}{N} - {\pi i}} \\
& = & e^{-\pi i} e^{\frac{-2\pi i k}{N}} \\
& = & -e^{\frac{-2\pi i k}{N}}
\end{matrix}
$$

This allows us to cut the number of "twiddle factor" calculations in half also. For  $0 \leq k < \frac{N}{2}$, we have

$$
\begin{matrix}
X_k & =
& E_k + e^{-\frac{2\pi i}{N}k} O_k \\
X_{k+\frac{N}{2}} & =
& E_k - e^{-\frac{2\pi i}{N}k} O_k
\end{matrix}
$$

This result, expressing the DFT of length $N$ recursively in terms of two DFTs of size $N/2$, is the core of the radix-2 DIT fast Fourier transform. The algorithm gains its speed by re-using the results of intermediate computations to compute multiple DFT outputs. Note that final outputs are obtained by a +/− combination of $E_k$ and $O_k \exp(-2\pi i k/N)$, which is simply a size-2 DFT (sometimes called a butterfly in this context); when this is generalized to larger radices below, the size-2 DFT is replaced by a larger DFT (which itself can be evaluated with an FFT).

This process is an example of the general technique of divide and conquer algorithms; in many traditional implementations, however, the explicit recursion is avoided, and instead one traverses the computational tree in breadth-first fashion.

