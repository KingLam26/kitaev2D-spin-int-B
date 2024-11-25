# Integer-spin-S Kitaev 2D Honeycomb Model

We consider the arbitary integer-spin-S Kitaev 2D Honeycomb Model defined by 

$$H_0=-J_z \sum_{\langle i, j\rangle_z} S_i^z S_j^z, \quad V=-J_x \sum_{\langle i, j\rangle_x} S_i^x S_j^x-J_y \sum_{\langle i, j\rangle_y} S_i^y S_j^y,$$

in the anisotropic limit. The lowest-order ($4S^{\textrm{th}}$) non-trivial effective Hamitonian is given by 

$$H_{\text {eff }} = J_{\text {eff }} \sum_i\sigma_i^x,$$

where the sum is taken over all $z$ bonds in the 2D honeycomb lattice, and $\sigma_i^x$ is a pseudo-half-spin with two degrees-of-freedom.

A key problem we are trying to solve is whether $J_{\text {eff }}$ is protected against fine-tuning across the spectrum. We hope to achieve this without an explicit computation of $J_{\text {eff }}$, which is computationally intractable for large arbitary spins. Using purely symmetry arguments, we showed that $J_{\text {eff }}$ can be expressed as a symbolic palindromic polynomial, whose coefficients represent the amplitudes of the diagrams that represent groups of terms contributing to $J_{\text {eff}}$. Then, the problem of showing that $J_{\text {eff }}$ is always non-vanishing is translated into a purely mathematical one of proving that the palindromic polynomial has no real roots.

This repository comprises of a number of separate scripts which were developed to manipulate and analyze these symbolic palindromic polynomials (PPs). A key feature we observed and discussed in our paper is the underlying structure relating the different diagram amplitudes. We exploit this in order to factorise the PPs and write them as a sum of squares of PPs.

As an example, in the case of even integer spins, the full PP $g(z)$ can be written as

$$g(z):=\tilde{A}_1(z,(a_1)_m)+z^2\tilde{A}_2(z,(\epsilon_1)_{mn}),$$

where 

$$\begin{aligned}&\tilde{A}_1(z,(a_1)_m)=A_1^2(z,(a_1)_m),\\&A_1(z,(a_1)_m):=z^S-(a_1)_1z^{S-1}+(a_1)_2z^{S-2}\cdots+(a_1)_{S/2}z^{S/2}\cdots+(a_1)_2z^2-(a_1)_1z^1+1.\end{aligned}$$

and $\tilde{A}_2(z,(\epsilon_1)_{mn})$ can also be factorised in a similar fashion. For odd integer spins, a similar factorisation process follows, which we discuss in greater detail in our paper.

**We stress that this repository is still under development, and is used mainly as a tool.**

## File Structure
1. ``main_poly.py``: Considers all possible diagrams that sum to $J_{\text {eff}}$ for a certain integer spin model, compute the full PP, and factorises it.
2. ``utils_poly.py``: Support for ``main_poly.py``.
3. ``quick_factor.py``: a quick and simple tool to factorise any polynomial.
4. ``coeff_anal.py``: Takes in explicit numerical diagram coefficients and analyzes them.

## Requirements
1. ``Python 3.x``
2. ``SymPy``