# SEECS: (ro)vibrational Schrödinger Equation solver with Exterior Complex Scaling

As above, so implemented with a B-spline basis.

## What ?

Numerically solves the 1D Schrödinger equation for a diatomic molecule,

$$\left[ -\frac{1}{2\mu} \frac{d^2}{dR^2} + V(R) + \frac{j(j+1)}{2\mu R^2} \right ] \psi(R) = E \psi(R)$$

given a tabulated internuclear potential $V(R)$, a reduced mass $\mu$, and a rotational quantum number $j$.
Returns the energies $E$, wavefunctions $\psi(R)$ evaluated on a grid, and the rotational constant $B_v$ for each state.

**Exterior Complex Scaling (ECS)** makes the solutions above the dissociation threshold meaningfully different than (half-)box states.
ECS is similar in goal and effect to other methods like Complex Absorbing Potentials (CAPs, AKA optical potentials).
CAPs get their name from their method.
Adding an imaginary potential to the tail end of the real-valued potential to produce a complex-valued potential beyond some cutoff distance $R_0$.
ECS differs in that the potential remains real and it is the R-grid that changes past a cutoff distance $R_0$: the R-grid past $R_0$ gets rotated into the complex plane by an angle $\theta$:

$$R \longrightarrow z(R) = R_0 + (R - R_0)\,e^{i\theta}, \qquad R > R_0,$$

Waves at $R>R_0$ decay, allowing the continuum to be represented in a finite basis of bound-state wavefunctions (hopefully) without reflection artifacts.

## How ?

- **Basis:** clamped B-splines of order $k$ on $[R_\mathrm{min}, R_\mathrm{max}]$ with `nelem` knot intervals. Boundary conditions are imposed by dropping basis functions at each edge (`ndropl`, `ndropr`).
- **Matrix elements:** Hamiltonian and overlap matrices are assembled element-by-element with Gauss–Legendre quadrature.
- **Eigenproblem:** the generalized problem $Hc = ESc$ is solved with a real symmetric solver (no ECS) or a complex-symmetric solver (ECS). True bound states in the ECS case will be mostly real, while resonances will have a significant imaginary component associated with the resonance width Γ.
- **Rotational constant:** $B_v = \langle v|1/(2\mu R^2)|v\rangle$, computed as a quadratic form on the eigenvectors.

## Units

If called as a library, everything is in atomic units.
If the `program` is called, then the namelist parameters can be used to selected unit conversions.
The scaling angle `theta` is given in degrees; `R0`, `rmin`, `rmax` are in the same units as the input $R$ grid (`runits_in`).

## Building

SEECS builds with [fpm](https://fpm.fortran-lang.org/):

```sh
fpm build
fpm run < input.namelist
```

LAPACK is required.
By default SEECS declares the LAPACK interfaces itself and links a system LAPACK.

## Input

TODO

## Output

TODO

## Status / limitations

TODO

## License

MIT
