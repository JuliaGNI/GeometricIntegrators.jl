
# Discontinuous Galerkin Variational Integrators

\\[
\delta \sum \limits_{n=0}^{N-1} \Bigg\{
    \sum \limits_{i=1}^{s} b_i \, L \big( q_h (t_n + c_i h), \, \dot{q}_h (t_n + c_i h) \big)
    + \dfrac{1}{2} \sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta \big( \phi (\gamma_i; q_{n} , q_{n}^+) \big) \, \dfrac{d\phi}{d\tau} (\gamma_i; q_{n} , q_{n}^+)
    + \dfrac{1}{2} \sum \limits_{i=1}^{\sigma} \beta_i \, \vartheta \big( \phi (\gamma_i; q_{n+1}^- , q_{n+1}) \big) \, \dfrac{d\phi}{d\tau} (\gamma_i; q_{n+1}^- , q_{n+1})
\Bigg\} = 0
\\]
