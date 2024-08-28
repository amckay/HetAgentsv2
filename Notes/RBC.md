Warming up with the RBC Model
=============================

Consider an RBC model in which preferences are given by
$$ \sum_{t=0}^\infty \beta^t \frac{C_t^{1-\gamma}}{1-\gamma},$$
production follows
$$Y_t = Z_t K_{t-1}^\alpha \bar L^{1-\alpha},$$
capital evolves according to
$$K_t = (1-\delta)K_{t-1}+Y_t - C_t,$$
and productivity evolves according to
$$\log Z_t = \rho \log Z_{t-1} + \varepsilon_{t},$$
where $\varepsilon$ is an exogenous, mean-zero innovation to TFP.  I have adopted a convention of dating the capital stock selected in
period $t$ and used in production in $t+1$ as $K_t$ so that a variable dated $t$ is measurable with respect
to date-$t$ information.

The model can be summarized by the following expectational difference equations
$$\begin{aligned}C_t^{-\gamma} &= \beta  R_{t+1} C_{t+1}^{-\gamma}  \\

   R_t &= \alpha Z_t \left( K_{t-1} / \bar L \right)^{\alpha-1} + 1 - \delta \\
 
   K_t &= (1-\delta)K_{t-1}+Y_t - C_t \\

   Y_t &= Z_t K_{t-1}^\alpha \bar L^{1-\alpha} \\

   \log Z_t &= \rho \log Z_{t-1} + \varepsilon_{t}. \end{aligned}
   $$

We will fix $\bar L = 1$.  The steady state is then
$$
\begin{aligned}
   \bar Z &= 1 \\
   \bar R &= 1/ \beta \\
   \bar K &= \left( \frac{\bar R-1+\delta}{\alpha} \right)^{1/(\alpha - 1)} \\
   \bar Y &= \bar K^\alpha \\
   \bar C &= \bar Y - \delta \bar K.
\end{aligned}
$$   
where we use the convention that bars denote steady state values.

We are going to solve for a perfect foresight transition path: the economy is at steady state and expected to remain there when at $t=0$ a realization of $\varepsilon_0 \neq 0$ occurs. The economy assumes no further $\varepsilon$'s will occur. This experiment sounds strange---why would the economy not expect these shocks even when they just saw one? It turns out that this experiment gives us the same impulse response to the $\varepsilon$ shock as a first-order perturbation of a stochastic environment.  The first-order perturbation solution features certainty equivalence, which means the agents behave as if future random variables are replaced by their mean values. In the case of $\varepsilon$, that means $\mathbb E_0[\varepsilon_t] = 0$ for all $t > 0$.

Let $X_t=\left\{C_t,R_t,K_t,Y_t,Z_t\right\}$ be the endogenous variables at date $t$. The unknown is a sequence $X \equiv \left\{X_t\right\}_{t=0}^T$.  In our computations, we will assume that after some large $T$ the economy has returned to steady state. We also assume the economy was in steady state before the shock occurs. Therefore we use the conventions $x_{T+1} = x_{-1} = \bar x$ for any variable $x$. We will use $E_t$ for endogenous variables, which in this case is $E_t = \varepsilon_t$. $E$ without a subscript is the sequence $E \equiv \{E_t\}_{t=0}^T.$  Finally, let $\bar X$ and $\bar E$ be the sequences in which all variables remain at their steady state values.


At a given date, the equations of the model can be written as 
$$f(X_{t-1},X_t,X_{t+1},E_t) = 0.$$
We can stack these equations for all dates to write 
$$f(X,E) = 0.$$
Our goal is to simply solve this equation for $X$ given a value of $E$.
We will consider two stategies. First, we can ``linearize'' the model. Let $f_X(\bar X, \bar E)$ and $f_E(\bar X, \bar E)$ be the Jacobians of $f$ with respect to $X$ and $E$, respectively, both evaluated at steady state. We can then solve for a first-order approximation as
$$\begin{aligned}
&f_X(\bar X,\bar E)(X - \bar X) + f_E(\bar X,\bar E)(E - \bar E) = 0 \\
&X = \bar X -  \left[ f_X(\bar X,\bar E) \right]^{-1}f_E(\bar X,\bar E)(E - \bar E).
\end{aligned}
$$
Alternatively, we can seek a non-linear solution using Newton's method. Suppose we have a candidate solution in iteration $j$ denoted $X^{(j)}$, we then form a new solution $X^{(j+1)}$ as
$$
X^{(j+1)} =  X^{(j)} -  \left[ f_X( X^{(j)}, E) \right]^{-1}f( X^{(j)}, E).
$$
We iterate on this equation until $f( X^{(j)}, E)$ is approximately zero.

## Code

The RBC example is solved by the file `RBC.jl`, which makes use of `ModelUtils.jl`. We will be using `ModelUtils.jl` again later in this course so reading the [documentation](https://github.com/amckay/ModelUtils.jl) would be a good investment. 

## Next step

We will now turn to heterogeneous agent models starting with the [endogenous grid method](EGM.md)