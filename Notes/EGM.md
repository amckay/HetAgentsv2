

Endogenous Grid Method
==================================


Consider a consumer facing fluctuating income $e$. We will assume that $e$ takes $N_e$ values and follows a Markov chain with transition probabilities $\Pi(e'|e)$ where a prime
denotes next period's value of a variable.  The consumer can save in a risk-free bond at return $R$. The budget constraint is therefore
$$  a' + c = R a + e,$$

where $a$ is assets before interest.  Let's suppose that there is a borrowing constraint $a \geq 0.$ The consumer wishes to maximize a utility function of
$$\mathbb E_0 \sum_{t=0}^\infty \beta^t u(c_t)$$
where $u(c) = c^{1-\gamma}/(1-\gamma)$.

The states of the consumer's problem are $(a,e)$ and so a solution to this problem can be represented as a function $g(a,e)$ that gives the choice of $a'$.

The Bellman equation is
$$V(a,e) = \max_{a'\geq 0} \left\{ u(Ra  + e - a') + \beta \mathbb E V(a',e') \right\},$$
where the expectation is taken over possible realizations of $e'$ conditional on $e$.  Notice that we are including uncertainty over idiosyncratic shocks even though we will (later) assume certainty equivalence (i.e. perfect foresight) over aggregate variables. For now we assume aggregate variables are constant.


The first-order and envelope conditions are
$$\begin{aligned}
u'(c) &\geq \beta  \mathbb E V_a(a',e') \\
V_a(a,e) &= R u'(c)
\end{aligned}$$


The heart of the endogenous grid method is as follows.  Suppose this period we have states $(a,e)$ where we know $e$, but as of 
yet we don't know $a$. Suppose we save an amount $a'$ so we have $g(a,e) = a'$ and we will have to determine $a$. Given  the function $V_a$ that will prevail in the future, we use the first-order condition to solve for $c$ and then the budget constraint to solve for $a$
$$a = (a' +c - e)/R$$
(recall we know $a'$). We then use the envelope condition to calculate $V_a(a,e)$ (i.e. the marginal value of assets today for a single $(a,e)$ pair. Doing this for many $a'$ values, we can map out the whole function $V_a(a,e)$. 

This procedure started with a guess of the $V_a$ function that prevails next period, and then produced the $V_a$ function that prevails this period. We use the new function as our guess and apply the same steps again, repeating until the guess reproduces itself.



## Implementation concepts


We now discuss some of the details of how we implement these ideas.

### Approximating the decision rules


A function cannot itself be represented in a computer so the first decision we have to make is how to approximate the savings policy rule.
Let $A$ be a grid of $N_A$ points on values of $a$ with the first gridpoint at the borrowing constraint, zero in this case. 
Now let $G$ be a $N_A \times N_e$ matrix where each column is a vector of length $N_A$ that represents the policy rule as a function of $a$ for one of the values of $e$.
The interpretation now is that if you have states $(A_{i},e_j)$ then the policy rule calls for saving $G_{ij}$.  For values of $a$ between two values of $A_{i}$ and $A_{i+1}$ we will use linear interpolation to fill in the function.

###  Implementing the endogenous grid method


The endogenous grid method imposes a grid on $a'$. We will use the same one that we use to approximate the savings policy. When we apply steps above, we find values of $a$ that map into the savings choices on the grid $A$. This is almost what we want,
but not quite. We want to know the values $G$ that correspond to the savings choices at current states $A$, but what we have calculated is the 
states $a$ that map into the savings choices $A$. 

For a given current $e_j$, the steps above yield a set of points $\tilde A$ such that $g(\tilde A_i,e_j) = A_i$. What we would like is a set of savings levels $G_{\cdot j}$ such that $g(A_i,e_j) = G_{ij}$.  We can use linear interpolation to find  $G_{\cdot j}$.

## Code 

An example application of the endogenous grid method can be found in `EGM.jl`

## Next step

We next discuss methods to analyze a population of households facing this type of consumption-savings problem, which brings us to the [Aiyagari (1994) model](Aiyagari.html).