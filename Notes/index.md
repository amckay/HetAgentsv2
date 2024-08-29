

Computational Notes on Heterogeneous-Agent Macroeconomics
===========================================================


Table of contents

 * [Warm-up with RBC Model](RBC.html)
 * [Endogenous Grid Method](EGM.html)
 * [Stationary Equilibrium](Aiyagari.html)
 * [Perfect-Foresight Transitions](MIT.html)
 * [Sequence-Space Methods](SSJac.html)
 * [Heterogeneous Agent New Keynesian Model](HANK.html)
 

Overview
=================================================

These notes provide a crash course on solving heterogeneous-agent macro models.   As a warm up, we will solve a simple representative-agent [RBC model](RBC.html). We then turn to a partial equilibrium consumption-savings problem to introduce the [endogenous grid method](EGM.html) for solving such problems. We then solve for a stationary equilibrium as in the [Aiyagari (1994) model](Aiyagari.html) for which we will  discuss non-stochastic simulation techniques.
It is relatively straightforward to add [perfect-foresight transitions](MIT.html) to a simple model like Aiygari (1994). However, in more complicated models we will find it useful to use richer [sequence-space methods](SSJac.html).  Finally, we use these ideas to solve a [heterogeneous agent New Keynesian Model](HANK.html).

These notes will present the ideas without getting too deep into the code that implements them. The methods are implemented by a set of of [Julia programs](https://github.com/amckay/HetAgentsv2), which are described throughout the notes.

These notes were heavily updated in the Summer of 2024. [The previous version is available here](https://alisdairmckay.com/Notes/HetAgentsV1/).


