# Weighted CPT

Implementation of weighted CPT, an efficient version of conditional permutation test, in R. 

## Summary

This is a follow-up project on the Conditional Permutation Test (CPT), proposed in [https://arxiv.org/abs/1807.05405](https://arxiv.org/abs/1807.05405). 

CPT is a general method for testing the conditional independence of two variables $X$ and $Y$ given a potentially high-dimensional random vector $Z$. Entries of $X$ are permuted non-uniformly to account for the dependence of $X$ on $Z$. However, implementation of CPT is challenging and usually relies on a Markov Chain Monte Carlo sampler. In this work, we propose an Importance-Sampling based version of CPT, where we permute $X$ within groups and define a weighted p-value. Theoretical analysis shows that our proposed test controls the Type I error rate and is robust to model misspecification. We show in simulations that the new test is computationally more efficient than the original CPT while maintaining almost comparable power. 
