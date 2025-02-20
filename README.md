# Replication code of "Computationally Efficient Methods for Solving Discrete-time Dynamic models with Continuous Actions" (Fukasawa, 2025)

This repository contains replication code of the paper titled "Computationally Efficient Methods for Solving Discrete-time Dynamic models with Continuous Actions" (previously entitled "Simple method for efficiently solving dynamic models with continuous actions using policy gradient") by Takeshi Fukasawa, investigating  computationally efficient algorithms for solving single-agent/multi-agent dynamic models with continuous actions. The code is written in MATLAB.

Note that the code relies on my other repository takeshi-fukasawa/spectral. If you run the replication code, please also download the spectral repository, and modify the path specified in the main code.

## growth_model_elastic_labor_supply (Section 4.1; Single-agent dynamic model)
The code is based on the replication code of Maliar and Maliar (2013), and it is designed to solve the single-agent neoclassical growth model with elastic labor supply. We can obtain the results of the experiments shown in Section 4.1 of the paper by running gen_results.m. 

## investment_competition_continuous_states (Section 4.2; Multi-agent dynamic game)
The code is designed to solve a dynamic investment competition model with continuous states. We can obtain the results of the experiments shown in Section 4.2 of the paper by running Main_code.m. Note that the code utilizes the Smolyak method. Functions related to the Smolyak method are based on Judd et al. (2014).

## References
* Arellano, C., Maliar, L., Maliar, S., & Tsyrennikov, V. (2016). Envelope condition method with an application to default risk models. Journal of Economic Dynamics and Control, 69, 436-459.
* Coleman, C., Lyon, S., Maliar, L. et al. (2021). Matlab, Python, Julia: What to Choose in Economics?. Comput Economics, 58, 1263–1288.
* Fukasawa, T. (2025). Computationally Efficient Methods for Solving Discrete-time Dynamic models with Continuous Actions, arXiv:2407.04227.
* Judd, K. L., Maliar, L., Maliar, S., & Valero, R. (2014). Smolyak method for solving dynamic economic models: Lagrange interpolation, anisotropic grid and adaptive domain. Journal of Economic Dynamics and Control, 44, 92-123.  
* Pakes, A. & McGuire, P. (1994). Computing Markov-Perfect Nash Equilibria: Numerical Implications of a Dynamic Differentiated Product Model. RAND Journal of Economics, 25, 555-589.
* Maliar, L., & Maliar, S. (2013). Envelope condition method versus endogenous grid method for solving dynamic programming problems. Economics Letters, 120(2), 262-266.
