# Replication code of "Simple method for efficiently solving dynamic models with continuous actions using policy gradient" (Fukasawa, 2024)

This repository contains replication codes of the paper titled "Simple method for efficiently solving dynamic models with continuous actions using policy gradient" by Takeshi Fukasawa, proposing VF-PGI-Spectral (Value function-Policy Gradient Iteration-Spectral) algorithm for computationally efficiently solving dynamic models with continuous actions. These codes are written in MATLAB.

Note that the codes rely on my other repository takeshi-fukasawa/spectral. If you run the replication codes, please also download the spectral repository, and modify the path specified in the main code.

## growth_model_elastic_labor_supply (Section 5.1, Appendix C.1)
The code is based on the replication code of Maliar and Maliar (2013), and it is designed to solve the single-agent neoclassical growth model with elastic labor supply. We can obtain the results of the experiments shown in Section 5.1 and Appendix C.1 of the paper by running gen_results.m and gen_results_VF_PGI_test.m.

## investment_competition_continuous_states (Section 5.2)
The code is designed to solve a dynamic investment competition model with continuous states. We can obtain the results of the experiments shown in Section 5.2 of the paper by running Main_code.m. Note that the code utilizes the Smolyak method. Functions related to the Smolyak method are based on Judd et al. (2014).

## Pakes_McGuire_model (Section 5.3)
The code is based on the replication code of Pakes and McGuire (1994), downloaded from https://scholar.harvard.edu/pakes/pages/pakes-maguire-algorithm-0. It is designed to solve dynamic investment competition model with discrete states considered by Pakes and McGuire (1994). We can obtain the results of the experiments shown in Section 5.3 of the paper by running runit.m.

## growth_model_inelastic_labor_supply (Appendix C.2)
The code is based on the replication code of Coleman (2021) and Arellano et al. (2016), and it is designed to solve the single-agent neoclassical growth model with inelastic labor supply. We can obtain the results of the experiments shown in Appendix C.2 of the paper by running gen_results.m.


## References
* Arellano, C., Maliar, L., Maliar, S., & Tsyrennikov, V. (2016). Envelope condition method with an application to default risk models. Journal of Economic Dynamics and Control, 69, 436-459.
* Coleman, C., Lyon, S., Maliar, L. et al. (2021). Matlab, Python, Julia: What to Choose in Economics?. Comput Economics, 58, 1263–1288.
* Fukasawa, T. (2024). Simple method for efficiently solving dynamic models with continuous actions using policy gradient, mimeo.  
* Judd, K. L., Maliar, L., Maliar, S., & Valero, R. (2014). Smolyak method for solving dynamic economic models: Lagrange interpolation, anisotropic grid and adaptive domain. Journal of Economic Dynamics and Control, 44, 92-123.  
* Pakes, A. & McGuire, P. (1994). Computing Markov-Perfect Nash Equilibria: Numerical Implications of a Dynamic Differentiated Product Model. RAND Journal of Economics, 25, 555-589.
* Maliar, L., & Maliar, S. (2013). Envelope condition method versus endogenous grid method for solving dynamic programming problems. Economics Letters, 120(2), 262-266.
