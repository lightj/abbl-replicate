# Arellano et al. 2021 - Replication Files

Replication codes for: "Heterogeneity of Consumption Responses to Income Shocks in the Presence of Nonlinear Persistence", by Arellano, Blundell, Bonhomme, Light

If you have any questions please e-mail: jdlight@uchicago.edu

# Overview

There are three sets of replication codes:

- Stata:
  - STATA codes to replicate the PSID sample:
    - 1_Master.do: runs everything
    - 2_Houehold_Panel.do: constructs the household panel
    - 3_Individual_Panel.do: constructs the household heads
    - 3b_Individual_Panel_Parent_Link.do: gets parental linkages
    - 4_Sample_Selection.do: does the main sample selection
    - 5_Aggregate_Measures.do: creates aggregates
    - 5b_Aggregate_Measures_Parent.do: creates aggregates for parents
    - 6_Get_Residuals_JOE_vF.do: generates the main samples
    
- Data:
  - data that you should be able to use to run the codes in folders Inc and Cons
  
- Inc:
  - Codes for labor and disposable income, codes to replicate the main figures
    - To generate the main income results you should use SMC_l_vF.m (labor income) or SMC_n_vF.m (disposable income)
    - The code will call the various other scripts
     
- Cons:
  - Quantile Model: codes to estimate baseline model, model with asset rule, model with no heterogeneity, model with heterogeneity and no latent income components
    - PMCMC_simple_cons_vF.m: this estimates the baseline model
    - PMCMC_assets_cons_vF.m: this estimates the baseline model with assets
    - PMCMC_y_cons_vF.m: this estimates the model with heterogeneity and no filtering
    - SMC_simple_cons_vF.m: this estimates the model with filtering and no heterogeneity
    - The codes will call the various other scripts
  - Plot Results: codes to replicate the main figures
    - Figures_baseline_vF.m: replicates the figures for the baseline model
    - Figures_assets_vF.m: replicates the figures for the model with assets
    - Figures_irf_baseline_vF.m: replicates the figures for the baseline model IRFs
    - Figures_smc_vF.m: replicates the figures for the model without heterogeneity
    - Figures_y_vF.m: replicates the figures for the model with heterogeneity and no filtering
    
- Bootstrap:
  - Codes that should let you run the paramteric and non-parametric bootstrap
