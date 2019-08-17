Homework for 10/31: Prepping to build the model
================
Kelsey Elwood
10/31/2018

Background Information
======================

The data I am using for the individual project looks at the first flowering date (FFD) of various alpine species in response to increases in temperature, snowpack, and N. Data was collected for 2007 and 2008, though for now I am only looking at 2008 data. The data is structured with three blocks of 16 plots each. In each block, a snowfence bisects the block to form sub-blocks, with half of the plots on the windward (normal to reduced snowpack) and half of the plots on the leeward (increased snowpack) side of the fence. In each sub-block, 2 plots have warming treatments, 2 plots have increased nitrogen, 2 plots have both warming and increased nitrogen, and 2 plots are left as controls. Figure 1 (below) captures the conceptual layout of the plots.

1. Scales of variables
======================

**Identify the scales of variation in the sampling design. What grouping variables are needed?**

For this project, the scales of variation could be both spatial and temporal, but for now, I plan to only explore the spatial variation. There are 2 grouping variables that are nested: 3 blocks and 2 snow sub-blocks (separated by a snow fence) within each block. As I gain confidence with using multi-level models, I could include the 2 temporal groupings of 2007 and 2008. Initially, I plan to only use data from 2008, so I would not have any temporal groupings.

*I'm still a little confused about the most appropriate way to address snow as a grouping variable versus a predictor variable. The intention of the experiment was to have snow be a predictor. The snow fence that creates the snow sub-block is purely a product of logistical constraints on how to produce plots with increased snow. Could you help me better understand the appropriate way to include snow in the model?*

**Identify the scale of the response variable - this will be the smallest scale of sampling unit that you can estimate.**

The response variable in this experiment is first flowering date (FFD) in 2008. Though many species were observed for FFD, only 5 had any data recorded. Of those 5, all had data missing from at least one plot:

| Species | Number of plots without data |
|---------|------------------------------|
| DESCAE  | 1                            |
| ACOROS  | 4                            |
| CALLEP  | 8                            |
| BISBIS  | 12                           |
| CARSCO  | 33                           |

I think I will limit my analysis to only DESCAE or ACOROS, but have also considered using 4-5 of these species. *Do you have any thoughts?*

**Identify the scales of measurement for the predictor variables.**

The predictor variables are (1) nitrogen and (2) temperature. The predictor variables are

**Make a sketch of your sample design, showing the different grouping scales and predictor scales.**

![Figure 1. Schematic diagram of sampling design](figures/indiv-proj-method-schematic_kelsey-e.png)

2. Model
========

**Write down the statistical model for your design using mathematical notation.**

The first equation describes the distribution of FFD between plots within, where ![\\mu\_i](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_i "\mu_i") is the average FFD in plots of three types: (1) Block 1, (2) Block 2, (3) Block 3.
![y\_i \\sim Binomial(\\mu\_i)](http://chart.apis.google.com/chart?cht=tx&chl=y_i%20%5Csim%20Binomial%28%5Cmu_i%29 "y_i \sim Binomial(\mu_i)")

The second equation describes the distribution of FFD between plots within sub-blocks, where ![\\mu\_j](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_j "\mu_j") is the average FFD in plots of two types: (1) more snow, (2) less snow.
![\\mu\_i \\sim Binomial(\\mu\_{j})](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_i%20%5Csim%20Binomial%28%5Cmu_%7Bj%7D%29 "\mu_i \sim Binomial(\mu_{j})")

The average FFD (![\\mu\_j](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_j "\mu_j")) in more/less snow plots (![\\mu\_i](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_i "\mu_i")) is a function of whether the plot has extra N (![x\_1](http://chart.apis.google.com/chart?cht=tx&chl=x_1 "x_1")) or warmer temperatures (![x\_2](http://chart.apis.google.com/chart?cht=tx&chl=x_2 "x_2")):
![\\mu\_j = \\alpha\_j + \\beta\_1x\_1 + \\beta\_1x\_2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cmu_j%20%3D%20%5Calpha_j%20%2B%20%5Cbeta_1x_1%20%2B%20%5Cbeta_1x_2 "\mu_j = \alpha_j + \beta_1x_1 + \beta_1x_2")

The fourth model equation describes the relationship BETWEEN sub-blocks: ![\\alpha\_j \\sim Normal(\\bar{\\mu}\_a,\\sigma\_{\\alpha})](http://chart.apis.google.com/chart?cht=tx&chl=%5Calpha_j%20%5Csim%20Normal%28%5Cbar%7B%5Cmu%7D_a%2C%5Csigma_%7B%5Calpha%7D%29 "\alpha_j \sim Normal(\bar{\mu}_a,\sigma_{\alpha})")

where...
![\\alpha\_j](http://chart.apis.google.com/chart?cht=tx&chl=%5Calpha_j "\alpha_j") is the mean FFD in a sub-block for a specific type of plot (N, temp)
![\\bar{\\mu}\_a](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbar%7B%5Cmu%7D_a "\bar{\mu}_a") is the sub-block overall mean for a specific type of plot (N, temp)
![\\sigma\_{\\alpha}](http://chart.apis.google.com/chart?cht=tx&chl=%5Csigma_%7B%5Calpha%7D "\sigma_{\alpha}") is the variation in FFD between sub-blocks

**Write down the corresponding linear model formula for stan\_lmer/stan\_glmer.**

`ppfit <- stan_glmer(FFD ~ snow + nitrogen + temperature + (1|block), REML=FALSE, data = nwt_ffd_data, family = gaussian)`
