Installation and usage
----------------------

The BhGLM package can be accessed at <https://github.com/nyiuab/BhGLM>.
It can be installed with `devtools` package as

    if(!require(devtools)) install.packages("devtools")
    library(devtools)
    install_github("nyiuab/BhGLM")

Or you can download BhGLM\_1.1.0.zip, and the source code can be
installed on your local machine.

    setwd()
    install.packages("BhGLM_1.1.0.zip", repos=NULL)
    library(BhGLM)

The Bayesian hierarchical negative binomial (HNB) model utilized in this
tutorial depends on the function bglm.

    bglm(formula, family = gaussian, data, offset, weights, subset, na.action,
        start = NULL, etastart, mustart, control = glm.control(epsilon = 1e-04, maxit = 50),
        prior = c("de", "t", "mde", "mt"), group = NULL, method.coef,
        dispersion = 1, prior.sd = 0.5, prior.scale = 0.5, prior.df = 1, prior.mean = 0, ss = c(0.04, 0.5),
        Warning = FALSE, verbose = TRUE, ...)

The arguments of bglm are explained as follows.

-   **formula, family, data, offset, weights, subset, na.action, start,
    etastart, mustart, control**: Functions same arguments in glm.
-   **family**: Identifies one of any of the standard families defined
    in glm or to the Negative Binomial distribution as NegBin or
    "NegBin".
-   **prior**: Distinguishes between several types of priors for the
    coefficents including double-exponetial ("de"), Student-t ("t"),
    spike-and-slab mixture double-exponential ("mde"), and
    spike-and-slab mixture t ("mt"). The mixture priors are only used
    for grouped predictors (defined by group). For ungrouped predictors,
    the priors are double-exponential or Student-t with scale =
    prior.scale.
-   **group**: Represents a numeric vector, an integer, or a list
    indicating the groups of predictors. If group=NULL, all the
    predictors form a single group. If group=K, the predictors are
    evenly divided into groups each with K predictors. If group is a
    numeric vector, it defines groups as follows: Group 1:
    (group\[1\]+1):group\[2\], Group 2: (group\[2\]+1):group\[3\], Group
    3: (group\[3\]+1):group\[4\], ..... If group is a list of variable
    names, group\[\[k\]\] includes variables in the k-th group.
-   **method.coef**: Determines whether to jointly update all
    coefficients or update coefficients group by group. The default is
    jointly updating. If method.coef=NULL or method.coef is missing, the
    function will jointly updating. If method.coef=K, the method will
    update K coefficients at a time. This argument can be a numeric
    vector or a list of variable names (as defined by group) that
    defines groups. If the number of coefficients is large, the
    group-by-group updating method can be much faster than the jointly
    updating.
-   **dispersion**: Identifies the dispersion parameter. If a value is
    provided, it is the starting value of dispersion. For Poisson and
    Binomial models, dispersion equals 1.
-   **prior.sd**: Represents the prior standard deviations in the normal
    priors of the coefficients. If provided, they are starting values of
    the prior standard deviations for the iterative algorithms.
-   **prior.scale**: Identifies the scale parameters in the priors, de
    and t. prior.scale = c(a1,a2,...,ak); if k &lt; J (J is the total
    number of predictors, not counting the intercept), it is internally
    expanded to c(a1,a2,...,ak, rep(ak,J-k)). For both t and de prior,
    smaller scale induces stronger shrinkage (weak shrinkage:
    prior.scale &gt; 1, strong shrinkage: prior.scale &lt; 0.3).
-   **prior.df**: Sets the prior degrees of freedom in the t and mt
    priors: default is 1 (leading to Cauchy prior)
-   **prior.mean**: Sets the prior mean for each coefficient.
-   **ss**: Represents a vector of two positive scale values for the
    mixture priors mde and mt, allowing for different scales for
    different predictors, leading to different amount of shrinkage.
    Smaller scale values give stronger shrinkage.
-   **warning**: A logical operator that shows the error messages of
    convergence and identifiability is set to TRUE.
-   **verbose**: A logical operator that prints the number of iterations
    and computational time when set to TRUE.
-   **...** All further arguments are derived from glm.

Simulation
----------

In order to assess the performance and the behavior of the HNB model,
the following R code was employed. The function glm.nb implemented in
the R package MASS was utilized as a means of comparison.

    ################  Simulation of Human Microbiome Count Data ####################
    ########### Comparison of glm.nb() and bglm() Modeling Techniques ##############

    #Amanda H. Pendegraft
    #R 3.5.0

    ################################################################################
    #Purpose of this R script:
    #(1) Simulate counts using NB distribution to represent a human microbiome dataset
    #    for multivariable analyses.  
    #(2) Apply glm.nb() and bglm() to compare performance.
    ################################################################################

    #Install packages.
    library(BhGLM)
    library(MASS)

    ################################################################################

    #Install functions.
    sim.out.updated <- function(coefs.p, coefs.est, alpha = c(0.05, 0.01))
    {
        p.if <- coefs.p
        b <- coefs.est
        power <- NULL
        for (j in 1:length(alpha)) {
            for (s in 1:ncol(coefs.p)) p.if[, s] <- (coefs.p[, s] <=
                alpha[j])
            power <- cbind(power, apply(p.if, 1, mean, na.rm=TRUE))
        }
        colnames(power) <- alpha
        b.mean <- apply(b, 1, mean, na.rm=TRUE)
        b.median <- apply(b, 1, median, na.rm=TRUE)
        b.l <- apply(b, 1, quantile, prob = 0.025, na.rm=TRUE)
        b.h <- apply(b, 1, quantile, prob = 0.975, na.rm=TRUE)
        est <- cbind(b.mean, b.median, b.l, b.h)
        colnames(est) <- c("mean", "median", "2.5%", "97.5%")
        out <- list(power = power, est = est)
        out
    }

    ################################################################################

    #Define simulation starting time.
    start.time=Sys.time()

    #Set working directory.
    setwd('/data/user/alhall91/simulation_HNB/')

    #Define parameters.
    niters=1000 #number of iterations
    R=rbind(c(-0.8, -0.5), c(-0.1, 0.1), c(0.5, 0.8)) #correlation
    N=c(50, 100, 200, 500) #number of samples
    M=c(1, 50, 100, 200) #number of covariates
    E=rbind(c(0.01, 0.15), c(0.20, 0.35), c(0.40, 0.55), c(0.00, 0.00)) #effect size
    S=c(0.01, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 2.00) #prior scale
    A=c(0.05) #significance level

    #Select parameter combination.
    comb=expand.grid(E=seq(1,4,1), M=seq(1,4,1),
        N=seq(1,4,1), R=seq(1,3,1))
    comb=comb[-which(comb$M!=1 & comb$E!=1),] #remove redundant combinations
    comb=comb[order(comb$E, decreasing=F),] #add combinations for k=1 type 1 error, seeds follow
    select=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

    #Starting seed.
    seed=(select-1)*1000+1
    #####seed=361001

    #Establish empty list to save results.
    results=vector("list", 33)
    names(results)=paste("res.", 1:33, sep="")

    #Run simulation.
    #####for (i in (1:25)){
    for (i in (1:niters)){

    #Set seed.
    set.seed(seed)

    #Generate coefficients.
    if (comb$M[select]==1) idx=c(1) else idx=c(15, 30, 45)
    if (comb$M[select]==1) B=runif(1, E[comb$E[select],1], E[comb$E[select],2]) else {
        B1=runif(1, E[1,1], E[1,2])
        B2=runif(1, E[2,1], E[2,2])
        B3=runif(1, E[3,1], E[3,2]) }
    Bvec=rep(0, M[comb$M[select]])
    if (comb$M[select]==1) Bvec[idx]=B else {
        Bvec[idx]=c(B1, B2, B3) }

    #Generate data.
    x=sim.x(n=N[comb$N[select]], m=M[comb$M[select]],
        corr=runif(1, R[comb$R[select],1], R[comb$R[select],2])); x=as.matrix(x)
    theta=runif(1, 0.1, 5); off=runif(N[comb$N[select]], 7.1, 10.5); T=exp(off)
    y=sim.y(x=x, mu=off-7,
        coefs=Bvec, p.neg=0.5, theta=theta); y=y$y.nb
    results[[31]]=cbind(results[[31]], c(Bvec))
    results[[32]]=cbind(results[[32]], c(theta))

    #Save data.
    setwd('/data/user/alhall91/data_HNB/')
    write.csv(cbind(y,x,T), paste(seed, "_HNB_data.csv", sep=""))
    setwd('/data/user/alhall91/simulation_HNB/')

    #Application of bglm() - nine prior scale values.
    adj_AIC=NULL
    for (j in (1:length(S))){
        f=NA; out.f=NA
        f=tryCatch(bglm(y~., data=data.frame(x), offset=off, family="NegBin", prior="t", prior.scale=S[j], verbose=T), error=function(e){NA})
        if (!is.na(f)) {
            out.f=summary(f)
            results[[j]]=cbind(results[[j]], c(out.f$coefficients[-1,1]))
            results[[j+9]]=cbind(results[[j+9]], c(out.f$coefficients[-1,4]))
            results[[j+18]]=cbind(results[[j+18]], c(f$nb.theta[1]))
            df=df.adj(f); AIC=2*df.adj(f)-2*f$loglik; adj_AIC[j]=AIC} else {
            results[[j]]=cbind(results[[j]], c(rep(NA, M[comb$M[select]])))
            results[[j+9]]=cbind(results[[j+9]], c(rep(NA, M[comb$M[select]])))
            results[[j+18]]=cbind(results[[j+18]], c(NA))
            adj_AIC[j]="NA"}}

    #Identify scale with best fit.
    adj_AIC=cbind(S, adj_AIC); best=adj_AIC[order(adj_AIC[,2]),1][1]
    results[[33]]=cbind(results[[33]], best)

    #Application of glm.nb(); prior scale not applicable.
    g=NA; out.g=NA
    g=tryCatch(glm.nb(y~offset(off)+., data=data.frame(x)), error=function(e){NA})
    if (!is.na(g)) {
        out.g=summary(g)
        if (length(which(!is.na(g$coefficients)))==M[comb$M[select]]+1) {
            results[[28]]=cbind(results[[28]], c(out.g$coefficients[-c(1),1]))
            results[[29]]=cbind(results[[29]], c(out.g$coefficients[-c(1),4]))
            results[[30]]=cbind(results[[30]], c(g$theta))} else {
            results[[28]]=cbind(results[[28]], c(rep(NA, M[comb$M[select]])))
            results[[29]]=cbind(results[[29]], c(rep(NA, M[comb$M[select]])))
            results[[30]]=cbind(results[[30]], c(NA))}} else {
            results[[28]]=cbind(results[[28]], c(rep(NA, M[comb$M[select]])))
            results[[29]]=cbind(results[[29]], c(rep(NA, M[comb$M[select]])))
            results[[30]]=cbind(results[[30]], c(NA))}

    #Update seed for subsequent iteration.  
    seed=seed+1

    print(i); print(seed)

    }

    #Summarization of frequency (power and false positive).
    out.1=cbind(sim.out.updated(results$res.10, results$res.1-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.11, results$res.2-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.12, results$res.3-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.13, results$res.4-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.14, results$res.5-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.15, results$res.6-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.16, results$res.7-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.17, results$res.8-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.18, results$res.9-results$res.31, alpha=A)$power,
        sim.out.updated(results$res.29, results$res.28-results$res.31, alpha=A)$power)
    colnames(out.1)=c("f0.01", "f0.05", "f0.10", "f0.15", "f0.25", "f0.50", "f0.75", "f1.00", "f2.00", "g")
    write.csv(out.1, file=paste(select, "_HNB_Freq_Coeff.csv", sep=""))

    #Summarization of coefficient estimation (mean only).
    out.2=cbind(sim.out.updated(results$res.10, results$res.1-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.11, results$res.2-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.12, results$res.3-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.13, results$res.4-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.14, results$res.5-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.15, results$res.6-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.16, results$res.7-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.17, results$res.8-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.18, results$res.9-results$res.31, alpha=A)$est[,1],
        sim.out.updated(results$res.29, results$res.28-results$res.31, alpha=A)$est[,1])
    colnames(out.2)=c("f0.01", "f0.05", "f0.10", "f0.15", "f0.25", "f0.50", "f0.75", "f1.00", "f2.00", "g")
    write.csv(out.2, file=paste(select, "_HNB_Bias_Coeff.csv", sep=""))

    #Summarization of dispersion estimation.
    out.3=cbind(t(sim.out.updated(results$res.19, results$res.19-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.20, results$res.20-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.21, results$res.21-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.22, results$res.22-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.23, results$res.23-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.24, results$res.24-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.25, results$res.25-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.26, results$res.26-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.27, results$res.27-results$res.32, alpha=A)$est),
        t(sim.out.updated(results$res.30, results$res.30-results$res.32, alpha=A)$est))
    colnames(out.3)=c("f0.01", "f0.05", "f0.10", "f0.15", "f0.25", "f0.50", "f0.75", "f1.00", "f2.00", "g")
    out.3=out.3[1,] #(mean only)
    write.csv(out.3, file=paste(select, "_HNB_Bias_Overdisp.csv", sep=""))

    #Summarization of AIC selection.
    write.csv(results$res.33, file=paste(select, "_HNB_adjAIC.csv", sep=""))

    #Define simulation stopping time; compute total simulation time.
    stop.time=Sys.time()
    time=round(difftime(stop.time, start.time, units="min"), 8)
    cat("The total simulation time was ", time, " minutes.")
    write.csv(time, paste(select, "_HNB_comp_time.csv", sep=""))

    #Erase workspace.
    rm(list=ls())

    ################################################################################

Application
-----------

An application of the HNB model was completed using the open-source,
open-access data published by the American Gut Project at
<ftp://ftp.microbio.me/AmericanGut/latest/11-packaged.zip>. The
following R code was employed to analysis those patients who
self-reported diagnosis of inflammatory bowel disease (IBD) by a medical
professional compared to those patients who self-reported no diagnosis
of inflammatory bowel disease (IBD) while adjusting for a mutlivariable
set of demographic, dietary, and systemic covariates.

    ###################  Application to American Gut Project #######################
    ################################ IBD Subset ####################################

    #Amanda H. Pendegraft
    #R 3.5.0

    ################################################################################

    #Purpose of this R script:
    #(1) To analyze subsets of AGP participants using the Bayesian HNB model as
    #    implemented in the bglm() function of the R package BhGLM in contrast to
    #    the classical NB model as implemented in the glm.nb() function of the R
    #    packages MASS.
    #(2) The FULL subsets focused on in this R script includes those patients diagnosed
    #    by a medical professional with IBD compared to those not diagnosed by a
    #    medical professional with IBD while adjusting for a high-dimensional set
    #    of covariates.  
    #(3) The REDUCED subsets focused on in this R scripts decreased the sample size of
    #    the FULL subsets using propensity score matching as implemented in the R
    #    package MatchIt.
    #(3) For the consideration of IBS, change IBD to IBS in the following R script.

    ################################################################################

    #Load packages.
    library(phyloseq)
    library(MatchIt)
    library(tableone)
    library(BhGLM)
    library(bmzim)
    library(ggplot2)

    ###############################################################################

    #Set working directory.
    setwd()

    #Read 16S rRNA data.
    #Selected dataset is fecal, untrimmed, unrarefied.
    counts=read.csv('ag_May_18_2017_fecal_notrim_unrarefied_species.csv', header=T, stringsAsFactors=F)
    row.names(counts)=counts[,1]; counts=counts[,-1]

    #Read metadata.
    clinic=read.csv('ag_May_18_2017_fecal_notrim_unrarefied_metadata_CLEAN.csv', header=T, stringsAsFactors=F)
    row.names(clinic)=clinic[,1]; clinic=clinic[, -1]

    #Split counts from taxonomy.
    otu=t(counts[,1:(dim(counts)[2]-7)]); otu=otu[row.names(clinic),]
    tax=counts[,(dim(counts)[2]-6):(dim(counts)[2])]

    #Define as phyloseq object.
    otu=otu_table(otu, taxa_are_rows=F)
    tax=tax_table(as.matrix(tax))
    clinic=sample_data(data.frame(clinic))
    agp1=phyloseq(otu, tax, clinic); agp1

    #Observe unique taxa at each phylogenentic level.
    unlist(lapply(apply(tax_table(agp1), 2, unique), length))

    ###############################################################################

    #Select samples meeting inclusion/exclusion criteria.
    agp2=subset_samples(agp1,
        IBD!='Diagnosed by an alternative medicine practitioner' &
        IBD!='Self-diagnosed' &
        IBD!='Unspecified' &
        IBD!='Unknown' &
        SEX!='unspecified' &
        SEX!='other' &
        RACE!='Unspecified' &
        ANTIBIOTIC_HISTORY=='I have not taken antibiotics in the past year.' &
        AGE_CORRECTED>=18 & AGE_CORRECTED<=69 &
        COUNTRY_RESIDENCE=='United States' &
        CANCER=='I do not have this condition' &
        DIABETES=='I do not have this condition'); agp2

    #Observe the distribution of IBD.
    table(sample_data(agp2)$IBD)

    #Redefine IBD.
    #Yes=Diagnosed by a medical professional (doctor, physician assistant).
    #No=I do not have this condition.
    idx=which(sample_data(agp2)$IBD=="I do not have this condition")
    sample_data(agp2)$IBD=1; sample_data(agp2)$IBD[idx]=0

    #Propensity score matching to control sample size.
    data=as(sample_data(agp2), "data.frame")
    data=data[,which(colnames(data) %in% c("IBD", "AGE_CORRECTED", "RACE", "SEX", "BMI_CORRECTED"))]
    data=data[-which(is.na(data$BMI_CORRECTED) || is.na(data$SEX) || is.na(data$RACE) ||
      is.na(data$AGE_CORRECTED)),]
    set.seed(0505)
    mdata=matchit(IBD~BMI_CORRECTED+SEX+RACE+AGE_CORRECTED, data, method="nearest", ratio=1); summary(mdata)

    #Extract REDUCED subset using mdata.
    idx=NULL
    idx[1:nrow(mdata$match.matrix)]=row.names(mdata$match.matrix)
    idx[(nrow(mdata$match.matrix)+1):(2*nrow(mdata$match.matrix))]=mdata$match.matrix[,1]
    REDUCED=subset_samples(agp2, row.names(sample_data(agp2)) %in% idx); REDUCED

    #And, FULL subset.
    FULL=agp2

    #Select FULL or reduced before proceeding.
    #agp3=FULL
    agp3=REDUCED

    #Select covariates.
    #Comment out 2nd, 17th, and 18th lines below if using REDUCED.
    covs=c('IBD',
    #   'AGE_CORRECTED', 'SEX', 'BMI_CORRECTED', 'RACE',
        'ALCOHOL_FREQUENCY','ALCOHOL_TYPES_RED_WINE','ALCOHOL_TYPES_WHITE_WINE',
        'ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL','ALCOHOL_TYPES_BEERCIDER',
        'ALCOHOL_TYPES_SOUR_BEERS', 'DRINKING_WATER_SOURCE',
        'ONE_LITER_OF_WATER_A_DAY_FREQUENCY','SUGAR_SWEETENED_DRINK_FREQUENCY',
        'POULTRY_FREQUENCY','SEAFOOD_FREQUENCY','RED_MEAT_FREQUENCY',
        'HIGH_FAT_RED_MEAT_FREQUENCY','CONSUME_ANIMAL_PRODUCTS_ABX', 'WHOLE_EGGS',
        'FRUIT_FREQUENCY','VEGETABLE_FREQUENCY', 'FERMENTED_PLANT_FREQUENCY',
        'MILK_CHEESE_FREQUENCY','LACTOSE','MILK_SUBSTITUTE_FREQUENCY','FROZEN_DESSERT_FREQUENCY',
        'WHOLE_GRAIN_FREQUENCY', 'GLUTEN', 'DIET_TYPE',
        'HOMECOOKED_MEALS_FREQUENCY','READY_TO_EAT_MEALS_FREQUENCY','PREPARED_MEALS_FREQUENCY',
        'OLIVE_OIL', 'SALTED_SNACKS_FREQUENCY','ARTIFICIAL_SWEETENERS','SUGARY_SWEETS_FREQUENCY',
        'MULTIVITAMIN', 'VITAMIN_B_SUPPLEMENT_FREQUENCY', 'VITAMIN_D_SUPPLEMENT_FREQUENCY',
        'OTHER_SUPPLEMENT_FREQUENCY','SMOKING_FREQUENCY', 'PROBIOTIC_FREQUENCY' )
    agp4=subset_samples(agp3, select=covs); agp4
    #sample_data(agp4)$AGE_CORRECTED=as.numeric(sample_data(agp4)$AGE_CORRECTED)
    #sample_data(agp4)$BMI_CORRECTED=as.numeric(sample_data(agp4)$BMI_CORRECTED)

    #Observe unique taxa at each phylogenentic level.
    unlist(lapply(apply(tax_table(agp4), 2, unique), length))

    #Observe descriptive statistics.
    table1=CreateTableOne(vars=colnames(sample_data(agp4))[c(2:length(covs))], strata=c('IBD'),
        testApprox=chisq.test, testNormal=oneway.test,
        data=data.frame(sample_data(agp4))); summary(table1); table1

    ###############################################################################

    #Agglomerate at other levels for analysis.
    #Comment out all except level of interest.
    agp5=tax_glom(agp4, taxrank="Phylum"); agp5
    #agp5=tax_glom(agp4, taxrank="Class"); agp5
    #agp5=tax_glom(agp4, taxrank="Order"); agp5
    #agp5=tax_glom(agp4, taxrank="Family"); agp5
    #agp5=tax_glom(agp4, taxrank="Genus"); agp5
    #agp5=tax_glom(agp4, taxrank="Species"); agp5

    #Analyze only those taxa with mean relative abundance >0.1%.
    #Comment out all except level of interest.
    agp7=transform_sample_counts(agp5, function(x) x / sum(x))
    agp8=filter_taxa(agp7, function(x) mean(x) > 1e-3, TRUE)
    agp7=prune_taxa(colnames(otu_table(agp8)), agp5)
    agp7=subset_taxa(agp7, Phylum!='p__'); agp7 #Exclude non-identified taxa.
    #agp7=subset_taxa(agp7, Class!='c__'); agp7
    #agp7=subset_taxa(agp7, Order!='o__'); agp7
    #agp7=subset_taxa(agp7, Family!='f__'); agp7
    #agp7=subset_taxa(agp7, Genus!='g__'); agp7
    #agp7=subset_taxa(agp7, Species!='s__'); agp7

    #Observe unique taxa at each phylogenentic level.
    unlist(lapply(apply(tax_table(agp7), 2, unique), length))

    #Establish indices mapping to levels of interest.
    #Comment out all except level of interest.
    idx=2
    #idx=3
    #idx=4
    #idx=5
    #idx=6
    #idx=7

    #Bayesian HNB model.
    #Select best prior scale using adjusted AIC.  
    #Comment out 3rd and 4th lines below if using REDUCED.
    S=c(0.01, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 2.00)
    x=covariates(x.cat=as(sample_data(agp7)[,1:39], 'data.frame'), cat.center=F, fill.missing=T)
    #x=covariates(x.con=as(sample_data(agp7)[,c(2,4)], 'data.frame'),
    #  x.cat=as(sample_data(agp7)[,c(1,3,5:43)], 'data.frame'), con.rescale=T, cat.center=F, fill.missing=T)
    off=log(rowSums(otu_table(agp1))); off=off[which(names(off) %in% row.names(sample_data(agp5)))]
    fplot=as.data.frame(matrix(0, nrow=ntaxa(agp7), ncol=5))
    colnames(fplot)=c('Scale', 'Estimate', 'SE', 'Stat', 'Pval')
    AIC=matrix(NA, nrow=ntaxa(agp7), ncol=9)
    for (i in (1:ntaxa(agp7))) {
        adjAIC=NULL
        y0=as.vector(otu_table(agp7)[,i])
        for (s in (1:length(S))){
            f=NA
            tryCatch( {
              f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[s],verbose=F) },
              error=function(e){NA} )
            if (!is.na(f)[1]) { if (f$iter<50 && length(which(is.na(f$coefficients)==T))==0) {
              adjAIC=c(adjAIC, 2*df.adj(f)-2*f$loglik) } else {
                  adjAIC=c(adjAIC, NA) } } else {
                    adjAIC=c(adjAIC, NA) } }
        AIC[i,1:9]=adjAIC
        if (length(which(is.na(adjAIC)))==9) { fplot[i,1:5]=rep(NA, 5) } else {
            s=order(adjAIC)[1]
            tryCatch( {
                f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[s], verbose=F)},
                  error=function(e){NA} )
            name=as.vector(tax_table(agp7)[which(row.names(tax_table(agp7))==colnames(otu_table(agp7)[,i])),idx])
            if (name%in%row.names(plot)) {name=paste(name, '__', i, sep='')}
            stats=summary(f)$coefficients[which(row.names(summary(f)$coefficients)=="IBD_1"),]
            row.names(fplot)[i]=name; fplot[i,1]=S[s]; fplot[i,2:5]=stats } }

    #Print adjust AIC statistics and Bayesian HNB modeling results.
    row.names(AIC)=tax_table(agp7)[,idx]; print(AIC)
    fplot$FDR=p.adjust(fplot$Pval, "fdr"); print(fplot)

    #Classical NB model.
    #Does not converge given when covariates exceed samples (i.e. REDUCED subset for IBD).
    gplot=as.data.frame(matrix(0, nrow=ntaxa(agp7), ncol=5))
    colnames(gplot)=c('Scale', 'Estimate', 'SE', 'Stat', 'Pval')
    AIC=vector('numeric', length=ntaxa(agp7))
    for (i in (1:ntaxa(agp7))){
        y0=as.vector(otu_table(agp7)[,i]); head(y0)
        g=NA
        g=tryCatch(glm.nb(y0~offset(off)+., data=x), error=function(e){NA})
        if (!is.na(g)[1]) { if (g$iter<50 && length(which(is.na(g$coefficients)==T))==0) {
            AIC[i]=g$aic
            name=as.vector(tax_table(agp7)[which(row.names(tax_table(agp7))==colnames(otu_table(agp7)[,i])),idx])
            if (name%in%row.names(gplot)) { name=paste(name, '__', i, sep='') }
            stats=summary(g)$coefficients[which(row.names(summary(g)$coefficients)=="IBD_1"),]
            row.names(gplot)[i]=name; gplot[i,1]=NA; gplot[i,2:5]=stats
            print(i) } else {
                gplot[i,1:5]=rep(NA, 5)
                AIC[i]=NA}} else {
                    gplot[i,1:5]=rep(NA, 5)
                    AIC[i]=NA } }

    #Print adjust AIC statistics and classical NB modeling results.
    row.names(AIC)=tax_table(agp7)[,idx]; print(AIC)
    gplot$FDR=p.adjust(gplot$Pval, "fdr"); print(gplot)

    ###############################################################################

    #After executing at all levels, save each set of results and plot those taxa
    #found to be significantly differentially abundant.

    #Identify significant taxa (i.e. P-value<0.10).
    #Rename files read in 4th line to reflect naming system.
    data2=NULL
    level=c('species', 'genus', 'family', 'order', 'class', 'phylum')
    for (i in (1:length(level))) {
        data1=read.csv(paste('HNB_IBD_application_', level[i], '_reduced.csv', sep=''), header=T,
          stringsAsFactors=F)
        data2=rbind(data2, data1[which(data1[,6]<0.10),]) }
    data2=data2[nrow(data2):1,]
    data1=data2[which(is.na(data2$Scale)),]; row.names(data1)=data1[,c(1)]; data1=data1[,-c(1)]
    data2=data2[which(!is.na(data2$Scale)),]; row.names(data2)=data2[,c(1)]; data2=data2[,-c(1)]

    #And, plot.
    par(mar=c(5,10,4,5)+.1)
    print(plot.coefs(res=as.matrix(data2[,c(2,3,5)]), cex.var=0.75, gap=1.25, cex.pts=0.75))

    ###############################################################################

Contact information
-------------------

Please contact Amanda H Pendegraft at <alhall91@uab.edu> with any
questions regarding the simulation and application code of the Bayesian
HNB tutorial.
