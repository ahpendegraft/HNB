Installation and usage
----------------------

The BhGLM package can be accessed at <https://github.com/nyiuab/BhGLM>.
After downloading BhGLM\_1.1.0.zip, the source code can be installed on
your local machine.

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

    ################  Simulation of High Dimensional Microbiome ####################
    ########### Comparison of glm.nb() and bglm() Modeling Techniques ##############

    #Amanda H. Pendegraft
    #R 3.5.0

    ################################################################################

    #Purpose of this R script: 
    #(1) Simulate counts using NB distribution to represent a high-dimensional microbiome 
    #    dataset.  Four small, moderate, or large non-zero effects embedded into dataset.
    #(2) Apply glm.nb() and bglm() to compare performance.

    ################################################################################

    #Install packages.
    library(BhGLM)
    library(MASS)

    ################################################################################

    #Define simulation starting time. 
    start.time=Sys.time()

    #Set working directory. 
    setwd('/data/scratch/alhall91/simulation_NB/small')
    #setwd('/data/scratch/alhall91/simulation_NB/moderate')
    #setwd('/data/scratch/alhall91/simulation_NB/large')

    #Define parameters.
    niters=1000 #number of iterations given each combination
    N=c(50, 100, 200, 500) #number of samples
    P=c(1, 50, 100, 200) #number of covariates
    S=c(0.01, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 2.00) #prior scale values
    C=rbind(c(0.5, 0.8), c(-0.1, 0.1), c(-0.8, -0.5)) #correlation among covariates

    #Establish object to link parameters to results files.
    key=NULL

    #Generate combinatios of prior scale, number of covariates, and sample size.
    comb=expand.grid(S=S, P=P, N=N)

    #Select a combination based on Slurm array task id.
    select=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    n=comb$N[select]; p=comb$P[select]; s=comb$S[select]

    #Define a seed.
    seed=1+300*(select-1)+432000

    #Run simulation looping around strong negative, weak, and strong positive correlation. 
    for (c in (1:nrow(C))){

    res1=NULL; res2=NULL; res3=NULL; res4=NULL; res5=NULL

    for (i in (1:niters)){

    #Set seed.
    set.seed(seed) 

    #Simulate data.
    r=runif(1, C[c,1], C[c,2]) 
    x=sim.x(n=n, m=p, corr=r) 
    x=as.matrix(x)
    if (p==1) {idx=c(1)} 
    if (p!=1) {idx=c(12,24,36,48)}
    if (p==1) {B=runif(1, 0.01, 0.15)} #small non-zero effects
    #if (p==1) {B=runif(1, 0.20, 0.35)} #moderate non-zero effects
    #if (p==1) {B=runif(1, 0.40, 0.55)} #large non-zero effects
    if (p!=1) {
        B=runif(4, 0.01, 0.15) #small non-zero effects
    #   B=runif(4, 0.20, 0.35) #moderate non-zero effects
    #   B=runif(4, 0.40, 0.55) #large non-zero effects
        B=B*c(-1,1,-1,1)}
    Bvec=rep(0,p); Bvec[idx]=B
    theta=runif(1, 0.1, 5) 
    off=runif(n, 7.1, 10.5) 
    mu=-7 
    mun=off+mu
    y=sim.y(x=x, mu=mun, coefs=Bvec, p.neg=0.5, theta=theta) 
    y0=y$y.nb 

    #Fit standard NB GLM with glm.nb() from MASS.
    g=tryCatch(glm.nb(y0~offset(off)+., data=data.frame(x)),
        error=function(e){NA})
    if (is.na(g)==F) {
        if(is.na(g$coefficients[p+1]==T)) {g=NA}}
        
    #Fit Bayesian HNB model with bglm() from BhGLM.
    f=tryCatch(bglm(y0~., data=data.frame(x), offset=off, family="NegBin", prior="t", prior.scale=s, verbose=F),
        error=function(e){NA})  
    if (is.na(f)==F) {
        if(is.na(f$coefficients[p+1]==T)) {f=NA}}
        
    #Save results. 
    if (is.na(g)==T & is.na(f)==T){
        res1=cbind(res1, c(theta, NA, NA))}
    if (is.na(g)==F & is.na(f)==T){
        res1=cbind(res1, c(theta, g$theta, NA))
        res4=cbind(res4, summary(g)$coefficients[-c(1),4])
        res2=cbind(res2, summary(g)$coefficients[-c(1),1]-Bvec)}
    if (is.na(g)==T & is.na(f)==F){
        res1=cbind(res1, c(theta, NA, f$nb.theta))
        res5=cbind(res5, summary.bh(f)[-c(1),3])
        res3=cbind(res3, summary.bh(f)[-c(1),1]-Bvec)}
    if (is.na(g)==F & is.na(f)==F){
        res1=cbind(res1, c(theta, g$theta, f$nb.theta))
        res4=cbind(res4, summary(g)$coefficients[-c(1),4])
        res5=cbind(res5, summary.bh(f)[-c(1),3])
        res2=cbind(res2, summary(g)$coefficients[-c(1),1]-Bvec)
        res3=cbind(res3, summary.bh(f)[-c(1),1]-Bvec)}
        
    #Increase seed for next combination of parameters.
    seed=seed+1 

    } 

    #Compute power, type 1 error rate, and coefficient estimate bias.
    if (length(res2)==0) {
        results1=list(c('Classical NB model does not converges for this parameter combination.'), 
            c('Classical NB model does not converges for this parameter combination.'))
        names(results1)=c('power', 'est')
        results2=sim.out(coefs.p=res5, coefs.est=res3, alpha=c(0.05, 0.01, 0.005, 0.001))
        res2=c('Classical NB model does not converges for this parameter combination.')
        res4=c('Classical NB model does not converges for this parameter combination.')
    } else {
        results1=sim.out(coefs.p=res4, coefs.est=res2, alpha=c(0.05, 0.01, 0.005, 0.001))
        results2=sim.out(coefs.p=res5, coefs.est=res3, alpha=c(0.05, 0.01, 0.005, 0.001))}

    #And, save.
    write.csv(res1, paste(select, '_', c, '_theta.csv', sep=''))
    write.csv(res2, paste(select, '_', c, '_standardized_estimates_glmnb.csv', sep=''))
    write.csv(res3, paste(select, '_', c, '_standardized_estimate_bglm.csv', sep=''))
    write.csv(res4, paste(select, '_', c, '_pvalue_glmnb.csv', sep=''))
    write.csv(res5, paste(select, '_', c, '_pvalue_bglm.csv', sep=''))
    write.csv(results1$power, paste(select, '_', c, '_power_glmnb.csv', sep=''))
    write.csv(results2$power, paste(select, '_', c, '_power_bglm.csv', sep=''))
    write.csv(results1$est, paste(select, '_', c, '_bias_glmnb.csv', sep=''))
    write.csv(results2$est, paste(select, '_', c, '_bias_bglm.csv', sep=''))

    #Also, save a key associated with this select file.
    if (res2=='Classical NB model does not converges for this parameter combination.'){
        key=data.frame(t(c(select, n, p, s, 0, ncol(res3))))
    } else {
        key=data.frame(t(c(select, n, p, s, ncol(res2), ncol(res3))))}
    colnames(key)=c('select', 'sample size', 'number of coefficients', 
        'scale', 'valid glm.nb', 'valid bglm')
    write.csv(key, paste(select, '_', c, '_key.csv', sep=''))

    }

    #Define simulation stopping time; compute total simulation time.
    stop.time=Sys.time()
    time=round(difftime(stop.time, start.time, units="min"), 8)
    cat("The total simulation time was ", time, " minutes.")

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
of inflammatory bowel disease (IBD) while adjusting for a
high-dimensional set of demographic, dietary, and systemic covariates.

    ###################  Application to American Gut Project #######################
    ################################ IBD Subset ####################################

    #Amanda H. Pendegraft
    #R 3.5.0

    ################################################################################

    #Purpose of this R script: 
    #(1) To analyze a subset of AGP participants using the Bayesian HNB model as 
    #    implemented in the bglm() function of the R package BhGLM.  
    #(2) The subset focused on in this R scripts includes those patients diagnosed
    #    by a medical professional with IBD compared to those not diagnosed by a 
    #    medical professional with IBD while adjusting for a high-dimensional set 
    #    of covariates.  
    #(3) For the consideration of IBS, change IBD to IBS in the following R script. 

    ################################################################################

    #Load packages.
    library(phyloseq)
    library(tableone)
    library(BhGLM)
    library(bmzim)

    ###############################################################################

    #Set working directory.
    setwd()

    #Microbiome count data.
    #Selected dataset is fecal, trimmed at 100 bp, unrarefied.
    #This dataset has already been agglomerated at the Genus level.
    counts=read.csv('ag_May_18_2017_fecal_100nt_unrarefied_genus.csv', header=T, stringsAsFactors=F)
    row.names(counts)=counts[,1]; counts=counts[,-1]

    #Clinical data.
    #This dataset has been cleaned. 
    clinic=read.csv('ag_May_18_2017_fecal_100nt_unrarefied_metadata_CLEAN.csv', header=T, stringsAsFactors=F)
    row.names(clinic)=clinic[,1]; clinic=clinic[, -1]

    #Split microbiome count data from taxonomy.
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

    #Select dietary variables to be included as covariates.  
    covs=c('AGE_CORRECTED','BMI_CORRECTED','SEX','RACE', 'CENSUS_REGION', 'IBD', 
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
    agp3=subset_samples(agp2, select=covs); agp3
    sample_data(agp3)$AGE_CORRECTED=as.numeric(sample_data(agp3)$AGE_CORRECTED)
    sample_data(agp3)$BMI_CORRECTED=as.numeric(sample_data(agp3)$BMI_CORRECTED)
    sample_data(agp3)$RACE[which(sample_data(agp3)$RACE=='Asian or Pacific Islander')]='Other'
    sample_data(agp3)$RACE[which(sample_data(agp3)$RACE=='Hispanic')]='Other'

    #Observe descriptive statistics. 
    table1=CreateTableOne(vars=colnames(sample_data(agp3))[c(1:5, 7:44)], strata=c('IBD'), 
        testApprox=chisq.test, testNormal=oneway.test,
        data=data.frame(sample_data(agp3)))
    #print(table1)

    ###############################################################################

    #Analyze only those taxa with mean relative abundance >0.1%.
    agp4=transform_sample_counts(agp3, function(x) x / sum(x))
    agp5=filter_taxa(agp4, function(x) mean(x) > 1e-3, TRUE)
    agp4=prune_taxa(colnames(otu_table(agp5)), agp3)

    #Observed unique taxa at each phylogenentic level. 
    unlist(lapply(apply(tax_table(agp5), 2, unique), length))

    #Application of bglm. 
    S=c(0.01, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 2.00)
    x=covariates(x.con=as(sample_data(agp4)[,1:2], 'data.frame'), 
        x.cat=as(sample_data(agp4)[,3:44], 'data.frame'), 
        con.rescale=F, cat.center=F, fill.missing=T)
    off=log(rowSums(otu_table(agp2)))
    plot=as.data.frame(matrix(0, nrow=ntaxa(agp4), ncol=5))
    colnames(plot)=c('Scale', 'Estimate', 'SE', 'Stat', 'Pval')
    for (i in (1:ntaxa(agp4))){
    adjAIC=NULL
    for (s in (1:length(S))){
    tryCatch({
        y0=as.vector(otu_table(agp4)[,i]); head(y0)
        f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[s], verbose=F)
        },error=function(e){NA})
    adjAIC=c(adjAIC, 2*df.adj(f)-2*f$loglik)}
    idx=order(adjAIC)[1]
    tryCatch({
        y0=as.vector(otu_table(agp4)[,i])
        f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[idx], verbose=F)
        },error=function(e){NA})
    name=as.vector(tax_table(agp4)[which(row.names(tax_table(agp4))==colnames(otu_table(agp4)[,i])),6])
    if (name%in%row.names(plot)) {name=paste(name, '__', i, sep='')}
    stats=summary(f)$coefficients[which(row.names(summary(f)$coefficients)=="`IBD_I do not have this condition`"),]
    row.names(plot)[i]=name; plot[i,1]=S[idx]; plot[i,2:5]=stats
    print(i)}

    #Print the results.
    #print(plot)

    ###############################################################################

    #Agglomerate at other phylogenetic levels for analysis.
    agp5=tax_glom(agp3, taxrank="Phylum"); agp5
    #agp5=tax_glom(agp3, taxrank="Class"); agp5
    #agp5=tax_glom(agp3, taxrank="Order"); agp5
    #agp5=tax_glom(agp3, taxrank="Family"); agp5

    #Analyze only those taxa with mean relative abundance >0.1%.
    agp6=transform_sample_counts(agp5, function(x) x / sum(x))
    agp7=filter_taxa(agp6, function(x) mean(x) > 1e-3, TRUE)
    agp6=prune_taxa(colnames(otu_table(agp6)), agp5)

    #Observe unique taxa at each phylogenentic level. 
    unlist(lapply(apply(tax_table(agp6), 2, unique), length))

    #Application of bglm.
    S=c(0.01, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 1.00, 2.00)
    x=covariates(x.con=as(sample_data(agp6)[,1:2], 'data.frame'), 
        x.cat=as(sample_data(agp7)[,3:44], 'data.frame'), 
        con.rescale=F, cat.center=F, fill.missing=T); dim(x)
    off=log(rowSums(otu_table(agp2)))
    plot=as.data.frame(matrix(0, nrow=ntaxa(agp6), ncol=5))
    colnames(plot)=c('Scale', 'Estimate', 'SE', 'Stat', 'Pval')
    for (i in (1:ntaxa(agp6))){
    adjAIC=NULL
    for (s in (1:length(S))){
    tryCatch({
        y0=as.vector(otu_table(agp6)[,i]); head(y0)
        f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[s], verbose=F)
        },error=function(e){NA})
    adjAIC=c(adjAIC, 2*df.adj(f)-2*f$loglik)}
    idx=order(adjAIC)[1]
    tryCatch({
        y0=as.vector(otu_table(agp7)[,i]); head(y0)
        f=bglm(y0~., data=x, offset=off, family="NegBin", prior="t", prior.scale=S[idx], verbose=F)
        },error=function(e){NA})
    name=as.vector(tax_table(agp6)[which(row.names(tax_table(agp6))==colnames(otu_table(agp6)[,i])),2])
    #name=as.vector(tax_table(agp6)[which(row.names(tax_table(agp6))==colnames(otu_table(agp6)[,i])),3])
    #name=as.vector(tax_table(agp6)[which(row.names(tax_table(agp6))==colnames(otu_table(agp6)[,i])),4])
    #name=as.vector(tax_table(agp6)[which(row.names(tax_table(agp6))==colnames(otu_table(agp6)[,i])),5])
    if (name%in%row.names(plot)) {name=paste(name, '__', i, sep='')}
    stats=summary(f)$coefficients[which(row.names(summary(f)$coefficients)=="`IBD_I do not have this condition`"),]
    row.names(plot)[i]=name; plot[i,1]=S[idx]; plot[i,2:5]=stats
    print(i)}

    #Print the results.
    #print(plot)

    ###############################################################################

    #After executing at all phylogenetic levels, save each set of results and plot those taxa
    #found to be significantly differentially abundant.

    #Identify significant taxa.
    data2=NULL
    level=c('Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum')
    for (i in (1:length(level))){
      data1=read.csv(paste(level[i], ' Level IBD Statistics.csv', sep=''), header=T, stringsAsFactors=F)
      data2=rbind(data2, data1[which(data1[,6]<0.05),])}
    row.names(data2)=data2[,1]; data2=data2[,-1]

    #And, plot.
    #pdf('IBD Plot.pdf', width=6.693, height=4.676)
    par(mar=c(5,10,4,5)+.1)
    print(plot.coefs(res=as.matrix(data2[,c(2,3,5)]), 
        cex.var=0.75, gap=1.25, cex.pts=0.75))
    #dev.off()

    ###############################################################################

Contact information
-------------------

Please contact Amanda H Pendegraft at <alhall91@uab.edu> with any
questions regarding the simulation and the application of the HNB
tutorial.
