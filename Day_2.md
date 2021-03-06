When we say "the analysis", we usually mean "the statistics". The
visualization can tell us about the trends, but nothing about their
relevance on a level of the general population. Ok, let's start.

### prepare the data

    library(vegan)
    phenotypes = read.table("./Phenotypes.txt",sep="\t",header=T,row.names=1)
    microbiome = read.table("./Microbiome.txt",sep="\t",header=T,row.names=1)
    microbiome = t(microbiome)
    microbiome = microbiome[rownames(phenotypes),]

Module 1. Alpha diversity
=========================

I hope you already understand the concept of alpha and beta-diversity.
let's calculate these metrics on species level.

calculate alpha diversity
-------------------------

    library(vegan)
    microbiome_species = microbiome[,grep("t__",colnames(microbiome),invert = T)]
    microbiome_species = microbiome_species[,grep("s__",colnames(microbiome_species))]
    colnames(microbiome_species) = sub(".*[|]s__","",colnames(microbiome_species))
    alpha_diversity_species = diversity(microbiome_species,index = "shannon")
    summary(alpha_diversity_species)

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.827   2.380   2.580   2.548   2.758   3.061

**Q1. get summary statistics of alpha diversity on genus level**

Are the phenotypes associated with alpha diversity? If alpha diversity
is normally (or at least symmetrically) distributed, we can use simple
parametric methods to explore that.

    shapiro.test(alpha_diversity_species)

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  alpha_diversity_species
    ## W = 0.97697, p-value = 0.8892

    hist(alpha_diversity_species)

![](Day_2_files/figure-markdown_strict/unnamed-chunk-3-1.png)

The distribution OK, and Shapiro test said it's normal. Let's run linear
regression for the subset of our phenotypes

    summary(lm(alpha_diversity_species ~ Age + BMI + PPI + DiagnosisCurrent,data = phenotypes))

    ## 
    ## Call:
    ## lm(formula = alpha_diversity_species ~ Age + BMI + PPI + DiagnosisCurrent, 
    ##     data = phenotypes)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.63553 -0.17488  0.07349  0.16341  0.45529 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          2.921120   0.547869   5.332 8.38e-05 ***
    ## Age                  0.005500   0.005705   0.964    0.350    
    ## BMI                 -0.021678   0.023880  -0.908    0.378    
    ## PPI                 -0.058463   0.209732  -0.279    0.784    
    ## DiagnosisCurrentIBD -0.133727   0.151460  -0.883    0.391    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3132 on 15 degrees of freedom
    ## Multiple R-squared:  0.1736, Adjusted R-squared:  -0.04674 
    ## F-statistic: 0.7879 on 4 and 15 DF,  p-value: 0.5508

**Q2. Make a boxplot of alpha diversity on different diagnosis**

Module 2. Beta diversity
========================

The concept of beta diversity is some kind of distance or dissimilarity
between samples. We can represent beta diversity as a distance matrix.
Being just a square matrix, it cannot be visualised directly, but we can
use RDA (redundancy analysis) to project the distances in two dimentions

    dist_euclidean = dist(microbiome_species)
    rda_euclid = dbrda(dist_euclidean ~ 0)
    plot(rda_euclid)

![](Day_2_files/figure-markdown_strict/unnamed-chunk-5-1.png)

However, but default, **dist** function uses euclidean distance metric.
It's known to be not appropriate for the data which have a lot of zeros
(zero-inflated). **Vegan** package contains a special function called
**vegdist** which contains several dissimilarity measures, some of them
are much more appropriate for microbiome data than euclidean distance.
By default, it uses Bray-Curtis dissimilarity.

    dist_bray = vegdist(microbiome_species,method = "bray")

**Q3. Make RDA plot for Bray-Curtis dissimilarity matrix**

Now let's do the statistical analysis on beta-diversity. We want to test
if phenotype has an effect on between-sample distances: the samples with
closer phenotype values might have more similar bacterial composition,
which means that they will be closer by distance. Let's do it for gender
phenotype.

    gender_test = adonis(dist_bray ~ phenotypes$Sex,permutations = 1000)
    gender_test

    ## 
    ## Call:
    ## adonis(formula = dist_bray ~ phenotypes$Sex, permutations = 1000) 
    ## 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## phenotypes$Sex  1    0.2116 0.21156 0.60884 0.03272 0.8412
    ## Residuals      18    6.2545 0.34747         0.96728       
    ## Total          19    6.4661                 1.00000

Now let's do it on gender and age simultaneously

    gender_age_test = adonis(dist_bray ~ phenotypes$Sex + phenotypes$Age,permutations = 1000)
    gender_age_test

    ## 
    ## Call:
    ## adonis(formula = dist_bray ~ phenotypes$Sex + phenotypes$Age,      permutations = 1000) 
    ## 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
    ## phenotypes$Sex  1    0.2116 0.21156 0.62819 0.03272 0.8472  
    ## phenotypes$Age  1    0.5294 0.52944 1.57210 0.08188 0.0969 .
    ## Residuals      17    5.7251 0.33677         0.88540         
    ## Total          19    6.4661                 1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

for *adonis* function, the order matters! try to run formulas: MB ~ age
+ gender and MB ~ gender + age. The variance explained by each variable
will be different between these two options. The idea is the variables
are correlated, the one which comes first takes as much variance of
distance matrix as possible. So, if we particularly interested in one
phenotype, correcting for all others, the phenotype of interest should
come last.

*Adonis* results are standard R objects, so you can directly extract the
values from summary tables. The parameter which contains matrix with
results is *aov.tab*

    gender_age_test$aov.tab[1,]

    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                Df SumsOfSqs MeanSqs F.Model       R2 Pr(>F)
    ## phenotypes$Sex  1   0.21155 0.21155 0.62819 0.032718 0.8472

**Q4 run beta-diveristy analysis for euclidean distance instead of
Bray-Curtis for age and gender**

**Q5 run beta-diversity analysis for disease state, correcting for all
other phenotypes**

Module 3. Statistical tests on particular bacterial features.
=============================================================

Both alpha- and beta diversity based associations give us only a general
picture: we test if phenotype has a large impact on overall bacteria
composition. But for interpretation it might be much more interesting if
disease is associated with particular bacterial taxa, such as species or
families. To find these bacteria, we should run some kind of statistical
test on every bacterial taxa present in the dataset. Let's do this
analysis for age phenotype. We'll use Spearman correlation, because
bacterial traits are poorly distributed, so parametic methods might be
affected by outliers, ties and other distribution artefacts. However, we
too many bacterial species in the dataset, so it make sense to remove
low abundant ones. let's take only species which are present in more
than 5 samples. Don't forget about adjustment for multiple testing
correction!

    library(foreach)
    microbiome_species_filtered = microbiome_species[,colSums(microbiome_species>0) > 5 ]
    correlation_with_age = foreach(i=1:ncol(microbiome_species_filtered),.combine = rbind)%do%{
      correlation_result = cor.test(phenotypes$Age,microbiome_species_filtered[,i],method = "spearman")
      data.frame(species = colnames(microbiome_species_filtered)[i],cor = correlation_result$estimate,pvalue = correlation_result$p.value)
    }
    correlation_with_age$adjP = p.adjust(correlation_with_age$pvalue,method = "BH")

**Q6 make the association with BMI. Use the same cutoff. Report 6
top-significant results **

We are mainly interested in IBD. It's two-level factor variable, so we
can use Wilcoxon test for that

**Q7 run the association between disease and species. Report the results
with FDR&lt;0.1 in a data.frame comprised of 1) name of the species;
2)association p-value; 3) adjusted p-value. Sort the results by adjusted
p-value, from smallest to largest. **

with ranked methods like Spearman correlation and Wilcoxon test, we
can't correct testing variables for covariates. How can we deal with it?
One of the ways to go is to use a proper transformation to make
variables distribution more normal. One of the commonly used
transformations for compositional data is arcsin square root
transformation. It only allows values from 0 to 1 as an input. Our data
is in percentages, so let's fix it and do the transformation.

    microbiome_species_filtered = microbiome_species_filtered/100
    microbiome_species_transformed = asin(sqrt(microbiome_species_filtered))

**Q8 Run test for normality (shapiro.test) for transformed and
untransformed datasets. Make a plot of -log(pvalues)**

Now we can try to use linear regression to correct for covariates. To
make it convenient, let's write a function for that.

    get_p_value = function(x) {
      lm1 = lm(x ~ .,data = phenotypes)
      aov.tab1 = summary(lm1)$coef
      RowIndex = grep("DiagnosisCurrent",rownames(aov.tab1))
      return(aov.tab1[RowIndex,4])
    }

    pvalues_transformed_species = apply(microbiome_species_transformed,2,get_p_value)

**Q9 run this analysis on untransformed species. Make a plot of
-log(p-values) for transformed and untransformed datasets.**

**Q10 modify the function to get not only p-values but also
coefficients. Plot coefficients for transformed data against
untransformed.**

**Q11 perform FDR for transformed results. Report a list of species
which are associated with disease at FDR&lt;0.01**
