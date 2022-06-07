# iPGS-construction
We present the codes used in [the GIGASTROKE project](https://doi.org/10.21203/rs.3.rs-1175817/v1) to construct integrative polygenic score (iPGS) models.


## Derivation and evaluation of iPGS models for Europeans and East Asians
In the GIGASTROKE project, iPGS models for Europeans and East Asians were derived from 10 GIGASTROKE GWAS and 12 GWAS of vascular risk traits. 

First, from the genome-wide summary statistics for each GWAS and a linkage disequilibrium (LD) reference panel of the European subjects (n=503) and East-Asian subjects (n=504) from the 1000 Genomes Project, 37 candidate PGS models were computed using P+T, [LDpred](https://github.com/bvilhjal/ldpred), and [PRScs](https://github.com/getian107/PRScs) algorithms. Then, the best PGS model was selected for each GWAS, where the best model was defined as the model that showed the maximal area under the curve (AUC) in the Europeans (European case-control data with 1,003 any ischemic stroke [AIS] cases and 8,997 controls) and East Asian (East-Asian case-control data with 577 any ischemic stroke [AIS] cases and 9,232 controls) model training dataset. 

Second, among the 22 selected PGS models derived from the 22 GWAS, 11 and 7 were significantly associated with AIS in the Europans and East Asian model training dataset respectively (Bonferroni-corrected P<0.05). The significant PGS models were used as the variables for elastic-net logistic regression and the weights for the variables were trained using the model training dataset. The European iPGS model consisting of 1,213,574 variants and an East Asian iPGS model consisting 6,010,730 variants were constructed by combining the 11 and 7 significant PGS models using the elastic-net derived weights respectively. 

Lastly, the European and East Asian iPGS models were evaluated in the European (a European prospective cohort data with 102,099 subjects including 1,128 incident IS cases) and East Asian (an East-Asian case-control data with 1,470 AIS cases and 40,459 controls) model evaluation datasets.

<img width="1082" alt="Derivation and evaluation of iPGS models for Europeans and East Asians" src="https://user-images.githubusercontent.com/31917903/172354128-bbfb3dcf-bcee-4415-a9ec-29cbca8364c8.png">
AS indicates any stroke; AIS, any ischemic stroke; LAS, large artery stroke; SVS, small vessel stroke; CES, cardioembolic stroke; AF, atrial fibrillation; CAD, coronary artery disease; T2D, type 2 diabetes; SBP, systolic blood pressure; DBP, diastolic blood pressure; TC, total cholesterol; LDL-C, low-density lipoprotein cholesterol; HDL-C, high-density lipoprotein cholesterol; TG, triglyceride; BMI, body mass index; AUC indicates area under the curve; EUR, Europeans; EAS, East-Asian; GWAS, genome-wide association study; LD, linkage disequilibrium; PGS, polygenic score.

## About this repository
This repository presents the sample code used in the second step (i.e., the construction of iPGS models). 
Note that the code requires the access to the individual-level data, and therefore, was modified file paths for each research institute. 

## `training-weights.R`
This code is written in R. 
This code loads a binary outcome, polygenic scores, covariates, and principal components from files. 
Multiple polygenic scores may be loaded to combine them into a iPGS model. 
Our code calculates a weight for each input polygenic model to construct iPGS model. 
The results will be output to four files:
- `scaling-factors.txt`
- `metaGRS.auc-best.model.beta`
- `metaGRS.auc-best.model.param`
- `glm.beta`

### Input score file format
The code assumes that polygenic score values are stored in the [sscore file format (produced by PLINK2 --score command)](https://www.cog-genomics.org/plink/2.0/formats#sscore). Our code recognizes the 4-th column (`SCORE1_AVG`) as polygenic scores. 

```text:example-of-sscore-format
#IID	NMISS_ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG
TEST_NA18939	2890	807.933	0.00270411
TEST_NA18940	2890	790.845	0.00315385
TEST_NA18941	2890	843.059	0.0027847
TEST_NA18942	2890	815.158	0.00216979
TEST_NA18943	2890	822.601	0.00232627
TEST_NA18944	2890	768.808	0.00397438
TEST_NA18945	2890	810.82	0.00217487
TEST_NA18946	2890	794.639	0.00308176
TEST_NA18947	2890	846.24	0.00218907
```

## `construct-iPGS.pl`
This Perl code combines multiple input PGS models into an iPGS model using the outputs of `training-weights.R`. 



