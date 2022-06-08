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
Note that the code requires the access to individual-level data, and therefore, file paths should be modified for each research institute. 

## `training-weights.R`
This code is written in R. 
This code loads a binary outcome, polygenic scores, covariates, and principal components from files. 
Multiple polygenic scores may be loaded to combine them into a iPGS model. 
Our code calculates a weight for each input polygenic model to construct iPGS model. 
The results will be output to four files:
- `scaling-factors.txt`
- `iPGS.auc-best.model.beta`
- `iPGS.auc-best.model.param`
- `glm.beta`

### Input score file format
The code assumes that polygenic score values are stored in the [sscore file format (produced by PLINK2 --score command)](https://www.cog-genomics.org/plink/2.0/formats#sscore). Our code recognizes the 4-th column (`SCORE1_AVG`) as polygenic scores. 

The below shows an example of the sscore fromat.
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

The below shows the usage of the program:
```
./construct-iPGS.pl \
  -i modified-iPGS.auc-best.model.beta \
  -s modified-scaling-factors.txt \
  -p list-of-variant-IDs.txt \
  -o output-model-file.txt.gz
```

### modified-iPGS.auc-best.model.beta
The below is an example of `iPGS.auc-best.model.beta`, which is one of the outputs of `training-weights.R`.
```text:iPGS.auc-best.model.beta
GIGASTROKE_AS_ASN	0.0444241457543297
GIGASTROKE_IS_ASN	0.0727897407227377
GIGASTROKE_AS_TRANS	0.0386965622551286
GIGASTROKE_IS_TRANS	0.0764857062983009
GIGASTROKE_SVS_TRANS	0.0325457745331939
bbj_SBP	0.0739281603318368
bbj_DBP	0.0500664791243444
sex	-0.136504634115316
age	0.757487766844012
PC1_AVG	-0.131346375829411
PC2_AVG	-0.00586400049684333
PC3_AVG	-0.0387839816683539
PC4_AVG	0
PC5_AVG	0
```

This file should be manually modified to use the program as follows:
- Remove lines for covariates (such as sex, age, and PCs)
- Substitute PGS model names by file paths for the PGS models
```text:modified-iPGS.auc-best.model.beta
/path/to/GIGASTROKE_AS_ASN.txt.gz	0.0444241457543297
/path/to/GIGASTROKE_IS_ASN.txt.gz	0.0727897407227377
/path/to/GIGASTROKE_AS_TRANS.txt.gz	0.0386965622551286
/path/to/GIGASTROKE_IS_TRANS.txt.gz	0.0764857062983009
/path/to/GIGASTROKE_SVS_TRANS.txt.gz	0.0325457745331939
/path/to/bbj_SBP.txt.gz	0.0739281603318368
/path/to/bbj_DBP.txt.gz	0.0500664791243444
```

### modified-scaling-factors.txt
The below is an example of `scaling-factors.txt`, which is one of the outputs of `training-weights.R`.
```text:scaling-factors.txt
GIGASTROKE_AS_ASN	0.0712152989452889
GIGASTROKE_IS_ASN	0.0714588270767169
GIGASTROKE_AS_TRANS	0.0486802253903617
GIGASTROKE_IS_TRANS	0.783356528059269
GIGASTROKE_SVS_TRANS	0.0131431025100374
bbj_SBP	0.178284845396525
bbj_DBP	0.161732484927453
sex	0.497816506122276
age	16.8341154160173
PC1_AVG	0.0162684577597879
PC2_AVG	0.0092715521245179
PC3_AVG	0.00733092631571973
PC4_AVG	0.00723126676335968
PC5_AVG	0.00693940113154949
```

This file should be manually modified to use the program as follows:
- Remove lines for covariates (such as sex, age, and PCs)
- Substitute PGS model names by file paths for the PGS models
```text:modified-scaling-factors.txt
/path/to/GIGASTROKE_AS_ASN.txt.gz	0.0712152989452889
/path/to/GIGASTROKE_IS_ASN.txt.gz	0.0714588270767169
/path/to/GIGASTROKE_AS_TRANS.txt.gz	0.0486802253903617
/path/to/GIGASTROKE_IS_TRANS.txt.gz	0.783356528059269
/path/to/GIGASTROKE_SVS_TRANS.txt.gz	0.0131431025100374
/path/to/bbj_SBP.txt.gz	0.178284845396525
/path/to/bbj_DBP.txt.gz	0.161732484927453
```

### list-of-variant-IDs.txt
The code assumes an variant ID format as `<chr_name>:<chr_position>:<other_allele>:<effect_allele>`.
The below shows an example of the variant IDs list file:
```text:example-of-variant-ID-list-file.txt
1:10539:C:A
1:10642:G:A
1:11008:C:G
1:11012:C:G
1:11063:T:G
1:13011:T:G
1:13110:G:A
1:13116:T:G
1:13118:A:G
1:13259:G:A
```

All variant IDs that are appeared in input iPGS models should be listed in this file. 

### output-model-file.txt.gz
The output file format is compatible with [the PGSCatalog format specification](https://www.pgscatalog.org/downloads/).
The below shows an example of the `output-model-file.txt.gz`. 

```text:example-of-output-model-file.txt.gz
chr_name	chr_position	reference_allele	effect_allele	effect_weight
1	108506	C	T	0.000380215245804818
1	109503	G	A	-0.000122842666960445
1	133855	C	T	-0.000520302949248368
1	135195	A	G	0.000455288216452582
1	136113	C	T	-0.000583210515841391
1	138041	G	A	-0.000130481848325896
1	138396	G	A	0.000441454621042695
1	138484	A	G	0.0013867688822883
1	138817	T	C	-0.000875354754123417
```
