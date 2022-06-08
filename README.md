# PhD Project 1: Data‐driven staging of genetic frontotemporal dementia using multi‐modal MRI 

In this project, we applied an unsupervised machine learning algorithm, the contrastive trajectory inference (cTI), to MRI data from individuals with genetic forms of frontotemporal dementia (FTD) to obtain individual scores of disease stage.

Brain changes begins many years prior to symptom onset in genetic FTD, and disease presentation and pathology varies widely. There are treatments under development for specific genetic forms. Asymptomatic carriers of FTD-causing genes will almost certainly develop the disease in future, meaning these individuals could be included in clinical trials. However, since genetic FTD is relatively rare, clinical trails will likely need to combine participants at different disease stages, including asymptmatic and symptomatic individuals. Methods to stage the disease during both the presymptomatic and symptomatic phases are needed for the development of clinical trials outcomes. Measures from MRI can detect brain changes in asymptomatic gene-carriers on a group level, but high variance exists across individuals. Combining these measures into a unified disease staging system may provide a more accurate staging method.

The cTI analyzes temporal patterns in multi-dimentional, large scale population datasets and, from this, calculates individual scores of disease stage. Here, we applied the cTI to five measures of brain structure and function from MRI data, derived from individuals who have FTD causing genetic mutations, including those with an FTD diagnosis (symptomatic), those who do not yet have symptoms (presymptomatic), and a control group of non-carriers. We compared the cTI-obtained disease scores in all gene carriers (presymptomatic and symptomatic) to clinical and neuropsychological test scores (measuring behavioral symptoms, attention, memory, language, and executive functions) as well as the estimated years to symptom onset (EYO: age - mean age of onset in relatives). 

The cTI based disease scores were significantly correlated with all clinical and neuropsychological tests and with the EYO. Mean cTI scores were significantly higher in the presymptomatic carriers than in the control group, indicating that the method may capture subtle pre-dementia cerebral changes, although this change was not replicated in a subset of subjects with complete MRI data. This study is a proof of concept that the cTI can identify data-driven disease stages in a heterogeneous sample of genetic FTD combining different genetic mutations and disease stages using only MRI meeasures.


[You can find the paper here](https://onlinelibrary.wiley.com/doi/10.1002/hbm.25727)

**In this repository you will find code for extracting the cTI scores from the model's output, conducting the statistical analyses, and creating the figures seen in the paper.**

