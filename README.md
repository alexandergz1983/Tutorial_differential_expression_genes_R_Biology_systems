## Tutorial para un proyecto de bioinformática oncológica: Análisis de Expresión Diferencial de Genes y Análisis para el descubrimiento de biomarcadores utilizando un conjunto de datos de cáncer colorrectal.

Este tutorial aborda el uso práctico de la bioinformática en un conjunto de datos reales, con el objetivo de proporcionar una experiencia realista en el análisis de expresión diferencial utilizando decenas de miles de transcriptos genéticas de tejidos normales y de cáncer colorrectal.

El tutorial se basa en R a tráves ipynb y/o RStudio como software para su implementación.


**origen de los datos:** https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Experiment%20Design

**Antecedentes:**

**A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients**

**https://doi.org/10.1016/j.molonc.2014.06.016**

Colorectal cancer (CRC) patients frequently experience disease recurrence and distant metastasis. This study aimed to identify prognostic indicators, including individual re- sponses to chemotherapy, in CRC patients. RNA-seq data was generated using 54 samples (normal colon, primary CRC, and liver metastases) from 18 CRC patients and genes associ- ated with CRC aggressiveness were identified. A risk score based on these genes was devel- oped and validated in four independent CRC patient cohorts (n 1⁄4 1063). Diverse statistical methods were applied to validate the risk scoring system, including a generalized linear model likelihood ratio test, KaplaneMeier curves, a log-rank test, and the Cox model. TREM1 and CTGF were identified as two activated regulators associated with CRC aggres- siveness. A risk score based on 19 genes regulated by TREM1 or CTGF activation (TCA19) was a significant prognostic indicator. In multivariate and subset analyses based on path- ological staging, TCA19 was an independent risk factor (HR 1⁄4 1.894, 95% CI 1⁄4 1.227e2.809, P 1⁄4 0.002). Subset stratification in stage III patients revealed that TCA19 had prognostic po- tential and identified patients who would benefit from adjuvant chemotherapy, regardless of age. The TCA19 predictor represents a novel diagnostic tool for identifying high-risk CRC patients and possibly predicting the response to adjuvant chemotherapy.

