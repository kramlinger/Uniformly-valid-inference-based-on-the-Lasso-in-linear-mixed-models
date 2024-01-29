# Uniformly Valid Inference Based on the Lasso in Linear Mixed Models

This repository contains the code that generates figures in the article ['Uniformly Valid Inference Based on the Lasso in Linear Mixed Models'](https://www.sciencedirect.com/science/article/pii/S0047259X23000763). 

---

### Real Data Study 

(R version 4.3.1, file: `realdata.R`)

The findings in Section 6 replicate the study of [Opsomer et. al](https://academic.oup.com/jrsssb/article/70/1/265/7109365). Our study re-scaled the acid-neutralizing capacity parameter to harmonize both covariates to a similar scale. 

I can not share the underlying data (`use2.csv`). 

---

### Simulation Study

(Julia version 1.6.0, file: `simulation.jl`)

Although the setting allows parallel implementation, this has not been done. The file generates three .csv files containing the data displayed in Section 5. 