# Introduction
This is the github for storing related functions, Shiny APP and codes for methodologies on Living Longer Probability (LLP).

# User Guide for Shiny APP of LLP Estimation and Graphical Display
This Shiny application performs LLP estimation using Kaplan Meier-based nonparametric and Cox model-based parametric methods with graphical display.

Under nonparametric method it provides:

- KM-based two-sample survival probability curves
- Half-tie LLP estimates
- Greenwood variance and logit-transformed 95% confidence intervals
- Kaplan-Meier survival curves with median survival times
- Unadjusted hazard ratios with 95% confidence intervals
- Log-rank p-values

Under semiparametric method it provides:

- Cox-based survival probability curves
- Category-specific Cox LLP and hazard ratios
- Overall Cox LLP for a binary effect modifier
- Confidence interval for the overall Cox LLP
- Optional Cox-model covariate adjustment
- Optional merging of two or three original categories into two binary groups

**Install the required R packages:**
install.packages(c("shiny", "dplyr", "ggplot2", "survival"))

The application loads:
library(shiny)
library(dplyr)
library(ggplot2)
library(grid)
library(survival)

**In RStudio:**

1. Open `app.R`.
2. Click **Run App**.
3. The application will open in the RStudio Viewer or a web browser.

Input file requirements:

The application accepts a CSV file.

The maximum upload size of the CSV file is approximately 500 MB

**The CSV must contain:**

1. A treatment-arm variable (Variable type: Categorical, At least 2 distinct non-missing categories)
2. A survival-time variable (Variable type: Numeric, non-negative)
3. A binary event/censoring variable (Variable type: Binary)

Optional variables include:

4. A categorical effect modifier (Variable type: categorical)
5. One or more Cox-model covariates (Variable type: Numeric or categorical)

**Treatment arm variable**

After uploading the CSV, select:

- the treatment-arm column
- the control-arm label
- the treatment-arm label


Only observations belonging to the selected control and treatment labels are analyzed.


**Survival time variable**

The time variable must be numeric.

Rows with missing, nonnumeric, or negative survival times are removed.

The unit may be days, months, or years, but it must be consistent across all observations.



**Event/censoring variable**

The status variable must use values `0` and `1`.

The app supports either convention.

Option 1:

1 = censored, 0 = event

Choose:

1 = censored, 0 = event


Option 2:

1 = event, 0 = censored

Choose:

1 = event, 0 = censored

Rows with status values other than `0` or `1` are removed.


**Analysis-method selection**

The user may select either or both:

- Non-parametric
- Semiparametric

**Nonparametric output**
The app generates:

- KM-based two-sample survival probability curves
- LLP estimates with 95% confidence intervals
- Kaplan-Meier survival curves
- Median survival times
- Unadjusted hazard ratios with 95% confidence intervals
- Log-rank p-values

**Semiparametric output**

The app generates:

- Cox-based survival probability curves
- Cox LLP estimates with 95% confidence intervals
- Hazard ratios with 95% confidence intervals
- Category-specific results when an effect modifier is selected
- Overall Cox LLP for a binary effect modifier
- Optional covariate-adjusted Cox results

**Effect modifier**

Check:

“Include categorical effect modifier”

Then:
1.	Select the effect-modifier column.
2.	Automatically all the categories will be selected 
3.	Select combine two categories to compute a combined LLP for the two categories


**Covariate adjustment**

Covariate adjustment applies to the semiparametric Cox analysis.
Check:

“Adjust Cox model for covariate(s)”
Then select one or more covariates.
