# BCA — BCKM 1980-2014 Replication Reference (Parts III & IV)

This file is the replication-target portion of the paper summary —
cross-country findings (§3.B) and the US 1980-2014 empirical tables
(Technical Appendix). It is the regression-test ground truth for the
Layer-1 BCKM replication and is **not** the live methodology spec.

For the live methodology and calibration that apply to **any**
country/period the toolkit is run on (Parts I and II — the prototype
economy, the accounting procedure, the application details, data
sources, and the calibration constants), see the root
[`BCA_info.md`](../BCA_info.md).

## Contents

- **Part III — Findings: Cross-Country Application** (main paper §3.B)
  - Cross-country evidence on the Great Recession and the 1982 recession
- **Part IV — Empirical Results: United States Replication**
  (Technical Appendix)
  - Tables 3-7: data construction, parameter choices
  - Tables 8, 9, 10: P matrix, P₀ vector, Q Cholesky factor
  - Table 11: f-statistics — the headline regression target
  - Tables 12-16: peak-trough decomposition, Great Recession
    counterfactuals, second-moment statistics

The numerical pins from these tables — Table 11 f-statistics, Table 12
peak-trough decomposition, Tables 8/9/10 (P, P₀, Q) — are the values
the BCKM-replication tests
([`tests/test_bckm_table12.py`](../tests/test_bckm_table12.py),
parts of [`tests/test_bckm_reference.py`](../tests/test_bckm_reference.py))
hold the pipeline accountable to.

---

# Part III — Findings: Cross-Country Application

*From Brinca, Chari, Kehoe, and McGrattan — "Accounting for Business Cycles," NBER Working Paper No. 22663, §3.B. The setup of §3 in the paper (functional forms, calibration constants, the consumer-durables and sales-tax treatments, and the OECD trend-removal step) duplicates Part I §3 above and is not repeated here. The paragraphs below pick up where §3.A in the paper leaves off.*

Now we describe the results of applying our procedure to OECD countries for the Great Recession and the 1982 recession. Here we focus primarily on the fluctuations due to the efficiency, labor, and investment wedges.[^1]

[^1]: We alert the reader that the quantitative results for Spain should be treated with caution. In some robustness analysis for Spain, we found that the nonlinear labor wedge computed directly from the consumer's first-order condition (4) moved substantially more than the labor wedge computed using our log-linearization procedure.

## The Great Recession

Here we discuss our findings for the 24 OECD countries. The main finding is that in terms of accounting for the downturn, in the United States the labor wedge is by far the most important, in Spain, Ireland, and Iceland the investment wedge is the most important, and in the rest of the countries, the efficiency wedge is the most important.

### Three Illustrative Recessions

Here we illustrate our findings for one country for which the efficiency wedge, labor wedge, and investment wedge, respectively, is the most important. In reporting our findings, we remove a country-specific trend from output, investment, and the government consumption wedge. Both output and labor are normalized to equal 100 in the base period 2008:1. Here we focus primarily on the fluctuations due to the efficiency, labor, and investment wedges. We discuss the government consumption wedge and its components in our Appendix.

#### France: Primarily an Efficiency Wedge Recession

Removed since it's not relevant to the exercise.

#### United States: Primarily a Labor Wedge Recession

Next consider the United States. In Figure 2, panel A, we see that output and labor both fell about 7% from 2008:1 to 2009:3 while investment fell about 23%. In Figure 2, panel B, we see that the efficiency wedge fell very modestly by only about 1%, while the labor wedge and the investment wedge both worsened dramatically, by about 8% and 9%, respectively. In Figure 2, panels C, D, and E, we see that the labor and investment wedges play the most important role in accounting for the downturn in output and labor, while the investment wedge accounts for the bulk of the downturn in investment.

Overall, considering the period from 2008 until the end of 2011, these results imply that the Great Recession in the United States should be thought of as primarily a labor wedge recession, with an important secondary role for the investment wedge. This finding implies that the most promising models must yield significant fluctuations in the labor wedge, with some role for the investment wedge. Models that emphasize the efficiency wedge are less promising.[^2]

[^2]: In the Appendix we show that if we estimated the stochastic process for the wedges from 1948 to 2015, the contribution of the labor wedge rises and that of the investment wedge falls. A similar change occurs if we decrease the investment adjustment cost parameter.

#### Ireland: Primarily an Investment Wedge Recession

Removed since it's not relevant to the exercise.

### Summary Statistics for our OECD Countries

So far we have described the Great Recession in three countries. Here we describe useful summary statistics over the period 2008:1 to 2011:3. One such statistic, referred to as the φ statistic, is intended to capture how closely a particular component, say, the output component due to the efficiency wedge, tracks the underlying variable, say, output. For our decomposition of output, we let

$$\phi^Y_i = \frac{1 / \sum_t (y_t - y_{it})^2}{\sum_j \left(1 / \sum_t (y_t - y_{jt})^2 \right)},$$

where y_{it} is the output component due to wedge *i* = (A, τ_l, τ_x, g). We compute similar statistics for labor and investment. The φ statistic has the desirable feature that it lies in [0, 1], sums to one across the four wedges, and when a particular output component tracks output perfectly — in that if (y_t − y_{it}) = 0 for all *t* — then φ^Y_i = 1, that is, the φ statistic for the wedge reaches its maximum value of 1. Note that this statistic is the inverse of the mean-square error for each wedge appropriately scaled so that the sum across wedges adds to one.

Now consider our main finding. In Figure 4, panel A, we display the φ statistic for the efficiency wedge and labor wedge components of output. The downward-sloping lines represent combinations for which the sum of the labor wedge and efficiency wedge components is constant at 70% and 90%, respectively. This figure shows that the United States stands out from the other countries in that the labor wedge accounts for a much greater fraction of the movements in output than it does in any other country. Specifically, the labor wedge accounts for about 46% of the movements in output in the United States but no more than 22% in any other country. In all other countries except Iceland, Ireland and Spain, the efficiency wedge accounts for roughly 50% or more of the movements in output. In Table 1, we report the decompositions of output, labor, and investment for all countries. There we see that for Iceland, Ireland, New Zealand, and Spain, the investment wedge accounts for 51%, 48%, 42%, and 82% of the movements in output, respectively. In the other panels of Figure 4, we display the φ statistics for the components of labor and investment.

Our main finding is also apparent if we use other ways to measure how important a given wedge is for the movements in output, labor, and investment. When we discussed the three illustrative recessions earlier, we compared simple peak-to-trough measures of output, labor, and investment to the corresponding measures for each of the components. In Tables 2A, 2B, and 2C, we report such measures for all of our countries. A quick perusal of these measures shows that they give the same message as the φ statistics do. Consider France, for example. The φ statistic indicates that the efficiency wedge accounts for the bulk of the movements in output, namely, about 92% of its decline. The peak-to-trough measure indicates that the efficiency wedge also accounts for the bulk of the peak-to-trough decline, namely about 5.9% of the 6.5% decline, or about 91% of the decline.

### Comparing the Great Recession with Recessions of the Early 1980s

The postwar era had essentially two periods during which most developed economies experienced recessions at roughly the same time: the early 1980s and the Great Recession of 2008. Here we compare the recessions of the early 1980s with the Great Recession. For the United States, we use the NBER business cycle dates; for the OECD countries, we use the business cycle dates as estimated by ECRI when available and otherwise use the CEPR Euro Area Business Cycle Dates. We use the stochastic process for wedges estimated over the 1980–2014 period for both episodes. (See the Appendix for details.)

In Figure 5, panel A, we compare the φ statistics for the efficiency wedge component of output for the two recessions. This panel shows that for most of the countries, the efficiency wedge in the Great Recession played a more important role than it did during the recessions of the 1980s. In Figure 5, panel B, we compare the φ statistics for the labor wedge component of output for the two recessions. This panel shows that in the Great Recession, the labor wedge accounts for over 40% of the fluctuations in output only in the United States, while in the 1982 recession it does so only in Belgium, the United Kingdom, and France. In Figure 5, panel C, we compare the φ statistics for the investment wedge component of output for the two recessions. This panel shows that in most of the countries, the investment wedge played a larger role in the recessions of the 1980s than it did during the Great Recession.

In Table 3 we report the φ statistics for the 1982 recessions. This table shows that the efficiency wedge played the most important role for ten countries, the labor wedge for three countries, and the investment wedge for seven countries. Together with Table 1, this table broadly reinforces our two main findings for the comparison. First, the labor wedge played an important role for output in the Great Recession only for the United States, and in the 1982 recession it played a dominant role only in Belgium, France, the United Kingdom, and New Zealand. Second, for most countries, in the Great Recession the efficiency wedge played a more important role and the investment wedge played a less important role than they did in the recessions of the 1980s.

In Table 4, panels A, B, and C, we report peak-to-trough results for the 1982 recession. Comparing Table 3 with the panels of Table 4, we see that the peak-to-trough results present the same overall picture as our φ statistics do. If we compare the classification of the most important wedge for each country using φ statistics for output to that using the peak-to-trough decline for output, we see that they agree in all but three cases.

### Summary Statistics for the Entire Period

In Tables 5A, 5B, and 5C, we present some summary statistics for the entire period 1980:1–2014:3 about the importance of the various wedges in accounting for the movements in output, labor, and investment. In Table 5A, for example, we report the standard deviation of the output component due to each wedge relative to the standard deviation of output during the entire period, along with the correlation of each such output component with output. In Tables 5B and 5C, we report similar statistics for labor and its components and for investment and its components.

Using these statistics to infer the importance of various wedges is more subtle than using the φ statistics. The φ statistic captures in one statistic how much the component due to a wedge moves, as well as how closely this component tracks the underlying variable. Instead, to evaluate the importance of a wedge using the statistics in this table, we need to jointly consider the relative standard deviations and the correlations.

Consider France, for example. Viewing the relative standard deviations alone suggests that the labor and investment wedges play roughly the same role in accounting for the movement in output. Indeed, the relative standard deviations of the labor and investment components of output are 93% and 92%, respectively. But the correlations of these variables with output suggests that the investment wedge plays a much more important role. Indeed, the labor component of output comoves negatively with output, whereas the investment component of output comoves positively with output.

With this perspective in mind, the averages across countries show that the efficiency wedge plays the most important role in accounting for output. The standard deviation of the efficiency component of output is 92% of output, and its correlation with output is 0.77. Even though the labor component of output is 89% as variable as output itself, it is essentially uncorrelated with output. In this sense, the labor wedge does not account for much of the movements in output.

### The Importance of the Classification of Consumer Durables

Macroeconomists have long argued that theory implies it is appropriate to treat the expenditures on consumer durables as a form of investment that yields a flow of consumption services. This treatment requires adjustments to the national income accounts classification of consumption and investment to make them consistent with the theory.

Here we show that while this adjustment is quantitatively important for some countries, for most countries it does not change the overall findings. In Figure 6, panel A, we contrast the φ statistic for the efficiency wedge component of output when this consistent adjustment is made and when it is not. Clearly, the countries with statistics most affected by this adjustment are Iceland and Spain. In Iceland, for example, the contribution of the efficiency wedge falls from 26% when durables are correctly accounted for to 12% when they are not. In Spain, the contribution of the efficiency wedge increases from 11% when durables are correctly accounted for to 29% when they are not.

In Figure 6, panels B and C, we contrast the analogous φ statistics for the labor wedge component of output and for the investment wedge component of output. In panel C we see that in Iceland and Spain, the contribution of the investment wedge to output is 51% and 82% when durables are correctly accounted for and 65% and 35% when they are not.

### Comparing our Procedure with a Perfect Foresight Procedure

Some authors implement a perfect foresight version of our procedure in which agents have perfect foresight about the future evolution of the wedges. The equilibrium conditions for the deterministic version of our prototype model are

$$(68)\quad c_t + x_t + g_t = y_t,$$

$$(69)\quad y_t = A_t F(k_t, (1 + \gamma)^t l_t),$$

$$(70)\quad \frac{-U_{lt}}{U_{ct}} = [1 - \tau_{lt}]\, A_t (1 + \gamma)^t F_{lt},\text{ and}$$

$$(71)\quad U_{ct}[1 + \tau_{xt}] = \beta U_{ct+1}\{A_{t+1} F_{kt+1} + (1 - \delta)[1 + \tau_{xt+1}]\}.$$

Clearly, the efficiency wedge, the labor wedge, and the government consumption wedge can be recovered from the static relationships in (68), (69), and (70). Recovering the investment wedge, however, requires solving the difference equation implied by the Euler equation (71). To do so we need to impose either an initial condition or a terminal condition. In practice, we imposed an initial condition that the investment wedge begins at zero.

In Figure 7, panels A, B, and C, we plot the φ statistics for the perfect foresight procedure against the same statistics for our procedure. These panels show that for a significant number of the countries, the φ statistics are very different. In particular, the perfect foresight procedure greatly exaggerates the importance of the labor wedge for the United States and Spain. Under perfect foresight, the labor wedge accounts for 92% and 72% of the movements in output for the United States and Spain, while under the standard business cycle accounting procedure, the labor wedge accounts for only 46% and 5%, respectively.

We highlight two important sources for these differences. One is that in the perfect foresight procedure, private agents anticipate the evolution of future wedges perfectly and thus react in the current period to actual future worsening or improvement of the wedges. In this sense, the perfect foresight procedure brings with it all the undesirable properties of the simple "news" models by which an anticipated worsening of, say, the labor wedge leads to a current boom as households choose to increase labor supply before times worsen. The other is that, as we noted earlier, the perfect foresight procedure uses the nonlinear version of the first-order conditions (68)–(71) to compute the wedges, whereas our procedure uses log-linearized versions of these conditions.

---

# Part IV — Empirical Results: United States Replication

*The following sections come from the Technical Appendix and apply the procedure of Parts I–III to U.S. data 1980Q1–2014Q4. They are the calibration target the codebase must replicate (CLAUDE.md → "Goal"). Section numbering 6/7/8 is preserved from the source.*

## 6. Whole-Sample Tables

**Table 3: Properties of the wedges**

| Country        | σA/σy | στl/σy | στx/σy | σg/σy | ρA,y  | ρτl,y | ρτx,y | ρg,y  |
|----------------|-------|--------|--------|-------|-------|-------|-------|-------|
| Australia      | 0.79  | 1.53   | 1.77   | 5.34  | 0.60  | 0.04  | 0.72  | -0.40 |
| Austria        | 0.92  | 1.10   | 1.52   | 2.83  | 0.83  | -0.01 | 0.35  | 0.12  |
| Belgium        | 0.68  | 1.51   | 1.01   | 3.23  | 0.68  | 0.38  | -0.35 | -0.37 |
| Canada         | 0.57  | 0.95   | 1.38   | 2.72  | 0.78  | 0.24  | 0.79  | -0.03 |
| Denmark        | 1.02  | 1.72   | 1.97   | 2.88  | 0.60  | 0.10  | 0.72  | -0.44 |
| Finland        | 0.67  | 0.88   | 1.62   | 2.40  | 0.81  | 0.00  | 0.68  | -0.11 |
| France         | 0.85  | 1.02   | 1.51   | 2.46  | 0.86  | 0.04  | 0.65  | -0.45 |
| Germany        | 0.61  | 0.69   | 1.05   | 2.09  | 0.88  | 0.36  | 0.67  | 0.29  |
| Iceland        | 0.89  | 2.18   | 3.00   | 10.13 | 0.78  | 0.09  | 0.25  | 0.14  |
| Ireland        | 0.79  | 1.80   | 2.25   | 2.98  | 0.59  | 0.08  | 0.52  | -0.17 |
| Israel         | 0.83  | 1.03   | 1.35   | 2.92  | 0.94  | -0.06 | 0.44  | -0.00 |
| Italy          | 0.84  | 1.14   | 1.82   | 3.29  | 0.87  | 0.16  | 0.47  | -0.25 |
| Japan          | 0.78  | 0.90   | 0.84   | 2.27  | 0.89  | 0.22  | 0.35  | -0.07 |
| Korea          | 0.77  | 2.03   | 1.80   | 9.43  | 0.75  | 0.53  | 0.57  | -0.53 |
| Luxembourg     | 1.08  | 1.57   | 2.17   | 3.00  | 0.95  | -0.61 | -0.23 | 0.52  |
| Mexico         | 0.82  | 1.10   | 1.42   | 4.14  | 0.92  | 0.43  | 0.25  | -0.65 |
| Netherlands    | 0.86  | 1.39   | 2.15   | 2.31  | 0.78  | 0.04  | 0.50  | 0.13  |
| New Zealand    | 0.87  | 1.47   | 1.96   | 4.66  | 0.70  | 0.10  | 0.58  | -0.23 |
| Norway         | 0.96  | 2.84   | 2.60   | 8.12  | 0.76  | -0.25 | 0.24  | 0.17  |
| Spain          | 0.53  | 1.84   | 2.37   | 5.76  | 0.28  | 0.70  | 0.34  | -0.36 |
| Sweden         | 0.76  | 1.04   | 0.92   | 1.56  | 0.80  | -0.15 | 0.84  | 0.10  |
| Switzerland    | 0.97  | 1.12   | 2.88   | 10.98 | 0.89  | -0.37 | 0.37  | 0.21  |
| United Kingdom | 0.66  | 1.35   | 1.08   | 2.79  | 0.64  | 0.62  | 0.45  | -0.41 |
| United States  | 0.52  | 1.12   | 1.37   | 3.14  | 0.72  | 0.69  | 0.74  | -0.61 |

*Series are first logged and detrended using the HP filter.*

**Table 4: Properties of output components**

| Country        | σyA/σy | σyτl/σy | σyτx/σy | σyg/σy | ρyA,y | ρyτl,y | ρyτx,y | ρyg,y |
|----------------|--------|---------|---------|--------|-------|--------|--------|-------|
| Australia      | 0.92   | 0.94    | 0.85    | 0.34   | 0.67  | -0.10  | 0.71   | -0.39 |
| Austria        | 1.06   | 0.98    | 1.05    | 0.34   | 0.82  | -0.32  | 0.37   | 0.16  |
| Belgium        | 0.77   | 1.00    | 0.44    | 0.37   | 0.72  | 0.68   | -0.34  | -0.29 |
| Canada         | 0.67   | 0.42    | 0.63    | 0.54   | 0.89  | -0.03  | 0.79   | -0.15 |
| Denmark        | 1.18   | 0.95    | 0.89    | 0.50   | 0.58  | -0.15  | 0.72   | -0.38 |
| Finland        | 0.74   | 0.72    | 0.89    | 0.23   | 0.80  | -0.33  | 0.71   | 0.05  |
| France         | 1.11   | 0.93    | 0.92    | 0.36   | 0.88  | -0.45  | 0.64   | -0.44 |
| Germany        | 0.74   | 0.34    | 0.61    | 0.24   | 0.87  | 0.02   | 0.69   | -0.24 |
| Iceland        | 0.97   | 1.19    | 1.44    | 0.64   | 0.75  | -0.15  | 0.27   | 0.06  |
| Ireland        | 0.84   | 0.92    | 0.92    | 0.41   | 0.62  | -0.02  | 0.53   | -0.09 |
| Israel         | 0.83   | 0.58    | 0.59    | 0.28   | 0.92  | 0.08   | 0.40   | -0.19 |
| Italy          | 0.99   | 1.03    | 1.39    | 0.37   | 0.85  | -0.32  | 0.51   | -0.59 |
| Japan          | 0.97   | 0.48    | 0.46    | 0.25   | 0.85  | 0.01   | 0.35   | 0.04  |
| Korea          | 1.04   | 0.99    | 0.90    | 0.32   | 0.69  | -0.12  | 0.58   | -0.49 |
| Luxembourg     | 1.14   | 1.01    | 1.14    | 0.49   | 0.95  | -0.18  | -0.20  | 0.71  |
| Mexico         | 0.97   | 0.69    | 0.68    | 0.20   | 0.91  | 0.15   | 0.21   | -0.74 |
| Netherlands    | 0.99   | 0.87    | 1.06    | 0.32   | 0.72  | -0.27  | 0.50   | -0.09 |
| New Zealand    | 1.06   | 0.83    | 0.88    | 0.29   | 0.66  | -0.14  | 0.58   | -0.37 |
| Norway         | 1.08   | 2.15    | 1.35    | 1.10   | 0.71  | -0.21  | 0.24   | 0.31  |
| Spain          | 0.72   | 1.15    | 1.29    | 0.35   | 0.34  | 0.35   | 0.35   | -0.32 |
| Sweden         | 0.93   | 0.53    | 0.40    | 0.26   | 0.93  | -0.28  | 0.84   | -0.21 |
| Switzerland    | 1.13   | 1.15    | 1.32    | 0.48   | 0.90  | -0.25  | 0.35   | -0.41 |
| United Kingdom | 0.73   | 0.85    | 0.55    | 0.25   | 0.61  | 0.50   | 0.43   | -0.38 |
| United States  | 0.60   | 0.58    | 0.61    | 0.37   | 0.76  | 0.64   | 0.74   | -0.66 |

*Series are first logged and detrended using the HP filter.*

**Table 5: Properties of labor components**

| Country        | σlA/σl | σlτl/σl | σlτx/σl | σlg/σl | ρlA,l | ρlτl,l | ρlτx,l | ρlg,l |
|----------------|--------|---------|---------|--------|-------|--------|--------|-------|
| Australia      | 0.27   | 1.20    | 1.08    | 0.43   | 0.39  | 0.42   | 0.50   | -0.33 |
| Austria        | 0.28   | 1.77    | 1.90    | 0.61   | -0.14 | 0.36   | 0.20   | 0.04  |
| Belgium        | 0.26   | 1.40    | 0.61    | 0.51   | 0.36  | 0.95   | -0.50  | -0.29 |
| Canada         | 0.39   | 0.66    | 0.99    | 0.84   | 0.75  | 0.36   | 0.82   | -0.40 |
| Denmark        | 0.23   | 1.10    | 1.03    | 0.58   | -0.44 | 0.73   | 0.53   | -0.44 |
| Finland        | 0.16   | 1.25    | 1.56    | 0.41   | 0.19  | 0.05   | 0.61   | -0.06 |
| France         | 0.63   | 1.90    | 1.87    | 0.74   | 0.25  | 0.20   | 0.38   | -0.34 |
| Germany        | 0.27   | 0.63    | 1.13    | 0.45   | 0.40  | 0.31   | 0.78   | -0.36 |
| Iceland        | 0.22   | 2.05    | 2.47    | 1.10   | -0.33 | 0.29   | 0.37   | -0.32 |
| Ireland        | 0.21   | 1.23    | 1.24    | 0.55   | 0.30  | 0.53   | 0.39   | -0.35 |
| Israel         | 0.09   | 1.69    | 1.74    | 0.81   | -0.88 | 0.38   | 0.33   | -0.13 |
| Italy          | 0.55   | 2.15    | 2.90    | 0.78   | 0.07  | 0.15   | 0.29   | -0.26 |
| Japan          | 0.49   | 1.06    | 1.02    | 0.55   | -0.05 | 0.46   | 0.51   | -0.04 |
| Korea          | 0.45   | 1.48    | 1.35    | 0.47   | -0.28 | 0.49   | 0.34   | -0.33 |
| Luxembourg     | 0.46   | 3.22    | 3.63    | 1.57   | -0.18 | 0.39   | 0.08   | -0.31 |
| Mexico         | 0.38   | 1.64    | 1.62    | 0.47   | 0.17  | 0.39   | 0.29   | -0.40 |
| Netherlands    | 0.39   | 1.45    | 1.76    | 0.53   | -0.35 | 0.39   | 0.41   | -0.28 |
| New Zealand    | 0.28   | 1.16    | 1.23    | 0.40   | -0.43 | 0.47   | 0.55   | -0.25 |
| Norway         | 0.58   | 3.49    | 2.20    | 1.79   | -0.13 | 0.31   | 0.25   | -0.32 |
| Spain          | 0.31   | 1.19    | 1.33    | 0.36   | 0.10  | 0.49   | 0.42   | -0.43 |
| Sweden         | 0.75   | 0.93    | 0.69    | 0.45   | 0.83  | 0.16   | 0.70   | -0.57 |
| Switzerland    | 0.38   | 2.62    | 3.00    | 1.09   | -0.03 | 0.30   | 0.13   | -0.16 |
| United Kingdom | 0.12   | 1.16    | 0.75    | 0.34   | -0.27 | 0.81   | 0.29   | -0.28 |
| United States  | 0.14   | 0.84    | 0.89    | 0.54   | 0.64  | 0.83   | 0.75   | -0.68 |

*Series are first logged and detrended using the HP filter.*

**Table 6: Properties of investment components**

| Country        | σxA/σx | σxτl/σx | σxτx/σx | σxg/σx | ρxA,x | ρxτl,x | ρxτx,x | ρxg,x |
|----------------|--------|---------|---------|--------|-------|--------|--------|-------|
| Australia      | 0.38   | 0.38    | 0.77    | 0.20   | 0.78  | -0.31  | 0.87   | 0.71  |
| Austria        | 0.62   | 0.71    | 1.35    | 0.10   | 0.53  | -0.71  | 0.89   | -0.31 |
| Belgium        | 0.39   | 0.76    | 0.47    | 0.27   | 0.84  | 0.91   | -0.69  | 0.83  |
| Canada         | 0.43   | 0.16    | 0.75    | 0.24   | 0.89  | -0.28  | 0.97   | -0.26 |
| Denmark        | 0.54   | 0.42    | 0.86    | 0.07   | 0.44  | -0.32  | 0.97   | 0.84  |
| Finland        | 0.34   | 0.39    | 0.95    | 0.10   | 0.73  | -0.66  | 0.98   | 0.81  |
| France         | 0.63   | 0.58    | 0.97    | 0.06   | 0.90  | -0.72  | 0.91   | -0.55 |
| Germany        | 0.53   | 0.22    | 0.93    | 0.22   | 0.58  | -0.12  | 0.96   | -0.68 |
| Iceland        | 0.29   | 0.40    | 1.12    | 0.32   | -0.17 | -0.36  | 0.93   | 0.64  |
| Ireland        | 0.34   | 0.40    | 0.92    | 0.08   | 0.49  | -0.36  | 0.95   | 0.92  |
| Israel         | 0.39   | 0.33    | 0.79    | 0.25   | 0.69  | -0.03  | 0.83   | 0.35  |
| Italy          | 0.47   | 0.48    | 1.37    | 0.24   | 0.54  | -0.73  | 0.90   | -0.68 |
| Japan          | 0.70   | 0.37    | 0.74    | 0.14   | 0.65  | -0.01  | 0.71   | 0.13  |
| Korea          | 0.56   | 0.50    | 1.01    | 0.13   | 0.57  | -0.59  | 0.93   | -0.01 |
| Luxembourg     | 0.58   | 0.58    | 1.33    | 0.22   | 0.23  | -0.92  | 0.87   | 0.90  |
| Mexico         | 0.50   | 0.36    | 0.92    | 0.09   | 0.67  | -0.12  | 0.72   | 0.30  |
| Netherlands    | 0.60   | 0.54    | 1.30    | 0.12   | 0.20  | -0.70  | 0.96   | -0.20 |
| New Zealand    | 0.51   | 0.40    | 0.96    | 0.20   | 0.36  | -0.47  | 0.94   | 0.52  |
| Norway         | 0.48   | 0.88    | 1.05    | 0.57   | -0.06 | 0.20   | 0.44   | 0.64  |
| Spain          | 0.38   | 0.49    | 1.24    | 0.11   | 0.13  | -0.36  | 0.90   | 0.21  |
| Sweden         | 0.74   | 0.32    | 0.51    | 0.20   | 0.94  | -0.36  | 0.97   | -0.43 |
| Switzerland    | 0.35   | 0.41    | 1.10    | 0.37   | 0.27  | -0.81  | 0.99   | 0.36  |
| United Kingdom | 0.39   | 0.50    | 0.73    | 0.16   | 0.42  | 0.23   | 0.84   | 0.83  |
| United States  | 0.35   | 0.29    | 0.92    | 0.12   | 0.79  | 0.15   | 0.94   | -0.86 |

*Series are first logged and detrended using the HP filter.*

**Table 7: Properties of consumption components**

| Country        | σcA/σc | σcτl/σc | σcτx/σc | σcg/σc | ρcA,c | ρcτl,c | ρcτx,c | ρcg,c |
|----------------|--------|---------|---------|--------|-------|--------|--------|-------|
| Australia      | 0.63   | 0.62    | 0.92    | 0.36   | 0.09  | 0.65   | 0.62   | 0.30  |
| Austria        | 1.07   | 0.73    | 0.46    | 0.75   | 0.53  | 0.29   | 0.18   | 0.12  |
| Belgium        | 0.56   | 0.98    | 1.02    | 0.49   | 0.09  | 0.64   | -0.26  | 0.56  |
| Canada         | 0.98   | 0.93    | 0.85    | 0.79   | -0.14 | 0.78   | -0.05  | 0.61  |
| Denmark        | 0.93   | 0.88    | 1.03    | 1.02   | 0.37  | 0.57   | -0.44  | 0.53  |
| Finland        | 1.32   | 0.93    | 1.11    | 1.21   | 0.12  | 0.42   | 0.14   | 0.07  |
| France         | 1.00   | 0.54    | 0.58    | 0.61   | 0.50  | 0.17   | 0.05   | 0.38  |
| Germany        | 1.07   | 0.52    | 0.42    | 0.71   | 0.62  | 0.44   | -0.18  | 0.21  |
| Iceland        | 0.51   | 0.44    | 0.86    | 0.89   | 0.13  | 0.47   | 0.35   | 0.40  |
| Ireland        | 0.57   | 0.47    | 0.97    | 0.99   | 0.06  | 0.75   | 0.28   | 0.05  |
| Israel         | 0.86   | 0.59    | 0.50    | 0.46   | 0.46  | 0.49   | 0.33   | 0.28  |
| Italy          | 0.85   | 0.82    | 0.57    | 0.61   | 0.49  | 0.28   | 0.49   | 0.24  |
| Japan          | 0.93   | 0.43    | 0.36    | 0.49   | 0.63  | 0.39   | 0.20   | 0.33  |
| Korea          | 0.40   | 0.38    | 0.57    | 1.43   | 0.41  | 0.19   | -0.61  | 0.85  |
| Luxembourg     | 1.46   | 1.07    | 0.88    | 2.45   | 0.02  | 0.37   | -0.07  | 0.22  |
| Mexico         | 0.51   | 0.38    | 0.28    | 0.32   | 0.72  | 0.37   | 0.27   | 0.67  |
| Netherlands    | 1.02   | 0.87    | 0.57    | 0.81   | 0.25  | 0.61   | 0.14   | 0.08  |
| New Zealand    | 0.66   | 0.56    | 0.72    | 0.57   | 0.32  | 0.58   | 0.28   | 0.32  |
| Norway         | 0.43   | 0.54    | 0.24    | 0.79   | -0.02 | 0.80   | 0.45   | 0.25  |
| Spain          | 0.40   | 0.62    | 0.89    | 0.86   | 0.26  | 0.79   | 0.17   | 0.22  |
| Sweden         | 0.92   | 0.46    | 0.41    | 0.41   | 0.51  | 0.71   | 0.36   | 0.04  |
| Switzerland    | 1.51   | 1.14    | 1.41    | 2.65   | 0.39  | 0.23   | 0.05   | 0.03  |
| United Kingdom | 0.72   | 0.80    | 0.42    | 0.30   | 0.52  | 0.67   | 0.07   | 0.35  |
| United States  | 0.67   | 0.85    | 0.74    | 0.25   | 0.45  | 0.81   | -0.28  | 0.60  |

*Series are first logged and detrended using the HP filter.*

## 7. United States — MLE Estimates

*MLE estimates reported for the period 1980Q1 to 2014Q4.*

**Table 8: P matrix**

|         |         |         |         |
|---------|---------|---------|---------|
|  0.9887 |  0.0307 | -0.0089 | -0.0407 |
| -0.0012 |  1.0011 | -0.0275 |  0.0175 |
| -0.0045 |  0.0449 |  0.9675 | -0.0426 |
|  0.0063 |  0.0017 |  0.0016 |  0.9945 |

> ⚠️ **Convention note (read before transcribing).** Table 8 above is
> printed in the paper's **"rows = drivers, columns = receivers"**
> narrative convention: row 0 reads "what z does" — its self-persistence
> in column 0, its outgoing spillover to τ_l in column 1, etc. This is
> the **transpose** of the textbook VAR convention
> `state_{t+1} = P · state_t` that BCKM's matlab code (`mleqadj.m:222`)
> and our codebase actually use. In the code convention, row 0 reads
> "what determines z's update" — coefficients on lags of [z, τ_l, τ_x, g].
> The canonical CODE-convention matrix lives in
> [`bca_core/constants.py`](bca_core/constants.py) as
> `P_BCKM_TABLE8`; **import from there, never re-transcribe Table 8**.
> A transposed-P bug across nine independent transcription sites cost us
> +501 nats of LL at BCKM-published θ; see Diary 2026-04-30.

**Table 9: P₀ vector**

|        |        |        |         |
|--------|--------|--------|---------|
| 0.0140 | 0.0008 | 0.0129 | -0.0137 |

**Table 10: Q matrix** (lower-triangular factor of V = QQ′)

|         |        |         |        |
|---------|--------|---------|--------|
|  0.0077 | 0.0024 | -0.0041 | 0.0003 |
|  0.0024 | 0.0043 |  0.0023 | 0.0153 |
| -0.0041 | 0.0023 |  0.0088 | 0.0121 |
|  0.0003 | 0.0153 |  0.0121 | 0.0139 |

## 8. Great Recession — United States

### 8.1 Summary Statistics

**Table 11: f statistics for observable's components (2008), United States**

| Component                             | United States |
|---------------------------------------|---------------|
| fYA (Output, efficiency wedge)        | 0.16          |
| fYτL (Output, labor wedge)            | 0.46          |
| fYτx (Output, investment wedge)       | 0.32          |
| fLA (Labor, efficiency wedge)         | 0.04          |
| fLτL (Labor, labor wedge)             | 0.70          |
| fLτx (Labor, investment wedge)        | 0.25          |
| fxA (Investment, efficiency wedge)    | 0.05          |
| fxτL (Investment, labor wedge)        | 0.05          |
| fxτx (Investment, investment wedge)   | 0.88          |
| fCA (Consumption, efficiency wedge)   | 0.07          |
| fCτL (Consumption, labor wedge)       | 0.82          |
| fCτx (Consumption, investment wedge)  | 0.03          |

*For reference — Avg all countries: fYA=0.64, fYτL=0.09, fYτx=0.20; fLA=0.33, fLτL=0.19, fLτx=0.31; fxA=0.36, fxτL=0.07, fxτx=0.48; fCA=0.47, fCτL=0.20, fCτx=0.12.*

**Table 12: Peak-trough declines — United States**

|                    | Trough   | Total  | Eff. (A) | Labor (L) | Invest. (X) | Govt. (G) |
|--------------------|----------|--------|----------|-----------|-------------|-----------|
| Output (%)         | 2009.625 | -7.0   | -1.9     | -3.4      | -4.5        | 2.7       |
| Labor (hours) (%)  | 2009.625 | -7.5   | -0.9     | -5.0      | -6.7        | 4.1       |
| Investment (%)     | 2009.625 | -23.2  | -4.9     | -3.0      | -21.6       | 3.2       |

### 8.2 Properties of the Output Components

**Table 13A: Summary Statistics**

| Output Components       | S/Y  | ρ(y_ω,y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|-------------------------|------|-----------------|--------|-------|-------|-------|
| Efficiency              | 0.60 | 0.57            | 0.66   | 0.76  | 0.48  | 0.25  |
| Labor                   | 0.58 | 0.44            | 0.56   | 0.64  | 0.69  | 0.66  |
| Investment              | 0.61 | 0.53            | 0.65   | 0.74  | 0.63  | 0.45  |
| Government Consumption  | 0.37 | -0.55           | -0.63  | -0.66 | -0.56 | -0.38 |

**Table 13B: Cross Correlations**

| Output components (X, Y) | ρ(X,Y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|--------------------------|---------------|--------|-------|-------|-------|
| Efficiency, labor        | 0.34          | 0.29   | 0.12  | 0.06  | -0.03 |
| Efficiency, investment   | 0.49          | 0.55   | 0.54  | 0.30  | 0.11  |
| Efficiency, gc           | -0.45         | -0.51  | -0.46 | -0.32 | -0.15 |
| Labor, investment        | 0.09          | 0.21   | 0.30  | 0.45  | 0.54  |
| Labor, gc                | -0.04         | -0.25  | -0.46 | -0.52 | -0.57 |
| Investment, gc           | -0.55         | -0.70  | -0.82 | -0.75 | -0.62 |

*Series are first logged and detrended using the HP filter.*

### 8.3 Properties of the Labor Components

**Table 14A: Summary Statistics**

| Labor Components        | S/L  | ρ(h_ω,h) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|-------------------------|------|-----------------|--------|-------|-------|-------|
| Efficiency              | 0.14 | 0.70            | 0.70   | 0.64  | 0.44  | 0.23  |
| Labor                   | 0.84 | 0.55            | 0.70   | 0.83  | 0.81  | 0.74  |
| Investment              | 0.89 | 0.71            | 0.77   | 0.75  | 0.60  | 0.38  |
| Government Consumption  | 0.54 | -0.70           | -0.72  | -0.68 | -0.51 | -0.30 |

**Table 14B: Cross Correlations**

| Labor components (X, Y)  | ρ(X,Y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|--------------------------|---------------|--------|-------|-------|-------|
| Efficiency, labor        | 0.61          | 0.59   | 0.49  | 0.30  | 0.10  |
| Efficiency, investment   | 0.62          | 0.73   | 0.76  | 0.59  | 0.41  |
| Efficiency, gc           | -0.63         | -0.79  | -0.89 | -0.70 | -0.50 |
| Labor, investment        | 0.09          | 0.21   | 0.30  | 0.45  | 0.54  |
| Labor, gc                | -0.04         | -0.25  | -0.46 | -0.52 | -0.57 |
| Investment, gc           | -0.55         | -0.70  | -0.82 | -0.75 | -0.62 |

*Series are first logged and detrended using the HP filter.*

### 8.4 Properties of the Investment Components

**Table 15A: Summary Statistics**

| Investment Components   | S/X  | ρ(x_ω,x) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|-------------------------|------|-----------------|--------|-------|-------|-------|
| Efficiency              | 0.35 | 0.63            | 0.73   | 0.79  | 0.53  | 0.30  |
| Labor                   | 0.29 | 0.01            | 0.09   | 0.15  | 0.31  | 0.42  |
| Investment              | 0.92 | 0.60            | 0.79   | 0.94  | 0.81  | 0.62  |
| Government Consumption  | 0.12 | -0.60           | -0.76  | -0.86 | -0.77 | -0.63 |

**Table 15B: Cross Correlations**

| Investment components (X, Y) | ρ(X,Y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|------------------------------|---------------|--------|-------|-------|-------|
| Efficiency, labor            | 0.43          | 0.39   | 0.23  | 0.13  | 0.01  |
| Efficiency, investment       | 0.55          | 0.62   | 0.63  | 0.40  | 0.20  |
| Efficiency, gc               | -0.52         | -0.62  | -0.61 | -0.45 | -0.26 |
| Labor, investment            | -0.17         | -0.13  | -0.10 | 0.11  | 0.28  |
| Labor, gc                    | 0.22          | 0.07   | -0.08 | -0.19 | -0.31 |
| Investment, gc               | -0.55         | -0.70  | -0.82 | -0.75 | -0.62 |

*Series are first logged and detrended using the HP filter.*

### 8.5 Properties of Model Data

**Table 16A: Summary Statistics**

| Observables (X / y)     | σ(X)/σ(y) | ρ(X,y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|-------------------------|-----------|---------------|--------|-------|-------|-------|
| Output                  | 1.00      | 0.67          | 0.85   | 1.00  | 0.85  | 0.67  |
| Hours                   | 1.03      | 0.60          | 0.76   | 0.88  | 0.85  | 0.73  |
| Investment              | 3.47      | 0.66          | 0.81   | 0.92  | 0.76  | 0.56  |
| Government Consumption  | 3.14      | -0.52         | -0.59  | -0.61 | -0.51 | -0.33 |

**Table 16B: Cross Correlations**

| Observables (X, Y)    | ρ(X,Y) lag -2 | lag -1 | lag 0 | lag 1 | lag 2 |
|-----------------------|---------------|--------|-------|-------|-------|
| Output, hours         | 0.73          | 0.85   | 0.88  | 0.76  | 0.60  |
| Output, investment    | 0.56          | 0.76   | 0.92  | 0.81  | 0.66  |
| Output, gc            | -0.33         | -0.51  | -0.61 | -0.59 | -0.52 |
| Hours, investment     | 0.47          | 0.68   | 0.85  | 0.85  | 0.78  |
| Hours, gc             | -0.25         | -0.46  | -0.64 | -0.68 | -0.68 |
| Investment, gc        | -0.52         | -0.69  | -0.80 | -0.73 | -0.60 |

*Series are first logged and detrended using the HP filter.*