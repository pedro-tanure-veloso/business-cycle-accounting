# Accounting for Business Cycles

*Based on Brinca, Chari, Kehoe, and McGrattan — NBER Working Paper No. 22663 (September 2016) and the accompanying Technical Appendix.*

## Contents

- **Part I — Methodology** (from the main paper)
  - 1. The Benchmark Prototype Economy
  - 2. The Accounting Procedure
  - 3. Details of the Application
- **Part II — Data and Calibration** (from the Technical Appendix)
  - 4. Data and Sources
  - 5. Parametrization and Calibration
- **Part III — Findings: Cross-Country Application** (from the main paper §3.B)
  - The Great Recession
    - Three Illustrative Recessions (United States detail kept; France and Ireland omitted)
    - Summary Statistics for OECD Countries
    - Comparing the Great Recession with the 1980s Recessions
    - Summary Statistics for the Entire Period
    - The Importance of the Classification of Consumer Durables
    - Comparing our Procedure with a Perfect Foresight Procedure
- **Part IV — Empirical Results: United States Replication** (from the Technical Appendix)
  - 6. Whole-Sample Tables
  - 7. United States — MLE Estimates
  - 8. Great Recession — United States

---

# Part I — Methodology

## 1. The Benchmark Prototype Economy

The benchmark prototype economy used later in our accounting procedure is a stochastic growth model. In each period *t*, the economy experiences one of finitely many events *s_t*, which index the shocks. We denote by *s^t* = (*s₀, ..., s_t*) the history of events up through and including period *t* and often refer to *s^t* as the state. The probability, as of period 0, of any particular history *s^t* is π_t(*s^t*). The initial realization *s₀* is given. The economy has four exogenous stochastic variables, all of which are functions of the underlying random variable *s^t*: the efficiency wedge A_t(*s^t*), the labor wedge 1 − τ_lt(*s^t*), the investment wedge 1/[1 + τ_xt(*s^t*)], and the government consumption wedge g_t(*s^t*).

In the model, consumers maximize expected utility over per capita consumption c_t and per capita labor l_t:

$$\sum_{t=0}^{\infty} \sum_{s^t} \beta^t \pi_t(s^t) U(c_t(s^t), l_t(s^t)) N_t,$$

subject to the budget constraint

$$c_t + [1 + \tau_{xt}(s^t)]\, x_t(s^t) = [1 - \tau_{lt}(s^t)]\, w_t(s^t)\, l_t(s^t) + r_t(s^t)\, k_t(s^{t-1}) + T_t(s^t)$$

and the capital accumulation law

$$(1)\quad (1 + n)\, k_{t+1}(s^t) = (1 - \delta)\, k_t(s^{t-1}) + x_t(s^t),$$

where k_t(*s^{t-1}*) denotes the per capita capital stock, x_t(*s^t*) per capita investment, w_t(*s^t*) the wage rate, r_t(*s^t*) the rental rate on capital, β the discount factor, δ the depreciation rate of capital, N_t the population with growth rate equal to 1 + n, and T_t(*s^t*) per capita lump-sum transfers.

The production function is A(*s^t*) F(k_t(*s^{t-1}*), (1 + γ)^t l_t(*s^t*)), where 1 + γ is the rate of labor-augmenting technical progress, which is assumed to be a constant. Firms maximize profits given by A_t(*s^t*) F(k_t(*s^{t-1}*), (1 + γ)^t l_t(*s^t*)) − r_t(*s^t*) k_t(*s^{t-1}*) − w_t(*s^t*) l_t(*s^t*).

The equilibrium of this benchmark prototype economy is summarized by the resource constraint,

$$(2)\quad c_t(s^t) + x_t(s^t) + g_t(s^t) = y_t(s^t),$$

where y_t(*s^t*) denotes per capita output, together with

$$(3)\quad y_t(s^t) = A_t(s^t)\, F(k_t(s^{t-1}), (1 + \gamma)^t\, l_t(s^t)),$$

$$(4)\quad \frac{-U_{lt}(s^t)}{U_{ct}(s^t)} = [1 - \tau_{lt}(s^t)]\, A_t(s^t)\, (1 + \gamma)^t\, F_{lt},\text{ and}$$

$$(5)\quad U_{ct}(s^t)[1 + \tau_{xt}(s^t)] = \beta \sum_{s^{t+1}} \pi_t(s^{t+1}|s^t)\, U_{ct+1}(s^{t+1}) \{ A_{t+1}(s^{t+1})\, F_{kt+1}(s^{t+1}) + (1 - \delta)[1 + \tau_{xt+1}(s^{t+1})] \},$$

where notations such as U_ct, U_lt, F_lt, and F_kt denote the derivatives of the utility function and the production function with respect to their arguments, and π_t(*s^{t+1}*|*s^t*) denotes the conditional probability π_t(*s^{t+1}*)/π_t(*s^t*). We assume that g_t(*s^t*) fluctuates around a trend of (1 + γ)^t.

Notice that in this benchmark prototype economy, the efficiency wedge resembles a blueprint technology parameter, and the labor and investment wedges resemble tax rates on labor income and investment. Other more elaborate models could be considered, such as models with other kinds of frictions that look like taxes on consumption or capital income. Consumption taxes induce a wedge between the consumption–leisure marginal rate of substitution and the marginal product of labor in the same way as do labor income taxes. Such taxes, if time-varying, also distort the intertemporal margins in (5). Capital income taxes induce a wedge between the intertemporal marginal rate of substitution and the marginal product of capital, which is only slightly different from the distortion induced by a tax on investment. We experimented with intertemporal distortions that resemble capital income taxes rather than investment taxes and found that our substantive conclusions are unaffected.

We emphasize that each of the wedges represents the overall distortion to the relevant equilibrium condition of the model. For example, distortions to labor supply affecting consumers and to labor demand affecting firms both distort the static first-order condition (4). Our labor wedge represents the sum of these distortions. Thus, our method identifies the overall wedge induced by both distortions and does not identify each separately. Likewise, liquidity constraints on consumers distort the consumer's intertemporal Euler equation, whereas investment financing frictions on firms distort the firm's intertemporal Euler equation. Our method combines the Euler equations for the consumer and the firm and therefore identifies only the overall wedge in the combined Euler equation given by (5). We focus on the overall wedges because what matters in determining business cycle fluctuations is the overall wedges, not each distortion separately.

For the equivalence results that follow, it is notationally convenient to work with the prototype model just described. For our quantitative results, we add investment adjustment costs by replacing the capital accumulation law (1) with

$$(6)\quad (1 + n)\, k_{t+1}(s^t) = (1 - \delta)\, k_t(s^{t-1}) + x_t(s^t) - \phi\!\left(\frac{x_t(s^t)}{k_t(s^{t-1})}\right),$$

where φ represents the per unit cost of adjusting the capital stock. We follow the macroeconomic literature in assuming that the adjustment costs are parameterized by the function

$$\phi\!\left(\frac{x}{k}\right) = \frac{a}{2}\left(\frac{x}{k} - b\right)^2,$$

where *b* = γ + δ + n is the steady-state value of the investment–capital ratio.

## 2. The Accounting Procedure

Having established our equivalence result, we now describe our accounting procedure at a conceptual level, discuss a Markovian implementation of it, and distinguish our procedure from others.

Our procedure is designed to answer questions of the following kind: How much would output fluctuate if the only wedge that fluctuated is the efficiency wedge and the probability distribution of the efficiency wedge is the same as in the prototype economy? Critically, our procedure ensures that agents' expectations of how the efficiency wedge will evolve are the same as in the prototype economy. For each experiment, we compare the properties of the resulting equilibria to those of the prototype economy. These comparisons, together with our equivalence results, allow us to identify promising classes of detailed economies.

### 2.A. The Accounting Procedure at a Conceptual Level

Recall that the state *s^t* is the history of the underlying abstract events *s_t*. Suppose for now that the stochastic process π_t(*s^t*) and the realizations of the state *s^t* in some particular episode are known. Recall that the prototype economy has one underlying (vector-valued) random variable, the state *s^t*, which has a probability of π_t(*s^t*). All of the other stochastic variables, including the four wedges — the efficiency wedge A_t(*s^t*), the labor wedge 1 − τ_lt(*s^t*), the investment wedge 1/[1 + τ_xt(*s^t*)], and the government consumption wedge g_t(*s^t*) — are simply functions of this random variable. Hence, when the state *s^t* is known, so are the wedges.

To evaluate the effects of just the efficiency wedge, for example, we consider an economy, referred to as an *efficiency wedge alone economy*, with the same underlying state *s^t* and probability π_t(*s^t*) and the same function A_t(*s^t*) for the efficiency wedge as in the prototype economy, but in which the other three wedges are set to be constant functions of the state, in that τ_lt(*s^t*) = τ̄_l, τ_xt(*s^t*) = τ̄_x, and g_t(*s^t*) = ḡ. Note that this construction ensures that the probability distribution of the efficiency wedge in this economy is identical to that in the prototype economy.

We compute the decision rules for the efficiency wedge alone economy, denoted y^e(*s^t*), l^e(*s^t*), and x^e(*s^t*). For a given initial value k₀, for any given sequence *s^t*, we refer to the resulting values of output, labor, and investment as the *efficiency wedge components* of output, labor, and investment.

In a similar manner, we define the labor wedge alone economy, the investment wedge alone economy, and the government consumption wedge alone economy, as well as economies with a combination of wedges, such as the efficiency and labor wedge economy.

### 2.B. A Markovian Implementation

So far we have described our procedure assuming that we know the stochastic process π_t(*s^t*) and that we can observe the state *s^t*. In practice, of course, we need to either specify the stochastic process a priori or use data to estimate it, and we need to uncover the state *s^t* from the data. Here we describe a set of assumptions that makes these efforts easy. Then we describe in detail the three steps involved in implementing our procedure.

We assume that the state *s_t* follows a Markov process π(*s_t*|*s_{t−1}*) and that the wedges in period *t* can be used to uniquely uncover the event *s_t*, in the sense that the mapping from the event *s_t* to the wedges (A_t, τ_lt, τ_xt, g_t) is one to one and onto. Given this assumption, without loss of generality, let the underlying event *s_t* = (s_{At}, s_{lt}, s_{xt}, s_{gt}), and let A_t(*s^t*) = s_{At}, τ_lt(*s^t*) = s_{lt}, τ_xt(*s^t*) = s_{xt}, and g_t(*s^t*) = s_{gt}. Note that we have effectively assumed that agents use only past wedges to forecast future wedges and that the wedges in period *t* are sufficient statistics for the event in period *t*. This assumption is only to make our estimation easier, and it can be relaxed.

In practice, to estimate the stochastic process for the state, we first specify a vector autoregressive AR(1) process for the event *s_t* = (s_{At}, s_{lt}, s_{xt}, s_{gt}) of the form

$$(60)\quad s_{t+1} = P_0 + P s_t + \varepsilon_{t+1},$$

where the shock ε_t is i.i.d. over time and is distributed normally with mean zero and covariance matrix V. To ensure that our estimate of V is positive semidefinite, we estimate the lower triangular matrix Q, where V = QQ′. The matrix Q has no structural interpretation.

The **first step** in our procedure is to use data on y_t, l_t, x_t, and g_t from an actual economy to estimate the parameters of the Markov process π(*s_t*|*s_{t−1}*). We can do so using a variety of methods, including the maximum likelihood procedure described later.

The **second step** in our procedure is to uncover the event *s_t* by measuring the realized wedges. We measure the government consumption wedge directly from the data as the sum of government consumption and net exports. To obtain the values of the other three wedges, we use the data and the model's decision rules. With y^d_t, l^d_t, x^d_t, g^d_t, and k^d_0 denoting the data and y(*s_t*, k_t), l(*s_t*, k_t), and x(*s_t*, k_t) denoting the decision rules of the model, the realized wedge series *s^d_t* solves

$$(61)\quad y^d_t = y(s^d_t, k_t),\quad l^d_t = l(s^d_t, k_t),\quad \text{and}\quad x^d_t = x(s^d_t, k_t),$$

with k_{t+1} = (1 − δ) k_t + x^d_t, k_0 = k^d_0, and g_t = g^d_t. In effect, we solve for the three unknown elements of the vector *s_t* using the three equations (3)–(5) and thereby uncover the state. The four wedges account for all of the movement in output, labor, investment, and government consumption.

Note also that, in measuring the realized wedges, the estimated stochastic process plays a role only in measuring the investment wedge. The efficiency and labor wedges can equivalently be calculated directly from (3) and (4) without computing the equilibrium of the model. Calculating the investment wedge requires computing the equilibrium because the right side of (5) has expectations over future values of consumption, the capital stock, the wedges, and so on.

The **third step** in our procedure is to conduct experiments to isolate the marginal effects of the wedges. To do that, we allow a subset of the wedges to fluctuate as they do in the data while the others are set to constants. To evaluate the effects of the efficiency wedge, we compute the decision rules for the efficiency wedge alone economy, denoted y^e(*s_t*, k_t), l^e(*s_t*, k_t), and x^e(*s_t*, k_t), in which A_t(*s^t*) = s_{At}, τ_lt(*s^t*) = τ̄_l, τ_xt(*s^t*) = τ̄_x, and g_t(*s^t*) = ḡ. Starting from k^d_0, we use *s^d_t*, the decision rules, and the capital accumulation law to compute y^e_t, l^e_t, and x^e_t — the *efficiency wedge components* of output, labor, and investment — and compare them to the data. Other components are computed and compared similarly.

By distinguishing the events to which the wedges are indexed from the wedges themselves, we can separate the *direct effect* and the *forecasting effect* of fluctuations in wedges. When we hold a particular wedge constant, we eliminate the direct effect of that wedge but retain the forecasting effect of the underlying state on the future evolution of the wedge. By doing so, we ensure that expectations of the fluctuating wedges are identical to those in the prototype economy.

### 2.C. Distinguishing Our Procedure from Others

Since this way of separating the direct and forecasting effects of wedges is critical to our procedure, here we describe an alternative procedure that might, at first, seem like the intuitive way to proceed but does not answer the question that interests us.

Consider a simple example with just two wedges, an efficiency wedge and a labor wedge, denoted W_t = (A_t, τ_lt)′. Suppose that we used our prototype model to estimate the following vector process for them of the form W_{t+1} = P W_t + ε_{t+1} where E ε_t ε′_t = V:

$$(62)\quad \begin{bmatrix} A_{t+1} \\ \tau_{lt+1} \end{bmatrix} = \begin{bmatrix} P_{AA} & P_{Al} \\ P_{lA} & P_{ll} \end{bmatrix} \begin{bmatrix} A_t \\ \tau_{lt} \end{bmatrix} + \begin{bmatrix} \varepsilon_{At+1} \\ \varepsilon_{lt+1} \end{bmatrix},$$

where we have suppressed the constant terms. Suppose also that we have decision rules of the form

$$(63)\quad y_t = y(W_t, k_t),\quad l_t = l(W_t, k_t),\quad \text{and}\quad x_t = x(W_t, k_t)$$

and that we have recovered the realized wedge series W^d_t along with the realized innovation series ε^d_{t+1}.

Now suppose we want to answer the question: How much would output fluctuate under the following three conditions? First, only the efficiency wedge fluctuates. Second, the realized sequence of the efficiency wedges coincides with that in the data. Third, the probability distribution of the efficiency wedge is the same as in the prototype economy.

A first attempt is to simply feed a realized innovation series ε̂_{t+1} = (ε^d_{At+1}, 0) for the event and to simulate the resulting shocks using

$$(64)\quad \begin{bmatrix} \hat{A}_{t+1} \\ \hat{\tau}_{lt+1} \end{bmatrix} = \begin{bmatrix} P_{AA} & P_{Al} \\ P_{lA} & P_{ll} \end{bmatrix} \begin{bmatrix} \hat{A}_t \\ \hat{\tau}_{lt} \end{bmatrix} + \begin{bmatrix} \varepsilon^d_{At+1} \\ 0 \end{bmatrix}.$$

This attempt meets our first condition but not the second if P or V has nonzero off-diagonal elements, as we show they do in the data. It also does not meet the third condition.

For a second attempt, we choose the sequence of innovations so that the first two conditions are met — that is, we choose {ε̂_{t+1}} so that (Â_t, τ̂_lt) = (A^d_t, τ̄_l) in the event. The problem is that agents' forecasts about future efficiency wedges are different under this procedure from what they are in the prototype economy. The expected value of A_{t+1} under this procedure is

$$(65)\quad E_t A_{t+1} = P_{AA} A^d_t + P_{Al} \bar{\tau}_l \quad \text{and} \quad E_t \tau_{lt+1} = P_{lA} A^d_t + P_{ll} \bar{\tau}_l.$$

The expectation of the underlying state *s_{t+1}* in the prototype economy is calculated from

$$(66)\quad E_t \begin{bmatrix} s_{At+1} \\ s_{lt+1} \end{bmatrix} = \begin{bmatrix} P_{AA} & P_{Al} \\ P_{lA} & P_{ll} \end{bmatrix} \begin{bmatrix} s^d_{At} \\ s^d_{lt} \end{bmatrix}$$

to be

$$(67)\quad E_t s_{At+1} = P_{AA} s^d_{At} + P_{Al} s^d_{lt} \quad \text{and} \quad E_t s_{lt+1} = P_{lA} s^d_{At} + P_{ll} s^d_{lt}.$$

Clearly, (65) and (67) do not agree when P_{Al} is not zero, so the second attempt does not meet our third condition.

In contrast, our procedure meets all three conditions. In the efficiency wedge alone economy only the efficiency wedge fluctuates, its realized value coincides with the data, and from (60) the probability distribution over *s_{t+1}* (and therefore over A_{t+1}) is the same in both the prototype economy and the efficiency wedge alone economy.

## 3. Details of the Application

To apply our accounting procedure, we use functional forms and parameter values that are familiar from the business cycle literature. We assume that the production function has the form

$$F(k, l) = k^{\alpha} l^{1-\alpha}$$

and the utility function the form

$$U(c, l) = \log c + \psi \log(1 - l).$$

We choose the capital share α to be one-third and the time allocation parameter ψ = 2.5. We choose the depreciation rate δ, the discount factor β, and growth rates γ and n so that, on an annualized basis, depreciation is 5%, the rate of time preference 2.5%, and the population growth rate and the growth of technology are country-specific and computed using OECD data. The adjustment cost parameter *b* = γ + δ + n is pinned down by the previous parameters and varies across countries. For the adjustment cost parameter *a*, we follow Bernanke, Gertler, and Gilchrist (1999) in choosing this parameter so that the elasticity, η, of the price of capital with respect to the investment–capital ratio is 0.25. In this setup, the price of capital q = 1/(1 − φ′), so that, evaluated at the steady state, η = ab. Given η and b, we then set *a* accordingly.

Our prototype economy is a closed economy. When confronting the data, we let government consumption in the model correspond to the sum of government consumption and net exports in the data. The rationale for this choice is given in Chari, Kehoe, and McGrattan (2005), where we prove an equivalence result between an open economy model and a closed economy model in which government consumption is treated in this fashion. We then use a standard maximum likelihood procedure to estimate the parameters P_0, P, and V of the vector AR(1) process for the wedges. In doing so, we use the log-linear decision rules of the prototype economy and data on output, labor, investment, and the sum of government consumption and net exports.

In confronting the theory with the data, we need to decide how to treat consumer durables and sales taxes. At a conceptual level, we think of current expenditures on consumer durables as augmenting the stock of consumer durables, which in turn provides a service flow of consumption to consumers. Based on this idea, we reallocate current expenditures of consumer durables from consumption to investment. We then add the imputed service flow from the stock of consumer durables to consumption and output. This imputed service flow is the rental rate on capital times the stock of durables. We assume that the stock of consumer durables depreciates at the same rate as the stock of physical capital. We also adjust the data to account for sales taxes. We assume that sales taxes are levied solely on consumption — we therefore subtract sales tax revenues from both consumption and measured output.

At a practical level, while the U.S. NIPA accounts have quarterly data on consumer durable expenditures for the 1980:1–2014:4 sample we use, the OECD has more limited data. For countries for which we only have annual data, we fill in quarterly estimates using maximum likelihood estimates of a state space model. For countries for which we only have quarterly data for a subsample, we regress consumer durables on investment and output and use the coefficients to construct estimates of the missing data. Once we have the quarterly series on consumer durables, we construct estimates of the capital stock using the perpetual inventory method. The service flow of durables is assumed to be 4% of the stock of durables.

We express all variables in per capita form and deflate by the GDP deflator. We then estimate separate sets of parameters for the stochastic process for wedges (60) for each of the OECD countries after removing country-specific trends in output, investment, and government consumption. The other parameters are the same across countries. The stochastic process parameters for the Great Recession are estimated using quarterly data for 1980:1–2014:4. The stochastic process (60) with these values is used by agents in our economy to form their expectations about future wedges. Details of the estimated values of the stochastic processes for each country are given in Part III below.

---

# Part II — Data and Calibration

*The following sections provide additional detail on data sources and country-specific parameters, drawn from the Technical Appendix.*

## 4. Data and Sources

The data used for the business cycle accounting exercises come mainly from OECD (variable codes in parenthesis). The time span is from 1980 to the end of 2014 and, unless mentioned otherwise, at the quarterly frequency. For some countries (such as Germany, Ireland, Israel and Mexico), data for most series were only available starting later than 1980Q1 and thus the business cycle accounting exercises were performed for shorter samples.

- **Economic Outlook 98**
  - Gross domestic product, value, market prices (GDP)
  - GDP deflator, market prices (PGDP)
  - Gross capital formation, current prices (ITISK)
  - Government final consumption expenditures, value, expenditure approach (CG)
  - Exports of goods and services, value, national accounts basis (XGS)
  - Imports of goods and services, value, national accounts basis (MGS)
  - Hours worked per employee, total economy (HRS)
  - Total employment (ET)

- **System of Quarterly National Accounts**
  - Durable goods (sub-category of CQRSA: private final consumption expenditure by durability, national currency, current prices)

- **Tax on goods and services**
  - Taxes on goods and services as a share of GDP, annual (TAXGOODSERV, PCGDP)

- **Population and Labor Force**
  - Population 15–64, persons, annual

All data are deflated by the GDP deflator. Data on durables are available for different time spans and frequency. When data was available at quarterly frequency, the series of durables were computed by regressing durables on a constant, Gross Capital Formation (ITISK) and Gross Domestic Product (GDP) in logs, for the available time span, and then using the coefficient estimates to compute the series for durables from the beginning of sample. When data on durables were only available at annual frequency, quarterly observations were estimated using maximum likelihood estimates of a state space model. Population data are interpolated from annual to quarterly frequency using cubic splines. The key per-capita constructions are:

- **per capita output (y):** real GDP − sales taxes + services from consumer durables (with return = 4%) + depreciation of durables (at an annualized rate of 25%), deflated by the GDP deflator and divided by population 16–64.
- **per capita hours (h):** hours worked × total employment, divided by population 16–64.
- **per capita investment (x):** gross capital formation + personal consumption expenditures on durables net of sales taxes, all deflated by the GDP deflator and divided by population 16–64.
- **per capita government consumption (g):** government final consumption expenditures + Exports of goods and services − Imports of goods and services, all deflated by the GDP deflator and divided by population 16–64.

## 5. Parametrization and Calibration

**Table 1: Parameters held fixed across countries**

| β     | δ    | ψ   | θ | α    |
|-------|------|-----|---|------|
| 0.975 | 0.05 | 2.5 | 1 | 0.33 |

where β is the (annualized) discount factor and δ the (annualized) depreciation rate of capital.

**Table 2: Country-specific parameters** — *n* is the average growth rate of population, *γ* the growth rate of labor-augmenting technology, and *a* the adjustment-costs coefficient.

| Country        | n    | γ    | a     |
|----------------|------|------|-------|
| Australia      | 1.44 | 2.08 | 11.65 |
| Austria        | 0.51 | 2.16 | 12.90 |
| Belgium        | 0.35 | 2.02 | 13.43 |
| Canada         | 1.08 | 1.72 | 12.69 |
| Denmark        | 0.26 | 1.99 | 13.64 |
| Finland        | 0.19 | 3.04 | 12.06 |
| France         | 0.46 | 1.77 | 13.67 |
| Germany        | 0.11 | 1.94 | 14.45 |
| Iceland        | 1.22 | 2.39 | 11.52 |
| Ireland        | 1.37 | 4.65 | 9.07  |
| Israel         | 1.98 | 2.15 | 10.87 |
| Italy          | 0.20 | 1.51 | 14.72 |
| Japan          | 0.09 | 2.09 | 14.12 |
| Korea          | 1.27 | 5.07 | 8.82  |
| Luxembourg     | 1.32 | 3.58 | 10.06 |
| Mexico         | 1.79 | 0.73 | 13.14 |
| Netherlands    | 0.46 | 2.41 | 12.60 |
| New Zealand    | 1.18 | 1.62 | 12.70 |
| Norway         | 0.81 | 2.35 | 12.15 |
| Spain          | 0.75 | 2.25 | 12.39 |
| Sweden         | 0.40 | 2.18 | 13.06 |
| Switzerland    | 0.85 | 1.38 | 13.67 |
| United Kingdom | 0.26 | 2.46 | 12.82 |
| United States  | 0.98 | 1.90 | 12.55 |

*n and γ multiplied by 100.*

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