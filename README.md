# Integrated Energy Model (IEM) Version 1.1

This repository contains the code for the Integrated Energy Model (IEM), version 1.1. It includes components for data input, cost optimization, and data output.

## 1. Input Data Preparation
The input section prepares all necessary data for the IEM, including:

- Fixed and variable costs for various electricity technologies.
- Generation mixes for different electricity technologies.
- Future IAM-projected targets for wind and solar penetration under different climate scenarios.
- Hourly power demand and capacity factors for wind and solar across various global climate models and scenarios.

Due to memory space limitations, only hourly capacity factors for key countries are uploaded in this repository. However, we are open to additional data requests for other countries.

## 2. Cost Optimization
The cost optimization section provides code for electricity dispatch modeling, with the objective of minimizing the overall electricity system cost. This includes:

-  Cost minimization of capital, fixed O&M (Operation and Maintenance), and variable costs for different technologies, including nuclear, hydropower, oil, gas, coal, wind, solar, long-duration storage, and short-duration storage.
- The optimization process is subject to constraints such as supply-demand balance, energy output limits, power generation constraints, and energy storage.

The least-cost optimization should be performed using the **Gurobi solver** (version 10.0 and above) and the **cvxpy Python package** (version 1.6.0 and above), with **Python** (version 3.12 and above).

## 3. Data Output
The output files will be shown after running cost optimization, which includes:

- Annual overall system cost for a given country and climate scenario.
- Annualized installed capacities various technologies.
- Hourly dispatch patterns various technologies.


## Contact Information
For more details, please contact:

- **zhengds22@mails.tsinghua.edu.cn**
- **dantong@tsinghua.edu.cn**

