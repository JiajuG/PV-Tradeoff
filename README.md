# PV-Tradeoff

This repository contains the pseudocode for the Pareto-constrained ecological-development co-optimization framework used to generate continuous photovoltaic development pathways.

## Pareto-Constrained Optimization of PV Layouts

The framework combines vegetation-greening probabilities predicted by XGBoost with PV development suitability predicted by MaxEnt. Candidate pathways are generated through a systematic search over the full weight domain. Pareto non-dominated solutions and their distances to the ideal point are then used to derive a fixed compromise weight for each vegetation-response scenario.

Uncertainty in vegetation responses is represented by the 10th, 50th, and 90th percentile scenarios.

```text
INPUTS:
    PV development suitability for pixel i: S_maxent[i]
    Annual energy yield for pixel i: Y[i]
    Area of pixel i: A[i]
    Vegetation-greening probability under scenario q: P_green[i, q]
    Vegetation-response class under scenario q: L[i, q]
    Vegetation-response scenarios Q = {10th, 50th, 90th}
    Candidate weights Alpha = {0.00, 0.01, ..., 1.00}
    Cumulative annual energy-yield evaluation targets T

FOR each vegetation-response scenario q in Q:

    Retain pixels with valid suitability, energy-yield,
    vegetation-response, and area data.

    FOR each valid pixel i:
        EcologicalCost[i, q] = 1 - P_green[i, q]
        DevelopmentCost[i]  = 1 - S_maxent[i]

    EcoRank[:, q] = percentile_rank(EcologicalCost[:, q])
    DevRank[:]    = percentile_rank(DevelopmentCost[:])

    Assign average ranks to tied values.
    CandidateSolutions = empty table

    FOR each alpha in Alpha:

        FOR each valid pixel i:
            PriorityScore[i] =
                alpha * DevRank[i]
                + (1 - alpha) * EcoRank[i, q]

        Sort all valid pixels in ascending order of PriorityScore
        to obtain one complete development sequence.

        Along this sequence, calculate at every position k:

            CumulativeEnergy[k] = SUM(Y[i]) / 10^9

            CumulativeBrowningArea[k] =
                SUM(A[i] * Indicator(L[i, q] = browning))

            CumulativeDevelopmentCost[k] =
                SUM(DevelopmentCost[i])

        FOR each cumulative annual energy-yield target t in T:
            Find the sequence position k nearest to t.

            Add the following record to CandidateSolutions:
                target = t
                weight = alpha
                browning_area = CumulativeBrowningArea[k]
                development_cost = CumulativeDevelopmentCost[k]

    TargetSpecificWeights = empty list

    FOR each cumulative annual energy-yield target t in T:

        Extract all candidate solutions evaluated at target t.

        Identify the Pareto non-dominated solutions by jointly
        minimizing cumulative browning area and cumulative
        development cost.

        Normalize both objectives within the Pareto set using
        min-max normalization.

        FOR each solution j in the normalized Pareto set:
            IdealDistance[j] = SQRT(
                NormalizedBrowningArea[j]^2
                + NormalizedDevelopmentCost[j]^2
            )

        Select the Pareto solution with the minimum IdealDistance.
        Record its weight as the target-specific compromise weight.

    Average the target-specific compromise weights over the focal
    planning evaluation nodes.

    Map the average to the nearest value in Alpha to obtain the
    scenario-specific fixed compromise weight alpha_fixed[q].

    Generate the final pathways:

        Scenario A, ecological priority:
            alpha = 0

        Scenario B, development priority:
            alpha = 1

        Scenario C, Pareto compromise:
            alpha = alpha_fixed[q]

    FOR each of Scenarios A, B, and C:
        Calculate the priority score using its fixed weight.
        Sort the valid pixels once in ascending score order.
        Use successive prefixes of this single sequence to represent
        increasing cumulative annual energy-yield levels.

        The resulting layouts are spatially nested: each layout at a
        lower energy-yield level is retained within layouts at higher
        energy-yield levels.

OUTPUTS:
    Scenario-specific fixed compromise weights
    Continuous development sequences
    Cumulative annual energy-yield trajectories
    Cumulative browning-area trajectories
    Cumulative development-cost trajectories
    Spatial PV-layout rasters
```

This repository documents the optimization procedure and does not contain the original spatial datasets or model-training data.
