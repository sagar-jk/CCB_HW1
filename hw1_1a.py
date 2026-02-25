import random

def propensities(x1, x2, x3):
    a1 = 0.5 * x1 * (x1 - 1) * x2 if (x1 >= 2 and x2 >= 1) else 0.0
    a2 = x1 * x3 * (x3 - 1) if (x1 >= 1 and x3 >= 2) else 0.0
    a3 = 3.0 * x2 * x3 if (x2 >= 1 and x3 >= 1) else 0.0
    return a1, a2, a3

def choose_reaction(x1, x2, x3):
    a1, a2, a3 = propensities(x1, x2, x3)
    a0 = a1 + a2 + a3
    if a0 == 0:
        return None

    r = random.random()
    if r < a1 / a0:
        return 1
    elif r < (a1 + a2) / a0:
        return 2
    else:
        return 3

def simulate_once():
    x1, x2, x3 = 110, 26, 55

    while True:
        # Check stopping conditions
        if x1 >= 150:
            return "C1"
        if x2 < 10:
            return "C2"
        if x3 > 100:
            return "C3"

        rxn = choose_reaction(x1, x2, x3)
        if rxn is None:
            return "none"

        # Apply reaction updates
        if rxn == 1:      # R1
            x1 -= 2; x2 -= 1; x3 += 4
        elif rxn == 2:    # R2
            x1 -= 1; x2 += 3; x3 -= 2
        else:             # R3
            x1 += 2; x2 -= 1; x3 -= 1


def estimate_probabilities(n_runs=20000):
    counts = {"C1": 0, "C2": 0, "C3": 0}

    for _ in range(n_runs):
        outcome = simulate_once()
        if outcome in counts:
            counts[outcome] += 1

    Pr_C1 = counts["C1"] / n_runs
    Pr_C2 = counts["C2"] / n_runs
    Pr_C3 = counts["C3"] / n_runs

    # Round to 0 or 1 for display
    Pr_C1 = round(Pr_C1)
    Pr_C2 = round(Pr_C2)
    Pr_C3 = round(Pr_C3)

    print(f"Pr(C1) = {Pr_C1}")
    print(f"Pr(C2) = {Pr_C2}")
    print(f"Pr(C3) = {Pr_C3}")


if __name__ == "__main__":
    estimate_probabilities()