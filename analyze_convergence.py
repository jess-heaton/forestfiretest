#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

def main():
    # 1) Read the data file
    #    We have columns in the order:
    #       p_from_script, M_from_script, nproc, N, p_from_code, M_from_code, avgSteps, avgTime, bottomFraction
    #
    #    Example line:
    #       0.3 10 1 100 0.3 10 22.4 0.129 0.73
    #
    #    We'll assign column names:
    col_names = [
        "p_script",      # Column 0
        "M_script",      # Column 1
        "nproc",         # Column 2
        "N",             # Column 3
        "p_code",        # Column 4
        "M_code",        # Column 5
        "avgSteps",      # Column 6
        "avgTime",       # Column 7
        "bottomFraction" # Column 8
    ]

    # Read in with space delimiter
    df = pd.read_csv("convergence_data.txt", sep=" ", header=None, names=col_names)

    # 2) Quick sanity check: p_script == p_code, M_script == M_code?
    #    (They should match, but if you want to confirm, you can compare them or drop duplicates.)

    # 3) Plot: Average Steps vs. M for each p
    plt.figure(figsize=(7,5))
    for p_val, grp in df.groupby("p_script"):
        # Sort by M just to ensure lines are in ascending order
        grp_sorted = grp.sort_values(by="M_script")
        plt.plot(grp_sorted["M_script"], grp_sorted["avgSteps"], marker='o', label=f"p={p_val}")

    plt.xlabel("Number of Repeats (M)")
    plt.ylabel("Average Steps")
    plt.title("Forest Fire Convergence: Steps vs. M (N=100)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("convergence_steps.png", dpi=150)
    plt.show()

    # 4) Plot: Average Time vs. M for each p (optional but often insightful)
    plt.figure(figsize=(7,5))
    for p_val, grp in df.groupby("p_script"):
        grp_sorted = grp.sort_values(by="M_script")
        plt.plot(grp_sorted["M_script"], grp_sorted["avgTime"], marker='s', label=f"p={p_val}")

    plt.xlabel("Number of Repeats (M)")
    plt.ylabel("Average Time (s)")
    plt.title("Forest Fire Convergence: Time vs. M (N=100)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("convergence_time.png", dpi=150)
    plt.show()

    # 5) (Optional) Plot bottomFraction vs. M for each p
    #    If you want to see how often the fire reached the bottom
    plt.figure(figsize=(7,5))
    for p_val, grp in df.groupby("p_script"):
        grp_sorted = grp.sort_values(by="M_script")
        plt.plot(grp_sorted["M_script"], grp_sorted["bottomFraction"], marker='^', label=f"p={p_val}")

    plt.xlabel("Number of Repeats (M)")
    plt.ylabel("Fraction that reached bottom")
    plt.title("Forest Fire: Bottom Fraction vs. M (N=100)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("convergence_bottom_fraction.png", dpi=150)
    plt.show()

    print("Plots saved to: convergence_steps.png, convergence_time.png, convergence_bottom_fraction.png")

if __name__ == "__main__":
    main()
