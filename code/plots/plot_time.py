import matplotlib.pyplot as plt
import numpy as np

# Dati della tabella
nodes = np.array([100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 900000])
naive = np.array([0.001, 0.002, 0.003, 0.015, 0.053, 0.350, 1.250, 16.460, 34.200])
path = np.array([0.002, 0.004, 0.006, 0.028, 0.056, 0.310, 0.690, 4.700, 8.510])

nodes = np.array([100, 500, 1000, 2000, 3000, 5000, 8000, 10000])
naive = np.array([0.333, 8.618, 32.698, 190.387, 393.587, 1870.175, 126308.081, 225162.195])
path = np.array([0.338, 2.046, 5.574, 10.546, 38.354, 59.348, 187.951, 208.295])

# Plot comparativo in scala log-log
plt.figure(figsize=(8,6))

plt.plot(nodes, naive, marker="o", linestyle="-", color="blue", label="Base reduction")
plt.plot(nodes, path, marker="s", linestyle="--", color="orange", label="Improved reduction")

# plt.xscale("log")

plt.xlabel("# Nodes")
plt.ylabel("Execution Time (ms)")
plt.title("Execution Time Comparison: Base vs Improved reduction")
plt.legend()

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig("output_plots/execution_time_comparison_reduction.pdf", format="pdf", bbox_inches="tight")
# plt.show()

# Calculate speedup (naive time / path time)
speedup = naive / path

# Speedup plot in a separate figure
plt.figure(figsize=(8,6))

plt.plot(nodes, speedup, marker="^", linestyle="-", color="green", linewidth=2, markersize=8, label="Speedup (Base/Improved)")

# plt.xscale("log")
plt.xlabel("# Nodes")
plt.ylabel("Speedup Factor")
plt.title("Speedup: Base vs Improved reduction")
plt.legend()

# Add horizontal line at y=1 for reference
plt.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='No speedup (1x)')
plt.legend()

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig("output_plots/speedup_comparison_reduction.pdf", format="pdf", bbox_inches="tight")
# plt.show()