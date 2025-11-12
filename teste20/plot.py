import pandas as pd
import matplotlib.pyplot as plt

# Load your CSV
data = pd.read_csv("output_python.csv")

# Extract first two columns
time = data.iloc[:, 0]   # first column (Time)
viremia = data.iloc[:, 1]  # second column (V)

# Plot
plt.figure(figsize=(10,6))
plt.plot(time, viremia, label="Viremia (from output_python.csv)", color="purple")
plt.xlabel("Time (days)")
plt.ylabel("Viremia")
plt.title("Viremia vs Time")
plt.legend()
plt.grid(True)
plt.show()

