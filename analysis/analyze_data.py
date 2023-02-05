import pandas as pd
import matplotlib.pyplot as plt

# load the csv:
df = pd.read_csv("../build/test.csv")
print(df.columns)
print(df.shape)

# get the docking port
docking_port = df.filter(like="target_d_vector_I_")
chaser_relative_position = df.filter(like="chaser_relative_position_")
time = df.filter(like="time").values
in_CLVF = df.filter(like="in_CLVF")
target_C_BI = df[
    ["target_C_BI_00", "target_C_BI_01","target_C_BI_02",
    "target_C_BI_10","target_C_BI_11","target_C_BI_12",
    "target_C_BI_20","target_C_BI_21","target_C_BI_22"]
]

print(docking_port.head())
print(chaser_relative_position.head())
print(target_C_BI.head())

# plot the docking port, and the relative position over time:
fig, axs = plt.subplots(4,1)
for i in range(3):
    docking_port_i = docking_port.values[:,i]
    chaser_relative_position_i = chaser_relative_position.values[:,i]
    ax = axs[i]
    ax.plot(time, docking_port_i, label="d")
    ax.plot(time, chaser_relative_position_i, label="r")
    ax.grid()
    ax.legend()
axs[3].plot(time, in_CLVF.values)
axs[3].grid()

fig, axs = plt.subplots(3,3)
for i in range(3):
    for j in range(3):
        target_C_BI_ij = target_C_BI.values[:, 3*i+j]
        ax = axs[i,j]
        ax.plot(time, target_C_BI_ij)
        ax.grid()
        ax.legend()

plt.show()