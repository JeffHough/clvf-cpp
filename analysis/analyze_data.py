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
target_q_BI = df[
    ["target_q_BI_x", "target_q_BI_y", "target_q_BI_z", "target_q_BI_w"]
]

print(docking_port.head())
print(chaser_relative_position.head())
print(target_q_BI.head())

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

fig, axs = plt.subplots(4,1)
for i in range(4):
    target_q_BI_i = target_q_BI.values[:, i]
    ax = axs[i]
    ax.plot(time, target_q_BI_i)
    ax.grid()
    ax.legend()

plt.show()