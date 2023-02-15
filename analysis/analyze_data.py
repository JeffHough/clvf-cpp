import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

# load the csv:
df = pd.read_csv("../build/test.csv")
print(df.columns)
print(df.shape)

# extract desired data:
docking_port = df.filter(like="target_d_vector_I_")
chaser_relative_position = df[["chaser_relative_position_0","chaser_relative_position_1","chaser_relative_position_2"]]
chaser_relative_position_B = df[["chaser_relative_position_B_0","chaser_relative_position_B_1","chaser_relative_position_B_2"]]
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

# 3d plot for the relative position:
fig = plt.figure()
ax0 = fig.add_subplot(1, 2, 1, projection='3d')
ax0.plot(
    chaser_relative_position["chaser_relative_position_0"],
    chaser_relative_position["chaser_relative_position_1"],
    chaser_relative_position["chaser_relative_position_2"],
)
ax0.grid()
ax1 = fig.add_subplot(1, 2, 2, projection='3d')
ax1.plot(
    chaser_relative_position_B["chaser_relative_position_B_0"],
    chaser_relative_position_B["chaser_relative_position_B_1"],
    chaser_relative_position_B["chaser_relative_position_B_2"],
)
ax1.grid()


#########################################################################################
####################### ANIMATION #######################################################
#########################################################################################
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')

points, = ax.plot([], [], [], 'o')

def update(index):
    points.set_data(chaser_relative_position_B.values[:index,:2].T)
    points.set_3d_properties(chaser_relative_position_B.values[:index,:2])
    return points,

ani=animation.FuncAnimation(
    fig, 
    update, 
    fargs=range(time.size), 
    interval=10, 
    blit=True, 
    repeat=True
)

plt.show()