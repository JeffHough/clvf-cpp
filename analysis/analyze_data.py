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
fig.suptitle("Inertial Position")

# plot the docking port, and the relative position over time:
fig, axs = plt.subplots(4,1)
for i in range(3):
    chaser_relative_position_i = chaser_relative_position_B.values[:,i]
    ax = axs[i]
    ax.plot(time, chaser_relative_position_i, label="r")
    ax.grid()
    ax.legend()
axs[3].plot(time, in_CLVF.values)
axs[3].grid()
fig.suptitle("Body-Fixed Position")

fig, axs = plt.subplots(4,1)
for i in range(4):
    target_q_BI_i = target_q_BI.values[:, i]
    ax = axs[i]
    ax.plot(time, target_q_BI_i)
    ax.grid()
    ax.legend()
fig.suptitle("Quaternion Values")

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
fig.suptitle("3D path")


#########################################################################################
####################### ANIMATION #######################################################
#########################################################################################
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')

# # plot the docking port position and direction:
# ax.quiver(
#     docking_port_B_x,
#     docking_port_B_y,
#     docking_port_B_z,
#     o_vec_body_x,
#     o_vec_body_y,
#     o_vec_body_z
# )

points, = ax.plot([], [], [], 'o')
fig.suptitle("Body-Fixed Animation")

def update(index: int):
    x = chaser_relative_position_B["chaser_relative_position_B_0"].values[:index]
    y = chaser_relative_position_B["chaser_relative_position_B_1"].values[:index]
    z = chaser_relative_position_B["chaser_relative_position_B_2"].values[:index]

    points.set_data(
        x,
        y,
    )

    points.set_3d_properties(z)
    
    if index > 0:
        ax.set_xbound(lower=np.min(x),upper=np.max(x))
        ax.set_ybound(lower=np.min(y), upper=np.max(y))
        ax.set_zlim(np.min(z), np.max(z))
    return points,

ani=animation.FuncAnimation(
    fig, 
    update, 
    frames=[x for x in range(time.shape[0]) if np.mod(x,100)==0],
    interval=0

)

plt.show()