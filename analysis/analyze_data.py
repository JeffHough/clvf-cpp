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

# o_hat_vector for the CLVF:
target_d_vector_B = df[["target_d_vector_B_0", "target_d_vector_B_1","target_d_vector_B_2"]]
target_o_hat_vector_B_CLVF = df[["target_o_hat_vector_B_CLVF_0", "target_o_hat_vector_B_CLVF_1","target_o_hat_vector_B_CLVF_2"]]
target_o_hat_vector_B_LVF = df[["target_o_hat_vector_B_LVF_0","target_o_hat_vector_B_LVF_1","target_o_hat_vector_B_LVF_2"]]

# o_have_vector for the inertial frame:
target_o_hat_vector_I_CLVF = df[["target_o_hat_vector_I_CLVF_0","target_o_hat_vector_I_CLVF_1","target_o_hat_vector_I_CLVF_2"]]

# alpha value for the CLVF:
target_alpha_CLVF = df[["alpha_CLVF"]]

# plot the docking port, and the relative position over time [INERTIAL]:
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

# plot the docking port, and the relative position over time [BODY-FIXED]:
fig, axs = plt.subplots(4,1)
for i in range(3):
    o_hat_target_B_i = target_o_hat_vector_B_CLVF.values[:,i]*target_alpha_CLVF.values[i]
    chaser_relative_position_i = chaser_relative_position_B.values[:,i]
    ax = axs[i]
    ax.plot(time, o_hat_target_B_i, label="target_vector")
    ax.plot(time, chaser_relative_position_i, label="r")
    ax.grid()
    ax.legend()
axs[3].plot(time, in_CLVF.values)
axs[3].grid()
fig.suptitle("Body-Fixed Position")

# plot the o-hat-vector in the inertial frame over time:
fig, axs = plt.subplots(3,1)
for i in range(3):
    o_vector_i = target_o_hat_vector_I_CLVF.values[:,i]
    axs[i].plot(time, o_vector_i)
    axs[i].grid()
    axs[i].set_title("O-Hat in the Inertial Frame")

# plot the quaternion values over time:
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


# Are we closing all the way? Check the size of "r" versus "alpha":
r_norm = np.zeros(time.shape)
for (i, r_vector) in zip(range(r_norm.size), chaser_relative_position.values):
    r_norm[i] = np.sqrt(r_vector[0]**2 + r_vector[1]**2 + r_vector[2]**2)

fig, axs = plt.subplots(1,1)
axs.plot(time, target_alpha_CLVF, label="alpha")
axs.plot(time, r_norm, label="r_norm")
axs.grid()
axs.legend()
axs.set_title("radius and alpha over time")

#########################################################################################
####################### ANIMATION #######################################################
#########################################################################################
indices_to_skip = 200

fig = plt.figure()
ax = fig.add_subplot(1,2,1,projection='3d')

points, = ax.plot([], [], [], 'o')
ax.set_title("Body-Fixed Animation")

# plot the target vector in the body-fixed frame:
target_vector_B = target_o_hat_vector_B_CLVF.values[0,0:3]*target_alpha_CLVF.values[0]
print(f"Target vector is: {target_vector_B}")
ax.quiver([0],[0],[0],[target_vector_B[0]],[target_vector_B[1]],[target_vector_B[2]])

def update(index: int):
    x_B = chaser_relative_position_B["chaser_relative_position_B_0"].values[:index]
    y_B = chaser_relative_position_B["chaser_relative_position_B_1"].values[:index]
    z_B = chaser_relative_position_B["chaser_relative_position_B_2"].values[:index]

    points.set_data(
        x_B,
        y_B,
    )

    points.set_3d_properties(z_B)
    
    if index > 0:
        ax.set_xbound(lower=np.min(x_B),upper=np.max(x_B))
        ax.set_ybound(lower=np.min(y_B), upper=np.max(y_B))
        ax.set_zlim(np.min(z_B), np.max(z_B))
    return points,

# animation for the inertial frame:
ax_2 = fig.add_subplot(1,2,2,projection='3d')
points_2, = ax_2.plot([], [], [], 'o')
ax_2.set_title("Inertial Animation")

target_vector_I = target_o_hat_vector_I_CLVF.values[0,0:3]*target_alpha_CLVF.values[0]
q = ax_2.quiver([0],[0],[0],[target_vector_I[0]],[target_vector_I[1]],[target_vector_I[2]])

def update_inertial(index: int):

    # update the quiver:
    global q
    target_vector_I = target_o_hat_vector_I_CLVF.values[index,0:3]*target_alpha_CLVF.values[0]
    q.remove()
    q = ax_2.quiver([0],[0],[0],[target_vector_I[0]],[target_vector_I[1]],[target_vector_I[2]])

    x = chaser_relative_position["chaser_relative_position_0"].values[:index]
    y = chaser_relative_position["chaser_relative_position_1"].values[:index]
    z = chaser_relative_position["chaser_relative_position_2"].values[:index]

    points_2.set_data(
        x,
        y,
    )

    points_2.set_3d_properties(z)
    
    if index > 0:
        ax_2.set_xbound(lower=np.min(x),upper=np.max(x))
        ax_2.set_ybound(lower=np.min(y), upper=np.max(y))
        ax_2.set_zlim(np.min(z), np.max(z))
    
    return points_2,

ani=animation.FuncAnimation(
    fig, 
    update, 
    frames=[x for x in range(time.shape[0]) if np.mod(x,indices_to_skip)==0],
    interval=0
)
ani_2=animation.FuncAnimation(
    fig, 
    update_inertial, 
    frames=[x for x in range(time.shape[0]) if np.mod(x,indices_to_skip)==0],
    interval=0
)

fig.suptitle("Animations of motion.")

plt.show()