import numpy as np
import pandas as pd
import mph

# Import data
df = pd.read_csv('2023-11-04_human_eyes_set_measures.csv')

# Start Comsol and import model
client = mph.start(cores=1)
model = client.load('../../Comsol_files/mfpt_model_splitting_prob_set_human_eyes.mph')
# Change of parameters to model IgG molecule
model.parameter('D', str(0.64e-10))
model.parameter('hya_perm', str(0.874e-7))
model.parameter('ilm_perm', str(1.19e-9))
model.build()
print(model.parameters())

# Define fixed parameters
# a = 1.1275
# b = 0.8885
lens_diam = 0.939
lens_thick = 0.3909
h_va = 0.251

# Assign values to varying parameters
a_list = df['Semi-axis a (cm)'].values
b_list = df['Semi-axis b (cm)'].values
vol_list = df['Vitreous volume (mL)'].values
al_list = df['Axial length (cm)'].values

# Solve for varying parameters
for i in range(0, len(a_list)):
    a = a_list[i]
    b = b_list[i]

    # Calculating params for different initial positions
    y_lens = np.sqrt(b**2*(1 - (lens_diam/2)**2/a**2))
    y_lens_hollow = lens_thick/2 - y_lens
    h_lens = b-y_lens
    h_os = b - (h_lens + h_va)
    y_mid = (y_lens_hollow + b)/2
    y_q1 = b - (b - y_lens_hollow)/4
    y_q2 = b - 3*(b - y_lens_hollow)/4
    x_os = np.sqrt(a**2*(1 - (h_os)**2/b**2))

    # Assigning parameter values in Comsol model
    model.parameter('a', str(a) + '[cm]')
    model.parameter('b', str(b) + '[cm]')
    model.parameter('lens_diam', str(lens_diam) + '[cm]')
    model.parameter('lens_thick', str(lens_thick) + '[cm]')
    model.parameter('h_va', str(h_va) + '[cm]')
    # model.build()
    # print(model.parameters())
    # Solver
    model.solve()
    [x, y, z, u] = model.evaluate(['x', 'y', 'z', 'u'])

    # Finding the arguments for mfpt at initial locations:
    x0_arg = np.nonzero(x == 0)
    y0_arg = np.nonzero(y == 0)
    z0_arg = np.nonzero(z == 0)
    z_mid = np.nonzero(np.isclose(z, y_mid))
    z_mid_q1 = np.nonzero(np.isclose(z, y_q1))
    z_mid_q2 = np.nonzero(np.isclose(z, y_q2))
    y_lens_diam = np.nonzero(np.isclose(y, lens_diam/2))
    z_lens_diam_q1 = np.nonzero(np.isclose(z, -y_lens/2))
    z_lens_diam_q2 = np.nonzero(np.isclose(z, y_lens/2))
    y_os_arg = np.nonzero(np.isclose(y, x_os))

    xy_intersect_0 = np.intersect1d(x0_arg, y0_arg)
    xz_intersect_0 = np.intersect1d(x0_arg, z0_arg)
    xz_intersect_0_zmid = np.intersect1d(x0_arg, z_mid)
    xz_intersect_0_zmidq1 = np.intersect1d(x0_arg, z_mid_q1)
    xz_intersect_0_zmidq2 = np.intersect1d(x0_arg, z_mid_q2)
    xy_intersect_0_ylens = np.intersect1d(x0_arg, y_lens_diam)
    xz_intersect_0_zlens_q1 = np.intersect1d(x0_arg, z_lens_diam_q1)
    xz_intersect_0_zlens_q2 = np.intersect1d(x0_arg, z_lens_diam_q2)
    xy_intersect_0_yos = np.intersect1d(x0_arg, y_os_arg)

    arg_mid = np.intersect1d(xy_intersect_0, xz_intersect_0_zmid)
    arg_mid_q1 = np.intersect1d(xy_intersect_0, xz_intersect_0_zmidq1)
    arg_mid_q2 = np.intersect1d(xy_intersect_0, xz_intersect_0_zmidq2)
    arg_lens_diam_mid = np.intersect1d(xz_intersect_0, xy_intersect_0_ylens)
    arg_lens_diam_q1 = np.intersect1d(xz_intersect_0_zlens_q1,
                                      xy_intersect_0_ylens)
    arg_lens_diam_q2 = np.intersect1d(xz_intersect_0_zlens_q2,
                                      xy_intersect_0_ylens)
    arg_xos_pt = np.intersect1d(xz_intersect_0, xy_intersect_0_yos)

    # Evaluating solution these points:
    mfpt_pt1 = u[arg_mid]
    print(mfpt_pt1)
    mfpt_pt2 = u[arg_mid_q1]
    mfpt_pt3 = u[arg_mid_q2]
    mfpt_pt4 = u[arg_lens_diam_mid]
    mfpt_pt5 = u[arg_lens_diam_q1]
    mfpt_pt6 = u[arg_lens_diam_q2]
    mfpt_pt7 = u[arg_xos_pt]
    mfpt = np.array([mfpt_pt1[0], mfpt_pt2[0], mfpt_pt3[0], mfpt_pt4[0],
                    mfpt_pt5[0], mfpt_pt6[0], mfpt_pt7[0]])

    df = pd.DataFrame({'a': [a],
                       'b': [b],
                       'AL': [al_list[i]],
                       'volume': [vol_list[i]],
                       '#Prop P1': mfpt[0],
                       '#Prop P2': mfpt[1],
                       '#Prop P3': mfpt[2],
                       '#Prop P4': mfpt[3],
                       '#Prop P5': mfpt[4],
                       '#Prop P6': mfpt[5],
                       '#Prop P7': mfpt[6]
                       })
    df.to_csv("2023-11-13_solved_igg_prob_exit_retina_human_eyes_set.csv",
              index=False,
              header=(i == 0),
              mode='a')
