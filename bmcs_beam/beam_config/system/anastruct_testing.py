from anastruct import SystemElements
import matplotlib.pyplot as plt
from bmcs_beam.beam_config.system.anastruct_custom_plotter import CustomPlotter
L = 1000
L_F = 260
n_x = 20
F = -2000

n_x_to_F = round((L_F/L) * n_x)
n_x_left = n_x - n_x_to_F

print(n_x_to_F)
print(n_x_left)

struct = SystemElements()
struct.plotter = CustomPlotter(struct, mesh=50)

n_x_to_F = round(0.5* n_x)
print(n_x_to_F)
n_x_remains = n_x - n_x_to_F
struct.add_multiple_elements([[0, 0], [L / 2, 0]], n_x_to_F)
struct.add_multiple_elements([[L / 2, 0], [L, 0]], n_x_remains)

struct.add_support_hinged(1)
struct.add_support_roll(n_x + 1)

struct.point_load(n_x_to_F + 1, Fy=F)
struct.q_load(100, 1)
struct.solve()

struct.show_structure()
# struct.show_displacement()
# axes = struct.plotter.plot_structure(figsize=(12, 8),
#         verbosity=0, show=False).get_axes()
# struct.show_shear_force()
struct.show_bending_moment()
# struct.show_reaction_force()


# shear_xy = struct.show_shear_force(values_only=True, factor=1)
#
# # Getting real internal forces values for 1 element
# normal_force = struct.get_element_results(element_id=1, verbose=True)['N']
# shear = struct.get_element_results(element_id=1, verbose=True)['Q']
# moment = struct.get_element_results(element_id=1, verbose=True)['M']

# moment = struct.get_element_result_range(unit = "moment")
# print(moment)

fig1, ax1 = plt.subplots()

# fig1 = plt.figure(figsize=(12, 8), dpi=300)
# fig1.set_figsize((12, 8))
# ax1 = fig1.add_subplot(111)

# plt.tight_layout()
struct.plotter.plot_structure(figsize=(12, 8),
        verbosity=1, show=False, ax=ax1)
fig1.show()