from anastruct import SystemElements
import matplotlib.pyplot as plt

L = 1000
L_F = 260
n_x = 20
F = -2000

n_x_to_F = round((L_F/L) * n_x)
n_x_left = n_x - n_x_to_F

print(n_x_to_F)
print(n_x_left)

struct = SystemElements()

n_x_to_F = round(0.5* n_x)
print(n_x_to_F)
n_x_remains = n_x - n_x_to_F
struct.add_multiple_elements([[0, 0], [L / 2, 0]], n_x_to_F)
struct.add_multiple_elements([[L / 2, 0], [L, 0]], n_x_remains)

struct.add_support_hinged(1)
struct.add_support_roll(n_x + 1)

struct.point_load(n_x_to_F + 1, Fy=F)
struct.solve()

# struct.show_structure()
# struct.show_displacement()
struct.show_shear_force()
struct.show_bending_moment()
# struct.show_reaction_force()


# shear_xy = struct.show_shear_force(values_only=True, factor=1)
#
# # Getting real internal forces values for 1 element
# normal_force = struct.get_element_results(element_id=1, verbose=True)['N']
# shear = struct.get_element_results(element_id=1, verbose=True)['Q']
# moment = struct.get_element_results(element_id=1, verbose=True)['M']

moment = struct.get_element_result_range(unit = "moment")
print(moment)