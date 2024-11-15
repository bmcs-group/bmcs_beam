from anastruct.fem.plotter.mpl import Plotter
import numpy as np
import math
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.patches as mpatches  # type: ignore
from anastruct.basic import find_nearest, rotate_xy
from anastruct.fem.plotter.values import (
    PlottingValues,
    det_scaling_factor,
    plot_values_axial_force,
    plot_values_bending_moment,
    plot_values_deflection,
    plot_values_element,
    plot_values_shear_force,
)

PATCH_SIZE = 0.03

class CustomPlotter(Plotter):
    def __init__(self, system, mesh):
        super().__init__(system, mesh)

    def __start_plot(self, figsize):
        plt.close("all")
        self.fig = plt.figure(figsize=figsize)
        self.one_fig = self.fig.add_subplot(111)
        plt.tight_layout()

    def __fixed_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        width = height = PATCH_SIZE * max_val
        for node in self.system.supports_fixed:
            support_patch = mpatches.Rectangle(
                (node.vertex.x - width * 0.5, -node.vertex.z - width * 0.5),
                width,
                height,
                color="r",
                zorder=9,
            )
            ax.add_patch(support_patch)

    def __hinged_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        radius = PATCH_SIZE * max_val
        for node in self.system.supports_hinged:
            support_patch = mpatches.RegularPolygon(
                (node.vertex.x, node.vertex.y - radius),
                numVertices=3,
                radius=radius,
                color="r",
                zorder=9,
            )
            ax.add_patch(support_patch)

    def __rotational_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        width = height = PATCH_SIZE * max_val
        for node in self.system.supports_rotational:
            support_patch = mpatches.Rectangle(
                (node.vertex.x - width * 0.5, -node.vertex.z - width * 0.5),
                width,
                height,
                color="r",
                zorder=9,
                fill=False,
            )
            ax.add_patch(support_patch)

    def __roll_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        radius = PATCH_SIZE * max_val
        count = 0
        for node in self.system.supports_roll:
            direction = self.system.supports_roll_direction[count]
            rotate = self.system.supports_roll_rotate[count]
            x1 = np.cos(np.pi) * radius + node.vertex.x + radius  # vertex.x
            z1 = np.sin(np.pi) * radius + node.vertex.y  # vertex.y
            x2 = (
                np.cos(np.radians(90)) * radius + node.vertex.x + radius
            )  # vertex.x + radius
            z2 = np.sin(np.radians(90)) * radius + node.vertex.y  # vertex.y + radius
            x3 = (
                np.cos(np.radians(270)) * radius + node.vertex.x + radius
            )  # vertex.x + radius
            z3 = np.sin(np.radians(270)) * radius + node.vertex.y  # vertex.y - radius

            triangle = np.array([[x1, z1], [x2, z2], [x3, z3]])

            if node.id in self.system.inclined_roll:
                angle = self.system.inclined_roll[node.id]
                triangle = rotate_xy(triangle, angle + np.pi * 0.5)
                support_patch = plt.Polygon(triangle, color="r", zorder=9)
                ax.add_patch(support_patch)
                ax.plot(
                    triangle[1:, 0] - 0.5 * radius * np.sin(angle),
                    triangle[1:, 1] - 0.5 * radius * np.cos(angle),
                    color="r",
                )
                if not rotate:
                    rect_patch = mpatches.RegularPolygon(
                        (node.vertex.x, node - node.vertex.y),
                        numVertices=4,
                        radius=radius,
                        orientation=angle,
                        color="r",
                        zorder=9,
                        fill=False,
                    )
                    ax.add_patch(rect_patch)

            elif direction == 2:  # horizontal roll
                support_patch = mpatches.RegularPolygon(
                    (node.vertex.x, node.vertex.y - radius),
                    numVertices=3,
                    radius=radius,
                    color="r",
                    zorder=9,
                )
                ax.add_patch(support_patch)
                y = -node.vertex.z - 2 * radius
                ax.plot(
                    [node.vertex.x - radius, node.vertex.x + radius], [y, y], color="r"
                )
                if not rotate:
                    rect_patch = mpatches.Rectangle(
                        (node.vertex.x - radius / 2, -node.vertex.z - radius / 2),
                        radius,
                        radius,
                        color="r",
                        zorder=9,
                        fill=False,
                    )
                    ax.add_patch(rect_patch)
            elif direction == 1:  # vertical roll
                # translate the support to the node

                support_patch = mpatches.Polygon(triangle, color="r", zorder=9)
                ax.add_patch(support_patch)

                y = node.vertex.y - radius
                ax.plot(
                    [node.vertex.x + radius * 1.5, node.vertex.x + radius * 1.5],
                    [y, y + 2 * radius],
                    color="r",
                )
                if not rotate:
                    rect_patch = mpatches.Rectangle(
                        (node.vertex.x - radius / 2, -node.vertex.z - radius / 2),
                        radius,
                        radius,
                        color="r",
                        zorder=9,
                        fill=False,
                    )
                    ax.add_patch(rect_patch)
            count += 1

    def __rotating_spring_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        radius = PATCH_SIZE * max_val

        for node, _ in self.system.supports_spring_y:
            r = np.arange(0, radius, 0.001)
            theta = 25 * np.pi * r / (0.2 * max_val)
            x = np.cos(theta) * r + node.vertex.x
            y = np.sin(theta) * r - radius + node.vertex.y
            ax.plot(x, y, color="r", zorder=9)

            # Triangle
            support_patch = mpatches.RegularPolygon(
                (node.vertex.x, node.vertex.y - radius * 3),
                numVertices=3,
                radius=radius * 0.9,
                color="r",
                zorder=9,
            )
            ax.add_patch(support_patch)

    def __spring_support_patch(self, max_val, ax=None):
        """
        :param max_val: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        h = PATCH_SIZE * max_val
        left = -0.5 * h
        right = 0.5 * h
        dh = 0.2 * h

        for node, _ in self.system.supports_spring_z:
            yval = np.arange(0, -9, -1) * dh + node.vertex.y
            xval = (
                np.array([0, 0, left, right, left, right, left, 0, 0]) + node.vertex.x
            )

            ax.plot(xval, yval, color="r", zorder=10)

            # Triangle
            support_patch = mpatches.RegularPolygon(
                (node.vertex.x, -node.vertex.z - h * 2.6),
                numVertices=3,
                radius=h * 0.9,
                color="r",
                zorder=10,
            )
            ax.add_patch(support_patch)

        for node, _ in self.system.supports_spring_x:
            xval = np.arange(0, 9, 1) * dh + node.vertex.x
            yval = (
                np.array([0, 0, left, right, left, right, left, 0, 0]) + node.vertex.y
            )

            ax.plot(xval, yval, color="r", zorder=10)

            # Triangle
            support_patch = mpatches.RegularPolygon(
                (node.vertex.x + h * 1.7, -node.vertex.z - h),
                numVertices=3,
                radius=h * 0.9,
                color="r",
                zorder=10,
            )
            ax.add_patch(support_patch)

    def __q_load_patch(self, max_val, verbosity, ax=None):
        """
        :param max_val: max scale of the plot

        xn1;yn1  q-load   xn1;yn1
        -------------------
        |__________________|
        x1;y1  element    x2;y2
        """
        if ax is None:
            ax = self.one_fig
        for q_id in self.system.loads_q.keys():
            el = self.system.element_map[q_id]
            if el.q_load[0] > 0:
                direction = 1
            else:
                direction = -1

            h1 = 0.05 * max_val * abs(el.q_load[0]) / self.max_q
            h2 = 0.05 * max_val * abs(el.q_load[1]) / self.max_q
            x1 = el.vertex_1.x
            y1 = el.vertex_1.y
            x2 = el.vertex_2.x
            y2 = el.vertex_2.y

            if el.q_direction == "y":
                ai = 0
            elif el.q_direction == "x":
                ai = 0.5 * np.pi
            else:
                ai = -el.angle

            # - value, because the positive z of the beam is opposite of positive y of the plotter
            xn1 = x1 + np.sin(ai) * h1 * direction
            yn1 = y1 + np.cos(ai) * h1 * direction
            xn2 = x2 + np.sin(ai) * h2 * direction
            yn2 = y2 + np.cos(ai) * h2 * direction
            qi = el.q_load[0]
            q = el.q_load[1]
            coordinates = ([x1, xn1, xn2, x2], [y1, yn1, yn2, y2])
            ax.plot(*coordinates, color="g")
            rec = plt.Polygon(np.vstack(coordinates).T, color="g", alpha=0.3)
            ax.add_patch(rec)

            ax.text(xn1, yn1, f"q={qi}", color="b", fontsize=9, zorder=10)
            ax.text(xn2, yn2, f"q={q}", color="b", fontsize=9, zorder=10)

            if verbosity == 0:
                # arrow
                pos = np.sqrt(((y1 - y2) ** 2) + ((x1 - x2) ** 2))
                cg = ((pos / 3) * (el.q_load[0] + 2 * el.q_load[1])) / (
                    el.q_load[0] + el.q_load[1]
                )
                height = math.sin(el.angle) * cg
                base = math.cos(el.angle) * cg

                len_x1 = np.sin(ai - np.pi) * 0.6 * h1 * direction
                len_x2 = np.sin(ai - np.pi) * 0.6 * h2 * direction
                len_y1 = np.cos(ai - np.pi) * 0.6 * h1 * direction
                len_y2 = np.cos(ai - np.pi) * 0.6 * h2 * direction
                step_x = np.linspace(xn1, xn2, 11)
                step_y = np.linspace(yn1, yn2, 11)
                step_len_x = np.linspace(len_x1, len_x2, 11)
                step_len_y = np.linspace(len_y1, len_y2, 11)
                average_h = (h1 + h2) / 2
                # fc = face color, ec = edge color

                # add multiple arrows to fill load
                for counter in range(len(step_x)):
                    if q + qi >= 0:
                        if counter == 0:
                            shape = "right"
                        elif counter == 10:
                            shape = "left"
                        else:
                            shape = "full"
                    else:
                        if counter == 0:
                            shape = "left"
                        elif counter == 10:
                            shape = "right"
                        else:
                            shape = "full"

                    ax.arrow(
                        step_x[counter],
                        step_y[counter],
                        step_len_x[counter],
                        step_len_y[counter],
                        head_width=average_h * 0.25,
                        head_length=0.4
                        * np.sqrt(step_len_y[counter] ** 2 + step_len_x[counter] ** 2),
                        ec="k",
                        fc="k",
                        shape=shape,
                    )

    @staticmethod
    def __arrow_patch_values(Fx, Fz, node, h):
        """
        :param Fx: (float)
        :param Fz: (float)
        :param node: (Node object)
        :param h: (float) Is a scale variable
        :return: Variables for the matplotlib plotter
        """
        F = (Fx ** 2 + Fz ** 2) ** 0.5
        len_x = Fx / F * h
        len_y = -Fz / F * h
        x = node.vertex.x - len_x * 1.2
        y = node.vertex.y - len_y * 1.2

        return x, y, len_x, len_y, F

    def __point_load_patch(self, max_plot_range, verbosity=0, ax=None):
        """
        :param max_plot_range: max scale of the plot
        """
        if ax is None:
            ax = self.one_fig
        for k in self.system.loads_point:
            Fx, Fz = self.system.loads_point[k]
            F = (Fx ** 2 + Fz ** 2) ** 0.5
            node = self.system.node_map[k]
            h = 0.1 * max_plot_range * F / self.max_system_point_load
            x, y, len_x, len_y, F = self.__arrow_patch_values(Fx, Fz, node, h)

            ax.arrow(
                x,
                y,
                len_x,
                len_y,
                head_width=h * 0.15,
                head_length=0.2 * h,
                ec="b",
                fc="orange",
                zorder=11,
            )
            ax.text(x, y, "F=%d" % F, color="k", fontsize=9, zorder=10)

    def __moment_load_patch(self, max_val, ax=None):
        if ax is None:
            ax = self.one_fig
        h = 0.2 * max_val
        for k, v in self.system.loads_moment.items():

            node = self.system.node_map[k]
            if v > 0:
                ax.plot(
                    node.vertex.x,
                    -node.vertex.z,
                    marker=r"$\circlearrowleft$",
                    ms=25,
                    color="orange",
                )
            else:
                ax.plot(
                    node.vertex.x,
                    -node.vertex.z,
                    marker=r"$\circlearrowright$",
                    ms=25,
                    color="orange",
                )
            ax.text(
                node.vertex.x + h * 0.2,
                -node.vertex.z + h * 0.2,
                "T=%d" % v,
                color="k",
                fontsize=9,
                zorder=10,
            )

    def plot_structure(
        self,
        figsize,
        verbosity,
        show=False,
        supports=True,
        scale=1,
        offset=(0, 0),
        gridplot=False,
        annotations=True,
        ax=None,
        loads=False,
        adjust_ax_width=True,
        adjust_ax_height=True,
    ):
        """
        :param show: (boolean) if True, plt.figure will plot.
        :param supports: (boolean) if True, supports are plotted.
        :param annotations: (boolean) if True, structure annotations are plotted. It includes section name.
        :return:
        """
        ax_is_None = ax is None

        if ax_is_None and not gridplot:
            self.__start_plot(figsize)

        if ax_is_None:
            ax = self.one_fig

        x, y = super().structure()

        max_x = np.max(x)
        min_x = np.min(x)
        max_y = np.max(y)
        min_y = np.min(y)

        center_x = (max_x - min_x) / 2 + min_x + -offset[0]
        center_y = (max_y - min_y) / 2 + min_x + -offset[1]

        max_plot_range = max(max_x, max_y)
        ax_range = max_plot_range * scale

        plusxrange = center_x + ax_range
        plusyrange = center_y + ax_range * figsize[1] / figsize[0]
        minxrange = center_x - ax_range
        minyrange = center_y - ax_range * figsize[1] / figsize[0]

        if adjust_ax_width:
            ax.set(xlim=(minxrange, plusxrange))
        if adjust_ax_height:
            ax.set(ylim=(minyrange, plusyrange))

        for el in self.system.element_map.values():
            x_val, y_val = plot_values_element(el)
            if verbosity == 0:
                ax.plot(x_val, y_val, color="black", marker="s")
                # add node ID to plot
                ax_range = max_plot_range * 0.015
                ax.text(
                    x_val[0] + ax_range,
                    y_val[0] + ax_range,
                    "%d" % el.node_id1,
                    color="g",
                    fontsize=9,
                    zorder=10,
                )
                ax.text(
                    x_val[-1] + ax_range,
                    y_val[-1] + ax_range,
                    "%d" % el.node_id2,
                    color="g",
                    fontsize=9,
                    zorder=10,
                )

                # add element ID to plot
                factor = 0.02 * self.max_val_structure
                x_val = (x_val[0] + x_val[-1]) / 2 - np.sin(el.angle) * factor
                y_val = (y_val[0] + y_val[-1]) / 2 + np.cos(el.angle) * factor

                ax.text(
                    x_val, y_val, str(el.id), color="r", fontsize=9, zorder=10
                )

                # add element annotation to plot
                # TODO: check how this holds with multiple structure scales.
                if annotations:
                    x_val += +np.sin(el.angle) * factor * 2.3
                    y_val += -np.cos(el.angle) * factor * 2.3
                    ax.text(
                        x_val, y_val, el.section_name, color="b", fontsize=9, zorder=10
                    )
            else:
                ax.plot(x_val, y_val, color="black")
        # add supports
        if supports:
            self.__fixed_support_patch(max_plot_range * scale, ax=ax)
            self.__hinged_support_patch(max_plot_range * scale, ax=ax)
            self.__rotational_support_patch(max_plot_range * scale, ax=ax)
            self.__roll_support_patch(max_plot_range * scale, ax=ax)
            self.__rotating_spring_support_patch(max_plot_range * scale, ax=ax)
            self.__spring_support_patch(max_plot_range * scale, ax=ax)

        # add_loads
        if loads or verbosity == 0:
            self.__q_load_patch(max_plot_range, verbosity, ax=ax)
            self.__point_load_patch(max_plot_range, verbosity, ax=ax)
            self.__moment_load_patch(max_plot_range, ax=ax)

        if ax_is_None:
            if show:
                self.plot()
            else:
                return self.fig

    def _add_node_values(self, x_val, y_val, value_1, value_2, digits, ax=None):
        if ax is None:
            ax = self.one_fig
        offset = self.max_val_structure * 0.015

        # add value to plot
        ax.text(
            x_val[1] - offset,
            y_val[1] + offset,
            "%s" % round(value_1, digits),
            fontsize=9,
            ha="center",
            va="center",
        )
        ax.text(
            x_val[-2] - offset,
            y_val[-2] + offset,
            "%s" % round(value_2, digits),
            fontsize=9,
            ha="center",
            va="center",
        )

    def _add_element_values(self, x_val, y_val, value, index, digits=2, ax=None):
        if ax is None:
            ax = self.one_fig

        ax.text(
            x_val[index],
            y_val[index],
            "%s" % round(value, digits),
            fontsize=9,
            ha="center",
            va="center",
        )

    def plot_result(
        self,
        axis_values,
        force_1=None,
        force_2=None,
        digits=2,
        node_results=True,
        fill_polygon=True,
        color=0,
        ax=None,
    ):
        if ax is None:
            ax = self.one_fig
        if fill_polygon:
            rec = plt.Polygon(
                np.vstack(axis_values).T, color="C{}".format(color), alpha=0.3
            )
            ax.add_patch(rec)
        # plot force
        x_val = axis_values[0]
        y_val = axis_values[1]

        ax.plot(x_val, y_val, color="C{}".format(color))

        if node_results:
            self._add_node_values(x_val, y_val, force_1, force_2, digits, ax=ax)

    def plot(self):
        plt.show()

    def axial_force(
        self,
        factor=None,
        figsize=None,
        verbosity=0,
        scale=1,
        offset=(0, 0),
        show=True,
        gridplot=False,
    ):
        self.plot_structure(figsize, 1, scale=scale, offset=offset, gridplot=gridplot)

        if factor is None:
            max_force = max(
                map(
                    lambda el: max(abs(el.N_1), abs(el.N_2)),
                    self.system.element_map.values(),
                )
            )
            factor = det_scaling_factor(max_force, self.max_val_structure)

        for el in self.system.element_map.values():
            if math.isclose(el.N_1, 0, rel_tol=1e-5, abs_tol=1e-9):
                continue
            else:
                axis_values = plot_values_axial_force(el, factor)
                color = 1 if el.N_1 < 0 else 0
                self.plot_result(
                    axis_values,
                    el.N_1,
                    el.N_2,
                    node_results=not bool(verbosity),
                    color=color,
                )

                point = (el.vertex_2 - el.vertex_1) / 2 + el.vertex_1
                if el.N_1 < 0:
                    point.displace_polar(
                        alpha=el.angle + 0.5 * np.pi,
                        radius=0.5 * el.N_1 * factor,
                        inverse_z_axis=True,
                    )

                    if verbosity == 0:
                        self.one_fig.text(
                            point.x,
                            point.y,
                            "-",
                            ha="center",
                            va="center",
                            fontsize=20,
                            color="b",
                        )
                if el.N_1 > 0:
                    point.displace_polar(
                        alpha=el.angle + 0.5 * np.pi,
                        radius=0.5 * el.N_1 * factor,
                        inverse_z_axis=True,
                    )

                    if verbosity == 0:
                        self.one_fig.text(
                            point.x,
                            point.y,
                            "+",
                            ha="center",
                            va="center",
                            fontsize=14,
                            color="b",
                        )

        if show:
            self.plot()
        else:
            return self.fig

    def bending_moment(
        self,
        factor=None,
        figsize=None,
        verbosity=0,
        scale=1,
        offset=(0, 0),
        show=True,
        gridplot=False,
        ax=None,
        adjust_ax_height=True,
    ):
        self.plot_structure(figsize, 1, supports=False, adjust_ax_height=adjust_ax_height,
                            scale=scale, offset=offset, gridplot=gridplot, ax=ax)

        ax_is_None = ax is None
        if ax_is_None:
            ax = self.one_fig

        con = len(self.system.element_map[1].bending_moment)
        if factor is None:
            # maximum moment determined by comparing the node's moments and the sagging moments.
            max_moment = max(
                map(
                    lambda el: max(
                        abs(el.node_1.Ty),
                        abs(((el.all_q_load[0] + el.all_q_load[1]) / 16) * el.l ** 2),
                    )
                    if el.type == "general"
                    else 0,
                    self.system.element_map.values(),
                )
            )
            factor = det_scaling_factor(max_moment, self.max_val_structure)

        # determine the axis values
        for el in self.system.element_map.values():
            if (
                math.isclose(el.node_1.Ty, 0, rel_tol=1e-5, abs_tol=1e-9)
                and math.isclose(el.node_2.Ty, 0, rel_tol=1e-5, abs_tol=1e-9)
                and not el.all_q_load[0]
                and not el.all_q_load[1]
            ):
                # If True there is no bending moment, so no need for plotting.
                continue
            axis_values = plot_values_bending_moment(el, factor, con)
            if verbosity == 0:
                node_results = True
            else:
                node_results = False
            self.plot_result(
                axis_values,
                abs(el.node_1.Ty),
                abs(el.node_2.Ty),
                node_results=node_results,
                ax=ax
            )

            if el.all_q_load:
                m_sag = min(el.bending_moment)
                index = find_nearest(el.bending_moment, m_sag)[1]
                offset = -self.max_val_structure * 0.05

                if verbosity == 0:
                    x = axis_values[0][index] + np.sin(-el.angle) * offset
                    y = axis_values[1][index] + np.cos(-el.angle) * offset
                    ax.text(x, y, "%s" % round(m_sag, 1), fontsize=9)

        # Apply equal padding up and below the plot
        self._add_vertical_padding(ax)

        if ax_is_None:
            if show:
                self.plot()
            else:
                return self.fig

    def _add_vertical_padding(self, ax):
        ymin, ymax = ax.get_ylim()
        max_diff = 0.5 * max(abs(ymin), abs(ymax))
        ax.set(ylim=(ymin + np.sign(ymin) * max_diff, ymax + ymax + np.sign(ymax) * max_diff))

    def shear_force(
        self,
        factor=None,
        figsize=None,
        verbosity=0,
        scale=1,
        offset=(0, 0),
        show=True,
        gridplot=False,
        include_structure=True,
        adjust_ax_height=True,
        ax=None,
    ):
        ax_is_None = ax is None
        if ax_is_None:
            ax = self.one_fig

        if include_structure:
            self.plot_structure(
                figsize, 1, scale=scale, offset=offset, gridplot=gridplot, ax=ax,
                adjust_ax_height=adjust_ax_height, supports=False
            )
        if factor is None:
            max_force = max(
                map(
                    lambda el: np.max(np.abs(el.shear_force)),
                    self.system.element_map.values(),
                )
            )
            factor = det_scaling_factor(max_force, self.max_val_structure)

        for el in self.system.element_map.values():
            if (
                math.isclose(el.node_1.Ty, 0, rel_tol=1e-5, abs_tol=1e-9)
                and math.isclose(el.node_2.Ty, 0, rel_tol=1e-5, abs_tol=1e-9)
                and el.q_load is None
            ):
                # If True there is no bending moment and no shear, thus no shear force, so no need for plotting.
                continue
            axis_values = plot_values_shear_force(el, factor)
            shear_1 = el.shear_force[0]
            shear_2 = el.shear_force[-1]

            self.plot_result(
                axis_values, shear_1, shear_2, node_results=not bool(verbosity), ax=ax
            )

        # Apply equal padding up and below the plot
        self._add_vertical_padding(ax)

        if ax_is_None:
            if show:
                self.plot()
            else:
                return self.fig

    def reaction_force(self, figsize, verbosity, scale, offset, show, gridplot=False):
        self.plot_structure(
            figsize, 1, supports=False, scale=scale, offset=offset, gridplot=gridplot
        )

        h = 0.2 * self.max_val_structure
        max_force = max(
            map(
                lambda node: max(abs(node.Fx), abs(node.Fz)),
                self.system.reaction_forces.values(),
            )
        )

        for node in self.system.reaction_forces.values():
            if not math.isclose(node.Fx, 0, rel_tol=1e-5, abs_tol=1e-9):
                # x direction
                scale = abs(node.Fx) / max_force * h
                sol = self.__arrow_patch_values(node.Fx, 0, node, scale)
                x = sol[0]
                y = sol[1]
                len_x = sol[2]
                len_y = sol[3]

                self.one_fig.arrow(
                    x,
                    y,
                    len_x,
                    len_y,
                    head_width=h * 0.15,
                    head_length=0.2 * scale,
                    ec="b",
                    fc="orange",
                    zorder=11,
                )

                if verbosity == 0:
                    self.one_fig.text(
                        x,
                        y,
                        "R=%s" % round(node.Fx, 2),
                        color="k",
                        fontsize=9,
                        zorder=10,
                    )

            if not math.isclose(node.Fz, 0, rel_tol=1e-5, abs_tol=1e-9):
                # z direction
                scale = abs(node.Fz) / max_force * h
                sol = self.__arrow_patch_values(0, node.Fz, node, scale)
                x = sol[0]
                y = sol[1]
                len_x = sol[2]
                len_y = sol[3]

                self.one_fig.arrow(
                    x,
                    y,
                    len_x,
                    len_y,
                    head_width=h * 0.15,
                    head_length=0.2 * scale,
                    ec="b",
                    fc="orange",
                    zorder=11,
                )

                if verbosity == 0:
                    self.one_fig.text(
                        x,
                        y,
                        "R=%s" % round(node.Fz, 2),
                        color="k",
                        fontsize=9,
                        zorder=10,
                    )

            if not math.isclose(node.Ty, 0, rel_tol=1e-5, abs_tol=1e-9):
                """
                '$...$': render the strings using mathtext
                """
                if node.Ty > 0:
                    self.one_fig.plot(
                        node.vertex.x,
                        -node.vertex.z,
                        marker=r"$\circlearrowleft$",
                        ms=25,
                        color="orange",
                    )
                if node.Ty < 0:
                    self.one_fig.plot(
                        node.vertex.x,
                        -node.vertex.z,
                        marker=r"$\circlearrowright$",
                        ms=25,
                        color="orange",
                    )

                if verbosity == 0:
                    self.one_fig.text(
                        node.vertex.x + h * 0.2,
                        -node.vertex.z + h * 0.2,
                        "T=%s" % round(node.Ty, 2),
                        color="k",
                        fontsize=9,
                        zorder=10,
                    )
        if show:
            self.plot()
        else:
            return self.fig

    def displacements(
        self,
        factor=None,
        figsize=None,
        verbosity=0,
        scale=1,
        offset=(0, 0),
        show=True,
        linear=False,
        gridplot=False,
    ):
        self.plot_structure(figsize, 1, scale=scale, offset=offset, gridplot=gridplot)
        if factor is None:
            # needed to determine the scaling factor
            max_displacement = max(
                map(
                    lambda el: max(
                        abs(el.node_1.ux), abs(el.node_1.uz), el.max_deflection
                    )
                    if el.type == "general"
                    else 0,
                    self.system.element_map.values(),
                )
            )
            factor = det_scaling_factor(max_displacement, self.max_val_structure)

        for el in self.system.element_map.values():
            axis_values = plot_values_deflection(el, factor, linear)
            self.plot_result(axis_values, node_results=False, fill_polygon=False)

            if el.type == "general":
                # index of the max deflection
                x = np.linspace(el.vertex_1.x, el.vertex_2.x, el.deflection.size)
                y = np.linspace(el.vertex_1.y, el.vertex_2.y, el.deflection.size)
                xd, yd = plot_values_deflection(el, 1.0, linear)
                deflection = ((xd - x) ** 2 + (yd - y) ** 2) ** 0.5
                index = np.argmax(np.abs(deflection))

                if verbosity == 0:

                    if index != 0 or index != el.deflection.size:
                        self._add_element_values(
                            axis_values[0],
                            axis_values[1],
                            deflection[index],
                            index,
                            3,
                        )
        if show:
            self.plot()
        else:
            return self.fig

    def results_plot(self, figsize, verbosity, scale, offset, show):
        """
        Aggregate all the plots in one grid plot.

        :param figsize: (tpl)
        :param verbosity: (int)
        :param scale: (flt)
        :param offset: (tpl)
        :param show: (bool)
        :return: Figure or None
        """
        plt.close("all")
        self.fig = plt.figure(figsize=figsize)
        a = 320
        self.one_fig = self.fig.add_subplot(a + 1)
        plt.title("structure")
        self.plot_structure(
            figsize, verbosity, show=False, scale=scale, offset=offset, gridplot=True
        )
        self.one_fig = self.fig.add_subplot(a + 2)
        plt.title("bending moment")
        self.bending_moment(None, figsize, verbosity, scale, offset, False, True)
        self.one_fig = self.fig.add_subplot(a + 3)
        plt.title("shear force")
        self.shear_force(None, figsize, verbosity, scale, offset, False, True)
        self.one_fig = self.fig.add_subplot(a + 4)
        plt.title("axial force")
        self.axial_force(None, figsize, verbosity, scale, offset, False, True)
        self.one_fig = self.fig.add_subplot(a + 5)
        plt.title("displacements")
        self.displacements(None, figsize, verbosity, scale, offset, False, False, True)
        self.one_fig = self.fig.add_subplot(a + 6)
        plt.title("reaction force")
        self.reaction_force(figsize, verbosity, scale, offset, False, True)

        if show:
            self.plot()
        else:
            return self.fig