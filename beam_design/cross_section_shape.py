import traits.api as tr
from bmcs_utils.api import InteractiveModel
import numpy as np
import sympy as sp
from bmcs_utils.models.interactive_window import View, Item
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


# TODO: decide whether to use this or just remove it, (maybe delete)
class ICrossSectionShape(tr.Interface):
    """This interface lists the functions need to be implemented by cross section classes."""

    def get_cs_area(self):
        """This function should return the area of the cross section."""

    def get_b(self, z_positions_array):
        """This function should return b values that correspond to the positions z_positions_array."""


class CrossSectionShapeBase(InteractiveModel):
    # TODO: handle the initialization of this base class because it doesn't have b to draw the cross section
    name = 'Cross section shape'

    """"This class describes the geometry of the cross section."""
    H = tr.Float(250)

    ipw_view = View(
        Item('H', minmax=(10, 3000), latex='H [mm]')
    )

    def get_area(self):
        """Calculate the integral of b over the height"""
        return self.get_cs_area()

    def get_b(self, z):
        """Width of the cross section must be specified by the subclasses"""
        raise NotImplementedError()

    def subplots(self, fig):
        return fig.subplots(1, 1)

    # TODO: decide on this update_plot if it will be done in subclasses for efficiency or here as a brief implementation
    def update_plot(self, ax):
        raise NotImplementedError()


@tr.provides(ICrossSectionShape)
class Rectangle(CrossSectionShapeBase):
    B = tr.Float(250)

    ipw_view = View(
        *CrossSectionShapeBase.ipw_view.content,        # this will add View Items of the base class CrossSectionShapeBase
        Item('B', minmax=(10, 500), latex='B [mm]')
    )

    def get_cs_area(self):
        return self.B * self.H

    def get_b(self, z_positions_array):
        return np.full_like(z_positions_array, self.B)

    def update_plot(self, ax):
        ax.axis('equal')
        ax.fill([0, self.B, self.B, 0, 0], [0, 0, self.H, self.H, 0], color='gray')
        ax.plot([0, self.B, self.B, 0, 0], [0, 0, self.H, self.H, 0], color='black')


@tr.provides(ICrossSectionShape)
class Circle(CrossSectionShapeBase):
    # TODO->Rostia: provide input field instead minmax range
    # H from the base class is used as the R

    ipw_view = View(
        Item('H', minmax=(10, 3000), latex='R [mm]'),
    )

    def get_cs_area(self):
        return np.pi * self.H * self.H

    def get_b(self, z_positions_array):
        # TODO->Saeed: complete this
        pass

    def update_plot(self, ax):
        # TODO->Saeed: fix this
        ax.axis('equal')
        ax.Circle((0, 0), self.R, color='gray', linewidth=4, fill=True)


@tr.provides(ICrossSectionShape)
class TShape(CrossSectionShapeBase):
    name = 'T-shape'

    B_f = tr.Float(250, input=True)
    B_w = tr.Float(100, input=True)
    H_w = tr.Float(100, input=True)

    ipw_view = View(
        *CrossSectionShapeBase.ipw_view.content,
        Item('B_f', minmax=(10, 3000), latex='B_f [mm]'),
        Item('B_w', minmax=(10, 3000), latex='B_w [mm]'),
        Item('H_w', minmax=(10, 3000), latex='H_w [mm]'),
    )

    def get_cs_area(self):
        return self.B_w * self.H_w + self.B_f * (self.H - self.H_w)

    get_b = tr.Property(tr.Callable, depends_on='+input')

    @tr.cached_property
    def _get_get_b(self):
        z_ = sp.Symbol('z')
        b_p = sp.Piecewise(
            (self.B_w, z_ < self.H_w),
            (self.B_f, True)
        )
        return sp.lambdify(z_, b_p, 'numpy')

    def update_plot(self, ax):
        # TODO->Saeed: fix this
        # Start drawing from bottom center of the cross section
        cs_points = np.array([  [self.B_w/2, 0],
                                [self.B_w/2, self.H_w],
                                [self.B_f/2, self.H_w],
                                [self.B_f/2, self.H],
                                [-self.B_f/2, self.H],
                                [-self.B_f/2, self.H_w],
                                [-self.B_w/2, self.H_w],
                                [-self.B_w/2, 0]])
        cs = Polygon(cs_points, True)
        patch_collection = PatchCollection([cs])

        ax.add_collection(patch_collection)


# TODO->Saeed: maybe complete this
@tr.provides(ICrossSectionShape)
class CustomShape(CrossSectionShapeBase):

    def get_cs_area(self):
        pass

    def get_b(self, z_positions_array):
        pass

    def update_plot(self, ax):
        # TODO->Saeed: fix this to use it for the CustomShape
        # self.update_plot(ax)
        z = np.linspace(0, self.H, 100)
        b = self.get_b(z)
        ax.axis([0, np.max(b), 0, self.H])
        ax.axis('equal')
        ax.fill(b, z, color='gray')
        ax.plot(b, z, color='black')


if __name__ == '__main__':
    cs = CrossSectionShapeBase()

    tv = cs.trait_view('ipw_view')
    cs.configure_traits()