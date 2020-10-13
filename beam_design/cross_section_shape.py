import traits.api as tr
from bmcs_utils.api import InteractiveModel
import numpy as np
import sympy as sp
from bmcs_utils.models.interactive_window import View, Item


class XICrossSectionShape(tr.Interface):
    """This interface lists the functions need to be implemented by cross section classes."""

    def get_cs_area(self):
        """This function should return the area of the cross section."""

    def get_b(self, z_positions_array):
        """This function should return b values that correspond to the positions z_positions_array."""


class CrossSectionShape(tr.HasTraits):
    """"This class describes the geometry of the cross section."""
    H = tr.Float(250)

    def get_area(self):
        '''Calculate the integral of b over the height'''

    def get_b(self, z):
        '''Width of the cross section must
        be specified by the subclasses'''
        raise NotImplementedError()

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        z = np.linspace(0, self.H, 100)
        b = self.get_b(z)
        ax.axis([0, np.max(b), 0, self.H])
        ax.axis('equal')
        ax.fill(b, z, color='gray')
        ax.plot(b, z, color='black')

class Rectangle(CrossSectionShape):
    B = tr.Float(250)

    ipw_view = View(
        Item('H', minmax=(1, 3000), latex='h'),
        Item('B', minmax=(1, 500), latex='b')
    )

    def get_cs_area(self):
        return self.b * self.h

    def get_b(self, z_positions_array):
        return np.full_like(z_positions_array, self.b)

    def subplots(self, fig):
        return super().subplots(fig)

    def xupdate_plot(self, ax):
        ax.axis([0, self.B, 0, self.h])
        ax.axis('equal')
        ax.fill([0, self.B, self.B, 0, 0], [0, 0, self.H, self.H, 0], color='gray')
        ax.plot([0, self.B, self.B, 0, 0], [0, 0, self.H, self.H, 0], color='black')


class Circle(CrossSectionShape):
    R = tr.Float(100, param=True, minmax=(1, 1000), auto_set=False, enter_set=True)

    def get_cs_area(self):
        pass

    def get_b(self, z_positions_array):
        pass

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        pass

class TShape(CrossSectionShape):
    name = 'T-shape'

    H = tr.Float(250, input=True)
    B_f = tr.Float(250, input=True)
    B_w = tr.Float(100, input=True)
    H_w = tr.Float(100, input=True)

    def get_cs_area(self):
        pass

    get_b = tr.Property(tr.Callable, depends_on='+input')
    @tr.cached_property
    def _get_get_b(self):
        z_ = sp.Symbol('z')
        b_p = sp.Piecewise(
            (self.B_w, z_ < self.H_w),
            (self.B_f, True)
        )
        return sp.lambdify(z_, b_p, 'numpy')

class CustomShape(CrossSectionShape):
    # The width b can be a sympy expression describing a variable B along the height
    # b = tr.Any(100, param=True, minmax=(1, 1000), auto_set=False, enter_set=True)

    def get_cs_area(self):
        pass

    def get_b(self, z_positions_array):
        pass

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        pass


if __name__ == '__main__':
    cs = CrossSectionShape()

    tv = cs.trait_view('ipw_view')
    print(tv.content.content)