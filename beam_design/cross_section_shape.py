import traits.api as tr
from bmcs_utils.api import InteractiveModel
import numpy as np
from bmcs_utils.models.interactive_window import View, Item


class ICrossSectionShape(tr.Interface):
    """This interface lists the functions need to be implemented by cross section classes."""

    def get_cs_area(self):
        """This function should return the area of the cross section."""

    def get_b(self, z_positions_array):
        """This function should return b values that correspond to the positions z_positions_array."""


class CrossSectionShape(InteractiveModel):
    """"This class describes the geometry of the cross section."""
    h = tr.Float(600, minmax=(1, 1000))
    b = tr.Float(250, minmax=(1, 1000)) # To be deleted later
    # cross_section_type = tr.Enum([tr.Instance(Rectangle), ])

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        pass


@tr.provides(ICrossSectionShape)
class Rectangle(CrossSectionShape):
    h = tr.Float(500)
    b = tr.Float(250)

    ipw_view = View(
        Item('h', minmax=(1, 3000), latex='h'),
        Item('b', minmax=(1, 500), latex='b')
    )

    def get_cs_area(self):
        return self.b * self.h

    def get_b(self, z_positions_array):
        return np.full_like(z_positions_array, self.b)

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        ax.axis([0, self.b, 0, self.h])
        ax.axis('equal')
        ax.fill([0, self.b, self.b, 0, 0], [0, 0, self.h, self.h, 0], color='gray')
        ax.plot([0, self.b, self.b, 0, 0], [0, 0, self.h, self.h, 0], color='black')


@tr.provides(ICrossSectionShape)
class Circle(CrossSectionShape):
    r = tr.Float(100, param=True, minmax=(1, 1000), auto_set=False, enter_set=True)

    def get_cs_area(self):
        pass

    def get_b(self, z_positions_array):
        pass

    def subplots(self, fig):
        return super().subplots(fig)

    def update_plot(self, ax):
        pass


@tr.provides(ICrossSectionShape)
class CustomShape(CrossSectionShape):
    h = tr.Float(100, param=True, minmax=(1, 1000), auto_set=False, enter_set=True)
    # The width b can be a sympy expression describing a variable B along the height
    b = tr.Any(100, param=True, minmax=(1, 1000), auto_set=False, enter_set=True)

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