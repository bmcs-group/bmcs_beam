from bmcs_utils.api import Model, Int, Item, View, Float


class System(Model):

    name = 'System'

    L = Float(3000)
    n_x = Int(100)

    tree = []

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}'),
        Item('n_x', latex='n_x'),
    )
