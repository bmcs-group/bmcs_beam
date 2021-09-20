from bmcs_utils.api import Model, Int, Item, View, Float
from bmcs_beam.beam_config.system.system import System

class ThreePBSystem(System):

    name = 'ThreePBSystem'

    F = Float(-1000)

    tree = []

    ipw_view = View(
        *System.ipw_view.content,
        Item('F', latex='F \mathrm{[N]}'),
    )
