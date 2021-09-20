from bmcs_utils.api import Model, Int, Item, View, Float
from bmcs_beam.beam_config.system.system import System

class FourPBSystem(System):

    name = 'FourPBSystem'

    F = Float(-1000)
    L_F = Float(1000, desc='length from support to the force')

    tree = []

    ipw_view = View(
        *System.ipw_view.content,
        Item('F', latex='F \mathrm{[N]}'),
        Item('L_F', latex='L_F \mathrm{[N]}'),
    )
