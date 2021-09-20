from bmcs_utils.api import Model, Int, Item, View, Float
from bmcs_beam.beam_config.system.system import System

class SimpleDistLoadSystem(System):

    name = 'SimpleDistLoadSystem'

    q = Float(-5000)

    tree = []

    ipw_view = View(
        *System.ipw_view.content,
        Item('q', latex='q \mathrm{[kN/m]}'), # kN/m = N/mm
    )
