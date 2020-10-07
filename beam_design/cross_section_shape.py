
import traits.api as tr
import traitsui.api as ui
from bmcs_utils.api import InteractiveModel, IPWInteract

class CrossSection(InteractiveModel):

    H = tr.Float(100, param=True, minmax=(1,1000),
                 auto_set=False, enter_set=True)
    B = tr.Float(100, param=True, minmax=(1,1000),
                 auto_set=False, enter_set=True)

    param_names = ['H', 'B']

if __name__ == '__main__':
    cs = CrossSection()

    tv = cs.trait_view('ipw_view')
    print(tv.content.content)