import numpy as np
import ipywidgets as widgets
from tpatools.plot import tpaplot_multi

def widgetplot(
        tabledict, 
        extraplotparams={}, 
        colours=None,
        widthslide=[0.001, 0.5, 0.001, 0.1],
        xoffslide=[-1,1,0.01,0],
        yoffslide=[0,1000,10,100],
        xminslide=[0,5,0.1,3.5],
        xmaxslide=[0.1,5,0.1,4.8],
        label_offset_y = 20,
        label_offset_x = 0,
        monocolour = None,
    ):
    """
    Create a jupyter notebook widget to interactively display tpa plots from
    a dictionary of table values

    :param dict extraplotparams: dictionary describing extra parameters supplied to the plot function (linewidth, etc.)
    :param colours: colours used for the plotting
    :type colours: list or None
    :param list widthslide: list to determine linewidth slider range - format [min, max, step, value]
    :param list xoffslide: x axis offset (in eV) - format [min, max, step, value]
    :param list yoffslide: y axis offset (in GM) - format [min, max, step, value]
    """
    widthwidget = widgets.FloatSlider(
        min=widthslide[0],
        max=widthslide[1], 
        step=widthslide[2],
        value=widthslide[3],
        description="Broadening:"
    )
    x_offsetwidget = widgets.FloatSlider(
        min=xoffslide[0],
        max=xoffslide[1], 
        step=xoffslide[2],
        value=xoffslide[3],
        description="x offset:"
    )
    y_offsetwidget = widgets.FloatSlider(
        min=yoffslide[0],
        max=yoffslide[1], 
        step=yoffslide[2],
        value=yoffslide[3],
        description="y offset:"
    )
    xminslider = widgets.FloatSlider(
        min=xminslide[0],
        max=xminslide[1], 
        step=xminslide[2],
        value=xminslide[3],
    )
    xmaxslider = widgets.FloatSlider(
        min=xmaxslide[0],
        max=xmaxslide[1], 
        step=xmaxslide[2],
        value=xmaxslide[3],
    )
    showybox = widgets.Checkbox(value=False, description='Show y axis')
    showlabelsbox=widgets.Checkbox(value=True, description='Show Labels')

    if len(tabledict) > 1:
        fromentryslider = widgets.IntSlider(min=1,max=len(tabledict), value=1)
        toentryslider = widgets.IntSlider(min=1,max=len(tabledict),value=len(tabledict))
    else:
        fromentryslider = widgets.fixed(1)
        toentryslider = widgets.fixed(1)

    widgets.interact(
        tpaplot_multi, 
        tabledict=widgets.fixed(tabledict),
        width=widthwidget, 
        #width = widgets.fixed(0.1),
        x_offset = x_offsetwidget,
        y_offset=y_offsetwidget,
        xmin=xminslider,
        xmax=xmaxslider,
        colours = widgets.fixed(colours),
        fromentry=fromentryslider,
        toentry=toentryslider,
        show_y = showybox,
        show_labels = showlabelsbox,
        save=widgets.fixed(None),
        monocolour=widgets.fixed(monocolour),
        figure_size=widgets.fixed((6,6)),
        extraplotparams=widgets.fixed(extraplotparams),
        label_offset_y = widgets.fixed(label_offset_y),
        label_offset_x = widgets.fixed(label_offset_x),
    )

def showtab(entryno, tabledict, roundto=3):
    print(f'{list(tabledict.keys())[entryno -1]}:\n')
    df = tabledict[list(tabledict.keys())[entryno - 1]]
    print(
        df.to_string(
            index = False,
            float_format=f'%.{roundto}f',
        )
    )

def tpatabs(tabledict, roundto=3):
    entryslider = widgets.IntSlider(min=1,max=len(tabledict))
    widgets.interact(
        showtab, 
        entryno=entryslider, 
        tabledict=widgets.fixed(tabledict),
        roundto=widgets.fixed(roundto),
    )

