import numpy as np
from tpatools.tools import eV_to_nm
import matplotlib.pyplot as plt
from tpatools.plot import tpabroaden
import ipywidgets as widgets

def tpaplot_multi(
        tabledict, 
        width, 
        x_offset, 
        y_offset, 
        xmin, 
        xmax, 
        colours = None,
        fromentry=0, 
        toentry=1000, 
        show_y=False, 
        show_labels=True, 
        nm=False, 
        justone=False, 
        showlegend=False,
        save=None, 
        monocolour=None, 
        figure_size=(6,6),
        extraplotparams={},
        label_offset_y=20,
        label_offset_x=0,
    ):

    fig, ax = plt.subplots(figsize=figure_size)
    rng = (xmin, xmax)

    if toentry is None:
        toentry = len(tabledict)

    def condition_statement(i, fromentry, toentry, justone):
        if justone:
            return i == (fromentry - 1)
        else:
            return i in range(fromentry-1,toentry)


    for i, (entryno, tab) in enumerate(tabledict.items()):
        #if i in range(fromentry-1, toentry):
        if condition_statement(i, fromentry, toentry, justone):
            cross_section = tab['Cross Section /GM'].values
            ex_energy = tab['Excitation Energy /eV'].values


            x, y = tpabroaden(
                ex_energy, 
                cross_section,
                rng=rng,
                width=width,
            )
            #print(i - fromentry + 1)
            y = y + y_offset * (i - fromentry + 1)
            x = x + x_offset * (i - fromentry + 1)

            if nm:
                x = eV_to_nm(x) * 2

            if monocolour is None and colours is None:
                current_colour = None
            elif monocolour is None and colours is not None:
                current_colour = colours[i]
            else:
                current_colour = monocolour

            if current_colour is not None:
                ax.plot(
                    x,
                    y,
                    color=current_colour,
                    label=entryno,
                    **extraplotparams,
                )
            else:
                ax.plot(
                    x,
                    y,
                    color=current_colour,
                    label=entryno,
                    **extraplotparams,
                )

            if show_labels:
                ax.annotate(
                    entryno,
                    (np.min(x) + label_offset_x, y[0] + label_offset_y),
                    color=current_colour
                )


    ax.set_ylabel('$\\sigma^{2PA}$')
    if nm:
        ax.set_xlabel('$\\lambda$ /nm')
    else:
        ax.set_xlabel('Energy /eV')
    
    if show_y == False:
        ax.set_yticks([])
    #ax.set_ylim(0,10)
    fig.tight_layout()
    if showlegend:
        ax.legend()

    if save is None:    
        plt.show()
    else:
        plt.savefig(save)
        plt.show()


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

