import numpy as np
from tpatools.tools import eV_to_nm
import matplotlib.pyplot as plt
from tpatools.plot import tpabroaden

def tpaplot(
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
                    (np.min(x), y[0] + 20),
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
