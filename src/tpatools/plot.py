import numpy as np
import matplotlib.pyplot as plt
from tpatools.tools import eV_to_nm

def voight(wavelength, intensity, width=25, shape=0.5, yscale=None, rng=None):
    if rng is None:
        l = np.min(wavelength) - 0.25*(np.max(wavelength) - np.min(wavelength))
        r = np.min(wavelength) + 0.25*(np.max(wavelength) - np.min(wavelength))
    else:
        l = rng[0]
        r = rng[1]

    #n = int(r - l + 1)
    x = np.linspace(l,r,5000)

    yvals=[]
    for i, wav in enumerate(wavelength):
        yvals.append(
             intensity[i] * (
            (
                shape * 0.832555/(width*1.772454)
                * np.exp(-2.772589*((x-wav)/width)**2)
                + (1-shape)/(3.1415927*width * (1+4*((x-wav)/width)**2) )
            ) /
            (
                shape*0.832555/(width*1.772454)
                + (1-shape)/(3.1415927*width)
            )
        )
        )
        y = sum(yvals)
        if yscale is None:
            y = y / np.max(y)
        else:
            y = y / yscale
    
    return x, y

def lorentzian(wavelength, intensity, width=25, rng=None, yscale=None):
    if rng is None:
        l = np.min(wavelength) - 0.25*(np.max(wavelength) - np.min(wavelength))
        r = np.min(wavelength) + 0.25*(np.max(wavelength) - np.min(wavelength))
    else:
        l = rng[0]
        r = rng[1]

    x = np.linspace(l,r,5000)
    yvals = []

    for i, wav in enumerate(wavelength):
        sq_term = (x - wav) / (width / 2)
        yvals.append(
            intensity[i] / (1 + sq_term**2)
        )
   
    y = sum(yvals)
    if yscale is None:
        y = y / np.max(y)
        y = y / yscale

    return x, y

def tpabroaden(energy, cross_section, width=0.1, rng=(3,5)):
    """
        Line-broadening for simulation of a 2PA spectrum
        Provide lifetime broadening parameter in eV (default 0.1)
    """
    x = np.linspace(rng[0],rng[1],5000)
    yvals = []
    cross_removelineshape = np.array(cross_section) * (np.pi * width)

    for i, en in enumerate(energy):
        yvals.append(
            cross_removelineshape[i] * width / (
                np.pi * (x - en)**2 + np.pi*width**2
            )
        )
    y = sum(yvals)
    return x, y


def tpaplot(
        tab, 
        width=0.1, 
        xmin=None, 
        xmax=None, 
        colour = None,
        nm=False, 
        save=None, 
        figure_size=(6,6),
        extraplotparams={},
        default_buffer_x_edge = 1,
    ):
    """
    Function to plot a simulated 2PA spectrum, given a table containing excitation energies and cross sections, as output by the tpaplot.parse.escf_table function
    """

    fig, ax = plt.subplots(figsize=figure_size)

    cross_section = tab['Cross Section /GM'].values
    ex_energy = tab['Excitation Energy /eV'].values

    if xmin is None:
        xmin = np.min(ex_energy) - default_buffer_x_edge

    if xmax is None:
        xmax = np.max(ex_energy) + default_buffer_x_edge

    rng = (xmin, xmax)



    x, y = tpabroaden(
        ex_energy, 
        cross_section,
        rng=rng,
        width=width,
    )

    if nm:
        x = eV_to_nm(x) * 2


    if colour is not None:
        ax.plot(
            x,
            y,
            color=colour,
            **extraplotparams,
        )
    else:
        ax.plot(
            x,
            y,
            color=colour,
            **extraplotparams,
        )

    ax.set_ylabel('$\\sigma^{2PA}$ /GM')
    if nm:
        ax.set_xlabel('$\\lambda$ /nm')
    else:
        ax.set_xlabel('Energy /eV')
    
    fig.tight_layout()

    if save is None:    
        plt.show()
    else:
        plt.savefig(save)
        plt.show()



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


    ax.set_ylabel('$\\sigma^{2PA}$ /GM')
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

