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


