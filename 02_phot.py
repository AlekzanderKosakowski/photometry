import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from photutils.centroids import centroid_com, centroid_sources, centroid_2dg, centroid_com, centroid_quadratic
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from photutils.aperture import aperture_photometry
from photutils.utils import calc_total_error

# For background statistic estimation
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint

from astropy.io import fits

import os
import sys


# Disable warning messages. (This is a bad idea in the long term. Address the warnings and remove this later.)
import warnings
warnings.filterwarnings('ignore') 

def update_plot():
    '''
    Update the image plot with new markers.
    '''
    plt.cla()
    ax.imshow(data, origin='upper', vmin=np.percentile(data, 0.25), vmax=np.percentile(data, 99.75))
    if len(x_init)>0:
        ax.scatter(x_init[0], y_init[0], marker='*', s=3*80, color="red", facecolor="None")
        ax.scatter(x_init[1:], y_init[1:], marker='+', s=80, color="red")
    plt.draw()

def on_click(event):
    '''
    Save locations of left-mouse clicks as initial guesses for centroiding algorithm.
    Remove entries from the list with right-mouse clicks.
    '''
    if event.button == 3:
        #
        # Right-click: delete nearest data point.
        #
        distances = [ ((event.xdata - x_init[i])**2 + (event.ydata - y_init[i])**2)**0.5 for i in range(len(x_init)) ]

        if len(x_init) > 0:
            x_init.pop(np.argmin(distances))
            y_init.pop(np.argmin(distances))


    elif event.button == 1:
        #
        # Left-click: Add initial centroid guess at clicked location as long as no other point exists at that location within 5 px
        #

        # Auto-centroid at the clicked location
        x_autocenter, y_autocenter = get_centroid(event.xdata, event.ydata)
        x_autocenter = x_autocenter[0] if type(x_autocenter) in [list, np.ndarray] else x_autocenter
        y_autocenter = y_autocenter[0] if type(y_autocenter) in [list, np.ndarray] else y_autocenter

        distances = [ ((event.xdata - x_init[i])**2 + (event.ydata - y_init[i])**2)**0.5 for i in range(len(x_init)) ] # Distance from existing plotted regions to the clicked/centroided regions.

        if len(distances) > 0:
            if min(distances)>5:
                x_init.append(x_autocenter)
                y_init.append(y_autocenter)
        else:
            x_init.append(x_autocenter)
            y_init.append(y_autocenter)    

    update_plot()
    return

def get_bkg(data):
    '''
    https://photutils.readthedocs.io/en/stable/background.html#background
    Currently uses the simple constant background assumption.
    However, bright stars on/near the CCD will create a gradient background which can be dealt with using the same tools.
    '''

    # Copy/Paste the code from photutils URL to find and remove sources from the data image, then perform 3-sigma clipping of cosmic rays, then calculate basic statistics.
    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
    threshold = detect_threshold(data, nsigma=2.0, sigma_clip=sigma_clip)
    segment_img = detect_sources(data, threshold, npixels=10)
    footprint = circular_footprint(radius=10)
    mask = segment_img.make_source_mask(footprint=footprint)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)

    # Create a new image with random Gaussian noise based on the estimated background noise parameters from above. We feed this into the photometry later to estimate errorbars.
    bkg_image = np.random.normal(mean, std, data.shape)

    if False:
        # Plot the image
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111)
        ax.imshow(data, origin='upper', vmin=np.percentile(data, 0.25), vmax=np.percentile(data, 99.75))
        plt.show()

        # Plot the image with the source regions cut out
        data2 = data*~mask
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111)
        ax.imshow(data2, origin='upper', vmin=np.percentile(data, 0.25), vmax=np.percentile(data, 99.75))
        plt.show()

        # Plot the estimated background image
        plt.imshow(bkg_image)
        plt.show()

    # Return an image with constant pixel value based on the mean of the Gaussian noise across the data image.
    # We return a constant array to allow for pixel masking if necessary, which would result in a non-constant array.
    return(mean*np.ones_like(data), std*np.ones_like(data))


def show_centroid(xcenter, ycenter):
    '''
    Create an imshow plot displaying the CCD image with markers at locations of current centroids.
    '''
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    ax.imshow(data, origin='upper', vmin=np.percentile(data, 0.25), vmax=np.percentile(data, 99.75))
    ax.scatter(xcenter[0], ycenter[0], marker='*', s=3*80, color="red", label="Centroids", facecolor="None")
    ax.scatter(xcenter[1:], ycenter[1:], marker='x', s=80, color="red", label="Centroids")
    plt.tight_layout()
    plt.show()

def get_centroid(x_init, y_init, box=11):
    #
    # Centroid the clicked locations and plot the best centers.
    #
    x, y = centroid_sources(data, x_init, y_init, box_size=box, centroid_func=centroid_2dg)
    return(x, y)

if __name__ == "__main__":
    #
    #  
    #
    global data
    global ax

    # Create list of image files
    path = "./"
    files = [path+k for k in os.listdir(path) if k[-5:]==".fits" and "._" not in k]

    # Obtain CCD parameters
    filename = np.random.choice(sorted(files[0:10]))
    file = fits.open(filename)
    try:
        data = file[0].data
        assert len(data.shape) == 2
    except AssertionError:
        data = file[0].data[0]
        assert len(data.shape) == 2, "Data must be 2-dimensional"

    header = file[0].header
    gain = header.get("gain", 1.0) # Use gain=1.0 if header value is not found.

    # Display the data in the first image. Click the image to select stars for photometry.
    x_init = []
    y_init = []

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.imshow(data, origin='upper', vmin=np.percentile(data, 0.25), vmax=np.percentile(data, 99.75))
    plt.connect('button_press_event', on_click)
    plt.tight_layout()
    plt.show()

    xcenter, ycenter = get_centroid(x_init, y_init, 11)
    # show_centroid(xcenter, ycenter)

    df = pd.DataFrame({})

    n_centroid = 30 # Re-centroid every n_centroid images
    for i,filename in enumerate(sorted(files)):
        #
        # Use astropy for photometry of each source.
        # Aperture size fixed for now, but later based on average PSF
        #
        if i%10 == 0:
            print(filename)
        aperture_r = 12.0

        # Read in each data file individually
        file = fits.open(filename)

        data = file[0].data
        if data.shape[0] == 1: # TODO: this line will break if data shape is something like (N>1, 256, 256), with more than one image saved per HDU
            data = data[0]
        
        if i % n_centroid == 0:
            try:
                # Re-center apertures based on previous centroid locations.
                xcenter, ycenter = get_centroid(xcenter, ycenter, 11)
                centers = [(k,l) for k,l in zip(xcenter, ycenter)]
            except ValueError:
                continue

        # Create circular source apertures located at each centroid.
        src_apertures = CircularAperture([k for k in centers], r=aperture_r)
        src_stats = ApertureStats(data, src_apertures)

        # Create annular background apertures surrounding source apertures to estimate local background.
        bkg_apertures = CircularAnnulus([k for k in centers], r_in=1.5*aperture_r, r_out=2.5*aperture_r)
        bkg_stats = ApertureStats(data, bkg_apertures, sigma_clip=SigmaClip(sigma=3.0, maxiters=10)) # 10 iterations of 3-sigma clipping to remove potential cosmic ray hits.

        iraf = False # Use equation from IRAF phot documentation to estimate flux uncertainty. This is nearly identical to the photutils background estimation method.
        if iraf:
            # Estimate source flux within each aperture.
            src_flux = src_stats.sum - src_apertures.area*bkg_stats.median
            src_flux_err = (src_flux / gain + src_apertures.area * bkg_stats.std**2 + src_apertures.area**2 * bkg_stats.std**2 / bkg_apertures.area)**0.5 # "error" equation taken from IRAF "phot" task (https://iraf.readthedocs.io/en/latest/tasks/noao/digiphot/apphot/phot.html)

        try:
            bkg_image, bkg_std = get_bkg(data)
            error = calc_total_error(data, bkg_std, gain)
        except AttributeError:
            print(f"File failed background estimation: {filename}. Skipping file.")
            continue

        phot_table = aperture_photometry(data-bkg_image, src_apertures, error=error)
        phot_table["bjd_tdb"] = file[0].header["BJD_TDB"]
        df = pd.concat((df, phot_table.to_pandas()))

    # Scale each of the light curves to their median values.
    # This ensures that we divide by the ratio of the dips caused by clouds/etc when applying the calibration light curve later.
    for i in np.unique(df["id"].values):
        df.loc[df["id"]==i, "aperture_sum_err"] /= np.nanmedian(df[df["id"]==i]["aperture_sum"])
        df.loc[df["id"]==i, "aperture_sum"] /= np.nanmedian(df[df["id"]==i]["aperture_sum"])

    # Create a weighted-mean light curve using the light curves of the field stars
    fluxes = np.array([df[df["id"]==i]["aperture_sum"].values for i in np.unique(df["id"].values)[1:]])
    weights = np.array([(1/df[df["id"]==i]["aperture_sum_err"].values)**2 for i in np.unique(df["id"].values)[1:]])

    # Weighted average calibration light curve using all stars except the first selected star.
    weighted, weighted_error = np.average(fluxes, axis=0, weights=weights, returned=True)
    weighted_error = weighted_error**-0.5

    time = df[df["id"]==1]["bjd_tdb"].values
    df.insert(6, "rel_flux", 1.0)
    df.insert(7, "rel_flux_err", 1.0)

    flux1 = df[df["id"]==1]["aperture_sum"]
    flux1_err = df[df["id"]==1]["aperture_sum_err"]

    for i in np.unique(df["id"].values):
        df.loc[df["id"]==i, "rel_flux"] = df.loc[df["id"]==i, "aperture_sum"]/weighted
        df.loc[df["id"]==i, "rel_flux_err"] = ((df.loc[df["id"]==i, "aperture_sum_err"]/weighted)**2 + (df.loc[df["id"]==i, "aperture_sum"]*weighted_error/weighted**2)**2)**0.5

    # Save the data to a CSV file. Keep only the best quality data based on the errorbars.

    df_clean = df[(df["id"]==1) &
                    (df["rel_flux_err"] < np.percentile(df[df["id"]==1]["rel_flux_err"],95)) &
                    (df["rel_flux_err"]<0.075)][["bjd_tdb","rel_flux","rel_flux_err"]]

    ofilename = "output_lc.csv"
    df_clean.to_csv(ofilename,index=False)

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    ax2 = ax.secondary_xaxis('top',functions=(lambda x: 24*(x-min(df_clean['bjd_tdb'].values)), lambda x: x/24+min(df_clean['bjd_tdb'].values)))
    ax2.set_xlabel('Time (hours)', fontsize='large')

    ax.errorbar(df_clean["bjd_tdb"], df_clean["rel_flux"], yerr=df_clean["rel_flux_err"], color='black', marker=".", markersize=3, elinewidth=0.5, capsize=1, linestyle='None', alpha=0.8)
    ax.set_xlabel("BJD_TDB (days)", fontsize="large")
    ax.set_ylabel("Relative Flux", fontsize="large")
    plt.grid(visible=True,which='both',linestyle="--", alpha=0.5, color='grey', linewidth=0.5)
    plt.savefig("output_lc.jpg", dpi=100)
    plt.show()


    # TODO: Airmass correction using a similar-colored non-variable star in the field as the main object of interest.
    #        This may require astroquery to check colors of nearby stars, or just have the user select the star(s) to use.


    # Fit and remove the light curve trend due to changing airmass.
    # mintime = np.min(time)
    # time -= mintime # PolyFit does not perform well with large numbers. JD are too large.
    # time2 = np.copy(time)
    # for k in range(5):
    #     std = np.std(weighted)
    #     med = np.median(weighted)
    #     ind = np.where(abs(weighted - med) < 3*std)[0]
    #     weighted = weighted[ind]
    #     weighted_error = weighted_error[ind]
    #     time2 = time2[ind]

    # order = 7
    # # Fit a polynomial to this light curve.
    # polfit = np.polyfit(x=time2, y=weighted, w=1/weighted_error, deg=order, rcond=None, full=False)
    # fit = np.zeros_like(time)
    # fit2 = np.zeros_like(time2)
    # for k in range(0, order+1):
    #     fit = fit + np.multiply(np.power(time, (order-k)), polfit[k])
    #     fit2 = fit2 + np.multiply(np.power(time2, (order-k)), polfit[k])

    # plt.errorbar(time2, weighted, yerr=weighted_error, color='black', marker=".", markersize=5, elinewidth=1, capsize=2, linestyle='None')
    # plt.plot(time2, fit2, color='red')
    # plt.show()

    # time += mintime
    # flux = df[df["id"]==1]["aperture_sum"]/fit
    # flux_err = df[df["id"]==1]["aperture_sum_err"]/fit

    # med = np.median(flux)
    # flux /= med
    # flux_err /= med

    # plt.errorbar(df[df["id"]==1]["bjd_tdb"]+mintime, df[df["id"]==1]["aperture_sum"]/fit, yerr=((df[df["id"]==1]["aperture_sum_err"]/fit)**2)**0.5, color='black', marker=".", markersize=5, elinewidth=1, capsize=2, linestyle='None')
    # plt.show()


    # plt.errorbar(time, flux, yerr=flux_err, color='black', marker=".", markersize=5, elinewidth=1, capsize=2, linestyle='None')
    # plt.show()

    # sys.exit()

        
