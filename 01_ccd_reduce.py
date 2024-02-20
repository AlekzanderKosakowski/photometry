from astropy.io import fits
import numpy as np
import os
# from numba import njit
from datetime import datetime

def imcombine(files, method="median"):
    '''
    Combine images (median or sum)
    '''
    data = np.array([fits.open(k)[0].data[0] for k in files]) # Collection of data matrices
    xdim = len(data[0][0]) # x-dimension length (pixels)
    ydim = len(data[0])    # y-dimension length (pixels)

    if method=="median":
        master = np.array([np.median([k[i][j] for k in data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)
    elif method=="sum":
        master = np.array([np.sum([k[i][j] for k in data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)

    return(master)

# @njit
def sub_bias(master_bias, fits_image):
    '''
    Subtract master bias values from each of the files in a list
    '''
    return(fits_image - master_bias)


if __name__ == "__main__":
    #
    # Basic CCD reduction pipeline for CCD image data from the McDonald Observatory 2.1-meter telescope.
    # TODO:
    #    Generalize by adjusting header values (EXPTIME is currently assumed to be milliseconds)
    #    Require user-provided directory or file names for bias, dark, and flat images. argparse would work well here.
    #    Require user-provided directory or file names for science images. Handle separately from calibration images.
    #
    # Bias file preparation.
    print("Bias")
    bias_files = ["./bias/"+k for k in os.listdir("./bias/") if "bias" in k and ".fits" in k]
    master_bias = imcombine(bias_files, "median")

    xdim = fits.open(bias_files[0])[0].header["NAXIS1"]
    ydim = fits.open(bias_files[0])[0].header["NAXIS2"]

    # Dark file preparation.
    print("Dark")
    dark_files = ["./dark/"+k for k in os.listdir("./dark/") if "dark" in k and ".fits" in k]
    darks = {str(k):[] for k in np.unique([round(float(fits.open(filename)[0].header["EXPTIME"])/1000) for filename in dark_files])}
    for exptime in darks:
        dark_files = ["./dark/"+k for k in os.listdir("./dark/") if f"dark{exptime}-" in k and ".fits" in k]
        dark_data = [fits.open(filename)[0].data[0] - master_bias for filename in dark_files]
        darks[exptime] = np.array([np.median([k[i][j] for k in dark_data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)

    # DomeFlats
    print("Flat")
    flat_files = ["./flat/"+k for k in os.listdir("./flat/") if "domeflat" in k and ".fits" in k]
    flats = {k:[] for k in np.unique([j.split("-")[0].split("_")[-1] for j in flat_files])}
    for filter in flats:
        flat_files = ["./flat/"+k for k in os.listdir("./flat/") if f"domeflat_{filter}" in k and ".fits" in k]
        exptime = str(round(float(fits.open(flat_files[0])[0].header["EXPTIME"])/1000))
        flat_data = [fits.open(filename)[0].data[0] - master_bias - darks[exptime] for filename in flat_files]
        # flat_data = [k/np.median(k) for k in flat_data]
        flats[filter] = np.array([np.median([k[i][j] for k in flat_data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)
        flats[filter] /= np.median(flats[filter])

    # Science
    print("Science")
    objects = [k for k in os.listdir() if "fits" not in k and k not in ["dark", "bias", "flat", ".DS_Store"]]
    for object in objects:
        print(object)

        sci_files = [f"./{object}/"+k for k in os.listdir(f"./{object}/") if f"{object}" in k and k[-5:]==".fits"]

        filter = fits.open(sci_files[0])[0].header.get("filter",sci_files[0].split("-")[0].split("_")[-1])
        exptime = round(float(fits.open(sci_files[0])[0].header["exptime"]))

        master_dark = darks[str(round(exptime/1000))]
        master_flat = flats[filter]

        for i, filename in enumerate(sorted(sci_files)):
            with fits.open(filename, mode="update") as ifile:
                if ifile[0].header.get("ccdred",False):
                    print(f"Skipping {filename}; file already modified.")
                    continue
                print(f"Updating {filename}.")
                ifile[0].header["ccdred"] = str(datetime.now())

                data0 = ifile[0].data[0]

                data = data0 - master_bias
                data -= master_dark
                data /= master_flat

                ifile[0].data = data

