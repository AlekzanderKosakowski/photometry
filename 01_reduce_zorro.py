import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from datetime import datetime

def rebin_image(data, bin=2):
    '''
    Bin data down to a fraction of its original size

    data: 2D NumPy array of CCD image
    bin:  factor by which to bin

    '''
    bdata_rows = []
    for row in data:
        bdata_rows.append([k.sum() for k in row.reshape(-1,new_bin)])
    bdata_rows = np.array(bdata_rows, dtype=np.float32)

    bdata = []
    for col in bdata_rows.T:
      bdata.append([k.sum() for k in col.reshape(-1,new_bin)])
    bdata = np.array(bdata, dtype=np.float32).T

    return(bdata)
    
def plot_image(image, iqr):
    '''
    Plot the image using the inner-quantile range as the [vmin,vmax] values.
    '''    
    min = np.nanpercentile(image, 50-iqr/2)
    max = np.nanpercentile(image, 50+iqr/2)
    plt.imshow(image, vmin=min, vmax=max, cmap="binary_r", origin="lower")
    plt.show()


def imcombine(files, method="median"):
    '''
    Combine images (median or sum)
    '''
    data = np.array([fits.open(k)[0].data for k in files]) # Collection of data matrices

    xdim = len(data[0][0]) # x-dimension length (pixels)
    ydim = len(data[0])    # y-dimension length (pixels)

    if method=="median":
        master = np.array([np.median([k[i][j] for k in data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)
    elif method=="sum":
        master = np.array([np.sum([k[i][j] for k in data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)

    return(master.astype(np.float64))


if __name__ == "__main__":
    #
    # Modified "01_ccd_reduce.py" for data from Zorro (and Alopeke) on Gemini South (North)
    #
    color = sys.argv[1] # r/b


    # Bias file preparation.
    print("Bias")
    bias_files = ["./bias/"+k for k in os.listdir("./bias/") if f"{color}.fits" in k]
    master_bias = imcombine(bias_files, "median")

    xdim = fits.open(bias_files[0])[0].header["NAXIS1"]
    ydim = fits.open(bias_files[0])[0].header["NAXIS2"]

    # Flats
    print("Flat")
    flat_files = ["./flat/"+k for k in os.listdir("./flat/") if f"{color}.fits" in k]
    flat_data0 = [fits.open(filename)[0].data for filename in flat_files]

    flat_data = [np.array(k, dtype=np.float32) for k in flat_data0 if np.max(k)<65000 and np.median(k)>3*np.median(master_bias)] # Ignore saturated flats and skip the bias frames taken with flats.
    
    # Build a mask which contains only the in-field-of-view pixels using the flats. One mask per flat, then combine where all masks are TRUE for the master mask
    # Replace all values outside of the field of view with np.nans so we can ignore them when calculating medians
    masks = []
    for flat in flat_data:
      masks.append(flat < 5*np.median(master_bias))
    masks = np.array(masks)

    master_mask = np.full(flat_data0[0].shape, True) 
    for i in range(len(master_mask)):
      for j in range(len(master_mask[i])):
        master_mask[i][j] = all(masks[:,i,j])

    # Apply the master_mask to each flat, replacing out of view pixels with nan
    for flat in flat_data:
      flat[master_mask] = np.nan
      flat -= master_bias

    master_flat = np.array([np.nanmedian([k[i][j] for k in flat_data]) for i in range(xdim) for j in range(ydim)]).reshape(xdim, ydim)
    master_flat /= np.nanmedian(master_flat)
    # plot_image(master_flat, 95)

    # Create another mask which includes only response>% within the master flat.
    master_mask2 = master_flat < 0.20
    master_flat[master_mask2] = np.nan
    # plot_image(master_flat, 95)



    # Science
    print("Science")
    sci_files = [f"./science/"+k for k in os.listdir(f"./science/") if k[-6:]==f"{color}.fits" and "._" not in k]

    for i, filename in enumerate(sorted(sci_files)):
        if (i+1)%100==0:
            print(f"  Progress: {(i+1)/len(sci_files)*100:>6.2f}%  {str(datetime.now().strftime('%H:%M:%S'))}")
        with fits.open(filename, mode="update") as ifile:
            if ifile[0].header.get("ccdred",False):
                print(f"Skipping {filename}; file already modified.")
                continue

            data = ifile[0].data.astype(np.float32) # Read in the data
            data[master_mask + master_mask2] = np.nan # Apply the masks

            data -= master_bias # Subtract bias
            data /= master_flat # Divide flat

            ifile[0].data = data # Replace data
            ifile[0].header["ccdred"] = str(datetime.now()) # Update the header to show the time of change
            print(f"  Updated {filename}")
