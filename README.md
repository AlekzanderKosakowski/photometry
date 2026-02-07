Contains pipelines for reducing and extracting photometry from CCD imaging data.

01_ccd_reduce.py
  Apply bias and flat correction to CCD image data.

01_reduce_zorro.py
  Apply bias and flat correction to CCD image data, modified for data with field of view smaller than the CCD image itself, such as Zorro on Gemini South

02_phot.py
  Perform simple flux extraction with a fixed-radius aperture. Handles image drifting through moderate box-size centroiding.

03_plot_both_zorro.py
  Plot dual-filter light curve data from the Zorro instrument.

update_header82.py
  Update the FITS header metadata for CCD images obtained with the 82-inch Otto Struve telescope at McDonald observatory. Includes new data for BJD_TDB timing, filter, RA, and Dec, read noise, and gain.
