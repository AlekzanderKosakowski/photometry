Contains pipelines for reducing and extracting photometry from CCD imaging data.

01_ccd_reduce.py
<ul><li>Apply bias and flat correction to CCD image data.</li></ul>

01_reduce_zorro.py
  <ul><li>Apply bias and flat correction to CCD image data, modified for data with field of view smaller than the CCD image itself, such as Zorro on Gemini South</li></ul>

02_phot.py
  <ul><li>Perform simple flux extraction with a fixed-radius aperture.</li><li>Handles image drifting through moderate box-size centroiding.</li></ul>

03_plot_both_zorro.py
  <ul><li>Plot the extracted dual-filter light curve data from the Zorro instrument.</li></ul>

update_header82.py
  <ul><li>Update the FITS header metadata for CCD images obtained with the 82-inch Otto Struve telescope at McDonald observatory.</li><li>Includes new data for BJD_TDB timing, filter, RA, and Dec, read noise, and gain.
</li></ul>
