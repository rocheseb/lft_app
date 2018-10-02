# README #

Linefit is an instrument line shape-calculating program written by Frank Hase at the Institute for Meteorology und Climate Research (IMK) in Karlsruhe. 

	Frank Hase, Thomas Blumenstock, and Clare Paton-Walsh, Analysis of the Instrumental Line Shape of High-Resolution Fourier Transform IR Spectrometers with Gas Cell Measurements and New Retrieval Software, Applied Optics, Vol. 38, Issue 15, pp. 3417-3422 (1999)

This app can be used to run linefit 14.5 and display its outputs.

This app requires python 2.7.x (not tested with python 3.x) with bokeh installed.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html

### Python ###

I suggest downloading python from https://www.anaconda.com/download/
Choose your operating system and get Python2.7

To install bokeh use the command (with windows you need to run the terminal as administator):

	conda install -c bokeh bokeh

To install a package:

	conda install PackageName

If the "conda" command does not find the package, you can also use pip:

	pip install PackageName

After installing anaconda, you should only need to install the "parse" and "re" packages

If you encounter error messages related to bokeh when running the app, you can try to revert to an earlier version of the package with:

	conda install bokeh=0.12.10

### How to use this app ###

This can read both OPUS and .dpt files.

Using .dpt files requires extra steps.

- Put the lft_app folder in the linefit/lft145/ directory
- lft145.exe (or lft145_ifort on linux) and the directories 'hit' and 'ergs' must be present in the linefit/lft145/ directory
- Spectrum file names need to follow this naming convention: YYMMDD_site_cell_MOPD_X_num.ext
	- YYMMDD year month day
	- site: two letter site abbreviation (e.g. eu for Eureka, oc for Lamont)
	- cell: one of 'hbr', 'n2o', 'hcl'
	- X: 'v' for vented instrument, 'e' for evacuated
	- MOPD: the maximum optical path difference in cm
	- num: an index number for the cell test (there might be more than one per day)
	- ext: either 'dpt' or a number
	
		e.g. 180308_eu_HCl_45_e_0.dpt or 180308_eu_HCl_45_e_0.0 
		For the first HCl cell test with 45 MOPD in an evacuated instrument at Eureka on March 8 2018
		For opus files the extension number does not matter

- For several tests in one day : 161122_eu_HCl_45_e_0.dpt, 161122_eu_HCl_45_e_1.dpt etc.
- Spectra must be placed in lft_app/spectra/cut
- For '.DPT' files:
	- .dpt (data point table) files can be generated in OPUS via the pop-up window generated with "Save as" 
	- they must have no headers and be cut e.g. between ~5200-5900 wavenumbers for HCl cells
	- In lft_app/spectra/cut/temp.dat write the spectrum filename, scanner temperature, and aperture size.
	
			spectrumfilename1,temperature1,apt_size1
			spectrumfilename2,temperature2,apt_size2
			etc.

	- in lft_app/lft_setup.py add the focal length of collimator of your instrument
- For opus files the parameters are read directly from the file
	
- In lft_app/lft_setup.py, add your cell information (follow the template)

- To run the app, navigate to the linefit/lft145/ directory in your terminal and use the command

	bokeh serve --show lft_app --args light

The --show option will pop up the browser.

By default the spectrum itself will be plotted in the browser and also saved in the data dictionary. This can lead to very large files and more loading time.

The --args light option will not display or save the whole spectrum.

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF documents with the plots are saved in lft_app/pdf/

There are two example spectra from Eureka in lft_app/spectra/cut/ so the app can be tested directly.

### Rationg of spectra ###

Spectra should be ratioed to ~1 to be used with the linefit extended mode:

- HCl cells, done in the code: 
	- no background
	- In the code I fit a 2nd order polynomial to the spectrum without the lines and use that to ratio the spectrum to normalize it to ~1 (seems more consistent than using a fixed numbers)

- HBr cells, done in the code:
	- background
	- the background file should be cut the same way as the spectrum if .dpt files are used, have the same file name but starting with 'ref_' (e.g. ref_180308_eu_HBr_180_e_0.dpt)
	- put the HBr background files in lft_app/spectra/background/
	- the spectra are ratioed with the background
	- the resulting ratioed spectrum is ratioed with its own average to normalize it to ~1

- N2O cells, done by you ! :
	- background, but different resolution from the spectrum
	- the rationg of spectrum with background shoul be done in OPUS
	- the resulting spectrum should be placed in lft_app/spectra/cut
	- it will be ratioed with its own average to normalize it to ~ 1

### User interface ###

##### Update dropdowns #####

This button allows you to update the list of files in the dropdown (e.g. if you added a file under spectra/cut after starting the app)

##### Previous session #####

Dropdown to select a .npy file from a previously saved session, it loads all the data from previous tests and allows you to add new runs.

##### Load session #####

Load the selected .npy file

##### Spectrum #####

Dropdown to select a spectrum file from the lft_app/spectra/cut/ directory

##### Regularisation factor #####

Textinput to quickly change the regularisation factor (same number for reg_mod and reg_phase)

##### Run linefit #####

Runs linefit !!

##### Save name #####

Textinput to specify a name for pdf and npy files that will be saved under lft_app/save/ and lft_app/pdf/, no need for an extension

##### Save #####

Button to save the current document in .npy file under lft_app/save/

Also saves static plots in a PDF (the last summary layout will break if there are a lot of tests)

##### Loop key #####

Textinput to provide a filename pattern and run linefit (with the specified regularisation factor) for all the spectra under lft_app/spectra/cut that match the pattern.

This accepts * notations 

e.g. \*_eu_\* will run linefit for all Eureka tests, 

\*_eu_HCl_\* will run all Eureka HCl tests

1803\*_eu_HCl_\* will run all Eureka HCl tests from march 2018

etc.

I would not advise using it for more than a few tens of spectra, the layout/browser may not handle it

Since it is not possible to have more than 20 contrasting colors, the colors will switch to a continuous viridis color palette if there are more spectra than that.

##### Test buttons #####

For each linefit run, a button will be added, to allow to switch between different runs, and a a red cross button will appear next to it to remove that test from the document (slow process as it needs to remake the whole document each time)

##### Legends #####

Legends are clickable to show / hide corresponding lines

##### Toolbars #####

Check out the different options of the toolbars ! (icons on the left of plots)

The save tool is clunky, it will save each plot in a separate png in your default download folder

##### ILS and fits tab #####

You can switch between the different microwindows with the 'MW' buttons

The average RMS from all the windows is also displayed at the bottom

##### Averaging kernels tab #####

You can click or select the scatter points to highlight the correspond averaging kernel row(s)

### Other info ###

N2O and HBr cell spectra are processed in a loop until the cell pressure converges; this usually take 2-3 linefit runs.

The python dictionaries saved in lft_app/saved_sessions/ can be merged with a utility program lft_app/utils/merge_sessions.py

The merged file can then be loaded from the browser.

### DISCLAIMER ###

If any warning or error message is given by linefit, this app will hang, you should then run linefit from the terminal to figure out what the problem is.

The app may hang if there is any convergence problem.

There will be more detailed outputs in the terminal than in the browser.

If a significant spectral detuning is detected, the app will use it to update the input file and re-run linefit.

### Contact ###

For this app:

	sebastien.roche@mail.utoronto.ca

For the linefit program:

	frank.hase@kit.edu