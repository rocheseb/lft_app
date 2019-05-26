# README #

Linefit is an instrument line shape-calculating program written by Frank Hase at the Institute for Meteorology und Climate Research (IMK) in Karlsruhe. 

	Frank Hase, Thomas Blumenstock, and Clare Paton-Walsh, Analysis of the Instrumental Line Shape of High-Resolution Fourier Transform IR Spectrometers with Gas Cell Measurements and New Retrieval Software, Applied Optics, Vol. 38, Issue 15, pp. 3417-3422 (1999)

This app can be used to run linefit 14.7 and display its outputs.

This app works with python 3.x and with bokeh > 1.1.0 installed.

	Bokeh: https://bokeh.pydata.org/en/latest/docs/installation.html

There are issues with the latest versions of Tornado, the apps works fine with Tornado 4.5.2

To use this app with linefit 14.5, you need to revert to the commit 2c7acbb with bokeh 0.12.10 and python 2.7

### Python ###

I suggest downloading python from https://www.anaconda.com/download/
Choose your operating system and get Python3.7

To install bokeh use the command (with windows you need to run the terminal as administator):

	conda install -c bokeh bokeh

To install a package:

	conda install PackageName

If the "conda" command does not find the package, you can also use pip:

	pip install PackageName

After installing anaconda, you should only need to install the "parse" and "re" packages

If you encounter error messages related to bokeh when running the app, you can try to revert to an earlier version of the package with:

	conda install bokeh=1.1.0

### How to use this app ###

This can read both OPUS and .dpt files.

Using .dpt files requires extra steps.

- Put the lft_app folder in the linefit/lft147/ directory
- lft147.exe (or lft147_ifort on linux) and the directories 'hit' and 'ergs' must be present in the linefit/lft147/ directory
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

- To run the app, navigate to the linefit/lft147/ directory in your terminal and use the command

	bokeh serve --show lft_app

The --show option will pop up the browser.

By default the spectrum itself will not be plotted in the browser and also not saved in the data dictionary as it can lead to very large files and more loading time.

The --args spec option will display and save the whole spectrum.

While the server is running, the app will be available in the browser at localhost:5006/lft_app

- Python dictionaries of the data are saved in lft_app/saved_sessions/
- PDF documents with the plots are saved in lft_app/pdf/

There are two example spectra from Eureka in lft_app/spectra/cut/ so the app can be tested directly.

### Rationg of spectra ###

Spectra should be ratioed to ~1 to be used with the linefit extended mode:

- HCl cells, done in linefit or in the code: 
	- if the TCCON 14.7 mode is used for HCl cells and linefit will take care of the ratioing
	- for other modes, a 2nd order polynomial is used to fit the spectrum without its lines and is then used for the ratioing

- HBr cells, done in the code:
	- background
	- the background file should be cut the same way as the spectrum if .dpt files are used, have the same file name but starting with 'ref_' (e.g. ref_180308_eu_HBr_180_e_0.dpt)
	- put the HBr background files in lft_app/spectra/background/
	- the spectra are ratioed with the background
	- the resulting ratioed spectrum is ratioed with its own average to normalize it to ~1
	- For dpt spectra I typically cut them between 2500-2700 cm-1

- N2O cells, done by you ! :
	- background, but different resolution from the spectrum
	- the rationg of spectrum with background shoul be done in OPUS
	- the resulting spectrum should be placed in lft_app/spectra/cut
	- it will be ratioed with its own average to normalize it to ~ 1
	- For dpt spectra I typically cut them between 2100-2300 cm-1

### User interface ###

##### Update dropdowns #####

This button allows you to update the list of files in the dropdown (e.g. if you added a file under spectra/cut after starting the app)

##### Previous session #####

Dropdown to select a .npy file from a previously saved session, it loads all the data from previous tests and allows you to add new runs.

##### Load session #####

Load the selected .npy file

##### Spectrum #####

Dropdown to select a spectrum file from the lft_app/spectra/cut/ directory

##### ILS model #####

Can be used to select the ILS model used by linefit.

If one of the TCCON modes is chosen with a HBr or N2O test, the Extended mode is used instead.

The ILS model will be indicated in test buttons and the legend with "m="

##### Regularisation factor #####

Textinput to quickly change the regularisation factor (same number for reg_mod and reg_phase)

When selected a HCl spectrum, this widget is disabled and its value set to 'TCCON'

If the HCl cell tests are processed with the TCCON mode of linefit 14.7, the regularisation factor is adjusted during the processing and will be output in /ergs/TCCON.dat

The regularisation factor will be indicated in test buttons and the legend with "r=", with the TCCON ILS model of linefit 14.7 this will be "r=T" instead

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

For HCl cell spectra this is done directly in linefit with the TCCON 14.7 mode, so it takes longer to process these spectra compared to the extended mode

The python dictionaries saved in lft_app/saved_sessions/ can be merged with a utility program lft_app/utils/merge_sessions.py

The merged file can then be loaded from the browser.

### DISCLAIMER ###

This app handles a limited number of warning/error messages that can be given by linefit.

If any other warning or error message is given by linefit, this app will hang, you should then run linefit from the terminal to figure out what the problem is.

The app may hang if there is any convergence problem.

There will be more detailed outputs in the terminal than in the browser, the actual linefit messages will be in the terminal.

If a significant spectral detuning is detected, the app will use it to update the input file and re-run linefit.

For HCl cells the spectral detuning is corrected directly by linefit for the TCCON 14.7 mode before the fitting.

### Contact ###

For this app:

	sebastien.roche@mail.utoronto.ca

For the linefit program:

	frank.hase@kit.edu