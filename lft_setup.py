def site_data():
	"""
	Return a dictionary with information specific to each site and cell
	"""

	data = {}

	data['hcl'] = hcl_cells()
	data['n2o'] = n2o_cells()
	data['hbr'] = hbr_cells()

	# focal length of collimator (mm), add the one for your site with the appropriate two letter site abbreviation
	data['FOC'] = {
					'eu':418.0, # Eureka
					'et':418.0, # East Trout Lake
					'ra':418.0, # St Denis
					'lr':418.0, # new Lauder
					'll':418.0, # old Lauder
					}

	return data

def hcl_cells():
	"""
	Returns a dictionary of the HCl cells data
	"""

	cell = {}

	template = {
				'ID':None,
				'length':None,
				'batch':None,
				'owner':None,
				'text_on_cell_body':None,
				'effp_h35cl_296k':None,
				'effp_h37cl_296k':None,
				'calibration_run':None,
				'h35cl_column':None,
				'h37cl_column':None,
				'location':None,
				}
	"""
	# Use the format below to add a cell
	# site
	cell['site'] = {key:template[key] for key in template}
	cell['site']['location'] =
	cell['site']['owner'] = 
	cell['site']['ID'] =
	cell['site']['length'] =
	cell['site']['batch'] =
	cell['site']['text_on_cell_body'] =
	cell['site']['effp_h35cl_296k'] = 
	cell['site']['effp_h37cl_296k'] = 
	cell['site']['h35cl_column'] = 
	cell['site']['h37cl_column'] = 
	cell['site']['calibration_run'] = 
	"""
	
	# Eureka
	cell['eu'] = {key:template[key] for key in template}
	cell['eu']['location'] = 'Eureka'
	cell['eu']['owner'] = 'University of Toronto'
	cell['eu']['ID'] = 30
	cell['eu']['length'] = 100
	cell['eu']['batch'] = 'CT1303'
	cell['eu']['text_on_cell_body'] = "cell #28"
	cell['eu']['effp_h35cl_296k'] = 4.78
	cell['eu']['effp_h37cl_296k'] = 4.817
	cell['eu']['h35cl_column'] = 1.2909e+22
	cell['eu']['h37cl_column'] = 1.2836e+22
	cell['eu']['calibration_run'] = 1312

	# East Trout Lake
	cell['et'] = {key:template[key] for key in template}
	cell['et']['location'] = 'East Trout Lake'
	cell['et']['owner'] = 'University of Toronto'
	cell['et']['ID'] = 62
	cell['et']['length'] = 100
	cell['et']['batch'] = 'CT1510'
	cell['et']['text_on_cell_body'] = 'cell #56'
	cell['et']['effp_h35cl_296k'] = 4.673
	cell['et']['effp_h37cl_296k'] = 4.683
	cell['et']['h35cl_column'] = 1.3310e+22
	cell['et']['h37cl_column'] = 1.3260e+22
	cell['et']['calibration_run'] = 1510

	# St Denis
	cell['ra'] = {key:template[key] for key in template}
	cell['ra']['location'] = 'St. Denis'
	cell['ra']['owner'] = 'Royal Belgian Institute for Space Aeronomy'
	cell['ra']['ID'] = 4
	cell['ra']['length'] = 100
	cell['ra']['batch'] = 'CT1303'
	cell['ra']['text_on_cell_body'] = 'cell #4'
	cell['ra']['effp_h35cl_296k'] = 5.10
	cell['ra']['effp_h37cl_296k'] = 5.10
	cell['ra']['h35cl_column'] = 1.268e+22
	cell['ra']['h37cl_column'] = 1.268e+22
	cell['ra']['calibration_run'] = 1303

	# New Lauder
	cell['lr'] = {key:template[key] for key in template}
	cell['lr']['location'] = 'Lauder'
	cell['lr']['owner'] = 'NIWA'
	cell['lr']['ID'] = 54
	cell['lr']['length'] = 100
	cell['lr']['batch'] = 'CT1510'
	cell['lr']['text_on_cell_body'] = 'cell #48'
	cell['lr']['effp_h35cl_296k'] = 4.706
	cell['lr']['effp_h37cl_296k'] = 4.847
	cell['lr']['h35cl_column'] = 1.3442e+22
	cell['lr']['h37cl_column'] = 1.3419e+22
	cell['lr']['calibration_run'] = 1510

	# Old Lauder
	cell['ll'] = {key:template[key] for key in template}
	cell['ll']['location'] = 'Lauder'
	cell['ll']['owner'] = 'NIWA'
	cell['ll']['ID'] = 33
	cell['ll']['length'] = 100
	cell['ll']['batch'] = 'Special'
	cell['ll']['text_on_cell_body'] = 'Caltech'
	cell['ll']['effp_h35cl_296k'] = 4.793
	cell['ll']['effp_h37cl_296k'] = 4.759
	cell['ll']['h35cl_column'] = 1.3028e+22
	cell['ll']['h37cl_column'] = 1.2947e+22
	cell['ll']['calibration_run'] = 1312

	# Wollongong
	cell['wg'] = {key:template[key] for key in template}
	cell['wg']['location'] = 'Wollongong'
	cell['wg']['owner'] = 'UoW'
	cell['wg']['ID'] = 36
	cell['wg']['length'] = 100
	cell['wg']['batch'] = 'CT0805'
	cell['wg']['text_on_cell_body'] = 'cell #4'
	cell['wg']['effp_h35cl_296k'] = 4.816
	cell['wg']['effp_h37cl_296k'] = 4.818
	cell['wg']['h35cl_column'] = 1.2247e22
	cell['wg']['h37cl_column'] = 1.2164e22
	cell['wg']['calibration_run'] = 1308

	return cell

def hbr_cells():
	"""
	Returns a dictionary of the HBr cells data
	"""
	cell = {}

	template = {
				'ID':None,
				'length':None,
				'batch':None,
				'owner':None,
				'text_on_cell_body':None,
				'calibration_run':None,
				'location':None,
				'pressure':None,
				'column':None,
				}

	# Eureka
	cell['eu'] = dict(template)
	cell['eu']['location'] = 'Eureka'
	cell['eu']['pressure'] = 1.52
	cell['eu']['column'] = 7.60e+20

	# Toronto
	cell['to'] = dict(template)
	cell['to']['location'] = 'Toronto'
	cell['to']['pressure'] = 2.4
	cell['to']['column'] = 1.08e+21

	return cell

def n2o_cells():
	"""
	Returns a dictionary of the N2O cells data
	"""

	cell = {}

	template = {
				'ID':None,
				'length':None,
				'batch':None,
				'owner':None,
				'text_on_cell_body':None,
				'calibration_run':None,
				'location':None,
				'pressure':None,
				'column':None,
				}

	# Eureka
	cell['eu'] = dict(template)
	cell['eu']['location'] = 'Eureka'
	cell['eu']['pressure'] = 0.8734
	cell['eu']['column'] = 4.373e+20

	# Toronto
	cell['to'] = dict(template)
	cell['to']['location'] = 'Toronto'
	cell['to']['pressure'] = 0.8734
	cell['to']['column'] = 4.478e+20

	return cell

def window_data():
	"""
	List of windows to be used by linefit for each cell
	"""

	window_dict = {}

	# NIR window for HCl with linefit 14.7; instead of multiple microwindows there is one big region that will have deweighted regions
	window_dict['hcl'] = ['(5675.00,5800.20)',]

	# MIR microwindows for N2O
	window_dict['n2o'] = [
					'(2167.03,2185.25)', # 1
					'(2222.825,2223.019)', # 2	
					'(2224.457,2224.715)', # 3
					]

	# MIR microwindows for HBr
	window_dict['hbr'] = [
					'(2590.32, 2590.72)', # 1
					'(2590.71, 2591.11)', # 2
					'(2605.60, 2606.00)', # 3
					'(2606.00, 2606.40)', # 4
					'(2620.39, 2620.79)', # 5
					'(2620.80, 2621.20)', # 6
					'(2634.70, 2635.10)', # 7
					'(2635.10, 2635.50)', # 8
					'(2648.50, 2648.90)', # 9
					'(2648.90, 2649.30)', # 10
					'(2661.76, 2662.16)', # 11
					'(2662.18, 2662.58)', # 12
					'(2674.52, 2674.92)', # 13
					'(2674.94, 2675.34)', # 14
					]

	return window_dict