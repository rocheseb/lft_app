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

	For each cell the main code uses by default:

	window_dict['hcl'] for hcl
	window_dict['n2o'] for hcl
	window_dict['hbr'] for hcl

	The use of the custom window lists is not yet supported
	"""

	window_dict = {}

	# custom MIR microwindows for HCl
	window_dict['hcl_test'] = [
					'(2843.1,2844.1)', # 1
					#'(2596.8,2597.8)',
					#'(2625.2,2626.2)',
					'(2862.5,2863.5)', # 2
					'(2864.6,2865.6)', # 3
					'(2903.6,2904.6)', # 4
					'(2905.7,2906.7)', # 5
					'(2923.2,2924.2)', # 6
					'(2923.2,2924.2)', # 7
					'(2923.2,2924.2)', # 8
					'(2925.4,2926.4)', # 9
					'(2960.6,2961.6)', # 10
					#'(2962.8,2963.8)', # 11
					'(3011.6,3012.6)', # 12
					'(3013.9,3014.9)', # 13
					'(3027.3,3028.3)', # 14
					#'(3029.6,3030.6)', # 15
					]

	# NIR microwindows for HCl
	window_dict['hcl'] = [	
					'(5683.0,5684.0)', # 1 
					'(5687.1,5688.1)', # 2
					'(5701.5,5702.5)', # 3
					'(5705.6,5706.6)', # 4
					'(5718.7,5719.7)', # 5 
					'(5734.6,5735.6)', # 6 
					'(5738.8,5739.8)', # 7
					'(5749.3,5750.3)', # 8
					'(5753.5,5754.5)', # 9
					'(5762.7,5763.7)', # 10
					'(5766.9,5767.9)', # 11
					'(5774.8,5775.8)', # 12
					'(5779.0,5780.0)', # 13
					]
					
	# NIR wide windows for HCl
	window_dict['hcl_1mw'] = [
					(5675.00,5799.99),	# 1
					(5712.0,5782.0),	# 2
					]

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