import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.table import Table,vstack
import matplotlib.font_manager as mfm

#Read Star Catalog
Bright = Table.read('../naked_eye.fits')
Bright.sort('mag')

#Remove the Sun
Bright = Bright[1:]

#Reverse the Right Ascension (unit:hour) for outside-in perspective
Bright['ra_rv'] = 24 - Bright['ra']

#Plot the Sky Map, the ratio of length to width = 2
fig = plt.figure(figsize=(18, 9))

pieces = 12
piece_ang = 360/pieces
#12-piece is the most appropriate

for i in range(pieces):
	rect_scatter = [1./pieces*i, 0.0, 1./pieces, 1]
	ax = fig.add_axes(rect_scatter)
	max_ra, min_ra = piece_ang*(i+1), piece_ang*i
	Const = Bright[(Bright['ra_rv']*15 < max_ra) & (Bright['ra_rv']*15 >= min_ra)]
	if (min_ra == 0) or (max_ra == 360):
		rot = 12
		max_ra, min_ra = (piece_ang*(i+1) + 15*rot) % 360, (piece_ang*i + 15*rot) % 360
	else:
		rot = 0
	m = Basemap(projection = 'sinu', lon_0 = (piece_ang*i + 15*rot) % 360 + piece_ang/2., lat_0 = 0)
	Const['ra_rv'] = (Const['ra_rv'] + rot) % 24
	Bright['ra_rv_t'] = (Bright['ra_rv'] + rot) % 24
	conste_names = [j for j in set(Const['con'])]
  #Plot ecliptic
	m.drawgreatcircle(-(270 + 15*rot) % 360, -23.44, -(90.001 + 15*rot) % 360, +23.44, linewidth = 2, color = '#907910')
	m.drawgreatcircle(-(270.001 + 15*rot) % 360, -23.44, -(90 + 15*rot) % 360, +23.44, linewidth = 2, color = '#907910')
  m.drawparallels(np.arange(-90,90,30), color = '#c6c2b6')
	m.drawmeridians(np.arange((piece_ang*i + 15*rot) % 360,(piece_ang*(i+1) + 15*rot) % 360 + 0.01, piece_ang/2.), color = '#c6c2b6')
	ax.set_xlim(-2e7-0.166666e7, -2e7+0.166666e7)
	ax.set_ylim(0, 2e7)
	for s in range(len(conste_names)):
		lines = Table.read('../Lines_Art/' + conste_names[s] + '.csv', format='csv')
		Conste = Bright[Bright['con'] == conste_names[s]]
		Conste['ra_rv'] = (Conste['ra_rv'] + rot) % 24
		if len(lines) > 0:
			for k in range(int(len(lines)/2)):
				idx = np.argmax(Bright['hip'] == int(lines['hip'][2*k]))
				idy = np.argmax(Bright['hip'] == int(lines['hip'][2*k+1]))
				lons_0 = Bright['ra_rv_t'][idx]*15.
				lats_0 = Bright['dec'][idx]
				lons_1 = Bright['ra_rv_t'][idy]*15.
				lats_1 = Bright['dec'][idy]
				lons = [lons_0, lons_1]
				lats = [lats_0, lats_1]
				if (min_ra < lons_0 < max_ra) or (min_ra < lons_1 < max_ra):
					x, y = m(lons, lats)
					m.plot(x, y, marker = None, ms = 10, lw = 1, color = '#913228', alpha=1.)
				else:
					pass
		if conste_names[s] == 'Ser':
			continue
		else:
			pass
		bound = Table.read('../bound/' + conste_names[s] + '_Bound', format='ascii')
		bound.add_row([bound['col1'][0], bound['col2'][-1], conste_names[s]])
		bound.add_row([bound['col1'][0], bound['col2'][0], conste_names[s]])
		bound['ra'] = np.zeros(len(bound)) - 999.
		bound['ra_rv'] = np.zeros(len(bound)) - 999.
		bound['dec'] = np.zeros(len(bound)) - 999.
		for k in range(len(bound)):
			bound['ra'][k] = (int(bound['col1'][k][:2]) + int(bound['col1'][k][3:5])/60. + float(bound['col1'][k][6:])/3600.)*15
			bound['dec'][k]= float(bound['col2'][k])
			bound['ra'][k] = (bound['ra'][k] + 15*rot) % 360
			bound['ra_rv'][k] = 360 - bound['ra'][k]
			if (bound['ra_rv'][k] > max_ra):
				bound['ra_rv'][k] = max_ra
			elif (bound['ra_rv'][k] < min_ra):
				bound['ra_rv'][k] = min_ra
			else:
				pass
		for k in range(len(bound)-1):
			lons = [bound['ra_rv'][k],bound['ra_rv'][k+1]]
			lats = [bound['dec'][k],bound['dec'][k+1]]
			x, y = m(lons,lats)
			m.plot(x, y, marker=None, color='#AAAAAA', lw=1, linestyle=':', alpha=1.)
	for k in range(len(Const)):
		x, y = m(Const['ra_rv'][k]*15., Const['dec'][k])
    #Modify the marksize of the stars with different mag
		if Const['mag'][k] > 2.:
			m.plot(x, y, marker = 'o', color = 'black', ms = 20*(10**(-0.4*0.5*max(Const['mag'][k],0.5))), mec = 'white', mew = 3*(10**(-0.4*0.5*max(Const['mag'][k], 0.5))))
		else:
			m.plot(x, y, marker = '*', color = 'black', ms = 1.25*20*(10**(-0.4*0.5*max(Const['mag'][k],0.5))), mec = 'white', mew = 0.8*(10**(-0.4*0.5*max(Const['mag'][k],0.5))))
	y_l = np.pi * np.linspace(-0.5, 0, 100)
	y_u = np.pi * np.linspace(0.5, 0, 100)
	x_r = 2e7 * (1./12 * np.cos(y_u) + 1)
	x_l = 2e7 * (-1./12 * np.cos(y_u) + 1)

plt.savefig('Skymap.pdf')
