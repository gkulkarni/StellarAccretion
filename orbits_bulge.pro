PRO orbits

  datafile = '/data/shens/Eris_star_prop/stellar_particle_coordinates_bulge.dat'
  openr, 1, datafile

  nstars = 0
  nfiles = 0 
  readf, 1, nstars, nfiles 

  kpcunit = 0.0
  readf, 1, kpcunit

  x = fltarr(nstars, nfiles)
  y = fltarr(nstars, nfiles)
  z = fltarr(nstars, nfiles)
  readf, 1, x
  readf, 1, y
  readf, 1, z

  xcen = fltarr(nfiles)
  ycen = fltarr(nfiles)
  zcen = fltarr(nfiles)
  readf, 1, xcen
  readf, 1, ycen
  readf, 1, zcen

  redshift = fltarr(nfiles)
  readf, 1, redshift

  close, 1

  x_physical = fltarr(nstars, nfiles)
  y_physical = fltarr(nstars, nfiles)
  z_physical = fltarr(nstars, nfiles)

  help, x_physical
  print, nfiles
  help, redshift

  for i = 0, nstars-1 do begin
     for j = 0, nfiles-1 do begin 
        x_physical(i, j) = (x(i, j)-xcen(j))*kpcunit/(1.0+redshift(j))
        if (x(i, j) EQ 0.0) then x_physical(i, j) = 1.0e10 
        y_physical(i, j) = (y(i, j)-ycen(j))*kpcunit/(1.0+redshift(j))
        if (y(i, j) EQ 0.0) then y_physical(i, j) = 1.0e10 
        z_physical(i, j) = (z(i, j)-zcen(j))*kpcunit/(1.0+redshift(j))
        if (z(i, j) EQ 0.0) then z_physical(i, j) = 1.0e10 
     endfor 
  endfor

  openw, 1, 'xs_bulge.txt'
  for i = 0, nstars-1 do begin
     for j = 0, nfiles-1 do begin 
        printf, 1, x_physical(i, j) 
     endfor
  endfor
  close, 1 

  openw, 1, 'ys_bulge.txt'
  for i = 0, nstars-1 do begin
     for j = 0, nfiles-1 do begin 
        printf, 1, y_physical(i, j) 
     endfor
  endfor
  close, 1 

  openw, 1, 'zs_bulge.txt'
  for i = 0, nstars-1 do begin
     for j = 0, nfiles-1 do begin 
        printf, 1, z_physical(i, j) 
     endfor
  endfor
  close, 1 
  
END
