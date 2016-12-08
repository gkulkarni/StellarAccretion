PRO rate

  datafile = '/data/shens/Eris_star_prop/accretion_rate_halo_agecut.dat'
  openr, 1, datafile

  nfiles = 400 
  
  ndisk = 0L
  readf, 1, ndisk
  print, ndisk

  x = fltarr(ndisk)
  y = fltarr(ndisk)
  z = fltarr(ndisk)
  readf, 1, x
  readf, 1, y
  readf, 1, z

  close, 1

  ;; kpcunit, xcen, ycen, and zen obtained using orbits.pro at 
  ;; z = 1.44329e-15.
  kpcunit = 90000.0
  xcen = 0.0861581
  ycen = -0.108580
  zcen = 0.0861788

  x_physical = (x-xcen)*kpcunit
  y_physical = (y-ycen)*kpcunit
  z_physical = (z-zcen)*kpcunit

  openw, 1, 'xs_disk.txt'
  for i = 0, ndisk-1 do begin
        printf, 1, x_physical(i)
  endfor
  close, 1 

  openw, 1, 'ys_disk.txt'
  for i = 0, ndisk-1 do begin
        printf, 1, y_physical(i)
  endfor
  close, 1 

  openw, 1, 'zs_disk.txt'
  for i = 0, ndisk-1 do begin
        printf, 1, z_physical(i)
  endfor
  close, 1 


end
