; radialdensityplot.pro

pro radialdensityplot, dens, rmax, $
	coordinate = co, momentum = mo, $
	col = col, key = key, legend = legend

; normalize to nuclear matter density
dens = dens/0.17

sz = size(dens)
npoints = sz(1)
r = findgen(npoints)*rmax/(npoints-1)

; set colors

tvlct, [0,255,0,0], [0,0,255,0], [0,0,0,255]
white = 255
black = 0
red = 1
green = 2 
blue = 3


; set names for x- and y-axes

if keyword_set(mo) then begin
	xtit = 'p [fm!u-1!n]'	
	ytit = '!mr!x(p) [fm!u3!n]'
	yr = [0.0005,4]
end else begin
	xtit = 'r [fm]'	
	ytit = '!mr!x(r) [!mr!x!d0!n]'
	yr = [0.00002,5]
end

xtick = ''
ytick = ''

; plot axes for setting plotting window (in white)

plot_io, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
	xtitle=xtit, ytitle=ytit, xtickname=xtick, $
	xstyle=1, ystyle=1, $
	ytickv = [0.0001, 0.001, 0.01, 0.1, 1], $
	ytickname = ['10!U-4!N', '10!U-3!N', '10!U-2!N', '10!U-1!N', '1'], $
	/nodata, color=white

; white the background

corner = convert_coord(!x.window, !y.window, /normal, /to_normal)
tv, [white,white], $
	corner(0,0), corner(1,0), $
	xsize = corner(0,1) - corner(0,0), $
	ysize = corner(1,1) - corner(1,0), $
	/normal

; plot the densities

if keyword_set(col) then begin
	plot_io, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, color=black, /noerase
	plot_io, r, dens(*,1), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, color=red, /noerase
	plot_io, r, dens(*,2), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, color=blue, /noerase
end else begin
	plot_io, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, linestyle=0, /noerase
	plot_io, r, dens(*,1), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, linestyle=2, /noerase
	plot_io, r, dens(*,2), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, linestyle=1, /noerase
end

plot_io, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
	xtitle=xtit, ytitle=ytit, xtickname=xtick, $
	xstyle=1, ystyle=1, $
	ytickv = [0.0001, 0.001, 0.01, 0.1, 1], $
	ytickname = ['10!U-4!N', '10!U-3!N', '10!U-2!N', '10!U-1!N', '1'], $
	/noerase, /nodata, color=black

pos = !p.position
xsize = pos(2)-pos(0)
ysize = pos(3)-pos(1)
if keyword_set(key) then begin
	xyouts, pos(0)+0.92*xsize, pos(1)+0.87*ysize, $	
		key, align=1.0, charsize=1.2, color=black, /norm 
end	

; legend

if keyword_set(legend) then begin
	if keyword_set(col) then begin
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.25*ysize, pos(1)+0.25*ysize], $	
			color=black, /norm
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.18*ysize, pos(1)+0.18*ysize], $	
			color=blue, /norm
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.11*ysize, pos(1)+0.11*ysize], $	
			color=red, /norm
	end else begin
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.25*ysize, pos(1)+0.25*ysize], $	
			linestyle=0, /norm
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.18*ysize, pos(1)+0.18*ysize], $	
			linestyle=1, /norm
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.11*ysize, pos(1)+0.11*ysize], $	
			linestyle=2, /norm
	end
	xyouts, pos(0)+0.17*(pos(2)-pos(0)), $
		pos(1)+0.24*(pos(3)-pos(1)), $	
		"total", color=black, /norm
	xyouts, pos(0)+0.17*(pos(2)-pos(0)), $
		pos(1)+0.17*(pos(3)-pos(1)), $	
		"neutron", color=black, /norm
	xyouts, pos(0)+0.17*(pos(2)-pos(0)), $
		pos(1)+0.1*(pos(3)-pos(1)), $	
		"proton", color=black, /norm 
end

end
