; chargedensityplot.pro

pro chargedensityplot, dens, rmax, $
	col = col, key = key

; normalize to nuclear matter density
; dens = dens/0.17

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

xtit = 'r [fm]'	
ytit = '!mr!x(r) [fm!u-3!n]'
rmax=8.0
yr = [0,0.12]

; plot axes for setting plotting window (in white)

plot, r, dens(*), xrange=[0,rmax], yrange=yr, $
	xtitle=xtit, ytitle=ytit, $
	xstyle=1, ystyle=1, $
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
	plot, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, color=black, /noerase
end else begin
	plot, r, dens(*,0), xrange=[0,rmax], yrange=yr, $
		xstyle=5, ystyle=5, linestyle=0, /noerase
end

plot, r, dens(*), xrange=[0,rmax], yrange=yr, $
	xtitle=xtit, ytitle=ytit, $
	xstyle=1, ystyle=1, $
	/noerase, /nodata, color=black

pos = !p.position
xsize = pos(2)-pos(0)
ysize = pos(3)-pos(1)
if keyword_set(key) then begin
	xyouts, pos(0)+0.92*xsize, pos(1)+0.87*ysize, $	
		key, align=1.0, charsize=1.2, color=black, /norm 
end	

end
