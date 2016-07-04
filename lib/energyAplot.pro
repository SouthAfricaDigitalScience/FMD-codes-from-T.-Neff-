; energyyplot.pro
; plot Binding energy over nucleon number

pro energyAplot, Afmd, energiesfmd, Aexp, energiesexp, $
	Amin = Amin, Amax = Amax, $
	Emin = Emin, Emax = Emax, $ 
	col = col, key = key, legend = legend

; set colors

tvlct, [0,255,0,0], [0,0,255,0], [0,0,0,255]
white = 255
black = 0
red = 1
green = 2 
blue = 3

xr = [Amin, Amax]
yr = [Emin, Emax]

; set names for x- and y-axes

xtit = 'A'	
ytit = 'E!dB!n/A [MeV]'

xtick = ''
ytick = ''

; plot axes for setting plotting window (in white)

plot, Afmd, energiesfmd, xrange=xr, yrange=yr, $
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

; plot the energies

; if keyword_set(col) then begin
	plot, Afmd, energiesfmd, xrange=xr, yrange=yr, $
		xstyle=5, ystyle=5, color=blue, /noerase
	plot, Aexp, energiesexp, xrange=xr, yrange=yr, $
		xstyle=5, ystyle=5, color=red, /noerase
; end

plot, Afmd, energiesfmd, xrange=xr, yrange=yr, $
	xtitle=xtit, ytitle=ytit, xtickname=xtick, $
	xstyle=1, ystyle=1, $
	/noerase, /nodata, color=black

pos = !p.position
xsize = pos(2)-pos(0)
ysize = pos(3)-pos(1)
if keyword_set(key) then begin
	xyouts, pos(0)+0.92*xsize, pos(1)+0.87*ysize, $	
		key, align=1.0, charsize=1.0, color=black, /norm 
end	

; legend

if keyword_set(legend) then begin
;	if keyword_set(col) then begin
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.2*ysize, pos(1)+0.2*ysize], $	
			color=blue, /norm
		plots, [pos(0)+0.05*xsize, pos(0)+0.15*xsize] , $
		          	[pos(1)+0.13*ysize, pos(1)+0.13*ysize], $	
			color=red, /norm
;	end
	xyouts, pos(0)+0.17*(pos(2)-pos(0)), $
		pos(1)+0.18*(pos(3)-pos(1)), $	
		"FMD", charsize=1.0, color=black, /norm
	xyouts, pos(0)+0.17*(pos(2)-pos(0)), $
		pos(1)+0.11*(pos(3)-pos(1)), $	
		"Experiment", charsize=1.0, color=black, /norm
end

end
