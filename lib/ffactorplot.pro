; ffactorplot.pro

pro ffactorplot, dens, qmax, qplotmax, n, $
	col = col, key = key, legend = legend

; calculate total density
if n gt 1 then totdens = total(dens, 2) else totdens = dens

sz = size(dens)
npoints = sz(1)
r = findgen(npoints)*qmax/(npoints-1)

; set colors

tvlct, [0,255,0,0], [0,0,255,0], [0,0,0,255]
white = 255
black = 0
red = 1
green = 2 
blue = 3


; set names for x- and y-axes

xtit = 'q [fm!u-1!n]'
ytit = '|F(q)|!u2!n'
xr = [0, qplotmax]	
yr = [0.8e-5, 1.4]

xtick = ''
ytick = ''

ytickss = 5
ytickvs = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
yticknames = ['10!U-5!N', '10!U-4!N', '10!U-3!N', '10!U-2!N', '10!U-1!N', '1']

; plot axes for setting plotting window (in white)

plot_io, r, dens(*,0), xrange=xr, yrange=yr, $
	xtitle=xtit, ytitle=ytit, xtickname=xtick, $
	xstyle=1, ystyle=1, $
        yticks = ytickss, ytickv = ytickvs, ytickname = yticknames, $
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
        if n gt 1 then begin
            for i=0,n-1 do begin
                plot_io, r, dens(*,i), xrange=xr, yrange=yr, $
		        thick=0.5, xstyle=5, ystyle=5, color=gray, /noerase
            endfor
        endif
	plot_io, r, totdens, xrange=xr, yrange=yr, $
		xstyle=5, ystyle=5, color=red, /noerase
end else begin
	plot_io, r, totdens, xrange=xr, yrange=yr, $
		xstyle=5, ystyle=5, linestyle=0, /noerase
end

plot_io, r, dens(*,0), xrange=xr, yrange=yr, $
	xtitle=xtit, ytitle=ytit, xtickname=xtick, $
	xstyle=1, ystyle=1, $
        yticks = ytickss, ytickv = ytickvs, ytickname = yticknames, $
	/noerase, /nodata, color=black

pos = !p.position
xsize = pos(2)-pos(0)
ysize = pos(3)-pos(1)
if keyword_set(key) then begin
	xyouts, pos(0)+0.92*xsize, pos(1)+0.87*ysize, $	
		key, align=1.0, charsize=1.2, color=black, /norm 
end	

; legend


end
