; densityplot.pro

pro densityplot, rho, view, min1, max1, min2, max2, $
	coordinate = co, momentum = mo, $
	cut = cut, integrated = int, $
	contour = cont, density = dens, $
	noxaxes = nox, noyaxes = noy, noframe = nof, $
        noannot = noan, key = key, logo = logo

; set colors

white = 255
black = 0

; x- and y-points

x = min1 + findgen(100)* (max1-min1)/99
y = min2 + findgen(100)* (max2-min2)/99

; set levels for contours

if keyword_set(cut) then begin
	if keyword_set(co) then begin
		maxrho = 5.0
		lvls = [0.001, 0.01, 0.1, 0.5, 1.0, 1.5]
		annt = ['0.001', '0.01', '0.1', '0.5', '1.0', '1.5']
;		maxrho = 20.0
;		lvls = [0.01, 0.1, 1.0, 4.0, 8.0]
;		annt = ['0.01', '0.1', '1.0', '4.0', '8.0']
	end else begin
		maxrho = 0.5
		lvls = [0.001, 0.01, 0.1, 0.2, 0.4]
		annt = ['0.001', '0.01', '0.1', '0.2',  '0.4']
	end
end else begin
	lvls = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
	annt = ['0.01', '0.1', '0.5', '1.0', '1.5', '2.0']
end	

; set names for x- and y-axes

if keyword_set(mo) then begin
	case view of
		0: begin
			xtit = 'k!dy!n [fm!u-1!n]'	
			ytit = 'k!dz!n [fm!u-1!n]'
		   end
		1: begin
			xtit = 'k!dx!n [fm!u-1!n]'	
			ytit = 'k!dz!n [fm!u-1!n]'
	           end
		2: begin
			xtit = 'k!dx!n [fm!u-1!n]'	
			ytit = 'k!dy!n [fm!u-1!n]'
	           end
	endcase
end else begin
	case view of
		0: begin
			xtit = 'y [fm]'	
			ytit = 'z [fm]'
		   end
		1: begin
			xtit = 'x [fm]'	
			ytit = 'z [fm]'
	           end
		2: begin
			xtit = 'x [fm]'	
			ytit = 'y [fm]'
	           end
	endcase
end

if keyword_set(nox) then begin
	xtit = ''
	xtick = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
end else begin
	xtick = ''
end
if keyword_set(noy) then begin
	ytit = ''
	ytick = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
end else begin
	ytick = ''
end


; plot axes for setting plotting window (in white)

plot, x, y, xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
	xtickname = xtick, ytickname = ytick, /nodata, color=white

corner = convert_coord(!x.window, !y.window, /normal, /to_device)

dims = size(rho)
rho = rebin(rho,3*dims(1),3*dims(2))

if keyword_set(dens) then begin
	if keyword_set(cut) then begin
		srho = bytscl(- rho, min=-maxrho, max=0.0)
	end else begin
		srho = bytscl(- rho, min=-maxrho, max=0.0)
	end
	tv, srho, $
		corner(0,0), corner(1,0), $
		xsize = corner(0,1) - corner(0,0), $
		ysize = corner(1,1) - corner(1,0), $
		/device
end

if keyword_set(cont) then begin
	if not keyword_set(noan) then $
		contour, rho, $	
			levels = lvls, c_annotation = annt, $
			xstyle = 5, ystyle = 5, $
			/follow, /noerase, color=black $
	else    contour, rho, $
			levels = lvls, c_labels = 0, $
			xstyle = 5, ystyle = 5, $
			/follow, /noerase, color=black
end

if not keyword_set(nof) then begin
    plot, x, y, xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
          xtickname = xtick, ytickname = ytick, /nodata, /noerase, color=black
end

if keyword_set(key) then begin
    xyouts, x(92), y(87), key, charsize=1.2, align=1.0, color=black
end	

device, /helvetica, /oblique
if keyword_set(logo) then begin
    xyouts, x(94), y(4), logo, charsize=0.9, align=1.0, color=black
end
device, /helvetica

end
