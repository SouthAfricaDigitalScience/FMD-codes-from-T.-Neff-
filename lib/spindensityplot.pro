; spindensityplot.pro


pro vecplot, u, v, x, y, maxrho = maxrho, color = col

; u(nx,ny), v(nx,ny), x(nx), y(ny) assumed
; plot every second density

s = size(u)
nx = s(1)
ny = s(2)

; minimum strength

minnorm = 0.05*maxrho

; factor to multiply vector length with

fac = 20.0/nx*(x(nx-1)-x(0))*0.25/maxrho

norm = sqrt(u*u + v*v)

for j=0,ny-1,2 do begin
	for i=0,nx-1,2 do begin
		if (norm(i,j)+norm(i+1,j+1))/2.0 gt minnorm then begin
			x0 = (x(i)+x(i+1))/2.0
			y0 = (y(j)+y(j+1))/2.0
			u0 = (u(i,j)+u(i+1,j+1))/4.0*fac
			v0 = (v(i,j)+v(i+1,j+1))/4.0*fac
			arrow, x0-u0, y0-v0, x0+u0, y0+v0, $
				/data, color=col, /solid, $
				hsize=-0.33
		end
	end
end

end


pro spindensityplot, rho, view, min1, max1, min2, max2, $
	coordinate = co, momentum = mo, $
	cut = cut, integrated = int, $
	vector = vec, density = dens, $
	noxaxes = nox, noyaxes = noy, key = key

; set colors

white = 255
black = 0

dims = size(rho)
dimx = dims(1)
dimy = dims(2)

; x- and y-points

x = min1 + findgen(dimx)* (max1-min1)/(dimx-1)
y = min2 + findgen(dimy)* (max2-min2)/(dimy-1)

; which directions to project out

case view of
	0: begin
		nx = 1
		ny = 2
	   end
	1: begin
		nx = 0
		ny = 2
	   end
	2: begin
		nx = 0
		ny = 1
	   end
endcase

; set names for x- and y-axes

if keyword_set(mo) then begin
	maxrho = 0.05
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
	maxrho = 0.1 
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

arho = sqrt(rho(*,*,0)*rho(*,*,0) + $ 
            rho(*,*,1)*rho(*,*,1) + $ 
	    rho(*,*,2)*rho(*,*,2)) 

arho = rebin(arho, 2*dimx, 2*dimy)

if keyword_set(dens) then begin

	if keyword_set(cut) then begin
		srho = bytscl(- arho, min=-maxrho, max=0.0)
	end else begin
		srho = bytscl(- arho, min=-maxrho, max=0.0)
	end
	tv, srho, $
		corner(0,0), corner(1,0), $
		xsize = corner(0,1) - corner(0,0), $
		ysize = corner(1,1) - corner(1,0), $
		/device
end

if keyword_set(vec) then begin
	vecplot, rho(*,*,nx), rho(*,*,ny), x, y, maxrho = maxrho, color=black
end

plot, x, y, xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
	xtickname = xtick, ytickname = ytick, /nodata, /noerase, color=black

if keyword_set(key) then begin
	xyouts, x(fix(0.92*dimx)), y(fix(0.87*dimy)), $
		key, align=1.0, charsize=1.2, color=black
end	

end
