; initialize multiplot postscript environment


function initpost, file, nx, ny, $
	gapx = gapx, gapy = gapy, aspect = aspect, color = col

cmtodev = 1000.0

lborder = 2.00*cmtodev
rborder = 0.25*cmtodev
bborder = 2.00*cmtodev
tborder = 0.25*cmtodev

if not keyword_set(gapx) then gapx=0.4
gapx = gapx*cmtodev
if not keyword_set(gapy) then gapy=0.4
gapy = gapy*cmtodev
if not keyword_set(aspect) then aspect=1.0 

if nx*ny gt 2 then !p.charsize = 2.0

if nx eq 1 then begin
	totx=7.9*cmtodev 
	if keyword_set(col) then fontsz = 12 else fontsz = 10
end else if nx eq 2 then begin 
	totx=15.0*cmtodev
	if keyword_set(col) then fontsz = 12 else fontsz = 10
end else begin
	totx=17.5*cmtodev
	if keyword_set(col) then fontsz = 10 else fontsz = 8
end

plotx = (totx - lborder - rborder - (nx-1)*gapx)/nx
ploty = aspect*plotx

toty = ny*ploty + (ny-1)*gapy + bborder + tborder


set_plot, 'ps'
!p.font = 0

if keyword_set(col) then begin
!x.thick = 3.0
!y.thick = 3.0
!p.thick = 2.0
end

if keyword_set(col) then begin
	device, /encapsul, /color, bits=8, $
		font_size = fontsz, filename = file, $
		xsize = totx/1000.0, ysize = toty/1000.0
end else begin
 	device, /encapsul, font_size = fontsz, filename = file, $
		xsize = totx/1000.0, ysize = toty/1000.0
end

!p.multi = [0, nx, ny]

pos = fltarr(nx*ny, 4)

for i = 0, nx-1 do begin
	for j = 0, ny-1 do begin
		pos(j+i*ny, 0) = (lborder + i*(plotx+gapx))/totx
		pos(j+i*ny, 1) = (bborder + (ny-j-1)*(ploty+gapy))/toty
		pos(j+i*ny, 2) = pos(j+i*ny, 0) + plotx/totx
		pos(j+i*ny, 3) = pos(j+i*ny, 1) + ploty/toty		
	end
end

return, pos

end



pro closepost

device, /close
set_plot, 'x'

end


