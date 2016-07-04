pro nucchart, nucnn, nucnp, data, scale = scale

if not keyword_set(scale) then scale = [ -4.0, 4.0 ]

demin = scale[0]
demax = scale[1]

imin = 16
imax = 255

i0 = byte(imin+(imax-imin)/(demax-demin)*(-demin))

; max neutron and proton number

maxn = max(nucnn)
maxp = max(nucnp)

; magic numbers for annotations

magic = [2,8,20,28,50]

xtickv = magic[ where(magic le maxn) ]
ytickv = magic[ where(magic le maxp) ]

xtickv = [ xtickv, xtickv-1 ]
ytickv = [ ytickv, ytickv-1 ]

xticks = n_elements(xtickv)-1
yticks = n_elements(ytickv)-1

nnuc = n_elements(data)

; put scaled data in byte-array 

datab = make_array(maxn+2, maxp+2, /byte, value = 0)
for i = 0,nnuc-1 do $
	if (data[i] > demin and data[i] < demax) then $
		datab[nucnn[i], nucnp[i]] = byte(imin+(imax-imin)/(demax-demin)*(data[i]-demin)) $
	else $
		datab[nucnn[i], nucnp[i]] = 3

; empty white plot to get display coordinates 

plot, [-1, maxn+1], [-1,maxp+1], $
        xticks = xticks, yticks = yticks, $
        xtickv = xtickv, ytickv = ytickv, $
        xtitle = 'neutrons', ytitle = 'protons', $
	xminor = -1, yminor = -1, $	
        xstyle = 1, ystyle = 1, $
	xtickname = replicate(' ', xticks+1), $
	ytickname = replicate(' ', yticks+1), $
        /isotropic, /nodata, color=white,$
	position = [ 0.15, 0.20, 0.80, 0.98 ]

corner = convert_coord(!x.window, !y.window, /normal, /to_device)

; set color encoding

hue = fltarr(240)
sat = fltarr(240)
light = fltarr(240)

; hue[ indgen(i0-16) ] = 60.0
; sat[ indgen(i0-16) ] = 1.0 - 1.0*indgen(i0-16)/(i0-16)
; light[ indgen(i0-16) ] = 1.0

; hue[ i0-16+indgen(255-i0) ] = 0.0
; sat[ i0-16+indgen(255-i0) ] = 1.0*indgen(255-i0)/(255-i0)
; light[ i0-16+indgen(255-i0) ] = 1.0

; hue[ indgen(240) ] = 180.0*(1.0-indgen(240)/240.0)
; sat[ indgen(240) ] = 0.5
; light[ indgen(240) ] = 1.0

hue[ indgen(i0-16) ] = 120.0+60.0*(1.0-1.0*indgen(i0-16)/(i0-16))
sat[ indgen(i0-16) ] = 0.5
light[ indgen(i0-16) ] = 0.8

hue[ i0-16+indgen(255-i0) ] = 120.0*(1.0-1.0*indgen(255-i0)/(255-i0))
sat[ i0-16+indgen(255-i0) ] = 0.5
light[ i0-16+indgen(255-i0) ] = 0.8

tvlct, 255,255,255
tvlct, hue, sat, light, 16, /hls

; color coded data

tv, datab, corner(0,0), corner(1,0), $
	xsize = corner(0,1) - corner(0,0), $
     	ysize = corner(1,1) - corner(1,0), $
       	/device

plot, [-1, maxn+1], [-1,maxp+1], $
        xticks = xticks, yticks = yticks, $
        xtickv = xtickv, ytickv = ytickv, $
	xminor = -1, yminor = -1, $	
        xstyle = 1, ystyle = 1, $
        ticklen=1.0,$
        xtitle = 'neutrons', ytitle = 'protons', $
	xtickname = replicate(' ', xticks+1), $
	ytickname = replicate(' ', yticks+1), $
        /isotropic, /noerase, /nodata, color=black,$
	position = [ 0.15, 0.20, 0.80, 0.98 ]

magic = [2,8,20,28,50]

xtickv = magic[ where(magic le maxn) ]
ytickv = magic[ where(magic le maxp) ]

xticks = n_elements(xtickv)-1
yticks = n_elements(ytickv)-1

plot, [-0.5, maxn+1.5], [-0.5,maxp+1.5],$
	xticks = xticks, yticks = yticks, $
        xtickv = xtickv, ytickv = ytickv, $
	xstyle = 1, ystyle = 1, $
	ticklen = 0.0, $
	/isotropic, /noerase, /nodata, color=black,$
	position = [ 0.15, 0.20, 0.80, 0.98 ]

colorbar, bottom = 16, ncolors=240, /vertical, /right, $
	format = '(F5.1)', divisions=6, range = scale, $
	position = [ 0.85, 0.20, 0.90, 0.92 ]  

end



