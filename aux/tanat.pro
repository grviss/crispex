;-------------------------------------------------------------
;+
; NAME:
;      	TANAT: Time slice ANAlysis Tool
;
; PURPOSE:
;	Analyse the output timeslices from CRISPEX in order to
;	determine velocities.
;
; CATEGORY:
;      	Data analysis
;
; CALLING SEQUENCE:
;	TANAT, filename, LINE_CENTER=line_center, ASECPIX=asecpix, DT=dt
;
; INPUTS:
;	filename	= filename of the timeslice
;
; KEYWORDS:
;	LINE_CENTER	= Not set: linecentre is determined from the data or from SPECTFILE.
;			  Integer scalar: linecentre is set to position specified by LINE_CENTER.
;			  2-element array (format: [WAVELENGTH, DELTA_LAMBDA]): linecentre is
;				determined from the data or from SPECTFILE and set to WAVELENGTH. 
;				The distance in wavelength between the linepositions is specified 
;				by DELTA LAMBDA.
;			  3-element array (format: [Integer scalar, WAVELENGTH, DELTA_LAMBDA]):
;				combination of the two above, as linecentre is specified by the
;				integer scalar and set to WAVELENGTH with tickmarks at DELTA_LAMBDA
;				distances from eachother.
;	ASECPIX		= Specifies the arcseconds per pixel. Defaults to current SST default of 0.0592.
;	DT		= Specifies the elapsed time in seconds per time step. Defaults to 1.
;
; OUTPUTS:
;	Measurements can be saved to a file, specified by the text in the filename input text box in
;	the upper right part of the Set up tab. Each measurement is written on a single line,
;	containing: the filename (i.e. of the saved timeslice, a *csav file), the spectral position at
;	which the measurement is taken, a measurement ID, the first spatial coordinate, the second
;	spatial coordinate, the actual distance in pixels between them, the first temporal coordinate,
;	the second temporal coordinate, the difference in temporal coordinates, the absolute speed
;	(in km/s) and an optional flag.
;	For fitted measurements the order is slightly altered: after the measurement ID follow first
;	the three spatial coordinates, the three temporal coordinates, the accelaration (in m/s^2) and 
;	finally an optional flag.
;
; COMMON BLOCKS:
;     	None.
;
; SIDE EFFECTS:
;     	None.
;
; RESTRICTIONS:
;     	Input must be in CRISPEX timeslice output format, i.e. certain
;	variable naming in the savefile is assumed.
;
; PROCEDURE:
;	Simple call (i.e. without any keywords) will yield a window containing
;	controls, the timeslice itself and a plot of the detailed spectrum. Spectral position as well
;	as conversion (i.e. arcseconds per pixel and seconds per time step) can be modified. Browsing
;	with the cursor over the timeslice will show the current frame number and pixel along the slice,
;	as well as the detailed spectrum at the current pixel (if such information is available of course). 
;
;	Left-clicking with the mouse will store and lock the current cursor position. Subsequently a
;	straight line will be drawn from the locked position to the then current cursor position enabling
;	the determination of the speed of a (constant velocity) feature. Re-clicking with the left mouse
;	button at a different location will relocate the locked position. Clicking with the middle mouse
;	button will fix the other end of the line. Right clicking will release the cursor lock.
;
; MODIFICATION HISTORY:
;	27 Mar 2009 GV:	setup of original version
;	30 Mar 2009 GV: release of beta version						(v0.9)
;	06 Apr 2009 GV: corrected errors in read-in and implemented calculation of	(v0.9.1)
;			correct distances
;	14 Apr 2009 GV: adapted to accomodate change in CRISPEX output save file	(v0.9.2)
;	22 May 2009 GV: implemented save measurement to file option			(v0.9.3)
;	24 Aug 2009 GV: implemented visualisation options (scaling and combination	(v0.9.4)
;			of spectral positions) and distributed options to tabs
;	22 Aug 2011 GV: extensive makeover, implemented parabolic fit, overlay of 	(v1.0)
;			measurements, option to open other files in-program, fixed
;			several display and saving bugs
;
; AUTHOR:
;	Gregal Vissers (g.j.m.vissers@astro.uio.no)
;	@ Institute for Theoretical Astrophysics, University of Oslo
;-
;-------------------------------------------------------------

;================================================================================= ABOUT WINDOW PROCEDURES
PRO TANAT_ABOUT_WINDOW, event 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_ABOUT_WINDOW'
	abouttlb = WIDGET_BASE(TITLE = 'TANAT: ABOUT', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS, TLB_SIZE_EVENTS = 0)
	disp = WIDGET_BASE(abouttlb, /COLUMN)
	aboutdrawid = WIDGET_DRAW(disp, XSIZE = (*(*info).winsizes).aboutwinx, YSIZE = (*(*info).winsizes).aboutwiny, RETAIN = 2)
	WIDGET_CONTROL, abouttlb, /REALIZE, TLB_SET_XOFFSET=(*(*info).winsizes).aboutwinxoffset, TLB_SET_YOFFSET=(*(*info).winsizes).aboutwinyoffset
	WIDGET_CONTROL, aboutdrawid, GET_VALUE = aboutwid
	TANAT_UPDATE_STARTUP_FEEDBACK, (*(*info).feedbparams).startup_im, (*(*info).feedbparams).xout, (*(*info).feedbparams).yout, $
		['Running TANAT version '+(*(*info).versioninfo).version_number+' ('+(*(*info).versioninfo).revision_number+')','',$
		'Developed by: Gregal Vissers', $
		'               Institute of Theoretical Astrophysics,',$
		'               University of Oslo',$
		'               2009-2011']
	WIDGET_CONTROL, aboutdrawid, EVENT_PRO = 'TANAT_ABOUT_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS, /DRAW_BUTTON_EVENTS
	WIDGET_CONTROL, abouttlb, SET_UVALUE = info
	XMANAGER, 'TANAT', abouttlb,/NO_BLOCK
	(*(*info).winids).abouttlb = abouttlb
END

PRO TANAT_ABOUT_CURSOR, event
; Handles cursor actions on the about window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_ABOUT_CURSOR'
	IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_DRAW' THEN BEGIN
		CASE event.TYPE OF
		0:	CASE event.PRESS OF
			1:	BEGIN	; left mouse button press
					WIDGET_CONTROL, (*(*info).winids).abouttlb, /DESTROY
					(*(*info).winids).abouttlb = 0
				END
			ELSE: BREAK
			ENDCASE
		ELSE: RETURN
		ENDCASE
	ENDIF
END

PRO TANAT_ABOUT_WINDOW_CLOSE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_ABOUT_WINDOW_CLOSE'
	WIDGET_CONTROL, (*(*info).winids).abouttlb, /DESTROY
END

;================================================================================= CALCULATE PROCEDURES
PRO TANAT_CALCULATE_DELTA, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_CALCULATE_DELTA'
	(*(*info).measparams).delta_lx = FLOAT((*(*info).dataparams).lx) - FLOAT((*(*info).curs).lxlock)
	(*(*info).measparams).delta_t = FLOAT((*(*info).dataparams).t) - FLOAT((*(*info).curs).tlock)
	IF (ABS((*(*info).measparams).delta_lx) GT 0.) THEN BEGIN
		IF ((*(*info).measparams).delta_lx GT 0.) THEN BEGIN
			lower = (*(*info).curs).lxlock
			upper = lower + (*(*info).measparams).delta_lx - 1
		ENDIF ELSE BEGIN
			upper = (*(*info).curs).lxlock
			lower = upper + (*(*info).measparams).delta_lx + 1
		ENDELSE
		(*(*info).measparams).act_delta_lx = FLOAT(TOTAL((*(*(*info).dataparams).lxdist)[lower:upper])) 
	ENDIF ELSE (*(*info).measparams).act_delta_lx = 0.
	time = FLOAT((*(*info).measparams).delta_t) * FLOAT((*(*info).measparams).secondststep)
	distance = 149597870. * TAN(FLOAT((*(*info).measparams).arcsecpix) * FLOAT((*(*info).measparams).act_delta_lx) / 3600. * !DTOR)
	IF ((time EQ 0.) OR ((time EQ 0.) AND (distance EQ 0.))) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_text, SET_VALUE = 'UNDEF'
	ENDIF ELSE BEGIN
		(*(*info).measparams).speed = FLOAT(distance) / ABS(FLOAT(time))
		WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_text, SET_VALUE = STRTRIM((*(*info).measparams).speed, 2)
	ENDELSE

	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_lx_text, SET_VALUE = STRTRIM((*(*info).measparams).act_delta_lx,2)
	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_t_text, SET_VALUE = STRTRIM(ABS((*(*info).measparams).delta_t),2)
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_params_text, SET_VALUE = '('+STRTRIM(FIX((*(*info).curs).lxlock),2)+','+STRTRIM(FIX((*(*info).dataparams).lx),2)+')'
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_params_text, SET_VALUE = '('+STRTRIM(FIX((*(*info).curs).tlock),2)+','+STRTRIM(FIX((*(*info).dataparams).t),2)+')'
END

;================================================================================= CLOSING AND PROGRAM EXIT PROCEDURES
PRO TANAT_CLEANUP, base
	WIDGET_CONTROL, base, GET_UVALUE = info
	PTR_FREE, info
END

PRO TANAT_CLOSE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info, /DESTROY
	PTR_FREE, info
END

PRO TANAT_CLOSE_EVENT_WINDOW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_CLOSE_EVENT_WINDOW'
	IF (event.TOP EQ (*(*info).winids).errtlb) THEN (*(*info).winids).errtlb = 0
	IF (event.TOP EQ (*(*info).winids).savnamenulltlb) THEN (*(*info).winids).savnamenulltlb = 0
	IF (event.TOP EQ (*(*info).winids).savnamedoubletlb) THEN (*(*info).winids).savnamedoubletlb = 0
	WIDGET_CONTROL, event.TOP, /DESTROY
END

;================================================================================= CURSOR PROCEDURES
PRO TANAT_CURSOR, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_CURSOR'
	IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_TRACKING' THEN BEGIN
		IF event.ENTER THEN BEGIN
			WIDGET_CONTROL, event.HANDLER, get_value = wid
			WSET, wid
			ci = UINTARR(16) & cim = ci & cim[8] = 1
			DEVICE, CURSOR_IMAGE = ci, CURSOR_MASK = cim, CURSOR_XY = [8,8]
		ENDIF ELSE BEGIN
			IF ((*(*info).measparams).parabolic_fit AND ((*(*info).measparams).np GE 1)) THEN BEGIN
				*(*(*info).measparams).lx_array = (*(*(*info).measparams).lx_array)[0:(*(*info).measparams).np-1]
				*(*(*info).measparams).t_array = (*(*(*info).measparams).t_array)[0:(*(*info).measparams).np-1]
				*(*(*info).measparams).slx_array = (*(*(*info).measparams).slx_array)[0:(*(*info).measparams).np-1]
				*(*(*info).measparams).st_array = (*(*(*info).measparams).st_array)[0:(*(*info).measparams).np-1]
				TANAT_PARABOLIC_FIT, event
				TANAT_DRAW, event
			ENDIF 
			DEVICE, /CURSOR_CROSSHAIR
		ENDELSE
	ENDIF ELSE IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_DRAW' THEN BEGIN
		CASE event.TYPE OF
		0:	CASE event.PRESS OF
			1:	BEGIN	; left mouse button
					(*(*info).curs).lockset = 1
					IF ((*(*info).measparams).parabolic_fit EQ 0) THEN BEGIN
						WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_lx_label, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_t_label, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_lx_text, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_t_text, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_label, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_text, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_params_label, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_params_text, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).t_params_label, SENSITIVE = 1
						WIDGET_CONTROL, (*(*info).ctrlsmeas).t_params_text, SENSITIVE = 1
					ENDIF
					(*(*info).curs).slxlock = event.X
					(*(*info).curs).stlock = event.Y
					(*(*info).curs).lxlock = FLOAT((*(*info).curs).slxlock) * (*(*info).dataparams).nlx / (*(*info).winsizes).windowx
					IF ((*(*info).dataparams).d_nt NE (*(*info).dataparams).nt) THEN $
						(*(*info).curs).tlock = FLOAT((*(*info).curs).stlock) * (*(*info).dataparams).d_nt / (*(*info).winsizes).windowy + (*(*info).dataparams).t_low ELSE $
						(*(*info).curs).tlock = FLOAT((*(*info).curs).stlock) * (*(*info).dataparams).nt / (*(*info).winsizes).windowy
					(*(*info).dataparams).lx = (*(*info).curs).lxlock
					(*(*info).dataparams).t = (*(*info).curs).tlock
					TANAT_CURSOR_SET_COORDSLIDERS, 0, 0, event
					IF (*(*info).measparams).parabolic_fit THEN BEGIN
						(*(*info).measparams).np += 1
						IF ((*(*info).measparams).np LT 2) THEN BEGIN
							*(*(*info).measparams).lx_array = (*(*info).curs).lxlock
							*(*(*info).measparams).t_array = (*(*info).curs).tlock
							*(*(*info).measparams).slx_array = (*(*info).curs).slxlock
							*(*(*info).measparams).st_array = (*(*info).curs).stlock
						ENDIF ELSE BEGIN
							*(*(*info).measparams).lx_array = [*(*(*info).measparams).lx_array, (*(*info).curs).lxlock]
							*(*(*info).measparams).t_array = [*(*(*info).measparams).t_array, (*(*info).curs).tlock]
							*(*(*info).measparams).slx_array = [*(*(*info).measparams).slx_array, (*(*info).curs).slxlock]
							*(*(*info).measparams).st_array = [*(*(*info).measparams).st_array, (*(*info).curs).stlock]
						ENDELSE
						IF ((*(*info).measparams).np GE 2) THEN BEGIN
							WIDGET_CONTROL, (*(*info).ctrlsmeas).rem_but, /SENSITIVE
							WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_label, /SENSITIVE
						ENDIF
					ENDIF
				END
			2:	BEGIN	; middle mouse button
					IF (*(*info).curs).lockset THEN BEGIN
						(*(*info).curs).lockset = 2
						(*(*info).curs).slx = event.X
						(*(*info).curs).st = event.Y
						(*(*info).dataparams).lx = FLOAT((*(*info).curs).slx) * (*(*info).dataparams).nlx / (*(*info).winsizes).windowx
						IF ((*(*info).dataparams).d_nt NE (*(*info).dataparams).nt) THEN $
							(*(*info).dataparams).t = FLOAT((*(*info).curs).st) * (*(*info).dataparams).d_nt / (*(*info).winsizes).windowy + (*(*info).dataparams).t_low ELSE $
							(*(*info).dataparams).t = FLOAT((*(*info).curs).st) * (*(*info).dataparams).nt / (*(*info).winsizes).windowy
						IF (*(*info).measparams).parabolic_fit THEN BEGIN
							(*(*info).measparams).np += 1
							*(*(*info).measparams).lx_array = [*(*(*info).measparams).lx_array, (*(*info).dataparams).lx]
							*(*(*info).measparams).t_array = [*(*(*info).measparams).t_array, (*(*info).dataparams).t]
							*(*(*info).measparams).slx_array = [*(*(*info).measparams).slx_array, (*(*info).curs).slx]
							*(*(*info).measparams).st_array = [*(*(*info).measparams).st_array, (*(*info).curs).st]
							IF ((*(*info).measparams).np EQ 3) THEN WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, /SENSITIVE ELSE WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, SENSITIVE = 0
						ENDIF ELSE WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, /SENSITIVE
						WIDGET_CONTROL, (*(*info).ctrlsmeas).flag_label, /SENSITIVE
						WIDGET_CONTROL, (*(*info).ctrlsmeas).flag_text, /SENSITIVE
						TANAT_CURSOR_SET_COORDSLIDERS, 0, 0, event
						IF (*(*info).measparams).parabolic_fit THEN TANAT_PARABOLIC_FIT, event
					ENDIF ELSE RETURN
				END
			4:	BEGIN	; right mouse button
					(*(*info).curs).lockset = 0
					TANAT_CURSOR_SET_COORDSLIDERS, 1, 1, event
					TANAT_RESET_OUTPUTS, event
				END
			ELSE: BREAK
			ENDCASE
		2:	BEGIN
				IF ((*(*info).curs).lockset NE 2) THEN BEGIN
					(*(*info).curs).slx = event.X
					(*(*info).curs).st = event.Y
					(*(*info).dataparams).lx = FLOAT((*(*info).curs).slx) * (*(*info).dataparams).nlx / (*(*info).winsizes).windowx
					IF ((*(*info).dataparams).d_nt NE (*(*info).dataparams).nt) THEN $
						(*(*info).dataparams).t = FLOAT((*(*info).curs).st) * (*(*info).dataparams).d_nt / (*(*info).winsizes).windowy+(*(*info).dataparams).t_low ELSE $
						(*(*info).dataparams).t = FLOAT((*(*info).curs).st) * (*(*info).dataparams).nt / (*(*info).winsizes).windowy
					IF ((*(*info).curs).lockset EQ 1) THEN BEGIN
						TANAT_CURSOR_SET_COORDSLIDERS, 0, 0, event
						IF ((*(*info).measparams).parabolic_fit EQ 0) THEN TANAT_CALCULATE_DELTA, event
						IF ((*(*info).measparams).parabolic_fit AND ((*(*info).measparams).np GE 1)) THEN BEGIN
							*(*(*info).measparams).lx_array = [(*(*(*info).measparams).lx_array)[0:(*(*info).measparams).np-1], (*(*info).dataparams).lx]
							*(*(*info).measparams).t_array = [(*(*(*info).measparams).t_array)[0:(*(*info).measparams).np-1], (*(*info).dataparams).t]
							*(*(*info).measparams).slx_array = [(*(*(*info).measparams).slx_array)[0:(*(*info).measparams).np-1], (*(*info).curs).slx]
							*(*(*info).measparams).st_array = [(*(*(*info).measparams).st_array)[0:(*(*info).measparams).np-1], (*(*info).curs).st]
							TANAT_PARABOLIC_FIT, event
						ENDIF
					ENDIF ELSE BEGIN
						TANAT_CURSOR_SET_COORDSLIDERS, 1, 1, event
					ENDELSE
				ENDIF ELSE RETURN
			END
		ELSE: RETURN
		ENDCASE
	ENDIF
	TANAT_DRAW, event
END

PRO TANAT_CURSOR_SET_COORDSLIDERS, xsensitive, ysensitive, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_CURSOR_SET_COORDSLIDERS'
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_slider, SET_VALUE = (*(*info).dataparams).lx, SENSITIVE = xsensitive
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_slider, SET_VALUE = (*(*info).dataparams).t, SENSITIVE = ysensitive
END
 
;================================================================================= BINARY CONVERSION FUNCTIONS
FUNCTION TANAT_DEC2BIN, decimal_number
; Handles the conversion of decimal to binary number
	IF (decimal_number NE 0) THEN BEGIN
		coarse_p = ALOG10(DOUBLE(decimal_number))/ALOG10(2.D)
		p = FLOOR(coarse_p)
		firstpass = 1
		i=0
		binary_array = INTARR(FLOOR(coarse_p)+1)
		b = REVERSE((2^FINDGEN(FLOOR(coarse_p)+1)))
		WHILE (p GE 0) DO BEGIN
			IF firstpass THEN BEGIN
				IF (2^p LE decimal_number) THEN BEGIN
					binary_array[i] = 1
					decimal_number = decimal_number - 2^p
				ENDIF ELSE binary_array[i] = 0
			ENDIF ELSE BEGIN
				IF (2^p LE decimal_number) THEN BEGIN
					binary_array[i] = 1
					decimal_number = decimal_number - 2^p
				ENDIF ELSE binary_array[i] = 0
				firstpass = 0
			ENDELSE
			p = p - 1
			i = i + 1
		ENDWHILE
		binary_array = REVERSE(binary_array)
		IF (N_ELEMENTS(binary_array) LT 5) THEN binary_array = [binary_array, REPLICATE(0,5-N_ELEMENTS(binary_array))]
	ENDIF ELSE binary_array = [0,0,0,0,0]
	RETURN, binary_array
END

;================================================================================= DRAW PROCEDURES
PRO TANAT_DRAW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_DRAW'
	WSET, (*(*info).winids).wid
	IF (*(*info).dispswitch).man_scale THEN BEGIN
		minimum = (*(*info).scaling).range / 100. * (*(*info).scaling).min_val + (*(*info).scaling).minimum
		maximum = (*(*info).scaling).range / 100. * (*(*info).scaling).max_val + (*(*info).scaling).minimum
	ENDIF ELSE BEGIN
		minimum = (*(*info).scaling).minimum
		maximum = (*(*info).scaling).range + (*(*info).scaling).minimum
	ENDELSE
	IF (*(*info).dispswitch).smoothed THEN TV,BYTSCL(CONGRID(*(*(*info).data).loopslice,(*(*info).winsizes).windowx,(*(*info).winsizes).windowy,/INTERP), MIN = minimum, MAX = maximum) $
	ELSE TV,BYTSCL(CONGRID(*(*(*info).data).loopslice,(*(*info).winsizes).windowx,(*(*info).winsizes).windowy), MIN = minimum, MAX = maximum)
	IF (*(*info).dispswitch).overlay_measurements THEN BEGIN
		FOR k=0,999 DO BEGIN
			IF ((*(*(*info).overlays).lxbpoint)[k] NE (*(*(*info).overlays).lxepoint)[k]) AND ((*(*(*info).overlays).tbpoint)[k] NE (*(*(*info).overlays).tepoint)[k]) THEN BEGIN
				PLOTS,[(*(*(*info).overlays).lxbpoint)[k], (*(*(*info).overlays).lxepoint)[k]],[((*(*(*info).overlays).tbpoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt,$
					((*(*(*info).overlays).tepoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt],/DEVICE
				PLOTS,[(*(*(*info).overlays).lxbpoint)[k], (*(*(*info).overlays).lxepoint)[k]],[((*(*(*info).overlays).tbpoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt,$
					((*(*(*info).overlays).tepoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt],/DEVICE, PSYM = 6
				PLOTS,[(*(*(*info).overlays).lxbpoint)[k], (*(*(*info).overlays).lxepoint)[k]],[((*(*(*info).overlays).tbpoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt,$
					((*(*(*info).overlays).tepoint)[k]-(*(*info).dataparams).t_low)*(*(*info).winsizes).windowy/(*(*info).dataparams).d_nt],/DEVICE, PSYM = 7
				IF ((*(*(*info).overlays).lxepoint)[k] GT (*(*(*info).overlays).lxbpoint)[k]) THEN BEGIN
					upp_x = (*(*(*info).overlays).lxepoint)[k]
					low_x = (*(*(*info).overlays).lxbpoint)[k]
				ENDIF ELSE BEGIN
					upp_x = (*(*(*info).overlays).lxbpoint)[k]
					low_x = (*(*(*info).overlays).lxepoint)[k]
				END
				IF ((*(*(*info).overlays).tepoint)[k] GT (*(*(*info).overlays).tbpoint)[k]) THEN BEGIN
					upp_t = ((*(*(*info).overlays).tepoint)[k]-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
					low_t = ((*(*(*info).overlays).tbpoint)[k]-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
				ENDIF ELSE BEGIN
					upp_t = ((*(*(*info).overlays).tbpoint)[k]-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
					low_t = ((*(*(*info).overlays).tepoint)[k]-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
				END
				XYOUTS, (upp_x-low_x)/2.+low_x, (upp_t-low_t)/2.+low_t+5, (*(*(*info).measparams).meas_id)[k], /DEVICE
			ENDIF
		ENDFOR
	ENDIF
	TANAT_DRAW_CURSCROSS_PLOT, event
	IF (*(*info).dispswitch).slab_set THEN BEGIN
		WSET, (*(*info).winids).ls_drawid
		PLOT, *(*(*info).dataparams).lps, *(*(*info).dataparams).spec, POS=[(*(*info).plotpos).lsx0,(*(*info).plotpos).lsy0,(*(*info).plotpos).lsx1,(*(*info).plotpos).lsy1],/NORM, CHARSIZE=1,YS=1,$
			YR=[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y],XSTYLE = (*(*info).dispswitch).v_dop_set * 8 + 1, XTITLE = (*(*info).plotaxes).spxtitle, YTITLE = 'Scaled intensity', $
			BACKGROUND = (*(*info).plotparams).bgplotcol,COLOR = (*(*info).plotparams).plotcol, /NODATA
		IF ( (*(*info).dispswitch).v_dop_set EQ 1) THEN BEGIN
			AXIS, XAXIS=1, XRANGE = [((*(*info).plotaxes).v_dop)[0], ((*(*info).plotaxes).v_dop)[(*(*info).dataparams).lp_last]], XSTYLE=1, XTITLE = 'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol
			XYOUTS,(*(*(*info).dataparams).lps)[FLOOR((*(*info).dataparams).lp_last-6/23.*(*(*info).dataparams).lp_last)],0.9*(*(*info).plotaxes).ls_yrange + (*(*info).plotaxes).ls_low_y, $
				STRTRIM(STRING((*(*info).plotaxes).v_dop[(*(*info).dataparams).lp],FORMAT='(3(F9.2,x))'),2)+' km/s', COLOR = (*(*info).plotparams).plotcol
		ENDIF
		IF ((*(*info).dispswitch).singlepos EQ 0) THEN BEGIN
			PLOTS, [1,1] * (*(*(*info).dataparams).lps)[(*(*info).dataparams).low_low_val],[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y], COLOR = (*(*info).plotparams).plotcol, LINESTYLE = 5
			PLOTS, [1,1] * (*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_upp_val],[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y], COLOR = (*(*info).plotparams).plotcol, LINESTYLE = 5
			IF (*(*info).dispswitch).multpos THEN BEGIN
				POLYFILL,[(*(*(*info).dataparams).lps)[(*(*info).dataparams).low_low_val],(*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_low_val],(*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_low_val],$
					(*(*(*info).dataparams).lps)[(*(*info).dataparams).low_low_val]], [(*(*info).plotaxes).ls_low_y, (*(*info).plotaxes).ls_low_y, (*(*info).plotaxes).ls_upp_y, (*(*info).plotaxes).ls_upp_y], $
					COLOR = 200
				POLYFILL,[(*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_upp_val],(*(*(*info).dataparams).lps)[(*(*info).dataparams).low_upp_val],(*(*(*info).dataparams).lps)[(*(*info).dataparams).low_upp_val],$
					(*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_upp_val]], [(*(*info).plotaxes).ls_low_y, (*(*info).plotaxes).ls_low_y, (*(*info).plotaxes).ls_upp_y, (*(*info).plotaxes).ls_upp_y], $
					COLOR = 200
				PLOTS, [1,1] * (*(*(*info).dataparams).lps)[(*(*info).dataparams).low_upp_val],[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y], COLOR = (*(*info).plotparams).plotcol, LINESTYLE = 5
				PLOTS, [1,1] * (*(*(*info).dataparams).lps)[(*(*info).dataparams).upp_low_val],[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y], COLOR = (*(*info).plotparams).plotcol, LINESTYLE = 5
			ENDIF
		ENDIF
		ssp = (  *(*(*info).data).loopdata)[FIX((*(*info).dataparams).lx),FIX((*(*info).dataparams).t),*]/REPLICATE((*(*info).dataparams).ms,(*(*info).dataparams).nlp)
		OPLOT, (*(*(*info).dataparams).lps)+(*(*info).dataparams).lp_first, ssp, LINE=2, COLOR = (*(*info).plotparams).plotcol
		OPLOT, (*(*(*info).dataparams).lps),*(*(*info).dataparams).spec,COLOR = (*(*info).plotparams).plotcol
		IF (*(*info).dispswitch).subtract THEN BEGIN
			OPLOT, (*(*(*info).dataparams).lps), *(*(*info).dataparams).spec-ssp, COLOR = (*(*info).plotparams).plotcol
			OPLOT, (*(*(*info).dataparams).lps), *(*(*info).dataparams).spec-ssp, PSYM = 4, COLOR = (*(*info).plotparams).plotcol
		ENDIF
		IF ((*(*info).plotaxes).ls_low_y LT 0.) THEN PLOTS, [((*(*(*info).dataparams).lps))[0],((*(*(*info).dataparams).lps))[(*(*info).dataparams).nlp-1]],[0.,0.], COLOR = (*(*info).plotparams).plotcol
		PLOTS, [1,1] * (*(*(*info).dataparams).lps)[(*(*info).dataparams).lp],[(*(*info).plotaxes).ls_low_y,(*(*info).plotaxes).ls_upp_y], COLOR = (*(*info).plotparams).plotcol
	ENDIF 
END

PRO TANAT_DRAW_CURSCROSS_PLOT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_DRAW_CURSCROSS_PLOT'
	PLOTS, (*(*info).curs).slx,(*(*info).curs).st, /DEVICE, COLOR = !P.COLOR, PSYM = 1
	IF ((*(*info).curs).lockset GE 1) THEN BEGIN
		IF (*(*info).measparams).parabolic_fit THEN BEGIN
			IF ((*(*info).measparams).np LT 3) THEN PLOTS, *(*(*info).measparams).slx_array, *(*(*info).measparams).st_array, /DEVICE, COLOR = !P.COLOR ELSE $
				PLOTS, *(*(*info).measparams).x_vals, *(*(*info).measparams).y_vals, /DEVICE, COLOR = !P.COLOR
			PLOTS, *(*(*info).measparams).slx_array, *(*(*info).measparams).st_array, /DEVICE, COLOR = !P.COLOR, PSYM = 1
			PLOTS, *(*(*info).measparams).slx_array, *(*(*info).measparams).st_array, /DEVICE, COLOR = !P.COLOR, PSYM = 4
		ENDIF ELSE BEGIN
			PLOTS, (*(*info).curs).slxlock, (*(*info).curs).stlock, /DEVICE, COLOR = !P.COLOR, PSYM = 1
			PLOTS, (*(*info).curs).slxlock, (*(*info).curs).stlock, /DEVICE, COLOR = !P.COLOR, PSYM = 4
			PLOTS, [(*(*info).curs).slx,(*(*info).curs).slxlock],[(*(*info).curs).st,(*(*info).curs).stlock], /DEVICE, COLOR = !P.COLOR
			IF ((*(*info).curs).lockset EQ 2) THEN PLOTS, (*(*info).curs).slx, (*(*info).curs).st, /DEVICE, COLOR = !P.COLOR, PSYM = 4
		ENDELSE
	ENDIF
END

PRO TANAT_DRAW_OVERLAY_SAVED_MEASUREMENTS, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_DRAW_OVERLAY_SAVED_MEASUREMENTS'
	(*(*info).dispswitch).overlay_measurements = event.SELECT
	IF (*(*info).dispswitch).overlay_measurements THEN BEGIN
		IF (FILE_INFO((*(*info).measparams).savefilename)).EXISTS THEN BEGIN
			nlines = FILE_LINES((*(*info).measparams).savefilename)
			datarr = STRARR(1,nlines)
			OPENR,unit1,(*(*info).measparams).savefilename,/GET_LUN
			READF,unit1,datarr
			FREE_LUN,unit1
			k=0
			cutfilename = STRMID((*(*info).dataparams).filename,STRPOS((*(*info).dataparams).filename,'/')+1,STRLEN((*(*info).dataparams).filename))
			FOR i=0,nlines-1 DO BEGIN
				lastline_1 = STRSPLIT(datarr(0,i), ' ', /EXTRACT)
				lastline_2 = STRSPLIT(datarr(0,i), '	', /EXTRACT)
				IF (SIZE(lastline_1, /DIMENSIONS) GT 1) THEN splitline = lastline_1 ELSE IF (SIZE(lastline_2, /DIMENSIONS) GT 1) THEN splitline = lastline_2
				IF ((splitline[0] EQ (*(*info).dataparams).filename) OR (splitline[0] EQ cutfilename)) THEN BEGIN
					(*(*(*info).overlays).lxbpoint)[k] = (FLOAT(splitline[3])) * (*(*info).winsizes).windowx / FLOAT((*(*info).dataparams).nlx)
					(*(*(*info).overlays).lxepoint)[k] = (FLOAT(splitline[4])) * (*(*info).winsizes).windowx / FLOAT((*(*info).dataparams).nlx)
					(*(*(*info).overlays).tbpoint)[k] = (FLOAT(splitline[6])); * (*(*info).winsizes).windowy / FLOAT((*(*info).dataparams).nt)
					(*(*(*info).overlays).tepoint)[k] = (FLOAT(splitline[7])); * (*(*info).winsizes).windowy / FLOAT((*(*info).dataparams).nt)
					(*(*(*info).measparams).meas_id)[k] = (STRSPLIT(splitline[2],'.',/EXTRACT))[1]
					k = k+1
				ENDIF
			ENDFOR
		ENDIF ELSE BEGIN
			CD,CURRENT=curpath
			TANAT_WINDOW_OK, event, 'TANAT: ERROR!','TANAT could not overlay saved measurements as '+(*(*info).measparams).savefilename,' does not exist in the current working directory','('+curpath+PATH_SEP()+').',$
				OK_EVENT='TANAT_CLOSE_EVENT_WINDOW', BASE=tlb
			(*(*info).winids).errtlb = tlb
			WIDGET_CONTROL, (*(*info).ctrlsmeas).overlay_button, SET_BUTTON = 0
		ENDELSE
	ENDIF
	TANAT_DRAW, event
END

;================================================================================= TAB EVENT PROCEDURES
PRO TANAT_EVENT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_EVENT'
END

;================================================================================= FILE OPENING PROCEDURES
PRO TANAT_FILE_OPEN, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_FILE_OPEN'
	new_filename = DIALOG_PICKFILE(TITLE='TANAT: Choose a CSAV file', /READ, /MUST_EXIST, FILTER='*.csav',/FIX_FILTER)
	IF (STRCOMPRESS(new_filename) NE '') THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		TANAT_FILE_RESTORE, new_filename, SPECTRUM=spectrum, MS=ms, SPECT_POS=spect_pos, X_LOOP_PTS=x_loop_pts, Y_LOOP_PTS=y_loop_pts, NLX=nlx, NT=nt, NLP=nlp, LOOP_DATA=loop_data, SLAB_SET=slab_set
		*(*(*info).dataparams).spec = spectrum		&	(*(*info).dataparams).ms = ms
		(*(*info).dataparams).nlx = nlx			&	(*(*info).dataparams).nt = nt
		(*(*info).dataparams).nlp = nlp			&	(*(*info).dispswitch).slab_set = slab_set
		*(*(*info).data).loopdata = loop_data		&	*(*(*info).data).loopslice = FLTARR(nlx,nt)
		(*(*info).dataparams).filename = new_filename	&	(*(*info).dataparams).lp = spect_pos
		TANAT_SET_SPECTPARAMS, spectrum, nlp, spect_pos, LINE_CENTER=line_center, LPS=lps, LC=lc, SPXTITLE=spxtitle, V_DOP_VALS=v_dop_vals, V_DOP_SET=v_dop_set, lp_first=lp_first, lp_last=lp_last
		*(*(*info).dataparams).lps = lps		&	(*(*info).dataparams).lc = lc
		(*(*info).plotaxes).spxtitle = spxtitle		&	(*(*info).plotaxes).v_dop = v_dop_vals
		(*(*info).dispswitch).v_dop_set = v_dop_set	&	(*(*info).dataparams).lp_first = lp_first
		(*(*info).dataparams).lp_last = lp_last
		TANAT_SET_TIMESLICE_PARAMS, nlx, nt, x_loop_pts, y_loop_pts, LXDIST=lxdist, LX_FIRST=lx_first, LX_LAST=lx_last, T_FIRST=t_first, T_LAST=t_last
		print,t_first,t_last
		*(*(*info).dataparams).lxdist = lxdist		&	(*(*info).dataparams).lx = FLOAT(lx_first)
		(*(*info).dataparams).t = FLOAT(t_first)	&	(*(*info).curs).slx = FLOAT(lx_first)
		(*(*info).dataparams).t_low = FLOAT(t_first)	&	(*(*info).dataparams).t_upp = FLOAT(t_last)
		(*(*info).dataparams).t_first = FLOAT(t_first)	&	(*(*info).dataparams).t_last = FLOAT(t_last)
		(*(*info).curs).st = FLOAT(t_first)		&	(*(*info).curs).lxlock = FLOAT(lx_first)
		(*(*info).curs).tlock = FLOAT(t_first)		&	(*(*info).curs).slxlock	= FLOAT(lx_first)
		(*(*info).curs).stlock = FLOAT(t_first)
		TANAT_SET_CONTROLS, event
		TANAT_T_RANGE, event
;		TANAT_UPDATE_LP, event
;		TANAT_DRAW, event
		IF (slab_set NE 1) THEN BEGIN
			WSET, (*(*info).winids).ls_drawid
			TV,CONGRID(REPLICATE(200,10,10),(*(*info).winsizes).lswinx,(*(*info).winsizes).lswiny)
			XYOUTS,(*(*info).winsizes).lswinx/2.,(*(*info).winsizes).lswiny/2.,'Could not display spectral information as!Conly one spectral position is available.', COLOR = 0, ALIGNMENT = 0.5, /DEVICE
		ENDIF
		TANAT_SET_SAVEFILENAME, event
	ENDIF
END

PRO TANAT_FILE_RESTORE, filename, SPECTRUM=spectrum, MS=ms, SPECT_POS=spect_pos, X_LOOP_PTS=x_loop_pts, Y_LOOP_PTS=y_loop_pts, NLX=nlx, NT=nt, NLP=nlp, LOOP_DATA=loop_data, SLAB_SET=slab_set
	RESTORE, filename
	spectrum = average_spectrum
	IF (N_ELEMENTS(spectrum) GT 1) THEN BEGIN
		ms	= scaling_factor
		nlp	= (SIZE(spectrum))[1]
	ENDIF ELSE BEGIN
		ms	= 1
		nlp 	= 1
	ENDELSE
	slab_set = 0
	IF (N_ELEMENTS(loop_slab) GT 0) THEN BEGIN
		slab_set = ((SIZE(loop_slab))[0] EQ 3) 
		loop_data = loop_slab
	ENDIF ELSE IF (N_ELEMENTS(loop_slice) GT 0) THEN BEGIN
		slab_set = ((SIZE(loop_slice))[0] EQ 3)
		loop_data = loop_slice
	ENDIF
	nlx = loop_size
	nt = (SIZE(loop_data))[2] 
END

;================================================================================= DETAILED SPECTRUM PROCEDURES
PRO TANAT_LS_LOW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_LS_LOW'
	WIDGET_CONTROL, (*(*info).ctrlsgen).lower_y_text, GET_VALUE = textvalue
	(*(*info).plotaxes).ls_low_y = FLOAT(textvalue[0])
	TANAT_LS_RANGE, event
END

PRO TANAT_LS_UPP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_LS_UPP'
	WIDGET_CONTROL, (*(*info).ctrlsgen).upper_y_text, GET_VALUE = textvalue
	(*(*info).plotaxes).ls_upp_y = FLOAT(textvalue[0])
	TANAT_LS_RANGE, event
END

PRO TANAT_LS_RANGE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_LS_RANGE'
	(*(*info).plotaxes).ls_yrange = (*(*info).plotaxes).ls_upp_y - (*(*info).plotaxes).ls_low_y
	TANAT_DRAW, event
END

PRO TANAT_LS_SUBTRACT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_LS_SUBTRACT'
	(*(*info).dispswitch).subtract = event.SELECT
	TANAT_DRAW, event
END

;================================================================================= PARABOLIC FIT PROCEDURES
PRO TANAT_PARABOLIC_FIT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_PARABOLIC_FIT'
	degree = ((*(*info).measparams).np-1) < 2
	degree = degree > 0
;	s_result = POLY_FIT(FLOAT(*(*(*info).measparams).slx_array), FLOAT(*(*(*info).measparams).st_array), degree)
	s_result = POLY_FIT(FLOAT(*(*(*info).measparams).st_array), FLOAT(*(*(*info).measparams).slx_array), degree)
	act_t_array = *(*(*info).measparams).t_array * (*(*info).measparams).secondststep
	act_lx_array =  149597870. * TAN(FLOAT((*(*info).measparams).arcsecpix) * FLOAT(*(*(*info).measparams).lx_array) / 3600. * !DTOR)
	act_s_result = POLY_FIT(act_t_array,act_lx_array, degree)
	IF ((*(*info).measparams).np GT 2) THEN BEGIN	
		IF ((*(*(*info).measparams).slx_array)[(*(*info).measparams).np-1] LT (*(*(*info).measparams).slx_array)[0]) THEN sign = -1. ELSE sign = 1.
;		*(*(*info).measparams).x_vals = FINDGEN(1000)/1001 * ( (*(*(*info).measparams).slx_array)[(*(*info).measparams).np-1] - (*(*(*info).measparams).slx_array)[0] ) + (*(*(*info).measparams).slx_array)[0]
;		*(*(*info).measparams).y_vals = s_result[2] * (FLOAT(*(*(*info).measparams).x_vals))^2 + s_result[1] * FLOAT(*(*(*info).measparams).x_vals) + s_result[0]
		*(*(*info).measparams).y_vals = FINDGEN(1000)/1001 * ( (*(*(*info).measparams).st_array)[(*(*info).measparams).np-1] - (*(*(*info).measparams).st_array)[0] ) + (*(*(*info).measparams).st_array)[0]
		*(*(*info).measparams).x_vals = s_result[2] * (FLOAT(*(*(*info).measparams).y_vals))^2 + s_result[1] * FLOAT(*(*(*info).measparams).y_vals) + s_result[0]
		(*(*info).measparams).acceleration = 2000. * act_s_result[2] * sign
		WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_text, SET_VALUE = STRTRIM((*(*info).measparams).acceleration,2), /SENSITIVE
	ENDIF
END

PRO TANAT_PARABOLIC_FIT_SET, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_PARABOLIC_FIT_SET'
	(*(*info).measparams).parabolic_fit = event.SELECT
	IF ((*(*info).curs).lockset GT 0) THEN TANAT_RESET_OUTPUTS, event
	IF ((*(*info).measparams).parabolic_fit EQ 0) THEN BEGIN
		(*(*info).measparams).np = 0
		*(*(*info).measparams).lx_array = 0.
		*(*(*info).measparams).t_array = 0.
		*(*(*info).measparams).slx_array = 0.
		*(*(*info).measparams).st_array = 0.
		WIDGET_CONTROL, (*(*info).ctrlsmeas).rem_but, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_label, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_text, SET_VALUE = '0', SENSITIVE = 0
		(*(*info).curs).lockset = 0
		TANAT_DRAW, event
	ENDIF
END

PRO TANAT_PARABOLIC_REM_POINT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_PARABOLIC_REM_POINT'
	(*(*info).measparams).np = (*(*info).measparams).np - 1
	*(*(*info).measparams).lx_array = (*(*(*info).measparams).lx_array)[0:(*(*info).measparams).np-1]
	*(*(*info).measparams).t_array = (*(*(*info).measparams).t_array)[0:(*(*info).measparams).np-1]
	*(*(*info).measparams).slx_array = (*(*(*info).measparams).slx_array)[0:(*(*info).measparams).np-1]
	*(*(*info).measparams).st_array = (*(*(*info).measparams).st_array)[0:(*(*info).measparams).np-1]
	IF ((*(*info).measparams).np LE 2) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsmeas).rem_but, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_label, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_text, SET_VALUE = '0', SENSITIVE = 0
	ENDIF
	IF ((*(*info).measparams).np EQ 3) THEN WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, /SENSITIVE ELSE WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, SENSITIVE = 0
	TANAT_PARABOLIC_FIT, event
	TANAT_DRAW, event
END

;================================================================================= PLAYBACK PROCEDURES
PRO TANAT_PB_BG, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_PB_BG'
	CASE (*(*info).pbparams).spmode OF
		1	: BEGIN
				(*(*info).dataparams).lp += (*(*info).pbparams).spdirection * (*(*info).pbparams).lp_step
				(*(*info).pbparams).spdirection *= -1
				IF ((*(*info).dataparams).lp GT (*(*info).dataparams).lp_last) THEN (*(*info).dataparams).lp -= (*(*info).dataparams).nlp ELSE $
					IF ((*(*info).dataparams).lp LT (*(*info).dataparams).lp_first) THEN (*(*info).dataparams).lp += (*(*info).dataparams).nlp
				WIDGET_CONTROL,(*(*info).ctrlsmeas).lp_slider, SET_VALUE = (*(*info).dataparams).lp
			  END
		ELSE: RETURN
	ENDCASE 
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 1. / (*(*info).pbparams).lp_speed
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_PB_SPECTBLINK, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_PB_SPECTBLINK'
	(*(*info).pbparams).spmode = event.SELECT
	IF (*(*info).pbparams).spmode THEN BEGIN
		(*(*info).pbparams).spmode = 1
		(*(*info).pbparams).spdirection = 1
		WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	ENDIF
END

;================================================================================= COMBINED SPECTRAL POSITIONS PROCEDURES
PRO TANAT_POS_SINGLE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_POS_SINGLE'
	(*(*info).dispswitch).singlepos = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_slider, SENSITIVE = (*(*info).dispswitch).singlepos
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_speed_slider, SENSITIVE = (*(*info).dispswitch).singlepos
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_slider, SENSITIVE = (*(*info).dispswitch).singlepos
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_button, SENSITIVE = (*(*info).dispswitch).singlepos
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_POS_DOUBLE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_POS_DOUBLE'
	(*(*info).dispswitch).doublepos = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_SLIDER_MAX = (*(*info).dataparams).nlp-2, SENSITIVE = (*(*info).dispswitch).doublepos
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_SLIDER_MIN = 2, SENSITIVE = (*(*info).dispswitch).doublepos
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_POS_MULT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_POS_MULT'
	(*(*info).dispswitch).multpos = event.SELECT
	IF ((*(*info).dataparams).low_low_val GE (*(*info).dataparams).upp_low_val) THEN (*(*info).dataparams).upp_low_val = (*(*info).dataparams).low_low_val + 1
	IF ((*(*info).dataparams).upp_upp_val LE (*(*info).dataparams).low_upp_val) THEN (*(*info).dataparams).low_upp_val = (*(*info).dataparams).upp_upp_val - 1
	WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_SLIDER_MAX = (*(*info).dataparams).nlp-4, SENSITIVE = (*(*info).dispswitch).multpos
	WIDGET_CONTROL, (*(*info).ctrlsview).low_upp_slider, SENSITIVE = (*(*info).dispswitch).multpos, SET_VALUE = (*(*info).dataparams).upp_low_val
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_low_slider, SENSITIVE = (*(*info).dispswitch).multpos, SET_VALUE = (*(*info).dataparams).low_upp_val
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_SLIDER_MIN = 3, SENSITIVE = (*(*info).dispswitch).multpos
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

;================================================================================= BMP PROCEDURES
FUNCTION TANAT_READ_BMP, filename, srcdir
; Handles the reading of (button) BMP files
	bmp_dummy = READ_BMP(srcdir+filename)  
	bmp_dummy = TRANSPOSE(bmp_dummy, [1,2,0])
	RETURN, bmp_dummy
END

;================================================================================= RESET PROCEDURES
PRO TANAT_RESET_OUTPUTS, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_RESET_OUTPUTS'
	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_lx_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_t_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_label,  SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_lx_text, SET_VALUE = '0', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).delta_t_text, SET_VALUE = '0', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).speed_text, SET_VALUE = '0', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_params_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_params_text, SET_VALUE = '(0,0)', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_params_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_params_text, SET_VALUE = '(0,0)', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).flag_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).flag_text, SET_VALUE = '0', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).save_button, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_label, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).acc_text, SET_VALUE = '0', SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).rem_but, SENSITIVE = 0
	(*(*info).measparams).flag = 0
	(*(*info).measparams).np = 0
	*(*(*info).measparams).lx_array = 0.
	*(*(*info).measparams).t_array = 0.
	*(*(*info).measparams).slx_array = 0.
	*(*(*info).measparams).st_array = 0.
END

;================================================================================= SAVE MEASUREMENT PROCEDURES
PRO TANAT_SAVE_MEASUREMENT, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SAVE_MEASUREMENT'
	savefiles = FILE_SEARCH((*(*info).measparams).savefilename, COUNT = savefilecount)
	IF savefilecount THEN BEGIN
		nlines = FILE_LINES((*(*info).measparams).savefilename)
		datarr = STRARR(1,nlines)
		OPENR,unit1,(*(*info).measparams).savefilename,/GET_LUN
		READF,unit1,datarr
		FREE_LUN,unit1
		lastline_1 = STRSPLIT(datarr(0,nlines-1), ' ', /EXTRACT)
		lastline_2 = STRSPLIT(datarr(0,nlines-1), '	', /EXTRACT)
		IF (SIZE(lastline_1, /DIMENSIONS) GT 1) THEN lastline = lastline_1 ELSE IF (SIZE(lastline_2, /DIMENSIONS) GT 1) THEN lastline = lastline_2
		lst_id = lastline[2]
		lst_id_parts = STRSPLIT(lst_id,'.',/EXTRACT)
		IF ((*(*info).dataparams).filename EQ lastline[0]) THEN BEGIN
			new_id = '	'+STRING(FLOAT(lst_id)+0.001, FORMAT = '(3(F9.3,x))')
		ENDIF ELSE BEGIN
			prev_main_id = lst_id_parts[0]
			new_id = '	'+STRING(FLOAT(prev_main_id)+1, FORMAT = '(3(F9.3,x))')
		ENDELSE
		OPENU, unit2,(*(*info).measparams).savefilename, WIDTH = 360, /GET_LUN, /APPEND
		IF (*(*info).measparams).parabolic_fit THEN $
			PRINTF, unit2, (*(*info).dataparams).filename, FIX((*(*info).dataparams).lp), new_id, (*(*(*info).measparams).lx_array)[0], (*(*(*info).measparams).lx_array)[1], (*(*(*info).measparams).lx_array)[2], $
			(*(*(*info).measparams).t_array)[0], (*(*(*info).measparams).t_array)[1], (*(*(*info).measparams).t_array)[2], (*(*info).measparams).acceleration, (*(*info).measparams).flag $
		ELSE PRINTF, unit2, (*(*info).dataparams).filename, FIX((*(*info).dataparams).lp), new_id, (*(*info).curs).lxlock, (*(*info).dataparams).lx, (*(*info).measparams).act_delta_lx, (*(*info).curs).tlock, $
			(*(*info).dataparams).t, ABS((*(*info).measparams).delta_t), (*(*info).measparams).speed, (*(*info).measparams).flag
		FREE_LUN,unit2
	ENDIF ELSE BEGIN
		new_id = '	'+STRING(0.0, FORMAT = '(3(F9.3,x))')
		OPENW, unit3, (*(*info).measparams).savefilename, WIDTH = 360, /GET_LUN
		PRINTF, unit3, '#	Filename				lp	id	lx_start		lx_end	delta_lx	t_start	t_end	delta_t	'+$
				'v_abs [km/s]	flag'
		IF (*(*info).measparams).parabolic_fit THEN $
			PRINTF, unit3, (*(*info).dataparams).filename, FIX((*(*info).dataparams).lp), new_id, (*(*(*info).measparams).lx_array)[0], (*(*(*info).measparams).lx_array)[1], (*(*(*info).measparams).lx_array)[2], $
			(*(*(*info).measparams).t_array)[0], (*(*(*info).measparams).t_array)[1], (*(*(*info).measparams).t_array)[2], (*(*info).measparams).acceleration, (*(*info).measparams).flag $
		ELSE PRINTF, unit3, (*(*info).dataparams).filename, FIX((*(*info).dataparams).lp), new_id, (*(*info).curs).lxlock, (*(*info).dataparams).lx, (*(*info).measparams).act_delta_lx, (*(*info).curs).tlock, $
			(*(*info).dataparams).t, ABS((*(*info).measparams).delta_t), (*(*info).measparams).speed, (*(*info).measparams).flag
		FREE_LUN, unit3
	ENDELSE
	PRINT, 'Measurement saved to "'+STRTRIM((*(*info).measparams).savefilename,2)+'" with ID: '+STRTRIM(new_id,2)+', flagged '+STRTRIM((*(*info).measparams).flag,2)
	IF (*(*info).dispswitch).overlay_measurements THEN TANAT_DRAW_OVERLAY_SAVED_MEASUREMENTS, event
END

PRO TANAT_SAVE_MEAUSUREMENT_DEFINE_FLAG, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SAVE_MEASUREMENT_DEFINE_FLAG'
	WIDGET_CONTROL, (*(*info).ctrlsmeas).flag_text, GET_VALUE = textvalue
	(*(*info).measparams).flag = FIX(textvalue[0])
END

;================================================================================= SCALING PROCEDURES
PRO TANAT_SCALING_MAN, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SCALING_MAN'
	(*(*info).dispswitch).man_scale = event.SELECT
	IF (*(*info).dispswitch).man_scale THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsview).max_scale_slider, SET_VALUE = (*(*info).scaling).max_val, /SENSITIVE
		WIDGET_CONTROL, (*(*info).ctrlsview).min_scale_slider, SET_VALUE = (*(*info).scaling).min_val, /SENSITIVE
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsview).max_scale_slider, SET_VALUE = 100, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsview).min_scale_slider, SET_VALUE = 0, SENSITIVE = 0
	ENDELSE
	TANAT_DRAW, event
END		

PRO TANAT_SCALING_RANGE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SCALING_RANGE'
	(*(*info).scaling).minimum = MIN( *(*(*info).data).loopslice )
	maximum = MAX( *(*(*info).data).loopslice )
	(*(*info).scaling).range = maximum - (*(*info).scaling).minimum
END

;================================================================================= SET PROCEDURES
PRO TANAT_SET_ARCSEC, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SET_ARCSEC'
	WIDGET_CONTROL, (*(*info).ctrlssetup).arcsec_text, GET_VALUE = textvalue
	(*(*info).measparams).arcsecpix = FLOAT(textvalue[0])
	IF ((*(*info).curs).lockset GE 1) THEN TANAT_CALCULATE_DELTA, event
END

PRO TANAT_SET_SECONDS, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SET_SECONDS'
	WIDGET_CONTROL, (*(*info).ctrlssetup).seconds_text, GET_VALUE = textvalue
	(*(*info).measparams).secondststep = FLOAT(textvalue[0])
	IF ((*(*info).curs).lockset GE 1) THEN TANAT_CALCULATE_DELTA, event
END

PRO TANAT_SET_SAVEFILENAME, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SET_SAVEFILENAME'
	WIDGET_CONTROL, (*(*info).ctrlssetup).savefile_text, GET_VALUE = textvalue
	(*(*info).measparams).savefilename = textvalue[0]
	savefiles = FILE_SEARCH((*(*info).measparams).savefilename, COUNT = savefilecount)
	compressedfilename = STRCOMPRESS((*(*info).measparams).savefilename, /REMOVE_ALL)
	IF (compressedfilename NE (*(*info).measparams).savefilename) OR ((*(*info).measparams).savefilename EQ '') THEN BEGIN 
		TANAT_WINDOW_OK, event, 'TANAT: ERROR!','Invalid filename. Please enter a filename of','at least one character and without any white spaces.', OK_EVENT='TANAT_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).savnamenulltlb = tlb
	ENDIF ELSE IF savefilecount THEN BEGIN
		TANAT_WINDOW_OK, event,'TANAT: WARNING!', 'The file already exists. New measurements', 'will be added to the existing file.', OK_EVENT='TANAT_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).savnamedoubletlb = tlb
	ENDIF
END

PRO TANAT_SET_CONTROLS, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SET_CONTROLS'
	; First set consequences of slab_set = {0,1}
	IF (*(*info).dispswitch).slab_set THEN BEGIN
		(*(*info).dataparams).low_low_val = 0				&	(*(*info).dataparams).upp_low_val = 1
		(*(*info).dataparams).low_upp_val = (*(*info).dataparams).nlp-2	&	(*(*info).dataparams).upp_upp_val = (*(*info).dataparams).nlp-1
		upp_low_slider_min = 1						&	upp_upp_slider_min = 2
		lp_slider_max = (*(*info).dataparams).lp_last
	ENDIF ELSE BEGIN
		(*(*info).dataparams).low_low_val = 0				& 	(*(*info).dataparams).upp_low_val = 0
		(*(*info).dataparams).low_upp_val = 1				&	(*(*info).dataparams).upp_upp_val = 1
		upp_low_slider_min = 0						&	upp_upp_slider_min = 0
		lp_slider_max = (*(*info).dataparams).lp
	ENDELSE
	(*(*info).dispswitch).overlay_measurements = 0
	(*(*info).dispswitch).smoothed = 0
	; Set controls on tab 1 (Measurements)
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_slider, SET_SLIDER_MIN = (*(*info).dataparams).lp_first, SET_SLIDER_MAX = lp_slider_max, SET_VALUE = (*(*info).dataparams).lp, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_speed_slider, SET_VALUE = 10, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lx_slider, SET_SLIDER_MIN = (*(*info).dataparams).lx_first, SET_SLIDER_MAX = (*(*info).dataparams).lx_last, SET_VALUE = (*(*info).dataparams).lx_first
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_slider, SET_SLIDER_MIN = (*(*info).dataparams).lp_first+1, SET_SLIDER_MAX = (*(*info).dataparams).lp_last, SET_VALUE = (*(*info).dataparams).lp_first+1, $
		SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_button, SET_BUTTON = 0, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_slider, SET_SLIDER_MIN = (*(*info).dataparams).t_first, SET_SLIDER_MAX = (*(*info).dataparams).t_last, SET_VALUE = (*(*info).dataparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlsmeas).overlay_button, SET_BUTTON = (*(*info).dispswitch).overlay_measurements
	; Set controls on tab 2 (View)
	WIDGET_CONTROL, (*(*info).ctrlsview).lower_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_first,2)
	WIDGET_CONTROL, (*(*info).ctrlsview).upper_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_last,2)
	WIDGET_CONTROL, (*(*info).ctrlsview).reset_trange_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).scale_button, SET_BUTTON = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).min_scale_slider, SET_VALUE = 0, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).max_scale_slider, SET_VALUE = 100, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).single_pos, /SET_BUTTON, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsview).double_pos, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsview).mult_pos, SENSITIVE = (*(*info).dispswitch).slab_set
	(*(*info).dispswitch).singlepos = 1
	(*(*info).dispswitch).doublepos = 0
	(*(*info).dispswitch).multpos = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_SLIDER_MIN = (*(*info).dataparams).low_low_val, SET_SLIDER_MAX = (*(*info).dataparams).nlp-2, SET_VALUE = (*(*info).dataparams).low_low_val, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).low_upp_slider, SET_SLIDER_MIN = (*(*info).dataparams).upp_low_val, SET_SLIDER_MAX = (*(*info).dataparams).nlp-3, SET_VALUE = (*(*info).dataparams).upp_low_val, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_low_slider, SET_SLIDER_MIN = upp_low_slider_min, SET_SLIDER_MAX = (*(*info).dataparams).low_upp_val, SET_VALUE = (*(*info).dataparams).low_upp_val, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_SLIDER_MIN = upp_upp_slider_min, SET_SLIDER_MAX = (*(*info).dataparams).upp_upp_val, SET_VALUE = (*(*info).dataparams).upp_upp_val, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsview).non_smooth_button, SET_BUTTON = ABS((*(*info).dispswitch).smoothed-1)
	WIDGET_CONTROL, (*(*info).ctrlsview).smooth_button, SET_BUTTON = (*(*info).dispswitch).smoothed
	; Set general controls
	WIDGET_CONTROL, (*(*info).ctrlsgen).subtract_but, SET_BUTTON = 0, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsgen).lower_y_label, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsgen).lower_y_text, SET_VALUE = '0', SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsgen).upper_y_label, SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsgen).upper_y_text, SET_VALUE = '1', SENSITIVE = (*(*info).dispswitch).slab_set
	WIDGET_CONTROL, (*(*info).ctrlsgen).det_spect, SENSITIVE = (*(*info).dispswitch).slab_set
	; Set to tab 2
	WIDGET_CONTROL, (*(*info).winids).tab_tlb, SET_TAB_CURRENT = 2
END

PRO TANAT_SET_SPECTPARAMS, spectrum, nlp, spect_pos, LINE_CENTER=line_center, lps=lps, LC=lc, SPXTITLE=spxtitle, V_DOP_VALS=v_dop_vals, V_DOP_SET=v_dop_set, lp_first=lp_first, lp_last=lp_last
	lps = 0
	spxtitle = ''
	v_dop_vals = 0
	v_dop_set = 0
	lp_first = 0
;	IF (N_ELEMENTS(spectrum) GT 1) THEN BEGIN
	IF (nlp GT 1) THEN BEGIN
		c_speed	= 2.99792458D5						; Speed of light indeed in km/s
		lps	= FINDGEN(nlp)
		IF (N_ELEMENTS(LINE_CENTER) EQ 0) THEN BEGIN			; If the LINE_CENTER keyword is not set:
			lc	= ( WHERE( spectrum EQ MIN(spectrum) ) )[0]	; autodetermine linecentre value from spectrum
			lps	= ( lps - lps[FIX(lc)] ) 				; reset scale to have lps=0 at lc
			spxtitle= 'Spectral position'
		ENDIF ELSE IF (N_ELEMENTS(LINE_CENTER) EQ 1) THEN BEGIN		; else, if the position is supplied
			lc 	= line_center					; use that given value
			IF (lc GE nlp) OR (lc LT 0) THEN BEGIN			; Check whether 0 LE lc LT nlp
				PRINT,'ERROR: Linecentre index value '+STRTRIM(FIX(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
				RETURN
			ENDIF
			lps	= ( lps - lps[FIX(lc)] ) 				; reset scale to have lps=0 at lc
			spxtitle= 'Spectral position'
		ENDIF ELSE IF (N_ELEMENTS(LINE_CENTER) EQ 2) THEN BEGIN		; else, if also the wavelength is supplied
			lc	 = ( WHERE( spectrum EQ MIN(spectrum) ) )[0]	; autodetermine linecentre value from spectrum
			lambda_c= line_center[0]				; get the linecentre wavelength
			dlambda	= line_center[1]				; and get delta lambda per lineposition
			lps	= ( lps - lps[FIX(lc)] ) 				; reset scale to have lps=0 at lc
			lps	= lps*dlambda + lambda_c				; reset scale to wavelength scale
			spxtitle= 'Wavelength'					; and adapt the xtitle accordingly
			v_dop_set = 1
		ENDIF ELSE BEGIN						; else, if linecentre and wavelength supplied
			lc 	= line_center[0]				; get the linecentre position
			IF (lc GE nlp) OR (lc LT 0) THEN BEGIN			; Check whether 0 LE lc LT nlp
				PRINT,'ERROR: Linecentre index value '+STRTRIM(FIX(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
				RETURN
			ENDIF
			lambda_c= line_center[1]				; and the linecentre wavelength
			dlambda	= line_center[2]				; and get delta(wavelength) per lineposition
			lps	= ( lps - lps[FIX(lc)] ) 				; reset scale to have lps=0 at lc
			lps	= lps*dlambda + lambda_c				; reset scale to wavelength scale
			spxtitle= 'Wavelength'					; and adapt the xtitle accordingly
			v_dop_set = 1
		ENDELSE
	
		IF v_dop_set THEN v_dop_vals = c_speed*(lps/lps[FIX(lc)]-1)		; array with Doppler velocities in km/s
		lp_last = nlp-1
	ENDIF ELSE BEGIN
		lc = 0
		IF (spect_pos GE 2) THEN lp_last = spect_pos ELSE lp_last = spect_pos + 2
	ENDELSE
END

PRO TANAT_SET_TIMESLICE_PARAMS, nlx, nt, x_loop_pts, y_loop_pts, LXDIST=lxdist, LX_FIRST=lx_first, LX_LAST=lx_last, T_FIRST=t_first, T_LAST=t_last
	lxdist = DBLARR(nlx)
	FOR i=0,nlx-2 DO BEGIN
		x_up = x_loop_pts[i+1]
		x_dn = x_loop_pts[i]
		y_up = y_loop_pts[i+1]
		y_dn = y_loop_pts[i]
		lxdist[i] = SQRT( (x_up - x_dn)^2 + (y_up - y_dn)^2 )
	ENDFOR
	lxdist[nlx-1] = MEAN(lxdist[0:nlx-2])
	lx_first	= 0						; Set number of first x-coordinate
	lx_last		= nlx-1						; Set number of last x-coordinate
	t_first		= 0						; Set number of first frame		
	t_last		= nt-1						; Set number of last frame
END

;================================================================================= SLIDER PROCEDURES
PRO TANAT_SLIDER_LOW_LOW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LOW_LOW'
	(*(*info).dataparams).low_low_val = event.VALUE
	IF (*(*info).dispswitch).doublepos THEN BEGIN
		IF ((*(*info).dataparams).low_low_val GE (*(*info).dataparams).upp_upp_val) THEN (*(*info).dataparams).upp_upp_val = (*(*info).dataparams).low_low_val + 1
		WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_VALUE = (*(*info).dataparams).upp_upp_val
	ENDIF ELSE IF (*(*info).dispswitch).multpos THEN BEGIN
		IF ((*(*info).dataparams).low_low_val GE (*(*info).dataparams).upp_low_val) THEN (*(*info).dataparams).upp_low_val = (*(*info).dataparams).low_low_val + 1
		IF ((*(*info).dataparams).upp_low_val GT (*(*info).dataparams).low_upp_val) THEN (*(*info).dataparams).low_upp_val = (*(*info).dataparams).low_upp_val + 1
		IF ((*(*info).dataparams).low_upp_val GE (*(*info).dataparams).upp_upp_val) THEN (*(*info).dataparams).upp_upp_val = (*(*info).dataparams).low_upp_val + 1
		WIDGET_CONTROL, (*(*info).ctrlsview).low_upp_slider, SET_VALUE = (*(*info).dataparams).upp_low_val
		WIDGET_CONTROL, (*(*info).ctrlsview).upp_low_slider, SET_VALUE = (*(*info).dataparams).low_upp_val
		WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_VALUE = (*(*info).dataparams).upp_upp_val
	ENDIF
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_LOW_UPP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LOW_UPP'
	(*(*info).dataparams).upp_low_val = event.VALUE
	IF ((*(*info).dataparams).upp_low_val GT (*(*info).dataparams).low_upp_val) THEN (*(*info).dataparams).low_upp_val = (*(*info).dataparams).low_upp_val + 1
	IF ((*(*info).dataparams).low_upp_val GE (*(*info).dataparams).upp_upp_val) THEN (*(*info).dataparams).upp_upp_val = (*(*info).dataparams).low_upp_val + 1
	IF ((*(*info).dataparams).upp_low_val LE (*(*info).dataparams).low_low_val) THEN (*(*info).dataparams).low_low_val = (*(*info).dataparams).upp_low_val - 1
	WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_VALUE = (*(*info).dataparams).low_low_val
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_VALUE = (*(*info).dataparams).upp_upp_val
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_low_slider, SET_VALUE = (*(*info).dataparams).low_upp_val
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_LP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LP'
	(*(*info).dataparams).lp = event.VALUE
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_LP_INCR, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LP_INCR'
	(*(*info).dataparams).lp += (*(*info).pbparams).lp_step
	IF ((*(*info).dataparams).lp GT (*(*info).dataparams).lp_last) THEN (*(*info).dataparams).lp = (*(*info).dataparams).lp_first
	WIDGET_CONTROL,(*(*info).ctrlsmeas).lp_slider, SET_VALUE = (*(*info).dataparams).lp
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_LP_DECR, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LP_DECR'
	(*(*info).dataparams).lp -= (*(*info).pbparams).lp_step
	IF ((*(*info).dataparams).lp LT (*(*info).dataparams).lp_first) THEN (*(*info).dataparams).lp = (*(*info).dataparams).lp_last
	WIDGET_CONTROL,(*(*info).ctrlsmeas).lp_slider, SET_VALUE = (*(*info).dataparams).lp
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_LX, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_LX'
	(*(*info).dataparams).lx = event.VALUE
	(*(*info).curs).slx = (*(*info).dataparams).lx * (*(*info).winsizes).windowx / (*(*info).dataparams).nlx
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_SCALING_MIN, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_SCALING_MIN'
	(*(*info).scaling).min_val = event.VALUE
	IF ((*(*info).scaling).min_val GE (*(*info).scaling).max_val) THEN BEGIN
		(*(*info).scaling).max_val = (*(*info).scaling).min_val + 1
		WIDGET_CONTROL, (*(*info).ctrlsview).max_scale_slider, SET_VALUE = (*(*info).scaling).max_val
	ENDIF
	TANAT_SCALING_RANGE, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_SCALING_MAX, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_SCALING_MAX'
	(*(*info).scaling).max_val = event.VALUE
	IF ((*(*info).scaling).max_val LE (*(*info).scaling).min_val) THEN BEGIN
		(*(*info).scaling).min_val = (*(*info).scaling).max_val - 1
		WIDGET_CONTROL, (*(*info).ctrlsview).min_scale_slider, SET_VALUE = (*(*info).scaling).min_val
	ENDIF
	TANAT_SCALING_RANGE, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_SPECTSTEP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_SPECTSTEP'
	(*(*info).pbparams).lp_step = event.VALUE
END

PRO TANAT_SLIDER_SPEED, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_SPEED'
	(*(*info).pbparams).lp_speed = event.VALUE
END

PRO TANAT_SLIDER_T, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_T'
	(*(*info).dataparams).t = event.VALUE
	(*(*info).curs).st = (*(*info).dataparams).t * (*(*info).winsizes).windowy / (*(*info).dataparams).nt
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_UPP_LOW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_UPP_LOW'
	(*(*info).dataparams).low_upp_val = event.VALUE
	IF ((*(*info).dataparams).low_upp_val LT (*(*info).dataparams).upp_low_val) THEN (*(*info).dataparams).upp_low_val = (*(*info).dataparams).upp_low_val - 1
	IF ((*(*info).dataparams).upp_low_val LE (*(*info).dataparams).low_low_val) THEN (*(*info).dataparams).low_low_val = (*(*info).dataparams).upp_low_val - 1
	IF ((*(*info).dataparams).low_upp_val GE (*(*info).dataparams).upp_upp_val) THEN (*(*info).dataparams).upp_upp_val = (*(*info).dataparams).low_upp_val + 1
	WIDGET_CONTROL, (*(*info).ctrlsview).low_upp_slider, SET_VALUE = (*(*info).dataparams).upp_low_val
	WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_VALUE = (*(*info).dataparams).low_low_val
	WIDGET_CONTROL, (*(*info).ctrlsview).upp_upp_slider, SET_VALUE = (*(*info).dataparams).upp_upp_val
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END

PRO TANAT_SLIDER_UPP_UPP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SLIDER_UPP_UPP'
	(*(*info).dataparams).upp_upp_val = event.VALUE
	IF (*(*info).dispswitch).doublepos THEN BEGIN
		IF ((*(*info).dataparams).upp_upp_val LE (*(*info).dataparams).low_low_val) THEN (*(*info).dataparams).low_low_val = (*(*info).dataparams).upp_upp_val - 1
		WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_VALUE = (*(*info).dataparams).low_low_val
	ENDIF ELSE IF (*(*info).dispswitch).multpos THEN BEGIN
		IF ((*(*info).dataparams).upp_upp_val LE (*(*info).dataparams).low_upp_val) THEN (*(*info).dataparams).low_upp_val = (*(*info).dataparams).upp_upp_val - 1
		IF ((*(*info).dataparams).low_upp_val LT (*(*info).dataparams).upp_low_val) THEN (*(*info).dataparams).upp_low_val = (*(*info).dataparams).upp_low_val - 1
		IF ((*(*info).dataparams).upp_low_val LE (*(*info).dataparams).low_low_val) THEN (*(*info).dataparams).low_low_val = (*(*info).dataparams).upp_low_val - 1
		WIDGET_CONTROL, (*(*info).ctrlsview).upp_low_slider, SET_VALUE = (*(*info).dataparams).low_upp_val
		WIDGET_CONTROL, (*(*info).ctrlsview).low_upp_slider, SET_VALUE = (*(*info).dataparams).upp_low_val
		WIDGET_CONTROL, (*(*info).ctrlsview).low_low_slider, SET_VALUE = (*(*info).dataparams).low_low_val
	ENDIF
	TANAT_UPDATE_LP, event
	TANAT_DRAW, event
END
		
;================================================================================= DISPLAY PROCEDURES
PRO TANAT_SMOOTH_VIEW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_SMOOTH_VIEW'
	(*(*info).dispswitch).smoothed = event.SELECT
	TANAT_DRAW, event
END

;================================================================================= DISPLAY T PROCEDURES
PRO TANAT_T_LOW, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_T_LOW'
	WIDGET_CONTROL, (*(*info).ctrlsview).lower_t_text, GET_VALUE = textvalue
	(*(*info).dataparams).t_low = FLOAT(textvalue[0])
	IF ((*(*info).dataparams).t_low GE (*(*info).dataparams).t_upp) THEN BEGIN
		(*(*info).dataparams).t_low = (*(*info).dataparams).t_upp - 1
		WIDGET_CONTROL, (*(*info).ctrlsview).lower_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_low,2)
	ENDIF
	IF ((*(*info).dataparams).t_low LT (*(*info).dataparams).t_first) THEN BEGIN
		(*(*info).dataparams).t_low = (*(*info).dataparams).t_first
		WIDGET_CONTROL, (*(*info).ctrlsview).lower_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_low,2)
	ENDIF
	TANAT_T_RANGE, event
END

PRO TANAT_T_UPP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_T_UPP'
	WIDGET_CONTROL, (*(*info).ctrlsview).upper_t_text, GET_VALUE = textvalue
	(*(*info).dataparams).t_upp = FLOAT(textvalue[0])
	IF ((*(*info).dataparams).t_upp LE (*(*info).dataparams).t_low) THEN BEGIN
		(*(*info).dataparams).t_upp = (*(*info).dataparams).t_low + 1
		WIDGET_CONTROL, (*(*info).ctrlsview).upper_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_upp,2)
	ENDIF
	IF ((*(*info).dataparams).t_upp GT (*(*info).dataparams).t_last) THEN BEGIN
		(*(*info).dataparams).t_upp = (*(*info).dataparams).t_last
		WIDGET_CONTROL, (*(*info).ctrlsview).upper_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_upp,2)
	ENDIF
	TANAT_T_RANGE, event
END

PRO TANAT_T_RANGE, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_T_RANGE'
	(*(*info).dataparams).t_range = (*(*info).dataparams).t_upp - (*(*info).dataparams).t_low
	IF ((*(*info).dataparams).t_range+1 NE (*(*info).dataparams).nt) THEN WIDGET_CONTROL, (*(*info).ctrlsview).reset_trange_but, /SENSITIVE ELSE WIDGET_CONTROL, (*(*info).ctrlsview).reset_trange_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_slider, SET_SLIDER_MIN = (*(*info).dataparams).t_low
	WIDGET_CONTROL, (*(*info).ctrlsmeas).t_slider, SET_SLIDER_MAX = (*(*info).dataparams).t_upp
	(*(*info).dataparams).d_nt = (*(*info).dataparams).t_range+1
	(*(*info).curs).st = ((*(*info).dataparams).t-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
	(*(*info).curs).stlock = ((*(*info).curs).tlock-(*(*info).dataparams).t_low) * (*(*info).winsizes).windowy / (*(*info).dataparams).d_nt
	TANAT_UPDATE_LP, event
	IF ((*(*info).curs).lockset GE 1) THEN TANAT_CALCULATE_DELTA, event
	TANAT_DRAW, event
END

PRO TANAT_T_RANGE_RESET, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_T_RANGE_RESET'
	(*(*info).dataparams).t_upp = (*(*info).dataparams).t_last
	(*(*info).dataparams).t_low = (*(*info).dataparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlsview).upper_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_upp,2)
	WIDGET_CONTROL, (*(*info).ctrlsview).lower_t_text, SET_VALUE = STRTRIM((*(*info).dataparams).t_low,2)
	WIDGET_CONTROL, (*(*info).ctrlsview).reset_trange_but, SENSITIVE = 0
	TANAT_T_RANGE, event
END

;================================================================================= UPDATE PROCEDURES
PRO TANAT_UPDATE_LP, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_UPDATE_LP'
	IF (*(*info).dispswitch).singlepos THEN BEGIN
		IF (*(*info).dispswitch).slab_set THEN BEGIN 
			*(*(*info).data).loopslice = REFORM((*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp,(*(*info).dataparams).lp-(*(*info).dataparams).lp_first]) 
		ENDIF ELSE BEGIN
			*(*(*info).data).loopslice = (*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp]
		ENDELSE
	ENDIF ELSE IF (*(*info).dispswitch).doublepos THEN BEGIN
		IF (*(*info).dispswitch).slab_set THEN BEGIN
			loopslice_1 = ((*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp,(*(*info).dataparams).low_low_val-(*(*info).dataparams).lp_first]) 
			loopslice_2 = ((*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp,(*(*info).dataparams).upp_upp_val-(*(*info).dataparams).lp_first])
			*(*(*info).data).loopslice = (loopslice_1 + loopslice_2)/2.
		ENDIF
	ENDIF ELSE IF (*(*info).dispswitch).multpos THEN BEGIN
		loopslice_1 = TOTAL(((*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp,((*(*info).dataparams).low_low_val-(*(*info).dataparams).lp_first):((*(*info).dataparams).upp_low_val-$
			(*(*info).dataparams).lp_first)]),3) / ((*(*info).dataparams).upp_low_val - (*(*info).dataparams).low_low_val + 1)
		loopslice_2 = TOTAL(((*(*(*info).data).loopdata)[*,(*(*info).dataparams).t_low:(*(*info).dataparams).t_upp,((*(*info).dataparams).low_upp_val-(*(*info).dataparams).lp_first):((*(*info).dataparams).upp_upp_val-$
			(*(*info).dataparams).lp_first)]),3) / ((*(*info).dataparams).upp_upp_val - (*(*info).dataparams).low_upp_val + 1)
		*(*(*info).data).loopslice = (loopslice_1 + loopslice_2)/2.
	ENDIF
	TANAT_SCALING_RANGE, event
END

PRO TANAT_UPDATE_STARTUP_FEEDBACK, bgim, xout, yout, feedback_text
	LOADCT,3,/SILENT
	TVSCL,bgim
	LOADCT,0,/SILENT
	FOR i=0,N_ELEMENTS(feedback_text)-1 DO XYOUTS, xout[i], yout[i], feedback_text[i], COLOR=255, /DEVICE, CHARSIZE=1.125
END

;================================================================================= VERBOSE PROCEDURES
PRO TANAT_VERBOSE_GET_ROUTINE, event, rname, IGNORE_LAST=ignore_last
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	prespace = STRJOIN(REPLICATE('  ',SCOPE_LEVEL()-2))
	IF KEYWORD_SET(IGNORE_LAST) THEN (*(*info).feedbparams).last_routine = ''
	IF ((rname NE (*(*info).feedbparams).last_routine) AND ((*(*info).feedbparams).last_routine_count GT 0)) THEN PRINT,''
	IF (rname EQ (*(*info).feedbparams).last_routine) THEN (*(*info).feedbparams).last_routine_count += 1 ELSE (*(*info).feedbparams).last_routine_count = 0
	IF ((*(*info).feedbparams).last_routine_count GT 0) THEN rcount = ' x '+STRTRIM((*(*info).feedbparams).last_routine_count,2)+'.' ELSE rcount = '.'
	IF (rname NE (*(*info).feedbparams).last_routine) THEN PRINT,prespace+'TANAT RUN: Called '+rname+'.' ELSE $
		WRITEU,-1,STRING(FORMAT='(%"\r'+prespace+'TANAT RUN: Called ",a'+STRTRIM(STRLEN(rname),2)+',a'+STRTRIM(STRLEN(rcount),2)+')',rname,rcount) 
	(*(*info).feedbparams).last_routine = rname
END
;================================================================================= WINDOW PROCEDURES
PRO TANAT_WINDOW, xsize, ysize, leader, title, base, wid, xoffset, yoffset, DRAWID = drawid, DRAWBASE =disp, XSCROLL = xscrollsize, YSCROLL = yscrollsize, SCROLL = scroll
	base = WIDGET_BASE(TITLE = STRTRIM(title), GROUP_LEADER = leader, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	drawid = WIDGET_DRAW(disp, XSIZE = xsize, YSIZE = ysize, RETAIN = 2)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET=xoffset, TLB_SET_YOFFSET=yoffset
	WIDGET_CONTROL, drawid, GET_VALUE = wid
END

PRO TANAT_WINDOW_OK, event, title, message1, message2, message3, message4, OK_EVENT=ok_event, BASE=base, BLOCK=block
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF ((((*(*info).feedbparams).verbosity)[2] EQ 1) OR (((*(*info).feedbparams).verbosity)[3] EQ 1)) THEN TANAT_VERBOSE_GET_ROUTINE, event, 'TANAT_WINDOW_OK'
	base = WIDGET_BASE(TITLE = title, GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	text_label1 = WIDGET_LABEL(message_base, VALUE = message1)
	IF (N_ELEMENTS(message2) GT 0) THEN text_label2 = WIDGET_LABEL(message_base, VALUE = message2)
	IF (N_ELEMENTS(message3) GT 0) THEN text_label3 = WIDGET_LABEL(message_base, VALUE = message3)
	IF (N_ELEMENTS(message4) GT 0) THEN text_label4 = WIDGET_LABEL(message_base, VALUE = message4)
	button_base = WIDGET_BASE(disp,/ROW,/ALIGN_CENTER)
	ok_but = WIDGET_BUTTON(button_base, VALUE = 'OK' , EVENT_PRO = ok_event)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	IF (N_ELEMENTS(BLOCK) NE 1) THEN block = 0
	XMANAGER, 'TANAT', base, NO_BLOCK=ABS(block-1)
END

;========================================================================= MAIN PROGRAM CODE
PRO TANAT,$							; call program
	filename,$						; name of file containing time slab
	LINE_CENTER=line_center, $				; line center keyword with spectral information
	ASECPIX=asecpix, $					; spatial resolution in arcseconds per pixel
	DT=dt, $						; time step
	VERBOSE=verbose
	
;========================================================================= PROGRAM-INFO ON CALL W/O PARAMS
	IF N_PARAMS() LT 1 THEN BEGIN
		PRINT,'TANAT, filename, LINE_CENTER=line_center, ASECPIX=asecpix, DT=dt, VERBOSE=verbose'
		RETURN
	ENDIF

;================================================================================= VERSION AND REVISION NUMBER
	version_number = '1.0'
	revision_number = '59'

;================================================================================= PROGRAM VERBOSITY CHECK
	IF (N_ELEMENTS(VERBOSE) NE 1) THEN BEGIN			
		IF (N_ELEMENTS(VERBOSE) GT 1) THEN PRINT,'ERROR: The VERBOSE keyword may only be set to a single integer number. Reverting to default verbosity level 0.'
		verbose = 0
		verbosity = [0,0,0,0,0]
	ENDIF ELSE BEGIN
		verbose >= 0	&	verbose <= 26
	ENDELSE
	verbosity = TANAT_DEC2BIN(verbose)


;================================================================================= TANAT DIRECTORY CHECK
	file_tanat		= (ROUTINE_INFO('TANAT',/SOURCE)).PATH
	dir_aux 		= FILE_DIRNAME(file_tanat,/MARK_DIRECTORY)
	dir_resources		= STRMID(dir_aux,0,STRPOS(STRMID(dir_aux,0,STRLEN(dir_aux)-1),'/',/REVERSE_SEARCH)+1)+'resources'+PATH_SEP()

;================================================================================= START-UP WINDOW
	screensize 	= GET_SCREEN_SIZE()											; Get the user screensize
	x_screen_mid	= screensize[0]/2.
	y_screen_mid	= screensize[1]/2.
	startup_im 	= REBIN(REFORM(TOTAL((TANAT_READ_BMP('tanat_startup.bmp',dir_resources))[*,*,1:2],3)),400,300)
	startup_nx 	= (SIZE(startup_im))[1]
	startup_ny 	= (SIZE(startup_im))[2]
	startup_xpos 	= FIX(x_screen_mid-startup_nx/2.)
	startup_ypos 	= FIX(y_screen_mid-startup_ny/2.)
	xout 		= REPLICATE(24,9)
	yout 		= REPLICATE(FIX(startup_ny/2.5)+10,9)-INDGEN(9)*15

;========================================================================= READ-IN OF FILE
	TANAT_FILE_RESTORE, filename, SPECTRUM=spectrum, MS=ms, SPECT_POS=spect_pos, X_LOOP_PTS=x_loop_pts, Y_LOOP_PTS=y_loop_pts, NLX=nlx, NT=nt, NLP=nlp, LOOP_DATA=loop_data, SLAB_SET=slab_set
	filename = STRMID(filename, STRPOS(filename,'/',/REVERSE_SEARCH)+1,STRLEN(filename)-1)		; Failsafe to remove leading /

	IF (N_ELEMENTS(type) GT 0) THEN BEGIN
		IF (type EQ 1) THEN PRINT,'WARNING: you are currently analysing data from an approximated loop slice!'
	ENDIF

;========================================================================= SETTING START-UP OPTIONS 
;--------------------------------------------------------------------------------- INITIAL SPECTRAL PARAMETERS
	TANAT_SET_SPECTPARAMS, spectrum, nlp, spect_pos, LINE_CENTER=line_center, lps=lps, LC=lc, SPXTITLE=spxtitle, V_DOP_VALS=v_dop_vals, V_DOP_SET=v_dop_set, lp_first=lp_first, lp_last=lp_last
	lp_start = spect_pos

;--------------------------------------------------------------------------------- INITIAL TIMESLICE PARAMETERS
	TANAT_SET_TIMESLICE_PARAMS, nlx, nt, x_loop_pts, y_loop_pts, LXDIST=lxdist, LX_FIRST=lx_first, LX_LAST=lx_last, T_FIRST=t_first, T_LAST=t_last
	IF (N_ELEMENTS(ASECPIX) EQ 1) THEN arcsecpix = asecpix ELSE arcsecpix = 0.0592
	IF (N_ELEMENTS(DT) EQ 1) THEN secondststep = dt ELSE secondststep = 1.

;--------------------------------------------------------------------------------- WINDOW SIZES (CHANGE ONLY
;--------------------------------------------------------------------------------- NUMERICAL VALUES!)
	imwinx 		= 0.5 * screensize[0]				; Set maximum x-extent of image window
	imwiny 		= 0.85 * screensize[1]				; Set maximum y-extent of image window

	windowx		= 0.2 * screensize[0]				; Set maximum x-extent of spectral win
	windowy		= imwiny					; Set maximum y-extent of spectral win
	lswinx 		= 0.25 * screensize[0]				; Set maximum x-extent of loc spec win
	lswiny 		= 0.2 * screensize[0]				; Set maximum y-extent of loc spec win

	lsx0		= 0.1						; x0 coordinate of the detailed spectrum plot
	lsy0		= 0.1						; y0 coordinate of the plot
	lsx1		= 0.95						; x1 coordinate of the plot
	IF (v_dop_set EQ 1) THEN lsy1 = 0.9 ELSE lsy1 = 0.95		; y1 coordinate of the plot

	spx0		= 0.2						; x0 coordinate of the temporal spectrum plot
	spy0		= 0.1						; y0 coordinate of the plot
	spx1		= 0.9						; x1 coordinate of the plot
	spy1		= 0.95						; y1 coordinate of the plot
	xplspw		= spx1 - spx0					; x-extent of the plot
	yplspw		= spy1 - spy0					; y-extent of the plot

	ntreb		= yplspw * windowy				; actual nt rebinning factor
	spwiny		= ntreb / yplspw				; extent and rebinning factor
	nlxreb		= xplspw * windowx 				; actual nl rebinning factor
	spwinx		= nlxreb / xplspw				; determine new window x- and y-size from plot

	xdelta		= 20						; Extra xoffset for positioning of windows
	ydelta		= 40						; Extra yoffset for positioning of windows

	ls_low_y	= 0.
	ls_upp_y	= 1.
	ls_yrange	= ls_upp_y - ls_low_y

	loopdata = PTR_NEW(loop_data, /NO_COPY)
	loopslice = PTR_NEW(FLTARR(nlx,nt))

	lxbpoint = PTR_NEW(FLTARR(1000))
	lxepoint = PTR_NEW(FLTARR(1000))
	tbpoint = PTR_NEW(FLTARR(1000))
	tepoint = PTR_NEW(FLTARR(1000))
	meas_id = PTR_NEW(STRARR(1000))
	lx_array = PTR_NEW(FLTARR(30))
	t_array = PTR_NEW(FLTARR(30))
	slx_array = PTR_NEW(FLTARR(30))
	st_array = PTR_NEW(FLTARR(30))
	x_vals = PTR_NEW(FLTARR(1000))
	y_vals = PTR_NEW(FLTARR(1000))

;========================================================================= SETTING UP WIDGET
;--------------------------------------------------------------------------------- INITIALISE CONTROL PANEL
	control_panel	= WIDGET_BASE(TITLE = 'TANAT: Timeslice Analysis Tool', TLB_FRAME_ATTR = 1, /ROW, KILL_NOTIFY = 'TANAT_CLEANUP', MBAR = menubar)
	filemenu	= WIDGET_BUTTON(menubar, VALUE = 'File', /MENU, UVALUE = 'file')
	about		= WIDGET_BUTTON(filemenu, VALUE = 'About', EVENT_PRO = 'TANAT_ABOUT_WINDOW', ACCELERATOR='Ctrl+A')
	open_file	= WIDGET_BUTTON(filemenu, VALUE = 'Open...', EVENT_PRO = 'TANAT_FILE_OPEN', /SEPARATOR, ACCELERATOR = 'Ctrl+O')
	exitmenu	= WIDGET_BUTTON(filemenu, VALUE = 'Quit', EVENT_PRO = 'TANAT_CLOSE', /SEPARATOR, ACCELERATOR='Ctrl+Q')

	shortcutmenu		= WIDGET_BUTTON(menubar, VALUE = 'Control shortcuts', /MENU, UVALUE = 'shortcut')
	sh_spectralmenu 	= WIDGET_BUTTON(shortcutmenu, VALUE = 'Spectral options', /MENU, UVALUE = 'spectral')
	sh_lp_incr_button 	= WIDGET_BUTTON(sh_spectralmenu, VALUE = 'Spectral position +', EVENT_PRO = 'TANAT_SLIDER_LP_INCR', ACCELERATOR = 'Shift+S')
	sh_lp_decr_button 	= WIDGET_BUTTON(sh_spectralmenu, VALUE = 'Spectral position -', EVENT_PRO = 'TANAT_SLIDER_LP_DECR', ACCELERATOR = 'Shift+A')

	buttons_base	= WIDGET_BASE(control_panel, /COLUMN)
	tab_tlb 	= WIDGET_TAB(buttons_base, LOCATION=location)

	measure_tab	= WIDGET_BASE(tab_tlb, TITLE = 'Measurements',/COLUMN)
	view_tab	= WIDGET_BASE(tab_tlb, TITLE = 'View',/COLUMN)
	setup_tab	= WIDGET_BASE(tab_tlb, TITLE='Set up',/COLUMN)
	WIDGET_CONTROL, tab_tlb, SET_TAB_CURRENT = 2

	settings	= WIDGET_BASE(setup_tab, /ROW, /FRAME)
	dimensions	= WIDGET_BASE(settings ,/COLUMN)
	arcsec_field	= WIDGET_BASE(dimensions, /ROW)
	arcsec_label	= WIDGET_LABEL(arcsec_field, VALUE = 'Arcseconds per pixel:', /ALIGN_LEFT)
	arcsec_text	= WIDGET_TEXT(arcsec_field, VALUE = STRTRIM(arcsecpix,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_SET_ARCSEC')
	seconds_field	= WIDGET_BASE(dimensions, /ROW)
	seconds_label	= WIDGET_LABEL(seconds_field, VALUE = 'Seconds per timestep:', /ALIGN_LEFT)
	seconds_text	= WIDGET_TEXT(seconds_field, VALUE = STRTRIM(secondststep,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_SET_SECONDS')
	savefile	= WIDGET_BASE(settings, /COLUMN, /FRAME)
	savefile_label	= WIDGET_LABEL(savefile, VALUE = 'Measurements saved to file:')
	savefile_text	= WIDGET_TEXT(savefile, VALUE = 'tanat_measurements.dat', /EDITABLE, XSIZE = 25, EVENT_PRO = 'TANAT_SET_SAVEFILENAME')

	t_range_field	= WIDGET_BASE(view_tab, /ROW, /FRAME)
	lower_t_label	= WIDGET_LABEL(t_range_field, VALUE = 'Lower t-value:', /ALIGN_LEFT)
	lower_t_text	= WIDGET_TEXT(t_range_field, VALUE = STRTRIM(t_first,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_T_LOW')
	upper_t_label	= WIDGET_LABEL(t_range_field, VALUE = 'Upper t-value:', /ALIGN_LEFT)
	upper_t_text	= WIDGET_TEXT(t_range_field, VALUE = STRTRIM(t_last,2),  /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_T_UPP')
	reset_trange_but= WIDGET_BUTTON(t_range_field, VALUE = 'Reset temporal boundaries', EVENT_PRO = 'TANAT_T_RANGE_RESET', SENSITIVE = 0)

	sliders		= WIDGET_BASE(measure_tab, /ROW, /FRAME)
	pos_sliders	= WIDGET_BASE(sliders, /COLUMN)
	lp_slid		= WIDGET_SLIDER(pos_sliders, TITLE = 'Spectral position', MIN = lp_first, MAX = lp_last, VALUE = lp_start, EVENT_PRO = 'TANAT_SLIDER_LP', /DRAG, SENSITIVE = slab_set)
	lx_slider	= WIDGET_SLIDER(pos_sliders, TITLE = 'Pixel position along the loop', MIN = lx_first, MAX = lx_last, VALUE = lx_first, EVENT_PRO = 'TANAT_SLIDER_LX', /DRAG)
	t_slider	= WIDGET_SLIDER(pos_sliders, TITLE = 'Frame number', MIN = t_first, MAX = t_last, VALUE = t_first, EVENT_PRO = 'TANAT_SLIDER_T', /DRAG)

	blink_sliders	= WIDGET_BASE(sliders,/COLUMN)
	lp_speed_slid	= WIDGET_SLIDER(blink_sliders, TITLE = 'Animation speed [blink/s]', MIN = 1, MAX = 100, VALUE = 10, EVENT_PRO = 'TANAT_SLIDER_SPEED', /DRAG, SENSITIVE = slab_set)
	lp_blink_slid	= WIDGET_SLIDER(blink_sliders, TITLE = 'Spectral increment', MIN = lp_first+1, MAX = lp_last, EVENT_PRO = 'TANAT_SLIDER_SPECTSTEP',/DRAG, SENSITIVE = slab_set)
	lp_blink_field	= WIDGET_BASE(blink_sliders, /ROW,/NONEXCLUSIVE)
	lp_blink_but	= WIDGET_BUTTON(lp_blink_field, VALUE = 'Blink between spectral positions', EVENT_PRO = 'TANAT_PB_SPECTBLINK', SENSITIVE = slab_set)

	parameters	= WIDGET_BASE(measure_tab, /ROW, /FRAME)
	speed		= WIDGET_BASE(parameters, /COLUMN)
	delta_lx	= WIDGET_BASE(speed, /ROW)
	delta_t		= WIDGET_BASE(speed, /ROW)
	delta_speed	= WIDGET_BASE(speed, /ROW)
	acc		= WIDGET_BASE(speed, /ROW)
	delta_lx_label	= WIDGET_LABEL(delta_lx, VALUE = 'Travelled distance [pixel]:', /ALIGN_LEFT, SENSITIVE = 0)
	delta_lx_text	= WIDGET_LABEL(delta_lx, VALUE = '0', SENSITIVE = 0, /DYNAMIC_RESIZE)
	delta_t_label	= WIDGET_LABEL(delta_t, VALUE = 'Elapsed time [frame]:', /ALIGN_LEFT, SENSITIVE = 0)
	delta_t_text	= WIDGET_LABEL(delta_t, VALUE = '0',  SENSITIVE = 0, /DYNAMIC_RESIZE)
	speed_label	= WIDGET_LABEL(delta_speed, VALUE = 'Absolute speed [km/s]:', /ALIGN_LEFT, SENSITIVE = 0)
	speed_text	= WIDGET_LABEL(delta_speed, VALUE = '0', SENSITIVE = 0, /DYNAMIC_RESIZE)
	acc_label	= WIDGET_LABEL(acc, VALUE = 'Acceleration [m/s^2]:', /ALIGN_LEFT, SENSITIVE = 0)
	acc_text	= WIDGET_LABEL(acc, VALUE = '0', SENSITIVe = 0, /DYNAMIC_RESIZE)

	params_and_save = WIDGET_BASE(parameters, /COLUMN)
	lx_params	= WIDGET_BASE(params_and_save, /ROW)
	t_params	= WIDGET_BASE(params_and_save, /ROW)
	save_params	= WIDGET_BASE(params_and_save, /ROW)
	lx_params_label	= WIDGET_LABEL(lx_params, VALUE = 'Position coordinates:', /ALIGN_LEFT, SENSITIVE = 0)
	lx_params_text	= WIDGET_LABEL(lx_params, VALUE = '(0,0)', SENSITIVE = 0, /DYNAMIC_RESIZE)
	t_params_label	= WIDGET_LABEL(t_params, VALUE = 'Time coordinates:', /ALIGN_LEFT, SENSITIVE = 0)
	t_params_text	= WIDGET_LABEL(t_params, VALUE = '(0,0)', SENSITIVE = 0, /DYNAMIC_RESIZE)
	flag_label	= WIDGET_LABEL(save_params, VALUE = 'Flag:', SENSITIVE = 0, /ALIGN_LEFT)
	flag_text	= WIDGET_TEXT(save_params, VALUE = '0', SENSITIVE = 0, XSIZE = 3, /EDITABLE, EVENT_PRO = 'TANAT_SAVE_MEASUREMENT_DEFINE_FLAG')
	save_button	= WIDGET_BUTTON(save_params, VALUE = 'Save measurement', EVENT_PRO = 'TANAT_SAVE_MEASUREMENT', SENSITIVE = 0)

	parabolic_fit	= WIDGET_BASE(measure_tab, /ROW, /FRAME)
	fit_but_base	= WIDGET_BASE(parabolic_fit, /ROW, /NONEXCLUSIVE)
	set_fit_but	= WIDGET_BUTTON(fit_but_base, VALUE = 'Store points for parabolic fit', EVENT_PRO = 'TANAT_PARABOLIC_FIT_SET')
	rem_but		= WIDGET_BUTTON(parabolic_fit, VALUE = 'Remove last point', EVENT_PRO = 'TANAT_PARABOLIC_REM_POINT', SENSITIVE = 0)

	overlay		= WIDGET_BASE(measure_tab, /ROW, /FRAME, /NONEXCLUSIVE)
	overlay_but	= WIDGET_BUTTON(overlay, VALUE = 'Overlay saved measurements for this timeslice', EVENT_PRO = 'TANAT_DRAW_OVERLAY_SAVED_MEASUREMENTS')


	detspect_frame	= WIDGET_BASE(buttons_base, /FRAME, /COLUMN)
	detspect_buts	= WIDGET_BASE(detspect_frame, /ROW, /NONEXCLUSIVE)
	subtract_but	= WIDGET_BUTTON(detspect_buts, VALUE = 'Subtract average', EVENT_PRO = 'TANAT_LS_SUBTRACT', TOOLTIP = 'Subtract detailed spectrum from average spectrum', SENSITIVE = slab_set)
	detspect_range	= WIDGET_BASE(detspect_frame, /ROW)
	lower_y_label	= WIDGET_LABEL(detspect_range, VALUE = 'Lower y-value:', /ALIGN_LEFT, SENSITIVE = slab_set)
	lower_y_text	= WIDGET_TEXT(detspect_range, VALUE = '0',  /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_LS_LOW', SENSITIVE = slab_set)
	upper_y_label	= WIDGET_LABEL(detspect_range, VALUE = 'Upper y-value:', /ALIGN_LEFT, SENSITIVE = slab_set)
	upper_y_text	= WIDGET_TEXT(detspect_range, VALUE = '1',  /EDITABLE, XSIZE = 5, EVENT_PRO = 'TANAT_LS_UPP', SENSITIVE = slab_set)
	detsp_drawbase	= WIDGET_BASE(detspect_frame, /COLUMN)
	det_spect	= WIDGET_DRAW(detsp_drawbase, XSIZE = lswinx, YSIZE = lswiny, SENSITIVE = slab_set)
	
	scale_base	= WIDGET_BASE(view_tab, /ROW, /FRAME)
	scale_but_base	= WIDGET_BASE(scale_base, /ROW, /NONEXCLUSIVE)
	scale_but	= WIDGET_BUTTON(scale_but_base, VALUE = 'Manual scaling', EVENT_PRO = 'TANAT_SCALING_MAN')
	min_scale_slider= WIDGET_SLIDER(scale_base, TITLE = 'Minimum', VALUE = 0, MIN = 0, MAX = 99, EVENT_PRO = 'TANAT_SLIDER_SCALING_MIN', /DRAG, SENSITIVE = 0)
	max_scale_slider= WIDGET_SLIDER(scale_base, TITLE = 'Maximum', VALUE = 100, MIN = 1, MAX = 100,  EVENT_PRO = 'TANAT_SLIDER_SCALING_MAX', /DRAG, SENSITIVE = 0)

	scale_view_base = WIDGET_BASE(view_tab, /COLUMN, /FRAME)
	pos_buts_1	= WIDGET_BASE(scale_view_base, /COLUMN, /EXCLUSIVE)
	single_pos	= WIDGET_BUTTON(pos_buts_1, VALUE = 'Single spectral position', EVENT_PRO = 'TANAT_POS_SINGLE', SENSITIVE = slab_set)
	WIDGET_CONTROL, single_pos, /SET_BUTTON
	double_pos	= WIDGET_BUTTON(pos_buts_1, VALUE = 'Two combined spectral positions', EVENT_PRO = 'TANAT_POS_DOUBLE', SENSITIVE = slab_set)
	mult_pos	= WIDGET_BUTTON(pos_buts_1, VALUE = 'More combined spectral positions', EVENT_PRO = 'TANAT_POS_MULT', SENSITIVE = slab_set)

	lower_sliders	= WIDGET_BASE(scale_view_base, /ROW)
	low_low_slider	= WIDGET_SLIDER(lower_sliders, TITLE = 'Lower lower spectral position', VALUE = 0, MIN = 0, MAX = nlp-3, SENSITIVE = 0, EVENT_PRO = 'TANAT_SLIDER_LOW_LOW', /DRAG)
	low_upp_slider	= WIDGET_SLIDER(lower_sliders, TITLE = 'Upper lower spectral position', VALUE = 1, MIN = 1, MAX = nlp-2, SENSITIVE = 0, EVENT_PRO = 'TANAT_SLIDER_LOW_UPP', /DRAG)
	upper_sliders	= WIDGET_BASE(scale_view_base, /ROW)
	upp_low_slider	= WIDGET_SLIDER(upper_sliders, TITLE = 'Lower upper spectral position', VALUE = nlp-2, MIN = 2, MAX = nlp-2, SENSITIVE = 0, EVENT_PRO = 'TANAT_SLIDER_UPP_LOW', /DRAG)
	upp_upp_slider	= WIDGET_SLIDER(upper_sliders, TITLE = 'Upper upper spectral position', VALUE = nlp-1, MIN = 3, MAX = nlp-1, SENSITIVE = 0, EVENT_PRO = 'TANAT_SLIDER_UPP_UPP', /DRAG)

	smooth_buttons	= WIDGET_BASE(view_tab, /ROW, /EXCLUSIVE, /FRAME)
	non_smooth_but	= WIDGET_BUTTON(smooth_buttons, VALUE = 'Non-smoothed view')
	smooth_but	= WIDGET_BUTTON(smooth_buttons, VALUE = 'Smoothed view', EVENT_PRO = 'TANAT_SMOOTH_VIEW')
	WIDGET_CONTROL, non_smooth_but, /SET_BUTTON

	draw_base	= WIDGET_BASE(control_panel, /COLUMN)
	timeslice	= WIDGET_DRAW(draw_base, XSIZE = windowx, YSIZE = windowy)

	WIDGET_CONTROL, control_panel, /REALIZE
	bg = WIDGET_BASE(control_panel, EVENT_PRO = 'TANAT_PB_BG')

	WIDGET_CONTROL, timeslice, GET_VALUE = drawid
	WIDGET_CONTROL, det_spect, GET_VALUE = ls_drawid
	WIDGET_CONTROL, timeslice, EVENT_PRO = 'TANAT_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS

;--------------------------------------------------------------------------------- MEASUREMENT TAB CONTROLS
	ctrlsgen = { $
		subtract_but:subtract_but, lower_y_label:lower_y_label, lower_y_text:lower_y_text, upper_y_label:upper_y_label, $
		upper_y_text:upper_y_text, det_spect:det_spect $
	}
;--------------------------------------------------------------------------------- MEASUREMENT TAB CONTROLS
	ctrlsmeas = { $
		lp_slider:lp_slid, lx_slider:lx_slider, t_slider:t_slider, $
		lp_speed_slider:lp_speed_slid, lp_blink_slider:lp_blink_slid, lp_blink_button:lp_blink_but, $
		delta_lx_label:delta_lx_label, delta_t_label:delta_t_label, delta_lx_text:delta_lx_text, delta_t_text:delta_t_text, $
		speed_label:speed_label, speed_text:speed_text, acc_label:acc_label, acc_text:acc_text, $
		lx_params_text:lx_params_text, t_params_text:t_params_text, lx_params_label:lx_params_label, t_params_label:t_params_label, $
		flag_label:flag_label, flag_text:flag_text, save_button:save_button, rem_but:rem_but, overlay_button:overlay_but $
	}
;--------------------------------------------------------------------------------- SETUP CONTROLS
	ctrlssetup = { $
		arcsec_text:arcsec_text, seconds_text:seconds_text, $
		savefile_text:savefile_text $
	}
;--------------------------------------------------------------------------------- VIEW TAB CONTROLS
	ctrlsview = { $
		lower_t_text:lower_t_text, upper_t_text:upper_t_text, reset_trange_but:reset_trange_but, $
		min_scale_slider:min_scale_slider, max_scale_slider:max_scale_slider, scale_button:scale_but, $
		single_pos:single_pos, double_pos:double_pos, mult_pos:mult_pos, $
		low_low_slider:low_low_slider, low_upp_slider:low_upp_slider, upp_low_slider:upp_low_slider, upp_upp_slider:upp_upp_slider, $
		non_smooth_button:non_smooth_but, smooth_button:smooth_but $
	}
;--------------------------------------------------------------------------------- CURSOR
	curs = { $
		lxlock:FLOAT(lx_first), tlock:FLOAT(t_first), slxlock:FLOAT(lx_first), stlock:FLOAT(t_first), $
		slx:FLOAT(lx_first), st:FLOAT(t_first), lockset:0 $
	}
;--------------------------------------------------------------------------------- DATA 
	data = { $
		loopdata:loopdata, loopslice:loopslice $
	}
;--------------------------------------------------------------------------------- DATA PARAMETERS
	dataparams = { $
		filename:filename, spec:PTR_NEW(spectrum), lps:PTR_NEW(lps), ms:ms, nlx:nlx, nt:nt, nlp:nlp, $
		lxdist:PTR_NEW(lxdist), lx:FLOAT(lx_first), t:FLOAT(t_first), lx_first:lx_first, lx_last:lx_last, $
		t_first:t_first, t_last:t_last, lp:lp_start, lc:lc, lp_first:lp_first, lp_last:lp_last, $
		t_low:t_first, t_upp:t_last, t_range:nt, d_nt:nt, low_low_val:0, upp_low_val:1, low_upp_val:nlp-2, upp_upp_val:nlp-1$
	}
;--------------------------------------------------------------------------------- DISPLAY SWITCHES
	dispswitch = { $
		overlay_measurements:0, singlepos:1, doublepos:0, multpos:0, subtract:0, $
		man_scale:0, smoothed:0, v_dop_set:v_dop_set, slab_set:slab_set $
	}
;--------------------------------------------------------------------------------- FEEDBACK PARAMS
	feedbparams = { $
		xout:xout, yout:yout, startup_im:startup_im, verbosity:verbosity, last_routine:'', last_routine_count:0 $
	}
;--------------------------------------------------------------------------------- SAVING MEASUREMENT PARAMETERS
	measparams = { $
		lx_array:lx_array, t_array:t_array, slx_array:slx_array, st_array:st_array, np:0, speed:0., delta_lx:0, act_delta_lx:0., delta_t:0, $
		acceleration:0., parabolic_fit:0, x_vals:x_vals, y_vals:y_vals, meas_id:meas_id, secondststep:secondststep, arcsecpix:arcsecpix, $
		savefilename:'tanat_measurements.dat', flag:0 $
	}
;--------------------------------------------------------------------------------- OVERLAY PARAMETERS
	overlays = { $
		lxbpoint:lxbpoint, lxepoint:lxepoint, tbpoint:tbpoint, tepoint:tepoint $
	}
;--------------------------------------------------------------------------------- PLAYBACK PARAMS
	pbparams = { $
		spmode:0, spdirection:1, lp_step:1, lp_speed:10, bg:bg $
	}
;--------------------------------------------------------------------------------- PLOTAXES
	plotaxes = { $
		ls_low_y:ls_low_y, ls_upp_y:ls_upp_y, ls_yrange:ls_yrange, v_dop:v_dop_vals, $
		spxtitle:spxtitle, spwinx:spwinx $
	}
;--------------------------------------------------------------------------------- PLOT PARAMETERS
	plotparams = { $
		plotcol:0, bgplotcol:!P.COLOR $
	}
;--------------------------------------------------------------------------------- PLOT POSITION
	plotpos = { $
		lsx0:lsx0, lsy0:lsy0, lsx1:lsx1, lsy1:lsy1 $
	}
;--------------------------------------------------------------------------------- SCALING PARAMS
	scaling = { $
		min_val:0., max_val:100., range:100.,minimum:0. $
	}
;--------------------------------------------------------------------------------- VERSION INFO
	versioninfo = { $
		version_number:version_number, revision_number:revision_number $				
	}
;--------------------------------------------------------------------------------- WINDOW IDs
	winids = { $
		root:control_panel, wid:drawid, ls_drawid:ls_drawid, $
		abouttlb:0, savnamenulltlb:0, savnamedoubletlb:0, errtlb:0, tab_tlb:tab_tlb $
	}
;--------------------------------------------------------------------------------- WINDOW SIZES
	winsizes = { $
		aboutwinx:startup_nx, aboutwiny:startup_ny, aboutwinxoffset:startup_xpos, aboutwinyoffset:startup_ypos, $
		lswinx:lswinx, lswiny:lswiny, windowx:windowx, windowy:windowy $
	}
;--------------------------------------------------------------------------------- DEFINE INFO POINTER
	info = { $
		ctrlsgen:PTR_NEW(ctrlsgen,/NO_COPY),$
		ctrlsmeas:PTR_NEW(ctrlsmeas,/NO_COPY),$
		ctrlssetup:PTR_NEW(ctrlssetup,/NO_COPY),$
		ctrlsview:PTR_NEW(ctrlsview,/NO_COPY),$
		curs:PTR_NEW(curs,/NO_COPY),$
		data:PTR_NEW(data,/NO_COPY), $
		dataparams:PTR_NEW(dataparams,/NO_COPY), $
		dispswitch:PTR_NEW(dispswitch,/NO_COPY),$
		feedbparams:PTR_NEW(feedbparams,/NO_COPY),$
		measparams:PTR_NEW(measparams,/NO_COPY), $
		overlays:PTR_NEW(overlays,/NO_COPY), $
		pbparams:PTR_NEW(pbparams,/NO_COPY),$
		plotaxes:PTR_NEW(plotaxes,/NO_COPY),$
		plotparams:PTR_NEW(plotparams,/NO_COPY),$
		plotpos:PTR_NEW(plotpos,/NO_COPY),$
		scaling:PTR_NEW(scaling,/NO_COPY),$
		versioninfo:PTR_NEW(versioninfo,/NO_COPY), $
		winids:PTR_NEW(winids,/NO_COPY), $
		winsizes:PTR_NEW(winsizes,/NO_COPY) $
	}
		
	info = PTR_NEW(info, /NO_COPY)
	WIDGET_CONTROL, control_panel, SET_UVALUE = info

	pseudoevent = { WIDGET_BUTTON, id:control_panel, top:control_panel, handler:0L, select:1 }

	IF (*(*info).dispswitch).slab_set THEN TANAT_UPDATE_LP, pseudoevent ;ELSE BEGIN
;		WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_slider, SENSITIVE = 0
;		WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_speed_slider, SENSITIVE = 0
;		WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_slider, SENSITIVE = 0
;		WIDGET_CONTROL, (*(*info).ctrlsmeas).lp_blink_button, SENSITIVE = 0
;	ENDELSE
	TANAT_UPDATE_LP, pseudoevent
	TANAT_SCALING_RANGE, pseudoevent
	TANAT_DRAW, pseudoevent
	IF (slab_set NE 1) THEN BEGIN
		WSET, (*(*info).winids).ls_drawid
		TV,CONGRID(REPLICATE(200,10,10),(*(*info).winsizes).lswinx,(*(*info).winsizes).lswiny)
		XYOUTS,(*(*info).winsizes).lswinx/2.,(*(*info).winsizes).lswiny/2.,'Could not display spectral information as!Conly one spectral position is available.', COLOR = 0, ALIGNMENT = 0.5, /DEVICE
	ENDIF
	TANAT_SET_SAVEFILENAME, pseudoevent

;--------------------------------------------------------------------------------- START MANAGING
	XMANAGER, 'TANAT', control_panel, /NO_BLOCK
END
