;---------------------------------------------------------------------------------------------------------
;+
; NAME:
;      	CSAV2CLSAV
;
; PURPOSE:
;     	Conversion of CRISPEX output *.CSAV files into CRISPEX input/output *.CLSAV files
;
; CATEGORY:
;      	CRISPEX auxiliary routine
;
; CALLING SEQUENCE:
;	CSAV2CLSAV, INPUTFILENAME, NX, NY
;
; INPUTS:
;	INPUTFILENAME	= One (or an array of) *.CSAV file(s).
;	NX		= The number of pixels in the x-direction.
;	NY		= The number of pixels in the y-direction.
;
; KEYWORDS:
;	None.
;
; OUTPUTS:
;	For each input *.CSAV file one *.CLSAV file.
;
; COMMON BLOCKS:
;     	None.
;
; SIDE EFFECTS:
;     	None.
;
; RESTRICTIONS:
;     	None.
;
; PROCEDURE:
;	Calling the routine with:
;		
;	CSAV2CLSAV, 'halpha_11Jun2008_2011Apr12_140515.csav', 925, 948
;
;	Will return one file: 'halpha_11Jun2008_2011Apr12_140515.clsav', which can be read in
;	by CRISPEX again.
;	
;	One may also supply an array of input files:
;		
;	CSAV2CLSAV, ['halpha_11Jun2008_2011Apr12_140515.csav','halpha_11Jun2008_2011Apr12_140759.csav'], 925, 948
;
; MODIFICATION HISTORY:
;	15 Dec 2009 GV:	Initial routine
;	12 Apr 2011 GV: Extended program info, small revisions
;
; AUTHOR:
;	Gregal Vissers (g.j.m.vissers@astro.uio.no)
;	@ Institute of Theoretical Astrophysics, University of Oslo
;-
;---------------------------------------------------------------------------------------------------------

PRO CSAV2CLSAV, inputfilename, nx, ny

	IF (N_PARAMS() NE 3) THEN BEGIN
		PRINT,'CSAV2CLSAV, INPUTFILENAME, NX, NY'
		RETURN
	ENDIF

	FOR i=0,N_ELEMENTS(inputfilename)-1 DO BEGIN
		RESTORE, inputfilename[i]
		xp = x_coords
		yp = y_coords
		SPLINE_P,xp,yp,xr,yr,INTERVAL=1
		w_lpts = WHERE((xr GE 0) AND (xr LE nx-1) AND (yr GE 0) AND (yr LE ny-1), nw_lpts)
	
		x_coords = xp
		y_coords = yp
		x_loop_pts = xr
		y_loop_pts = yr
		w_loop_pts = w_lpts
		
		IF (N_ELEMENTS(t_saved) EQ 1) THEN tsav = 1 ELSE tsav = 0
		IF (N_ELEMENTS(crispex_version) EQ 2) THEN cvers = 1 ELSE cvers = 0
	
		dirfilename = FILE_DIRNAME(inputfilename[i],/MARK_DIRECTORY)
		splitfilename = STRMID(inputfilename[i],STRLEN(dirfilename),STRLEN(inputfilename[i]))
		outputfilename = STRMID(splitfilename,0,STRPOS(splitfilename,'.',/REVERSE_SEARCH))+'.clsav'
		IF (tsav AND cvers) THEN $
			SAVE, crispex_version, spect_pos, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, t_saved, FILENAME = dirfilename+outputfilename ELSE $
		IF tsav THEN $
			SAVE, spect_pos, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, t_saved, FILENAME = dirfilename+outputfilename ELSE $
		IF cvers THEN $
			SAVE, crispex_version, spect_pos, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, FILENAME = dirfilename+outputfilename ELSE $
			SAVE, spect_pos, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, FILENAME = dirfilename+outputfilename
		PRINT,'Written: '+dirfilename+outputfilename
	ENDFOR
END
