;---------------------------------------------------------------------------------------------------------
;+
; NAME:
;      	CRISPEX: CRIsp SPectral EXplorer
;
; PURPOSE:
;     	Initially, analysis of temporal spectropolarimetric data obtained with the CRISP instrument
;	at the Swedish Solar Telescope (SST) at La Palma (Spain), however any data formatted in a 
;	particular way (cf. the online reference pages) may be fed to CRISPEX for browsing and/on
;	analysis purposes
;
; CATEGORY:
;      	Image data browsing and analysis
;
; CALLING SEQUENCE:
;    	CRISPEX, IMCUBE, SPCUBE, REFCUBE=refcube, SPECTFILE=spectfile, 
;			LINE_CENTER=line_center, MNSPEC=mnspec, DT=dt, EXTS=exts, SINGLE_CUBE=single_cube,
;			SCALE_STOKES=scale_stokes, SCALE_CUBES=scale_cubes, VALS_IMG=vals_img, VALS_REF=vals_ref, 
;			NO_WARP=no_warp, XTITLE=xtitle, YTITLE=ytitle, MASKCUBE=maskcube, WINDOW_LARGE=window_large, VERBOSE=verbose
;
; INPUTS:
;	IMCUBE		= 3D image datacube (dimensions [nx,ny,nt*nlp*ns]) or, if SPCUBE is not provided, a scan (dimensions [nx,ny,nlp]). 
;			  If SINGLE_CUBE is specified a 3D datacube may be provided even if SPCUBE is not. Required input.
;	SPCUBE		= 3D spectral datacube (dimensions [nlp,nt,nx*ny*ns]). Required input if SINGLE_CUBE is not supplied with the number
;			  of linepositions. Otherwise optional input.
;
; KEYWORDS:
;	REFCUBE		= Reference data of same spatial dimensions as main data. REFCUBE may be supplied with either:
;			  2D data array: dimensions [nx,ny].
;			  3D data array: dimensions [nx,ny,refnt], where refnt must be equal to nt.
;			  Scalar string: pointing to a reference image cube of dimensions [nx,nx,refnt*refnlp] where refnt must be equal to
;				1 or to nt, but refnlp need not be equal to nlp.
;		 	  2-element string array: first element pointing to a reference image cube (dimensions [nx,ny,nt*refnlp]), 
;				the second element pointing to the corresponding reference spectral cube (dimensions 
;				[refnlp,refnt,nx*ny]). Again, refnt must be equal to 1 or nt, while any (positive integer) value of
;				refnlp is allowed.
;	SPECTFILE	= File containing the normalised spectrum as function of the linepositions or wavelength. Contains the required variables:
;			  NORM_SPECT (normalised spectrum), NORM_FACTOR (normalisation factor used to produce NORM_SPECT) and SPECT_POS
;			  (spectral wavelength positions in any desired unit). Optional variables: XTITLE and YTITLE (x- and y-title label 
;			  for plots, respectively). If not set, the mean spectrum is determined from the scan(s) determined by the MNSPEC 
;			  keyword (i.e. from all x,y pixels at all linepositions at t given by the setting of MNSPEC).
;			  Scalar string: spectral file corresponding to main data.
;			  2-element string array: first element is the spectral file corresponding to the main data,
;				second element corresponds to the reference data.
;			  Set first element to '' if you only want to set the reference spectral file, e.g. SPECTFILE=['','reference.spectfile'].
;	LINE_CENTER	= Specifies line centre or wavelength information. Not set: linecentre is determined from the data or from SPECTFILE.
;			  Integer scalar: linecentre is set to position specified by LINE_CENTER.
;			  1D 2-element array (format: [WAVELENGTH, DELTA_LAMBDA]): linecentre is determined from 
;				the data and set to WAVELENGTH. The distance in wavelength between the linepositions is specified by 
;				DELTA_LAMBDA.
;			  1D 3-element array (format: [Integer scalar, WAVELENGTH, DELTA_LAMBDA]): combination of
;				the two above, as linecentre is specified by the integer scalar and set to 
;				WAVELENGTH with tickmarks at DELTA_LAMBDA distances from eachother.
;			  2D 1-element integer array (format: [[Integer scalar], [Integer scalar]]): first element sets the linecentre
;				position for the main data, the second element that of the reference data.
;			  2D 2-element array (format: [[MAIN_WAVELENGTH, MAIN_DELTA_LAMBDA],[REF_WAVELENGTH, REF_DELTA_LAMBDA]]):
;				where the elements from the first subarray set the linecentre wavelength and wavelength spacing
;				for the main data, while those of the second subarray set those values for the reference data.
;			  2D 3-element array (format: [[Integer scalar, MAIN_WAVELENGTH, MAIN_DELTA_LAMBDA],[Integer scalar,
;				REF_WAVELENGTH, REF_DELTA_LAMBDA]]): combination of the two above as linecentre for the main data
;				is specified by the first integer scalar and set to MAIN_WAVELENGTH with tickmarks at MAIN_DELTA_LAMBDA 
;				distances from eachother, and correspondingly for the reference data using the values from the second
;				subarray.				
;	MNSPEC		= Not set: mean spectrum is determined from the t=0 scan.
;			  Integer scalar: mean spectrum is determined from the t=MNSPEC scan.
;			  2-element integer array: mean spectrum is determined from the t=MNSPEC[0] through t=MNSPEC[1] scans.
;	DT		= Specifies the elapsed time in seconds per time step. Defaults to not defined, 
;			  showing frame number instead of time on the vertical axes.
;	EXTS		= If set, the time slices/slabs displayed in the program will be exact timeslices, 
;			  obtained through linear interpolation, rather than approximated timeslices, obtained
;			  through nearest-neighbour interpolation (which is the default setting). Note that 
;			  setting this keyword may slow down the browsing of the spectral range (i.e. movement
;			  of the spectral slider) considerably.
;	SINGLE_CUBE	= Single integer value specifying the number of spectral positions of the datacube.
;			  Only to be used when IMCUBE is provided with a 3D spectrotemporal datacube and SPCUBE is not specified.
;	SCALE_STOKES	= If set, the detailed and average spectra of Stokes Q, U and/or V will be scaled to the
;			  maximum of Stokes I (i.e. I/I, Q/I, U/I and/or V/I). If not set, each Stokes component will be scaled 
;			  to its respective maximum.
;	SCALE_CUBES	= Specifies the value that the data should be multiplied with. Default set to 1.
;			  Integer scalar: value used for the main data.
;			  2-element integer array: first element is used for the main data, the second for the reference data.
;	VALS_IMG	= If set, the value of the pixel of the image image under the cursor will be returned in the parameters overview window.
;	VALS_REF	= If set, the value of the pixel of the reference image under the cursor will be returned in the parameters overview window.
;	NO_WARP		= Prevents the warping of the temporal spectrum when the wavelength spacing is non-equidistant. Applies equally
;			  to both sets of data if reference data is supplied. Defaults to not set.
;	XTITLE		= Sets the x-title of the temporal spectrum and the detailed spectrum.
;			  Scalar string: x-title labels for the main temporal spectrum and detailed spectrum.
;			  2-element string array: first element sets the x-title labels for the main temporal spectrum and detailed
;				spectrum, while the second element sets the labels for the reference data.
;			  Set first element to '' if you only want to set the reference x-title, e.g. XTITLE=['','Height [km]'].
;	YTITLE		= Sets the y-title of the detailed spectrum.
;			  Scalar string: y-title label for the detailed spectrum.
;			  2-element string array: first element sets the y-title label for the main detailed spectrum, while the 
;				second element sets the label for the reference data.
;			  Set first element to '' if you only want to set the reference y-title, e.g. YTITLE=['','Velocity [km/s]'].
;	MASKCUBE	= Mask data of same spatial dimensions as main data. Must be supplied with a scalar string pointing to a 
;			  mask image cube of dimensions [nx,nx,masknt] where masknt must be equal to nt or 1.
;	WINDOW_LARGE	= Override the "1:1 window scaling whenever possible" setting. Useful for data with small nx and/or ny, where 1:1 image 
;			  display is possible, but would yield small image display windows. Default set to 0.
;	VERBOSE		= Verbosity setting for program setup and running. Mainly for maintenance purposes. Verbosity levels can be set
;			  by supplying the keyword with the following values (add bitwise):
;				0:  No verbosity.
;				1:  Basic setup verbosity.
;				2:  Extended setup verbosity.
;				4:  Basic runtime verbosity.
;				8:  Extended runtime verbosity.
;				16: Playback statistics verbosity.
;			  In practice the values 3 and 7 are not useful and 1 has become obsolete for all practical purposes. All values 
;			  larger than 26 are reduced to 26, all values smaller than 0 are set to 0.
;
; OUTPUTS:
;	Window outputs:
;    		Depending on the call. In all cases at least three windows will appear, one of which being the
;		control panel, one the main image window and one the detailed spectrum window. If called with:
;			- one datacube, then an additional window will open, showing the spectrum along a slit;
;			- two datacubes, then an additional window will open, showing the temporal spectrum;
;			- three datacubes, then two additional windows will open, one showing the temporal
;			  spectrum and one showing the reference image from the third cube;
;			- four datacubes, then four additional windows will open, one showing the main temporal spectrum,
;			  one showing the reference temporal spectrum, one showing the reference detailed spectrum and
;			  one showing the reference image.
;		Additional data display windows may be accessed through the tabs in the control panel in all 
;		cases, although not all data display windows may be available, depending on the number of data 
;		cubes with which the program is called.
;	Saveable outputs:
;		- loop path points (*.CLSAV file);
;		- loop path or detection timeslice or timeslab (*.CSAV file);
;		- intensity versus time data for specific linepositions (*.CINT file);
;		- (selected) timeseries as MPEG movie;
;		- selected frame as JPEG snapshot;
;		- (selected) timeseries as JPEG files;
;		- current session (*.CSES file).
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
;     	In default setting, four windows are opened, one control panel and three subsidiary windows containing
;	the main image, the temporal spectrum and the detailed spectrum, respectively. Additional data browsing and
;	analysis options/windows may be obtained through control panel options.
;
;	Some example calling sequences are discussed below. The basic calling sequence is:
;
;		CRISPEX, 'main.imcube', 'main.spcube'
;
;	Spectral information may be supplied either through the use of a spectral save file or the use of the 
;	LINE_CENTER keyword, e.g.:
;
;		CRISPEX, 'main.imcube', 'main.spcube', SPECTFILE='main.spectfile'
;	or
;		CRISPEX, 'main.imcube', 'main.spcube', LINE_CENTER=[6562,0.1]
;
;	Reference data may be viewed by supplying such data to the REFCUBE keyword, either in the simple reference mode:
;
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE='ref.imcube'
;
;	or in the dual cube mode:
;
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube']
;
;	In addition, when running in dual cube mode, one may provide arrays to certain keywords in order to set
;	options for both the main and the reference data, e.g.:
;	
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			SPECTFILE=['main.spectfile','ref.specftile']
;	or
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			LINE_CENTER=[[6562,0.1],[8542,0.1]]
;	or
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			XTITLE=['Height [km]','Height [km]']
;	etc.
;
;	One may also set some of the options for the reference cube only (while retaining the default options for the main
;	data) by setting the element corresponding to the main data to an empty scalar string, e.g.:
;
;		CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			YTITLE=['','Velocity [km/s]']
;
;	Shortcut controls (requires focus on the control panel to work):
;		- File menu:
;			- 'About' (Ctrl+A);
;			- 'Preferences (Ctrl+P);
;			- 'Help' (Ctrl+H, redirecting to online reference pages);
;			- 'Quit' (Ctrl+Q).
;		- Control menu:
;			- 'Playback options':
;				- 'one frame backward' (Shift+B);
;				- 'play backwards' (Shift+Backspace);
;				- 'pause' (Shift+Spacebar);
;				- 'play forwards' (Shift+Tab);
;				- 'one frame forwards' (Shift+F).
;			- 'Spectral options':
;				- 'increase spectral position' (Shift+S); 
;				- 'decrease spectral position' (Shift+A).
;			- 'Zoom options':
;				- '1x' (Ctrl+Shift+1); 
;				- '2x' (Ctrl+Shift+2);
;				- '3x' (Ctrl+Shift+3);
;				- '4x' (Ctrl+Shift+4);
;				- '6x' (Ctrl+Shift+6);
;				- '8x' (Ctrl+Shift+8).
;			- 'Runtime options:
;				- 'Interrupt' (Ctrl+Shift+C).
;
;	Note on the time slice along a loop: for better user performance a quicker, but less exact method is 
;	employed in the in-program display of the time slice. Once saving the current timeslice or all the 
;	loop data, the computationally more expensive interpolation of the data is used to extract the 
;	timeslice (or slab). The differences in the time slice obtained with the quicker method and with the 
;	interpolation method are small (though non-negligible when considering velocities) and the quick 
;	method is good enough for the explorational purposes this code was written for. 
;	Additionally, when saving a time slice or slab only that part of the cube between the lower t-value 
;	and the upper t-value will be saved, i.e. if you wish to save it for the full range in time, be sure 
;	to reset the temporal boundaries.
;
;	Further information on the calling sequence and specific options can be found in the online reference pages,
;	which can be reached through the in-program help function or by going directly to:
;
;		http://folk.uio.no/~gregal/crispex
;
; MODIFICATION HISTORY:
;	11 Feb 2009 GV: start buildup of program, based on Øystein Langangen's August 2008 
;			version of CRISP_SPECTRAL_EXPLORER.PRO
;			further (structural) inspiration from XIMOVIE.PRO and XSLICE.PRO
;	21 Feb 2009 GV: incorporation of LP_HEADER.PRO functionality
;	24 Feb 2009 GV: release of beta version 							(v0.9)
;	02 Mar 2009 GV: implementation of x- and y-slice display 					(v0.9.1)
;	04 Mar 2009 GV: implementation of display of spectral slice along a slit, including controls	(v0.9.3)
;			to set slit angle and length
;	09 Mar 2009 GV: implementation of movement along the slit direction				(v0.9.4)
;	10 Mar 2009 GV: implementation of option to study spectral scan and show a reference image	(v0.9.5)
;			or cube
;	15 Mar 2009 GV: implementation of zoom option							(v0.9.6)
;	23 Mar 2009 GV: implementation of time slice along a segmented line and	option to save the	(v0.9.7)
;			resulting data
;	14 Apr 2009 GV: corrected extraction and saving of time slice, removed x- and y-slice display	(v0.9.8)
;			due to redundancy from spectral slice along a slit
;	15 May 2009 GV: implementation of extended save and retrieve options, extended image scaling	(v0.9.9)
;			options, save and restore session options, adjustable time and spectral range, 
;			different loop linestyles, selection menu for loop overlays, parameter overview 
;			window and display of saved timeslices
;	24 Aug 2009 GV: implementation of extra zoomfactors, save as MPEG and JPEG options, save from	(v1.0)
;			detection file, spectral and temporal range choice in saving timeslabs, exact 
;			timeslice display in-program, resizable display windows, read-in of full 
;			reference cubes and an option to calculate mean spectrum over a range in 
;			timesteps
;	05 Nov 2009 GV: implementation of loop path feedback, shortcut controls through keyboard,	(v1.1)
;			single full cube call and also fixed a number of bugs
;	21 Mar 2010 GV: enabled visualisation of Stokes cube data, display of reference and image cube	(v1.5)
;			data values, drawing of loop paths for 3D temporal image cube, retreival and
;			saving timeslices from reference cube, extended scaling options for reference
;			cube image, implemented spatial measurement tool and help function, moved user
;			feedback to pop-up windows, disposed of obsolete keywords and fixed a number
;			of bugs
;	01 Jul 2011 GV: enabled dual cube mode, in-program Doppler images, extended (Stokes) spectral	(v1.6)
;			options, extended plot options, blinking while playing, setting of preferences, 
;			extraction of intensity-time plots, and made aesthetic improvements (bitmap
;			play buttons, better cursor visibility and startup screen).
;	22 Aug 2011 GV: extended save as options, implemented save as PNG, saving of color MPEG, an	(v1.6.1)
;			option to open a restored loop in TANAT and fixed a number of bugs
;	23 Mar 2012 GV: extended reference and Stokes data input options, implemented mask overlays,	(v1.6.2)
;			enabled display of multiple space-time diagrams and fixed a number of bugs
;	04 Dec 2012 GV: extended in-program analysis options through display of reference space-time 	(v1.6.3)
;			diagram, extended image and space-time diagram scaling options and fixed a
;			number of bugs
;
; ACKNOWLEDGEMENTS:
;	This code would not be present in its current state and with the current functionalities without
;	the relentless practical testing and the valuable input and ideas of Luc Rouppe van der Voort,
;	Sven Wedemeyer-Böhm, Mats Carlsson, Patrick Antolin, Jorrit Leenaarts, Bart de Pontieu, 
;	Eamon Scullion and Jaime de la Cruz Rodriguez. 
;
; AUTHOR:
;	Gregal Vissers (g.j.m.vissers@astro.uio.no)
;	@ Institute of Theoretical Astrophysics, University of Oslo
; $Id$
;-
;---------------------------------------------------------------------------------------------------------

;================================================================================= ABOUT WINDOW PROCEDURES
PRO CRISPEX_ABOUT_WINDOW, event 							
; Creates an about-window displaying code name, version and revision number
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ABOUT_WINDOW'
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': ABOUT'
	CRISPEX_WINDOW, (*(*info).winsizes).aboutwinx, (*(*info).winsizes).aboutwiny, (*(*info).winids).root, title, abouttlb, aboutwid, (*(*info).winsizes).aboutxoffset, (*(*info).winsizes).aboutyoffset, $
		DRAWID = aboutdrawid, DRAWBASE = aboutdrawbase
	CRISPEX_UPDATE_STARTUP_FEEDBACK, (*(*info).feedbparams).startup_im, (*(*info).feedbparams).xout, (*(*info).feedbparams).yout, $
		['Running CRISPEX version '+(*(*info).versioninfo).version_number+' ('+(*(*info).versioninfo).revision_number+')','',$
		'Developed by: Gregal Vissers', $
		'               Institute of Theoretical Astrophysics,',$
		'               University of Oslo',$
		'               2009-2013']
	WIDGET_CONTROL, aboutdrawid, EVENT_PRO = 'CRISPEX_ABOUT_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS, /DRAW_BUTTON_EVENTS
	WIDGET_CONTROL, abouttlb, SET_UVALUE = info
	XMANAGER, 'CRISPEX', abouttlb,/NO_BLOCK
	(*(*info).winids).abouttlb = abouttlb
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).abouttlb], labels=['abouttlb']
END

PRO CRISPEX_ABOUT_CURSOR, event
; Handles cursor actions on the about window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ABOUT_CURSOR'
	IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_DRAW' THEN BEGIN
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [event.TYPE,event.PRESS], labels=['WIDGET_DRAW: event.TYPE','WIDGET_DRAW: event.PRESS']
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

;================================================================================= BINARY CONVERSION FUNCTIONS
FUNCTION CRISPEX_DEC2BIN, decimal_number
; Handles the conversion of decimal to binary number
	IF (decimal_number NE 0) THEN BEGIN
		coarse_p = ALOG10(DOUBLE(decimal_number))/ALOG10(2.D)
		p = FLOOR(coarse_p)
		i=0
		binary_array = INTARR(FLOOR(coarse_p)+1)
		b = REVERSE((2^FINDGEN(FLOOR(coarse_p)+1)))
		WHILE (p GE 0) DO BEGIN
			IF (2^p LE decimal_number) THEN BEGIN
				binary_array[i] = 1
				decimal_number -= 2^p
			ENDIF ELSE binary_array[i] = 0
			p -= 1
			i += 1
		ENDWHILE
		binary_array = REVERSE(binary_array)
		IF (N_ELEMENTS(binary_array) LT 5) THEN binary_array = [binary_array, REPLICATE(0,5-N_ELEMENTS(binary_array))]
	ENDIF ELSE binary_array = [0,0,0,0,0]
	RETURN, binary_array
END

;================================================================================= CLEAR ESTIMATE PROCEDURES
PRO CRISPEX_CLEAR_CURRENT_ESTIMATE, event								
; Clears current saving time estimate
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CLEAR_CURRENT_ESTIMATE'
	IF (*(*info).feedbparams).estimate_run THEN BEGIN
		(*(*info).feedbparams).estimate_lx = 0
		(*(*info).feedbparams).estimate_time = 0.
		(*(*info).feedbparams).estimate_run = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_estimate, SENSITIVE = 0
	ENDIF
END

PRO CRISPEX_CLEAR_CURRENT_CPFT, event								
; Clears current saving time estimate
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CLEAR_CURRENT_ESTIMATE'
	pftfiles = FILE_SEARCH((*(*info).paths).dir_cpft+'crispex.'+(*(*info).paths).hostname+'.cpft', COUNT = pftfilecount)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, ['crispex.'+(*(*info).paths).hostname+'.cpft',pftfilecount], labels=['File to be deleted','Filecount']
	IF pftfilecount THEN BEGIN
		SPAWN,'rm '+(*(*info).paths).dir_cpft+'crispex.'+(*(*info).paths).hostname+'.cpft'
		(*(*info).feedbparams).estimate_lx = 0
		(*(*info).feedbparams).estimate_time = 0.
		(*(*info).feedbparams).estimate_run = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_estimate, SENSITIVE = 0
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','Could not delete crispex.'+((*(*info).paths).hostname)[0]+'.cpft','from '+(*(*info).paths).dir_cpft+'.','File does not exist.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END

PRO CRISPEX_CLEAR_CURRENT_INST, event								
; Clears current saving time estimate
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CLEAR_CURRENT_ESTIMATE'
	instfiles = FILE_SEARCH((*(*info).paths).dir_inst+'crispex.'+(*(*info).paths).hostname+'.inst', COUNT = instfilecount)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, ['crispex.'+(*(*info).paths).hostname+'.inst',instfilecount], labels=['File to be deleted','Filecount']
	IF instfilecount THEN BEGIN
		SPAWN,'rm '+(*(*info).paths).dir_inst+'crispex.'+(*(*info).paths).hostname+'.inst'
		(*(*info).feedbparams).estimate_lx = 0
		(*(*info).feedbparams).estimate_time = 0.
		(*(*info).feedbparams).estimate_run = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_inst, SENSITIVE = 0
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','Could not delete crispex.'+((*(*info).paths).hostname)[0]+'.inst','from '+(*(*info).paths).dir_inst+'.','File does not exist.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END
;================================================================================= PROGRAM EXIT PROCEDURES
PRO CRISPEX_CLOSE, event								
; Called upon closing program, checks for existence of performance test file; if not present it is written
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CLOSE'
	IF ((*(*info).paths).dir_cpft_write EQ 1) THEN BEGIN
		pftfiles = FILE_SEARCH((*(*info).paths).dir_cpft+'crispex.'+(*(*info).paths).hostname+'.cpft', COUNT = pftfilecount)
		IF (pftfilecount EQ 0) AND ((*(*info).feedbparams).estimate_run EQ 1) THEN BEGIN
			estimate_lx = (*(*info).feedbparams).estimate_lx
			estimate_time = (*(*info).feedbparams).estimate_time
			estimate_run = (*(*info).feedbparams).estimate_run
			SAVE, estimate_lx, estimate_time, estimate_run, FILENAME = (*(*info).paths).dir_cpft+'crispex.'+(*(*info).paths).hostname+'.cpft'
			PRINT,'Written: '+(*(*info).paths).dir_cpft+'crispex.'+(*(*info).paths).hostname+'.cpft'
		ENDIF
	ENDIF ELSE BEGIN
		PRINT, 'ERROR: Could not write performance file crispex.'+(*(*info).paths).hostname+'.cpft '
		PRINT, '       to '+(*(*info).paths).dir_cpft+'. Permission denied.'
	ENDELSE
	FREE_LUN, (*(*info).data).lun
	IF (*(*info).dataswitch).spfile THEN FREE_LUN, (*(*info).data).lur
	IF ((*(*info).winswitch).showref AND ((*(*info).data).luf GT 0)) THEN FREE_LUN, (*(*info).data).luf
	IF ((*(*info).dataswitch).refspfile AND ((*(*info).data).lufs GT 0)) THEN FREE_LUN, (*(*info).data).lufs
	IF ((*(*info).dataswitch).maskfile AND ((*(*info).data).lum GT 0)) THEN FREE_LUN, (*(*info).data).lum
	WIDGET_CONTROL, (*(*info).winids).root, /DESTROY
	PTR_FREE, info
END

PRO CRISPEX_CLOSE_CLEANUP, base								
; Clean-up upon closing program
	WIDGET_CONTROL, base, GET_UVALUE = info
	FREE_LUN, (*(*info).data).lun
	IF (*(*info).dataswitch).spfile THEN FREE_LUN, (*(*info).data).lur
	IF ((*(*info).winswitch).showref AND ((*(*info).data).luf GT 0)) THEN FREE_LUN, (*(*info).data).luf
	IF ((*(*info).dataswitch).refspfile AND ((*(*info).data).lufs GT 0)) THEN FREE_LUN, (*(*info).data).lufs
	IF ((*(*info).dataswitch).maskfile AND ((*(*info).data).lum GT 0)) THEN FREE_LUN, (*(*info).data).lum
	CRISPEX_CLOSE_CLEAN_INSTANCE_FILE, (*(*info).paths).dir_inst_write, (*(*info).paths).dir_inst, (*(*info).paths).hostname, ((*(*info).sesparams).curr_instance_id)[0]
	PTR_FREE, info
END

PRO CRISPEX_CLOSE_CLEAN_INSTANCE_FILE, dir_inst_write, dir_inst, hostname, curr_instance_id
; Called upon closing program, checks for existence of performance test file; if not present it is written
	IF (dir_inst_write EQ 1) THEN BEGIN
		instfile = FILE_SEARCH(dir_inst+'crispex.'+hostname+'.inst', COUNT = instfilecount)
		IF instfilecount THEN BEGIN
			nlines = (FILE_LINES(instfile))[0]
			datarr = STRARR(1,nlines)
			OPENR,unit1,instfile,/GET_LUN
			READF,unit1,datarr
			FREE_LUN,unit1
			routine_name = STRARR(nlines)
			instance_id = LONARR(nlines)
			FOR i=1,nlines-1 DO BEGIN
				splitline = STRSPLIT(datarr[i],'	',/EXTRACT)
				routine_name[i] = splitline[0]
				instance_id[i] = splitline[3]
			ENDFOR
			where_crispex = WHERE(routine_name EQ 'CRISPEX')
			sel_instance_id = instance_id[where_crispex]
			clean_line = WHERE(sel_instance_id EQ curr_instance_id)+1
			first_part = datarr[0:(clean_line-1)]
			IF (clean_line NE (nlines-1)) THEN BEGIN
				last_part = datarr[(clean_line+1):(nlines-1)] 
				rewritten_arr = [first_part,last_part]
			ENDIF ELSE rewritten_arr = first_part
			OPENW, unit, instfile[0], WIDTH = 360, /GET_LUN
			FOR i=0,nlines-2 DO PRINTF, unit,rewritten_arr[i]
			FREE_LUN, unit
		ENDIF
	ENDIF
END

PRO CRISPEX_CLOSE_EVENT_WINDOW, event
; Called upon closing window when no extra processes need to be run
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CLOSE_EVENT_WINDOW'
	IF (event.TOP EQ (*(*info).winids).restsestlb) THEN (*(*info).winids).restsestlb = 0
	IF (event.TOP EQ (*(*info).winids).savewintlb) THEN BEGIN
		(*(*info).winids).savewintlb = 0
		IF ((*(*info).winids).saveoptwintlb GT 0) THEN BEGIN
			WIDGET_CONTROL, (*(*info).winids).saveoptwintlb,/DESTROY
			(*(*info).winids).saveoptwintlb = 0
		ENDIF
	ENDIF
	IF (event.TOP EQ (*(*info).winids).saveoptwintlb) THEN (*(*info).winids).saveoptwintlb = 0
	IF (event.TOP EQ (*(*info).winids).abouttlb) THEN (*(*info).winids).abouttlb = 0
	IF (event.TOP EQ (*(*info).winids).errtlb) THEN (*(*info).winids).errtlb = 0
	IF (event.TOP EQ (*(*info).winids).warntlb) THEN (*(*info).winids).warntlb = 0
	WIDGET_CONTROL, event.TOP, /DESTROY
END

;================================================================================= CURSOR PROCEDURES
PRO CRISPEX_CURSOR, event								
; Cursor handling procedure, tracks and handles events from the cursor on the main image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CURSOR', /IGNORE_LAST
	IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_TRACKING' THEN BEGIN
		IF event.ENTER THEN BEGIN
			WIDGET_CONTROL, event.HANDLER, get_value = wid
			WSET, wid
			ci = UINTARR(16) & cim = ci & cim[7] = 1
			DEVICE, CURSOR_IMAGE = ci, CURSOR_MASK = cim, CURSOR_XY = [8,8]
		ENDIF ELSE BEGIN
			IF (((*(*info).loopparams).np GE 1) AND (*(*info).overlayswitch).looppath_feedback AND ((*(*info).curs).lockset GT 0)) THEN BEGIN
				*(*(*info).loopparams).xp = (*(*(*info).loopparams).xp)[*,0:(*(*info).loopparams).np-1]
				*(*(*info).loopparams).yp = (*(*(*info).loopparams).yp)[*,0:(*(*info).loopparams).np-1]
				*(*(*info).overlayparams).sxp = (*(*(*info).overlayparams).sxp)[*,0:(*(*info).loopparams).np-1]
				*(*(*info).overlayparams).syp = (*(*(*info).overlayparams).syp)[*,0:(*(*info).loopparams).np-1]
				IF ((*(*info).loopparams).np GE 2) THEN CRISPEX_LOOP_GET_PATH, event ELSE BEGIN
					*(*(*info).loopparams).xr = *(*(*info).loopparams).xp
					*(*(*info).loopparams).yr = *(*(*info).loopparams).yp
				ENDELSE
				(*(*info).dataparams).x = (*(*(*info).loopparams).xp)[*,(*(*info).loopparams).np-1]
				(*(*info).dataparams).y = (*(*(*info).loopparams).yp)[*,(*(*info).loopparams).np-1]
				CRISPEX_COORDSLIDERS_SET, 0, 0, event
				IF (*(*info).winswitch).showphis THEN BEGIN
					(*(*info).curs).sx = (*(*info).curs).sxlock
					(*(*info).curs).sy = (*(*info).curs).sylock
					CRISPEX_PHISLIT_DIRECTION, event
					CRISPEX_UPDATE_PHIS, event
				ENDIF ELSE CRISPEX_DRAW, event
			ENDIF		
			IF (!D.WINDOW NE -1) THEN DEVICE, /CURSOR_CROSSHAIR 
		ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [event.ENTER], labels=['WIDGET_TRACKING: event.Enter']
	ENDIF ELSE IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_DRAW' THEN BEGIN
		CASE event.TYPE OF
		0:	CASE event.PRESS OF
			1:	BEGIN	; left mouse button press
					(*(*info).curs).lockset = 1
					WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, SET_BUTTON = (*(*info).curs).lockset-1
					WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = (*(*info).curs).lockset
					(*(*info).curs).sxlock = event.X
					(*(*info).curs).sylock = event.Y
					(*(*info).curs).sx = (*(*info).curs).sxlock
					(*(*info).curs).sy = (*(*info).curs).sylock
					IF ((*(*info).zooming).factor EQ 1) THEN BEGIN
						(*(*info).curs).xlock = FLOAT((*(*info).curs).sxlock * (*(*info).dataparams).nx) / (*(*info).winsizes).xywinx
						(*(*info).curs).ylock = FLOAT((*(*info).curs).sylock * (*(*info).dataparams).ny) / (*(*info).winsizes).xywiny
					ENDIF ELSE BEGIN
						(*(*info).curs).xlock = FLOAT((*(*info).curs).sxlock * ((*(*info).dataparams).d_nx+1)) / (*(*info).winsizes).xywinx + (*(*info).zooming).xpos
						(*(*info).curs).ylock = FLOAT((*(*info).curs).sylock * ((*(*info).dataparams).d_ny+1)) / (*(*info).winsizes).xywiny + (*(*info).zooming).ypos
					ENDELSE
					(*(*info).dataparams).x = (*(*info).curs).xlock
					(*(*info).dataparams).y = (*(*info).curs).ylock
					IF (*(*info).overlayswitch).loopslit THEN BEGIN
						(*(*info).loopparams).np += 1
						IF ((*(*info).loopparams).np EQ 2) THEN WIDGET_CONTROL, (*(*info).ctrlscp).loop_slit_but, SET_VALUE = 'Erase loop path'
						IF ((*(*info).loopparams).np GE 2) THEN BEGIN
							*(*(*info).loopparams).xp = [[*(*(*info).loopparams).xp],[(*(*info).curs).xlock]]
							*(*(*info).loopparams).yp = [[*(*(*info).loopparams).yp],[(*(*info).curs).ylock]]
							*(*(*info).overlayparams).sxp = [[*(*(*info).overlayparams).sxp],[(*(*info).curs).sxlock]]
							*(*(*info).overlayparams).syp = [[*(*(*info).overlayparams).syp],[(*(*info).curs).sylock]]
							IF ((*(*info).winids).looptlb EQ 0) THEN WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 1
							IF ((*(*info).loopparams).np EQ 3) THEN WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = 1
							CRISPEX_LOOP_GET, event
							CRISPEX_UPDATE_LP, event
							IF ((*(*info).zooming).factor NE 1) THEN BEGIN
								*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
								*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
							ENDIF
						ENDIF ELSE BEGIN
							(*(*(*info).loopparams).xp)[0] = (*(*info).curs).xlock
							(*(*(*info).loopparams).yp)[0] = (*(*info).curs).ylock
							(*(*(*info).overlayparams).sxp)[0] = (*(*info).curs).sxlock
							(*(*(*info).overlayparams).syp)[0] = (*(*info).curs).sylock
						ENDELSE
					ENDIF
					IF (*(*info).meas).spatial_measurement THEN BEGIN
						(*(*info).meas).np = 1
						(*(*(*info).meas).xp) = (*(*info).curs).xlock
						(*(*(*info).meas).yp) = (*(*info).curs).ylock
						(*(*(*info).meas).sxp) = (*(*info).curs).sxlock
						(*(*(*info).meas).syp) = (*(*info).curs).sylock 
					ENDIF
					CRISPEX_COORDSLIDERS_SET, 0, 0, event
				END
			2:	BEGIN	; middle mouse button press
					IF ((*(*info).meas).spatial_measurement AND ((*(*info).meas).np LT 2)) THEN BEGIN
						(*(*info).curs).sx = event.X
						(*(*info).curs).sy = event.Y
						CRISPEX_CURSOR_GET_XY, event
						(*(*info).meas).np = 2
						*(*(*info).meas).xp = [(*(*(*info).meas).xp)[0],(*(*info).dataparams).x]	
						*(*(*info).meas).yp = [(*(*(*info).meas).yp)[0],(*(*info).dataparams).y]	
						*(*(*info).meas).sxp = [(*(*(*info).meas).sxp)[0],(*(*info).curs).sx]	
						*(*(*info).meas).syp = [(*(*(*info).meas).syp)[0],(*(*info).curs).sy]
						CRISPEX_COORDSLIDERS_SET, 0, 0, event
						CRISPEX_MEASURE_CALC, event
					ENDIF
				END
			4:	BEGIN	; right mouse button press
					(*(*info).curs).lockset = 0
					(*(*info).meas).np = 0
					*(*(*info).meas).xp = 0
					*(*(*info).meas).yp = 0
					*(*(*info).meas).sxp = 0
					*(*(*info).meas).syp = 0
					WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SET_VALUE = '0'
					WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SET_VALUE = '0'
					WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, SET_BUTTON = 1
					WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = 0
					CRISPEX_COORDSLIDERS_SET, 1, 1, event
				END
			ELSE: BREAK
			ENDCASE
		2:	BEGIN	; mouse movement
				IF ((*(*info).curs).lockset EQ 0) THEN BEGIN
					(*(*info).curs).sx = event.X
					(*(*info).curs).sy = event.Y
					CRISPEX_CURSOR_GET_XY, event
					CRISPEX_COORDSLIDERS_SET, 1, 1, event 
				ENDIF ELSE IF ((*(*info).overlayswitch).loopslit AND (*(*info).overlayswitch).looppath_feedback AND ((*(*info).loopparams).np GE 1)) THEN BEGIN
					(*(*info).curs).sx = event.X
					(*(*info).curs).sy = event.Y
					CRISPEX_CURSOR_GET_XY, event
					*(*(*info).loopparams).xp = [[(*(*(*info).loopparams).xp)[*,0:(*(*info).loopparams).np-1]],[(*(*info).dataparams).x]]
					*(*(*info).loopparams).yp = [[(*(*(*info).loopparams).yp)[*,0:(*(*info).loopparams).np-1]],[(*(*info).dataparams).y]]
					*(*(*info).overlayparams).sxp = [[(*(*(*info).overlayparams).sxp)[*,0:(*(*info).loopparams).np-1]],[(*(*info).curs).sx]]
					*(*(*info).overlayparams).syp = [[(*(*(*info).overlayparams).syp)[*,0:(*(*info).loopparams).np-1]],[(*(*info).curs).sy]]
					CRISPEX_LOOP_GET_PATH, event
					IF ((*(*info).zooming).factor NE 1) THEN BEGIN
						*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
						*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
					ENDIF
					CRISPEX_COORDSLIDERS_SET, 0, 0, event
				ENDIF ELSE IF ((*(*info).meas).spatial_measurement AND ((*(*info).meas).np EQ 1)) THEN BEGIN
					(*(*info).curs).sx = event.X
					(*(*info).curs).sy = event.Y
					CRISPEX_CURSOR_GET_XY, event
					*(*(*info).meas).xp = [(*(*(*info).meas).xp)[0],(*(*info).dataparams).x]	
					*(*(*info).meas).yp = [(*(*(*info).meas).yp)[0],(*(*info).dataparams).y]	
					*(*(*info).meas).sxp = [(*(*(*info).meas).sxp)[0],(*(*info).curs).sx]	
					*(*(*info).meas).syp = [(*(*(*info).meas).syp)[0],(*(*info).curs).sy]
					CRISPEX_MEASURE_CALC, event
					CRISPEX_COORDSLIDERS_SET, 0, 0, event
				ENDIF ELSE RETURN
			END
		ELSE: RETURN
		ENDCASE
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [event.TYPE,event.PRESS,(*(*info).dataparams).x,(*(*info).dataparams).y,(*(*info).curs).sx,(*(*info).curs).sy], $
			labels=['WIDGET_DRAW: event.TYPE','WIDGET_DRAW: event.PRESS','x','y','sx','sy']
		IF (*(*info).winswitch).showphis THEN BEGIN
			CRISPEX_PHISLIT_DIRECTION, event
			CRISPEX_UPDATE_PHIS, event
		ENDIF ELSE CRISPEX_DRAW, event
	ENDIF
END

PRO CRISPEX_CURSOR_GET_XY, event
; Converts the window x and y coordinates to data x and y coordinates
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CURSOR_GET_XY'
	IF ((*(*info).zooming).factor EQ 1) THEN BEGIN
		(*(*info).dataparams).x = (*(*info).curs).sx * (*(*info).dataparams).nx / (*(*info).winsizes).xywinx
		(*(*info).dataparams).y = (*(*info).curs).sy * (*(*info).dataparams).ny / (*(*info).winsizes).xywiny
	ENDIF ELSE BEGIN
		(*(*info).dataparams).x = (*(*info).curs).sx * ((*(*info).dataparams).d_nx+1) / (*(*info).winsizes).xywinx + (*(*info).zooming).xpos
		(*(*info).dataparams).y = (*(*info).curs).sy * ((*(*info).dataparams).d_ny+1) / (*(*info).winsizes).xywiny + (*(*info).zooming).ypos
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x,(*(*info).dataparams).y], labels=['x','y']
END

PRO CRISPEX_CURSOR_LOCK, event								
; Called upon locking/unlocking cursor with 'lock cursor' or 'unlock cursor' button, handles cursor (un)locking
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_CURSOR_LOCK'
	(*(*info).curs).lockset = event.SELECT
	IF (*(*info).curs).lockset THEN BEGIN
		(*(*info).curs).xlock = (*(*info).dataparams).x	&	(*(*info).curs).ylock = (*(*info).dataparams).y
		(*(*info).curs).sxlock = (*(*info).curs).sx	&	(*(*info).curs).sylock = (*(*info).curs).sy
		CRISPEX_COORDSLIDERS_SET, 0, 0, event
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
			CRISPEX_VERBOSE_GET, event, [(*(*info).curs).xlock,(*(*info).curs).ylock,(*(*info).curs).sxlock,(*(*info).curs).sylock], labels=['xlock','ylock','sxlock','sylock']
	ENDIF ELSE CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_COORDSLIDERS_SET, xsensitive, ysensitive, event				
; Adjusts sliders according to change in cursor position or locked/unlocked state
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_COORDSLIDERS_SET'
	WIDGET_CONTROL, (*(*info).ctrlscp).x_slider, SET_VALUE = (*(*info).dataparams).x, SENSITIVE = xsensitive
	WIDGET_CONTROL, (*(*info).ctrlscp).y_slider, SET_VALUE = (*(*info).dataparams).y, SENSITIVE = ysensitive
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x,(*(*info).dataparams).y,xsensitive,ysensitive], labels=['x','y','xsensitive','ysensitive']
END

;================================================================================= DISPLAYS PROCEDURES
PRO CRISPEX_DISPLAYS_ALL_TO_FRONT, event
; Brings all opened session windows to front
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info	
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_ALL_TO_FRONT', /IGNORE_LAST
	; Data windows
	WSHOW, (*(*info).winids).imwid
	IF ((*(*info).winids).sptlb NE 0) THEN WSHOW, (*(*info).winids).spwid
	IF ((*(*info).winids).lstlb NE 0) THEN WSHOW, (*(*info).winids).lswid
	IF ((*(*info).winids).phistlb NE 0) THEN WSHOW, (*(*info).winids).phiswid
	IF ((*(*info).winids).reftlb NE 0) THEN WSHOW, (*(*info).winids).refwid
	IF ((*(*info).winids).doptlb NE 0) THEN WSHOW, (*(*info).winids).dopwid
	IF ((*(*info).winids).imreftlb NE 0) THEN WSHOW, (*(*info).winids).imrefwid
	IF (TOTAL(*(*(*info).winids).restlooptlb) NE 0) THEN FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO WSHOW, (*(*(*info).winids).restloopwid)[i]
	IF ((*(*info).winids).retrdettlb NE 0) THEN WSHOW, (*(*info).winids).retrdetwid
	IF ((*(*info).winids).looptlb NE 0) THEN WSHOW, (*(*info).winids).loopwid
	IF ((*(*info).winids).refsptlb NE 0) THEN WSHOW, (*(*info).winids).refspwid
	IF ((*(*info).winids).reflstlb NE 0) THEN WSHOW, (*(*info).winids).reflswid
	IF ((*(*info).winids).inttlb NE 0) THEN WSHOW, (*(*info).winids).intwid
	; Action windows
	IF ((*(*info).winids).savetlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).savetlb, /SHOW
	IF ((*(*info).winids).detsavetlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).detsavetlb, /SHOW
	IF ((*(*info).winids).restoretlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).restoretlb, /SHOW
	IF ((*(*info).winids).restsestlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).restsestlb, /SHOW
	IF ((*(*info).winids).savewintlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).savewintlb, /SHOW
	IF ((*(*info).winids).saveoptwintlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).saveoptwintlb, /SHOW
	IF ((*(*info).winids).intmenutlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).intmenutlb, /SHOW
	IF ((*(*info).winids).preftlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).preftlb, /SHOW
	; Warning and feedback windows
	IF ((*(*info).winids).paramtlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).paramtlb, /SHOW
	IF ((*(*info).winids).estimatetlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).estimatetlb, /SHOW
	IF ((*(*info).winids).feedbacktlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).feedbacktlb, /SHOW
	IF ((*(*info).winids).restsesfeedbtlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
	IF ((*(*info).winids).abouttlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).abouttlb, /SHOW
	IF ((*(*info).winids).errtlb NE 0) THEN WIDGET_CONTROL, (*(*info).winids).errtlb, /SHOW
END

PRO CRISPEX_DISPWIDS, event
; Brings all opened session windows to front
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info	
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPWIDS'
	(*(*info).winswitch).dispwids = ABS((*(*info).winswitch).dispwids-1)
	WIDGET_CONTROL,(*(*info).ctrlscp).dispwid, SET_BUTTON = (*(*info).winswitch).dispwids
	tlbarr = [(*(*info).winids).imtlb,(*(*info).winids).sptlb,(*(*info).winids).lstlb,(*(*info).winids).reftlb,(*(*info).winids).refsptlb,(*(*info).winids).reflstlb,(*(*info).winids).imreftlb,(*(*info).winids).doptlb, $
		(*(*info).winids).phistlb,*(*(*info).winids).restlooptlb,(*(*info).winids).retrdettlb,(*(*info).winids).looptlb,(*(*info).winids).reflooptlb,(*(*info).winids).inttlb]
	widarr = [(*(*info).winids).imwid,(*(*info).winids).spwid,(*(*info).winids).lswid,(*(*info).winids).refwid,(*(*info).winids).refspwid,(*(*info).winids).reflswid,(*(*info).winids).imrefwid,(*(*info).winids).dopwid, $
		(*(*info).winids).phiswid,*(*(*info).winids).restloopwid,(*(*info).winids).retrdetwid,(*(*info).winids).loopwid,(*(*info).winids).refloopwid,(*(*info).winids).intwid]
	title_arr = [(*(*info).winids).imwintitle,(*(*info).winids).spwintitle,(*(*info).winids).lswintitle,(*(*info).winids).refwintitle,(*(*info).winids).refspwintitle,(*(*info).winids).reflswintitle,$
		(*(*info).winids).imrefwintitle,(*(*info).winids).dopwintitle, (*(*info).winids).phiswintitle,*(*(*info).winids).restloopwintitle,(*(*info).winids).retrdetwintitle,(*(*info).winids).loopwintitle,$
		(*(*info).winids).refloopwintitle,(*(*info).winids).intwintitle]
	wherenot0 = WHERE(tlbarr NE 0)
	IF (*(*info).winswitch).dispwids THEN BEGIN
		FOR i=0,N_ELEMENTS(wherenot0)-1 DO WIDGET_CONTROL,tlbarr[wherenot0[i]], BASE_SET_TITLE = STRTRIM(widarr[wherenot0[i]],2)+' - '+title_arr[wherenot0[i]]
	ENDIF ELSE BEGIN
		FOR i=0,N_ELEMENTS(wherenot0)-1 DO WIDGET_CONTROL,tlbarr[wherenot0[i]], BASE_SET_TITLE = title_arr[wherenot0[i]]
	ENDELSE
END

PRO CRISPEX_DISPLAYS_DETSPECT_IM_SELECT, event
; Handles the selection of detspect options for the main image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_DETSPECT_IM_SELECT'
	(*(*info).ctrlsswitch).imrefdetspect = 0
	CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
END

PRO CRISPEX_DISPLAYS_DETSPECT_REF_SELECT, event
; Handles the selection of detspect options for the reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_DETSPECT_REF_SELECT'
	(*(*info).ctrlsswitch).imrefdetspect = 1
	CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
END

PRO CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
; Handles the setting of scaling buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS'
	IF (*(*info).ctrlsswitch).imrefdetspect THEN BEGIN		; If selected options for reference
		WIDGET_CONTROL, (*(*info).ctrlscp).ls_toggle_but, SET_BUTTON = (*(*info).winswitch).showrefls, SET_VALUE = 'Display '+STRLOWCASE(((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).refheightset])
		WIDGET_CONTROL, (*(*info).ctrlscp).subtract_but, SET_BUTTON = (*(*info).plotswitch).ref_subtract
		WIDGET_CONTROL, (*(*info).ctrlscp).scale_detspect_but, SET_BUTTON = (*(*info).dispswitch).ref_detspect_scale, $
			SET_VALUE = 'Scale '+STRLOWCASE(((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).refheightset])+' to maximum of average'
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_low_y_ref,2)
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_upp_y_ref,2)
		WIDGET_CONTROL, (*(*info).ctrlscp).detspect_label, SET_VALUE = ((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).refheightset]+':'
	ENDIF ELSE BEGIN				; If selected options for main
		WIDGET_CONTROL, (*(*info).ctrlscp).ls_toggle_but, SET_BUTTON = (*(*info).winswitch).showls, SET_VALUE = 'Display '+STRLOWCASE(((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).heightset])
		WIDGET_CONTROL, (*(*info).ctrlscp).subtract_but, SET_BUTTON = (*(*info).plotswitch).subtract
		WIDGET_CONTROL, (*(*info).ctrlscp).scale_detspect_but, SET_BUTTON = (*(*info).dispswitch).detspect_scale, $
			SET_VALUE = 'Scale '+STRLOWCASE(((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).heightset])+' to maximum of average'
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],2)
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],2)
		WIDGET_CONTROL, (*(*info).ctrlscp).detspect_label, SET_VALUE = ((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).heightset]+':'
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).ctrlsswitch).imrefdetspect], labels=['imrefdetspect']
END

PRO CRISPEX_DISPLAYS_DOPPLER_TOGGLE, event, NO_DRAW=no_draw
; Reference image window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_DOPPLER_TOGGLE', /IGNORE_LAST
	(*(*info).winswitch).showdop = event.SELECT
	IF (*(*info).winswitch).showdop THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Doppler image'
		CRISPEX_WINDOW, (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny, (*(*info).winids).root, title, doptlb, dopwid, (*(*info).winsizes).xdelta,(*(*info).winsizes).ydelta, $
			DRAWID = dopdrawid, DRAWBASE = dopdrawbase
		(*(*info).winids).doptlb = doptlb		&	(*(*info).winids).dopwid = dopwid	&	(*(*info).winids).dopdrawid = dopdrawid
		(*(*info).winids).dopdrawbase = dopdrawbase	&	(*(*info).winids).dopwintitle = title
		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
			CRISPEX_UPDATE_T, event
			CRISPEX_DRAW_DOPPLER, event
		ENDIF
		WIDGET_CONTROL, dopdrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS
		WIDGET_CONTROL, doptlb, SET_UVALUE = info
		XMANAGER, 'CRISPEX', doptlb, /NO_BLOCK
	ENDIF ELSE BEGIN
		(*(*info).dispswitch).drawdop = 0
		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
			CRISPEX_DRAW_SPECTRAL, event
			CRISPEX_DRAW_TIMESLICES, event
		ENDIF
		WIDGET_CONTROL, (*(*info).winids).doptlb, /DESTROY
		(*(*info).winids).doptlb = 0
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).dop_scaling_but, SENSITIVE = ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop)
	IF (*(*info).overlayswitch).mask THEN CRISPEX_MASK_BUTTONS_SET, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).doptlb,(*(*info).winids).dopwid,(*(*info).winids).dopdrawid], labels=['doptlb','dopwid','dopdrawid']
END

PRO CRISPEX_DISPLAYS_INT_MENU, event, set_but_array
; Sets up the intensity-time plot options menu
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_MENU'
	WIDGET_CONTROL,/HOURGLASS
	eventval = INDGEN((*(*info).dataparams).nlp)
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Intensity-time plot options', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	disp2 = WIDGET_BASE(disp, /COLUMN, /FRAME)
	sel_allnone = WIDGET_BASE(disp2, /ROW)
	sel_allnone_lab = WIDGET_LABEL(sel_allnone, VALUE = 'Plot selected diagnostics:', /ALIGN_LEFT)
	sel_allnone_buts = WIDGET_BASE(sel_allnone, /ROW, /EXCLUSIVE)
	(*(*info).ctrlsint).int_sel_all = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'All', EVENT_PRO = 'CRISPEX_DISPLAYS_INT_SEL_ALL')
	(*(*info).ctrlsint).int_sel_none = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'None', EVENT_PRO = 'CRISPEX_DISPLAYS_INT_SEL_NONE')
	IF ((*(*info).dataparams).nlp GT 10) THEN sel_opts = WIDGET_BASE(disp2, /COLUMN, Y_SCROLL_SIZE = 440) ELSE sel_opts = WIDGET_BASE(disp2, /COLUMN)
	uniq_cols = ((*(*info).intparams).collab_diagnostics)[UNIQ((*(*info).intparams).collab_diagnostics)]
	FOR i=0,(*(*info).dataparams).nlp-1 DO BEGIN
		sel_subopts = WIDGET_BASE(sel_opts, /ROW)
		sel_buts = WIDGET_BASE(sel_subopts, /NONEXCLUSIVE)
		name = 'int_sel_but_'+STRTRIM(i,2)
		but_val = ((*(*info).intparams).diagnostics)[i]
		sel_but = WIDGET_BUTTON(sel_buts, VALUE = but_val, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_DISPLAYS_INT_MENU_EVENT', UNAME = name)
		WIDGET_CONTROL, sel_but, SET_BUTTON = ((*(*(*info).intparams).sel_diagnostics)[i] EQ 1)
		lname = 'int_sel_line_'+STRTRIM(i,2)
		sel_line_list = WIDGET_COMBOBOX(sel_subopts, VALUE = ((*(*info).intparams).linlab_diagnostics)[0:5], UVALUE = eventval[i], /DYNAMIC_RESIZE, EVENT_PRO = 'CRISPEX_DISPLAYS_INT_SEL_LINE', UNAME = lname)
		WIDGET_CONTROL, sel_line_list, SET_COMBOBOX_SELECT = ((*(*(*info).intparams).lines_diagnostics)[i] MOD 6), SENSITIVE = ((*(*(*info).intparams).sel_diagnostics)[i] EQ 1) 
		cname = 'int_sel_cols_'+STRTRIM(i,2)
		sel_cols_list = WIDGET_COMBOBOX(sel_subopts, VALUE = uniq_cols, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_DISPLAYS_INT_SEL_COLS', UNAME = cname)
		WIDGET_CONTROL, sel_cols_list, SET_COMBOBOX_SELECT = WHERE(uniq_cols EQ ((*(*info).intparams).collab_diagnostics)[(*(*(*info).intparams).selcol_diagnostics)[i]]), $
			SENSITIVE = ((*(*(*info).intparams).sel_diagnostics)[i] EQ 1) 
	ENDFOR
	disp_label = WIDGET_LABEL(disp2, VALUE = 'Display options:', /ALIGN_LEFT)
	yrange_base = WIDGET_BASE(disp2, /ROW)
	lower_y_label = WIDGET_LABEL(yrange_base, VALUE = 'Lower y-value:')
	(*(*info).ctrlsint).lower_y_int_text = WIDGET_TEXT(yrange_base, VALUE = STRTRIM((*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_INT_LOW')
	upper_y_label = WIDGET_LABEL(yrange_base, VALUE = 'Upper y-value:')
	(*(*info).ctrlsint).upper_y_int_text = WIDGET_TEXT(yrange_base, VALUE = STRTRIM((*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s],2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_INT_UPP')
	trange_base = WIDGET_BASE(disp2, /ROW)
	lower_t_label = WIDGET_LABEL(trange_base, VALUE = 'Lower t-index:')
	(*(*info).ctrlsint).lower_t_int_text = WIDGET_TEXT(trange_base, VALUE = STRTRIM((*(*info).plotaxes).int_low_t,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_INT_T_LOW')
	upper_t_label = WIDGET_LABEL(trange_base, VALUE = 'Upper t-index:')
	(*(*info).ctrlsint).upper_t_int_text = WIDGET_TEXT(trange_base, VALUE = STRTRIM((*(*info).plotaxes).int_upp_t,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_INT_T_UPP')
	(*(*info).ctrlsint).reset_trange_but = WIDGET_BUTTON(disp2, VALUE = 'Reset temporal boundaries', EVENT_PRO = 'CRISPEX_DISPRANGE_INT_T_RESET', $
		SENSITIVE = (((*(*info).plotaxes).int_upp_t-(*(*info).plotaxes).int_low_t+1) NE (*(*info).dataparams).nt))
	lock_base = WIDGET_BASE(disp2, /ROW, /NONEXCLUSIVE)
	lock_t_range = WIDGET_BUTTON(lock_base, VALUE = 'Lock main to intensity-time temporal range', EVENT_PRO = 'CRISPEX_DISPRANGE_INT_LOCK_T')
	WIDGET_CONTROL, lock_t_range, SET_BUTTON = (*(*info).intparams).lock_t
	button_base = WIDGET_BASE(disp2, COLUMN=2, /GRID_LAYOUT, /ALIGN_CENTER)
	(*(*info).ctrlsint).int_sel_save = WIDGET_BUTTON(button_base, VALUE = 'Save selected', EVENT_PRO = 'CRISPEX_INT_SAVE')
	closebut = WIDGET_BUTTON(button_base, VALUE = 'Close', EVENT_PRO = 'CRISPEX_DISPLAYS_INT_MENU_CLOSE')
	IF (N_ELEMENTS(WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)) EQ (*(*info).dataparams).nlp) THEN WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_all, /SET_BUTTON ELSE $
		IF (TOTAL(WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)) EQ -1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_none, /SET_BUTTON
				WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_save, SENSITIVE = 0
		ENDIF
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).spxoffset, TLB_SET_YOFFSET = 0
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).intmenutlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).intmenutlb], labels=['intmenutlb']
END

PRO CRISPEX_DISPLAYS_INT_MENU_EVENT, event
; Handles the selection of diagnostics to be shown in the intensity-time plot
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_MENU_EVENT'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).intparams).sel_diagnostics)[eventval] = ( (*(*(*info).intparams).sel_diagnostics)[eventval] EQ 0) 
	lname = 'int_sel_line_'+STRTRIM(eventval,2)
	WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = lname), SENSITIVE = (*(*(*info).intparams).sel_diagnostics)[eventval]
	cname = 'int_sel_cols_'+STRTRIM(eventval,2)
	WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = cname), SENSITIVE = (*(*(*info).intparams).sel_diagnostics)[eventval]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).intparams).sel_diagnostics)[eventval]], labels=['Diagnostic ID','Diagnostic selected']
	CRISPEX_DISPLAYS_INT_BUTTON_CONDITION, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_DISPLAYS_INT_BUTTON_CONDITION, event
; Handles the update of buttons after selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_BUTTON_CONDITION'
	condition = WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_save, SENSITIVE = ((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))
	WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_none, SET_BUTTON = ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1)
	WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_all, SET_BUTTON = (((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).dataparams).nlp))
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).dataparams).nlp)),$
		ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1),((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))], labels=['All selected','None selected','Save enabled']
END

PRO CRISPEX_DISPLAYS_INT_MENU_CLOSE, event
; Handles the closing of the intensity versus time plot options menu and clean-up of display afterwards
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_MENU_CLOSE'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).intmenutlb,(*(*info).winids).inttlb], labels=['intmenutlb was','inttlb was']
	WIDGET_CONTROL, (*(*info).winids).intmenutlb, /DESTROY
	WIDGET_CONTROL, (*(*info).winids).inttlb, /DESTROY
	(*(*info).winids).intmenutlb = 0
	(*(*info).winids).inttlb = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).int_toggle_but, SET_BUTTON = 0
	(*(*info).winswitch).showint = 0
END

PRO CRISPEX_DISPLAYS_INT_SEL_ALL, event
; Handles selection of all intensity versus time plots
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_SEL_ALL'
	*(*(*info).intparams).sel_diagnostics = REPLICATE(1,(*(*info).dataparams).nlp)
	FOR i=0,(*(*info).dataparams).nlp-1 DO BEGIN
		name = 'int_sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 1
		lname = 'int_sel_line_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = lname), /SENSITIVE
		cname = 'int_sel_cols_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = cname), /SENSITIVE
	ENDFOR		
	WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_save, /SENSITIVE
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_SEL_COLS, event
; Handles selection of linestyle of intensity versus time plot
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_SEL_COLS'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	IF ( (*(*(*info).intparams).sel_diagnostics)[eventval] EQ 1) THEN (*(*(*info).intparams).selcol_diagnostics)[eventval] = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).intparams).selcol_diagnostics)[eventval]], labels=['Diagnostic ID','Color index selected']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_SEL_LINE, event
; Handles selection of linestyle of intensity versus time plot
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_SEL_LINE'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	IF ( (*(*(*info).intparams).sel_diagnostics)[eventval] EQ 1) THEN (*(*(*info).intparams).lines_diagnostics)[eventval] = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).intparams).lines_diagnostics)[eventval]], labels=['Diagnostic ID','Linestyle selected']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_SEL_NONE, event
; Handles selection of none intensity versus time plots
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_SEL_NONE'
	*(*(*info).intparams).sel_diagnostics = REPLICATE(0,(*(*info).dataparams).nlp)
	FOR i=0,(*(*info).dataparams).nlp-1 DO BEGIN
		name = 'int_sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 0
		lname = 'int_sel_line_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = lname), SENSITIVE = 0
		cname = 'int_sel_cols_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = cname), SENSITIVE = 0
	ENDFOR		
	WIDGET_CONTROL, (*(*info).ctrlsint).int_sel_save, SENSITIVE = 0
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_RESIZE, event						
; Intensity versus time window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).intxres, (*(*info).winsizes).intyres, (*(*info).plotpos).intxmargin_init, (*(*info).plotpos).intxwall_init, $
		intxres, intyres, intwidth, intheight, intx0, intx1, inty0, inty1, ERROR=error, /GOLDEN
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).intxres = intxres		& 	(*(*info).winsizes).intyres = intyres
		(*(*info).plotpos).intx0 = intx0		&	(*(*info).plotpos).intx1 = intx1
		(*(*info).plotpos).inty0 = inty0		&	(*(*info).plotpos).inty1 = inty1
		(*(*info).plotaxes).intxticklen = (*(*info).plotaxes).ticklen / intheight
		(*(*info).plotaxes).intyticklen = (*(*info).plotaxes).ticklen / intwidth
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).intxres,(*(*info).winsizes).intyres,(*(*info).plotpos).intx0,(*(*info).plotpos).intx1,$
		(*(*info).plotpos).inty0,(*(*info).plotpos).inty1], labels=['error','intxres','intyres','intx0','intx1','inty0','inty1']
	WIDGET_CONTROL, (*(*info).winids).intdrawid, DRAW_XSIZE = (*(*info).winsizes).intxres, DRAW_YSIZE = (*(*info).winsizes).intyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_TOGGLE, event, NO_DRAW=no_draw
; Intensity versus time window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_INT_TOGGLE'
	(*(*info).winswitch).showint = event.SELECT
	IF (*(*info).winswitch).showint THEN BEGIN
		CRISPEX_DISPLAYS_INT_MENU, event
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Intensity-time plot'
		CRISPEX_WINDOW, (*(*info).winsizes).intxres, (*(*info).winsizes).intyres, (*(*info).winids).root, title, inttlb, intwid, (*(*info).winsizes).lsxoffset, 0, DRAWID = intdrawid, $
			RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_INT_RESIZE'
		(*(*info).winids).inttlb = inttlb	&	(*(*info).winids).intwid = intwid	&	(*(*info).winids).intdrawid = intdrawid
		(*(*info).winids).intwintitle = title
		WIDGET_CONTROL, (*(*info).winids).inttlb, SET_UVALUE = info
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).inttlb,(*(*info).winids).intwid,(*(*info).winids).intdrawid], labels=['inttlb','intwid','intdrawid']
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_INT, event
	ENDIF ELSE CRISPEX_DISPLAYS_INT_MENU_CLOSE, event
END

PRO CRISPEX_DISPLAYS_LOOPSLAB_GET, event
; Loopslab window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_LOOPSLAB_GET'
	IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN
		CRISPEX_DISPLAYS_LOOPSLAB, event,/NO_DRAW
		CRISPEX_DISPLAYS_REFLOOPSLAB, event
	ENDIF ELSE CRISPEX_DISPLAYS_LOOPSLAB, event
END

PRO CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event					
; Updates loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES'
	WSET, (*(*info).winids).loopwid
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).loopx0,(*(*info).plotpos).loopy0, (*(*info).plotpos).loopx1,(*(*info).plotpos).loopy1], $
		YTICKLEN = (*(*info).plotaxes).loopyticklen, XTICKLEN = (*(*info).plotaxes).loopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).loopwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_LOOPSLAB_RESIZE, event
; Loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_LOOPSLAB_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).loopxres, (*(*info).winsizes).loopyres, (*(*info).plotpos).loopxmargin_init, (*(*info).plotpos).loopxwall_init, $
		loopxres, loopyres, loopwidth, loopheight, loopx0, loopx1, loopy0, loopy1, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).loopxres = loopxres			& 	(*(*info).winsizes).loopyres = loopyres
		(*(*info).plotpos).loopx0 = loopx0			&	(*(*info).plotpos).loopx1 = loopx1
		(*(*info).plotpos).loopy0 = loopy0			&	(*(*info).plotpos).loopy1 = loopy1
		(*(*info).plotpos).loopxplspw = loopx1 - loopx0		&	(*(*info).plotpos).loopyplspw = loopy1 - loopy0
		(*(*info).plotaxes).loopxticklen = -1 * (*(*info).plotaxes).ticklen / loopheight
		(*(*info).plotaxes).loopyticklen = -1 * (*(*info).plotaxes).ticklen / loopwidth
		(*(*info).dispparams).loopnlxreb = (*(*info).plotpos).loopxplspw * (*(*info).winsizes).loopxres 
		(*(*info).dispparams).loopntreb = (*(*info).plotpos).loopyplspw * (*(*info).winsizes).loopyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).loopxres,(*(*info).winsizes).loopyres,(*(*info).plotpos).loopx0,(*(*info).plotpos).loopx1,$
		(*(*info).plotpos).loopy0,(*(*info).plotpos).loopy1], labels=['error','loopxres','loopyres','loopx0','loopx1','loopy0','loopy1']
	WIDGET_CONTROL, (*(*info).winids).loopdrawid, DRAW_XSIZE = (*(*info).winsizes).loopxres, DRAW_YSIZE = (*(*info).winsizes).loopyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event
	CRISPEX_DRAW_LOOPSLAB, event
END

PRO CRISPEX_DISPLAYS_LOOPSLAB, event, NO_DRAW=no_draw
; Loopslab window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_LOOPSLAB'
	(*(*info).winswitch).showloop = 1
	WIDGET_CONTROL,/HOURGLASS
	CRISPEX_LOOP_GET_PATH, event
	CRISPEX_LOOP_GET_SLAB, event
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': T-slice along loop'
	CRISPEX_WINDOW, (*(*info).winsizes).loopxres, (*(*info).winsizes).loopyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta, $
		((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = loopdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_LOOPSLAB_RESIZE'
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).loopx0,(*(*info).plotpos).loopy0,(*(*info).plotpos).loopx1,(*(*info).plotpos).loopy1],$
		YTICKLEN = (*(*info).plotaxes).loopyticklen, XTICKLEN = (*(*info).plotaxes).loopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop',$
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	(*(*info).winids).looptlb = tlb		&	(*(*info).winids).loopwid = wid		&	(*(*info).winids).loopdrawid = loopdrawid
	(*(*info).winids).loopwintitle = title 
	WIDGET_CONTROL, (*(*info).winids).looptlb, SET_UVALUE = info
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		CRISPEX_UPDATE_LP, event
		CRISPEX_ZOOM_LOOP, event
		CRISPEX_DRAW, event
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).save_loop_pts, SENSITIVE = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).timeslicemenu, SENSITIVE = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).looptlb,(*(*info).winids).loopwid,(*(*info).winids).loopdrawid], labels=['looptlb','loopwid','loopdrawid']
END

PRO CRISPEX_DISPLAYS_REFLOOPSLAB_REPLOT_AXES, event					
; Updates reference loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFLOOPSLAB_REPLOT_AXES'
	WSET, (*(*info).winids).refloopwid
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopy0, (*(*info).plotpos).refloopx1,(*(*info).plotpos).refloopy1], $
		YTICKLEN = (*(*info).plotaxes).refloopyticklen, XTICKLEN = (*(*info).plotaxes).refloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refloopwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_REFLOOPSLAB_RESIZE, event
; Reference loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFLOOPSLAB_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).refloopxres, (*(*info).winsizes).refloopyres, (*(*info).plotpos).refloopxmargin_init, (*(*info).plotpos).refloopxwall_init, $
		refloopxres, refloopyres, refloopwidth, refloopheight, refloopx0, refloopx1, refloopy0, refloopy1, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).refloopxres = refloopxres			& 	(*(*info).winsizes).refloopyres = refloopyres
		(*(*info).plotpos).refloopx0 = refloopx0			&	(*(*info).plotpos).refloopx1 = refloopx1
		(*(*info).plotpos).refloopy0 = refloopy0			&	(*(*info).plotpos).refloopy1 = refloopy1
		(*(*info).plotpos).refloopxplspw = refloopx1 - refloopx0	&	(*(*info).plotpos).refloopyplspw = refloopy1 - refloopy0
		(*(*info).plotaxes).refloopxticklen = -1 * (*(*info).plotaxes).ticklen / refloopheight
		(*(*info).plotaxes).refloopyticklen = -1 * (*(*info).plotaxes).ticklen / refloopwidth
		(*(*info).dispparams).refloopnlxreb = (*(*info).plotpos).refloopxplspw * (*(*info).winsizes).refloopxres 
		(*(*info).dispparams).refloopntreb = (*(*info).plotpos).refloopyplspw * (*(*info).winsizes).refloopyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).refloopxres,(*(*info).winsizes).refloopyres,(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopx1,$
		(*(*info).plotpos).refloopy0,(*(*info).plotpos).refloopy1], labels=['error','refloopxres','refloopyres','refloopx0','refloopx1','refloopy0','refloopy1']
	WIDGET_CONTROL, (*(*info).winids).refloopdrawid, DRAW_XSIZE = (*(*info).winsizes).refloopxres, DRAW_YSIZE = (*(*info).winsizes).refloopyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_REFLOOPSLAB_REPLOT_AXES, event
	CRISPEX_DRAW_REFLOOPSLAB, event
END

PRO CRISPEX_DISPLAYS_REFLOOPSLAB, event, NO_DRAW=no_draw
; Reference, loopslab window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFLOOPSLAB'
	(*(*info).winswitch).showrefloop = 1
	WIDGET_CONTROL,/HOURGLASS
	CRISPEX_LOOP_GET_REFSLAB, event		
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Reference T-slice along loop'
	CRISPEX_WINDOW, (*(*info).winsizes).refloopxres, (*(*info).winsizes).refloopyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta, $
		((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = refloopdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFLOOPSLAB_RESIZE'
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopy0,(*(*info).plotpos).refloopx1,(*(*info).plotpos).refloopy1],$
		YTICKLEN = (*(*info).plotaxes).refloopyticklen, XTICKLEN = (*(*info).plotaxes).refloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop',$
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	(*(*info).winids).reflooptlb = tlb		&	(*(*info).winids).refloopwid = wid		&	(*(*info).winids).refloopdrawid = refloopdrawid
	(*(*info).winids).refloopwintitle = title 
	WIDGET_CONTROL, (*(*info).winids).reflooptlb, SET_UVALUE = info
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		CRISPEX_UPDATE_LP, event
		CRISPEX_ZOOM_LOOP, event
		CRISPEX_DRAW, event
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).reflooptlb,(*(*info).winids).refloopwid,(*(*info).winids).refloopdrawid], labels=['reflooptlb','refloopwid','refloopdrawid']
END

PRO CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE, event 
; Sets up the parameter overview window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE', /IGNORE_LAST
	(*(*info).winswitch).showparam = event.SELECT
	IF ((*(*info).winswitch).showparam EQ 1) THEN BEGIN
		base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Parameters overview', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
		disp_base = WIDGET_BASE(base, /ROW)
    ; Column 1 of parameters overview containing cursor x,y and zoomfactor
		disp = WIDGET_BASE(disp_base, /COLUMN)
    ; Cursor x info
		x_coord_base = WIDGET_BASE(disp, /ROW)
		x_coord_txt = WIDGET_LABEL(x_coord_base, VALUE = 'Cursor x:')
		(*(*info).ctrlsparam).x_coord_val = WIDGET_LABEL(x_coord_base, VALUE = STRTRIM((*(*info).dataparams).x,2), /DYNAMIC_RESIZE)
    ; Cursor y info
		y_coord_base = WIDGET_BASE(disp, /ROW)
		y_coord_txt = WIDGET_LABEL(y_coord_base, VALUE = 'Cursor y:')
		(*(*info).ctrlsparam).y_coord_val = WIDGET_LABEL(y_coord_base, VALUE = STRTRIM((*(*info).dataparams).y,2), /DYNAMIC_RESIZE)
    ; Zommfactor info
		zoom_base = WIDGET_BASE(disp, /ROW)
		zoom_txt = WIDGET_LABEL(zoom_base, VALUE = 'Zoomfactor:')
		(*(*info).ctrlsparam).zoom_val = WIDGET_LABEL(zoom_base, VALUE = STRTRIM(LONG((*(*info).zooming).factor),2), /DYNAMIC_RESIZE)
    ; Column 2 containing main spectral/height info, including Doppler
		disp2 = WIDGET_BASE(disp_base, /COLUMN)
    ; Main spectral info
		lp_coord_base = WIDGET_BASE(disp2, /ROW)
		lp_coord_txt = WIDGET_LABEL(lp_coord_base, VALUE = ((*(*info).paramparams).wav_h)[(*(*info).plotswitch).heightset]+' index:')
		(*(*info).ctrlsparam).lp_coord_val = WIDGET_LABEL(lp_coord_base, VALUE = STRTRIM(LONG((*(*info).dataparams).lp),2), /DYNAMIC_RESIZE)
		act_lp_base = WIDGET_BASE(disp2, /ROW)
		act_lp_txt = WIDGET_LABEL(act_lp_base, $
      VALUE = ((*(*info).paramparams).wav_h)[(*(*info).plotswitch).heightset]+' ['+$
              ((*(*info).dataparams).lpunit)[0]+']:')
    ; Main Doppler info
  	IF ((*(*info).plotswitch).v_dop_set OR (*(*info).plotswitch).heightset) THEN $
      (*(*info).ctrlsparam).act_lp_val = WIDGET_LABEL(act_lp_base, $
        VALUE = STRTRIM(STRING((*(*info).dataparams).lps[(*(*info).dataparams).lp],$
        FORMAT = '(3(F9.1,x))'),2),/DYNAMIC_RESIZE) $
    ELSE (*(*info).ctrlsparam).act_lp_val = WIDGET_LABEL(act_lp_base, VALUE = 'N/A')
    IF ((*(*info).plotswitch).heightset EQ 0) THEN BEGIN
  		v_dop_base = WIDGET_BASE(disp2, /ROW)
  		v_dop_txt = WIDGET_LABEL(v_dop_base, VALUE = 'Doppler velocity (km/s):')
  		IF (*(*info).plotswitch).v_dop_set THEN BEGIN
  			(*(*info).ctrlsparam).v_dop_val = WIDGET_LABEL(v_dop_base, VALUE = STRTRIM(STRING((*(*info).plotaxes).v_dop[(*(*info).dataparams).lp],FORMAT='(3(F9.2,x))'),2), /DYNAMIC_RESIZE)
  		ENDIF ELSE IF (*(*info).plotswitch).heightset THEN BEGIN
  			(*(*info).ctrlsparam).v_dop_val = WIDGET_LABEL(v_dop_base, VALUE = 'N/A')
  		ENDIF ELSE BEGIN
  			(*(*info).ctrlsparam).v_dop_val = WIDGET_LABEL(v_dop_base, VALUE = 'N/A')
  		ENDELSE
    ENDIF
    ; Column 3 (if reference present) containing reference spectral/height info, including Doppler
		IF ((*(*info).dataparams).refnlp GT 1) THEN BEGIN
			dispref = WIDGET_BASE(disp_base, /COLUMN)
      ; Reference spectral info
			lp_ref_coord_base = WIDGET_BASE(dispref, /ROW)
			lp_ref_coord_txt = WIDGET_LABEL(lp_ref_coord_base, VALUE = 'Reference '+STRLOWCASE(((*(*info).paramparams).wav_h)[(*(*info).plotswitch).refheightset])+' index:')
			(*(*info).ctrlsparam).lp_ref_coord_val = WIDGET_LABEL(lp_ref_coord_base, VALUE = STRTRIM(LONG((*(*info).dataparams).lp_ref),2), /DYNAMIC_RESIZE)
			act_lp_ref_base = WIDGET_BASE(dispref, /ROW)
			act_lp_ref_txt = WIDGET_LABEL(act_lp_ref_base, VALUE = 'Reference '+$
        STRLOWCASE(((*(*info).paramparams).wav_h)[(*(*info).plotswitch).refheightset])+' ['+$
        ((*(*info).dataparams).lpunit)[1]+']:')
      ; Reference Doppler info
    	IF ((*(*info).plotswitch).v_dop_set_ref OR (*(*info).plotswitch).refheightset) THEN $
        (*(*info).ctrlsparam).act_lp_ref_val = WIDGET_LABEL(act_lp_ref_base, $
          VALUE = STRTRIM(STRING((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref],$
          FORMAT = '(3(F9.1,x))'),2),/DYNAMIC_RESIZE) $
      ELSE (*(*info).ctrlsparam).act_lp_ref_val = WIDGET_LABEL(act_lp_ref_base, VALUE = 'N/A')
      IF ((*(*info).plotswitch).heightset EQ 0) THEN BEGIN
  			v_dop_ref_base = WIDGET_BASE(dispref, /ROW)
  			v_dop_ref_txt = WIDGET_LABEL(v_dop_ref_base, VALUE = 'Reference Doppler velocity (km/s):')
  			IF (*(*info).plotswitch).v_dop_set_ref THEN BEGIN
  				(*(*info).ctrlsparam).v_dop_ref_val = WIDGET_LABEL(v_dop_ref_base, VALUE = STRTRIM(STRING((*(*info).plotaxes).v_dop_ref[(*(*info).dataparams).lp_ref],FORMAT='(3(F9.2,x))'),2), /DYNAMIC_RESIZE)
  			ENDIF ELSE IF (*(*info).plotswitch).refheightset THEN BEGIN
  				(*(*info).ctrlsparam).v_dop_ref_val = WIDGET_LABEL(v_dop_ref_base, VALUE = 'N/A')
  			ENDIF ELSE BEGIN
  				(*(*info).ctrlsparam).v_dop_ref_val = WIDGET_LABEL(v_dop_ref_base, VALUE = 'N/A')
  			ENDELSE
      ENDIF
		ENDIF
    ; Column 4 (if reference present, else column 3) containing channel and time info
		disp3 = WIDGET_BASE(disp_base, /COLUMN)
    ; Channel info
		stokes_base = WIDGET_BASE(disp3, /ROW)
		stokes_txt = WIDGET_LABEL(stokes_base, VALUE = 'Stokes:')
		(*(*info).ctrlsparam).stokes_val = WIDGET_LABEL(stokes_base, VALUE = STRTRIM(((*(*info).stokesparams).labels)[(*(*info).dataparams).s],2), /DYNAMIC_RESIZE)
    ; Time info
		t_coord_base = WIDGET_BASE(disp3, /ROW)
		t_coord_txt = WIDGET_LABEL(t_coord_base, VALUE = 'Time index:')
		(*(*info).ctrlsparam).t_coord_val = WIDGET_LABEL(t_coord_base, VALUE = STRTRIM((*(*info).dataparams).t,2), /DYNAMIC_RESIZE)
		act_t_base = WIDGET_BASE(disp3, /ROW)
;		act_t_txt = WIDGET_LABEL(act_t_base, VALUE = 'Actual time (s):')
		act_t_txt = WIDGET_LABEL(act_t_base, VALUE = 'Actual '+(*(*info).plottitles).spytitle+':')
		IF (*(*info).paramswitch).dt_set THEN (*(*info).ctrlsparam).act_t_val = WIDGET_LABEL(act_t_base, VALUE = STRTRIM(STRING((*(*info).dataparams).t * (*(*info).plotaxes).dt, FORMAT='(3(F9.2,x))'),2), $
			/DYNAMIC_RESIZE) ELSE (*(*info).ctrlsparam).act_t_val = WIDGET_LABEL(act_t_base, VALUE = 'N/A')
    ; Column 5 (if reference present, else column 4) containing image values below cursor
    ; Main image value info
		IF (*(*info).paramswitch).img_get THEN BEGIN
			disp4 = WIDGET_BASE(disp_base, /COLUMN)
			img_base = WIDGET_BASE(disp4, /ROW)
			img_label = WIDGET_LABEL(img_base, VALUE = 'Main cube value ['+$
        ((*(*info).dataparams).bunit)[0]+'] :')
			(*(*info).ctrlsparam).img_val = WIDGET_LABEL(img_base, VALUE = STRTRIM(STRING((*(*info).paramparams).img_get_val,FORMAT='(3(F9.2,x))'),2),/DYNAMIC_RESIZE)
		ENDIF
    ; Reference image value info
		IF (*(*info).paramswitch).ref_get THEN BEGIN
			IF ((*(*info).paramswitch).img_get EQ 0) THEN disp4 = WIDGET_BASE(disp_base, /COLUMN)
			ref_base = WIDGET_BASE(disp4, /ROW)
			ref_label = WIDGET_LABEL(ref_base, VALUE = 'Reference cube value ['+$
        ((*(*info).dataparams).bunit)[1]+'] :')
      (*(*info).ctrlsparam).ref_val = WIDGET_LABEL(ref_base, VALUE = STRTRIM(STRING((*(*info).paramparams).ref_get_val,FORMAT='(3(F9.2,x))'),2),/DYNAMIC_RESIZE)
		ENDIF
		IF ((*(*info).paramswitch).img_get OR (*(*info).paramswitch).ref_get) THEN BEGIN
			IF ((*(*info).paramswitch).ref_get EQ 0) THEN BEGIN
				mainlabel = 'main cube value:' 
				extra_label = ''
			ENDIF ELSE BEGIN
				mainlabel = 'cube values, main:'
				extra_label = ','
			ENDELSE
			IF ((*(*info).paramswitch).img_get EQ 0) THEN mainlabel = 'reference cube value:' ELSE label = 'reference:'
			sc_base = WIDGET_BASE(disp4, /ROW)
			sc_label = WIDGET_LABEL(sc_base, VALUE = 'Scaled '+mainlabel)
			IF (*(*info).paramswitch).img_get THEN (*(*info).ctrlsparam).imgsc_val = WIDGET_LABEL(sc_base, VALUE = STRTRIM(STRING((*(*info).paramparams).imgsc_get_val,FORMAT='(3(F9.2,x))'),2)+extra_label,$
				/DYNAMIC_RESIZE)
			IF (*(*info).paramswitch).ref_get THEN BEGIN
				IF (*(*info).paramswitch).img_get THEN refsc_label = WIDGET_LABEL(sc_base, VALUE = label)
				(*(*info).ctrlsparam).refsc_val = WIDGET_LABEL(sc_base, VALUE = STRTRIM(STRING((*(*info).paramparams).refsc_get_val,FORMAT='(3(F9.2,x))'),2),/DYNAMIC_RESIZE)
			ENDIF
		ENDIF
		WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 0, TLB_SET_YOFFSET = (*(*info).winsizes).xywiny + (*(*info).winsizes).ydelta
		WIDGET_CONTROL, base, SET_UVALUE = info
		XMANAGER, 'CRISPEX', base, /NO_BLOCK
		(*(*info).winids).paramtlb = base
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).paramtlb, /DESTROY
		(*(*info).winids).paramtlb = 0
		(*(*info).winswitch).showparam = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).looptlb], labels=['looptlb']
END

PRO CRISPEX_DISPLAYS_PLOT_RESIZE, event, new_xres_tmp, new_yres_tmp, init_xres, init_yres, init_xmargin, init_xwall, new_xres, new_yres, new_width, new_height, $
	x0, x1, y0, y1, v_dop_set, INX0=inx0, INX1=inx1, INY0=iny0, INY1=iny1, ERROR=error, GOLDEN=golden, ACTUAL_RESIZE=actual_resize, DETSPECT=detspect, STOKES_SELECT=stokes_select
; Handles the display plot resizing
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_PLOT_RESIZE'
	IF (N_ELEMENTS(v_dop_set) NE 1) THEN v_dop_set = 0
	IF KEYWORD_SET(ACTUAL_RESIZE) THEN BEGIN												; Used for display of loopslices
		new_xres = new_xres_tmp
		new_yres = new_yres_tmp
		new_margin = init_xmargin/new_xres
		new_wall = init_xwall/new_xres
		new_width = (1. - (new_margin + new_wall))
		IF ((v_dop_set EQ 1) OR ((*(*info).dataparams).ns GT 1)) THEN new_height = (1. - (new_margin * 2.) * new_xres/new_yres) ELSE new_height = (1. - (new_margin + new_wall) * new_xres/new_yres)
		x0 = new_margin
		x1 = x0 + new_width
		y0 = new_margin * new_xres/new_yres
		y1 = y0 + new_height; * new_xres/new_yres 
	ENDIF ELSE BEGIN															; Used for regular plots
		IF KEYWORD_SET(DETSPECT) THEN BEGIN												; If considering main detailed spectrum
			curns = TOTAL((*(*info).stokesparams).select_sp)
			IF KEYWORD_SET(STOKES_SELECT) THEN BEGIN
				prevns = TOTAL((*(*info).stokesparams).prev_select_sp)
				new_xres_tmp = init_xres											; Default: no change in xsize
				new_yres_tmp = init_yres											; Default: no change in ysize
				IF (v_dop_set EQ 1) THEN offset = init_xmargin ELSE offset = init_xwall
				IF (prevns GT curns) THEN BEGIN											; If reducing the selected Stokes
					IF (curns EQ 2) THEN new_yres_tmp = (init_yres + offset - v_dop_set*init_xwall) / 2. ELSE $		; If new number is 2, prev was 3
						IF (curns EQ 1) THEN new_xres_tmp = (init_xres + init_xwall) / 2. 				;If new number is 1, prev was 2
				ENDIF ELSE BEGIN												; If increasing the selected Stokes
					IF (curns EQ 3) THEN new_yres_tmp = 2. * init_yres - offset + v_dop_set*init_xwall ELSE $		; If new number is 3, prev was 2
						IF (curns EQ 2) THEN new_xres_tmp = 2. * init_xres - init_xwall					; If new number is 2, prev was 1
				ENDELSE
			ENDIF
			IF (curns LE 2) THEN BEGIN
				npanels = curns	&	cols = curns
				rowarr = REPLICATE(0,curns)
			ENDIF ELSE BEGIN
				npanels = 4	&	cols = 2
				rowarr = [1,1,0,0]
			ENDELSE
			rows = CEIL(npanels / FLOAT(cols))
			x0 = FLTARR(npanels)
			x1 = FLTARR(npanels)
			y0 = FLTARR(npanels)
			y1 = FLTARR(npanels)
		ENDIF ELSE BEGIN														; All other plot windows
			rows = 1	&	cols = 1
		ENDELSE
		dx = ABS(new_xres_tmp - init_xres)
		dy = ABS(new_yres_tmp - init_yres)
		IF (dx GT dy) THEN BEGIN
			new_margin = init_xmargin/new_xres_tmp
			new_wall = init_xwall/new_xres_tmp
			new_xres = new_xres_tmp
			new_width = (1 - (cols*new_margin + new_wall))/FLOAT(cols)
			IF KEYWORD_SET(GOLDEN) THEN new_height = new_width * 2D / (1 + SQRT(5)) ELSE BEGIN
				IF (v_dop_set EQ 1) THEN new_height = (1. - (new_margin * 2.)) ELSE new_height = (1. - (new_margin + new_wall))
			ENDELSE
			IF (v_dop_set EQ 1) THEN new_yres = ((rows+1) * new_margin + rows*new_height + (rows-1)*new_wall) * new_xres ELSE new_yres = (new_wall + rows*new_height + rows*new_margin) * new_xres
			IF KEYWORD_SET(DETSPECT) THEN BEGIN
				x0 = new_margin * new_xres/new_xres + (INDGEN(npanels) MOD cols) * (new_width + new_margin) * new_xres/new_xres
				x1 = x0 + new_width * new_xres/new_xres
				y0 = new_margin * new_xres/new_yres + rowarr * (new_height + new_margin + v_dop_set*new_wall) * new_xres/new_yres
				y1 = y0 + new_height * new_xres/new_yres
			ENDIF ELSE BEGIN
				x0 = new_margin 
				x1 = x0 + new_width
				y0 = new_margin * new_xres/new_yres
				y1 = y0 + new_height * new_xres/new_yres 
			ENDELSE
		ENDIF ELSE IF (dx LT dy) THEN BEGIN
			new_margin = init_xmargin/new_yres_tmp
			new_wall = init_xwall/new_yres_tmp
			new_yres = new_yres_tmp
			IF (v_dop_set EQ 1) THEN new_height = (1. - ((rows+1)*new_margin + (rows-1)*new_wall))/FLOAT(rows) ELSE new_height = (1. - (rows*new_margin + new_wall))/FLOAT(rows)
			IF KEYWORD_SET(GOLDEN) THEN new_width = new_height / 2D * (1 + SQRT(5)) ELSE new_width = (1. - (new_margin + new_wall))
			new_xres = (cols*new_margin + cols*new_width + new_wall) * new_yres
			IF KEYWORD_SET(DETSPECT) THEN BEGIN
				x0 = new_margin * new_yres/new_xres + (INDGEN(npanels) MOD cols) * (new_width + new_margin) * new_yres/new_xres
				x1 = x0 + new_width * new_yres/new_xres
				y0 = new_margin + rowarr * (new_height + new_margin + v_dop_set*new_wall)
				y1 = y0 + new_height
			ENDIF ELSE BEGIN
				x0 = new_margin * new_yres/new_xres
				x1 = x0 + new_width * new_yres/new_xres
				y0 = new_margin
				y1 = y0 + new_height
			ENDELSE
		ENDIF ELSE BEGIN															; If no change in size (Stokes select)
			x0 = inx0	&	x1 = inx1
			y0 = iny0	&	y1 = iny1
			new_xres = new_xres_tmp
			new_yres = new_yres_tmp
			new_width = inx1[0] - inx0[0]
			new_height = iny1[0] - iny0[0]
		ENDELSE
	ENDELSE
	dxpl = x1[0] - x0[0]	&	dypl = y1[0] - y0[0]
	IF ((dxpl LE 0) OR (dypl LE 0) OR (x0[0] LE 0) OR (y0[0] LE 0)) THEN error = 1 ELSE error = 0
END

PRO CRISPEX_DISPLAYS_LS_RESIZE, event, STOKES_SELECT=stokes_select
; Detailed spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_LS_RESIZE'
	IF KEYWORD_SET(STOKES_SELECT) THEN BEGIN
		newlsxres = (*(*info).winsizes).lsxres	&	newlsyres = (*(*info).winsizes).lsyres
	ENDIF ELSE BEGIN
		newlsxres = event.X		&	newlsyres = event.Y
	ENDELSE
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, newlsxres, newlsyres, (*(*info).winsizes).lsxres, (*(*info).winsizes).lsyres, (*(*info).plotpos).lsxmargin_init, (*(*info).plotpos).lsxwall_init, lsxres, lsyres, lswidth, lsheight, $
		lsx0, lsx1, lsy0, lsy1, (*(*info).plotswitch).v_dop_set, INX0=(*(*info).plotpos).lsx0, INX1=(*(*info).plotpos).lsx1, INY0=(*(*info).plotpos).lsy0, INY1=(*(*info).plotpos).lsy1, ERROR=error, $
		/GOLDEN, /DETSPECT, STOKES_SELECT=stokes_select
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).lsxres = lsxres	& 	(*(*info).winsizes).lsyres = lsyres
		(*(*info).plotpos).lsx0 = lsx0		&	(*(*info).plotpos).lsx1 = lsx1
		(*(*info).plotpos).lsy0 = lsy0		&	(*(*info).plotpos).lsy1 = lsy1
		(*(*info).plotaxes).lsxticklen = (*(*info).plotaxes).ticklen / lsheight
		(*(*info).plotaxes).lsyticklen = (*(*info).plotaxes).ticklen / lswidth
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN BEGIN
		nstokes_sel = TOTAL((*(*info).stokesparams).select_sp)
		CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).lsxres,(*(*info).winsizes).lsyres,(*(*info).plotpos).lsx0,(*(*info).plotpos).lsx1,(*(*info).plotpos).lsy0,(*(*info).plotpos).lsy1], $
			labels=['error','lsxres','lsyres',REPLICATE('lsx0',nstokes_sel),REPLICATE('lsx1',nstokes_sel),REPLICATE('lsy0',nstokes_sel),REPLICATE('lsy1',nstokes_sel)]
	ENDIF
	WIDGET_CONTROL, (*(*info).winids).lsdrawid, DRAW_XSIZE = (*(*info).winsizes).lsxres, DRAW_YSIZE = (*(*info).winsizes).lsyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DRAW_LS, event
END

PRO CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE, event
; Sets the playback mode to blink of main and reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE', /IGNORE_LAST
	(*(*info).winswitch).showimref = event.SELECT
	(*(*info).pbparams).imrefmode = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, SENSITIVE=ABS((*(*info).pbparams).imrefmode-1)
	IF (*(*info).winswitch).showimref THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Main vs. reference image blink'
		CRISPEX_WINDOW, (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny, (*(*info).winids).root, title, imreftlb, imrefwid, $
			(*(*info).winsizes).xdelta,(*(*info).winsizes).ydelta, DRAWID = imrefdrawid, DRAWBASE = imrefdrawbase
		(*(*info).winids).imreftlb = imreftlb		&	(*(*info).winids).imrefwid = imrefwid	&	(*(*info).winids).imrefdrawid = imrefdrawid
		(*(*info).winids).imrefdrawbase = imrefdrawbase	&	(*(*info).winids).imrefwintitle = title
		IF ((*(*info).feedbparams).count_pbstats EQ 0) THEN (*(*info).feedbparams).pbstats = SYSTIME(/SECONDS)
		WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
		WIDGET_CONTROL, imrefdrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS
		WIDGET_CONTROL, imreftlb, SET_UVALUE = info
		XMANAGER, 'CRISPEX', imreftlb, /NO_BLOCK
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).imreftlb, /DESTROY
		(*(*info).scaling).imrefscaling = 0
		(*(*info).winids).imreftlb = 0
		IF ((*(*info).winids).feedbacktlb NE 0) THEN BEGIN
			(*(*info).feedbparams).count_pbstats = 0
			WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, /SENSITIVE
		ENDIF
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).imreftlb,(*(*info).winids).imrefwid,(*(*info).winids).imrefdrawid], labels=['imreftlb','imrefwid','imrefdrawid']
END

PRO CRISPEX_DISPLAYS_IMREF_LS_TOGGLE, event, NO_DRAW=no_draw
; Detailed spectrum window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_IMREF_LS_TOGGLE', /IGNORE_LAST
	IF (*(*info).ctrlsswitch).imrefdetspect THEN BEGIN	; For reference detailed spectrum window
		IF ((*(*info).winswitch).showrefls EQ 0) THEN BEGIN
			title = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).reflswintitle)[(*(*info).plotswitch).refheightset]
			CRISPEX_WINDOW, (*(*info).winsizes).reflsxres, (*(*info).winsizes).reflsyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).lsxoffset, $
				(*(*info).winswitch).showsp * (*(*info).winsizes).ydelta, DRAWID = lsdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFLS_RESIZE'
			(*(*info).winids).reflstlb = tlb		&	(*(*info).winids).reflswid = wid	&	(*(*info).winswitch).showrefls = 1
			(*(*info).winids).reflsdrawid = lsdrawid	&	(*(*info).winids).reflswintitle = title
			WIDGET_CONTROL, (*(*info).winids).reflstlb, SET_UVALUE = info
			IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_REFLS, event
		ENDIF ELSE BEGIN
			WIDGET_CONTROL, (*(*info).winids).reflstlb, /DESTROY
			(*(*info).winids).reflstlb = 0
			(*(*info).winswitch).showrefls = 0
		ENDELSE
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).reflstlb,(*(*info).winids).reflswid,(*(*info).winids).reflsdrawid], labels=['reflstlb','reflswid','reflsdrawid']
	ENDIF ELSE BEGIN
		IF ((*(*info).winswitch).showls EQ 0) THEN BEGIN
			title = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).lswintitle)[(*(*info).plotswitch).heightset]
			CRISPEX_WINDOW, (*(*info).winsizes).lsxres, (*(*info).winsizes).lsyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).lsxoffset, 0, DRAWID = lsdrawid, RESIZING = 1, $
				RES_HANDLER = 'CRISPEX_DISPLAYS_LS_RESIZE'
			(*(*info).winids).lstlb = tlb		&	(*(*info).winids).lswid = wid	&	(*(*info).winswitch).showls = 1
			(*(*info).winids).lsdrawid = lsdrawid	&	(*(*info).winids).lswintitle = title
			WIDGET_CONTROL, (*(*info).winids).lstlb, SET_UVALUE = info
			IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_LS, event
		ENDIF ELSE BEGIN
			WIDGET_CONTROL, (*(*info).winids).lstlb, /DESTROY
			(*(*info).winids).lstlb = 0
			(*(*info).winswitch).showls = 0
		ENDELSE
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).lstlb,(*(*info).winids).lswid,(*(*info).winids).lsdrawid], labels=['lstlb','lswid','lsdrawid']
	ENDELSE
END

PRO CRISPEX_DISPLAYS_PHIS_RESIZE, event							
; Spectral phi slice window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_PHIS_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).phisxres, (*(*info).winsizes).phisyres, (*(*info).plotpos).phisxmargin_init, (*(*info).plotpos).phisxwall_init, $
		phisxres, phisyres, phiswidth, phisheight, phisx0, phisx1, phisy0, phisy1, (*(*info).plotswitch).v_dop_set, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).phisxres = phisxres			& 	(*(*info).winsizes).phisyres = phisyres
		(*(*info).plotpos).phisx0 = phisx0			&	(*(*info).plotpos).phisx1 = phisx1
		(*(*info).plotpos).phisy0 = phisy0			&	(*(*info).plotpos).phisy1 = phisy1
		(*(*info).plotpos).phisxplspw = phisx1 - phisx0		&	(*(*info).plotpos).phisyplspw = phisy1 - phisy0
		(*(*info).plotaxes).phisxticklen = -1 * (*(*info).plotaxes).ticklen / phisheight
		(*(*info).plotaxes).phisyticklen = -1 * (*(*info).plotaxes).ticklen / phiswidth
		(*(*info).dispparams).phisnlpreb = (*(*info).plotpos).phisxplspw * (*(*info).winsizes).phisxres 
		(*(*info).dispparams).nphireb = (*(*info).plotpos).phisyplspw * (*(*info).winsizes).phisyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).phisxres,(*(*info).winsizes).phisyres,(*(*info).plotpos).phisx0,(*(*info).plotpos).phisx1,$
		(*(*info).plotpos).phisy0,(*(*info).plotpos).phisy1], labels=['error','phisxres','phisyres','phisx0','phisx1','phisy0','phisy1']
	WIDGET_CONTROL, (*(*info).winids).phisdrawid, DRAW_XSIZE = (*(*info).winsizes).phisxres, DRAW_YSIZE = (*(*info).winsizes).phisyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DRAW_PHIS, event
END

PRO CRISPEX_DISPLAYS_PHIS_TOGGLE, event, NO_DRAW=no_draw
; Spectral phi slice window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_PHIS_TOGGLE', /IGNORE_LAST
	IF ((*(*info).winswitch).showphis EQ 0) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		IF (*(*info).winswitch).showsp THEN BEGIN
			IF ((*(*info).dataparams).nt GT 1) THEN *(*(*info).data).sspscan = (*(*(*info).data).scan)[(*(*info).dataparams).t] ELSE *(*(*info).data).sspscan = (*(*(*info).data).scan)
			*(*(*info).data).phiscan = (*(*(*info).data).sspscan)[*,*,((*(*info).dataparams).s * (*(*info).dataparams).nlp):(((*(*info).dataparams).s+1)*(*(*info).dataparams).nlp-1)] 
		ENDIF
		IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
		IF (*(*info).plotswitch).stokesfile THEN title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle ELSE title = ''
		wintitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).phiswintitle)[(*(*info).plotswitch).heightset]
		CRISPEX_WINDOW, (*(*info).winsizes).phisxres, (*(*info).winsizes).phisyres, (*(*info).winids).root, wintitle, tlb, wid, $
			(*(*info).winsizes).spxoffset, (*(*info).winswitch).showsp * (*(*info).winsizes).ydelta, DRAWID = phisdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_PHIS_RESIZE'
		PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).phiparams).nphi), /NODATA, YRANGE = [-(*(*info).phiparams).nphi/2.,(*(*info).phiparams).nphi/2.], /YSTYLE, $
			XR = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], POS = [(*(*info).plotpos).phisx0,(*(*info).plotpos).phisy0,$
			(*(*info).plotpos).phisx1,(*(*info).plotpos).phisy1], YTICKLEN = (*(*info).plotaxes).phisyticklen, XTICKLEN = (*(*info).plotaxes).phisxticklen, XTITLE = (*(*info).plottitles).spxtitle, $
			YTITLE = 'Position along slit [pixel]', XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, XRANGE = [((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_low], ((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTITLE = 'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol ELSE $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, XRANGE = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTICKNAME = REPLICATE(' ',60), XTITLE = title, COLOR = (*(*info).plotparams).plotcol
		(*(*info).winids).phistlb = tlb			&	(*(*info).winids).phiswid = wid
		(*(*info).winids).phisdrawid = phisdrawid	&	(*(*info).winids).phiswintitle = wintitle 
		WIDGET_CONTROL, (*(*info).winids).phistlb, SET_UVALUE = info
		(*(*info).ctrlsswitch).bwd_insensitive = 0	
		(*(*info).ctrlsswitch).fwd_insensitive = 0
		(*(*info).winswitch).showphis = 1
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).phistlb, /DESTROY
		(*(*info).winids).phistlb = 0
		(*(*info).winswitch).showphis = 0
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).phi_slider, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).nphi_slider, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = (*(*info).winswitch).showphis
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
		CRISPEX_UPDATE_PHIS, event
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).phistlb,(*(*info).winids).phiswid,(*(*info).winids).phisdrawid], labels=['phistlb','phiswid','phisdrawid']
END

PRO CRISPEX_DISPLAYS_REF_TOGGLE, event, NO_DRAW=no_draw
; Reference image window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REF_TOGGLE', /IGNORE_LAST
	(*(*info).winswitch).showref = event.SELECT
	IF (*(*info).winswitch).showref THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Reference image'
		CRISPEX_WINDOW, (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny, (*(*info).winids).root, title, reftlb, refwid, (*(*info).winsizes).xdelta,(*(*info).winsizes).ydelta, $
			DRAWID = refdrawid, DRAWBASE = refdrawbase
		(*(*info).winids).reftlb = reftlb		&	(*(*info).winids).refwid = refwid	&	(*(*info).winids).refdrawid = refdrawid
		(*(*info).winids).refdrawbase = refdrawbase	&	(*(*info).winids).refwintitle = title
		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
			CRISPEX_UPDATE_T, event
			CRISPEX_DRAW_REF, event
		ENDIF
		WIDGET_CONTROL, refdrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS
		WIDGET_CONTROL, reftlb, SET_UVALUE = info
		XMANAGER, 'CRISPEX', reftlb, /NO_BLOCK
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).ref_scaling_but, SET_BUTTON=0
		WIDGET_CONTROL, (*(*info).ctrlscp).xy_scaling_but, /SET_BUTTON
		(*(*info).scaling).imrefscaling = 0
		WIDGET_CONTROL, (*(*info).winids).reftlb, /DESTROY
		(*(*info).winids).reftlb = 0
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).ref_scaling_but, SENSITIVE = (*(*info).winswitch).showref
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SENSITIVE = (((*(*info).dataparams).nlp EQ (*(*info).dataparams).refnlp) AND ((*(*info).dataparams).refnlp GT 1))
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SENSITIVE = ((*(*info).dataparams).refnlp GT 1)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).reftlb,(*(*info).winids).refwid,(*(*info).winids).refdrawid], labels=['reftlb','refwid','refdrawid']
END

PRO CRISPEX_DISPLAYS_RESIZE_ERROR, event
; Opens error window on resize error
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	CRISPEX_WINDOW_OK, event,'ERROR!','Window resize request cannot be completed:','resize values beyond boundaries.','Reverted to old window size.',$
		OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
	(*(*info).winids).errtlb = tlb
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).errtlb], labels=['errtlb']
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event				
; Updates restored loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES'
	FOR i=0,N_ELEMENTS(*(*(*info).winids).restloopwid)-1 DO BEGIN
		WSET, (*(*(*info).winids).restloopwid)[i]
		PLOT, FINDGEN((*(*(*info).loopsdata).rest_loopsize)[i]), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
			(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).restloopx0,(*(*info).plotpos).restloopy0,	(*(*info).plotpos).restloopx1,(*(*info).plotpos).restloopy1], $
			YTICKLEN = (*(*info).plotaxes).restloopyticklen, XTICKLEN = (*(*info).plotaxes).restloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
			BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).winids).restloopwid)[i]], labels=['Window ID for replot']
	ENDFOR
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_RESIZE, event				
; Restored loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).restloopxres, (*(*info).winsizes).restloopyres, (*(*info).plotpos).restloopxmargin_init, (*(*info).plotpos).restloopxwall_init, $
		restloopxres, restloopyres, restloopwidth, restloopheight, restloopx0, restloopx1, restloopy0, restloopy1, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).restloopxres = restloopxres			& 	(*(*info).winsizes).restloopyres = restloopyres
		(*(*info).plotpos).restloopx0 = restloopx0			&	(*(*info).plotpos).restloopx1 = restloopx1
		(*(*info).plotpos).restloopy0 = restloopy0			&	(*(*info).plotpos).restloopy1 = restloopy1
		(*(*info).plotpos).restloopxplspw = restloopx1 - restloopx0	&	(*(*info).plotpos).restloopyplspw = restloopy1 - restloopy0
		(*(*info).plotaxes).restloopxticklen = -1 * (*(*info).plotaxes).ticklen / restloopheight
		(*(*info).plotaxes).restloopyticklen = -1 * (*(*info).plotaxes).ticklen / restloopwidth
		(*(*info).dispparams).restloopnlxreb = (*(*info).plotpos).restloopxplspw * (*(*info).winsizes).restloopxres 
		(*(*info).dispparams).restloopntreb = (*(*info).plotpos).restloopyplspw * (*(*info).winsizes).restloopyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).loopxres,(*(*info).winsizes).loopyres,(*(*info).plotpos).loopx0,(*(*info).plotpos).loopx1,$
		(*(*info).plotpos).loopy0,(*(*info).plotpos).loopy1], labels=['error','loopxres','loopyres','loopx0','loopx1','loopy0','loopy1']
	FOR i=0,N_ELEMENTS(*(*(*info).winids).restloopdrawid)-1 DO WIDGET_CONTROL, (*(*(*info).winids).restloopdrawid)[i], DRAW_XSIZE = (*(*info).winsizes).restloopxres, DRAW_YSIZE = (*(*info).winsizes).restloopyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event
	CRISPEX_DRAW_REST_LOOP, event
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_SELECT, event						
; Restored loopslab display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_SELECT'
	WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
;	print,event.INDEX,(*(*info).winswitch).showrestloop
	IF (event.INDEX GT 0) THEN BEGIN
		(*(*info).winswitch).showrestloop = 1
		sel_disp_loop = WHERE(*(*(*info).restoreparams).disp_loopnr EQ (event.INDEX-1))
		IF (TOTAL(sel_disp_loop) GE 0) THEN BEGIN			; If the selected loop is being displayed, it should be destroyed
			WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[sel_disp_loop], /DESTROY
			list_values[event.INDEX] = 'Display time slice '+STRTRIM(event.INDEX-1,2)
			WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = event.INDEX
			(*(*(*info).winids).restlooptlb)[sel_disp_loop] = 0
			wherenot0 = WHERE(*(*(*info).winids).restlooptlb NE 0, COMPLEMENT=where0)
			where0 = where0[0]
			old_restricted_t_range = TOTAL(*(*(*info).dispswitch).restricted_t_range)
			old_restricted_lp_range = TOTAL(*(*(*info).dispswitch).restricted_lp_range)
			old_restricted_lp_ref_range = TOTAL(*(*(*info).dispswitch).restricted_lp_ref_range)
			IF (TOTAL(wherenot0) NE -1) THEN BEGIN
				*(*(*info).winids).restlooptlb = (*(*(*info).winids).restlooptlb)[wherenot0]
				*(*(*info).winids).restloopwid = (*(*(*info).winids).restloopwid)[wherenot0]
				*(*(*info).winids).restloopdrawid = (*(*(*info).winids).restloopdrawid)[wherenot0]
				*(*(*info).restoreparams).disp_loopnr = (*(*(*info).restoreparams).disp_loopnr)[wherenot0]
				*(*(*info).restoreparams).disp_imref = (*(*(*info).restoreparams).disp_imref)[wherenot0]
				*(*(*info).dispswitch).restricted_t_range = (*(*(*info).dispswitch).restricted_t_range)[wherenot0]
				*(*(*info).dispswitch).restricted_lp_range = (*(*(*info).dispswitch).restricted_lp_range)[wherenot0]
				*(*(*info).dispswitch).restricted_lp_ref_range = (*(*(*info).dispswitch).restricted_lp_ref_range)[wherenot0]
				*(*(*info).restoreparams).disp_slices = (*(*(*info).restoreparams).disp_slices)[wherenot0]
				*(*(*info).restoreparams).disp_ref_slices = (*(*(*info).restoreparams).disp_ref_slices)[wherenot0]
				sel_reorder = WHERE(wherenot0 GT where0)
				IF (sel_reorder[0] NE -1) THEN BEGIN
					FOR k=0,N_ELEMENTS(sel_reorder)-1 DO BEGIN
						*(*(*(*info).loopsdata).rest_loopslice[where0+k]) = *(*(*(*info).loopsdata).rest_loopslice[wherenot0[sel_reorder[k]]])
						*(*(*(*info).loopsdata).rest_loopslab[where0+k]) = *(*(*(*info).loopsdata).rest_loopslab[wherenot0[sel_reorder[k]]])
						*(*(*(*info).loopsdata).rest_crossloc[where0+k]) = *(*(*(*info).loopsdata).rest_crossloc[wherenot0[sel_reorder[k]]])
					ENDFOR
				ENDIF ELSE k=0
				*(*(*(*info).loopsdata).rest_loopslice[where0+k]) = 0
				*(*(*(*info).loopsdata).rest_loopslab[where0+k]) = 0
				*(*(*(*info).loopsdata).rest_crossloc[where0+k]) = 0
			ENDIF ELSE BEGIN
				*(*(*info).restoreparams).disp_loopnr = -1
				*(*(*info).restoreparams).disp_imref = -1
				(*(*info).restoreparams).disp_slices = PTR_NEW(0)
				(*(*info).restoreparams).disp_ref_slices = PTR_NEW(0)
				(*(*info).winids).restlooptlb = PTR_NEW(0)
				(*(*info).dispswitch).restricted_t_range = PTR_NEW(0)
				(*(*info).dispswitch).restricted_lp_range = PTR_NEW(0)
				(*(*info).dispswitch).restricted_lp_ref_range = PTR_NEW(0)
				(*(*info).winswitch).showrestloop = 0
				*(*(*(*info).loopsdata).rest_loopslice[0]) = 0 
				*(*(*(*info).loopsdata).rest_loopslab[0]) = 0
				*(*(*(*info).loopsdata).rest_crossloc[0]) = 0
			ENDELSE
			IF ((TOTAL(*(*(*info).dispswitch).restricted_t_range) EQ 0) AND old_restricted_t_range) THEN CRISPEX_DISPRANGE_T_RESET, event
			IF ((TOTAL(*(*(*info).dispswitch).restricted_lp_range) EQ 0) AND old_restricted_lp_range) THEN CRISPEX_DISPRANGE_LP_RESET, event
			IF ((TOTAL(*(*(*info).dispswitch).restricted_lp_ref_range) EQ 0) AND old_restricted_lp_ref_range) THEN CRISPEX_DISPRANGE_LP_REF_RESET, event
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SENSITIVE = (TOTAL(*(*(*info).restoreparams).disp_slices) EQ 0) 
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SENSITIVE = (TOTAL(*(*(*info).restoreparams).disp_ref_slices) EQ 0) 
		ENDIF ELSE BEGIN
			(*(*info).restoreparams).disp_loopfile = (*(*(*info).restoreparams).cfiles)[event.INDEX-1]
			refbase = FILE_BASENAME(STRMID((*(*info).dataparams).reffilename,0,STRPOS((*(*info).dataparams).reffilename,'.',/REVERSE_SEARCH)))
			IF (STRLEN(refbase) GT 0) THEN disp_imref = STRCMP(refbase,FILE_BASENAME((*(*info).restoreparams).disp_loopfile),STRLEN(refbase)) ELSE disp_imref = 0
			IF (TOTAL(*(*(*info).restoreparams).disp_loopnr) GE 0) THEN BEGIN
				*(*(*info).restoreparams).disp_loopnr = [*(*(*info).restoreparams).disp_loopnr,event.INDEX-1] 
				*(*(*info).restoreparams).disp_imref = [*(*(*info).restoreparams).disp_imref,disp_imref]
			ENDIF ELSE BEGIN
				*(*(*info).restoreparams).disp_loopnr = event.INDEX-1
				*(*(*info).restoreparams).disp_imref = disp_imref
			ENDELSE
			list_values[event.INDEX] = 'Hide time slice '+STRTRIM(event.INDEX-1,2)
			WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = event.INDEX
			CRISPEX_DISPLAYS_RESTORE_LOOPSLAB, event
		ENDELSE
	ENDIF ELSE BEGIN
		(*(*info).restoreparams).disp_loopfile = '0'
		IF (*(*info).winswitch).showrestloop THEN BEGIN
			(*(*info).winswitch).showrestloop = 0
			CRISPEX_DISPRANGE_T_RESET, event
			CRISPEX_DISPRANGE_LP_RESET, event
			CRISPEX_DISPRANGE_LP_REF_RESET, event
			FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO BEGIN
				WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[i], /DESTROY
				list_values[(*(*(*info).restoreparams).disp_loopnr)[i]+1] = 'Display time slice '+STRTRIM((*(*(*info).restoreparams).disp_loopnr)[i],2)
				*(*(*(*info).loopsdata).rest_loopslice[i]) = 0
				*(*(*(*info).loopsdata).rest_loopslab[i]) = 0
				*(*(*(*info).loopsdata).rest_crossloc[i]) = 0
			ENDFOR
			WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = event.INDEX
			(*(*info).winids).restlooptlb = PTR_NEW(0)
			(*(*info).dispswitch).restricted_t_range = PTR_NEW(0)
			(*(*info).dispswitch).restricted_lp_range = PTR_NEW(0)
			*(*(*info).restoreparams).disp_loopnr = -1
			(*(*info).restoreparams).disp_slices = PTR_NEW(0)
			(*(*info).restoreparams).disp_ref_slices = PTR_NEW(0)
		ENDIF
		WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider,/SENSITIVE
	ENDELSE
;	IF ((*(*(*info).restoreparams).disp_loopnr)[0] NE -1) THEN print,*(*(*info).restoreparams).disp_loopnr
;	FOR k=0,(*(*info).restoreparams).cfilecount-1 DO BEGIN
;		print,'------'
;		help,*(*(*(*info).loopsdata).rest_loopslice[k]),/str
;		help,*(*(*(*info).loopsdata).rest_loopslab[k]),/str
;		help,*(*(*(*info).loopsdata).rest_crossloc[k]),/str 
;	ENDFOR
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB, event, NO_DRAW=no_draw, INDEX=index
; Restored loopslab display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB'
		update_t_range = 0
		update_lp_range = 0
		restricted_t_range = 0
		restricted_lp_range = 0
		restricted_lp_ref_range = 0
		ref_slice_only = 0
		WIDGET_CONTROL,/HOURGLASS
		RESTORE, (*(*info).restoreparams).disp_loopfile
		IF (N_ELEMENTS(INDEX) EQ 1) THEN idx = index ELSE idx = N_ELEMENTS(*(*(*info).restoreparams).disp_loopnr)-1
		*(*(*(*info).loopsdata).rest_crossloc[idx]) = vertices
		IF (N_ELEMENTS(loop_slab) GT 0) THEN loopslab = loop_slab ELSE loopslab = loop_slice
		*(*(*(*info).loopsdata).rest_loopslab[idx]) = loopslab
		slice_only = (SIZE(*(*(*(*info).loopsdata).rest_loopslab[idx]),/N_DIMENSIONS) LT 3)
		IF slice_only THEN BEGIN			; Only a slice
			*(*(*(*info).loopsdata).rest_loopslice[idx]) = *(*(*(*info).loopsdata).rest_loopslab[idx])
			IF ((SIZE(*(*(*(*info).loopsdata).rest_loopslice[idx])))[2] NE (*(*info).dataparams).nt) THEN BEGIN
				IF (N_ELEMENTS(t_low) EQ 0) THEN BEGIN
					CRISPEX_WINDOW_OK, event, 'WARNING!', $
						'The slice to be loaded has a reduced temporal range,','however the format in which it was saved does not',$
						'allow for correct slice restoration. If display','is required, please save the slice again.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
					(*(*info).winids).warntlb = tlb
					RETURN
				ENDIF ELSE BEGIN
					IF ((t_low GE (*(*info).dispparams).t_upp) OR (t_upp LT (*(*info).dispparams).t_low)) THEN BEGIN
						CRISPEX_WINDOW_OK, event, 'WARNING!', $
							'The temporal range of the loaded slice falls outside','the range set by the currently loaded slices. If',$
							'display is required, please close all currently','loaded slices before proceeding.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
						(*(*info).winids).warntlb = tlb
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
						caseidx = (*(*(*info).restoreparams).disp_loopnr)[idx]
						list_values[caseidx+1] = 'Display time slice '+STRTRIM(caseidx,2)
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = caseidx+1
						*(*(*info).restoreparams).disp_loopnr = (*(*(*info).restoreparams).disp_loopnr)[0:idx-1]
						*(*(*info).restoreparams).disp_imref = (*(*(*info).restoreparams).disp_imref)[0:idx-1]
						RETURN
					ENDIF ELSE BEGIN
						update_t_range = 1
						restricted_t_range = 1
						(*(*info).dispparams).t_low = t_low
						(*(*info).dispparams).t_upp = t_upp
						IF (N_ELEMENTS(t_saved) EQ 0) THEN (*(*info).dataparams).t = (*(*info).dispparams).t_low ELSE (*(*info).dataparams).t = t_saved
						WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = 0
					ENDELSE
				ENDELSE
			ENDIF
			IF (*(*(*info).restoreparams).disp_imref)[idx] THEN BEGIN
				restricted_lp_ref_range = 1
				(*(*info).dataparams).lp_ref = spect_pos > (*(*info).dispparams).lp_ref_low < (*(*info).dispparams).lp_ref_upp 
				WIDGET_CONTROL,(*(*info).ctrlscp).lp_ref_slider,SENSITIVE = 0, SET_VALUE = (*(*info).dataparams).lp_ref
				ref_slice_only = 1	&	slice_only = 0
			ENDIF ELSE BEGIN
				restricted_lp_range = 1
				(*(*info).dataparams).lp = spect_pos > (*(*info).dispparams).lp_low < (*(*info).dispparams).lp_upp
				WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider,SENSITIVE = 0, SET_VALUE = (*(*info).dataparams).lp
				WIDGET_CONTROL,(*(*info).ctrlscp).lower_lp_text, SENSITIVE = 0
				WIDGET_CONTROL,(*(*info).ctrlscp).upper_lp_text, SENSITIVE = 0
			ENDELSE
			CRISPEX_UPDATE_T, event
		ENDIF ELSE BEGIN										; A full or partial slab
			IF ((SIZE(*(*(*(*info).loopsdata).rest_loopslab[idx])))[2] NE (*(*info).dataparams).nt) THEN BEGIN
				IF (N_ELEMENTS(t_saved) EQ 0) THEN (*(*info).dataparams).t = (*(*info).dispparams).t_low ELSE (*(*info).dataparams).t = t_saved
				IF (N_ELEMENTS(t_low) EQ 0) THEN BEGIN
					CRISPEX_WINDOW_OK, event, 'WARNING!', $
						'The slice to be loaded has a reduced temporal range,','however the format in which it was saved does not',$
						'allow for correct slice restoration. If display','is required, please save the slice again.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
					(*(*info).winids).warntlb = tlb
					RETURN
				ENDIF ELSE BEGIN
					IF ((t_low GE (*(*info).dispparams).t_upp) OR (t_upp LT (*(*info).dispparams).t_low)) THEN BEGIN
						CRISPEX_WINDOW_OK, event, 'WARNING!', $
							'The temporal range of the loaded slice falls outside','the range set by the currently loaded slices. If',$
							'display is required, please close all currently','loaded slices before proceeding.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
						(*(*info).winids).warntlb = tlb
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
						caseidx = (*(*(*info).restoreparams).disp_loopnr)[idx]
						list_values[caseidx+1] = 'Display time slice '+STRTRIM(caseidx,2)
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = caseidx+1
						*(*(*info).restoreparams).disp_loopnr = (*(*(*info).restoreparams).disp_loopnr)[0:idx-1]
						*(*(*info).restoreparams).disp_imref = (*(*(*info).restoreparams).disp_imref)[0:idx-1]
						RETURN
					ENDIF ELSE BEGIN
						update_t_range = 1
						restricted_t_range = 1
						(*(*info).dispparams).t_low = t_low
						(*(*info).dispparams).t_upp = t_upp
						IF (N_ELEMENTS(t) EQ 0) THEN (*(*info).dataparams).t = (*(*info).dispparams).t_low ELSE (*(*info).dataparams).t = t
						WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = 0
					ENDELSE
				ENDELSE
			ENDIF
			IF (*(*(*info).restoreparams).disp_imref)[idx] THEN nlp_comp = (*(*info).dataparams).refnlp ELSE nlp_comp = (*(*info).dataparams).nlp
			IF ((SIZE(*(*(*(*info).loopsdata).rest_loopslab[idx])))[3] NE nlp_comp) THEN BEGIN
				update_lp_range = 1
				IF (*(*(*info).restoreparams).disp_imref)[idx] THEN BEGIN	; Reference
					restricted_lp_ref_range = 1
					(*(*info).dispparams).lp_ref_low = spect_pos_low > (*(*info).dispparams).lp_ref_low
					(*(*info).dispparams).lp_ref_upp = spect_pos_upp < (*(*info).dispparams).lp_ref_upp
					(*(*info).dispparams).lp_ref_range = (*(*info).dispparams).lp_ref_upp - (*(*info).dispparams).lp_ref_low + 1
					(*(*info).dataparams).lp_ref = spect_pos > (*(*info).dispparams).lp_ref_low < (*(*info).dispparams).lp_ref_upp
					WIDGET_CONTROL,(*(*info).ctrlscp).lp_ref_slider, SET_SLIDER_MIN = (*(*info).dispparams).lp_ref_low, SET_SLIDER_MAX = (*(*info).dispparams).lp_ref_upp, $
						SET_VALUE = (*(*info).dataparams).lp_ref
				ENDIF ELSE BEGIN						; Main
					restricted_lp_range = 1
					(*(*info).dispparams).lp_low = spect_pos_low > (*(*info).dispparams).lp_low
					(*(*info).dispparams).lp_upp = spect_pos_upp < (*(*info).dispparams).lp_upp
					(*(*info).dataparams).lp = spect_pos > (*(*info).dispparams).lp_low < (*(*info).dispparams).lp_upp
					WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2), SENSITIVE = 0
					WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2), SENSITIVE = 0
					WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider,SET_SLIDER_MIN = (*(*info).dispparams).lp_low, SET_SLIDER_MAX = (*(*info).dispparams).lp_upp, SET_VALUE = (*(*info).dataparams).lp
				ENDELSE
			ENDIF
			CRISPEX_UPDATE_LP, event
		ENDELSE
		IF (loop_size GT 0) THEN (*(*(*info).loopsdata).rest_loopsize)[idx] = loop_size ELSE (*(*(*info).loopsdata).rest_loopsize)[idx] = (SIZE(loopslab))[1]
		wintitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': T-slice along loop '
		CRISPEX_WINDOW, (*(*info).winsizes).restloopxres, (*(*info).winsizes).restloopyres, (*(*info).winids).root, wintitle+STRTRIM((*(*(*info).restoreparams).disp_loopnr)[idx],2), tlb, wid, $
			(*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta,((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = disp_rest_loopdrawid, RESIZING = 1, $
			RES_HANDLER = 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_RESIZE'
		PLOT, FINDGEN((*(*(*info).loopsdata).rest_loopsize)[idx]), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
			(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).restloopx0,(*(*info).plotpos).restloopy0,(*(*info).plotpos).restloopx1,(*(*info).plotpos).restloopy1], $
			YTICKLEN = (*(*info).plotaxes).restloopyticklen, XTICKLEN = (*(*info).plotaxes).restloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop',$
			BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		IF ((*(*(*info).winids).restlooptlb)[0] NE 0) THEN BEGIN
			*(*(*info).winids).restlooptlb = [*(*(*info).winids).restlooptlb,tlb]
			*(*(*info).winids).restloopwid = [*(*(*info).winids).restloopwid,wid]
			*(*(*info).winids).restloopdrawid = [*(*(*info).winids).restloopdrawid,disp_rest_loopdrawid]
			*(*(*info).winids).restloopwintitle = [*(*(*info).winids).restloopwintitle,wintitle]
			*(*(*info).dispswitch).restricted_t_range = [*(*(*info).dispswitch).restricted_t_range,restricted_t_range]
			*(*(*info).dispswitch).restricted_lp_range = [*(*(*info).dispswitch).restricted_lp_range,restricted_lp_range]
			*(*(*info).dispswitch).restricted_lp_ref_range = [*(*(*info).dispswitch).restricted_lp_ref_range,restricted_lp_ref_range]
			*(*(*info).restoreparams).disp_slices = [*(*(*info).restoreparams).disp_slices,slice_only]
			*(*(*info).restoreparams).disp_ref_slices = [*(*(*info).restoreparams).disp_ref_slices,ref_slice_only]
		ENDIF ELSE BEGIN
			*(*(*info).winids).restlooptlb = tlb
			*(*(*info).winids).restloopwid = wid
			*(*(*info).winids).restloopdrawid = disp_rest_loopdrawid
			*(*(*info).winids).restloopwintitle = wintitle
			*(*(*info).dispswitch).restricted_t_range = restricted_t_range
			*(*(*info).dispswitch).restricted_lp_range = restricted_lp_range
			*(*(*info).dispswitch).restricted_lp_ref_range = restricted_lp_ref_range
			*(*(*info).restoreparams).disp_slices = slice_only
			*(*(*info).restoreparams).disp_ref_slices = ref_slice_only
		ENDELSE
		WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[N_ELEMENTS(*(*(*info).winids).restlooptlb)-1], SET_UVALUE = info
		IF update_t_range THEN BEGIN
			CRISPEX_DISPRANGE_T_RANGE, event
			WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
		ENDIF
		IF update_lp_range THEN BEGIN
			CRISPEX_DISPRANGE_LP_RANGE, event
			WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, SENSITIVE = 0
		ENDIF
		IF ((update_t_range EQ 0) AND (update_lp_range EQ 0) AND ~KEYWORD_SET(NO_DRAW)) THEN CRISPEX_DRAW, event
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [*(*(*info).winids).restlooptlb,*(*(*info).winids).restloopwid,*(*(*info).winids).restloopdrawid], $
			labels=['restlooptlb','restloopwid','restloopdrawid']
END

PRO CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
; Updates retrieved detection loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES'
	WSET, (*(*info).winids).retrdetwid
	PLOT, FINDGEN((*(*info).loopsdata).det_loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).retrdetx0,(*(*info).plotpos).retrdety0,(*(*info).plotpos).retrdetx1,(*(*info).plotpos).retrdety1], $
		YTICKLEN = (*(*info).plotaxes).retrdetyticklen, XTICKLEN = (*(*info).plotaxes).retrdetxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).retrdetwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_RESIZE, event	
; Retrieved detection loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).retrdetxres, (*(*info).winsizes).retrdetyres, (*(*info).plotpos).retrdetxmargin_init, (*(*info).plotpos).retrdetxwall_init, $
		retrdetxres, retrdetyres, retrdetwidth, retrdetheight, retrdetx0, retrdetx1, retrdety0, retrdety1, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).retrdetxres = retrdetxres			& 	(*(*info).winsizes).retrdetyres = retrdetyres
		(*(*info).plotpos).retrdetx0 = retrdetx0			&	(*(*info).plotpos).retrdetx1 = retrdetx1
		(*(*info).plotpos).retrdety0 = retrdety0			&	(*(*info).plotpos).retrdety1 = retrdety1
		(*(*info).plotpos).retrdetxplspw = retrdetx1 - retrdetx0	&	(*(*info).plotpos).retrdetyplspw = retrdety1 - retrdety0
		(*(*info).plotaxes).retrdetxticklen = -1 * (*(*info).plotaxes).ticklen / retrdetheight
		(*(*info).plotaxes).retrdetyticklen = -1 * (*(*info).plotaxes).ticklen / retrdetwidth
		(*(*info).dispparams).retrdetnlxreb = (*(*info).plotpos).retrdetxplspw * (*(*info).winsizes).retrdetxres 
		(*(*info).dispparams).retrdetntreb = (*(*info).plotpos).retrdetyplspw * (*(*info).winsizes).retrdetyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).loopxres,(*(*info).winsizes).loopyres,(*(*info).plotpos).loopx0,(*(*info).plotpos).loopx1,$
		(*(*info).plotpos).loopy0,(*(*info).plotpos).loopy1], labels=['error','loopxres','loopyres','loopx0','loopx1','loopy0','loopy1']
	WIDGET_CONTROL, (*(*info).winids).retrdetdrawid, DRAW_XSIZE = (*(*info).winsizes).retrdetxres, DRAW_YSIZE = (*(*info).winsizes).retrdetyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
	CRISPEX_DRAW_RETR_DET, event
END

PRO CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB, event, NO_DRAW=no_draw
; Retrieved detection loopslab display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB'
	IF ~KEYWORD_SET(NO_DRAW) THEN (*(*info).detparams).idx = event.INDEX-1
	IF ((*(*info).detparams).idx GE 0) THEN BEGIN
		IF ((*(*info).winids).retrdettlb GT 0) THEN WIDGET_CONTROL, (*(*info).winids).retrdettlb, /DESTROY
		(*(*info).winswitch).showretrdet = 1
		WIDGET_CONTROL,/HOURGLASS
		CRISPEX_RETRIEVE_DET_GET_SLICE, event
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_UPDATE_LP, event
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': T-slice along detection '
		CRISPEX_WINDOW, (*(*info).winsizes).retrdetxres, (*(*info).winsizes).retrdetyres, (*(*info).winids).root, title+STRTRIM((*(*info).detparams).idx,2), tlb, wid, $
			(*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta,((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = disp_retr_detdrawid, RESIZING = 1, $
			RES_HANDLER = 'CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_RESIZE'
		PLOT, FINDGEN((*(*info).loopsdata).det_loopsize), FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt, $
			(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).retrdetx0,(*(*info).plotpos).retrdety0,(*(*info).plotpos).retrdetx1,(*(*info).plotpos).retrdety1], $
			YTICKLEN = (*(*info).plotaxes).retrdetyticklen, XTICKLEN = (*(*info).plotaxes).retrdetxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
			BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		(*(*info).winids).retrdettlb = tlb	&	(*(*info).winids).retrdetwid = wid	&	(*(*info).winids).retrdetdrawid = disp_retr_detdrawid
		(*(*info).winids).retrdetwintitle = title
		WIDGET_CONTROL, (*(*info).winids).retrdettlb, SET_UVALUE = info
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DISPRANGE_T_RANGE, event
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = 0 
	ENDIF ELSE BEGIN
		IF (*(*info).winswitch).showretrdet THEN BEGIN
			(*(*info).winswitch).showretrdet = 0
			WIDGET_CONTROL, (*(*info).winids).retrdettlb, /DESTROY
			(*(*info).winids).retrdettlb = 0
		ENDIF
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).retrdettlb,(*(*info).winids).retrdetwid,(*(*info).winids).retrdetdrawid], $
		labels=['retrdettlb','retrdetwid','retrdetdrawid']
END
	
PRO CRISPEX_DISPLAYS_REFLS_RESIZE, event							
; Detailed spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFLS_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).reflsxres, (*(*info).winsizes).reflsyres, (*(*info).plotpos).reflsxmargin_init, (*(*info).plotpos).reflsxwall_init, $
		reflsxres, reflsyres, reflswidth, reflsheight, reflsx0, reflsx1, reflsy0, reflsy1, (*(*info).plotswitch).v_dop_set_ref, ERROR=error, /GOLDEN
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).reflsxres = reflsxres	& 	(*(*info).winsizes).reflsyres = reflsyres
		(*(*info).plotpos).reflsx0 = reflsx0		&	(*(*info).plotpos).reflsx1 = reflsx1
		(*(*info).plotpos).reflsy0 = reflsy0		&	(*(*info).plotpos).reflsy1 = reflsy1
		(*(*info).plotaxes).reflsxticklen = (*(*info).plotaxes).ticklen / reflsheight
		(*(*info).plotaxes).reflsyticklen = (*(*info).plotaxes).ticklen / reflswidth
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).reflsxres,(*(*info).winsizes).reflsyres,(*(*info).plotpos).reflsx0,(*(*info).plotpos).reflsx1,$
		(*(*info).plotpos).reflsy0,(*(*info).plotpos).reflsy1], labels=['error','reflsxres','reflsyres','reflsx0','reflsx1','reflsy0','reflsy1']
	WIDGET_CONTROL, (*(*info).winids).reflsdrawid, DRAW_XSIZE = (*(*info).winsizes).reflsxres, DRAW_YSIZE = (*(*info).winsizes).reflsyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DRAW_REFLS, event
END

PRO CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
; Updates reference temporal spectrum display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFSP_REPLOT_AXES'
	WSET, (*(*info).winids).refspwid
	PLOT, (*(*info).dataparams).reflps, FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).refspx0,(*(*info).plotpos).refspy0,(*(*info).plotpos).refspx1,(*(*info).plotpos).refspy1], $
		YTICKLEN = (*(*info).plotaxes).refspxticklen, XTICKLEN = (*(*info).plotaxes).refspyticklen, XTITLE = (*(*info).plottitles).refspxtitle, YTITLE = (*(*info).plottitles).spytitle, $
		XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
		XR = [((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_low],((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_upp]]
	IF ((*(*info).plotswitch).v_dop_set_ref EQ 1) THEN AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).refspxticklen, XRANGE = [((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_low], $
		((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_upp]], XSTYLE=1, XTITLE = 'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refspwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_REFSP_RESIZE, event
; Reference temporal spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFSP_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).refspxres, (*(*info).winsizes).refspyres, (*(*info).plotpos).refspxmargin_init, (*(*info).plotpos).refspxwall_init, $
		refspxres, refspyres, refspwidth, refspheight, refspx0, refspx1, refspy0, refspy1, (*(*info).plotswitch).v_dop_set_ref, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).refspxres = refspxres		& 	(*(*info).winsizes).refspyres = refspyres
		(*(*info).plotpos).refspx0 = refspx0			&	(*(*info).plotpos).refspx1 = refspx1
		(*(*info).plotpos).refspy0 = refspy0			&	(*(*info).plotpos).refspy1 = refspy1
		(*(*info).plotpos).refxplspw = refspx1 - refspx0	&	(*(*info).plotpos).refyplspw = refspy1 - refspy0
		(*(*info).plotaxes).refspxticklen = -1 * (*(*info).plotaxes).ticklen / refspheight
		(*(*info).plotaxes).refspyticklen = -1 * (*(*info).plotaxes).ticklen / refspwidth
		(*(*info).dispparams).refnlpreb = (*(*info).plotpos).refxplspw * (*(*info).winsizes).refspxres 
		(*(*info).dispparams).refntreb = (*(*info).plotpos).refyplspw * (*(*info).winsizes).refspyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error, (*(*info).winsizes).refspxres,(*(*info).winsizes).refspyres,(*(*info).plotpos).refspx0,(*(*info).plotpos).refspx1,$
		(*(*info).plotpos).refspy0,(*(*info).plotpos).refspy1], labels=['error','refspxres','refspyres','refspx0','refspx1','refspy0','refspy1']
	WIDGET_CONTROL, (*(*info).winids).refspdrawid, DRAW_XSIZE = (*(*info).winsizes).refspxres, DRAW_YSIZE = (*(*info).winsizes).refspyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
	CRISPEX_DRAW_REFSP, event
END

PRO CRISPEX_DISPLAYS_REFSP_TOGGLE, event, NO_DRAW=no_draw
; Reference temporal spectrum display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_REFSP_TOGGLE', /IGNORE_LAST
	IF ((*(*info).winswitch).showrefsp EQ 0) THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).refspwintitle)[(*(*info).plotswitch).refheightset]
		CRISPEX_WINDOW, (*(*info).winsizes).refspxres, (*(*info).winsizes).refspyres, (*(*info).winids).root, title, refsptlb, refspwid, $
			(*(*info).winsizes).spxoffset, (*(*info).winswitch).showsp * (*(*info).winsizes).ydelta, DRAWID = refspdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFSP_RESIZE'
		PLOT, (*(*info).dataparams).reflps, FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
			(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt],/YS, POS = [(*(*info).plotpos).refspx0,(*(*info).plotpos).refspy0,(*(*info).plotpos).refspx1,(*(*info).plotpos).refspy1], $
			YTICKLEN = (*(*info).plotaxes).refspyticklen, XTICKLEN = (*(*info).plotaxes).refspxticklen, XTITLE = (*(*info).plottitles).refspxtitle, YTITLE = (*(*info).plottitles).spytitle, $
			XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
			XR = [(*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low], (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_upp]]
		IF ((*(*info).plotswitch).v_dop_set_ref EQ 1) THEN $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).refspxticklen, XRANGE = [((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_low], $
				((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_upp]], XSTYLE=1,XTITLE = 'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol
		(*(*info).winids).refsptlb = refsptlb		&	(*(*info).winids).refspwid = refspwid	&	(*(*info).winswitch).showrefsp = 1
		(*(*info).winids).refspdrawid = refspdrawid	&	(*(*info).winids).refspwintitle = title 
		WIDGET_CONTROL, (*(*info).winids).refsptlb, SET_UVALUE = info
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_REFSP, event
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).refsptlb, /DESTROY
		(*(*info).winids).refsptlb = 0
		(*(*info).winswitch).showrefsp = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refsptlb,(*(*info).winids).refspwid,(*(*info).winids).refspdrawid], labels=['refsptlb','refspwid','refspdrawid']
END

PRO CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
; Updates temporal spectrum display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_SP_REPLOT_AXES'
	WSET, (*(*info).winids).spwid
	IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
	IF (*(*info).plotswitch).stokesfile THEN title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle ELSE title = ''
	PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
		(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt], /YS, POS = [(*(*info).plotpos).spx0,(*(*info).plotpos).spy0,(*(*info).plotpos).spx1,(*(*info).plotpos).spy1], $
		YTICKLEN = (*(*info).plotaxes).spyticklen, XTICKLEN = (*(*info).plotaxes).spxticklen, XTITLE = (*(*info).plottitles).spxtitle, YTITLE = (*(*info).plottitles).spytitle, $
		XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, XR = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low],$
		(*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]]
	IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN $
		AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, XRANGE = [((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_low], ((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
			XTITLE = title+'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol ELSE $
		AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, XRANGE = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
			XTICKNAME = REPLICATE(' ',60), XTITLE = title, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).spwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_SP_RESIZE, event
; Temporal spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_SP_RESIZE'
	CRISPEX_DISPLAYS_PLOT_RESIZE, event, event.X, event.Y, (*(*info).winsizes).spxres, (*(*info).winsizes).spyres, (*(*info).plotpos).spxmargin_init, (*(*info).plotpos).spxwall_init, $
		spxres, spyres, spwidth, spheight, spx0, spx1, spy0, spy1, (*(*info).plotswitch).v_dop_set, ERROR=error, /ACTUAL_RESIZE
	IF error THEN CRISPEX_DISPLAYS_RESIZE_ERROR, event ELSE BEGIN
		(*(*info).winsizes).spxres = spxres		& 	(*(*info).winsizes).spyres = spyres
		(*(*info).plotpos).spx0 = spx0			&	(*(*info).plotpos).spx1 = spx1
		(*(*info).plotpos).spy0 = spy0			&	(*(*info).plotpos).spy1 = spy1
		(*(*info).plotpos).xplspw = spx1 - spx0		&	(*(*info).plotpos).yplspw = spy1 - spy0
		(*(*info).plotaxes).spxticklen = -1 * (*(*info).plotaxes).ticklen / spheight
		(*(*info).plotaxes).spyticklen = -1 * (*(*info).plotaxes).ticklen / spwidth
		(*(*info).dispparams).nlpreb = (*(*info).plotpos).xplspw * (*(*info).winsizes).spxres 
		(*(*info).dispparams).ntreb = (*(*info).plotpos).yplspw * (*(*info).winsizes).spyres 
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [error,(*(*info).winsizes).spxres,(*(*info).winsizes).spyres,(*(*info).plotpos).spx0,(*(*info).plotpos).spx1,$
		(*(*info).plotpos).spy0,(*(*info).plotpos).spy1], labels=['error','spxres','spyres','spx0','spx1','spy0','spy1']
	WIDGET_CONTROL, (*(*info).winids).spdrawid, DRAW_XSIZE = (*(*info).winsizes).spxres, DRAW_YSIZE = (*(*info).winsizes).spyres
	WIDGET_CONTROL, event.TOP, SET_UVALUE = info
	CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
	CRISPEX_DRAW_SP, event
END

PRO CRISPEX_DISPLAYS_SP_TOGGLE, event, NO_DRAW=no_draw
; Temporal spectrum display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_SP_TOGGLE', /IGNORE_LAST
	IF ((*(*info).winswitch).showsp EQ 0) THEN BEGIN
		IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
		IF (*(*info).plotswitch).stokesfile THEN title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle ELSE title = ''
		wintitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).spwintitle)[(*(*info).plotswitch).heightset]
		CRISPEX_WINDOW, (*(*info).winsizes).spxres, (*(*info).winsizes).spyres, (*(*info).winids).root, wintitle, tlb, wid, $
			(*(*info).winsizes).spxoffset, 0, DRAWID = spdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_SP_RESIZE'
		PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).dispparams).t_range * (*(*info).plotaxes).dt), /NODATA, YR = [(*(*info).dispparams).t_low * (*(*info).plotaxes).dt,$
			(*(*info).dispparams).t_upp * (*(*info).plotaxes).dt],/YS, POS = [(*(*info).plotpos).spx0,(*(*info).plotpos).spy0,(*(*info).plotpos).spx1,(*(*info).plotpos).spy1], $
			YTICKLEN = (*(*info).plotaxes).spyticklen, XTICKLEN = (*(*info).plotaxes).spxticklen, XTITLE = (*(*info).plottitles).spxtitle, YTITLE = (*(*info).plottitles).spytitle, $
			XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
			XR = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]]
		IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, XRANGE = [((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_low], ((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTITLE = 'Doppler velocity [km/s]',COLOR = (*(*info).plotparams).plotcol ELSE $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, XRANGE = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTICKNAME = REPLICATE(' ',60), XTITLE = title, COLOR = (*(*info).plotparams).plotcol
		(*(*info).winids).sptlb = tlb		&	(*(*info).winids).spwid = wid	&	(*(*info).winswitch).showsp = 1
		(*(*info).winids).spdrawid = spdrawid	&	(*(*info).winids).spwintitle = wintitle
		WIDGET_CONTROL, (*(*info).winids).sptlb, SET_UVALUE = info
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_SP, event
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).sptlb, /DESTROY
		(*(*info).winids).sptlb = 0
		(*(*info).winswitch).showsp = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).sptlb,(*(*info).winids).spwid,(*(*info).winids).spdrawid], labels=['sptlb','spwid','spdrawid']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_I, event
; Stokes I for main image and temporal spectrum selection 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_I'
	(*(*info).dataparams).s = WHERE((*(*info).stokesparams).labels EQ 'I')
	CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT, event
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_Q, event
; Stokes Q for main image and temporal spectrum selection 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_Q'
	(*(*info).dataparams).s = WHERE((*(*info).stokesparams).labels EQ 'Q')
	CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT, event
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_U, event
; Stokes U for main image and temporal spectrum selection 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_U'
	(*(*info).dataparams).s = WHERE((*(*info).stokesparams).labels EQ 'U')
	CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT, event
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_V, event
; Stokes V for main image and temporal spectrum selection 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_V'
	(*(*info).dataparams).s = WHERE((*(*info).stokesparams).labels EQ 'V')
	CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT, event
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT, event
; Processing of Stokes component for main image and temporal spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_CONT'
	WIDGET_CONTROL, (*(*info).ctrlsparam).stokes_val, SET_VALUE = STRTRIM(((*(*info).stokesparams).labels)[(*(*info).dataparams).s],2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [((*(*info).stokesparams).labels)[(*(*info).dataparams).s]], labels=['Stokes image selected']
	CRISPEX_DISPLAYS_STOKES_SELECT_XY_RECOVER_YRANGE, event
	CRISPEX_UPDATE_T, event
	CRISPEX_UPDATE_SLICES, event
	IF ((*(*(*info).scaling).imagescale)[0] EQ 2) THEN CRISPEX_SCALING_MAN_FIRST, event
	IF ((*(*(*info).scaling).imagescale)[0] EQ 3) THEN CRISPEX_SCALING_MAN_CURR, event
	IF ((*(*info).winids).sptlb GT 0) THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_RECOVER_YRANGE, event
; Restores lower/upper y-values of specific Stokes component detailed spectrum plot for input
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_RECOVER_YRANGE'
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],2)
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s]], $
		labels=['Lower detspect y-value','Upper detspect y-value']
	IF (*(*info).winswitch).showint THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsint).lower_y_int_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],2)
		WIDGET_CONTROL, (*(*info).ctrlsint).upper_y_int_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s],2)
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]], $
			labels=['Lower int y-value','Upper int y-value']
	ENDIF
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_SP_I, event
; Stokes I for detailed spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_I'
	(*(*info).stokesparams).prev_select_sp = (*(*info).stokesparams).select_sp
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN BEGIN
		(*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'I')] = event.SELECT
		CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 1
	ENDIF ELSE (*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'I')] = event.SELECT
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 0
	CRISPEX_DISPLAYS_LS_RESIZE, event, /STOKES_SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRJOIN(((*(*info).stokesparams).labels)[WHERE((*(*info).stokesparams).select_sp EQ 1)],', ')],labels=['Stokes detspect selected']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_SP_Q, event
; Stokes Q for detailed spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_Q'
	(*(*info).stokesparams).prev_select_sp = (*(*info).stokesparams).select_sp
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN BEGIN
		(*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'Q')] = event.SELECT
		CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 1
	ENDIF ELSE (*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'Q')] = event.SELECT
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 0
	CRISPEX_DISPLAYS_LS_RESIZE, event, /STOKES_SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRJOIN(((*(*info).stokesparams).labels)[WHERE((*(*info).stokesparams).select_sp EQ 1)],', ')],labels=['Stokes detspect selected']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_SP_U, event
; Stokes U for detailed spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_U'
	(*(*info).stokesparams).prev_select_sp = (*(*info).stokesparams).select_sp
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN BEGIN
		(*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'U')] = event.SELECT
		CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 1
	ENDIF ELSE (*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'U')] = event.SELECT
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 0 
	CRISPEX_DISPLAYS_LS_RESIZE, event, /STOKES_SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRJOIN(((*(*info).stokesparams).labels)[WHERE((*(*info).stokesparams).select_sp EQ 1)],', ')],labels=['Stokes detspect selected']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_SP_V, event
; Stokes V for detailed spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_V'
	(*(*info).stokesparams).prev_select_sp = (*(*info).stokesparams).select_sp
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN BEGIN
		(*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'V')] = event.SELECT
		CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 1
	ENDIF ELSE (*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ 'V')] = event.SELECT
	IF (TOTAL((*(*info).stokesparams).select_sp) EQ 1) THEN CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, 0 
	CRISPEX_DISPLAYS_LS_RESIZE, event, /STOKES_SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRJOIN(((*(*info).stokesparams).labels)[WHERE((*(*info).stokesparams).select_sp EQ 1)],', ')],labels=['Stokes detspect selected']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON, event, sens
; Processing of button sensitivity of Stokes component for detailed spectrum selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_SENSBUTTON'
	FOR i=0,TOTAL((*(*info).stokesparams).select_sp)-1 DO BEGIN
		IF (((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))[i]] EQ 'I') THEN WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_i_but, SENSITIVE = sens, /SET_BUTTON 
		IF (((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))[i]] EQ 'Q') THEN WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_q_but, SENSITIVE = sens, /SET_BUTTON 
		IF (((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))[i]] EQ 'U') THEN WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_u_but, SENSITIVE = sens, /SET_BUTTON 
		IF (((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))[i]] EQ 'V') THEN WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_v_but, SENSITIVE = sens, /SET_BUTTON
	ENDFOR
END

;================================================================================= DISPLAY RANGE PROCEDURES
PRO CRISPEX_DISPRANGE_INT_LOW, event
; Handles change in lower y-value of intensity versus time display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_LOW'
	WIDGET_CONTROL, (*(*info).ctrlsint).lower_y_int_text, GET_VALUE = textvalue
	(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] = FLOAT(textvalue[0])
	IF ((*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] GE (*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]) THEN $
		(*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] + 1 
	WIDGET_CONTROL, (*(*info).ctrlsint).upper_y_int_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s],2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]], $
		labels=['Lower int y-value','Upper int y-value']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPRANGE_INT_UPP, event
; Handles change in upper y-value of intensity versus time display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_UPP'
	WIDGET_CONTROL, (*(*info).ctrlsint).upper_y_int_text, GET_VALUE = textvalue
	(*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s] = FLOAT(textvalue[0]) 
	IF ((*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s] LE (*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s]) THEN $
		(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s] - 1 
	WIDGET_CONTROL, (*(*info).ctrlsint).lower_y_int_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]], $
		labels=['Lower int y-value','Upper int y-value']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPRANGE_INT_T_LOW, event
; Handles change in lower t-value of intensity versus time display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_T_LOW'
	WIDGET_CONTROL, (*(*info).ctrlsint).lower_t_int_text, GET_VALUE = textvalue
	(*(*info).plotaxes).int_low_t = FLOAT(textvalue[0])
	IF ((*(*info).plotaxes).int_low_t GE (*(*info).plotaxes).int_upp_t) THEN (*(*info).plotaxes).int_low_t = (*(*info).plotaxes).int_upp_t - 1
	IF ((*(*info).plotaxes).int_low_t LT (*(*info).dispparams).t_first) THEN (*(*info).plotaxes).int_low_t = (*(*info).dispparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlsint).lower_t_int_text, SET_VALUE = STRTRIM((*(*info).plotaxes).int_low_t,2)
	CRISPEX_DISPRANGE_INT_T_RANGE, event
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPRANGE_INT_T_UPP, event
; Handles change in upper t-value of intensity versus time display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_T_UPP'
	WIDGET_CONTROL, (*(*info).ctrlsint).upper_t_int_text, GET_VALUE = textvalue
	(*(*info).plotaxes).int_upp_t = FLOAT(textvalue[0]) 
	IF ((*(*info).plotaxes).int_upp_t LE (*(*info).plotaxes).int_low_t) THEN (*(*info).plotaxes).int_upp_t = (*(*info).plotaxes).int_low_t + 1
	IF ((*(*info).plotaxes).int_upp_t GT (*(*info).dispparams).t_last) THEN (*(*info).plotaxes).int_upp_t = (*(*info).dispparams).t_last
	WIDGET_CONTROL, (*(*info).ctrlsint).upper_t_int_text, SET_VALUE = STRTRIM((*(*info).plotaxes).int_upp_t,2)
	CRISPEX_DISPRANGE_INT_T_RANGE, event
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPRANGE_INT_LOCK_T, event
; Locks main temporal range to intensity versus time display temporal range
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_LOCK_T'
	(*(*info).intparams).lock_t = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).intparams).lock_t], labels=['Lock main to int temporal range']
	CRISPEX_DISPRANGE_INT_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_INT_T_RANGE, event
; Locks main temporal range to intensity versus time display temporal range
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_T_RANGE'
	IF (*(*info).intparams).lock_t THEN BEGIN
		(*(*info).dispparams).t_low = (*(*info).plotaxes).int_low_t
		(*(*info).dispparams).t_upp = (*(*info).plotaxes).int_upp_t
	ENDIF ELSE BEGIN
		(*(*info).dispparams).t_low = (*(*info).dispparams).t_first
		(*(*info).dispparams).t_upp = (*(*info).dispparams).t_last
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).plotaxes).int_low_t,(*(*info).plotaxes).int_upp_t,(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp], $
		labels=['Lower int t-value','Upper int t-value','Lower t-value','Upper t-value']
	CRISPEX_DISPRANGE_T_RANGE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	IF (*(*info).winswitch).showint THEN WIDGET_CONTROL, (*(*info).ctrlsint).reset_trange_but, SENSITIVE = (((*(*info).plotaxes).int_upp_t-(*(*info).plotaxes).int_low_t+1) NE (*(*info).dataparams).nt)
END

PRO CRISPEX_DISPRANGE_INT_T_RESET, event
; Handles reset of temporal boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_INT_T_RESET'
	(*(*info).plotaxes).int_upp_t = (*(*info).dispparams).t_last
	(*(*info).plotaxes).int_low_t = (*(*info).dispparams).t_first
	CRISPEX_DISPRANGE_INT_T_RANGE, event
	WIDGET_CONTROL, (*(*info).ctrlsint).upper_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	WIDGET_CONTROL, (*(*info).ctrlsint).lower_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	WIDGET_CONTROL, (*(*info).ctrlsint).reset_trange_but, SENSITIVE = 0
END

PRO CRISPEX_DISPRANGE_LS_LOW, event
; Handles change in lower y-value of detailed spectrum display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_LOW'
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, GET_VALUE = textvalue
	IF (*(*info).ctrlsswitch).imrefdetspect THEN BEGIN
		(*(*info).plotaxes).ls_low_y_ref = FLOAT(textvalue[0]) 
		IF ((*(*info).plotaxes).ls_low_y_ref GE (*(*info).plotaxes).ls_upp_y_ref) THEN (*(*info).plotaxes).ls_upp_y_ref = (*(*info).plotaxes).ls_low_y_ref + 1 
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_upp_y_ref,2)
	ENDIF ELSE BEGIN
		(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] = FLOAT(textvalue[0]) 
		IF ((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] GE (*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s]) THEN $
			(*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] + 1 
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],2)
	ENDELSE
	CRISPEX_DISPRANGE_LS_RANGE, event
END

PRO CRISPEX_DISPRANGE_LS_UPP, event
; Handles change in upper y-value of detailed spectrum display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_UPP'
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, GET_VALUE = textvalue
	IF (*(*info).ctrlsswitch).imrefdetspect THEN BEGIN
		(*(*info).plotaxes).ls_upp_y_ref = FLOAT(textvalue[0])
		IF ((*(*info).plotaxes).ls_upp_y_ref LE (*(*info).plotaxes).ls_low_y_ref) THEN (*(*info).plotaxes).ls_low_y_ref = (*(*info).plotaxes).ls_upp_y_ref - 1 
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_low_y_ref,2)
	ENDIF ELSE BEGIN
		(*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] = FLOAT(textvalue[0])
		IF ((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] LE (*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s]) THEN $
			(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] - 1 
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],2)
	ENDELSE
	CRISPEX_DISPRANGE_LS_RANGE, event
END

PRO CRISPEX_DISPRANGE_LS_RANGE, event
; Determines range from change in lower or upper y-value of detailed spectrum display window and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_RANGE'
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).plotaxes).ls_yrange_ref = (*(*info).plotaxes).ls_upp_y_ref - (*(*info).plotaxes).ls_low_y_ref ELSE $
		(*(*(*info).plotaxes).ls_yrange)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] - (*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],$
		(*(*info).plotaxes).ls_low_y_ref,(*(*info).plotaxes).ls_upp_y_ref], labels=['Lower main y-value','Upper main y-value','Lower ref y-value','Upper ref y-value']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DISPRANGE_LS_SCALE_SELECT, event
; Handles the selection of scaling (or not) of the detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_SCALE'
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).dispswitch).ref_detspect_scale = event.SELECT ELSE (*(*info).dispswitch).detspect_scale = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dispswitch).detspect_scale,(*(*info).dispswitch).ref_detspect_scale], labels=['Scale main detspect','Scale ref detspect']
	CRISPEX_DISPRANGE_LS_SCALE, event
END


PRO CRISPEX_DISPRANGE_LS_SCALE, event
; Handles the scaling (or not) of the detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_SCALE'
	IF (*(*info).ctrlsswitch).imrefdetspect THEN CRISPEX_DISPRANGE_LS_SCALE_REF, event ELSE CRISPEX_DISPRANGE_LS_SCALE_MAIN, event
	CRISPEX_DISPRANGE_LS_RANGE, event
END

PRO CRISPEX_DISPRANGE_LS_SCALE_MAIN, event
; Handles the scaling (or not) of the main detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_SCALE_MAIN'
	IF (*(*info).dispswitch).detspect_scale THEN BEGIN
		IF (*(*info).plotswitch).scalestokes THEN ms = (*(*info).dataparams).ms * ((*(*info).paramparams).scale_cubes)[0] ELSE ms = ((*(*info).dataparams).ms) * ((*(*info).paramparams).scale_cubes)[0]
		(*(*(*info).plotaxes).ls_low_y) /= ms 
		(*(*(*info).plotaxes).ls_upp_y) /= ms
	ENDIF ELSE BEGIN
		IF (*(*info).plotswitch).scalestokes THEN ms = (*(*info).dataparams).ms ELSE ms = ((*(*info).dataparams).ms)
		(*(*(*info).plotaxes).ls_low_y) *= ms * ((*(*info).paramparams).scale_cubes)[0]
		(*(*(*info).plotaxes).ls_upp_y) *= ms * ((*(*info).paramparams).scale_cubes)[0]
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],2)
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],2)
END

PRO CRISPEX_DISPRANGE_LS_SCALE_REF, event
; Handles the scaling (or not) of the reference detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_SCALE_REF'
	IF (*(*info).dispswitch).ref_detspect_scale THEN BEGIN
		(*(*info).plotaxes).ls_low_y_ref /= ((*(*info).dataparams).refms * ((*(*info).paramparams).scale_cubes)[1]) 
		(*(*info).plotaxes).ls_upp_y_ref /= ((*(*info).dataparams).refms * ((*(*info).paramparams).scale_cubes)[1])
	ENDIF ELSE BEGIN
		(*(*info).plotaxes).ls_low_y_ref *= (*(*info).dataparams).refms * ((*(*info).paramparams).scale_cubes)[1]
		(*(*info).plotaxes).ls_upp_y_ref *= (*(*info).dataparams).refms * ((*(*info).paramparams).scale_cubes)[1]
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_low_y_ref,2)
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, SET_VALUE = STRTRIM((*(*info).plotaxes).ls_upp_y_ref,2)
END

PRO CRISPEX_DISPRANGE_LS_SUBTRACT, event
; Handles display of average subtracted spectrum in detailed spectrum window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LS_SUBTRACT'
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).plotswitch).ref_subtract = event.SELECT ELSE (*(*info).plotswitch).subtract = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).plotswitch).subtract,(*(*info).plotswitch).ref_subtract], labels=['Subtract from main','Subtract from ref']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DISPRANGE_T_LOW, event
; Handles change in lower t-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_T_LOW'
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, GET_VALUE = textvalue
	(*(*info).dispparams).t_low = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).t_low GE (*(*info).dispparams).t_upp) THEN (*(*info).dispparams).t_low = (*(*info).dispparams).t_upp - 1
	IF ((*(*info).dispparams).t_low LT (*(*info).dispparams).t_first) THEN (*(*info).dispparams).t_low = (*(*info).dispparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	IF (*(*info).intparams).lock_t THEN BEGIN
		(*(*info).plotaxes).int_low_t = (*(*info).dispparams).t_low
		CRISPEX_DISPRANGE_INT_T_RANGE, event
		IF (*(*info).winswitch).showint THEN WIDGET_CONTROL, (*(*info).ctrlsint).lower_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	ENDIF ELSE CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_T_UPP, event
; Handles change in upper t-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_T_UPP'
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, GET_VALUE = textvalue
	(*(*info).dispparams).t_upp = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).t_upp LE (*(*info).dispparams).t_low) THEN (*(*info).dispparams).t_upp = (*(*info).dispparams).t_low + 1
	IF ((*(*info).dispparams).t_upp GT (*(*info).dispparams).t_last) THEN (*(*info).dispparams).t_upp = (*(*info).dispparams).t_last
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	IF (*(*info).intparams).lock_t THEN BEGIN
		(*(*info).plotaxes).int_upp_t = (*(*info).dispparams).t_upp
		CRISPEX_DISPRANGE_INT_T_RANGE, event
		IF (*(*info).winswitch).showint THEN WIDGET_CONTROL, (*(*info).ctrlsint).upper_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	ENDIF ELSE CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_T_RANGE, event, NO_DRAW=no_draw
; Determines range from change in lower or upper t-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_T_RANGE'
	IF (*(*info).winswitch).showretrdet THEN BEGIN
		(*(*info).dispparams).t_low = (*(*(*info).detparams).t_restored)[(*(*info).detparams).idx] - (*(*info).detparams).delta_t_dn > (*(*info).dispparams).t_first
		(*(*info).dispparams).t_upp = (*(*(*info).detparams).t_restored)[(*(*info).detparams).idx] + (*(*info).detparams).delta_t_up < (*(*info).dispparams).t_last
	ENDIF
	(*(*info).dispparams).t_range = (*(*info).dispparams).t_upp - (*(*info).dispparams).t_low + 1
	IF ((*(*info).dispparams).t_range NE (*(*info).dataparams).nt) THEN WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, /SENSITIVE ELSE WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
	IF ((*(*info).winswitch).showretrdet EQ 0) THEN (*(*info).dataparams).t = (*(*info).dispparams).t_low
	WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_SLIDER_MIN = (*(*info).dispparams).t_low, SET_SLIDER_MAX = (*(*info).dispparams).t_upp, SET_VALUE = (*(*info).dataparams).t
	IF ((*(*info).dispparams).t_range - 1 EQ 1) THEN BEGIN
		t_step = 1
		t_sens = 0 
	ENDIF ELSE BEGIN
		t_step = (*(*info).pbparams).t_step
		t_sens = 1
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp], labels=['Lower t-value','Upper t-value']
	WIDGET_CONTROL, (*(*info).ctrlscp).t_step_slider, SET_SLIDER_MAX = (*(*info).dispparams).t_range - 1, SET_VALUE = t_step, SENSITIVE = t_sens
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		IF (*(*info).winswitch).showsp THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
		IF (*(*info).winswitch).showrefsp THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
		IF (*(*info).winswitch).showloop THEN CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event
		IF (*(*info).winswitch).showrestloop THEN CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event
		IF (*(*info).winswitch).showretrdet THEN CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
		CRISPEX_UPDATE_T, event
		CRISPEX_DRAW, event
	ENDIF
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF (*(*info).dataswitch).onecube THEN WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SET_VALUE = 'Update spectral windows', SENSITIVE = 1 ELSE WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1
	ENDIF
END

PRO CRISPEX_DISPRANGE_T_RESET, event, NO_DRAW=no_draw
; Handles reset of temporal boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_T_RESET'
	(*(*info).dispparams).t_upp = (*(*info).dispparams).t_last
	(*(*info).dispparams).t_low = (*(*info).dispparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
	IF ((*(*info).winswitch).showint AND (*(*info).intparams).lock_t) THEN CRISPEX_DISPRANGE_INT_T_RESET, event ELSE CRISPEX_DISPRANGE_T_RANGE, event, NO_DRAW=no_draw
END

PRO CRISPEX_DISPRANGE_LP_LOW, event
; Handles change in lower s-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_LOW'
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, GET_VALUE = textvalue
	(*(*info).dispparams).lp_low = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).lp_low GE (*(*info).dispparams).lp_upp) THEN (*(*info).dispparams).lp_low = (*(*info).dispparams).lp_upp - 1
	IF ((*(*info).dispparams).lp_low LT (*(*info).dispparams).lp_first) THEN (*(*info).dispparams).lp_low = (*(*info).dispparams).lp_first
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2)
	CRISPEX_DISPRANGE_LP_RANGE, event
END

PRO CRISPEX_DISPRANGE_LP_UPP, event
; Handles change in upper s-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_UPP'
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, GET_VALUE = textvalue
	(*(*info).dispparams).lp_upp = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).lp_upp LE (*(*info).dispparams).lp_low) THEN (*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_low + 1
	IF ((*(*info).dispparams).lp_upp GT (*(*info).dispparams).lp_last) THEN (*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_last
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2)
	CRISPEX_DISPRANGE_LP_RANGE, event
END

PRO CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
; Determines range from change in lower or upper s-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_RANGE'
	(*(*info).dispparams).lp_range = (*(*info).dispparams).lp_upp - (*(*info).dispparams).lp_low + 1
	IF ((*(*info).dispparams).lp_range NE (*(*info).dataparams).nlp) THEN WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, /SENSITIVE ELSE WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, SENSITIVE = 0
	IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN (*(*info).dataparams).lp = (*(*info).dispparams).lp_low ELSE $
		IF ((*(*info).dataparams).lp GT (*(*info).dispparams).lp_upp) THEN (*(*info).dataparams).lp = (*(*info).dispparams).lp_upp ELSE $
		(*(*info).dataparams).lp = (*(*info).dataparams).lp
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_SLIDER_MIN = (*(*info).dispparams).lp_low, SET_SLIDER_MAX = (*(*info).dispparams).lp_upp, SET_VALUE = (*(*info).dataparams).lp
	IF ((*(*info).dispparams).lp_range - 1 EQ 1) THEN BEGIN
		lp_step = 1
		lp_sens = 0 
	ENDIF ELSE BEGIN
		lp_step = (*(*info).pbparams).lp_step
		lp_sens = 1
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_slider, SET_SLIDER_MAX = (*(*info).dispparams).lp_range - 1, SET_VALUE = lp_step, SENSITIVE = lp_sens
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		IF (*(*info).winswitch).showsp THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
		CRISPEX_UPDATE_T, event
		CRISPEX_DRAW, event
	ENDIF
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF (*(*info).dataswitch).onecube THEN WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SET_VALUE = 'Update spectral windows', SENSITIVE = 1 ELSE WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp], labels=['Lower lp-value','Upper lp-value']
END

PRO CRISPEX_DISPRANGE_LP_RESET, event, NO_DRAW=no_draw
; Handles reset of spectral boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_RESET'
	(*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_last
	(*(*info).dispparams).lp_low = (*(*info).dispparams).lp_first
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, /SENSITIVE
	CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
END

PRO CRISPEX_DISPRANGE_LP_REF_RANGE, event, NO_DRAW=no_draw
; Determines range from change in lower or upper s-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_REF_RANGE'
	(*(*info).dispparams).lp_ref_range = (*(*info).dispparams).lp_ref_upp - (*(*info).dispparams).lp_ref_low + 1
	IF ((*(*info).dataparams).lp_ref LT (*(*info).dispparams).lp_ref_low) THEN (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_low ELSE $
		IF ((*(*info).dataparams).lp_ref GT (*(*info).dispparams).lp_ref_upp) THEN (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_upp 
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SET_SLIDER_MIN = (*(*info).dispparams).lp_ref_low, SET_SLIDER_MAX = (*(*info).dispparams).lp_ref_upp, SET_VALUE = (*(*info).dataparams).lp_ref
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		IF (*(*info).winswitch).showrefsp THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
		CRISPEX_UPDATE_T, event
		CRISPEX_DRAW, event
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).lp_ref_low,(*(*info).dispparams).lp_ref_upp], labels=['Lower lp-value','Upper lp-value']
END

PRO CRISPEX_DISPRANGE_LP_REF_RESET, event, NO_DRAW=no_draw
; Handles reset of reference spectral boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DISPRANGE_LP_REF_RESET'
	(*(*info).dispparams).lp_ref_upp = (*(*info).dispparams).lp_ref_last
	(*(*info).dispparams).lp_ref_low = (*(*info).dispparams).lp_ref_first
	CRISPEX_DISPRANGE_LP_REF_RANGE, event, NO_DRAW=no_draw
END
;================================================================================= DISPLAY DRAW PROCEDURES
PRO CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, no_cursor=no_cursor, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize, draw_mask=draw_mask
; Handles overplotting of the cursor, slits and loop paths
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_CURSCROSS_PLOT'
	(*(*info).phiparams).d_nphi_set = (*(*info).zooming).factor * (*(*info).phiparams).nphi_set
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF ((*(*info).zooming).factor EQ 1) THEN $
			PLOTS, ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywinx/(2.*(*(*info).dataparams).nx))*COS((*(*info).phiparams).angle*!DTOR) + (*(*info).curs).sx, $
				([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywiny/(2.*(*(*info).dataparams).ny))*SIN((*(*info).phiparams).angle*!DTOR) + (*(*info).curs).sy, /DEVICE, COLOR = !P.COLOR $
		ELSE PLOTS, ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywinx/(2.*(*(*info).dataparams).nx))*COS((*(*info).phiparams).angle*!DTOR) + (*(*info).curs).sx, $
			([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywiny/(2.*(*(*info).dataparams).ny))*SIN((*(*info).phiparams).angle*!DTOR) + (*(*info).curs).sy, /DEVICE, COLOR = !P.COLOR
	ENDIF
	IF (*(*info).overlayswitch).loopslit AND ((*(*info).loopparams).np GT 0) THEN BEGIN
		CRISPEX_ZOOM_LOOP, event
		IF ~KEYWORD_SET(NO_ENDPOINTS) THEN PLOTS, *(*(*info).overlayparams).sxp, *(*(*info).overlayparams).syp, /DEVICE, COLOR = !P.COLOR, PSYM = 1, THICK=thick, SYMSIZE=symsize
		IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN PLOTS,*(*(*info).overlayparams).sxr,*(*(*info).overlayparams).syr, /DEVICE, COLOR = !P.COLOR, PSYM = 3, THICK=thick $
		ELSE PLOTS,*(*(*info).overlayparams).sxr,*(*(*info).overlayparams).syr,/DEVICE, COLOR = !P.COLOR, LINESTYLE = (*(*info).overlayparams).loop_linestyle, THICK=thick
	ENDIF ELSE IF ((*(*info).meas).np GE 1) THEN BEGIN
		CRISPEX_ZOOM_MEAS, event
		IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
			PLOTS, *(*(*info).meas).sxp,*(*(*info).meas).syp, /DEVICE, COLOR = !P.COLOR, PSYM = 1, THICK=thick, SYMSIZE=symsize
			PLOTS, *(*(*info).meas).sxp,*(*(*info).meas).syp, /DEVICE, COLOR = !P.COLOR, LINESTYLE = 0, THICK=thick, SYMSIZE=symsize
		ENDIF
	ENDIF ELSE IF ~KEYWORD_SET(NO_CURSOR) THEN PLOTS, (*(*info).curs).sx,(*(*info).curs).sy, /DEVICE, PSYM = 1, COLOR = curscolor, THICK=thick, SYMSIZE=symsize
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).curs).sx,(*(*info).curs).sy,(*(*info).winsizes).xywinx,(*(*info).winsizes).xywiny], labels=['sx','sy','xywinx','xywiny']
	CRISPEX_DRAW_LOOP_OVERLAYS, event, NO_NUMBER=no_number, THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize
	IF draw_mask THEN CRISPEX_DRAW_MASK_OVERLAYS, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_0, event
; Handles setting of loop path linestyle to solid
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LOOP_LINESTYLE_0'
	(*(*info).overlayparams).loop_linestyle = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_1, event
; Handles setting of loop path linestyle to dotted (i.e. actual loop points, not dots between the vertices)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LOOP_LINESTYLE_1'
	(*(*info).overlayparams).loop_linestyle = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_2, event
; Handles setting of loop path linestyle to dashed
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LOOP_LINESTYLE_2'
	(*(*info).overlayparams).loop_linestyle = 2
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_OVERLAYS, event, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize
; Handles overplotting of loop paths from the restored and retrieved loops as well as from the retrieved detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LOOP_OVERLAYS'
	IF (((*(*info).loopswitch).restore_loops EQ 1) AND ((*(*info).restoreparams).cfilecount GT 0)) THEN BEGIN
		drawrestore = 0
		FOR i=0,(*(*info).restoreparams).cfilecount-1 DO BEGIN
			IF (((*(*info).overlayswitch).overlalways EQ 1) OR ((*(*(*info).restoreparams).lp_restored)[i] EQ (*(*info).dataparams).lp) AND ((*(*(*info).restoreparams).sel_loops)[i] EQ 1)) THEN BEGIN
				drawrestore += 1
				xp_orig = (*(*(*info).restoreparams).xp_restored)[(0+(*(*(*info).restoreparams).psizes)[i]):((*(*(*info).restoreparams).psizes)[i+1]-1)]
				yp_orig = (*(*(*info).restoreparams).yp_restored)[(0+(*(*(*info).restoreparams).psizes)[i]):((*(*(*info).restoreparams).psizes)[i+1]-1)]
				xr_orig = (*(*(*info).restoreparams).xr_restored)[(0+(*(*(*info).restoreparams).rsizes)[i]):((*(*(*info).restoreparams).rsizes)[i+1]-1)]
				yr_orig = (*(*(*info).restoreparams).yr_restored)[(0+(*(*(*info).restoreparams).rsizes)[i]):((*(*(*info).restoreparams).rsizes)[i+1]-1)]
				IF ((*(*info).zooming).factor NE 1) THEN BEGIN
					sxp = (xp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
					syp = (yp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
					sxr = (xr_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
					syr = (yr_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
				ENDIF ELSE BEGIN
					sxp = xp_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
					syp = yp_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
					sxr = xr_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx	
					syr = yr_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
				ENDELSE
				sxp_last = sxp[(SIZE(sxp))[1]-1]-1.5*(*(*info).zooming).factor
				syp_last = syp[(SIZE(sxp))[1]-1]+1.5*(*(*info).zooming).factor
				IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
					PLOTS, sxp, syp, PSYM = 1, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
					PLOTS, sxp, syp, PSYM = 4, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
				ENDIF
				IF ~KEYWORD_SET(NO_NUMBER) THEN XYOUTS,sxp_last,syp_last,STRTRIM(i,2), /DEVICE
				IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN PLOTS, sxr, syr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
				ELSE PLOTS, sxr, syr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR, /DEVICE, THICK=thick
			ENDIF
		ENDFOR
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).restoreparams).cfilecount,drawrestore], labels=['Loops restored','Loops drawn']
	ENDIF 
	IF (((*(*info).retrparams).clfilecount GT 0) AND ((*(*info).loopswitch).retrieve_loops EQ 1)) THEN BEGIN
		drawretr = 0
		FOR k=0,(*(*info).retrparams).clfilecount-1 DO BEGIN
			IF ((*(*(*info).retrparams).sel_loops)[k] EQ 1) THEN BEGIN
				drawretr += 1
				xlp_orig = (*(*(*info).retrparams).xlp_restored)[(0+(*(*(*info).retrparams).lpsizes)[k]):((*(*(*info).retrparams).lpsizes)[k+1]-1)]
				ylp_orig = (*(*(*info).retrparams).ylp_restored)[(0+(*(*(*info).retrparams).lpsizes)[k]):((*(*(*info).retrparams).lpsizes)[k+1]-1)]
				xlr_orig = (*(*(*info).retrparams).xlr_restored)[(0+(*(*(*info).retrparams).lrsizes)[k]):((*(*(*info).retrparams).lrsizes)[k+1]-1)]
				ylr_orig = (*(*(*info).retrparams).ylr_restored)[(0+(*(*(*info).retrparams).lrsizes)[k]):((*(*(*info).retrparams).lrsizes)[k+1]-1)]
				IF ((*(*info).zooming).factor NE 1) THEN BEGIN
					sxlp = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
					sylp = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
					sxlr = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
					sylr = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
				ENDIF ELSE BEGIN
					sxlp = xlp_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
					sylp = ylp_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
				sxlr = xlr_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
					sylr = ylr_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
				ENDELSE
				sxlp_last = sxlp[(SIZE(sxlp))[1]-1]-1.5*(*(*info).zooming).factor
				sylp_last = sylp[(SIZE(sxlp))[1]-1]+1.5*(*(*info).zooming).factor
				IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
					PLOTS, sxlp, sylp, PSYM = 1, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
					PLOTS, sxlp, sylp, PSYM = 4, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
				ENDIF
				IF ~KEYWORD_SET(NO_NUMBER) THEN XYOUTS,sxlp_last,sylp_last,'L'+STRTRIM(k,2), /DEVICE
				IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN PLOTS, sxlr, sylr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
				ELSE PLOTS, sxlr, sylr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR, /DEVICE, THICK=thick
			ENDIF
		ENDFOR
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).retrparams).clfilecount,drawretr], labels=['Loops retrieved','Retrieved loops drawn']
	ENDIF 
	IF ((*(*info).loopswitch).retrieve_detfile EQ 1) THEN BEGIN
		low = (*(*info).detparams).mid-FLOOR((*(*info).detparams).width/2.)
		upp = (*(*info).detparams).mid+FLOOR((*(*info).detparams).width/2.)
		condition = WHERE(*(*(*info).detparams).sel_dets EQ 1, conditioncount)
		IF (*(*info).overlayswitch).det_overlay_all THEN BEGIN
			FOR j=0,(*(*info).detparams).nr_dets-1 DO CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, j, low, upp, NO_NUMBER=no_number, THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize
			detdrawn = (*(*info).detparams).nr_dets
		ENDIF ELSE IF (((*(*info).overlayswitch).det_overlay_all EQ 0) AND (N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) THEN BEGIN
				indices = WHERE((*(*(*info).detparams).sel_dets) EQ 1)
				FOR m=0,N_ELEMENTS(condition)-1 DO CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, indices[m], low, upp, NO_NUMBER=no_number, THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize
				detdrawn = N_ELEMENTS(indices)
		ENDIF ELSE detdrawn = 0
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).detparams).nr_dets,conditioncount,detdrawn], labels=['Detections retrieved','Detections selected','Detections drawn']
	ENDIF ELSE RETURN
END

PRO CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, j, low, upp, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize
; Handles the actual overplotting of loop paths from the retrieved detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS'
	xlp_orig = (*(*(*info).detparams).xlp_restored)[(0+(*(*(*info).detparams).lpsizes)[j]):((*(*(*info).detparams).lpsizes)[j+1]-1)]
	ylp_orig = (*(*(*info).detparams).ylp_restored)[(0+(*(*(*info).detparams).lpsizes)[j]):((*(*(*info).detparams).lpsizes)[j+1]-1)]
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		sxlp = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		sylp = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
	ENDIF ELSE BEGIN
		sxlp = xlp_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
		sylp = ylp_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
	ENDELSE
	sxlp_last = sxlp[(SIZE(sxlp))[1]-1]-1.5*(*(*info).zooming).factor
	sylp_last = sylp[(SIZE(sxlp))[1]-1]+1.5*(*(*info).zooming).factor
	IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
		PLOTS, sxlp, sylp, PSYM = 1, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
		PLOTS, sxlp, sylp, PSYM = 4, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
	ENDIF
	IF ~KEYWORD_SET(NO_NUMBER) THEN XYOUTS,sxlp_last,sylp_last,'D'+STRTRIM(j,2), /DEVICE
	FOR i=low,upp DO BEGIN
		xlr_orig = (*(*(*info).detparams).xlr_restored)[(0+(*(*(*info).detparams).lrsizes)[j]):((*(*(*info).detparams).lrsizes)[j+1]-1),i]
		ylr_orig = (*(*(*info).detparams).ylr_restored)[(0+(*(*(*info).detparams).lrsizes)[j]):((*(*(*info).detparams).lrsizes)[j+1]-1),i]
		IF ((*(*info).zooming).factor NE 1) THEN BEGIN
			sxlr = (xlr_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
			sylr = (ylr_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
		ENDIF ELSE BEGIN
			sxlr = xlr_orig * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
			sylr = ylr_orig * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
		ENDELSE
		IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN PLOTS, sxlr, sylr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
		ELSE PLOTS, sxlr, sylr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR,/DEVICE, THICK=thick
	ENDFOR
END

PRO CRISPEX_DRAW_MASK_OVERLAYS, event
; Handles the overlay of a mask
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_MASK_OVERLAYS'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		x_low = (*(*info).zooming).xpos
		x_upp = (*(*info).zooming).xpos + (*(*info).dataparams).d_nx
		y_low = (*(*info).zooming).ypos
		y_upp = (*(*info).zooming).ypos + (*(*info).dataparams).d_ny
	ENDIF ELSE BEGIN
		x_low = 0
		x_upp = (*(*info).dataparams).nx-1
		y_low = 0
		y_upp = (*(*info).dataparams).ny-1
	ENDELSE
	LOADCT, (*(*info).overlayparams).maskct, /SILENT
	CONTOUR,*(*(*info).data).maskslice,COLOR=(*(*info).overlayparams).maskcolor, LEVELS = 1, /ISOTROPIC, XS=13,YS=13,POSITION=[0,0,1,1],/NORMAL, /NOERASE
	LOADCT, 0, /SILENT
END

PRO CRISPEX_DRAW, event
; Handles the actual drawing of the data into the respective open display windows
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW'
	IF (*(*info).winswitch).showparam THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsparam).x_coord_val, SET_VALUE = STRTRIM(FIX((*(*info).dataparams).x),2)
		WIDGET_CONTROL, (*(*info).ctrlsparam).y_coord_val, SET_VALUE = STRTRIM(FIX((*(*info).dataparams).y),2)
		WIDGET_CONTROL, (*(*info).ctrlsparam).t_coord_val, SET_VALUE = STRTRIM((*(*info).dataparams).t,2)
		WIDGET_CONTROL, (*(*info).ctrlsparam).lp_coord_val, SET_VALUE = STRTRIM(LONG((*(*info).dataparams).lp),2)
		WIDGET_CONTROL, (*(*info).ctrlsparam).zoom_val, SET_VALUE = STRTRIM(LONG((*(*info).zooming).factor),2)
		IF (*(*info).paramswitch).dt_set THEN WIDGET_CONTROL, (*(*info).ctrlsparam).act_t_val, SET_VALUE = STRTRIM(STRING((*(*info).dataparams).t * (*(*info).plotaxes).dt, FORMAT='(3(F9.2,x))'),2) 
		IF (*(*info).plotswitch).v_dop_set THEN BEGIN
			WIDGET_CONTROL, (*(*info).ctrlsparam).act_lp_val, SET_VALUE = STRTRIM(STRING((*(*info).dataparams).lps[(*(*info).dataparams).lp], FORMAT = '(3(F9.1,x))'),2)
			WIDGET_CONTROL, (*(*info).ctrlsparam).v_dop_val, SET_VALUE = STRTRIM(STRING((*(*info).plotaxes).v_dop[(*(*info).dataparams).lp],FORMAT = '(3(F9.2,x))'),2)
		ENDIF ELSE IF (*(*info).plotswitch).heightset THEN BEGIN
			WIDGET_CONTROL, (*(*info).ctrlsparam).act_lp_val, SET_VALUE = STRTRIM(STRING((*(*info).dataparams).lps[(*(*info).dataparams).lp], FORMAT = '(3(F9.1,x))'),2)
		ENDIF
		IF ((*(*info).dataparams).refnlp GT 1) THEN BEGIN
			WIDGET_CONTROL, (*(*info).ctrlsparam).lp_ref_coord_val, SET_VALUE = STRTRIM(LONG((*(*info).dataparams).lp_ref),2)
			IF (*(*info).plotswitch).v_dop_set_ref THEN BEGIN
				WIDGET_CONTROL, (*(*info).ctrlsparam).act_lp_ref_val, SET_VALUE = STRTRIM(STRING((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref], FORMAT = '(3(F9.1,x))'),2)
				WIDGET_CONTROL, (*(*info).ctrlsparam).v_dop_ref_val, SET_VALUE = STRTRIM(STRING((*(*info).plotaxes).v_dop_ref[(*(*info).dataparams).lp_ref],FORMAT = '(3(F9.2,x))'),2)
			ENDIF ELSE IF (*(*info).plotswitch).refheightset THEN BEGIN
				WIDGET_CONTROL, (*(*info).ctrlsparam).act_lp_ref_val, SET_VALUE = STRTRIM(STRING((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref], FORMAT = '(3(F9.1,x))'),2)
			ENDIF
		ENDIF		
		IF (*(*info).paramswitch).ref_get THEN BEGIN
			IF ((*(*info).dataparams).refnt EQ 0) THEN refdata = *(*(*info).data).refdata ELSE IF ((*(*info).dataparams).refnt EQ 1) THEN refdata = (*(*(*info).data).refdata)[0] ELSE $
					refdata = (*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref]
			(*(*info).paramparams).ref_get_val = refdata[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y)] * ((*(*info).paramparams).scale_cubes)[1]
			IF ((FLOOR(ALOG10(ABS((*(*info).paramparams).ref_get_val))) LE -2) OR (FLOOR(ALOG10(ABS((*(*info).paramparams).ref_get_val))) GE 3)) THEN ref_val_format = '(3(E10.4,x))' ELSE ref_val_format = '(3(F9.2,x))'
			(*(*info).paramparams).refsc_get_val = refdata[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y)]/((*(*info).dataparams).refms)
			IF ((FLOOR(ALOG10(ABS((*(*info).paramparams).refsc_get_val))) LE -2) OR (FLOOR(ALOG10(ABS((*(*info).paramparams).refsc_get_val))) GE 3)) THEN refsc_val_format = '(3(E10.4,x))' $
				ELSE refsc_val_format = '(3(F9.2,x))'
			WIDGET_CONTROL, (*(*info).ctrlsparam).ref_val, SET_VALUE = STRTRIM(STRING((*(*info).paramparams).ref_get_val,FORMAT=ref_val_format),2)
			WIDGET_CONTROL, (*(*info).ctrlsparam).refsc_val, SET_VALUE = STRTRIM(STRING((*(*info).paramparams).refsc_get_val,FORMAT=refsc_val_format),2)
		ENDIF
		IF (*(*info).paramswitch).img_get THEN BEGIN
			IF (*(*info).plotswitch).scalestokes THEN ms = (*(*info).dataparams).ms ELSE ms = ((*(*info).dataparams).ms)[(*(*info).dataparams).s]
			IF ((*(*info).dataswitch).spfile EQ 1) OR (*(*info).dataswitch).onecube THEN $
				imdata = (*(*(*info).data).imagedata)[(*(*info).dataparams).t*(*(*info).dataparams).nlp*(*(*info).dataparams).ns + (*(*info).dataparams).s*(*(*info).dataparams).nlp + (*(*info).dataparams).lp] $
				ELSE imdata = (*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp] 
			(*(*info).paramparams).img_get_val = imdata[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y)] * ((*(*info).paramparams).scale_cubes)[0]
			IF ((FLOOR(ALOG10(ABS((*(*info).paramparams).img_get_val))) LE -2) OR (FLOOR(ALOG10(ABS((*(*info).paramparams).img_get_val)) GE 3))) THEN img_val_format = '(3(E10.4,x))' ELSE img_val_format = '(3(F9.2,x))'
			(*(*info).paramparams).imgsc_get_val = imdata[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y)]/ms
			IF ((FLOOR(ALOG10(ABS((*(*info).paramparams).imgsc_get_val))) LE -2) OR (FLOOR(ALOG10(ABS((*(*info).paramparams).imgsc_get_val))) GE 3)) THEN imgsc_val_format = '(3(E10.4,x))' $
				ELSE imgsc_val_format = '(3(F9.2,x))'
			IF ((*(*info).paramswitch).ref_get EQ 0) THEN extra_label = '' ELSE extra_label = ','
			WIDGET_CONTROL, (*(*info).ctrlsparam).img_val, SET_VALUE = STRTRIM(STRING((*(*info).paramparams).img_get_val,FORMAT=img_val_format),2)
			WIDGET_CONTROL, (*(*info).ctrlsparam).imgsc_val, SET_VALUE = STRTRIM(STRING((*(*info).paramparams).imgsc_get_val,FORMAT=imgsc_val_format),2)+extra_label
		ENDIF
	ENDIF
	IF ((*(*info).curs).lockset AND ((*(*info).overlayswitch).loopslit NE 1)) THEN BEGIN
		(*(*info).curs).sx = (*(*info).curs).sxlock	&	(*(*info).curs).sy = (*(*info).curs).sylock
	ENDIF 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paramswitch).img_get,(*(*info).paramparams).img_get_val], labels=['Get image value','Image value']
	CRISPEX_DRAW_IMREF, event
	CRISPEX_DRAW_SPECTRAL, event
	CRISPEX_DRAW_TIMESLICES, event
	CRISPEX_DRAW_OTHER, event
END

PRO CRISPEX_DRAW_IMREF, event
; (Re)draw main and reference image window procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_IMREF'
	CRISPEX_DRAW_XY, event
	IF (*(*info).winswitch).showref THEN CRISPEX_DRAW_REF, event
	IF (*(*info).winswitch).showimref THEN CRISPEX_DRAW_IMREF_BLINK, event
	IF (*(*info).winswitch).showdop THEN CRISPEX_DRAW_DOPPLER, event
END

PRO CRISPEX_DRAW_SPECTRAL, event
; (Re)draw spectral windows procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_SPECTRAL'
	IF (*(*info).winswitch).showls THEN CRISPEX_DRAW_LS, event
	IF (*(*info).winswitch).showrefls THEN CRISPEX_DRAW_REFLS, event
	IF (*(*info).winswitch).showphis AND (((*(*info).pbparams).mode EQ 'PAUSE') OR (*(*info).dispparams).phislice_update) THEN CRISPEX_DRAW_PHIS, event		
END

PRO CRISPEX_DRAW_TIMESLICES, event
; (Re)draw timeslices procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_TIMESLICES'
	IF (*(*info).winswitch).showsp THEN CRISPEX_DRAW_SP, event
	IF (*(*info).winswitch).showrefsp THEN CRISPEX_DRAW_REFSP, event
	IF (*(*info).winswitch).showloop THEN CRISPEX_DRAW_LOOPSLAB, event 
	IF (*(*info).winswitch).showrefloop THEN CRISPEX_DRAW_REFLOOPSLAB, event 
	IF (*(*info).winswitch).showrestloop THEN CRISPEX_DRAW_REST_LOOP, event
	IF (*(*info).winswitch).showretrdet THEN CRISPEX_DRAW_RETR_DET, event
END

PRO CRISPEX_DRAW_OTHER, event
; (Re)draw other windows procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_OTHER'
	IF (*(*info).winswitch).showint THEN CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DRAW_SUBCOLOR, event, imref, subcolor, minimum, maximum, xyrange=xyrange
; Determines the color beneath the cursor to get the cursor color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_SUBCOLOR'
	IF (N_ELEMENTS(XYRANGE) EQ 0) THEN BEGIN
		x_upp = FIX((*(*info).dataparams).x) + 5 < FIX((*(*info).dispparams).x_last)
		x_low = FIX((*(*info).dataparams).x) - 5 > FIX((*(*info).dispparams).x_first)
		y_upp = FIX((*(*info).dataparams).y) + 5 < FIX((*(*info).dispparams).y_last)
		y_low = FIX((*(*info).dataparams).y) - 5 > FIX((*(*info).dispparams).y_first)
	ENDIF ELSE BEGIN
		x_low = xyrange[0]	&	x_upp = xyrange[1]
		y_low = xyrange[2]	&	y_upp = xyrange[3]
	ENDELSE
	IF (imref EQ 0) THEN $
		selected_data = (*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp] $
	ELSE IF (imref EQ 1) THEN BEGIN
		IF ((*(*info).dataparams).refnt GT 1) THEN selected_data = (*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref] ELSE $
			selected_data = (*(*(*info).data).refdata)[(*(*info).dataparams).lp_ref]
	ENDIF ELSE IF (imref EQ 2) THEN	selected_data = *(*(*info).data).dopslice
	scaled_data = BYTSCL(selected_data[x_low:x_upp, y_low:y_upp], MIN = minimum, MAX = maximum)
	subcolor = MEAN(scaled_data)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [x_low,x_upp,y_low,y_upp,subcolor], labels=['x_low','x_upp','y_low','y_upp','subcolor']
END

PRO CRISPEX_DRAW_XY, event, no_cursor=no_cursor, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize, asecbar=asecbar
; (Re)draw main image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_XY'
	WSET, (*(*info).winids).imwid
	CRISPEX_DRAW_XYSLICE_SCALING, event, imdisp, minimum, maximum
	IF ((*(*info).zooming).factor EQ 1) THEN TV, CONGRID(imdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE TV, imdisp 
	CRISPEX_DRAW_SUBCOLOR, event, 0, subcolor, minimum, maximum
	IF (subcolor GE 122) THEN curscolor = 0 ELSE curscolor = 255
	CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, NO_CURSOR=no_cursor, NO_NUMBER=no_number, THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize, $
		DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[0])
	IF KEYWORD_SET(ASECBAR) THEN BEGIN
		xlow = 25							& 	ylow = 25 
		xupp = xlow + (*(*info).savparams).overlays_asecbar_pix		&	yupp = ylow + 15 + (*(*info).savparams).overlays_thick
		asecbarcol = 255
		PLOTS, [xlow,xupp],[ylow,ylow], THICK=(*(*info).savparams).overlays_thick, COLOR=asecbarcol, /DEVICE
		XYOUTS, (xupp-xlow)/2.+xlow, ylow+5, STRTRIM((*(*info).savparams).overlays_asecbar_length,2)+'"', /DEVICE, COLOR=asecbarcol, ALIGN=0.5,	CHARSIZE=(*(*info).savparams).overlays_symsize, $
			CHARTHICK=(*(*info).savparams).overlays_thick
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).imwid,curscolor], labels=['Window ID for draw','Main curscolor']
END

PRO CRISPEX_DRAW_XYSLICE_SCALING, event, finalimage, minimum, maximum
; Determines the minimum and maximum value for the main image scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_XYSLICE_SCALING'
		IF ((*(*(*info).scaling).imagescale)[0] EQ 0) THEN BEGIN
			IF ((*(*info).scaling).scalestokes_max AND ((*(*info).dataparams).s GE 1)) THEN BEGIN
				minimum = MIN(((*(*info).scaling).imagemin)[*,(*(*info).dataparams).s])
				maximum = MAX(((*(*info).scaling).imagemax)[*,(*(*info).dataparams).s])
			ENDIF ELSE BEGIN
				minimum = ((*(*info).scaling).imagemin)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
				maximum = ((*(*info).scaling).imagemax)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
			ENDELSE
		ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[0] EQ 1) THEN BEGIN
			minimum = MIN( *(*(*info).data).xyslice, MAX=maximum)
		ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[0] EQ 2) OR ((*(*(*info).scaling).imagescale)[0] EQ 3) THEN BEGIN
			IF (*(*(*info).scaling).relative)[0] THEN BEGIN
				minimum = (*(*(*info).scaling).scale_range)[0] / 100. * (*(*(*info).scaling).rel_scale_min_val)[0] + (*(*(*info).scaling).scale_minimum)[0]
				maximum = (*(*(*info).scaling).scale_range)[0] / 100. * (*(*(*info).scaling).rel_scale_max_val)[0] + (*(*(*info).scaling).scale_minimum)[0]
			ENDIF ELSE BEGIN
				minimum = (*(*(*info).scaling).scale_range)[0] / 255. * (*(*(*info).scaling).scale_min_val)[0] + (*(*(*info).scaling).scale_minimum)[0]
				maximum = (*(*(*info).scaling).scale_range)[0] / 255. * (*(*(*info).scaling).scale_max_val)[0] + (*(*(*info).scaling).scale_minimum)[0]
			ENDELSE
		ENDIF 
		finalimage = BYTSCL(*(*(*info).data).xyslice, MIN = minimum, MAX = maximum) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [minimum,maximum], labels=['minimum','maximum']
END

PRO CRISPEX_DRAW_DOPPLER, event
; (Re)draw Doppler-image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_DOPPLER'
	WSET, (*(*info).winids).dopwid
	IF (*(*info).dispswitch).drawdop THEN BEGIN
		CRISPEX_DRAW_DOPSLICE_SCALING, event, dopdisp, dopminimum, dopmaximum
		IF ((*(*info).zooming).factor EQ 1) THEN TV, CONGRID(dopdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE TV, dopdisp
		CRISPEX_DRAW_SUBCOLOR, event, 2, dopsubcolor, dopminimum, dopmaximum
		IF (dopsubcolor GE 122) THEN dopcurscolor = 0 ELSE dopcurscolor = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, dopcurscolor, DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[2])
	ENDIF ELSE BEGIN
		(*(*info).scaling).imrefscaling = 0
		CRISPEX_SCALING_SET_BUTTONS, event
		CRISPEX_SCALING_SET_SLIDERS, event
		WIDGET_CONTROL, (*(*info).ctrlscp).xy_scaling_but, /SET_BUTTON
		WIDGET_CONTROL, (*(*info).ctrlscp).dop_scaling_but, SET_BUTTON = 0
		TV, CONGRID(*(*(*info).data).emptydopslice,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny)
		IF ((*(*info).dataparams).lp GT (*(*info).dataparams).lp_dop) THEN BEGIN
			lp_blue = (*(*info).dataparams).lp_dop		&	lp_red = (*(*info).dataparams).lp
		ENDIF ELSE BEGIN
			lp_blue = (*(*info).dataparams).lp		&	lp_red = (*(*info).dataparams).lp_dop
		ENDELSE
		IF (lp_red EQ lp_blue) THEN extramessage = 'Same spectral position' ELSE extramessage = 'Spectral position outside set spectral range'
		XYOUTS,(*(*info).winsizes).xywinx/2.,(*(*info).winsizes).xywiny/2.,'Could not create Doppler image for selected spectral positions.!C'+extramessage+': (lp_blue,lp_red)=('+$
			STRTRIM(lp_blue,2)+','+STRTRIM(lp_red,2)+')', CHARSIZE = 1.2, COLOR = 255, ALIGNMENT = 0.5, /DEVICE
		dopcurscolor = 255
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).dop_scaling_but, SENSITIVE = (*(*info).dispswitch).drawdop
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).dopwid,(*(*info).dispswitch).drawdop,dopcurscolor], $
		labels=['Window ID for draw','Drawing Doppler image','Doppler curscolor']
END

PRO CRISPEX_DRAW_DOPSLICE_SCALING, event, finalimage, minimum, maximum
; Determines the minimum and maximum value for the Doppler-image scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_DOPSLICE_SCALING'
		IF ((*(*(*info).scaling).imagescale)[2] EQ 0) THEN BEGIN
			IF ((*(*info).scaling).scalestokes_max AND ((*(*info).dataparams).s GE 1)) THEN BEGIN
				minimum = MIN(((*(*info).scaling).dopplermin)[*,(*(*info).dataparams).s])
				maximum = MAX(((*(*info).scaling).dopplermax)[*,(*(*info).dataparams).s])
			ENDIF ELSE BEGIN
				minimum = ((*(*info).scaling).dopplermin)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
				maximum = ((*(*info).scaling).dopplermax)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
			ENDELSE
		ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[2] EQ 1) THEN minimum = MIN( *(*(*info).data).dopslice, MAX=maximum) ELSE $
			IF ((*(*(*info).scaling).imagescale)[2] EQ 2) OR ((*(*(*info).scaling).imagescale)[2] EQ 3) THEN BEGIN
			IF (*(*(*info).scaling).relative)[2] THEN BEGIN
				minimum = (*(*(*info).scaling).scale_range)[2] / 100. * (*(*(*info).scaling).rel_scale_min_val)[2] + (*(*(*info).scaling).scale_minimum)[2]
				maximum = (*(*(*info).scaling).scale_range)[2] / 100. * (*(*(*info).scaling).rel_scale_max_val)[2] + (*(*(*info).scaling).scale_minimum)[2]
			ENDIF ELSE BEGIN
				minimum = (*(*(*info).scaling).scale_range)[2] / 255. * (*(*(*info).scaling).scale_min_val)[2] + (*(*(*info).scaling).scale_minimum)[2]
				maximum = (*(*(*info).scaling).scale_range)[2] / 255. * (*(*(*info).scaling).scale_max_val)[2] + (*(*(*info).scaling).scale_minimum)[2]
			ENDELSE
		ENDIF 
		finalimage = BYTSCL(*(*(*info).data).dopslice, MIN = minimum, MAX = maximum) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [minimum,maximum], labels=['minimum','maximum']
END

PRO CRISPEX_DRAW_REF, event
; (Re)draw reference image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REF'
	WSET, (*(*info).winids).refwid
	CRISPEX_DRAW_REFSLICE_SCALING, event, refdisp, refmin, refmax
	IF ((*(*info).zooming).factor EQ 1) THEN TV, CONGRID(refdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE TV, refdisp 
	CRISPEX_DRAW_SUBCOLOR, event, 1, subcolor_ref, refmin, refmax
	IF (subcolor_ref GE 122) THEN curscolor_ref = 0 ELSE curscolor_ref = 255
	CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor_ref, DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[1])
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refwid,curscolor_ref], labels=['Window ID for draw','Reference curscolor']
END

PRO CRISPEX_DRAW_REFSLICE_SCALING, event, refdisp, refmin, refmax
; Determines the minimum and maximum value for the reference image scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REFSLICE_SCALING'
	IF ((*(*(*info).scaling).imagescale)[1] EQ 0) THEN BEGIN
		IF ((*(*info).dataparams).refnlp EQ 1) THEN BEGIN
			refmin = (*(*info).scaling).refmin 
			refmax = (*(*info).scaling).refmax
		ENDIF ELSE BEGIN
			refmin = ((*(*info).scaling).refmin)[(*(*info).dataparams).lp_ref]
			refmax = ((*(*info).scaling).refmax)[(*(*info).dataparams).lp_ref]
		ENDELSE
	ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[1] EQ 1) THEN refmin = MIN( *(*(*info).data).refslice, MAX=refmax) ELSE $
		IF ((*(*(*info).scaling).imagescale)[1] EQ 2) OR ((*(*(*info).scaling).imagescale)[1] EQ 3) THEN BEGIN
		IF (*(*(*info).scaling).relative)[1] THEN BEGIN
			refmin = (*(*(*info).scaling).scale_range)[1] / 100. * (*(*(*info).scaling).rel_scale_min_val)[1] + (*(*(*info).scaling).scale_minimum)[1]
			refmax = (*(*(*info).scaling).scale_range)[1] / 100. * (*(*(*info).scaling).rel_scale_max_val)[1] + (*(*(*info).scaling).scale_minimum)[1]
		ENDIF ELSE BEGIN
			refmin = (*(*(*info).scaling).scale_range)[1] / 255. * (*(*(*info).scaling).scale_min_val)[1] + (*(*(*info).scaling).scale_minimum)[1]
			refmax = (*(*(*info).scaling).scale_range)[1] / 255. * (*(*(*info).scaling).scale_max_val)[1] + (*(*(*info).scaling).scale_minimum)[1]
		ENDELSE
	ENDIF 
	refdisp = BYTSCL(*(*(*info).data).refslice, MIN = refmin, MAX = refmax)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [refmin,refmax], labels=['minimum','maximum']
END

PRO CRISPEX_DRAW_IMREF_BLINK, event
; Handles the blinking of the main and reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_IMREF_BLINK'
	WSET,(*(*info).winids).imrefwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).imrefwid], labels=['Window ID for draw']
	IF (SIZE((*(*info).plotaxes).dt,/TYPE) NE 2) THEN time = STRING((*(*info).dataparams).t*(*(*info).plotaxes).dt,FORMAT='(F'+STRTRIM(FLOOR(ALOG10(((*(*info).dataparams).t+1)*(*(*info).plotaxes).dt))+4,2)+'.2)')+' s' ELSE $
		time = STRTRIM((*(*info).dataparams).t,2)
	IF (*(*info).winids).imrefdisp THEN BEGIN
		CRISPEX_DRAW_REFSLICE_SCALING, event, refdisp, refmin, refmax
		IF ((*(*info).zooming).factor EQ 1) THEN TV, CONGRID(refdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE TV, refdisp
		CRISPEX_DRAW_SUBCOLOR, event, 1, subcolor_ref, refmin, refmax
		IF (subcolor_ref GE 122) THEN curscolor_ref = 0 ELSE curscolor_ref = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor_ref, DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[1])
		CRISPEX_DRAW_SUBCOLOR, event, 1, color_reftxt, refmin, refmax, xyrange=[10,100,5,20]
		IF (color_reftxt GE 122) THEN reftxtcol = 0 ELSE reftxtcol = 255
		XYOUTS, 10, 10, 'Reference image, t='+time, /DEVICE, COLOR=reftxtcol
	ENDIF ELSE BEGIN
		CRISPEX_DRAW_XYSLICE_SCALING, event, imdisp, minimum, maximum
		IF ((*(*info).zooming).factor EQ 1) THEN TV, CONGRID(imdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE TV, imdisp 
		CRISPEX_DRAW_SUBCOLOR, event, 0, subcolor, minimum, maximum
		IF (subcolor GE 122) THEN curscolor = 0 ELSE curscolor = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[0])
		CRISPEX_DRAW_SUBCOLOR, event, 0, color_txt, minimum, maximum, xyrange=[10,70,5,20]
		IF (color_txt GE 122) THEN txtcol = 0 ELSE txtcol = 255
		XYOUTS, 10, 10, 'Main image, t='+time,/DEVICE, COLOR=txtcol
	ENDELSE
END

PRO CRISPEX_DRAW_REFLS, event
; (Re)draw reference detailed spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REFLS'
	order_corr = 0.
	WSET, (*(*info).winids).reflswid
	s = (WHERE((*(*info).stokesparams).select_sp EQ 1))[0]
	refspec = ((*(*info).dataparams).refspec)
	refms = (*(*info).dataparams).refms
	title = 'Stokes I'
	ls_low_y = (*(*info).plotaxes).ls_low_y_ref
	ls_upp_y = (*(*info).plotaxes).ls_upp_y_ref
	order_corr=0.
	IF (*(*info).dispswitch).ref_detspect_scale THEN reflsytitle = 'Scaled '+STRLOWCASE((*(*info).plottitles).reflsytitle) ELSE BEGIN
		IF ((FLOOR(ALOG10(ls_low_y)) LE -2) OR (FLOOR(ALOG10(ls_upp_y)) GE 3)) THEN BEGIN
			order_corr = FLOOR(ALOG10(ls_upp_y))
			reflsytitle = (*(*info).plottitles).reflsytitle+' (x10!U'+STRTRIM(order_corr,2)+'!N)'
		ENDIF ELSE reflsytitle = (*(*info).plottitles).reflsytitle
		ls_low_y /= (10.^(order_corr))
		ls_upp_y /= (10.^(order_corr))
	ENDELSE
	IF ((*(*info).dispswitch).ref_detspect_scale EQ 0) THEN BEGIN
		refspec *= refms * ((*(*info).paramparams).scale_cubes)[1] / (10.^(order_corr))
		refms = 1
	ENDIF
	PLOT, (*(*info).dataparams).reflps, refspec, /NORM, CHARSIZE=1, YS=1, YR=[ls_low_y,ls_upp_y], XR=[((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_low],$
		((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_upp]], $
		XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, XTITLE = (*(*info).plottitles).refspxtitle, YTITLE = reflsytitle, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
		POSITION = [(*(*info).plotpos).reflsx0,(*(*info).plotpos).reflsy0,(*(*info).plotpos).reflsx1,(*(*info).plotpos).reflsy1], XTICKLEN = (*(*info).plotaxes).reflsxticklen, YTICKLEN = (*(*info).plotaxes).reflsyticklen
	XYOUTS, (((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_upp]-((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_low])*0.1+((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_low], $
		(ls_upp_y-ls_low_y)*0.9+ls_low_y, title, COLOR = (*(*info).plotparams).plotcol
	IF ( (*(*info).plotswitch).v_dop_set_ref EQ 1) THEN BEGIN
		AXIS, XAXIS=1, XRANGE = [((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_low], ((*(*info).plotaxes).v_dop_ref)[(*(*info).dispparams).lp_ref_upp]], XSTYLE=1, XTITLE = 'Doppler velocity [km/s]', $
		COLOR = (*(*info).plotparams).plotcol
	ENDIF
	IF ((*(*info).dataswitch).refspfile EQ 1) THEN refssp = ( ( *(*(*info).data).refspdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx + FIX((*(*info).dataparams).x) ] )[*,(*(*info).dataparams).t]/refms $
		ELSE refssp = (*(*(*info).data).refsspscan)[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y),(s * (*(*info).dataparams).refnlp):((s+1) * (*(*info).dataparams).refnlp - 1)]/refms
	IF ((*(*info).dispswitch).ref_detspect_scale EQ 0) THEN refssp = refssp * ((*(*info).paramparams).scale_cubes)[1] / (10.^(order_corr))
	OPLOT, (*(*info).dataparams).reflps, refssp, LINE=2, COLOR = (*(*info).plotparams).plotcol
	IF (*(*info).plotswitch).ref_subtract THEN BEGIN
		OPLOT, (*(*info).dataparams).reflps, refspec-refssp, COLOR = (*(*info).plotparams).plotcol
		OPLOT, (*(*info).dataparams).reflps, refspec-refssp, PSYM = 4, COLOR = (*(*info).plotparams).plotcol
	ENDIF
	PLOTS, [1,1] * (*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref],[ls_low_y,ls_upp_y], COLOR = (*(*info).plotparams).plotcol
	IF ((ls_low_y LT 0.) AND (ls_upp_y GT 0.)) THEN PLOTS, [((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_low],((*(*info).dataparams).reflps)[(*(*info).dispparams).lp_ref_upp]], [0.,0.], $
		COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).reflswid,ls_low_y,ls_upp_y], labels=['Window ID for draw','Lower y-value','Upper y-value']
END

PRO CRISPEX_DRAW_REFSP, event	
; (Re)draw reference temporal spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REFSP'
	WSET,(*(*info).winids).refspwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refspwid], labels=['Window ID for draw']
	refspslice = ((*(*(*info).data).refspdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx + FIX((*(*info).dataparams).x) ])[(*(*info).dispparams).lp_ref_low:(*(*info).dispparams).lp_ref_upp,$
		(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
	IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_REFSLICE_SCALING,event,refdisp,refmin,refmax ELSE refmin = MIN(refspslice,MAX=refmax)
	IF (*(*info).dispswitch).warprefspslice THEN dispslice = WARP_TRI((*(*info).dispparams).xi_ref, (*(*info).dispparams).yi_ref, (*(*info).dispparams).xo_ref, (*(*info).dispparams).yo_ref, refspslice) $
		ELSE dispslice = refspslice
       	TV,(CONGRID( BYTSCL(dispslice,MIN=refmin,MAX=refmax), (*(*info).dispparams).refnlpreb, (*(*info).dispparams).refntreb, INTERP = (*(*info).dispparams).interpspslice, /CENTER) ), (*(*info).plotpos).refspx0, $
		(*(*info).plotpos).refspy0, /NORM
	PLOTS, [(*(*info).plotpos).refspx0, (*(*info).plotpos).refspx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).refyplspw + $
		(*(*info).plotpos).refspy0), /NORMAL, COLOR = 100
	PLOTS, [1,1] * ( FLOAT((*(*info).dataparams).lp_ref - (*(*info).dispparams).lp_ref_low) / FLOAT((*(*info).dispparams).lp_ref_range-1) * (*(*info).plotpos).refxplspw + $
		(*(*info).plotpos).refspx0 ), [(*(*info).plotpos).refspy0,(*(*info).plotpos).refspy1], /NORMAL, COLOR = 100
END

PRO CRISPEX_DRAW_INT, event
; (Re)draw intensity versus time plot procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_INT'
	WSET, (*(*info).winids).intwid
	title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]
	int_low_y = (*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] 
	int_upp_y = (*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]
	condition = WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)
	PLOT, FINDGEN((*(*info).dataparams).nt)*(*(*info).plotaxes).dt, FINDGEN((*(*info).dataparams).nt), /NODATA, /NORM, CHARSIZE=1, YR=[int_low_y, int_upp_y], XR = [(*(*info).plotaxes).int_low_t*(*(*info).plotaxes).dt,$
		(*(*info).plotaxes).int_upp_t*(*(*info).plotaxes).dt], /YS, /XS, XTITLE = (*(*info).plottitles).spytitle, YTITLE = 'Counts/Mean Counts', BACKGROUND = (*(*info).plotparams).bgplotcol, $
		COLOR = (*(*info).plotparams).plotcol, LINESTYLE = 0, POSITION = [(*(*info).plotpos).intx0,(*(*info).plotpos).inty0,(*(*info).plotpos).intx1,(*(*info).plotpos).inty1], XTICKLEN = (*(*info).plotaxes).intxticklen, $
		YTICKLEN = (*(*info).plotaxes).intyticklen
	IF (condition[0] NE -1) THEN BEGIN
		selcol = (*(*(*info).intparams).selcol_diagnostics)[condition]
		FOR i=0,N_ELEMENTS(condition)-1 DO BEGIN
			IF (*(*info).dataswitch).spfile THEN BEGIN
				ssp = REFORM( ( ( *(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + $
					(*(*info).dataparams).s ] )[condition[i],*] )
			ENDIF ELSE BEGIN
				ssp = FLTARR((*(*info).dataparams).nt)
				FOR t=0,(*(*info).dataparams).nt-1 DO BEGIN
					ssp[t] = ((*(*(*info).data).imagedata)[t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + (*(*info).dataparams).s * (*(*info).dataparams).nlp + $
						(*(*info).dataparams).lp])[(*(*info).dataparams).x,(*(*info).dataparams).y]
				ENDFOR
			ENDELSE
			avgdint = ssp / ABS(MEAN(ssp))
			LOADCT, 12, /SILENT
			plotcol = ((*(*info).intparams).colors_diagnostics)[selcol[i]]
			OPLOT, FINDGEN((*(*info).dataparams).nt)*(*(*info).plotaxes).dt, avgdint, COLOR=plotcol, LINESTYLE = (*(*(*info).intparams).lines_diagnostics)[condition[i]]
			LOADCT, 0, /SILENT
		ENDFOR
	ENDIF
	t_range = (*(*info).plotaxes).int_upp_t - (*(*info).plotaxes).int_low_t + 1
	XYOUTS, (t_range*(*(*info).plotaxes).dt)*0.1+(*(*info).plotaxes).int_low_t*(*(*info).plotaxes).dt, (int_upp_y-int_low_y)*0.9+int_low_y, title, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).dataparams).t GE (*(*info).plotaxes).int_low_t) AND ((*(*info).dataparams).t LE (*(*info).plotaxes).int_upp_t)) THEN PLOTS, [1,1] * (*(*info).dataparams).t * (*(*info).plotaxes).dt, $
		[int_upp_y, int_low_y], COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).intwid,int_low_y,int_upp_y], labels=['Window ID for draw','Lower y-value','Upper y-value']
END

PRO CRISPEX_DRAW_LS, event
; (Re)draw detailed spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LS'
	WSET, (*(*info).winids).lswid
	IF (*(*info).plotswitch).stokesfile THEN title = 'Stokes '+((*(*info).stokesparams).labels)[WHERE((*(*info).stokesparams).select_sp EQ 1)] ELSE title = ''
	FOR i=0,TOTAL((*(*info).stokesparams).select_sp)-1 DO BEGIN
		s = (WHERE((*(*info).stokesparams).select_sp EQ 1))[i]
		spec = ((*(*info).dataparams).spec)[*,s]
		IF (*(*info).plotswitch).scalestokes THEN ms = (*(*info).dataparams).ms ELSE ms = ((*(*info).dataparams).ms)[s]
		title = 'Stokes '+((*(*info).stokesparams).labels)[s]
		ls_low_y = (*(*(*info).plotaxes).ls_low_y)[s] 
		ls_upp_y = (*(*(*info).plotaxes).ls_upp_y)[s]
		order_corr=0.
		IF (*(*info).dispswitch).detspect_scale THEN lsytitle = 'Scaled '+STRLOWCASE((*(*info).plottitles).lsytitle) ELSE BEGIN
			IF ((FLOOR(ALOG10(ls_low_y)) LE -2) OR (FLOOR(ALOG10(ls_upp_y)) GE 3)) THEN BEGIN
				order_corr = FLOOR(ALOG10(ls_upp_y))
				lsytitle = (*(*info).plottitles).lsytitle+' (x10!U'+STRTRIM(order_corr,2)+'!N)'
			ENDIF ELSE lsytitle = (*(*info).plottitles).lsytitle
			ls_low_y /= (10.^(order_corr))
			ls_upp_y /= (10.^(order_corr))
		ENDELSE
		IF ((*(*info).dispswitch).detspect_scale EQ 0) THEN BEGIN
			spec *= ms * ((*(*info).paramparams).scale_cubes)[0] / (10.^(order_corr))
			ms = 1
		ENDIF
		PLOT, (*(*info).dataparams).lps, spec, /NORM, CHARSIZE=1, YS=1, YR=[ls_low_y,ls_upp_y], XR=[(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], $
			XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, XTITLE = (*(*info).plottitles).spxtitle, YTITLE = lsytitle, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
			POSITION = [((*(*info).plotpos).lsx0)[i],((*(*info).plotpos).lsy0)[i],((*(*info).plotpos).lsx1)[i],((*(*info).plotpos).lsy1)[i]], XTICKLEN = (*(*info).plotaxes).lsxticklen, $
			YTICKLEN = (*(*info).plotaxes).lsyticklen, NOERASE=(i GT 0)
		XYOUTS, ((*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]-(*(*info).dataparams).lps[(*(*info).dispparams).lp_low])*0.1+(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], $
			(ls_upp_y-ls_low_y)*0.9+ls_low_y, title, COLOR = (*(*info).plotparams).plotcol
		IF ( (*(*info).plotswitch).v_dop_set EQ 1) THEN BEGIN
			AXIS, XAXIS=1, XRANGE = [((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_low], ((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_upp]], XSTYLE=1, XTITLE = 'Doppler velocity [km/s]', $
				COLOR = (*(*info).plotparams).plotcol
		ENDIF
		IF ((*(*info).dataswitch).spfile EQ 1) THEN BEGIN
			ssp = ( ( *(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + $
				FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + s ] )[*,(*(*info).dataparams).t]/ms
		ENDIF ELSE BEGIN
			IF (*(*info).dataswitch).onecube THEN $
				ssp = (*(*(*info).data).sspscan)[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y),(s * (*(*info).dataparams).nlp):((s+1) * (*(*info).dataparams).nlp - 1)]/ms ELSE $
				ssp = (*(*(*info).data).scan)[FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y),(s * (*(*info).dataparams).nlp):((s+1) * (*(*info).dataparams).nlp -1)]/ms
		ENDELSE			
		IF ((*(*info).dispswitch).detspect_scale EQ 0) THEN ssp = ssp * ((*(*info).paramparams).scale_cubes)[0] / (10.^(order_corr))
		OPLOT, (*(*info).dataparams).lps, ssp, LINE=2, COLOR = (*(*info).plotparams).plotcol
		IF (*(*info).plotswitch).subtract THEN BEGIN
			OPLOT, (*(*info).dataparams).lps, spec-ssp, COLOR = (*(*info).plotparams).plotcol
			OPLOT, (*(*info).dataparams).lps, spec-ssp, PSYM = 4, COLOR = (*(*info).plotparams).plotcol
		ENDIF
		PLOTS, [1,1] * (*(*info).dataparams).lps[(*(*info).dataparams).lp],[ls_low_y,ls_upp_y], COLOR = (*(*info).plotparams).plotcol
		IF (*(*info).dispswitch).drawdop THEN PLOTS, [1,1] * (*(*info).dataparams).lps[(*(*info).dataparams).lp_dop],[ls_low_y,ls_upp_y], COLOR = (*(*info).plotparams).plotcol
		IF ((ls_low_y LT 0.) AND (ls_upp_y GT 0.)) THEN PLOTS, [((*(*info).dataparams).lps)[(*(*info).dispparams).lp_low],((*(*info).dataparams).lps)[(*(*info).dispparams).lp_upp]], [0.,0.], $
			COLOR = (*(*info).plotparams).plotcol
		IF ((*(*info).winswitch).showref AND ((*(*info).ctrlsswitch).lp_ref_lock EQ 0) AND ((*(*info).dataswitch).refspfile EQ 0) AND ((*(*info).dataparams).refnlp GT 1)) THEN BEGIN
			PLOTS, [1,1] * (*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref],[ls_low_y,ls_upp_y], COLOR = (*(*info).plotparams).plotcol
			IF ((ls_low_y LT 0.) AND (ls_upp_y GT 0.)) THEN PLOTS, [((*(*info).dataparams).lps)[(*(*info).dispparams).lp_low],((*(*info).dataparams).lps)[(*(*info).dispparams).lp_upp]], [0.,0.], $
				COLOR = (*(*info).plotparams).plotcol
		ENDIF
	ENDFOR
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).lswid,ls_low_y,ls_upp_y], labels=['Window ID for draw','Lower y-value','Upper y-value']
END
	
PRO CRISPEX_DRAW_SP, event	
; (Re)draw temporal spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_SP'
	WSET,(*(*info).winids).spwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).spwid], labels=['Window ID for draw']
	spslice = ((*(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns +$
		(*(*info).dataparams).s ])[(*(*info).dispparams).lp_low:(*(*info).dispparams).lp_upp,(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
	IF (*(*info).dispswitch).warpspslice THEN dispslice = WARP_TRI((*(*info).dispparams).xi, (*(*info).dispparams).yi, (*(*info).dispparams).xo, (*(*info).dispparams).yo, spslice) ELSE dispslice = spslice
	IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_XYSLICE_SCALING,event,imdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
       	TV,(CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum), (*(*info).dispparams).nlpreb, (*(*info).dispparams).ntreb, INTERP = (*(*info).dispparams).interpspslice, /CENTER) ), (*(*info).plotpos).spx0, $
		(*(*info).plotpos).spy0, /NORM
	PLOTS, [(*(*info).plotpos).spx0, (*(*info).plotpos).spx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).yplspw + (*(*info).plotpos).spy0), $
		/NORMAL, COLOR = 100
	PLOTS, [1,1] * ( FLOAT((*(*info).dataparams).lp-(*(*info).dispparams).lp_low) / FLOAT((*(*info).dispparams).lp_range-1) * (*(*info).plotpos).xplspw + (*(*info).plotpos).spx0 ), $
		[(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
	IF (*(*info).dispswitch).drawdop THEN PLOTS, [1,1] * ( FLOAT((*(*info).dataparams).lp_dop-(*(*info).dispparams).lp_low) / FLOAT((*(*info).dispparams).lp_range-1) * (*(*info).plotpos).xplspw + (*(*info).plotpos).spx0 ), $
		[(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
	IF ((*(*info).winswitch).showref AND ((*(*info).ctrlsswitch).lp_ref_lock EQ 0) AND ((*(*info).dataswitch).refspfile EQ 0) AND ((*(*info).dataparams).refnlp GT 1)) THEN $
		PLOTS, [1,1] * ( FLOAT((*(*info).dataparams).lp_ref-(*(*info).dispparams).lp_low) / FLOAT((*(*info).dispparams).lp_range-1) * (*(*info).plotpos).xplspw + (*(*info).plotpos).spx0 ), $
			[(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
END

PRO CRISPEX_DRAW_PHIS, event
; (Re)draw spectral slice along a slit procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_PHIS'
	WSET, (*(*info).winids).phiswid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).phiswid], labels=['Window ID for draw']
	IF ((*(*info).curs).lockset EQ 0) THEN BEGIN		
		IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
		IF (*(*info).plotswitch).stokesfile THEN title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle ELSE title = ''
		PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).phiparams).nw_cur), /NODATA, YRANGE = [-(((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + (*(*info).phiparams).sphi)-0.5, $
			(*(*info).phiparams).nw_cur-(((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + (*(*info).phiparams).sphi)-0.5], /YSTYLE, $
			POS = [(*(*info).plotpos).phisx0,(*(*info).plotpos).phisy0,(*(*info).plotpos).phisx1,(*(*info).plotpos).phisy1], YTICKLEN = (*(*info).plotaxes).phisxticklen, XTICKLEN = (*(*info).plotaxes).phisyticklen, $
			XTITLE = (*(*info).plottitles).spxtitle, XR = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low],(*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], YTITLE = 'Position along slit [pixel]', $
			XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, XRANGE = [((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_low], ((*(*info).plotaxes).v_dop)[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTITLE = title+'Doppler velocity [km/s]', COLOR = (*(*info).plotparams).plotcol ELSE $
			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, XRANGE = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]], XSTYLE=1, $
				XTICKNAME = REPLICATE(' ',60), XTITLE = title, COLOR = (*(*info).plotparams).plotcol
	ENDIF
	dispslice = (*(*(*info).data).phislice)[(*(*info).dispparams).lp_low:(*(*info).dispparams).lp_upp,*]
	IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_XYSLICE_SCALING,event,imdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
	TV, CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum), (*(*info).dispparams).phisnlpreb, (*(*info).dispparams).nphireb, $
		INTERP = (*(*info).dispparams).interpspslice), (*(*info).plotpos).phisx0, (*(*info).plotpos).phisy0, /NORM
	PLOTS, [(*(*info).plotpos).phisx0, (*(*info).plotpos).phisx1], [1,1] * ( (((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + $
		(*(*info).phiparams).sphi+0.5)/FLOAT((*(*info).phiparams).nw_cur)*(*(*info).plotpos).phisyplspw + (*(*info).plotpos).phisy0), /NORMAL, COLOR = 100
	PLOTS, [1,1] * ( FLOAT((*(*info).dataparams).lp-(*(*info).dispparams).lp_low) / FLOAT((*(*info).dispparams).lp_range-1) * (*(*info).plotpos).phisxplspw + (*(*info).plotpos).phisx0 ), $
		[(*(*info).plotpos).phisy0, (*(*info).plotpos).phisy1], /NORMAL, COLOR = 100
END

PRO CRISPEX_DRAW_LOOPSLAB, event
; (Re)draw loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_LOOPSLAB'
	IF ((SIZE(*(*(*info).loopsdata).crossloc))[1] NE (*(*info).loopparams).np) THEN BEGIN
		CRISPEX_LOOP_GET, event
		CRISPEX_UPDATE_LP, event
		IF ((SIZE(*(*(*info).loopsdata).crossloc))[1] NE (*(*info).loopparams).np) THEN RETURN
	ENDIF
	WSET, (*(*info).winids).loopwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).loopwid], labels=['Window ID for draw']
	dispslice = (*(*(*info).loopsdata).loopslice)[*,(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
	IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_XYSLICE_SCALING,event,imdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
	TV, CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum), (*(*info).dispparams).loopnlxreb, (*(*info).dispparams).loopntreb, /INTERP), $
		(*(*info).plotpos).loopx0, (*(*info).plotpos).loopy0,/NORM
	PLOTS, [(*(*info).plotpos).loopx0, (*(*info).plotpos).loopx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).loopyplspw + $
		(*(*info).plotpos).loopy0), /NORMAL, COLOR = 100
	FOR i=0,(*(*info).loopparams).np-1 DO BEGIN
		PLOTS, [1,1] * ( FLOAT((*(*(*info).loopsdata).crossloc)[i]) / FLOAT((*(*info).loopsdata).loopsize-1) * (*(*info).plotpos).loopxplspw + (*(*info).plotpos).loopx0 ), $
			[(*(*info).plotpos).loopy0, (*(*info).plotpos).loopy1], /NORMAL, COLOR = 100
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, (*(*(*info).loopsdata).crossloc)[i], labels='crossloc['+STRTRIM(i,2)+']'
	ENDFOR
END

PRO CRISPEX_DRAW_REFLOOPSLAB, event
; (Re)draw reference loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REFLOOPSLAB'
	WSET, (*(*info).winids).refloopwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refloopwid], labels=['Window ID for draw']
	dispslice = (*(*(*info).loopsdata).refloopslice)[*,(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
	IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_REFSLICE_SCALING,event,refdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
	TV, CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum), (*(*info).dispparams).refloopnlxreb, (*(*info).dispparams).refloopntreb, /INTERP), $
		(*(*info).plotpos).refloopx0, (*(*info).plotpos).refloopy0,/NORM
	PLOTS, [(*(*info).plotpos).refloopx0, (*(*info).plotpos).refloopx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).refloopyplspw + $
		(*(*info).plotpos).refloopy0), /NORMAL, COLOR = 100
	FOR i=0,(*(*info).loopparams).np-1 DO BEGIN
		PLOTS, [1,1] * ( FLOAT((*(*(*info).loopsdata).crossloc)[i]) / FLOAT((*(*info).loopsdata).loopsize-1) * (*(*info).plotpos).refloopxplspw + (*(*info).plotpos).refloopx0 ), $
			[(*(*info).plotpos).refloopy0, (*(*info).plotpos).refloopy1], /NORMAL, COLOR = 100
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, (*(*(*info).loopsdata).crossloc)[i], labels='crossloc['+STRTRIM(i,2)+']'
	ENDFOR
END

PRO CRISPEX_DRAW_REST_LOOP, event
; (Re)draw restored loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_REST_LOOP'
	FOR k=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO BEGIN
		WSET, (*(*(*info).winids).restloopwid)[k]
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).winids).restloopwid)[k]], labels=['Window ID for draw']
		IF (*(*(*info).dispswitch).restricted_t_range)[k] THEN BEGIN
			lower_t = 0
			upper_t = (*(*info).dispparams).t_upp - (*(*info).dispparams).t_low
		ENDIF ELSE BEGIN
			lower_t = (*(*info).dispparams).t_low
			upper_t = (*(*info).dispparams).t_upp
		ENDELSE
		dispslice = (*(*(*(*info).loopsdata).rest_loopslice[k]))[*,lower_t:upper_t]
		IF (*(*(*info).restoreparams).disp_imref)[k] THEN BEGIN
			IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_REFSLICE_SCALING,event,refdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
		ENDIF ELSE BEGIN
			IF (*(*info).dispparams).slices_imscale THEN CRISPEX_DRAW_XYSLICE_SCALING,event,imdisp,minimum,maximum ELSE minimum = MIN(dispslice,MAX=maximum)
		ENDELSE
		TV, CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum), (*(*info).dispparams).restloopnlxreb, (*(*info).dispparams).restloopntreb, /INTERP), $
			(*(*info).plotpos).restloopx0, (*(*info).plotpos).restloopy0, /NORM
		PLOTS, [(*(*info).plotpos).restloopx0, (*(*info).plotpos).restloopx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * $
			(*(*info).plotpos).restloopyplspw + (*(*info).plotpos).restloopy0), /NORMAL, COLOR = 100
		FOR i=0,(SIZE(*(*(*(*info).loopsdata).rest_crossloc[k])))[1]-1 DO BEGIN
			PLOTS, [1,1] * ( FLOAT((*(*(*(*info).loopsdata).rest_crossloc[k]))[i]) / FLOAT((*(*(*info).loopsdata).rest_loopsize)[k]-1) * (*(*info).plotpos).restloopxplspw + (*(*info).plotpos).restloopx0 ),$
				[(*(*info).plotpos).restloopy0, (*(*info).plotpos).restloopy1], /NORMAL, COLOR = 100
		ENDFOR
	ENDFOR
END

PRO CRISPEX_DRAW_RETR_DET, event
; (Re)draw restreived detection procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_RETR_DET'
	WSET, (*(*info).winids).retrdetwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).retrdetwid], labels=['Window ID for draw']
	TVSCL, CONGRID((*(*(*info).loopsdata).det_loopslice)[*,(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp], (*(*info).dispparams).retrdetnlxreb, (*(*info).dispparams).retrdetntreb, /INTERP), $
		(*(*info).plotpos).retrdetx0, (*(*info).plotpos).retrdety0, /NORM
	PLOTS, [(*(*info).plotpos).retrdetx0, (*(*info).plotpos).retrdetx1], [1,1] * ( ((*(*info).dataparams).t-(*(*info).dispparams).t_low) / FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).retrdetyplspw + $
		(*(*info).plotpos).retrdety0), /NORMAL, COLOR = 100
END

;================================================================================= ESTIMATE SAVING TIME PROCEDURES
PRO CRISPEX_ESTIMATE_TIME_WINDOW, event
; Opens the calculating saving time estimate window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_ESTIMATE_TIME_WINDOW'
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': CALCULATING...', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	(*(*info).ctrlsfeedb).estimate_label = WIDGET_LABEL(message_base, VALUE = 'Calculating time required for save procedure. Please wait...', /ALIGN_LEFT)
	label2	= WIDGET_LABEL(message_base, VALUE = ' ')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).estimatetlb = base
	(*(*info).winswitch).estimate_win = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).estimatetlb], labels=['estimatetlb']
END

PRO CRISPEX_ESTIMATE_TIME_CALCULATION, event
; Performs the actual saving time estimate calculation
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_DRAW_ESTIMATE_TIME_CALCULATION'
	IF (*(*info).loopswitch).retrieve_loops THEN BEGIN
		RESTORE,(*(*(*info).retrparams).retrieve_files)[0]
		*(*(*info).loopparams).xr = x_loop_pts		&	*(*(*info).loopparams).yr = y_loop_pts
		*(*(*info).loopparams).w_lpts = w_loop_pts	&	(*(*info).dataparams).lp = spect_pos
	ENDIF 
	IF ((*(*info).dataparams).nt LT 100) THEN BEGIN
		endt = (*(*info).dataparams).nt & modulus = FIX((*(*info).dataparams).nt/3.) & result2 = FIX((*(*info).dataparams).nt/9.) & result3 = FIX((*(*info).dataparams).nt/9.*2.)
	ENDIF ELSE BEGIN
		endt = 100 & modulus = 30 & result2 = 10 & result3 = 20
	ENDELSE
	t_0 = SYSTIME(/SECONDS)
	FOR t=0,endt-1 DO BEGIN
			IF (t EQ 0) THEN tmp = INTERPOLATE( (*(*(*info).data).imagedata)[t*(*(*info).dataparams).nlp + (*(*info).dataparams).lp], (*(*(*info).loopparams).xr)[*(*(*info).loopparams).w_lpts],$
				(*(*(*info).loopparams).yr)[*(*(*info).loopparams).w_lpts]) $
			ELSE tmp = [[tmp], [INTERPOLATE( (*(*(*info).data).imagedata)[t*(*(*info).dataparams).nlp + (*(*info).dataparams).lp], (*(*(*info).loopparams).xr)[*(*(*info).loopparams).w_lpts],$
				(*(*(*info).loopparams).yr)[*(*(*info).loopparams).w_lpts])]]
		IF (t mod modulus EQ 0) THEN WIDGET_CONTROL, (*(*info).ctrlsfeedb).estimate_label, SET_VALUE = 'Calculating time required for save procedure. Please wait.'
		IF (t mod modulus EQ result2) THEN WIDGET_CONTROL, (*(*info).ctrlsfeedb).estimate_label, SET_VALUE = 'Calculating time required for save procedure. Please wait..'
		IF (t mod modulus EQ result3) THEN WIDGET_CONTROL, (*(*info).ctrlsfeedb).estimate_label, SET_VALUE = 'Calculating time required for save procedure. Please wait...'
	ENDFOR
	(*(*info).feedbparams).estimate_lx = N_ELEMENTS(*(*(*info).loopparams).w_lpts)
	(*(*info).feedbparams).estimate_time = (SYSTIME(/SECONDS) - t_0) / FLOAT(endt)
	(*(*info).feedbparams).estimate_run = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_estimate, /SENSITIVE
	IF (*(*info).loopswitch).retrieve_loops THEN BEGIN
		*(*(*info).loopparams).xr = 0	&	*(*(*info).loopparams).yr = 0	&	*(*(*info).loopparams).w_lpts = 0
	ENDIF
END

PRO CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
; Gets the units and denominator for the saving time estimate calculation
	IF (time LT 60) THEN BEGIN
		units = ' seconds'
		denom = 1. 
	ENDIF ELSE IF (time LT 3600) THEN BEGIN
		units = ' minutes' 
		denom = 60. 
	ENDIF ELSE IF (time LT 86400) THEN BEGIN
		units = ' hours' 
		denom = 3600.
	ENDIF ELSE BEGIN
		units = ' days'
		denom = 86400.
	ENDELSE
END

PRO CRISPEX_ESTIMATE_FULL_TIME_RUNNING, pass, totalpasses, t0, t1, denom, units, accumsectime, totalsectime
; Gets the units and denominator for the saving time running estimate
	accumsectime = t1-t0
	totalsectime = accumsectime/pass*totalpasses
	IF (totalsectime GT 60) THEN BEGIN
		IF (totalsectime GT 3600) THEN BEGIN
			IF (totalsectime GT 86400) THEN BEGIN
				units = ' days.'
				denom = 86400.
			ENDIF ELSE BEGIN
				units = ' hours.'
				denom = 3600.
			ENDELSE
		ENDIF ELSE BEGIN
			units = ' minutes.'
			denom = 60.
		ENDELSE
	ENDIF ELSE BEGIN
		units = ' seconds.'
		denom = 1.
	ENDELSE
END

;================================================================================= TAB EVENT PROCEDURE
PRO CRISPEX_EVENT, event
; Handles tab events
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_EVENT', /IGNORE_LAST
END

;================================================================================= FIND CRISPEX OUTPUT FILE PROCEDURES
PRO CRISPEX_FIND_CLSAV, event
; Finds CLSAV output files (i.e. saved loop points files)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_FIND_CLSAV'
	imagefilename = STRMID((*(*info).dataparams).imfilename,STRPOS((*(*info).dataparams).imfilename,PATH_SEP(),/REVERSE_SEARCH)+1,STRLEN((*(*info).dataparams).imfilename))
	firstsplit = STRMID(imagefilename,0,STRPOS(imagefilename,'.',/REVERSE_SEARCH))
	fstr = STRSPLIT(firstsplit[0],'_',/EXTRACT)
	IF (N_ELEMENTS(fstr) GE 2) THEN filename = fstr[0]+'_'+fstr[1] ELSE filename = fstr[0]
	clfiles = FILE_SEARCH((*(*info).paths).ipath+filename+"*clsav", COUNT = clfilecount)
	*(*(*info).retrparams).clfiles  = clfiles
	(*(*info).retrparams).clfilecount = clfilecount
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [imagefilename, filename, STRTRIM((*(*info).retrparams).clfilecount,2)], labels=['filename','basename','clfilecount']
END

PRO CRISPEX_FIND_CSAV, event
; Finds CSAV output files (i.e. saved loopslice/slab files)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_FIND_CSAV'
	imagefilename = STRMID((*(*info).dataparams).imfilename,STRPOS((*(*info).dataparams).imfilename,PATH_SEP(),/REVERSE_SEARCH)+1,STRLEN((*(*info).dataparams).imfilename))
	firstsplit = STRMID(imagefilename,0,STRPOS(imagefilename,'.',/REVERSE_SEARCH))
	fstr = STRSPLIT(firstsplit[0],'_',/EXTRACT)
	IF (N_ELEMENTS(fstr) GE 2) THEN filename = fstr[0]+'_'+fstr[1] ELSE filename = fstr[0]
	cfiles = FILE_SEARCH((*(*info).paths).ipath+filename+"*csav", COUNT = cfilecount)
	IF (*(*info).dataswitch).reffile THEN BEGIN
		reffilename = STRMID((*(*info).dataparams).reffilename,STRPOS((*(*info).dataparams).reffilename,PATH_SEP(),/REVERSE_SEARCH)+1,STRLEN((*(*info).dataparams).reffilename))
		reffirstsplit = STRMID(reffilename,0,STRPOS(reffilename,'.',/REVERSE_SEARCH))
		reffstr = STRSPLIT(reffirstsplit[0],'_',/EXTRACT)
		IF (N_ELEMENTS(reffstr) GE 2) THEN reffilename = reffstr[0]+'_'+reffstr[1] ELSE reffilename = reffstr[0]
		refcfiles = FILE_SEARCH((*(*info).paths).ipath+reffilename+"*csav", COUNT = refcfilecount)
		*(*(*info).restoreparams).cfiles  = [cfiles,refcfiles]
		(*(*info).restoreparams).cfilecount = cfilecount+refcfilecount
	ENDIF ELSE BEGIN
		*(*(*info).restoreparams).cfiles  = cfiles
		(*(*info).restoreparams).cfilecount = cfilecount
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [imagefilename, filename, STRTRIM((*(*info).restoreparams).cfilecount,2)], labels=['filename','basename','cfilecount']
END

;================================================================================= HELP PROCEDURE
PRO CRISPEX_HELP, event
; Opens the CRISPEX Reference Pages
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_HELP'
	tempfile = FILEPATH('temp_crispex_redirect.html', /TMP)
	OPENW, lun, tempfile, /GET_LUN
	PRINTF, lun, '<HTML><HEADER><TITLE>CRISPEX help pages</TITLE><META HTTP-EQUIV="REFRESH" CONTENT="0;URL=http://folk.uio.no/gregal/crispex">'+$
		'</HEADER><BODY></BODY></HTML>'
	FREE_LUN, lun
	ONLINE_HELP, BOOK=tempfile
	WAIT, 10.0
	FILE_DELETE, tempfile
END

PRO CRISPEX_MAIL_BUG, event
; Opens a new message for bug reporting
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MAIL_BUG'
	tempfile = FILEPATH('temp_crispex_report_bug.html', /TMP)
	subject = 'CRISPEX:%20Reporting%20bug%20in%20v'+(*(*info).versioninfo).version_number+'%20(rev%20'+(*(*info).versioninfo).revision_number+')'
	OPENW, lun, tempfile, /GET_LUN
	PRINTF, lun, '<HTML><HEADER><TITLE>CRISPEX: Report a bug</TITLE><META HTTP-EQUIV="REFRESH" CONTENT="0;URL=mailto:g.j.m.vissers@astro.uio.no?SUBJECT='+subject+'">'+$
		'</HEADER><BODY>If your e-mail client does not automatically open up a new message window, you may also send your bug report manually to: g.j.m.vissers@astro.uio.no, '+$
		'preferrably with <i>"CRISPEX: Reporting bug in v'+(*(*info).versioninfo).version_number+' (rev '+(*(*info).versioninfo).revision_number+')"</i> as subject heading.</BODY></HTML>'
	FREE_LUN, lun
	ONLINE_HELP, BOOK=tempfile
	WAIT, 10.0
	FILE_DELETE, tempfile
END

PRO CRISPEX_MAIL_SUGGESTION, event
; Opens a new message for suggestion reporting
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MAIL_SUGGESTION'
	tempfile = FILEPATH('temp_crispex_report_suggestion.html', /TMP)
	subject = 'CRISPEX:%20Suggestion%20after%20v'+(*(*info).versioninfo).version_number+'%20(rev%20'+(*(*info).versioninfo).revision_number+')'
	OPENW, lun, tempfile, /GET_LUN
	PRINTF, lun, '<HTML><HEADER><TITLE>CRISPEX: Send a suggestion</TITLE><META HTTP-EQUIV="REFRESH" CONTENT="0;URL=mailto:g.j.m.vissers@astro.uio.no?SUBJECT='+subject+'">'+$
		'</HEADER><BODY>If your e-mail client does not automatically open up a new message window, you may also send your suggestion manually to: g.j.m.vissers@astro.uio.no, '+$
		'preferrably with <i>"CRISPEX: Suggestion after v'+(*(*info).versioninfo).version_number+' (rev '+(*(*info).versioninfo).revision_number+')"</i> as subject heading.</BODY></HTML>'
	FREE_LUN, lun
	ONLINE_HELP, BOOK=tempfile
	WAIT, 10.0
	FILE_DELETE, tempfile
END

;================================================================================= INTENSITY-TIME SAVE PROCEDURES
PRO CRISPEX_INT_SAVE, event
; Handles the actual saving of the intensity-time plots
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_INT_SAVE'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=intfilename, /tlab, ext='cint'
	condition = WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)
	intensities = FLTARR((*(*info).dispparams).t_range,N_ELEMENTS(condition))
	avg_intensity = FLTARR(N_ELEMENTS(condition))
	FOR i=0,N_ELEMENTS(condition)-1 DO BEGIN
		intensities[0,i] = REFORM( ( ( *(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + $
			(*(*info).dataparams).s ] )[condition[i],(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp] )
		avg_intensity[i] = MEAN(intensities[*,i])
	ENDFOR
	diagnostics = ((*(*info).intparams).diagnostics)[condition]
	x = (*(*info).dataparams).x		&	y = (*(*info).dataparams).y
	nt = (*(*info).dispparams).t_range	& 	dt = (*(*info).plotaxes).dt
	t_low = (*(*info).dispparams).t_low	&	t_upp = (*(*info).dispparams).t_upp
	t_saved = (*(*info).dataparams).t	&	crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
	SAVE, crispex_version, intensities, avg_intensity, diagnostics, nt, dt, t_low, t_upp, t_saved, x, y, FILENAME=(*(*info).paths).opath+intfilename
	PRINT, 'Written: '+(*(*info).paths).opath+intfilename+'.cint'
END


;================================================================================= INTERRUPT PROCEDURE
PRO CRISPEX_INTERRUPT, event
; Handles interrupting at runtime
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_INTERRUPT', /IGNORE_LAST
	PRINT,'Interrupted CRISPEX at runtime. Type [.c] to continue...'
	STOP
END

;================================================================================= LOOP PROCEDURES
PRO CRISPEX_LOOP_DEFINE, event
; Handles the start of loop definition procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_DEFINE'
	(*(*info).overlayswitch).loopslit = event.SELECT
	IF ((*(*info).overlayswitch).loopslit EQ 0) THEN BEGIN
		IF ((*(*info).winids).looptlb GT 0) THEN WIDGET_CONTROL, (*(*info).winids).looptlb, /DESTROY
		*(*(*info).loopparams).xp = 0.		&	*(*(*info).loopparams).yp = 0.		&	*(*(*info).overlayparams).sxp = 0.	
		*(*(*info).overlayparams).syp = 0.	&	*(*(*info).loopparams).xr = 0.		&	*(*(*info).loopparams).yr = 0.
		*(*(*info).overlayparams).sxr = 0.	&	*(*(*info).overlayparams).syr = 0.	&	(*(*info).loopparams).np = 0
		(*(*info).loopsdata).loopsize = 0	&	(*(*info).winids).looptlb = 0		&	(*(*info).overlayswitch).loopslit = 0
		(*(*info).winswitch).showloop = 0
		IF (*(*info).winswitch).showrefloop THEN BEGIN
			IF ((*(*info).winids).reflooptlb GT 0) THEN WIDGET_CONTROL, (*(*info).winids).reflooptlb, /DESTROY
			(*(*info).winids).reflooptlb = 0	&	(*(*info).winswitch).showrefloop = 0
		ENDIF
		WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).save_loop_pts, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).timeslicemenu, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).loop_slit_but, SET_VALUE = 'Draw loop path'
		WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, /SET_BUTTON
		(*(*info).curs).lockset = event.SELECT
		CRISPEX_COORDSLIDERS_SET, 1, 1, event
		CRISPEX_DRAW, event
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_but, SENSITIVE = ABS(event.SELECT-1)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayswitch).loopslit], labels=['Draw loop path']
END

PRO CRISPEX_LOOP_FEEDBACK, event
; Loop path feedback selection procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_FEEDBACK', /IGNORE_LAST
	(*(*info).overlayswitch).looppath_feedback = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayswitch).looppath_feedback], labels=['Loop path feedback']
END

PRO CRISPEX_LOOP_GET, event
; Gets the loop path and (if the loop display window is open) also the loopslab for display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET'
	CRISPEX_LOOP_GET_PATH, event
	IF (*(*info).winswitch).showloop THEN BEGIN
		CRISPEX_LOOP_GET_SLAB, event
		CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event					
	ENDIF
END

PRO CRISPEX_LOOP_GET_PATH, event
; Gets the actual loop path from spline interpolation
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET_PATH'
	np_local = (SIZE(*(*(*info).loopparams).xp))[2] 
	IF (np_local GE 2) THEN BEGIN
		IF ((*(*(*info).loopparams).xp)[np_local-1] EQ (*(*(*info).loopparams).xp)[np_local-2]) AND ((*(*(*info).loopparams).yp)[np_local-1] EQ (*(*(*info).loopparams).yp)[np_local-2]) THEN RETURN
	ENDIF
	SPLINE_P,*(*(*info).loopparams).xp,*(*(*info).loopparams).yp,xr,yr,INTERVAL=1
	*(*(*info).loopparams).xr = xr	&	*(*(*info).loopparams).yr = yr
	*(*(*info).loopparams).w_lpts = WHERE((*(*(*info).loopparams).xr GE 0) AND (*(*(*info).loopparams).xr LE (*(*info).dispparams).x_last) AND (*(*(*info).loopparams).yr GE 0) AND $
		(*(*(*info).loopparams).yr LE (*(*info).dispparams).y_last), nw_lpts)
	(*(*info).loopparams).nw_lpts = nw_lpts
	*(*(*info).overlayparams).sxr = *(*(*info).loopparams).xr * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
	*(*(*info).overlayparams).syr = *(*(*info).loopparams).yr * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
END

PRO CRISPEX_LOOP_GET_SLAB, event
; Handles the extraction of the loopslab for display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET_SLAB'
	IF (*(*info).dispswitch).exts THEN BEGIN
		WIDGET_CONTROL, /HOURGLASS
		lp_orig = (*(*info).dataparams).lp
		FOR i=(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp DO BEGIN
			(*(*info).dataparams).lp = i
			CRISPEX_LOOP_GET_EXACT_SLICE, event, *(*(*info).data).imagedata, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, $
				*(*(*info).loopparams).w_lpts, loopslice, crossloc, loopsize,/im
			IF (i EQ (*(*info).dispparams).lp_low) THEN loopslab = loopslice ELSE loopslab = [[[loopslab]], [[loopslice]]]
			PRINT,'> '+STRTRIM(i-(*(*info).dispparams).lp_low+1,2)+'/'+STRTRIM((*(*info).dispparams).lp_upp-(*(*info).dispparams).lp_low+1,2)+' slices extracted...'
		ENDFOR
		(*(*info).dataparams).lp = lp_orig
	ENDIF ELSE BEGIN
		CRISPEX_LOOP_GET_APPROX_SLAB, event, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, loopslab, crossloc, loopsize
	ENDELSE
	*(*(*info).loopsdata).loopslab = loopslab
	*(*(*info).loopsdata).crossloc = crossloc
	(*(*info).loopsdata).loopsize = loopsize
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [N_ELEMENTS(crossloc),loopsize], labels=['np','loopsize']
END

PRO CRISPEX_LOOP_GET_REFSLAB, event
; Handles the extraction of the reference loopslab for display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET_REFSLAB'
	IF (*(*info).dispswitch).refexts THEN BEGIN
		WIDGET_CONTROL, /HOURGLASS
		lp_orig = (*(*info).dataparams).lp_ref
		FOR i=(*(*info).dispparams).lp_ref_low,(*(*info).dispparams).lp_ref_upp DO BEGIN
			(*(*info).dataparams).lp_ref = i
			CRISPEX_LOOP_GET_EXACT_SLICE, event, *(*(*info).data).refdata, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, $
				*(*(*info).loopparams).w_lpts, refloopslice, crossloc, loopsize
			IF (i EQ (*(*info).dispparams).lp_ref_low) THEN refloopslab = refloopslice ELSE refloopslab = [[[refloopslab]], [[refloopslice]]]
			PRINT,'> '+STRTRIM(i-(*(*info).dispparams).lp_ref_low+1,2)+'/'+STRTRIM((*(*info).dispparams).lp_ref_upp-(*(*info).dispparams).lp_ref_low+1,2)+' slices extracted...'
		ENDFOR
		(*(*info).dataparams).lp_ref = lp_orig
	ENDIF ELSE BEGIN
		CRISPEX_LOOP_GET_APPROX_SLAB, event, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, refloopslab, crossloc, loopsize, /REFERENCE
	ENDELSE
	*(*(*info).loopsdata).refloopslab = refloopslab
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [N_ELEMENTS(crossloc),loopsize], labels=['np','loopsize']
END

PRO CRISPEX_LOOP_GET_APPROX_SLAB, event, xrs, yrs, xps, yps, ap_loopslab, ap_crossloc, ap_loopsize, REFERENCE=reference
; Gets the approximated loopslab for in-pragram display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET_APPROX_SLAB'
	IF KEYWORD_SET(REFERENCE) THEN BEGIN
		FOR i=0,(SIZE(xrs))[1]-1 DO BEGIN
			IF (i EQ 0) THEN tmp = ((*(*(*info).data).refspdata)[ROUND(yrs[i]) * (*(*info).dataparams).nx + ROUND(xrs[i])]) ELSE $
				tmp = [ [[tmp]] , [[ ((*(*(*info).data).refspdata)[ROUND(yrs[i]) * (*(*info).dataparams).nx + ROUND(xrs[i])]) ]] ]
		ENDFOR
	ENDIF ELSE BEGIN
		FOR i=0,(SIZE(xrs))[1]-1 DO BEGIN
			IF (i EQ 0) THEN BEGIN
				tmp = ((*(*(*info).data).spdata)[ROUND(yrs[i]) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + ROUND(xrs[i]) * (*(*info).dataparams).ns + (*(*info).dataparams).s])
			ENDIF ELSE BEGIN
				tmp = [[[tmp]],[[((*(*(*info).data).spdata)[ROUND(yrs[i]) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + ROUND(xrs[i]) * (*(*info).dataparams).ns + (*(*info).dataparams).s])]]]
			ENDELSE
		ENDFOR
	ENDELSE
	ap_loopslab = TRANSPOSE(tmp, [2,1,0])
	n = 1 + (ABS(((SIZE(xps))[0] EQ 1)-1))
	FOR k=0,(SIZE(xps))[n]-1 DO BEGIN
		IF (k EQ 0) THEN ap_crossloc = [0] ELSE ap_crossloc = [ap_crossloc,WHERE( (xrs EQ (xps[k])) AND (yrs EQ (yps[k])))]
	ENDFOR
	ap_loopsize = (SIZE(ap_loopslab))[1]
END

PRO CRISPEX_LOOP_GET_EXACT_SLICE, event, extractdata, xrs, yrs, xps, yps, w_lpts, ex_loopslice, ex_crossloc, ex_loopsize, NO_NLP=no_nlp, im=im
; Gets the exact loopslice for in-program display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info			
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_GET_EXACT_SLICE'
	IF KEYWORD_SET(NO_NLP) THEN BEGIN
		FOR t=(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp DO BEGIN
			IF (t EQ (*(*info).dispparams).t_low) THEN tmp = INTERPOLATE( extractdata[t], (xrs)[w_lpts],(yrs)[w_lpts]) ELSE tmp = [[tmp], [INTERPOLATE( extractdata[t], (xrs)[w_lpts],(yrs)[w_lpts])]]
		ENDFOR
	ENDIF ELSE BEGIN
		IF KEYWORD_SET(IM) THEN BEGIN		; Get space-time diagram from main cube
			FOR t=(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp DO BEGIN
				IF (t EQ (*(*info).dispparams).t_low) THEN tmp = INTERPOLATE( extractdata[t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + (*(*info).dataparams).s * (*(*info).dataparams).nlp + $
					(*(*info).dataparams).lp], (xrs)[w_lpts],(yrs)[w_lpts]) ELSE $
					tmp = [[tmp], [INTERPOLATE( extractdata[t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + (*(*info).dataparams).s * (*(*info).dataparams).nlp + $
						(*(*info).dataparams).lp], (xrs)[w_lpts],(yrs)[w_lpts])]] 
			ENDFOR
		ENDIF ELSE BEGIN			; Get space-time diagram from reference cube
			FOR t=(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp DO BEGIN
				IF (t EQ (*(*info).dispparams).t_low) THEN tmp = INTERPOLATE( extractdata[t*(*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref], (xrs)[w_lpts],(yrs)[w_lpts]) ELSE $
					tmp = [[tmp], [INTERPOLATE( extractdata[t*(*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref], (xrs)[w_lpts],(yrs)[w_lpts])]]
			ENDFOR
		ENDELSE
	ENDELSE
	ex_loopslice = tmp
	n = 1 + (ABS(((SIZE(xps))[0] EQ 1)-1))
	FOR k=0,(SIZE(xps))[n]-1 DO BEGIN
		IF (k EQ 0) THEN ex_crossloc = [0] ELSE ex_crossloc = [ex_crossloc,WHERE( (xrs EQ (xps)[k]) AND (yrs EQ (yps)[k]))]
	ENDFOR
	ex_loopsize = (SIZE(ex_loopslice))[1]
END


PRO CRISPEX_LOOP_REMOVE_POINT, event
; Handles the removal of a loop point
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_LOOP_REMOVE_POINT'
	(*(*info).loopparams).np -= 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopparams).np], labels=['np']
	*(*(*info).loopparams).xp = (*(*(*info).loopparams).xp)[*,0:(*(*info).loopparams).np-1]
	*(*(*info).loopparams).yp = (*(*(*info).loopparams).yp)[*,0:(*(*info).loopparams).np-1]
	*(*(*info).overlayparams).sxp = (*(*(*info).overlayparams).sxp)[*,0:(*(*info).loopparams).np-1]
	*(*(*info).overlayparams).syp = (*(*(*info).overlayparams).syp)[*,0:(*(*info).loopparams).np-1]
	(*(*info).dataparams).x = (*(*(*info).loopparams).xp)[(*(*info).loopparams).np-1]
	(*(*info).dataparams).y = (*(*(*info).loopparams).yp)[(*(*info).loopparams).np-1]
	(*(*info).curs).sxlock = (*(*(*info).overlayparams).sxp)[(*(*info).loopparams).np-1]
	(*(*info).curs).sylock = (*(*(*info).overlayparams).syp)[(*(*info).loopparams).np-1]
	CRISPEX_COORDSLIDERS_SET, 0, 0, event
	IF ((*(*info).loopparams).np EQ 2) THEN WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = 0
	IF ((*(*info).loopparams).np LT 2) THEN WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 0
	IF ((*(*info).zooming).factor NE 1) THEN CRISPEX_ZOOM_LOOP, event
	IF (*(*info).winswitch).showloop THEN BEGIN
		CRISPEX_UPDATE_LP, event
		WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 0
		IF ((*(*info).loopparams).np EQ 2) THEN WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = 0
		CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event					
	ENDIF ELSE CRISPEX_LOOP_GET, event
	IF (*(*info).winswitch).showphis THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
		CRISPEX_UPDATE_PHIS, event
	ENDIF ELSE CRISPEX_DRAW, event
END
;================================================================================= MASK PROCEDURES
FUNCTION CRISPEX_MASK_OVERLAY, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MASK_OVERLAY'
	maskim = (*(*info).overlayswitch).maskim
	maskim[event.value] = event.select
	(*(*info).overlayswitch).maskim = maskim
	CRISPEX_DRAW, event
END

PRO CRISPEX_MASK_OVERLAY_COLOR_SLIDER, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MASK_OVERLAY_COLOR_SLIDER'
	(*(*info).overlayparams).maskcolor = event.VALUE
	CRISPEX_DRAW, event
END

PRO CRISPEX_MASK_OVERLAY_SELECT_COLOR_TABLE, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MASK_OVERLAY_SELECT_COLOR_TABLE'
	(*(*info).overlayparams).maskct = event.index
	CRISPEX_DRAW, event
END

PRO CRISPEX_MASK_BUTTONS_SET, event
; Handles the setting of mask buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MASK_BUTTONS_SET'
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[0],SET_BUTTON = ((*(*info).overlayswitch).maskim)[0], SENSITIVE = (*(*info).overlayswitch).mask
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[1],SET_BUTTON = ((*(*info).overlayswitch).maskim)[1], SENSITIVE = ((*(*info).overlayswitch).mask AND (*(*info).dataswitch).reffile)
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[2],SET_BUTTON = ((*(*info).overlayswitch).maskim)[2], SENSITIVE = ((*(*info).overlayswitch).mask AND (*(*info).winswitch).showdop)
END

;================================================================================= MEASUREMENT PROCEDURES
PRO CRISPEX_MEASURE_ARCSEC, event
; Handles the change of arcseconds per pixel resolution
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MEASURE_ARCSEC', /IGNORE_LAST
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_text, GET_VALUE = textvalue
	(*(*info).meas).arcsecpix = FLOAT(textvalue[0])
	IF ((*(*info).meas).np GT 0) THEN CRISPEX_MEASURE_CALC, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).meas).arcsecpix], labels=['arcsecpix']
END

PRO CRISPEX_MEASURE_ENABLE, event
; Enables the spatial measurement tool
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MEASURE_ENABLE', /IGNORE_LAST
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_label, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_text, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_unit, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_lab, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_lab, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SENSITIVE = event.SELECT
	(*(*info).meas).spatial_measurement = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).meas).spatial_measurement], labels=['Enabled measurement']
	IF (event.SELECT EQ 0) THEN BEGIN
		(*(*info).meas).np = 0
		*(*(*info).meas).xp = 0
		*(*(*info).meas).yp = 0
		*(*(*info).meas).sxp = 0
		*(*(*info).meas).syp = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SET_VALUE = '0'
		WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SET_VALUE = '0'
		(*(*info).curs).lockset = event.SELECT
		WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, SET_BUTTON = ABS((*(*info).curs).lockset-1)
		WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = (*(*info).curs).lockset
		CRISPEX_COORDSLIDERS_SET, 1, 1, event
		CRISPEX_DRAW, event
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).loop_slit_but, SENSITIVE = ABS(event.SELECT-1)
	WIDGET_CONTROL, (*(*info).ctrlscp).loop_feedb_but, SENSITIVE = ABS(event.SELECT-1)
END

PRO CRISPEX_MEASURE_CALC, event
; Computes the actual distance from the measurement line length and resolution
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_MEASURE_CALC'
  au = 149697870. 
	delta_x = FLOAT((*(*(*info).meas).xp)[1]-(*(*(*info).meas).xp)[0]) * (*(*info).dataparams).dx
	delta_y = FLOAT((*(*(*info).meas).yp)[1]-(*(*(*info).meas).yp)[0]) * (*(*info).dataparams).dy
  ; Check whether dx/dy units are arcsec or (c/k/M)m
  IF (STRCMP((*(*info).dataparams).xunit,'arcsec') OR $
      STRCMP((*(*info).dataparams).xunit,'asec')) THEN BEGIN
  	delta_asec = SQRT(delta_x^2 + delta_y^2)
;	delta_asec = delta_pix * FLOAT((*(*info).meas).arcsecpix)
  	delta_km = au * TAN(delta_asec / 3600. * !DTOR)
  ENDIF ELSE BEGIN
    IF STRCMP((*(*info).dataparams).xunit,'m') THEN fact = 1E3 ELSE $
      IF STRCMP((*(*info).dataparams).xunit,'cm') THEN fact = 1E5 ELSE $
      IF STRCMP((*(*info).dataparams).xunit,'km') THEN fact = 1 ELSE $
      IF STRCMP((*(*info).dataparams).xunit,'Mm') THEN fact = 1E-3 
    delta_km = SQRT(delta_x^2 + delta_y^2) / fact
    delta_asec = ATAN(delta_km / au) / !DTOR * 3600.
  ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SET_VALUE = STRTRIM(delta_asec,2)
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SET_VALUE = STRTRIM(delta_km,2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [delta_asec,delta_km], labels=['Arcseconds','Kilometers']
END

;================================================================================= PLAYBACK PROCEDURES
PRO CRISPEX_PB_BG, event
; Handles the actual playback, given the mode of playback
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_BG'
	CASE (*(*info).pbparams).mode OF
		'PAUSE'	: BEGIN
				IF (*(*info).pbparams).spmode THEN BEGIN
					(*(*info).dataparams).lp += (*(*info).pbparams).spdirection * (*(*info).pbparams).lp_step
					(*(*info).pbparams).spdirection *= -1
					IF ((*(*info).dataparams).lp GT (*(*info).dispparams).lp_upp) THEN (*(*info).dataparams).lp -= (*(*info).dispparams).lp_range ELSE $
						IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN (*(*info).dataparams).lp += (*(*info).dispparams).lp_range
			  	ENDIF ELSE IF (*(*info).pbparams).imrefmode THEN BEGIN
					(*(*info).winids).imrefdisp = ABS((*(*info).winids).imrefdisp-1)
				ENDIF ELSE RETURN
				WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
			  END
		'PLAY'	: BEGIN
				IF ((*(*info).feedbparams).count_pbstats EQ 0) THEN (*(*info).feedbparams).pbstats = SYSTIME(/SECONDS)
				(*(*info).dataparams).t += (*(*info).pbparams).direction * (*(*info).pbparams).t_step
				IF (*(*info).pbparams).spmode THEN BEGIN
					(*(*info).dataparams).lp += (*(*info).pbparams).spdirection * (*(*info).pbparams).lp_step
					(*(*info).pbparams).spdirection *= -1
					IF ((*(*info).dataparams).lp GT (*(*info).dispparams).lp_upp) THEN (*(*info).dataparams).lp -= (*(*info).dispparams).lp_range ELSE $
						IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN (*(*info).dataparams).lp += (*(*info).dispparams).lp_range
					WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
			  	ENDIF ELSE IF (*(*info).pbparams).imrefmode THEN BEGIN
					(*(*info).winids).imrefdisp = ABS((*(*info).winids).imrefdisp-1)
				ENDIF
				CASE (*(*info).pbparams).lmode OF
					'LOOP'	: BEGIN
							IF (*(*info).dataparams).t GT (*(*info).dispparams).t_upp THEN (*(*info).dataparams).t -= (*(*info).dispparams).t_range
							IF (*(*info).dataparams).t LT (*(*info).dispparams).t_low THEN (*(*info).dataparams).t += (*(*info).dispparams).t_range
						  END
					'CYCLE'	: BEGIN
							IF (*(*info).dataparams).t GT (*(*info).dispparams).t_upp THEN BEGIN
								(*(*info).pbparams).direction = -1
								(*(*info).dataparams).t = (*(*info).dispparams).t_upp - ((*(*info).dataparams).t MOD (*(*info).dispparams).t_upp)
							ENDIF ELSE IF (*(*info).dataparams).t LT (*(*info).dispparams).t_low THEN BEGIN
								(*(*info).pbparams).direction = 1
								IF ((*(*info).dispparams).t_low EQ (*(*info).dispparams).t_first) THEN (*(*info).dataparams).t *= -1 ELSE $
									(*(*info).dataparams).t = -1 * ((*(*info).dataparams).t - (*(*info).dispparams).t_low) + (*(*info).dispparams).t_low
							ENDIF
							IF ((*(*info).pbparams).direction GT 0) THEN CRISPEX_PB_BUTTONS_SET, event, /FWD_SET ELSE CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
						  END
					'BLINK'	: BEGIN
							(*(*info).pbparams).direction *= -1
							IF (*(*info).dataparams).t GT (*(*info).dispparams).t_upp THEN (*(*info).dataparams).t -=(*(*info).dispparams).t_range ELSE $
								IF (*(*info).dataparams).t LT (*(*info).dispparams).t_low THEN (*(*info).dataparams).t += (*(*info).dispparams).t_range
							IF ((*(*info).pbparams).direction GT 0) THEN CRISPEX_PB_BUTTONS_SET, event, /FWD_SET ELSE CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
						  END
				ENDCASE
			  END
	ENDCASE
	WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dataparams).t
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 1. / (*(*info).pbparams).t_speed
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2),STRTRIM((*(*info).pbparams).direction,2),$
		STRTRIM((*(*info).pbparams).spmode,2),STRTRIM((*(*info).dataparams).lp,2),STRTRIM((*(*info).pbparams).spdirection,2)], $
		labels=['Play mode','Loop mode','t','Play direction','Spectral blink mode','lp','Blink direction']
	CRISPEX_UPDATE_T, event
	IF (*(*info).pbparams).spmode THEN CRISPEX_UPDATE_LP, event
	IF (*(*info).dispparams).phislice_update THEN CRISPEX_UPDATE_SLICES, event, /NO_DRAW
	CRISPEX_DRAW, event
	IF (((*(*info).feedbparams).verbosity)[4] EQ 1) THEN BEGIN
		(*(*info).feedbparams).count_pbstats += 1
		newtime = SYSTIME(/SECONDS)
		timediff = newtime-(*(*info).feedbparams).pbstats
		(*(*(*info).feedbparams).sum_pbstats)[(*(*info).feedbparams).count_pbstats MOD 10] = timediff
		(*(*info).feedbparams).av_pbstats = (((*(*info).feedbparams).av_pbstats)*((*(*info).feedbparams).count_pbstats-1) + timediff)/FLOAT((*(*info).feedbparams).count_pbstats)
		(*(*info).feedbparams).pbstats = newtime
		IF ((*(*info).feedbparams).count_pbstats GE 10) THEN average = STRING(MEAN(*(*(*info).feedbparams).sum_pbstats),FORMAT='(F6.4)')+' s' ELSE average = 'N/A'
		CRISPEX_UPDATE_USER_FEEDBACK, event, title='CRISPEX DEBUGGING: Playback statistics', var=((*(*info).feedbparams).count_pbstats+(*(*info).winids).feedbacktlb), minvar=1, /close_button, $
			feedback_text='Time elapsed since last update: '+STRING(timediff,FORMAT='(F6.4)')+' s, average (over last 10): '+average+', (over last '+$
			STRING((*(*info).feedbparams).count_pbstats,FORMAT='(I'+STRTRIM(FLOOR(ALOG10((*(*info).feedbparams).count_pbstats))+1,2)+')')+'): '+STRING((*(*info).feedbparams).av_pbstats,FORMAT='(F6.4)')+' s'
	ENDIF
END

PRO CRISPEX_PB_BUTTONS_SET, event, fbwd_set=fbwd_set, bwd_set=bwd_set, pause_set=pause_set, fwd_set=fwd_set, ffwd_set=ffwd_set, loop_set=loop_set, cycle_set=cycle_set, blink_set=blink_set
; Sets playback buttons according to actions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_BUTTONS_SET'
	IF (KEYWORD_SET(LOOP_SET) OR KEYWORD_SET(CYCLE_SET) OR KEYWORD_SET(BLINK_SET)) THEN BEGIN
		IF KEYWORD_SET(LOOP_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).loop_button, SET_VALUE = (*(*info).ctrlspbbut).loop_pressed, /SET_BUTTON ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_button, SET_VALUE = (*(*info).ctrlspbbut).loop_idle
		IF KEYWORD_SET(CYCLE_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).cycle_button, SET_VALUE = (*(*info).ctrlspbbut).cycle_pressed, /SET_BUTTON ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).cycle_button, SET_VALUE = (*(*info).ctrlspbbut).cycle_idle
		IF KEYWORD_SET(BLINK_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).blink_button, SET_VALUE = (*(*info).ctrlspbbut).blink_pressed, /SET_BUTTON ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).blink_button, SET_VALUE = (*(*info).ctrlspbbut).blink_idle
	ENDIF ELSE BEGIN
		IF KEYWORD_SET(FBWD_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).fbwd_button, SET_VALUE = (*(*info).ctrlspbbut).fbwd_pressed ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).fbwd_button, SET_VALUE = (*(*info).ctrlspbbut).fbwd_idle
		IF KEYWORD_SET(BWD_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).bwd_button, SET_VALUE = (*(*info).ctrlspbbut).bwd_pressed ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).bwd_button, SET_VALUE = (*(*info).ctrlspbbut).bwd_idle
		IF KEYWORD_SET(PAUSE_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).pause_button, SET_VALUE = (*(*info).ctrlspbbut).pause_pressed ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).pause_button, SET_VALUE = (*(*info).ctrlspbbut).pause_idle
		IF KEYWORD_SET(FWD_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).fwd_button, SET_VALUE = (*(*info).ctrlspbbut).fwd_pressed ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).fwd_button, SET_VALUE = (*(*info).ctrlspbbut).fwd_idle
		IF KEYWORD_SET(FFWD_SET) THEN WIDGET_CONTROL, (*(*info).ctrlscp).ffwd_button, SET_VALUE = (*(*info).ctrlspbbut).ffwd_pressed ELSE $
			WIDGET_CONTROL, (*(*info).ctrlscp).ffwd_button, SET_VALUE = (*(*info).ctrlspbbut).ffwd_idle
	ENDELSE
END

PRO CRISPEX_PB_BACKWARD, event
; Sets the playback mode to backward play
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_BACKWARD'
	IF (*(*info).pbparams).mode EQ 'PLAY' AND (*(*info).pbparams).direction EQ -1 THEN RETURN
	(*(*info).pbparams).direction = -1			&	(*(*info).pbparams).mode = 'PLAY'
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,STRTRIM((*(*info).pbparams).direction,2)], labels=['Play mode','Play direction']
END

PRO CRISPEX_PB_FORWARD, event
; Sets the playback mode to forward play
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_FORWARD'
	IF (*(*info).pbparams).mode EQ 'PLAY' AND (*(*info).pbparams).direction EQ 1 THEN RETURN
	(*(*info).pbparams).direction = 1			&	(*(*info).pbparams).mode = 'PLAY'
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	CRISPEX_PB_BUTTONS_SET, event, /FWD_SET
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,STRTRIM((*(*info).pbparams).direction,2)], labels=['Play mode','Play direction']
END

PRO CRISPEX_PB_FASTBACKWARD, event
; Sets the playback mode to single step backward
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_FASTBACKWARD'
	IF (*(*info).pbparams).mode EQ 'PAUSE' OR (*(*info).pbparams).lmode EQ 'BLINK' THEN BEGIN
		CRISPEX_PB_BUTTONS_SET, event, /FBWD_SET
		(*(*info).dataparams).t -= (*(*info).pbparams).t_step
		IF ((*(*info).dataparams).t LT (*(*info).dispparams).t_low) THEN (*(*info).dataparams).t = (*(*info).dispparams).t_upp
		CRISPEX_UPDATE_T, event
		WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dataparams).t
		CRISPEX_UPDATE_SLICES, event
		CRISPEX_DRAW, event
		CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
END

PRO CRISPEX_PB_FASTFORWARD, event
; Sets the playback mode to single step forward
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_FASTFORWARD'
	IF (*(*info).pbparams).mode EQ 'PAUSE' OR (*(*info).pbparams).lmode EQ 'BLINK' THEN BEGIN
		CRISPEX_PB_BUTTONS_SET, event, /FFWD_SET
		(*(*info).dataparams).t = (((*(*info).dataparams).t - (*(*info).dispparams).t_low + (*(*info).pbparams).t_step) MOD ((*(*info).dispparams).t_upp - (*(*info).dispparams).t_low + 1)) + (*(*info).dispparams).t_low
		CRISPEX_UPDATE_T, event
		WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dataparams).t
		CRISPEX_UPDATE_SLICES, event
		CRISPEX_DRAW, event
		CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
END

PRO CRISPEX_PB_PAUSE, event
; Sets the playback mode pause
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_PAUSE'
	IF (*(*info).pbparams).mode EQ 'PAUSE' THEN RETURN
	(*(*info).pbparams).mode = 'PAUSE'
	IF ((*(*info).winids).feedbacktlb NE 0) THEN BEGIN
		(*(*info).feedbparams).count_pbstats = 0
		WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, /SENSITIVE
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	IF ((*(*info).dispparams).phislice_update NE 1) THEN CRISPEX_UPDATE_SLICES, event		
END

PRO CRISPEX_PB_BLINK, event
; Sets the playback mode to temporal blink
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_BLINK'
	(*(*info).pbparams).lmode = 'BLINK'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /BLINK_SET
END

PRO CRISPEX_PB_CYCLE, event
; Sets the playback mode to cycle
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_CYCLE'
	(*(*info).pbparams).lmode = 'CYCLE'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /CYCLE_SET
END

PRO CRISPEX_PB_LOOP, event
; Sets the playback mode to loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_LOOP'
	(*(*info).pbparams).lmode = 'LOOP'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dataparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /LOOP_SET
END

PRO CRISPEX_PB_SPECTBLINK, event
; Sets the playback mode to spectral blink
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PB_SPECTBLINK'
	(*(*info).pbparams).spmode = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).imref_blink_but, SENSITIVE=ABS((*(*info).pbparams).spmode-1)
	IF (*(*info).pbparams).spmode THEN BEGIN
		(*(*info).pbparams).spmode = 1	&	(*(*info).pbparams).spdirection = 1
		IF ((*(*info).feedbparams).count_pbstats EQ 0) THEN (*(*info).feedbparams).pbstats = SYSTIME(/SECONDS)
		WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	ENDIF ELSE BEGIN
		IF ((*(*info).winids).feedbacktlb NE 0) THEN BEGIN
			(*(*info).feedbparams).count_pbstats = 0
			WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, /SENSITIVE
		ENDIF
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_speed_slider, SENSITIVE = (*(*info).pbparams).spmode
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).spmode,(*(*info).pbparams).spdirection], labels=['Spectral blink mode','Blink direction']
END

;================================================================================= SPECTRAL PHI SLIT PROCEDURES
PRO CRISPEX_PHISLIT_GET_SLICE, event
; Gets the spectral phi slit slice
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PHISLIT_GET_SLICE'
	x_pts = COS(!DTOR * (*(*info).phiparams).angle) * (FINDGEN( 2 * (*(*info).phiparams).nphi_set/2) - (*(*info).phiparams).nphi_set/2) + FIX((*(*info).dataparams).x)
	y_pts = SIN(!DTOR * (*(*info).phiparams).angle) * (FINDGEN( 2 * (*(*info).phiparams).nphi_set/2) - (*(*info).phiparams).nphi_set/2) + FIX((*(*info).dataparams).y)
	w = WHERE((x_pts GE 0 ) AND (x_pts LE (*(*info).dispparams).x_last) AND (y_pts GE 0) AND (y_pts LE (*(*info).dispparams).y_last), nw)
	(*(*info).phiparams).nw_cur = nw
	x_pts = REBIN(x_pts[w], nw, (*(*info).dataparams).nlp)
	y_pts = REBIN(y_pts[w], nw, (*(*info).dataparams).nlp)
	lp_pts = REBIN(FINDGEN(1,(*(*info).dataparams).nlp), nw, (*(*info).dataparams).nlp)
	tmp = INTERPOLATE( *(*(*info).data).phiscan, x_pts, y_pts, lp_pts)
	midphi = WHERE(w EQ (*(*info).phiparams).nphi_set/2)
	(*(*info).phiparams).sphi = midphi +  ((*(*info).phiparams).nphi - nw)/2
	delta = (*(*info).phiparams).sphi - ((*(*info).phiparams).nphi)/2
	*(*(*info).data).phislice = TRANSPOSE(tmp, [1,0])
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [nw,midphi], labels=['nw','midphi']
END

PRO CRISPEX_PHISLIT_DIRECTION, event
; Determines the direction of the spectral phi slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PHISLIT_DIRECTION'
	x_maxpts = COS(!DTOR * (*(*info).phiparams).angle) * (FINDGEN( 2 * (*(*info).phiparams).nphi/2) - (*(*info).phiparams).nphi/2) + FIX((*(*info).dataparams).x)
	y_maxpts = SIN(!DTOR * (*(*info).phiparams).angle) * (FINDGEN( 2 * (*(*info).phiparams).nphi/2) - (*(*info).phiparams).nphi/2) + FIX((*(*info).dataparams).y)
	wmax = WHERE((x_maxpts GE 0 ) AND (x_maxpts LE (*(*info).dispparams).x_last) AND (y_maxpts GE 0) AND (y_maxpts LE (*(*info).dispparams).y_last), nwmax)
	thex = REBIN(x_maxpts[wmax], nwmax, 2)
	they = REBIN(y_maxpts[wmax], nwmax, 2)
	thep = REBIN(FINDGEN(1,2), nwmax, 2)
	*(*(*info).data).indices = INTERPOLATE( *(*(*info).data).indexmap, thex, they, thep)
	search = WHERE(((*(*(*info).data).indices)[*,0] EQ FIX((*(*info).dataparams).x)) AND ((*(*(*info).data).indices)[*,1] EQ FIX((*(*info).dataparams).y)))
	(*(*info).phiparams).curindex = search[0]
	(*(*info).phiparams).maxindex = (SIZE(wmax))[1]
	IF ((*(*info).phiparams).curindex GT 1) OR ((*(*info).phiparams).curindex NE (*(*info).phiparams).maxindex-1) THEN BEGIN
		IF (*(*info).ctrlsswitch).bwd_insensitive THEN WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = 1
		IF (*(*info).ctrlsswitch).fwd_insensitive THEN WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = 1
	ENDIF ELSE IF ((*(*info).phiparams).curindex LT 1) OR ((*(*info).phiparams).curindex EQ (*(*info).phiparams).maxindex-1) THEN BEGIN
		IF ((*(*info).ctrlsswitch).bwd_insensitive EQ 0) THEN WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = 0
		IF ((*(*info).ctrlsswitch).fwd_insensitive EQ 0) THEN WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = 0
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x,(*(*info).dataparams).y,(*(*info).phiparams).curindex,(*(*info).phiparams).maxindex], $
		labels=['x','y','current index', 'maximum index']
END

PRO CRISPEX_PHISLIT_MOVE_BWD, event
; Handles the movement of the spectral phi slit backward along the slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PHISLIT_MOVE_BWD'
	newindex = (*(*info).phiparams).curindex-1
	(*(*info).phiparams).curindex = newindex
	(*(*info).dataparams).x = ((*(*(*info).data).indices)[newindex,*])[0]
	(*(*info).dataparams).y = ((*(*(*info).data).indices)[newindex,*])[1]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [newindex+1,newindex,(*(*info).dataparams).x,(*(*info).dataparams).y], labels=['Previous index','New index','x','y']
	CRISPEX_PHISLIT_MOVE_SET_COORDS, event
	IF ((*(*info).phiparams).curindex EQ 0) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = 0
		(*(*info).ctrlsswitch).bwd_insensitive = 1
	ENDIF
	IF (*(*info).ctrlsswitch).fwd_insensitive THEN WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = 1
	CRISPEX_UPDATE_PHIS, event 
	CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_PHISLIT_MOVE_FWD, event
; Handles the movement of the spectral phi slit forward along the slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PHISLIT_MOVE_FWD'
	newindex = (*(*info).phiparams).curindex+1
	(*(*info).phiparams).curindex = newindex
	(*(*info).dataparams).x = ((*(*(*info).data).indices)[newindex,*])[0]
	(*(*info).dataparams).y = ((*(*(*info).data).indices)[newindex,*])[1]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [newindex-1,newindex,(*(*info).dataparams).x,(*(*info).dataparams).y], labels=['Previous index','New index','x','y']
	CRISPEX_PHISLIT_MOVE_SET_COORDS, event
	IF ((*(*info).phiparams).curindex EQ (*(*info).phiparams).maxindex-1) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = 0
		(*(*info).ctrlsswitch).fwd_insensitive = 1
	ENDIF
	IF (*(*info).ctrlsswitch).bwd_insensitive THEN WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = 1
	CRISPEX_UPDATE_PHIS, event 
	CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_PHISLIT_MOVE_SET_COORDS, event
; Handles the coordinate update after the movement of the spectral phi slit along the slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PHISLIT_MOVE_SET_COORDS'
	IF ((*(*info).zooming).factor EQ 1) THEN BEGIN
		(*(*info).curs).sx = (*(*info).dataparams).x * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
		(*(*info).curs).sy = (*(*info).dataparams).y * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
	ENDIF ELSE BEGIN
		(*(*info).curs).sx = (*(*info).dataparams).x * ((*(*info).dataparams).d_nx+1) / (*(*info).dataparams).nx
		(*(*info).curs).sy = (*(*info).dataparams).y * ((*(*info).dataparams).d_ny+1) / (*(*info).dataparams).ny
	ENDELSE
	IF (*(*info).curs).lockset EQ 1 THEN BEGIN
		(*(*info).curs).xlock = (*(*info).dataparams).x	&	(*(*info).curs).ylock = (*(*info).dataparams).y
		(*(*info).curs).sxlock= (*(*info).curs).sx	&	(*(*info).curs).sylock= (*(*info).curs).sy
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).curs).sx,(*(*info).curs).sy,(*(*info).curs).xlock,(*(*info).curs).ylock,(*(*info).curs).sxlock,(*(*info).curs).sylock], $
		labels=['sx','sy','xlock','ylock','sxlock','sylock']
END

;================================================================================= PREFERENCES ROUTINES
PRO CRISPEX_PREFERENCES_WINDOW, event
; Opens up the preferences window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_WINDOW'
	base 		= WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Preferences', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	main		= WIDGET_BASE(base, /COLUMN)
	startup_base	= WIDGET_BASE(main, /COLUMN, /FRAME)
	startup_lab 	= WIDGET_LABEL(startup_base, VALUE = 'At start-up:', /ALIGN_LEFT)
	startup_buts 	= WIDGET_BASE(startup_base, /COLUMN, /NONEXCLUSIVE)
	(*(*info).ctrlspref).startup_win	= WIDGET_BUTTON(startup_buts, VALUE = 'Show start-up window', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_STARTUPWIN')
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_win, SET_BUTTON = (*(*info).prefs).startupwin
	(*(*info).ctrlspref).startup_autopl	= WIDGET_BUTTON(startup_buts, VALUE = 'Start playing automatically', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_AUTOPLAY')
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_autopl, SET_BUTTON = (*(*info).prefs).autoplay
	displays_base	= WIDGET_BASE(main,/COLUMN, /FRAME)
	displays_lab 	= WIDGET_LABEL(displays_base, VALUE = 'Display options:', /ALIGN_LEFT)
	displays_buts 	= WIDGET_BASE(displays_base, /GRID_LAYOUT, COLUMN=2)
	(*(*info).prefs).bgplotcol_old = (*(*info).plotparams).bgplotcol
	(*(*info).ctrlspref).displays_bgcols = WIDGET_SLIDER(displays_buts, TITLE = 'Default background plot color', MIN = 0, MAX = 255, VALUE = (*(*info).plotparams).bgplotcol, /DRAG, $
		EVENT_PRO = 'CRISPEX_PREFERENCES_SET_BGPLOTCOL')
	(*(*info).prefs).plotcol_old = (*(*info).plotparams).plotcol
	(*(*info).ctrlspref).displays_plcols= WIDGET_SLIDER(displays_buts, TITLE = 'Default line plot color', MIN = 0, MAX =255, VALUE = (*(*info).plotparams).plotcol, /DRAG, EVENT_PRO = 'CRISPEX_PREFERENCES_SET_PLOTCOL')
	displays_int_base= WIDGET_BASE(displays_base, /ROW, /NONEXCLUSIVE)
	(*(*info).prefs).interpspslice_old = (*(*info).dispparams).interpspslice
	(*(*info).ctrlspref).displays_interp = WIDGET_BUTTON(displays_int_base, VALUE = 'Interpolate spectral slices', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_INTERPOLATE')
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_interp, SET_BUTTON = (*(*info).dispparams).interpspslice
	displays_preview= WIDGET_BUTTON(displays_int_base, VALUE = 'Preview', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_PREVIEW')
	WIDGET_CONTROL, displays_preview, SET_BUTTON = (*(*info).prefs).preview
	displays_opts	= WIDGET_BASE(displays_base, /COLUMN, /NONEXCLUSIVE)
	(*(*info).ctrlspref).displays_slices = WIDGET_BUTTON(displays_opts, VALUE = 'Scale slices according main/reference image', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_SLICES_IMSCALE')
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_slices, SET_BUTTON = (*(*info).dispparams).slices_imscale
	(*(*info).ctrlspref).displays_phislice = WIDGET_BUTTON(displays_opts, VALUE = 'Automatically update spectral Phi-slice', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_PHISLICE_UPDATE')
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_phislice, SET_BUTTON = (*(*info).dispparams).phislice_update
	paths_base	= WIDGET_BASE(main, /COLUMN, /FRAME)
	paths_lab	= WIDGET_LABEL(paths_base, VALUE = 'Paths:', /ALIGN_LEFT)
	paths_i_base	= WIDGET_BASE(paths_base, /COLUMN)
	paths_i_labbuts= WIDGET_BASE(paths_i_base, /ROW)
	paths_i_lab	= WIDGET_LABEL(paths_i_labbuts, VALUE = 'Default input path:')
	paths_i_buts	= WIDGET_BASE(paths_i_labbuts, /ROW, /EXCLUSIVE)
	(*(*info).ctrlspref).paths_i_def_but = WIDGET_BUTTON(paths_i_buts, VALUE = 'Local working directory', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_IPATH_SEL_DEFAULT', /NO_RELEASE)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_def_but, SET_BUTTON = ABS((*(*info).prefs).defipath-1)
	(*(*info).ctrlspref).paths_i_sav_but = WIDGET_BUTTON(paths_i_buts, VALUE = 'Other directory', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_IPATH_SEL_OTHER', /NO_RELEASE)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_sav_but, SET_BUTTON = (*(*info).prefs).defipath
	(*(*info).ctrlspref).paths_ipath_text = WIDGET_TEXT(paths_i_base, VALUE = (*(*info).prefs).prefipath, /EDITABLE, EVENT_PRO = 'CRISPEX_PREFERENCES_SET_IPATH_OTHER', SENSITIVE = (*(*info).prefs).defipath)
	paths_o_base	= WIDGET_BASE(paths_base, /COLUMN)
	paths_o_labbuts= WIDGET_BASE(paths_o_base, /ROW)
	paths_o_lab	= WIDGET_LABEL(paths_o_labbuts, VALUE = 'Default output path:')
	paths_o_buts	= WIDGET_BASE(paths_o_labbuts, /ROW, /EXCLUSIVE)
	(*(*info).ctrlspref).paths_o_def_but = WIDGET_BUTTON(paths_o_buts, VALUE = 'Local working directory', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_OPATH_SEL_DEFAULT', /NO_RELEASE)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_def_but, SET_BUTTON = ABS((*(*info).prefs).defopath-1)
	(*(*info).ctrlspref).paths_o_sav_but = WIDGET_BUTTON(paths_o_buts, VALUE = 'Other directory', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_OPATH_SEL_OTHER', /NO_RELEASE)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_sav_but, SET_BUTTON = (*(*info).prefs).defopath
	(*(*info).ctrlspref).paths_opath_text = WIDGET_TEXT(paths_o_base, VALUE = (*(*info).prefs).prefopath, /EDITABLE, EVENT_PRO = 'CRISPEX_PREFERENCES_SET_OPATH_OTHER', SENSITIVE = (*(*info).prefs).defopath)
	(*(*info).ctrlspref).paths_iopath = WIDGET_BUTTON(paths_base, VALUE = 'Set output path to input path', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_IOPATH', SENSITIVE = ((*(*info).prefs).prefipath NE (*(*info).prefs).prefopath))
	save_base	= WIDGET_BASE(main, /COLUMN, /FRAME)
	save_lab	= WIDGET_LABEL(save_base, VALUE = 'Saving:', /ALIGN_LEFT)
	saveid_buts 	= WIDGET_BASE(save_base, /ROW)
	saveid_lab	= WIDGET_LABEL(saveid_buts, VALUE = 'Default unique file ID:', /ALIGN_LEFT)
	saveids		= ['YYYYMMMDD_hhmmss (default)','DDMMMYYYY_hhmmss', 'YYYYMMDD_hhmmss','DDMMYYYY_hhmmss']
	(*(*info).ctrlspref).save_defsaveid= WIDGET_COMBOBOX(saveid_buts, VALUE = saveids, /DYNAMIC_RESIZE, EVENT_PRO = 'CRISPEX_PREFERENCES_SET_SAVEID')
	CRISPEX_SAVE_DETERMINE_SAVEID, event, defsaveid_sample, /PREF
	saveid_sample	= WIDGET_BASE(save_base, /ROW)
	saveid_sample_lab = WIDGET_LABEL(saveid_sample, VALUE = 'Example:')
	(*(*info).ctrlspref).save_defsaveid_sample = WIDGET_LABEL(saveid_sample, VALUE = defsaveid_sample, /DYNAMIC_RESIZE)
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid, SET_COMBOBOX_SELECT = (*(*info).prefs).defsaveid
	dec_buts 	= WIDGET_BASE(main, /ALIGN_CENTER, /GRID_LAYOUT, COLUMN=3)
	(*(*info).ctrlspref).set_defaults 	= WIDGET_BUTTON(dec_buts, VALUE = 'Default settings', EVENT_PRO = 'CRISPEX_PREFERENCES_SET_DEFAULTS', SENSITIVE = nondefault)
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
	cancel 		= WIDGET_BUTTON(dec_buts, VALUE = 'Cancel', EVENT_PRO = 'CRISPEX_PREFERENCES_CANCEL')
	save_settings 	= WIDGET_BUTTON(dec_buts, VALUE = 'Save settings', EVENT_PRO = 'CRISPEX_PREFERENCES_SAVE_SETTINGS')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).spxoffset, TLB_SET_YOFFSET = 0
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).preftlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).preftlb], labels=['preftlb']
END

PRO CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
; Handles the checking whether preference buttons and values are in their default position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON'
	IF (((*(*info).prefs).tmp_autoplay NE (*(*info).prefs).default_autoplay) OR $
		((*(*info).prefs).tmp_startupwin NE (*(*info).prefs).default_startupwin) OR $
		((*(*info).prefs).tmp_bgplotcol NE (*(*info).prefs).default_bgplotcol) OR $
		((*(*info).prefs).tmp_plotcol NE (*(*info).prefs).default_plotcol) OR $
		((*(*info).prefs).tmp_interpspslice NE (*(*info).prefs).default_interpspslice) OR $
		((*(*info).prefs).tmp_slices_imscale NE (*(*info).prefs).default_slices_imscale) OR $			
		((*(*info).prefs).tmp_phislice_update NE (*(*info).prefs).default_phislice_update) OR $			
		((*(*info).prefs).tmp_defipath NE (*(*info).prefs).default_defipath) OR $
		((*(*info).prefs).tmp_prefipath NE (*(*info).prefs).default_prefipath) OR $
		((*(*info).prefs).tmp_defopath NE (*(*info).prefs).default_defopath) OR $
		((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).default_prefopath) OR $
		((*(*info).prefs).tmp_defsaveid NE (*(*info).prefs).default_defsaveid)) THEN nondefault = 1 ELSE nondefault = 0
	WIDGET_CONTROL, (*(*info).ctrlspref).set_defaults, SENSITIVE = nondefault
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [ABS(nondefault-1)], labels=['Buttons on default']
END

PRO CRISPEX_PREFERENCES_SET_STARTUPWIN, event
; Handles the toggle on/off setting of start-up window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_STARTUPWIN'
	(*(*info).prefs).tmp_startupwin = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_startupwin], labels=['Startup window']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_AUTOPLAY, event
; Handles the toggle on/off setting of autoplay at start-up
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_AUTOPLAY'
	(*(*info).prefs).tmp_autoplay = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_autoplay], labels=['Autoplay']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_BGPLOTCOL, event
; Handles the setting of the background plot color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_BGPLOTCOL'
	(*(*info).prefs).tmp_bgplotcol = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_bgplotcol], labels=['Background color']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).plotparams).bgplotcol = (*(*info).prefs).tmp_bgplotcol
		CRISPEX_PREFERENCES_REDRAW, event
	ENDIF
END

PRO CRISPEX_PREFERENCES_SET_PHISLICE_UPDATE, event									
; Handles the toggle on/off live update of spectral phi-slice
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_PHISLICE_UPDATE'
	(*(*info).prefs).tmp_phislice_update = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_phislice_update], labels=['Live update spectral Phi-slice']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_PLOTCOL, event
; Handles the setting of the line plot color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_PLOTCOL'
	(*(*info).prefs).tmp_plotcol = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_plotcol], labels=['Plot color']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).plotparams).plotcol = (*(*info).prefs).tmp_plotcol
		CRISPEX_PREFERENCES_REDRAW, event
	ENDIF
END

PRO CRISPEX_PREFERENCES_SET_SLICES_IMSCALE, event									
; Handles the toggle on/off scaling of slices according to main/reference image scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_SLICES_IMSCALE'
	(*(*info).prefs).tmp_slices_imscale = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_slices_imscale], labels=['Scale slices with main/reference image']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_INTERPOLATE, event
; Handles the toggle on/off setting of interpolating the spectral slices
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_INTERPOLATE'
	(*(*info).prefs).tmp_interpspslice = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_interpspslice], labels=['Interpolate']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).dispparams).interpspslice = (*(*info).prefs).tmp_interpspslice
		CRISPEX_PREFERENCES_REDRAW, event
	ENDIF
END

PRO CRISPEX_PREFERENCES_SET_IPATH_SEL_DEFAULT, event
; Handles the toggle on/off selection of the default input path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_IPATH_SEL_DEFAULT'
	(*(*info).prefs).tmp_defipath = ABS(event.SELECT-1)
	(*(*info).prefs).tmp_prefipath = (*(*info).prefs).default_prefipath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, SET_VALUE = (*(*info).prefs).default_prefipath, SENSITIVE = (*(*info).prefs).tmp_defipath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM(event.SELECT,2),(*(*info).prefs).tmp_prefipath], labels=['Default input path set','Input path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_IPATH_SEL_OTHER, event
; Handles the toggle on/off selection of different input path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_IPATH_SEL_OTHER'
	(*(*info).prefs).tmp_defipath = event.SELECT
	IF (*(*info).prefs).tmp_defipath THEN BEGIN
		IF ((*(*info).prefs).tmp_prefipath NE (*(*info).paths).ipath) THEN (*(*info).prefs).tmp_prefipath = (*(*info).paths).ipath
		WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, SET_VALUE = (*(*info).prefs).tmp_prefipath, SENSITIVE = (*(*info).prefs).tmp_defipath
		CRISPEX_SAVE_SET_IPATH, event
		WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, SET_VALUE = (*(*info).prefs).tmp_prefipath
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM(ABS(event.SELECT-1),2),(*(*info).prefs).tmp_prefipath], labels=['Default input path set','Input path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_IPATH_OTHER, event
; Handles the setting of the input path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_IPATH_OTHER'
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, GET_VALUE = textvalue
	(*(*info).prefs).tmp_prefipath = textvalue[0]
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_prefipath], labels=['Input path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_OPATH_SEL_DEFAULT, event
; Handles the toggle on/off selection of the default output path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_OPATH_SEL_DEFAULT'
	(*(*info).prefs).tmp_defopath = ABS(event.SELECT-1)
	(*(*info).prefs).tmp_prefopath = (*(*info).prefs).default_prefopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, SET_VALUE = (*(*info).prefs).default_prefopath, SENSITIVE = (*(*info).prefs).tmp_defopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM(event.SELECT,2),(*(*info).prefs).tmp_prefopath], labels=['Default output path set','Output path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_OPATH_SEL_OTHER, event
; Handles the toggle on/off selection of different output path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_OPATH_SEL_OTHER'
	(*(*info).prefs).tmp_defopath = event.SELECT
	IF (*(*info).prefs).tmp_defopath THEN BEGIN
		IF ((*(*info).prefs).tmp_prefopath NE (*(*info).paths).opath) THEN (*(*info).prefs).tmp_prefopath = (*(*info).paths).opath
		WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, SET_VALUE = (*(*info).prefs).tmp_prefopath, SENSITIVE = (*(*info).prefs).tmp_defopath
		CRISPEX_SAVE_SET_OPATH, event
		WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, SET_VALUE = (*(*info).prefs).tmp_prefopath
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM(ABS(event.SELECT-1),2),(*(*info).prefs).tmp_prefopath], labels=['Default output path set','Output path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_OPATH_OTHER, event
; Handles the setting of the output path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_OPATH_OTHER'
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, GET_VALUE = textvalue
	(*(*info).prefs).tmp_prefopath = textvalue[0]
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_prefopath], labels=['Output path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_IOPATH, event
; Handles the setting of the output path to the input path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_IOPATH'
	(*(*info).prefs).tmp_prefopath = (*(*info).prefs).tmp_prefipath
	(*(*info).prefs).tmp_defopath = (*(*info).prefs).tmp_defipath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_def_but, SET_BUTTON = ABS((*(*info).prefs).tmp_defopath-1)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_sav_but, SET_BUTTON = (*(*info).prefs).tmp_defopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, SET_VALUE = (*(*info).prefs).tmp_prefopath, SENSITIVE = (*(*info).prefs).tmp_defopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM(ABS((*(*info).prefs).tmp_defipath-1),2),(*(*info).prefs).tmp_prefipath,STRTRIM(ABS((*(*info).prefs).tmp_defopath-1),2),$
		(*(*info).prefs).tmp_prefopath], labels=['Default input path set','Input path','Default output path set','Output path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_PREVIEW, event
; Handles the enabling of preview mode
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_PREVIEW'
	(*(*info).prefs).preview = event.SELECT
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).plotparams).plotcol = (*(*info).prefs).tmp_plotcol		&	(*(*info).plotparams).bgplotcol = (*(*info).prefs).tmp_bgplotcol
		(*(*info).dispparams).interpspslice = (*(*info).prefs).tmp_interpspslice
	ENDIF ELSE BEGIN
		(*(*info).plotparams).plotcol = (*(*info).prefs).plotcol_old		& 	(*(*info).plotparams).bgplotcol = (*(*info).prefs).bgplotcol_old
		(*(*info).dispparams).interpspslice = (*(*info).prefs).interpspslice_old
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).preview], labels=['Preview']
	CRISPEX_PREFERENCES_REDRAW, event
END

PRO CRISPEX_PREFERENCES_SET_SAVEID, event
; Handles the setting of unique save ID
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_SAVEID'
	(*(*info).prefs).tmp_defsaveid = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_defsaveid], labels=['Default save ID']
	CRISPEX_SAVE_DETERMINE_SAVEID, event, defsaveid_sample, /PREF
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid_sample, SET_VALUE = defsaveid_sample
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_DEFAULTS, event
; Handles the setting of all defaults
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SET_DEFAULTS'
	(*(*info).prefs).tmp_startupwin = (*(*info).prefs).default_startupwin		&	(*(*info).prefs).tmp_interpspslice = (*(*info).prefs).default_interpspslice
	(*(*info).prefs).tmp_autoplay = (*(*info).prefs).default_autoplay		&	(*(*info).prefs).tmp_defsaveid = (*(*info).prefs).default_defsaveid
	(*(*info).prefs).tmp_defipath = (*(*info).prefs).default_defipath		&	(*(*info).prefs).tmp_prefipath = (*(*info).prefs).default_prefipath
	(*(*info).prefs).tmp_defopath = (*(*info).prefs).default_defopath		&	(*(*info).prefs).tmp_prefopath = (*(*info).prefs).default_prefopath
	(*(*info).prefs).tmp_bgplotcol = (*(*info).prefs).default_bgplotcol		&	(*(*info).prefs).tmp_plotcol = (*(*info).prefs).default_plotcol
	(*(*info).prefs).tmp_phislice_update = (*(*info).prefs).default_phislice_update	&	(*(*info).prefs).tmp_slices_imscale = (*(*info).prefs).default_slices_imscale		
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_win, SET_BUTTON = (*(*info).prefs).tmp_startupwin
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_autopl, SET_BUTTON = (*(*info).prefs).tmp_autoplay
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_bgcols, SET_VALUE = (*(*info).prefs).tmp_bgplotcol
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_plcols, SET_VALUE = (*(*info).prefs).tmp_plotcol
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_interp, SET_BUTTON = (*(*info).prefs).tmp_interpspslice
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_phislice, SET_BUTTON = (*(*info).prefs).tmp_phislice_update		
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_slices, SET_BUTTON = (*(*info).prefs).tmp_slices_imscale
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_def_but, SET_BUTTON = ABS((*(*info).prefs).tmp_defipath-1)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_sav_but, SET_BUTTON = (*(*info).prefs).tmp_defipath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, SET_VALUE = (*(*info).prefs).tmp_prefipath, SENSITIVE = (*(*info).prefs).tmp_defipath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_def_but, SET_BUTTON = ABS((*(*info).prefs).tmp_defopath-1)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_sav_but, SET_BUTTON = (*(*info).prefs).tmp_defopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, SET_VALUE = (*(*info).prefs).tmp_prefopath, SENSITIVE = (*(*info).prefs).tmp_defopath
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid, SET_COMBOBOX_SELECT = (*(*info).prefs).tmp_defsaveid
	CRISPEX_SAVE_DETERMINE_SAVEID, event, defsaveid_sample
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid_sample, SET_VALUE = defsaveid_sample
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).plotparams).plotcol = (*(*info).prefs).tmp_plotcol
		(*(*info).plotparams).bgplotcol = (*(*info).prefs).tmp_bgplotcol
		(*(*info).dispparams).interpspslice = (*(*info).prefs).tmp_interpspslice
		CRISPEX_PREFERENCES_REDRAW, event
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [STRTRIM((*(*info).prefs).tmp_startupwin,2),STRTRIM((*(*info).prefs).tmp_autoplay,2),STRTRIM((*(*info).prefs).tmp_bgplotcol,2),$
		STRTRIM((*(*info).prefs).tmp_plotcol,2),STRTRIM((*(*info).prefs).tmp_interpspslice,2),STRTRIM((*(*info).prefs).preview,2),STRTRIM((*(*info).prefs).tmp_phislice_update,2),$
		STRTRIM((*(*info).prefs).tmp_slices_imscale,2),STRTRIM(ABS((*(*info).prefs).tmp_defipath-1),2),(*(*info).prefs).tmp_prefipath,STRTRIM(ABS((*(*info).prefs).tmp_defopath-1),2),$
		STRTRIM((*(*info).prefs).tmp_prefopath,2),STRTRIM((*(*info).prefs).tmp_defsaveid,2)], labels=['Startup window','Autoplay','Background color','Plot color',$
		'Interpolate','Preview','Update phi-slice','Scale slices','Default input path set','Input path','Default output path set','Output path','Default save ID']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SAVE_SETTINGS, event, RESAVE=resave
; Handles the saving of settings
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_SAVE_SETTINGS'
	CD, CURRENT=curpath
	IF (*(*info).prefs).tmp_defipath THEN BEGIN
		IF ((*(*info).prefs).tmp_prefipath EQ curpath) THEN (*(*info).prefs).tmp_defipath = 0
	ENDIF
	IF (*(*info).prefs).tmp_defopath THEN BEGIN
		IF ((*(*info).prefs).tmp_prefopath EQ curpath) THEN (*(*info).prefs).tmp_defopath = 0
	ENDIF
	curpath = curpath+PATH_SEP()
	startupwin = (*(*info).prefs).tmp_startupwin		&	interpspslice = (*(*info).prefs).tmp_interpspslice
	autoplay = (*(*info).prefs).tmp_autoplay		&	defsaveid = (*(*info).prefs).tmp_defsaveid
	defipath = (*(*info).prefs).tmp_defipath		&	prefipath = (*(*info).prefs).tmp_prefipath
	defopath = (*(*info).prefs).tmp_defopath		&	prefopath = (*(*info).prefs).tmp_prefopath
	bgplotcol = (*(*info).prefs).tmp_bgplotcol		&	plotcol = (*(*info).prefs).tmp_plotcol
	phislice_update = (*(*info).prefs).tmp_phislice_update	&	slices_imscale = (*(*info).prefs).tmp_slices_imscale	
	crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
	SAVE, crispex_version, startupwin, interpspslice, phislice_update, slices_imscale, autoplay, defsaveid, defipath, defopath, bgplotcol, plotcol, prefipath, prefopath, FILENAME = (*(*info).paths).dir_settings+'crispex.cpref'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).dir_settings+'crispex.cpref'],labels=['Written']
	(*(*info).prefs).startupwin = startupwin		&	(*(*info).dispparams).interpspslice = interpspslice
	(*(*info).prefs).autoplay = autoplay			&	(*(*info).prefs).defsaveid = defsaveid
	(*(*info).prefs).defipath = defipath			&	(*(*info).prefs).defopath = defopath
	(*(*info).plotparams).bgplotcol = bgplotcol		&	(*(*info).plotparams).plotcol = plotcol
	(*(*info).prefs).prefipath = prefipath			&	(*(*info).paths).ipath = prefipath
	(*(*info).prefs).prefopath = prefopath			&	(*(*info).paths).opath = prefopath
	(*(*info).dispparams).phislice_update = phislice_update	&	(*(*info).dispparams).slices_imscale = slices_imscale	
	IF ~KEYWORD_SET(RESAVE) THEN BEGIN
		CRISPEX_PREFERENCES_REDRAW, event
		WIDGET_CONTROL, (*(*info).winids).preftlb, /DESTROY
		(*(*info).winids).preftlb = 0
	ENDIF
END

PRO CRISPEX_PREFERENCES_CANCEL, event
; Handles the exiting from the preferences window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_CANCEL'
	IF (*(*info).prefs).preview THEN BEGIN
		(*(*info).plotparams).plotcol = (*(*info).prefs).plotcol_old		& 	(*(*info).plotparams).bgplotcol = (*(*info).prefs).bgplotcol_old
		(*(*info).dispparams).interpspslice = (*(*info).prefs).interpspslice_old
		CRISPEX_PREFERENCES_REDRAW, event
	ENDIF
	WIDGET_CONTROL, (*(*info).winids).preftlb, /DESTROY
	(*(*info).winids).preftlb = 0
END

PRO CRISPEX_PREFERENCES_REDRAW, event
; Handles the redrawing of plot windows in case the preview mode is set
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_PREFERENCES_REDRAW'
	IF ((*(*info).winids).sptlb NE 0) THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
	IF (TOTAL(*(*(*info).winids).restlooptlb) NE 0) THEN CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).retrdettlb NE 0) THEN CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).looptlb NE 0) THEN CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).refsptlb NE 0) THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
	CRISPEX_DRAW_SPECTRAL, event
	CRISPEX_DRAW_TIMESLICES, event
	CRISPEX_DRAW_OTHER, event
END

;================================================================================= READ BMP BUTTONS FUNCTION
FUNCTION CRISPEX_READ_BMP_BUTTONS, filename, srcdir
; Handles the reading of (button) BMP files
	button_dummy = READ_BMP(srcdir+filename)  
	button_dummy = TRANSPOSE(button_dummy, [1,2,0])
	RETURN, button_dummy
END

;================================================================================= READ HEADER PROCEDURE
PRO CRISPEX_READ_FITSHEADER, filename, datatype=datatype, nx=nx, ny=ny, nlp=nlp, nt=nt, ns=ns, $
                             imnt=imnt, spnt=spnt, lps=lps, offset=offset, lptitle=lptitle, $
                             inttitle=inttitle, ttitle=ttitle, dt=dt, bunit=bunit, lpunit=lpunit, $
                             dx=dx, dy=dy, xunit=xunit, yunit=yunit, exten_no=exten_no
; Handles read in of the header of the fits input file
		offset = fitspointer(filename,exten_no=KEYWORD_SET(exten_no),hdr)
		parseheader,hdr,key
		datatype = key.datatype
    ; Spatial: dimensions 1 and 2
    nx = key.nx  &  dx = key.dx  &  xunit = key.xunit
    ny = key.ny  &  dy = key.dy  &  yunit = key.yunit
		nlp = key.nlp
		nt = key.nt
    dt = key.dt
    ns = key.ns
    imnt = key.nlp * key.nt
		spnt = key.nx*key.ny
    lps = key.lam
    bunit = key.bunit
    lpunit = key.lpunit
    lptitle=key.lplab+' ['+lpunit+']'
    inttitle=key.btype+' ['+bunit+']'
    ttitle = key.tlab+' ['+key.tunit+']'
END

PRO CRISPEX_READ_HEADER, filename, header=header, datatype=datatype, dims=dims, nx=nx, ny=ny, $
                         nt=nt, endian=endian, stokes=stokes, ns=ns, diagnostics=diagnostics
; Handles the read in of the header of the input files
	OPENR, lun, filename, /GET_LUN
	rec	= ASSOC(lun, BYTARR(512))	&	header	= STRING(rec[0])
	FREE_LUN, lun
	len	= STRLEN(header)

	search	= 'datatype='	&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown datatype!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	datatype= LONG( STRMID(header, pos + STRLEN(search), 1) )
	search	= 'dims='	&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown number of dimensions!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	dims	= LONG( STRMID(header, pos + STRLEN(search), 1) )
	IF (dims LT 2) OR (dims GT 3) THEN BEGIN
		MESSAGE, /INFO, 'Number of dimensions not supported!'
		PRINT, 'Dimensions: '+dims
		RETALL
	ENDIF
	search	= 'nx='		&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown number of pixels in x-direction!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	pos1	= STRPOS(header, ',', pos)
	nx	= LONG( STRMID(header, pos + STRLEN(search), pos1-pos) )
	search	= 'ny='		&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown number of pixels in y-direction!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	pos1	= STRPOS(header, ',', pos)
	ny	= LONG( STRMID(header, pos + STRLEN(search), pos1-pos) )
	search	= 'nt='		&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown number of pixels in t-direction!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	pos1	= STRPOS(header, ',', pos)
	nt	= LONG( STRMID(header, pos + STRLEN(search), pos1-pos) )
	search	= 'endian='	&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN BEGIN
		MESSAGE, /INFO, 'Unknown endianness!'
		PRINT, 'Header: '+header
		RETALL
	ENDIF
	endian	= STRMID(header, pos + STRLEN(search), 1)
	search	= 'ns='	&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN ns = 1 ELSE ns = STRMID(header, pos + STRLEN(search), 1)
	search	= 'stokes='	&	pos	= STRPOS(header, search)
	IF pos EQ -1 THEN stokes = ['I'] ELSE stokes = STRMID(header, pos + STRLEN(search), 2*ns+1)
	search	= 'diagnostics='	&	pos	= STRPOS(header, search)
	pos1	= STRPOS(header,']',pos)
	IF pos NE -1 THEN diagnostics = STRMID(header, pos + STRLEN(search), pos1-pos)
END
 
;================================================================================= RESTORE LOOPS PROCEDURES
PRO CRISPEX_RESTORE_LOOPS_MAIN, event
; Start the restore loops procedures, opens the menu if CLSAV files are present or otherwise returns an error message
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_MAIN'
	(*(*info).loopswitch).restore_loops = event.SELECT
	CRISPEX_FIND_CSAV, event
	IF ((*(*info).restoreparams).cfilecount GT 0) THEN BEGIN
		IF (*(*info).loopswitch).restore_loops THEN BEGIN
			CRISPEX_RESTORE_LOOPS_MENU, event
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_all, SENSITIVE = 1
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_sav, SENSITIVE = 1
		ENDIF ELSE CRISPEX_RESTORE_LOOPS_MENU_CLOSE, event
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','No saved time slice (*csav) files found corresponding','to the current data file. Unable to produce loop overlays.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
		WIDGET_CONTROL, (*(*info).ctrlscp).overlay_but, SET_BUTTON = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopswitch).restore_loops],labels=['Restoring loops']
END

PRO CRISPEX_RESTORE_LOOPS_MENU, event, set_but_array
; Sets up the restored loops menu and reads in the loop points
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_MENU'
	WIDGET_CONTROL,/HOURGLASS
	filenames = *(*(*info).restoreparams).cfiles
	filecount = (*(*info).restoreparams).cfilecount
	IF (N_ELEMENTS(set_but_array) NE filecount) THEN *(*(*info).restoreparams).sel_loops = INTARR((*(*info).restoreparams).cfilecount)
	eventval = INDGEN(filecount)
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Choose loop overlays', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	disp2 = WIDGET_BASE(disp, /COLUMN, /FRAME)
	sel_allnone = WIDGET_BASE(disp2, /ROW)
	sel_allnone_lab = WIDGET_LABEL(sel_allnone, VALUE = 'Select:', /ALIGN_LEFT)
	sel_allnone_buts = WIDGET_BASE(sel_allnone, /ROW, /EXCLUSIVE)
	(*(*info).ctrlsrestore).sel_all = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'All', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_SEL_ALL')
	(*(*info).ctrlsrestore).sel_none = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'None', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_SEL_NONE')
	disp_txt = REPLICATE('Display time slice ', filecount)
	list_values = ['Do not display any slice', disp_txt+STRTRIM(eventval,2)]
	(*(*info).ctrlsrestore).disp_list = WIDGET_COMBOBOX(sel_allnone, VALUE = list_values, /DYNAMIC_RESIZE, EVENT_PRO = 'CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_SELECT')
	WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_none, /SET_BUTTON
	IF (filecount GT 15) THEN sel_buts = WIDGET_BASE(disp2, /COLUMN, /NONEXCLUSIVE,Y_SCROLL_SIZE = 425) ELSE sel_buts = WIDGET_BASE(disp2, /COLUMN, /NONEXCLUSIVE)
	FOR i=0,filecount-1 DO BEGIN
		IF (N_ELEMENTS(set_but_array) NE filecount) THEN (*(*(*info).restoreparams).sel_loops)[i] = 0
		singlefilename = (STRSPLIT((*(*(*info).restoreparams).cfiles)[i], PATH_SEP(), /EXTRACT))[N_ELEMENTS(STRSPLIT((*(*(*info).restoreparams).cfiles)[i], PATH_SEP(), /EXTRACT))-1]
		name = 'sel_restore_but_'+STRTRIM(i,2)
		but_val = STRTRIM(i,2)+': '+singlefilename
		sel_but = WIDGET_BUTTON(sel_buts, VALUE = but_val, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_MENU_EVENT', UNAME = name) 
		IF (N_ELEMENTS(set_but_array) GT 0) THEN WIDGET_CONTROL, sel_but, SET_BUTTON = set_but_array[i]
		RESTORE,(*(*(*info).restoreparams).cfiles)[i]
		nel = N_ELEMENTS(SIZE(x_coords))						; Support for different saving methods of xp and yp points
		IF (nel EQ 5) THEN BEGIN
			xpcoords = REFORM(x_coords[0,*])
			ypcoords = REFORM(y_coords[0,*])
		ENDIF ELSE BEGIN
			xpcoords = x_coords
			ypcoords = y_coords
		ENDELSE
		IF (i EQ 0) THEN BEGIN
			*(*(*info).restoreparams).xp_restored = xpcoords
			*(*(*info).restoreparams).yp_restored = ypcoords
			*(*(*info).restoreparams).xr_restored = x_loop_pts
			*(*(*info).restoreparams).yr_restored = y_loop_pts
			*(*(*info).restoreparams).psizes = [0, LONG((SIZE(xpcoords))[1])]
			*(*(*info).restoreparams).rsizes = [0, LONG(loop_size)]
			*(*(*info).restoreparams).lp_restored = spect_pos
		ENDIF ELSE BEGIN
			*(*(*info).restoreparams).xp_restored = [*(*(*info).restoreparams).xp_restored, xpcoords]
			*(*(*info).restoreparams).yp_restored = [*(*(*info).restoreparams).yp_restored, ypcoords]
			*(*(*info).restoreparams).xr_restored = [*(*(*info).restoreparams).xr_restored, x_loop_pts]
			*(*(*info).restoreparams).yr_restored = [*(*(*info).restoreparams).yr_restored, y_loop_pts]
			singlepsize = (SIZE(xpcoords))[1]
			*(*(*info).restoreparams).psizes = [*(*(*info).restoreparams).psizes, (*(*(*info).restoreparams).psizes)[i]+LONG(singlepsize)]
			*(*(*info).restoreparams).rsizes = [*(*(*info).restoreparams).rsizes, (*(*(*info).restoreparams).rsizes)[i]+LONG(loop_size)]
			*(*(*info).restoreparams).lp_restored = [*(*(*info).restoreparams).lp_restored, spect_pos]
		ENDELSE
		initial_feedback = 'Processed save file '+STRTRIM(i,2)+': '+singlefilename
		CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring loops...', var=i, maxvar=filecount-1, feedback_text=initial_feedback
		IF (i EQ filecount-1) THEN CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event
	ENDFOR
	decision_base = WIDGET_BASE(disp2,COLUMN=3,/GRID_LAYOUT,/ALIGN_CENTER)
	update_filelist = WIDGET_BUTTON(decision_base, VALUE = 'Update filelist', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_UPDATE_FILELIST')
	(*(*info).ctrlsrestore).open_tanat = WIDGET_BUTTON(decision_base, VALUE = 'Open in TANAT...', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_OPEN_TANAT');, SENSITIVE = 0)
	closebut = WIDGET_BUTTON(decision_base, VALUE = 'Close', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_MENU_CLOSE')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).lsxoffset, TLB_SET_YOFFSET = (*(*info).winsizes).lswiny+1.5*(*(*info).winsizes).ydelta
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).restoretlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).restoretlb,filecount],labels=['restoretlb','Restored loops']
END

PRO CRISPEX_RESTORE_LOOPS_MENU_EVENT, event
; Handles the selection of loops to be restored
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_MENU_EVENT'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).restoreparams).sel_loops)[eventval] = ( (*(*(*info).restoreparams).sel_loops)[eventval] EQ 0) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).restoreparams).sel_loops)[eventval]], labels=['Loop ID','Loop selected']
	CRISPEX_RESTORE_LOOPS_BUTTON_CONDITION, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_RESTORE_LOOPS_BUTTON_CONDITION, event
; Handles the update of buttons after selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_BUTTON_CONDITION'
	condition = WHERE(*(*(*info).restoreparams).sel_loops EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_none, SET_BUTTON = ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1)
	WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_all, SET_BUTTON = (((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).restoreparams).cfilecount))
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).restoreparams).cfilecount)),$
		ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1),N_ELEMENTS(condition)-(TOTAL(condition) EQ -1)], labels=['All selected','None selected','Total selected']
END

PRO CRISPEX_RESTORE_LOOPS_MENU_CLOSE, event
; Handles the closing of the restored loops menu and clean-up of display afterwards
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_MENU_CLOSE'
	IF ((*(*info).dispparams).t_range NE (*(*info).dataparams).nt) THEN CRISPEX_DISPRANGE_T_RESET, event, /NO_DRAW
	IF ((*(*info).dispparams).lp_range NE (*(*info).dataparams).nlp) THEN CRISPEX_DISPRANGE_LP_RESET, event, /NO_DRAW
	IF ((*(*info).dispparams).lp_ref_range NE (*(*info).dataparams).refnlp) THEN CRISPEX_DISPRANGE_LP_REF_RESET, event, /NO_DRAW
	WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_all, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_sav, SENSITIVE = 0
	(*(*info).loopswitch).restore_loops = 0
	(*(*info).restoreparams).disp_loopfile = '0'
	*(*(*info).restoreparams).disp_loopnr = -1
	IF (*(*info).winswitch).showrestloop THEN BEGIN
		FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO BEGIN
			WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[i], /DESTROY
			*(*(*(*info).loopsdata).rest_loopslice[i]) = 0
			*(*(*(*info).loopsdata).rest_loopslab[i]) = 0
			*(*(*(*info).loopsdata).rest_crossloc[i]) = 0
		ENDFOR
		(*(*info).winswitch).showrestloop = 0
		*(*(*info).winids).restlooptlb = 0
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).overlay_but, SET_BUTTON = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, /SENSITIVE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopswitch).restore_loops,(*(*info).winids).restoretlb], labels=['Restoring loops','restoretlb was']
	CRISPEX_DRAW, event
	WIDGET_CONTROL, (*(*info).winids).restoretlb, /DESTROY
	(*(*info).winids).restoretlb = 0
END

PRO CRISPEX_RESTORE_LOOPS_SEL_ALL, event
; Handles selection of all restored loops
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_SEL_ALL'
	*(*(*info).restoreparams).sel_loops = REPLICATE(1,(*(*info).restoreparams).cfilecount)
	FOR i=0,(*(*info).restoreparams).cfilecount-1 DO BEGIN
		name = 'sel_restore_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 1
	ENDFOR		
	CRISPEX_DRAW, event
END

PRO CRISPEX_RESTORE_LOOPS_SEL_NONE, event
; Handles selection of none restored loops
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_SEL_NONE'
	*(*(*info).restoreparams).sel_loops = REPLICATE(0,(*(*info).restoreparams).cfilecount)
	FOR i=0,(*(*info).restoreparams).cfilecount-1 DO BEGIN
		name = 'sel_restore_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 0
	ENDFOR		
	CRISPEX_DRAW, event
END

PRO CRISPEX_RESTORE_LOOPS_UPDATE_FILELIST, event
; Handles the update of the restored file list
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_UPDATE_FILELIST'
	IF ((*(*info).dispparams).t_range NE (*(*info).dataparams).nt) THEN CRISPEX_DISPRANGE_T_RESET, event, /NO_DRAW
	IF ((*(*info).dispparams).lp_range NE (*(*info).dataparams).nlp) THEN CRISPEX_DISPRANGE_LP_RESET, event, /NO_DRAW
	IF ((*(*info).dispparams).lp_ref_range NE (*(*info).dataparams).refnlp) THEN CRISPEX_DISPRANGE_LP_REF_RESET, event, /NO_DRAW
	(*(*info).restoreparams).disp_loopfile = '0'
	*(*(*info).restoreparams).disp_loopnr = -1
	IF (*(*info).winswitch).showrestloop THEN BEGIN
		FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO BEGIN
			WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[i], /DESTROY
			*(*(*(*info).loopsdata).rest_loopslice[i]) = 0
			*(*(*(*info).loopsdata).rest_loopslab[i]) = 0
			*(*(*(*info).loopsdata).rest_crossloc[i]) = 0
		ENDFOR
		(*(*info).winswitch).showrestloop = 0
		*(*(*info).winids).restlooptlb = 0
	ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, /SENSITIVE
	WIDGET_CONTROL, (*(*info).winids).root, SET_UVALUE = info
	WIDGET_CONTROL, event.TOP, /DESTROY
	event.TOP = (*(*info).winids).root
	CRISPEX_RESTORE_LOOPS_MAIN, event
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RESTORE_LOOPS_ALWAYS, event
; Sets the display of restored loops to always (i.e. at all spectral positions) or selected (i.e. at saved positions)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_ALWAYS'
	(*(*info).overlayswitch).overlalways = event.SELECT
	IF ((*(*info).loopswitch).restore_loops EQ 0) THEN CRISPEX_RESTORE_LOOPS_MAIN, event ELSE CRISPEX_DRAW, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayswitch).overlalways], labels=['Overlay loops always']
END

PRO CRISPEX_RESTORE_LOOPS_OPEN_TANAT, event			
; Handles all prior to loading TANAT to analyse the selected loop slice/slab
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_OPEN_TANAT'
	path_tanat = (ROUTINE_INFO('TANAT',/SOURCE)).PATH								; CHeck whether TANAT has already been compiled
	tanat_file = FILE_SEARCH((*(*info).paths).dir_aux,'tanat.pro', /FULLY_QUALIFY_PATH, COUNT=tanat_count)	; Try to find TANAT where it should be
	IF (tanat_count EQ 1) THEN (*(*info).paths).dir_tanat = (*(*info).paths).dir_aux
	IF (path_tanat NE '') THEN BEGIN										; TANAT has been compiled, so check whether it is the one at the expected location
		IF (path_tanat NE (*(*info).paths).dir_tanat+'tanat.pro') THEN BEGIN					; Compiled TANAT is not the correct TANAT
			IF (tanat_count EQ 1) THEN BEGIN								; Correct TANAT exists, so compile that one and run
				dir_tanat = STRMID((*(*info).paths).dir_tanat,0,STRPOS((*(*info).paths).dir_tanat,'/',/REVERSE_SEARCH))
				IF (STRPOS(!PATH,dir_tanat) NE 0) THEN !PATH = dir_tanat+':'+!PATH
				RESOLVE_ROUTINE, 'TANAT'
				CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
			ENDIF ELSE IF ((*(*info).paths).tanat_repointed NE 1) THEN BEGIN				; Correct TANAT can not be found, run anyway but warn
				CRISPEX_WINDOW_OK, event,'WARNING!','CRISPEX could not find TANAT at the expected location','('+(*(*info).paths).dir_aux+'),','but a local copy has been compiled before from',$
					'('+STRMID(path_tanat,0,STRPOS(path_tanat,'/',/REVERSE_SEARCH))+').','Press OK to continue with that local copy of TANAT, or select a different one.',$
					OK_EVENT='CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN', CANCEL_EVENT='CRISPEX_RESTORE_LOOPS_OPEN_TANAT_REPOINT', CANCEL_LABEL='Select different', BASE=tlb
				(*(*info).winids).errtlb = tlb
				(*(*info).paths).tanat_repointed = 1
			ENDIF ELSE CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
		ENDIF ELSE CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
	ENDIF ELSE BEGIN												; TANAT has not been compiled, but perhaps resides where it should
		IF (tanat_count EQ 1) THEN BEGIN									; Correct TANAT can be found
			dir_tanat = STRMID((*(*info).paths).dir_tanat,0,STRPOS((*(*info).paths).dir_tanat,'/',/REVERSE_SEARCH))
			IF (STRPOS(!PATH,dir_tanat) NE 0) THEN !PATH = dir_tanat+':'+!PATH				; Make sure that that version will be compiled
			RESOLVE_ROUTINE, 'TANAT'
			CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
		ENDIF ELSE BEGIN											; TANAT has not been compiled, nor can be found where it should: have user point to local copy or abort	
			CRISPEX_WINDOW_OK, event,'ERROR!','CRISPEX could not find TANAT at the expected location','('+(*(*info).paths).dir_aux+').',$
				'Press OK to point CRISPEX to a local copy of TANAT, or cancel to abort.',OK_EVENT='CRISPEX_RESTORE_LOOPS_OPEN_TANAT_REPOINT', CANCEL_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
			(*(*info).winids).errtlb = tlb
		ENDELSE
	ENDELSE 
END

PRO CRISPEX_RESTORE_LOOPS_OPEN_TANAT_REPOINT, event				
; Handles (re)pointing CRISPEX to the correct copy of TANAT
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_OPEN_TANAT_REPOINT'
	file_tanat = DIALOG_PICKFILE(PATH='~/', FILTER='*.pro', /FIX_FILTER, /MUST_EXIST, TITLE='CRISPEX'+(*(*info).sesparams).instance_label+': Choose copy of TANAT to compile')
	(*(*info).paths).dir_tanat = STRMID(file_tanat,0,STRPOS(file_tanat,'/',/REVERSE_SEARCH))
	!PATH = (*(*info).paths).dir_tanat+':'+!PATH
	RESOLVE_ROUTINE,'TANAT'
	(*(*info).paths).tanat_repointed = 1
	WIDGET_CONTROL, (*(*info).winids).root, SET_UVALUE = info
	WIDGET_CONTROL, (*(*info).winids).errtlb,/DESTROY
	(*(*info).winids).errtlb = 0
	event.TOP = (*(*info).winids).root
	CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
END

PRO CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN, event
; Handles the actual loading TANAT to analyse the selected loop slice/slab
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RESTORE_LOOPS_OPEN_TANAT_OPEN'
	IF ((*(*info).winids).errtlb GT 0) THEN BEGIN
		WIDGET_CONTROL, (*(*info).winids).root, SET_UVALUE = info
		WIDGET_CONTROL, (*(*info).winids).errtlb,/DESTROY
		(*(*info).winids).errtlb = 0
		event.TOP = (*(*info).winids).root
	ENDIF
	file_csav = DIALOG_PICKFILE(PATH=(*(*info).paths).ipath, FILTER='*.csav', /FIX_FILTER, /MUST_EXIST, TITLE='CRISPEX'+(*(*info).sesparams).instance_label+': Choose CSAV file to load in TANAT')
	IF (STRLEN(file_csav) GT 0) THEN TANAT, file_csav, ASECPIX = (*(*info).meas).arcsecpix, DT = (*(*info).plotaxes).dt
END


;================================================================================= RETRIEVE FROM DETECTIONS PROCEDURES
PRO CRISPEX_RETRIEVE_DET_FILE_MENU, event, set_but_array, DETFILENAME=detfilename, NO_DRAW=no_draw
; Opens the retrieved detections menu and reads in the data from the detection file
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_FILE_MENU'
	IF (N_ELEMENTS(detfilename) NE 1) THEN (*(*info).detparams).detfilename = DIALOG_PICKFILE(/READ, /MUST_EXIST, TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Select detection file')
	IF ((*(*info).detparams).detfilename EQ '') THEN BEGIN
		(*(*info).loopswitch).retrieve_detfile = 0
		RETURN
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).det_file_but, SENSITIVE = 0
		WIDGET_CONTROL, /HOURGLASS
		(*(*info).loopswitch).retrieve_detfile = 1
		RESTORE,(*(*info).detparams).detfilename
		(*(*info).detparams).nr_dets = ntot
		IF (N_ELEMENTS(set_but_array) NE (*(*info).detparams).nr_dets) THEN *(*(*info).detparams).sel_dets = INTARR((*(*info).detparams).nr_dets)
		eventval = INDGEN((*(*info).detparams).nr_dets)
		base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Save from detection file', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
		disp = WIDGET_BASE(base, /COLUMN)
		disp2 = WIDGET_BASE(disp, /COLUMN, /FRAME)
		sel_allnone = WIDGET_BASE(disp2, /ROW)
		sel_allnone_lab = WIDGET_LABEL(sel_allnone, VALUE = 'Select:', /ALIGN_LEFT)
		sel_allnone_buts = WIDGET_BASE(sel_allnone, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsdet).sel_all = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'All', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_SEL_ALL')
		(*(*info).ctrlsdet).sel_none = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'None', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_SEL_NONE')
		WIDGET_CONTROL, (*(*info).ctrlsdet).sel_none, /SET_BUTTON
		disp_txt = REPLICATE('Display time slice D', (*(*info).detparams).nr_dets)
		list_values = ['Do not display any slice', disp_txt+STRTRIM(eventval,2)]
		(*(*info).ctrlsdet).disp_list = WIDGET_COMBOBOX(sel_allnone, VALUE = list_values, /DYNAMIC_RESIZE, EVENT_PRO = 'CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB')
		overlays = WIDGET_BASE(disp2,/ROW)
		overlay_lab = WIDGET_LABEL(overlays, VALUE = 'Overlay:', /ALIGN_LEFT)
		overlay_buts = WIDGET_BASE(overlays, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsdet).overlay_all = WIDGET_BUTTON(overlay_buts, VALUE = 'All detections', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_OVERLAY_ALL')
		WIDGET_CONTROL, (*(*info).ctrlsdet).overlay_all, /SET_BUTTON
		(*(*info).ctrlsdet).overlay_sel = WIDGET_BUTTON(overlay_buts, VALUE = 'Selected detections')
		IF ((*(*info).detparams).nr_dets GT 15) THEN sel_buts = WIDGET_BASE(disp2, /ROW, Y_SCROLL_SIZE = 425) ELSE sel_buts = WIDGET_BASE(disp2, /ROW)
		ncols = 5
		nlines = CEIL((*(*info).detparams).nr_dets / FLOAT(ncols))
		sel_buts_col_0 = WIDGET_BASE(sel_buts, /NONEXCLUSIVE, /COLUMN)
		sel_buts_col_1 = WIDGET_BASE(sel_buts, /NONEXCLUSIVE, /COLUMN)
		sel_buts_col_2 = WIDGET_BASE(sel_buts, /NONEXCLUSIVE, /COLUMN)
		sel_buts_col_3 = WIDGET_BASE(sel_buts, /NONEXCLUSIVE, /COLUMN)
		sel_buts_col_4 = WIDGET_BASE(sel_buts, /NONEXCLUSIVE, /COLUMN)
		FOR i=0,(*(*info).detparams).nr_dets-1 DO BEGIN
			IF (N_ELEMENTS(set_but_array) NE (*(*info).detparams).nr_dets) THEN (*(*(*info).detparams).sel_dets)[i] = 0
			col = FLOOR(i/nlines)
			name = 'det_sel_but_'+STRTRIM(i,2)
			but_val = 'D'+STRTRIM(i,2)
			IF (col EQ 0) THEN act_sel_buts_col = sel_buts_col_0
			IF (col EQ 1) THEN act_sel_buts_col = sel_buts_col_1
			IF (col EQ 2) THEN act_sel_buts_col = sel_buts_col_2
			IF (col EQ 3) THEN act_sel_buts_col = sel_buts_col_3
			IF (col EQ 4) THEN act_sel_buts_col = sel_buts_col_4
			sel_but = WIDGET_BUTTON(act_sel_buts_col, VALUE = but_val, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_RETRIEVE_DET_MENU_EVENT', UNAME = name) 
			IF (N_ELEMENTS(set_but_array) EQ (*(*info).detparams).nr_dets) THEN WIDGET_CONTROL, sel_but, SET_BUTTON = set_but_array[i]
			npoints = (SIZE((*detections[i]).x))[1]
			maxwidth = (SIZE((*detections[i]).x))[2]
			stdwidth = CEIL(maxwidth/2.)
			IF ~KEYWORD_SET(NO_DRAW) THEN (*(*info).detparams).width = stdwidth
			(*(*info).detparams).mid = CEIL(maxwidth/2.)-1
			IF (i EQ 0) THEN BEGIN
				*(*(*info).detparams).t_restored = [(*detections[i]).t]
				*(*(*info).detparams).xlp_restored = [(*detections[i]).x[0,2],(*detections[i]).x[npoints-1,2]]
				*(*(*info).detparams).ylp_restored = [(*detections[i]).y[0,2],(*detections[i]).y[npoints-1,2]]
				*(*(*info).detparams).xlr_restored = [(*detections[i]).x]
				*(*(*info).detparams).ylr_restored = [(*detections[i]).y]
				*(*(*info).detparams).lpsizes = [0, LONG(2)]
				*(*(*info).detparams).lrsizes = [0, LONG(npoints)]
			ENDIF ELSE BEGIN
				*(*(*info).detparams).t_restored = [*(*(*info).detparams).t_restored, (*detections[i]).t]
				*(*(*info).detparams).xlp_restored = [*(*(*info).detparams).xlp_restored, (*detections[i]).x[0,2],(*detections[i]).x[npoints-1,2]]
				*(*(*info).detparams).ylp_restored = [*(*(*info).detparams).ylp_restored, (*detections[i]).y[0,2],(*detections[i]).y[npoints-1,2]]
				*(*(*info).detparams).xlr_restored = [*(*(*info).detparams).xlr_restored, (*detections[i]).x]
				*(*(*info).detparams).ylr_restored = [*(*(*info).detparams).ylr_restored, (*detections[i]).y]
				*(*(*info).detparams).lpsizes = [*(*(*info).detparams).lpsizes, (*(*(*info).detparams).lpsizes)[i]+LONG(2)]
				*(*(*info).detparams).lrsizes = [*(*(*info).detparams).lrsizes, (*(*(*info).detparams).lrsizes)[i]+LONG(npoints)]
			ENDELSE
		ENDFOR			
		time_base = WIDGET_BASE(disp, /ROW, /FRAME)
		dtmin_label = WIDGET_LABEL(time_base, VALUE = 'Delta t down:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).dtmin_text = WIDGET_TEXT(time_base, VALUE = STRTRIM((*(*info).detparams).delta_t_dn,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_DELTA_T_DN')
		dtmax_label = WIDGET_LABEL(time_base, VALUE = 'Delta t up:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).dtmax_text = WIDGET_TEXT(time_base, VALUE = STRTRIM((*(*info).detparams).delta_t_up,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_DELTA_T_UP')
		(*(*info).ctrlsdet).width_slider = WIDGET_SLIDER(disp, TITLE = 'Detection loop width', MIN = 1, MAX = maxwidth, VALUE = (*(*info).detparams).width, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_WIDTH')
		spectral_pos = WIDGET_BASE(disp, /COLUMN, /FRAME)
		spectral_pos_base = WIDGET_BASE(spectral_pos, /ROW)
		spectral_pos_label = WIDGET_LABEL(spectral_pos_base, VALUE = 'Spectral positions:', /ALIGN_LEFT)
		spectral_pos_buts_base = WIDGET_BASE(spectral_pos_base, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsdet).all_pos = WIDGET_BUTTON(spectral_pos_buts_base, VALUE = 'All', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_ALL_POS', /NO_RELEASE)
		WIDGET_CONTROL, (*(*info).ctrlsdet).all_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 1)
		(*(*info).ctrlsdet).saved_pos = WIDGET_BUTTON(spectral_pos_buts_base, VALUE = 'Current', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_CUR_POS', /NO_RELEASE)
		WIDGET_CONTROL, (*(*info).ctrlsdet).saved_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 2)
		(*(*info).ctrlsdet).sel_range_pos = WIDGET_BUTTON(spectral_pos_buts_base, VALUE = 'Selected range', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_SEL_RANGE_POS', /NO_RELEASE)
		WIDGET_CONTROL, (*(*info).ctrlsdet).sel_range_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 3)
		range_pos_base = WIDGET_BASE(spectral_pos, /ROW)
		main_label = WIDGET_LABEL(range_pos_base, VALUE = 'Main:', /ALIGN_LEFT)
		dlpmin_label = WIDGET_LABEL(range_pos_base, VALUE = 'Lower lp-index:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).dlpmin_text = WIDGET_TEXT(range_pos_base, VALUE = STRTRIM((*(*info).detparams).lp_dn,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_LP_DN', SENSITIVE = 0)
		dlpmax_label = WIDGET_LABEL(range_pos_base, VALUE = 'Upper lp-index:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).dlpmax_text = WIDGET_TEXT(range_pos_base, VALUE = STRTRIM((*(*info).detparams).lp_up,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_LP_UP', SENSITIVE = 0)
		refrange_pos_base = WIDGET_BASE(spectral_pos, /ROW)
		ref_label = WIDGET_LABEL(refrange_pos_base, VALUE = 'Reference:', /ALIGN_LEFT)
		refdlpmin_label = WIDGET_LABEL(refrange_pos_base, VALUE = 'Lower lp-index:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).refdlpmin_text = WIDGET_TEXT(refrange_pos_base, VALUE = STRTRIM((*(*info).detparams).lp_ref_dn,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_LP_REF_DN', SENSITIVE = 0)
		refdlpmax_label = WIDGET_LABEL(refrange_pos_base, VALUE = 'Upper lp-index:', /ALIGN_LEFT)
		(*(*info).ctrlsdet).refdlpmax_text = WIDGET_TEXT(refrange_pos_base, VALUE = STRTRIM((*(*info).detparams).lp_ref_up,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_RETRIEVE_DET_LP_REF_UP', SENSITIVE = 0)
		save_cube = WIDGET_BASE(disp, /ROW, /FRAME)
		save_cube_lab = WIDGET_LABEL(save_cube, VALUE = 'Save from: ')
		save_cube_but = WIDGET_BASE(save_cube, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsdet).save_imonly = WIDGET_BUTTON(save_cube_but, VALUE = 'IMCUBE only', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_IMCUBE_ONLY', /NO_RELEASE)
		WIDGET_CONTROL, (*(*info).ctrlsdet).save_imonly, /SET_BUTTON
		refsens = (*(*info).winswitch).showref AND ((*(*info).dataparams).refnt GT 1) 
		(*(*info).ctrlsdet).save_refonly = WIDGET_BUTTON(save_cube_but, VALUE = 'REFCUBE only', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_REFCUBE_ONLY', SENSITIVE = refsens, /NO_RELEASE)
		(*(*info).ctrlsdet).save_imref = WIDGET_BUTTON(save_cube_but, VALUE = 'both cubes', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_IMREF', SENSITIVE = refsens, /NO_RELEASE)
		dec_buts = WIDGET_BASE(disp, /ROW, /ALIGN_CENTER)
		change_path_but = WIDGET_BUTTON(dec_buts, VALUE = 'Change path', EVENT_PRO = 'CRISPEX_SAVE_SET_OPATH')
		(*(*info).ctrlsdet).get_dets = WIDGET_BUTTON(dec_buts, VALUE = 'Save selected detection(s)', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_SEL_LOOPS', SENSITIVE = 0)
		cancel = WIDGET_BUTTON(dec_buts, VALUE = 'Cancel', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_CANCEL')
		WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).lsxoffset
		WIDGET_CONTROL, base, SET_UVALUE = info
		XMANAGER, 'CRISPEX', base, /NO_BLOCK
		(*(*info).winids).detsavetlb = base
	ENDELSE
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopswitch).retrieve_detfile,(*(*info).winids).detsavetlb,(*(*info).detparams).nr_dets],$
		labels=['Retrieving detections','detsavetlb','Retrieved detections']
END

PRO CRISPEX_RETRIEVE_DET_SEL_ALL, event
; Handles the selection of all detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_SEL_ALL'
	*(*(*info).detparams).sel_dets = REPLICATE(1,(*(*info).detparams).nr_dets)
	WIDGET_CONTROL, (*(*info).ctrlsdet).get_dets, SENSITIVE = 1
	FOR i=0,(*(*info).detparams).nr_dets-1 DO BEGIN
		name = 'det_sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 1
	ENDFOR		
	IF ((*(*info).overlayswitch).det_overlay_all EQ 0) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_RETRIEVE_DET_SEL_NONE, event
; Handles the selection of no detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_SEL_NONE'
	*(*(*info).detparams).sel_dets = REPLICATE(0,(*(*info).detparams).nr_dets)
	WIDGET_CONTROL, (*(*info).ctrlsdet).get_dets, SENSITIVE = 0
	FOR i=0,(*(*info).detparams).nr_dets-1 DO BEGIN
		name = 'det_sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 0
	ENDFOR		
	IF ((*(*info).overlayswitch).det_overlay_all EQ 0) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_RETRIEVE_DET_MENU_EVENT, event
; Handles the selection of a certain detection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_MENU_EVENT', /IGNORE_LAST
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).detparams).sel_dets)[eventval] = ( (*(*(*info).detparams).sel_dets)[eventval] EQ 0) 
	condition = WHERE(*(*(*info).detparams).sel_dets EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlsdet).get_dets, SENSITIVE = ((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))
	WIDGET_CONTROL, (*(*info).ctrlsdet).sel_all, SET_BUTTON = (((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).detparams).nr_dets))
	WIDGET_CONTROL, (*(*info).ctrlsdet).sel_none, SET_BUTTON = ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1)
	IF ((*(*info).overlayswitch).det_overlay_all EQ 0) THEN CRISPEX_DRAW, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).detparams).sel_dets)[eventval],$
		(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).restoreparams).cfilecount)),ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1),$
		N_ELEMENTS(condition)-(TOTAL(condition) EQ -1)], labels=['Detection ID','Detection selected','All selected','None selected','Total selected']

END

PRO CRISPEX_RETRIEVE_DET_OVERLAY_ALL, event
; Enables the overlay of all or no detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_OVERLAY_ALL'
	(*(*info).overlayswitch).det_overlay_all = event.SELECT
	CRISPEX_DRAW, event
END

PRO CRISPEX_RETRIEVE_DET_WIDTH, event
; Handles the width of the retrieved detection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_WIDTH'
	prev_width = (*(*info).detparams).width
	set_width = event.VALUE
	IF (set_width NE prev_width) THEN BEGIN
		WIDGET_CONTROL, /HOURGLASS
		IF (FLOOR(set_width/2.) EQ set_width/2.) THEN BEGIN
			IF (prev_width GT set_width) THEN set_width = set_width-1 ELSE set_width = set_width+1
		ENDIF
		(*(*info).detparams).width = set_width
		WIDGET_CONTROL, (*(*info).ctrlsdet).width_slider, SET_VALUE = (*(*info).detparams).width
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [prev_width,set_width],labels=['Previous width','New width']
		IF (*(*info).winswitch).showretrdet THEN BEGIN
			CRISPEX_RETRIEVE_DET_GET_SLICE, event 
			CRISPEX_UPDATE_LP, event
			CRISPEX_DISPRANGE_T_RANGE, event
			WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text,SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2),SENSITIVE = 0
			WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text,SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2),SENSITIVE = 0 
		ENDIF ELSE CRISPEX_DRAW, event
	ENDIF ELSE RETURN
END

PRO CRISPEX_RETRIEVE_DET_ALL_POS, event
; Enables the retreival of detections at all spectral positions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_ALL_POS'
	(*(*info).savswitch).pos_dets = 1
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = 0
	IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = 0
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).pos_dets],labels=['Saving spectral position setting']
END

PRO CRISPEX_RETRIEVE_DET_CUR_POS, event
; Enables the retreival of detections at the current spectral position only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_CUR_POS'
	(*(*info).savswitch).pos_dets = 2
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = 0
	IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = 0
		WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = 0
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).pos_dets],labels=['Saving spectral position setting']
END
	
PRO CRISPEX_RETRIEVE_DET_SEL_RANGE_POS, event
; Enables the retreival of detections at a selected range of spectral positions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_SEL_RANGE_POS'
	(*(*info).savswitch).pos_dets = 3
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = ((*(*info).savswitch).det_imref_only NE 2)
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = ((*(*info).savswitch).det_imref_only NE 2)
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = ((*(*info).savswitch).det_imref_only GE 2)
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = ((*(*info).savswitch).det_imref_only GE 2)
	CRISPEX_RETRIEVE_DET_LP_DN, event
	CRISPEX_RETRIEVE_DET_LP_UP, event
	CRISPEX_RETRIEVE_DET_LP_REF_DN, event
	CRISPEX_RETRIEVE_DET_LP_REF_UP, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).pos_dets],labels=['Saving spectral position setting']
END

PRO CRISPEX_RETRIEVE_DET_IMCUBE_ONLY, event
; Enables the retreival of detections from the image cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_IMCUBE_ONLY'
	(*(*info).savswitch).det_imref_only = 1
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).det_imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_DET_REFCUBE_ONLY, event
; Enables the retreival of detections from the reference cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_REFCUBE_ONLY'
	(*(*info).savswitch).det_imref_only = 2
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).det_imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_DET_IMREF, event
; Enables the retreival of detections from both the image and the reference cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_IMREF'
	(*(*info).savswitch).det_imref_only = 3
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).det_imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_DET_SEL_LOOPS, event
; Opens the warning windows giving the saving time estimates of the chosen saving procedure, intermediate step towards saving detection loopslabs
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_SEL_LOOPS'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
		*(*(*info).detparams).sel_loops = WHERE(*(*(*info).detparams).sel_dets EQ 1)
		(*(*info).detparams).nr_sel_loops = N_ELEMENTS(*(*(*info).detparams).sel_loops)
		IF ((*(*info).savswitch).pos_dets EQ 2) THEN BEGIN
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).detparams).nr_sel_loops-1 DO BEGIN
				w_loop_pts = (*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]+1]-(*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]]
				sub_time = (*(*info).feedbparams).estimate_time * w_loop_pts/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).detparams).width * ((*(*info).detparams).delta_t_up + $
					(*(*info).detparams).delta_t_dn + 1) * CEIL((*(*info).savswitch).det_imref_only/2.)
				time = time + sub_time
			ENDFOR
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			(*(*info).savswitch).cont = 4
			CRISPEX_SAVE_WARNING_YESNO, event,STRTRIM((*(*info).detparams).nr_sel_loops,2)+' detection(s) selected. Saving the exact loop slice(s) may take' ,$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+'. Do you wish to continue saving?', OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF ELSE IF ((*(*info).savswitch).pos_dets EQ 1) THEN BEGIN
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).detparams).nr_sel_loops-1 DO BEGIN
				w_loop_pts = (*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]+1]-(*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]]
				sub_im_time = (*(*info).feedbparams).estimate_time * w_loop_pts/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dataparams).nlp * (*(*info).detparams).width * $
					((*(*info).detparams).delta_t_up + (*(*info).detparams).delta_t_dn + 1)
				IF ((*(*info).dataparams).refnt EQ (*(*info).dataparams).nt) THEN sub_ref_time = sub_im_time/(*(*info).dataparams).nlp ELSE sub_ref_time = sub_im_time
				IF ((*(*info).savswitch).det_imref_only EQ 1) THEN sub_time = sub_im_time ELSE IF ((*(*info).savswitch).det_imref_only EQ 2) THEN sub_time = sub_ref_time ELSE sub_time = sub_im_time + sub_ref_time
				time += sub_time
			ENDFOR
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			(*(*info).savswitch).cont = 5
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).detparams).nr_sel_loops,2)+' detection(s) selected. Saving the exact loop slab(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+'. Do you wish to continue saving?', OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF ELSE IF ((*(*info).savswitch).pos_dets EQ 3) THEN BEGIN
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).detparams).nr_sel_loops-1 DO BEGIN
				w_loop_pts = (*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]+1]-(*(*(*info).detparams).lrsizes)[(*(*(*info).detparams).sel_loops)[i]]
				sub_im_time = (*(*info).feedbparams).estimate_time * w_loop_pts/FLOAT((*(*info).feedbparams).estimate_lx) * ((*(*info).detparams).lp_up-(*(*info).detparams).lp_dn+1) * $
					(*(*info).detparams).width * ((*(*info).detparams).delta_t_up + (*(*info).detparams).delta_t_dn + 1)
				IF ((*(*info).dataparams).refnt EQ (*(*info).dataparams).nt) THEN sub_ref_time = sub_im_time/((*(*info).detparams).lp_up-(*(*info).detparams).lp_dn+1) ELSE sub_ref_time = sub_im_time
				IF ((*(*info).savswitch).det_imref_only EQ 1) THEN sub_time = sub_im_time ELSE IF ((*(*info).savswitch).det_imref_only EQ 2) THEN sub_time = sub_ref_time ELSE sub_time = sub_im_time + sub_ref_time
				time += sub_time
			ENDFOR
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			(*(*info).savswitch).cont = 5
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).detparams).nr_sel_loops,2)+' detection(s) selected. Saving the exact loop slab(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+'. Do you wish to continue saving?', OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF ELSE RETURN
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).detparams).nr_sel_loops,(*(*info).savswitch).cont],labels=['Number of detections','Saving procedure']
	ENDIF
END

PRO CRISPEX_RETRIEVE_DET_DELTA_T_DN, event
; Handles the change in delta_t down from the detection framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_DELTA_T_DN'
	WIDGET_CONTROL, (*(*info).ctrlsdet).dtmin_text, GET_VALUE = textvalue
	(*(*info).detparams).delta_t_dn = FLOAT(LONG(FLOAT(textvalue[0])))
	IF (((*(*info).detparams).delta_t_dn GE 0) AND ((*(*info).detparams).delta_t_dn LT 1)) THEN (*(*info).detparams).delta_t_dn = 1 
	IF ((*(*info).detparams).delta_t_dn LT 0) THEN (*(*info).detparams).delta_t_dn = FLOAT(LONG(ABS((*(*info).detparams).delta_t_dn)))
	IF ((*(*info).detparams).delta_t_dn GT (*(*info).dispparams).t_last) THEN (*(*info).detparams).delta_t_dn = (*(*info).dispparams).t_last
	WIDGET_CONTROL, (*(*info).ctrlsdet).dtmin_text, SET_VALUE = STRTRIM((*(*info).detparams).delta_t_dn, 2)
	CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_RETRIEVE_DET_DELTA_T_UP, event
; Handles the change in delta_t up from the detection framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_DELTA_T_UP'
	WIDGET_CONTROL, (*(*info).ctrlsdet).dtmax_text, GET_VALUE = textvalue
	(*(*info).detparams).delta_t_up = FLOAT(LONG(FLOAT(textvalue[0])))
	IF (((*(*info).detparams).delta_t_up GE 0) AND ((*(*info).detparams).delta_t_up LT 1)) THEN (*(*info).detparams).delta_t_up = 1 
	IF ((*(*info).detparams).delta_t_up LT 0) THEN (*(*info).detparams).delta_t_up = FLOAT(LONG(ABS((*(*info).detparams).delta_t_up)))
	IF ((*(*info).detparams).delta_t_up GT (*(*info).dispparams).t_last) THEN (*(*info).detparams).delta_t_up = (*(*info).dispparams).t_last
	WIDGET_CONTROL, (*(*info).ctrlsdet).dtmax_text, SET_VALUE = STRTRIM((*(*info).detparams).delta_t_up, 2)
	CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_RETRIEVE_DET_LP_DN, event
; Handles the change in lower spectral boundary for the detection extraction
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_LP_DN'
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, GET_VALUE = textvalue
	(*(*info).detparams).lp_dn = FLOAT(textvalue[0])
	IF ((*(*info).detparams).lp_dn LT (*(*info).dispparams).lp_first) THEN (*(*info).detparams).lp_dn = (*(*info).dispparams).lp_first ELSE IF ((*(*info).detparams).lp_dn GE (*(*info).detparams).lp_up) THEN $
		(*(*info).detparams).lp_dn = (*(*info).detparams).lp_up-1
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SET_VALUE = STRTRIM((*(*info).detparams).lp_dn,2)
	(*(*info).dispparams).lp_low = (*(*info).detparams).lp_dn
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2)
	CRISPEX_DISPRANGE_LP_RANGE, event
END

PRO CRISPEX_RETRIEVE_DET_LP_UP, event
; Handles the change in upper spectral boundary for the detection extraction
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_LP_UP'
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, GET_VALUE = textvalue
	(*(*info).detparams).lp_up = FLOAT(textvalue[0])
	IF ((*(*info).detparams).lp_up GT (*(*info).dispparams).lp_last) THEN (*(*info).detparams).lp_up = (*(*info).dispparams).lp_last ELSE IF ((*(*info).detparams).lp_up LE (*(*info).detparams).lp_dn) THEN $
		(*(*info).detparams).lp_up = (*(*info).detparams).lp_dn+1
	WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SET_VALUE = STRTRIM((*(*info).detparams).lp_up,2)
	(*(*info).dispparams).lp_upp = (*(*info).detparams).lp_up
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2)
	CRISPEX_DISPRANGE_LP_RANGE, event
END

PRO CRISPEX_RETRIEVE_DET_LP_REF_DN, event
; Handles the change in lower spectral boundary for the detection extraction
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_LP_REF_DN'
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, GET_VALUE = textvalue
	(*(*info).detparams).lp_ref_dn = FLOAT(textvalue[0])
	IF ((*(*info).detparams).lp_ref_dn LT (*(*info).dispparams).lp_ref_first) THEN (*(*info).detparams).lp_ref_dn = (*(*info).dispparams).lp_ref_first ELSE $
		IF ((*(*info).detparams).lp_ref_dn GE (*(*info).detparams).lp_ref_up) THEN (*(*info).detparams).lp_ref_dn = (*(*info).detparams).lp_ref_up-1
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SET_VALUE = STRTRIM((*(*info).detparams).lp_ref_dn,2)
	(*(*info).dispparams).lp_ref_low = (*(*info).detparams).lp_ref_dn
END

PRO CRISPEX_RETRIEVE_DET_LP_REF_UP, event
; Handles the change in upper spectral boundary for the detection extraction
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_LP_REF_UP'
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, GET_VALUE = textvalue
	(*(*info).detparams).lp_ref_up = FLOAT(textvalue[0])
	IF ((*(*info).detparams).lp_ref_up GT (*(*info).dispparams).lp_ref_last) THEN (*(*info).detparams).lp_ref_up = (*(*info).dispparams).lp_ref_last ELSE $
		IF ((*(*info).detparams).lp_ref_up LE (*(*info).detparams).lp_ref_dn) THEN (*(*info).detparams).lp_ref_up = (*(*info).detparams).lp_ref_dn+1
	WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SET_VALUE = STRTRIM((*(*info).detparams).lp_ref_up,2)
	(*(*info).dispparams).lp_ref_upp = (*(*info).detparams).lp_ref_up
END

PRO CRISPEX_RETRIEVE_DET_GET_SLICE, event
; Handles the in-program display of the detection loopslab
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_GET_SLICE'
	tmp_loopslab = 0
	delta = FLOOR((*(*info).detparams).width/2.)
	FOR k=((*(*info).detparams).mid-delta),((*(*info).detparams).mid+delta) DO BEGIN
		xlp_det = (*(*(*info).detparams).xlp_restored)[(0+(*(*(*info).detparams).lpsizes)[(*(*info).detparams).idx]):((*(*(*info).detparams).lpsizes)[(*(*info).detparams).idx+1]-1)]
		ylp_det = (*(*(*info).detparams).ylp_restored)[(0+(*(*(*info).detparams).lpsizes)[(*(*info).detparams).idx]):((*(*(*info).detparams).lpsizes)[(*(*info).detparams).idx+1]-1)]
		xlr_det = (*(*(*info).detparams).xlr_restored)[(0+(*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx]):((*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx+1]-1), k]
		ylr_det = (*(*(*info).detparams).ylr_restored)[(0+(*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx]):((*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx+1]-1), k]
		IF (*(*info).dispswitch).exts THEN BEGIN
			WIDGET_CONTROL,/HOURGLASS
			w_lpts = INDGEN((*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx+1] - (*(*(*info).detparams).lrsizes)[(*(*info).detparams).idx])
			lp_orig = (*(*info).dataparams).lp
			FOR i=(*(*info).detparams).lp_dn,(*(*info).detparams).lp_up DO BEGIN
				(*(*info).dataparams).lp = i
				CRISPEX_LOOP_GET_EXACT_SLICE, event, *(*(*info).data).imagedata, xlr_det, ylr_det, xlp_det, ylp_det, w_lpts, det_loopslice, det_crossloc, det_loopsize, /im
				IF (i EQ (*(*info).detparams).lp_dn) THEN loopslab = det_loopslice ELSE loopslab = [[[loopslab]], [[det_loopslice]]]
				feedback_text = 'Detection '+STRTRIM((*(*info).detparams).idx,2)+', position '+STRTRIM(k,2)+'/'+STRTRIM((*(*info).detparams).width,2)+': '+$
					STRTRIM(i-(*(*info).detparams).lp_dn+1,2)+'/'+STRTRIM((*(*info).detparams).lp_up-(*(*info).detparams).lp_dn+1,2)+' slices extracted...'
				CRISPEX_UPDATE_USER_FEEDBACK, event, title='Retrieved detection loop slice', var=(i-(*(*info).detparams).lp_dn+1)+k, minvar=((*(*info).detparams).mid-delta)+1, $
					maxvar=((*(*info).detparams).lp_up-(*(*info).detparams).lp_dn+1), feedback_text=feedback_text
			ENDFOR
			(*(*info).dataparams).lp = lp_orig
			IF (k EQ (*(*info).detparams).mid) THEN BEGIN
				crossloc = det_crossloc
				loopsize = det_loopsize
			ENDIF
			CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event
		ENDIF ELSE BEGIN
			CRISPEX_LOOP_GET_APPROX_SLAB, event, xlr_det, ylr_det, xlp_det, xlp_det, loopslab, crossloc, loopsize
		ENDELSE
		tmp_loopslab = tmp_loopslab + loopslab
	ENDFOR
	*(*(*info).loopsdata).det_loopslab = tmp_loopslab/(*(*info).detparams).width
	*(*(*info).loopsdata).det_crossloc = crossloc
	(*(*info).loopsdata).det_loopsize = loopsize
	t_det = (*(*(*info).detparams).t_restored)[(*(*info).detparams).idx]
	(*(*info).dispparams).t_low = t_det - (*(*info).detparams).delta_t_dn > (*(*info).dispparams).t_first
	(*(*info).dispparams).t_upp = t_det + (*(*info).detparams).delta_t_up < (*(*info).dispparams).t_last
	(*(*info).dataparams).t = t_det
END

PRO CRISPEX_RETRIEVE_DET_CANCEL, event
; Handles the closing of the detection file menu
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_DET_CANCEL'
	(*(*info).loopswitch).retrieve_detfile = 0
	(*(*info).winswitch).showretrdet = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopswitch).retrieve_detfile,(*(*info).winids).detsavetlb], labels=['Retrieving detections','retrdettlb was']
	CRISPEX_DISPRANGE_T_RESET, event
	WIDGET_CONTROL, (*(*info).winids).detsavetlb, /DESTROY
	(*(*info).winids).detsavetlb = 0
	IF ((*(*info).winids).retrdettlb GT 0) THEN WIDGET_CONTROL, (*(*info).winids).retrdettlb, /DESTROY
	(*(*info).winids).retrdettlb = 0
	*(*(*info).detparams).sel_dets = REPLICATE(0,(*(*info).detparams).nr_dets)
	(*(*info).overlayswitch).det_overlay_all = 1
	(*(*info).savswitch).pos_dets = 1
	(*(*info).savswitch).det_imref_only = 1
	(*(*info).detparams).lp_dn = (*(*info).dispparams).lp_first		&	(*(*info).detparams).lp_up = (*(*info).dispparams).lp_last
	(*(*info).detparams).lp_ref_dn = (*(*info).dispparams).lp_ref_first	&	(*(*info).detparams).lp_ref_up = (*(*info).dispparams).lp_ref_last
	WIDGET_CONTROL, (*(*info).ctrlscp).det_file_but, /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SENSITIVE = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SENSITIVE = 1
END

;================================================================================= RETRIEVE FROM SAVED LOOP POINT PROCEDURES
PRO CRISPEX_RETRIEVE_LOOP_MENU, event, set_but_array
; Sets up the retrieved loop points menu and reads in the data from the respective files
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_MENU'
	CRISPEX_FIND_CLSAV, event
	IF ((*(*info).retrparams).clfilecount GT 0) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 0
		(*(*info).loopswitch).retrieve_loops = 1
		filenames = *(*(*info).retrparams).clfiles
		filecount = (*(*info).retrparams).clfilecount
		IF (N_ELEMENTS(set_but_array) NE filecount) THEN *(*(*info).retrparams).sel_loops = INTARR((*(*info).retrparams).clfilecount)
		eventval = INDGEN(filecount)
		base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Save from retrieved loop(s)', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
		disp = WIDGET_BASE(base, /COLUMN)
		disp2 = WIDGET_BASE(disp, /COLUMN, /FRAME)
		sel_allnone = WIDGET_BASE(disp2, /ROW)
		sel_allnone_lab = WIDGET_LABEL(sel_allnone, VALUE = 'Select:', /ALIGN_LEFT)
		sel_allnone_buts = WIDGET_BASE(sel_allnone, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsloop).sel_all = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'All', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_SEL_ALL')
		(*(*info).ctrlsloop).sel_none = WIDGET_BUTTON(sel_allnone_buts, VALUE = 'None', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_SEL_NONE')
		WIDGET_CONTROL, (*(*info).ctrlsloop).sel_none, /SET_BUTTON
		update_filelist = WIDGET_BUTTON(sel_allnone, VALUE = 'Update filelist', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_UPDATE_FILELIST')
		IF (filecount GT 16) THEN sel_buts = WIDGET_BASE(disp2, /COLUMN, /NONEXCLUSIVE, X_SCROLL_SIZE=400, Y_SCROLL_SIZE=425) ELSE sel_buts = WIDGET_BASE(disp2, /COLUMN, /NONEXCLUSIVE)
		FOR i=0,filecount-1 DO BEGIN
			IF (N_ELEMENTS(set_but_array) NE filecount) THEN (*(*(*info).retrparams).sel_loops)[i] = 0
			singlefilename = (STRSPLIT((*(*(*info).retrparams).clfiles)[i], PATH_SEP(), /EXTRACT))[N_ELEMENTS(STRSPLIT((*(*(*info).retrparams).clfiles)[i], PATH_SEP(), /EXTRACT))-1]
			name = 'sel_but_'+STRTRIM(i,2)
			but_val = 'L'+STRTRIM(i,2)+': '+singlefilename
			sel_but = WIDGET_BUTTON(sel_buts, VALUE = but_val, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_MENU_EVENT', UNAME = name) 
			IF (N_ELEMENTS(set_but_array) GT 0) THEN WIDGET_CONTROL, sel_but, SET_BUTTON = set_but_array[i]
			RESTORE,(*(*(*info).retrparams).clfiles)[i]
			IF (i EQ 0) THEN BEGIN
				*(*(*info).retrparams).xlp_restored = REFORM(x_coords[0,*])
				*(*(*info).retrparams).ylp_restored = REFORM(y_coords[0,*])
				*(*(*info).retrparams).xlr_restored = x_loop_pts
				*(*(*info).retrparams).ylr_restored = y_loop_pts
				*(*(*info).retrparams).lpsizes = [0, LONG((SIZE(REFORM(x_coords[0,*])))[1])]
				*(*(*info).retrparams).lrsizes = [0, LONG((SIZE(REFORM(x_loop_pts)))[1])]
			ENDIF ELSE BEGIN
				*(*(*info).retrparams).xlp_restored = [*(*(*info).retrparams).xlp_restored, REFORM(x_coords[0,*])]
				*(*(*info).retrparams).ylp_restored = [*(*(*info).retrparams).ylp_restored, REFORM(y_coords[0,*])]
				*(*(*info).retrparams).xlr_restored = [*(*(*info).retrparams).xlr_restored, x_loop_pts]
				*(*(*info).retrparams).ylr_restored = [*(*(*info).retrparams).ylr_restored, y_loop_pts]
				singlepsize = (SIZE(REFORM(x_coords[0,*])))[1]
				singlersize = (SIZE(REFORM(x_loop_pts)))[1]
				*(*(*info).retrparams).lpsizes = [*(*(*info).retrparams).lpsizes, (*(*(*info).retrparams).lpsizes)[i]+LONG(singlepsize)]
				*(*(*info).retrparams).lrsizes = [*(*(*info).retrparams).lrsizes, (*(*(*info).retrparams).lrsizes)[i]+LONG(singlersize)]
			ENDELSE
		ENDFOR
		spectral_pos = WIDGET_BASE(disp, /ROW, /FRAME, /EXCLUSIVE)
		(*(*info).ctrlsloop).all_pos = WIDGET_BUTTON(spectral_pos, VALUE = 'At all spectral positions', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_ALL_POS')
		WIDGET_CONTROL, (*(*info).ctrlsloop).all_pos, /SET_BUTTON
		(*(*info).ctrlsloop).saved_pos = WIDGET_BUTTON(spectral_pos, VALUE = 'At saved spectral position')
		del_clsav = WIDGET_BASE(disp, /ROW, /FRAME, /EXCLUSIVE)
		(*(*info).ctrlsloop).del_files = WIDGET_BUTTON(del_clsav, VALUE = 'Delete *clsav file(s) after saving', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_DELETE_CLSAV')
		WIDGET_CONTROL, (*(*info).ctrlsloop).del_files, /SET_BUTTON
		(*(*info).ctrlsloop).keep_files = WIDGET_BUTTON(del_clsav, VALUE = 'Keep *clsav file(s) after saving')
		save_cube = WIDGET_BASE(disp, /ROW, /FRAME)
		save_cube_lab = WIDGET_LABEL(save_cube, VALUE = 'Save from: ')
		save_cube_but = WIDGET_BASE(save_cube, /ROW, /EXCLUSIVE)
		(*(*info).ctrlsloop).save_imonly = WIDGET_BUTTON(save_cube_but, VALUE = 'IMCUBE only', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_IMCUBE_ONLY', /NO_RELEASE)
		WIDGET_CONTROL, (*(*info).ctrlsloop).save_imonly, /SET_BUTTON
		refsens = (*(*info).winswitch).showref AND ((*(*info).dataparams).refnt GT 1)
		(*(*info).ctrlsloop).save_refonly = WIDGET_BUTTON(save_cube_but, VALUE = 'REFCUBE only', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_REFCUBE_ONLY', SENSITIVE = refsens, /NO_RELEASE)
		(*(*info).ctrlsloop).save_imref = WIDGET_BUTTON(save_cube_but, VALUE = 'both cubes', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_IMREF', SENSITIVE = refsens, /NO_RELEASE)
		dec_buts = WIDGET_BASE(disp, COLUMN=3, /GRID_LAYOUT, /ALIGN_CENTER)
		change_path_but = WIDGET_BUTTON(dec_buts, VALUE = 'Change path', EVENT_PRO = 'CRISPEX_SAVE_SET_OPATH')
		cancel = WIDGET_BUTTON(dec_buts, VALUE = 'Cancel', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_MENU_CANCEL')
		get_loops = WIDGET_BUTTON(dec_buts, VALUE = 'Retrieve and save', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_MENU_CONTINUE', SENSITIVE = 0)
		(*(*info).ctrlsloop).get_loops = get_loops
		WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).spxoffset, TLB_SET_YOFFSET = 0
		WIDGET_CONTROL, base, SET_UVALUE = info
		XMANAGER, 'CRISPEX', base, /NO_BLOCK
		(*(*info).winids).savetlb = base
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).savetlb,filecount,(*(*info).loopswitch).restore_loops],labels=['savetlb','Retrieved loops','Was restoring loops']
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','No saved loop points (*.clsav) files found', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
	IF (*(*info).loopswitch).restore_loops THEN BEGIN
		(*(*info).loopswitch).was_restore_loops = (*(*info).loopswitch).restore_loops
		(*(*info).loopswitch).restore_loops = 0
	ENDIF
END

PRO CRISPEX_RETRIEVE_LOOP_MENU_EVENT, event
; Handles the selection of a certain retrieved loop point file
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_MENU_EVENT'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).retrparams).sel_loops)[eventval] = ( (*(*(*info).retrparams).sel_loops)[eventval] EQ 0)
	condition = WHERE(*(*(*info).retrparams).sel_loops EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlsloop).get_loops, SENSITIVE = ((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))
	WIDGET_CONTROL, (*(*info).ctrlsloop).sel_all, SET_BUTTON = (((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).retrparams).clfilecount)) 
	WIDGET_CONTROL, (*(*info).ctrlsloop).sel_none, SET_BUTTON = ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).retrparams).sel_loops)[eventval],$
		(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).retrparams).clfilecount)),ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1),$
		N_ELEMENTS(condition)-(TOTAL(condition) EQ -1)], labels=['Retrieved loop ID','Retreived loop selected','All selected','None selected','Total selected']
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RETRIEVE_LOOP_MENU_CANCEL, event
; Handles the closing of the retrieved loop menu
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_MENU_CANCEL'
	WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 1
	(*(*info).loopswitch).retrieve_loops = 0
	(*(*info).loopswitch).restore_loops = (*(*info).loopswitch).was_restore_loops
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopswitch).retrieve_loops,(*(*info).winids).savetlb,(*(*info).loopswitch).restore_loops], $
		labels=['Retrieving loops','savetlb was','Restoring loops']
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW_IMREF, event
	WIDGET_CONTROL, (*(*info).winids).savetlb, /DESTROY
	(*(*info).winids).savetlb = 0
END

PRO CRISPEX_RETRIEVE_LOOP_MENU_CONTINUE, event
; Opens the warning windows giving the saving time estimates of the chosen saving procedure, intermediate step towards saving selected retrieved loopslabs
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_MENU_CONTINUE'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
		*(*(*info).retrparams).retrieve_files = (*(*(*info).retrparams).clfiles)[WHERE(*(*(*info).retrparams).sel_loops EQ 1)]
		(*(*info).retrparams).retrieve_filecount = N_ELEMENTS(*(*(*info).retrparams).retrieve_files)
		IF ((*(*info).savswitch).all_pos_loops EQ 0) THEN BEGIN
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).retrparams).retrieve_filecount-1 DO BEGIN
				RESTORE,(*(*(*info).retrparams).retrieve_files)[i]
				sub_time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(w_loop_pts)/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dispparams).t_range * CEIL((*(*info).savswitch).imref_only/2.)
				time += sub_time
			ENDFOR
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			(*(*info).savswitch).cont = 1
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).retrparams).retrieve_filecount,2)+' *clsav file(s) selected. Saving the exact loop slice(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+' and equally named (approximated) loop','slice(s) may be overwritten. Do you wish to continue saving?', $
				OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF ELSE IF (*(*info).savswitch).all_pos_loops THEN BEGIN
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).retrparams).retrieve_filecount-1 DO BEGIN
				RESTORE,(*(*(*info).retrparams).retrieve_files)[i]
				sub_im_time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(w_loop_pts)/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dataparams).nlp * (*(*info).dispparams).t_range
				IF ((*(*info).dataparams).refnt EQ (*(*info).dataparams).nt) THEN sub_ref_time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(w_loop_pts)/FLOAT((*(*info).feedbparams).estimate_lx) * $
					(*(*info).dispparams).t_range ELSE sub_ref_time = sub_im_time
				IF ((*(*info).savswitch).imref_only EQ 1) THEN sub_time = sub_im_time ELSE IF ((*(*info).savswitch).imref_only EQ 2) THEN sub_time = sub_ref_time ELSE sub_time = sub_im_time + sub_ref_time
				time += sub_time
			ENDFOR
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			(*(*info).savswitch).cont = 2
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).retrparams).retrieve_filecount,2)+' *clsav file(s) selected. Saving the exact loop slab(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+' and equally named (approximated) loop','slab(s) may be overwritten. Do you wish to continue saving?',$
				OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF ELSE RETURN
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).retrparams).retrieve_filecount,(*(*info).savswitch).cont],labels=['Number of retrieved loops','Saving procedure']
	ENDIF
END

PRO CRISPEX_RETRIEVE_LOOP_SEL_ALL, event
; Handles the selection of all retrieved loop files
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_SEL_ALL'
	*(*(*info).retrparams).sel_loops = REPLICATE(1,(*(*info).retrparams).clfilecount)
	WIDGET_CONTROL, (*(*info).ctrlsloop).get_loops, SENSITIVE = 1
	FOR i=0,(*(*info).retrparams).clfilecount-1 DO BEGIN
		name = 'sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 1
	ENDFOR		
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RETRIEVE_LOOP_SEL_NONE, event
; Handles the selection of no retrieved loop files
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_SEL_NONE'
	*(*(*info).retrparams).sel_loops = REPLICATE(0,(*(*info).retrparams).clfilecount)
	WIDGET_CONTROL, (*(*info).ctrlsloop).get_loops, SENSITIVE = 0
	FOR i=0,(*(*info).retrparams).clfilecount-1 DO BEGIN
		name = 'sel_but_'+STRTRIM(i,2)
		WIDGET_CONTROL, WIDGET_INFO(event.TOP, FIND_BY_UNAME = name), SET_BUTTON = 0
	ENDFOR		
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RETRIEVE_LOOP_UPDATE_FILELIST, event
; Handles the update of the retrieved file list
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_UPDATE_FILELIST'
	WIDGET_CONTROL, (*(*info).winids).root, SET_UVALUE = info
	WIDGET_CONTROL, event.TOP, /DESTROY
	event.TOP = (*(*info).winids).root
	CRISPEX_RETRIEVE_LOOP_MENU, event
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RETRIEVE_LOOP_DELETE_CLSAV, event
; Enables or disables the deletion of CLSAV files after saving the respective loopslabs
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_DELETE_CLSAV'
	(*(*info).savswitch).delete_clsav = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).delete_clsav],labels=['Delete loop path file']
END

PRO CRISPEX_RETRIEVE_LOOP_ALL_POS, event
; Enables the extraction of the retrieved loop at all or only the saved spectral positions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_ALL_POS'
	(*(*info).savswitch).all_pos_loops = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).all_pos_loops],labels=['Saving spectral position setting']
END

PRO CRISPEX_RETRIEVE_LOOP_IMCUBE_ONLY, event
; Enables the retreival of loop paths from the image cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_IMCUBE_ONLY'
	(*(*info).savswitch).imref_only = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_REFCUBE_ONLY, event
; Enables the retreival of loop paths from the reference cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_REFCUBE_ONLY'
	(*(*info).savswitch).imref_only = 2
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_IMREF, event
; Enables the retreival of loop paths from both the image and reference cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_IMREF'
	(*(*info).savswitch).imref_only = 3
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLAB, event
; Opens the warning windows giving the saving time estimate of the procedure, intermediate step towards saving all retrieved loopslabs (i.e. at all spectral positions)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLAB'
	IF ((*(*info).retrparams).clfilecount GT 0) THEN BEGIN
		CRISPEX_SAVE_CHECK_PATH_WRITE, event
		IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
			(*(*info).savparams).lp_orig = (*(*info).dataparams).lp
			(*(*info).retrparams).retrieve_filecount = (*(*info).retrparams).clfilecount
			*(*(*info).retrparams).retrieve_files = *(*(*info).retrparams).clfiles
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).retrparams).retrieve_filecount-1 DO BEGIN
				RESTORE,(*(*(*info).retrparams).retrieve_files)[i]
				sub_time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(w_loop_pts)/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dataparams).nlp * (*(*info).dispparams).t_range
				time += sub_time
			ENDFOR
			(*(*info).savswitch).cont = 2
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).retrparams).retrieve_filecount,2)+' *clsav file(s) found. Saving the exact loop slab(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+' and equally named (approximated) loop','slab(s) may be overwritten. Do you wish to continue saving?',$
				OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
			(*(*info).dataparams).lp = (*(*info).savparams).lp_orig
		ENDIF
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).retrparams).retrieve_filecount,(*(*info).savswitch).cont],labels=['Number of retrieved loops','Saving procedure']
	ENDIF
END

PRO CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLICE, event
; Opens the warning windows giving the saving time estimate of the procedure, intermediate step towards saving all retrieved loopslices (i.e. at saved spectral positions)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLICE'
	CRISPEX_FIND_CLSAV, event
	IF ((*(*info).retrparams).clfilecount GT 0) THEN BEGIN
		CRISPEX_SAVE_CHECK_PATH_WRITE, event
		IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
			(*(*info).retrparams).retrieve_filecount = (*(*info).retrparams).clfilecount
			*(*(*info).retrparams).retrieve_files = *(*(*info).retrparams).clfiles
			IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
				CRISPEX_ESTIMATE_TIME_WINDOW, event
				CRISPEX_ESTIMATE_TIME_CALCULATION, event
			ENDIF
			time = 0
			FOR i=0,(*(*info).retrparams).retrieve_filecount-1 DO BEGIN
				RESTORE,(*(*(*info).retrparams).retrieve_files)[i]
				sub_time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(w_loop_pts)/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dispparams).t_range
				time += sub_time
			ENDFOR
			(*(*info).savswitch).cont = 1
			CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
			CRISPEX_SAVE_WARNING_YESNO, event, STRTRIM((*(*info).retrparams).retrieve_filecount,2)+' *clsav file(s) found. Saving the exact loop slice(s) may take',$
				'up to about '+STRTRIM(CEIL(time/denom),2)+units+' and equally named (approximated) loop','slice(s) may be overwritten. Do you wish to continue saving?',$
				OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
			IF ((*(*info).winswitch).estimate_win EQ 1) THEN BEGIN
				WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
				(*(*info).winids).estimatetlb = 0
				(*(*info).winswitch).estimate_win = 0
			ENDIF
		ENDIF
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).retrparams).retrieve_filecount,(*(*info).savswitch).cont],labels=['Number of retrieved loops','Saving procedure']
	ENDIF
END

;================================================================================= IMAGE SCALING PROCEDURES
PRO CRISPEX_SCALING_AUTO_CONST, event
; Enables the automatic constant scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_AUTO_CONST'
	(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SENSITIVE = 0, SET_VALUE = 255
	WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SENSITIVE = 0, SET_VALUE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).abs_scale_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).rel_scale_but, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling]],$
		labels=['Scaling image','Scaling setting']
	CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_AUTO_SING, event
; Enables the automatic scaling per timestep
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_AUTO_SING'
	(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SENSITIVE = 0, SET_VALUE = 255
	WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SENSITIVE = 0, SET_VALUE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).abs_scale_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).rel_scale_but, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling]],$
		labels=['Scaling image','Scaling setting']
	CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_MAN_FIRST, event										
; Enables the manual scaling based on the first timestep
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_MAN_FIRST'
	(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] = 2
	IF ((*(*info).scaling).imrefscaling EQ 0) THEN BEGIN
		IF ((*(*info).scaling).scalestokes_max AND ((*(*info).dataparams).s GE 1)) THEN BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = MIN(((*(*info).scaling).imagemin)[*,(*(*info).dataparams).s])
			maximum = MAX(((*(*info).scaling).imagemax)[*,(*(*info).dataparams).s])
		ENDIF ELSE BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = ((*(*info).scaling).imagemin)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
			maximum = ((*(*info).scaling).imagemax)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
		ENDELSE
	ENDIF ELSE IF ((*(*info).scaling).imrefscaling EQ 1) THEN BEGIN
		IF ((*(*info).dataparams).refnlp EQ 1) THEN BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = (*(*info).scaling).refmin
			maximum = (*(*info).scaling).refmax
		ENDIF ELSE BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = ((*(*info).scaling).refmin)[(*(*info).dataparams).lp_ref]
			maximum = ((*(*info).scaling).refmax)[(*(*info).dataparams).lp_ref]
		ENDELSE
	ENDIF ELSE IF ((*(*info).scaling).imrefscaling EQ 2) THEN BEGIN
		IF ((*(*info).scaling).scalestokes_max AND ((*(*info).dataparams).s GE 1)) THEN BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = MIN(((*(*info).scaling).dopplermin)[*,(*(*info).dataparams).s])
			maximum = MAX(((*(*info).scaling).dopplermax)[*,(*(*info).dataparams).s])
		ENDIF ELSE BEGIN
			(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = ((*(*info).scaling).dopplermin)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
			maximum = ((*(*info).scaling).dopplermax)[(*(*info).dataparams).lp,(*(*info).dataparams).s]
		ENDELSE
	ENDIF
	(*(*(*info).scaling).scale_range)[(*(*info).scaling).imrefscaling] = maximum - (*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling],maximum],labels=['Scaling image','Scaling setting','Scaling minimum','Scaling maximum']
	CRISPEX_SCALING_SET_SLIDERS, event
	CRISPEX_SCALING_REDRAW, event
	WIDGET_CONTROL, (*(*info).ctrlscp).abs_scale_but, /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).rel_scale_but, /SENSITIVE
END

PRO CRISPEX_SCALING_MAN_CURR, event					
; Enables the manual scaling based on the current timestep
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_MAN_CURR'
	(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] = 3
	IF ((*(*info).scaling).imrefscaling EQ 0) THEN BEGIN
		IF (*(*info).dataswitch).spfile OR (*(*info).dataswitch).onecube THEN imdata = (*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
			(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp] ELSE imdata = (*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp]
		(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = MIN(imdata, MAX=maximum)
	ENDIF ELSE IF ((*(*info).scaling).imrefscaling EQ 1) THEN BEGIN
		IF ((*(*info).dataparams).refnt EQ 0) THEN refdata = *(*(*info).data).refdata ELSE IF ((*(*info).dataparams).refnt EQ 1) THEN refdata = (*(*(*info).data).refdata)[0] ELSE $
			IF ((*(*info).dataparams).refnlp EQ 1) THEN refdata = (*(*(*info).data).refdata)[(*(*info).dataparams).t] ELSE $
			refdata = (*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref]
		(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = MIN(refdata, MAX=maximum)
	ENDIF ELSE IF ((*(*info).scaling).imrefscaling EQ 2) THEN BEGIN
		(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling] = MIN(*(*(*info).data).dopslice, MAX=maximum)
	ENDIF
	(*(*(*info).scaling).scale_range)[(*(*info).scaling).imrefscaling] = maximum - (*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).scale_minimum)[(*(*info).scaling).imrefscaling],maximum],labels=['Scaling image','Scaling setting','Scaling minimum','Scaling maximum']
	CRISPEX_SCALING_SET_SLIDERS, event
	WIDGET_CONTROL, (*(*info).ctrlscp).abs_scale_but, /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).rel_scale_but, /SENSITIVE
	CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_REL, event
; Enables the relative or absolute scaling in the range 0-100% or 0-255 respectively
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_MAN_REL'
	(*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] = event.SELECT
	IF (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] THEN BEGIN
		(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] / 2.55
		(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] / 2.55
		WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling], /SENSITIVE, SET_SLIDER_MIN = 1, SET_SLIDER_MAX = 100
		WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling], /SENSITIVE, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = 99
	ENDIF ELSE BEGIN
		(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] * 2.55
		(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] * 2.55
		WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling], /SENSITIVE, SET_SLIDER_MIN = 1, SET_SLIDER_MAX = 255
		WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling], /SENSITIVE, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = 254
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling],(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling]],labels=['Scaling image','Absolute minimum','Absolute maximum','Relative minimum','Relative maximum']
	CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_MAX_SLIDER, event
; Handles events from the maximum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_MAX_SLIDER'
	IF (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] THEN BEGIN
		(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] = event.VALUE
		IF ((*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] LE (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling]) THEN BEGIN
			(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] - 1
			WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling]
		ENDIF
	ENDIF ELSE BEGIN		
		(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] = event.VALUE
		IF ((*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] LE (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling]) THEN BEGIN
			(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] - 1
			WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling]
		ENDIF
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling],(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling]],labels=['Scaling image','Absolute minimum','Absolute maximum','Relative minimum','Relative maximum']
	IF ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 3) THEN CRISPEX_SCALING_MAN_CURR, event ELSE CRISPEX_SCALING_REDRAW, event 
END

PRO CRISPEX_SCALING_MIN_SLIDER, event
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_MIN_SLIDER'
	IF (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] THEN BEGIN
		(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] = event.VALUE
		IF ((*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] GE (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling]) THEN BEGIN
			(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling] + 1
			WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling]
		ENDIF
	ENDIF ELSE BEGIN
		(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] = event.VALUE
		IF ((*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] GE (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling]) THEN BEGIN
			(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling] = (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling] + 1
			WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling]
		ENDIF
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling],(*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling],$
		(*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling]],labels=['Scaling image','Absolute minimum','Absolute maximum','Relative minimum','Relative maximum']
	IF ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 3) THEN CRISPEX_SCALING_MAN_CURR, event ELSE CRISPEX_SCALING_REDRAW, event 
END

PRO CRISPEX_SCALING_XY_SELECT, event
; Handles the selection of scaling the main image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_XY_SELECT'
	(*(*info).scaling).imrefscaling = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [((*(*info).scaling).imrefscaling EQ 0)],labels=['Image scaling select']
	CRISPEX_SCALING_SET_BUTTONS, event
	CRISPEX_SCALING_SET_SLIDERS, event
END

PRO CRISPEX_SCALING_REF_SELECT, event
; Handles the selection of scaling the reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_REF_SELECT'
	(*(*info).scaling).imrefscaling = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [((*(*info).scaling).imrefscaling EQ 1)],labels=['Reference scaling select']
	CRISPEX_SCALING_SET_BUTTONS, event
	CRISPEX_SCALING_SET_SLIDERS, event
END

PRO CRISPEX_SCALING_DOP_SELECT, event
; Handles the selection of scaling the reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_DOP_SELECT'
	(*(*info).scaling).imrefscaling = 2
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [((*(*info).scaling).imrefscaling EQ 2)],labels=['Doppler scaling select']
	CRISPEX_SCALING_SET_BUTTONS, event
	CRISPEX_SCALING_SET_SLIDERS, event
END

PRO CRISPEX_SCALING_REDRAW, event
; Handles the redrawing of window contents after adjustment of scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_REDRAW'
	IF ((*(*info).scaling).imrefscaling EQ 0) THEN BEGIN
		CRISPEX_DRAW_XY, event 
		IF (*(*info).winswitch).showsp THEN CRISPEX_DRAW_SP, event 
	ENDIF ELSE BEGIN
		IF (*(*info).winswitch).showref THEN CRISPEX_DRAW_REF, event
		IF (*(*info).winswitch).showdop THEN CRISPEX_DRAW_DOPPLER, event
	ENDELSE
	IF ((*(*info).dispparams).slices_imscale AND (*(*info).winswitch).showrestloop) THEN CRISPEX_DRAW_REST_LOOP, event
END

PRO CRISPEX_SCALING_SET_BUTTONS, event
; Handles the setting of scaling buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_SET_BUTTONS'
	WIDGET_CONTROL, (*(*info).ctrlscp).auto_const, SET_BUTTON = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 0) 
	WIDGET_CONTROL, (*(*info).ctrlscp).auto_sing, SET_BUTTON = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlscp).man_first, SET_BUTTON = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 2) 
	WIDGET_CONTROL, (*(*info).ctrlscp).man_curr, SET_BUTTON = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 3)
	sens = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] GE 2)
	WIDGET_CONTROL, (*(*info).ctrlscp).abs_scale_but, SENSITIVE = sens, SET_BUTTON = ABS( (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] - 1 )
	WIDGET_CONTROL, (*(*info).ctrlscp).rel_scale_but, SENSITIVE = sens, SET_BUTTON = (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling],$
		ABS((*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling]-1),(*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling]],$
		labels=['Image scaling','Scaling button set','Absolute scaling set','Relative scaling set']
END

PRO CRISPEX_SCALING_SET_SLIDERS, event
; Handles the setting of scaling sliders
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SCALING_SET_SLIDERS'
	sens = ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] GE 2)
	IF (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_max_val)[(*(*info).scaling).imrefscaling], SENSITIVE = sens, SET_SLIDER_MIN = 1, SET_SLIDER_MAX = 100
		WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).rel_scale_min_val)[(*(*info).scaling).imrefscaling], SENSITIVE = sens, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = 99
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).max_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_max_val)[(*(*info).scaling).imrefscaling], SENSITIVE = sens, SET_SLIDER_MIN = 1, SET_SLIDER_MAX = 255
		WIDGET_CONTROL, (*(*info).ctrlscp).min_scale_slider, SET_VALUE = (*(*(*info).scaling).scale_min_val)[(*(*info).scaling).imrefscaling], SENSITIVE = sens, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = 254
	ENDELSE
END

;================================================================================= SESSION SAVE/RESTORE PROCEDURES
PRO CRISPEX_SESSION_SAVE_WINDOW, event
; Gets the filename for the session save routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_SAVE_WINDOW'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_GET_FILENAME, event, 'Save session', 'crispex_session','CRISPEX_SESSION_SAVE_CONTINUE', /SESSION_SAVE
END

PRO CRISPEX_SESSION_SAVE_CONTINUE, event
; Checks the session save filename for validity and overwrite problems
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_SAVE_CONTINUE'
	(*(*info).savparams).savpro = 'SESSION'
	CRISPEX_SAVE_CHECK_FILENAME, event, 'cses', 'CRISPEX_SESSION_SAVE_OVER_CONTINUE'
END

PRO CRISPEX_SESSION_SAVE_OVER_CONTINUE, event
; Handles the overwriting and activates the subsequent saving of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_SAVE_OVER_CONTINUE'
	WIDGET_CONTROL, (*(*info).savparams).filename_text, GET_VALUE = session_filename
	CRISPEX_SESSION_SAVE, event, session_filename
	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SESSION_SAVE, event, sesfilename
; Handles the actual saving of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	WIDGET_CONTROL, /HOURGLASS
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_SAVE'
	ctrlsswitch = *(*info).ctrlsswitch	&	curs = *(*info).curs
	dataparams = *(*info).dataparams	&	dataswitch = *(*info).dataswitch		&	detparams = *(*info).detparams
	dispparams = *(*info).dispparams	&	dispswitch = *(*info).dispswitch		&	intparams = *(*info).intparams
	loopparams = *(*info).loopparams	&	loopswitch = *(*info).loopswitch		&	meas = *(*info).meas
	overlayparams = *(*info).overlayparams	&	overlayswitch = *(*info).overlayswitch		&	paramswitch = *(*info).paramswitch
	pbparams = *(*info).pbparams		&	phiparams = *(*info).phiparams
	plotaxes = *(*info).plotaxes		&	plotparams = *(*info).plotparams
	plotpos = *(*info).plotpos		&	plotswitch = *(*info).plotswitch		&	plottitles = *(*info).plottitles
	restoreparams = *(*info).restoreparams	&	retrparams = *(*info).retrparams		&	savswitch = *(*info).savswitch
	scaling = *(*info).scaling		&	stokesparams = *(*info).stokesparams		&	versioninfo = *(*info).versioninfo
	winsizes = *(*info).winsizes		&	winswitch = *(*info).winswitch			&	zooming = *(*info).zooming
	SAVE, ctrlsswitch, curs, dataparams, dataswitch, detparams, dispparams, dispswitch, intparams, loopparams, loopswitch, meas, overlayparams, overlayswitch, paramswitch, pbparams, $
		phiparams, plotaxes, plotparams, plotpos, plotswitch, plottitles, restoreparams, retrparams, savswitch, scaling, stokesparams, versioninfo, winsizes, winswitch, zooming, $
		FILENAME = (*(*info).paths).opath+sesfilename+'.cses'
	PRINT,'Written: '+(*(*info).paths).opath+sesfilename+'.cses'
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SESSION_RESTORE_WINDOW, event
; Opens the session restore window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_RESTORE_WINDOW'
	*(*(*info).sesparams).csesfiles = FILE_SEARCH((*(*info).paths).ipath+"*cses", COUNT = csesfilecount)
	IF (csesfilecount GT 0) THEN BEGIN
		eventval = INDGEN(csesfilecount)
		base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Restore session', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
		disp = WIDGET_BASE(base, /COLUMN)
		buts_disp = WIDGET_BASE(disp, /COLUMN, /FRAME)
		choose_text = WIDGET_LABEL(buts_disp, VALUE = 'Choose the session file that is to be restored:', /ALIGN_LEFT)
		IF (csesfilecount GT 15) THEN sel_buts = WIDGET_BASE(buts_disp, /COLUMN, /EXCLUSIVE, Y_SCROLL_SIZE = 425) ELSE sel_buts = WIDGET_BASE(buts_disp, /COLUMN, /EXCLUSIVE)
		FOR i=0,csesfilecount-1 DO BEGIN
			(*(*(*info).sesparams).sessions)[i] = 0
			singlefilename = (STRSPLIT((*(*(*info).sesparams).csesfiles)[i], PATH_SEP(), /EXTRACT))[N_ELEMENTS(STRSPLIT((*(*(*info).sesparams).csesfiles)[i], PATH_SEP(), /EXTRACT))-1]
			but_val = singlefilename
			sel_but = WIDGET_BUTTON(sel_buts, VALUE = but_val, UVALUE = eventval[i], EVENT_PRO = 'CRISPEX_SESSION_RESTORE_EVENT') 
		ENDFOR
		dec_buts = WIDGET_BASE(disp, COLUMN=2, /GRID_LAYOUT, /ALIGN_CENTER)
		cancel = WIDGET_BUTTON(dec_buts, VALUE = 'Cancel', EVENT_PRO = 'CRISPEX_CLOSE_EVENT_WINDOW')
		(*(*info).sesparams).rest_sessions = WIDGET_BUTTON(dec_buts, VALUE = 'Restore session', EVENT_PRO = 'CRISPEX_SESSION_RESTORE', SENSITIVE = 0)
		WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
		WIDGET_CONTROL, base, SET_UVALUE = info
		XMANAGER, 'CRISPEX', base, /NO_BLOCK
		(*(*info).winids).restsestlb = base
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [csesfilecount,(*(*info).winids).restsestlb], labels=['Restored sessions','restsestlb']
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','No stored session (*.cses) files found.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END

PRO CRISPEX_SESSION_RESTORE_EVENT, event
; Handles the selection of a session to be restored
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_RESTORE_EVENT'
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).sesparams).sessions)[eventval] = ( (*(*(*info).sesparams).sessions)[eventval] EQ 0) 
	condition = WHERE(*(*(*info).sesparams).sessions EQ 1)
	WIDGET_CONTROL, (*(*info).sesparams).rest_sessions, SENSITIVE = ((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).retrparams).sel_loops)[eventval]], labels=['Restored session ID','Restored session selected']
END

PRO CRISPEX_SESSION_RESTORE_READ_POINTER, event, currpointer, restpointer, NO_RESTORE=no_restore
; Handles the actual restoration of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_RESTORE_READ_POINTER'
	currtags = TAG_NAMES(currpointer)	&	resttags = TAG_NAMES(restpointer)
	ncurr = N_ELEMENTS(currtags)		&	nrest = N_ELEMENTS(resttags)
	no_rest = N_ELEMENTS(NO_RESTORE)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [ncurr,nrest,no_rest-1],labels=['Current tags','Restored tags','Prevent replace tag']
	IF (ncurr EQ nrest) THEN BEGIN											; Both the restored and current pointer have the same number of tags
		IF (no_rest EQ 0) THEN BEGIN										; Pointer to be restored without skipping tags
			IF (WHERE((resttags EQ currtags) EQ 0) EQ -1) THEN currpointer = restpointer ELSE BEGIN		; If all tagnames are equal, then just replace the pointer
				FOR i=0,nrest-1 DO (currpointer).(WHERE(currtags EQ resttags[i])) = restpointer.(i)	; Else go through all tags and replace only where tagnames are the same
			ENDELSE
		ENDIF ELSE BEGIN											; Pointer to be restored while skipping tags
			no_replace = WHERE((STRLOWCASE(resttags) EQ STRLOWCASE(no_restore)) EQ 1)
			FOR i=0,nrest-1 DO BEGIN
				IF (no_replace NE i) THEN (currpointer).(i) = restpointer.(i)				; If pointer tag is not the one to be skipped, replace its value
			ENDFOR
		ENDELSE
	ENDIF ELSE BEGIN												; Unequal number of tags between current and restored pointer
		FOR i=0,nrest-1 DO BEGIN
			IF (WHERE(currtags EQ resttags[i]) NE -1) THEN (currpointer).(WHERE(currtags EQ resttags[i])) = restpointer.(i)			; Go through all tags and replace only where tagnames are the same
		ENDFOR
	ENDELSE
END

PRO CRISPEX_SESSION_RESTORE, event
; Handles the actual restoration of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SESSION_RESTORE'
	WIDGET_CONTROL, /HOURGLASS
	CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=0, feedback_text='Restoring session file and checking version...', /SESSION
	restore_session = (*(*(*info).sesparams).csesfiles)[WHERE(*(*(*info).sesparams).sessions EQ 1)]
	RESTORE, restore_session
	; Check revision number
	IF (N_ELEMENTS(versioninfo) GT 0) THEN cont = (versioninfo.revision_number GE '546') ELSE cont = 0
	IF cont THEN BEGIN
		IF (((*(*info).dataparams).imfilename EQ dataparams.imfilename) AND ((*(*info).dataparams).spfilename EQ dataparams.spfilename) AND $
			((*(*info).dataparams).reffilename EQ dataparams.reffilename) AND ((*(*info).dataparams).refspfilename EQ dataparams.refspfilename)) THEN BEGIN
			; Kill all open windows (except control panel and main image window)
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=1, feedback_text='Closing open windows...'
			IF ((*(*info).winids).sptlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).sptlb, /DESTROY 
				(*(*info).winids).sptlb = 0 
			ENDIF
			IF ((*(*info).winids).lstlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).lstlb, /DESTROY
				(*(*info).winids).lstlb = 0 
			ENDIF
			IF ((*(*info).winids).phistlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).phistlb, /DESTROY 
				(*(*info).winids).phistlb = 0 
			ENDIF
			IF ((*(*info).winids).reftlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).reftlb, /DESTROY 
				(*(*info).winids).reftlb = 0 
			ENDIF
			IF ((*(*info).winids).doptlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).doptlb, /DESTROY 
				(*(*info).winids).doptlb = 0 
			ENDIF
			IF ((*(*info).winids).imreftlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).imreftlb, /DESTROY 
				(*(*info).winids).imreftlb = 0 
			ENDIF
			IF (TOTAL(*(*(*info).winids).restlooptlb) NE 0) THEN BEGIN 
				FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO WIDGET_CONTROL, (*(*(*info).winids).restlooptlb)[i], /DESTROY
				(*(*info).winids).restlooptlb = PTR_NEW(0) 
			ENDIF
			IF ((*(*info).winids).retrdettlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).retrdettlb, /DESTROY 
				(*(*info).winids).retrdettlb = 0 
			ENDIF
			IF ((*(*info).winids).looptlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).looptlb, /DESTROY 
				(*(*info).winids).looptlb = 0 
			ENDIF
			IF ((*(*info).winids).reflooptlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).reflooptlb, /DESTROY 
				(*(*info).winids).reflooptlb = 0 
			ENDIF
			IF ((*(*info).winids).refsptlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).refsptlb, /DESTROY 
				(*(*info).winids).refsptlb = 0 
			ENDIF
			IF ((*(*info).winids).reflstlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).reflstlb, /DESTROY 
				(*(*info).winids).reflstlb = 0 
			ENDIF
			IF ((*(*info).winids).inttlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).inttlb, /DESTROY 
				(*(*info).winids).inttlb = 0 
			ENDIF
			IF ((*(*info).winids).paramtlb NE 0) THEN BEGIN 
				WIDGET_CONTROL, (*(*info).winids).paramtlb, /DESTROY 
				(*(*info).winids).paramtlb = 0 
			ENDIF
			; Restore all parameters, for each parameter checking whether the number of subparameters has changed (if revision_number on session file differs from current CRISPEX)
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=1, feedback_text='Restoring session variables...'
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).ctrlsswitch, ctrlsswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).curs, curs
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).dataparams, dataparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).dataswitch, dataswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).detparams, detparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).dispparams, dispparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).dispswitch, dispswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).intparams, intparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).loopparams, loopparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).loopswitch, loopswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).meas, meas
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).overlayparams, overlayparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).overlayswitch, overlayswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).paramswitch, paramswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).pbparams, pbparams, NO_RESTORE='BG'
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).phiparams, phiparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).plotaxes, plotaxes
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).plotparams, plotparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).plotpos, plotpos
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).plotswitch, plotswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).plottitles, plottitles
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).restoreparams, restoreparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).retrparams, retrparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).scaling, scaling
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).savswitch, savswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).stokesparams, stokesparams
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).winsizes, winsizes
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).winswitch, winswitch
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).zooming, zooming
			; Reset controls on control panel given the control switches and other parameters
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=1, feedback_text='Resetting control panel controls...'
			; Temporal
			WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dataparams).t, SET_SLIDER_MIN = (*(*info).dispparams).t_low, SET_SLIDER_MAX = (*(*info).dispparams).t_upp, $
				SENSITIVE = ((*(*info).dataparams).nt GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = ((*(*info).dataparams).nt GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = ((*(*info).dataparams).nt GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = ((*(*info).dispparams).t_range NE (*(*info).dataparams).nt)
			WIDGET_CONTROL, (*(*info).ctrlscp).t_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed, SENSITIVE = ((*(*info).dataparams).nt GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).t_step_slider, SET_SLIDER_MAX = (*(*info).dispparams).t_range - 1, SET_VALUE = (*(*info).pbparams).t_step, $
				SENSITIVE = (((*(*info).dispparams).t_range - 1 NE 1) AND (*(*info).dataparams).nt GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).imref_blink_but, SET_BUTTON = (*(*info).winswitch).showimref, SENSITIVE = (*(*info).dataswitch).reffile
			; Spectral
			WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2), SENSITIVE = ((*(*info).dataparams).nlp GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2), SENSITIVE = ((*(*info).dataparams).nlp GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, SENSITIVE = ((*(*info).dispparams).lp_range NE (*(*info).dataparams).nlp)
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp, SET_SLIDER_MIN = (*(*info).dispparams).lp_low, SET_SLIDER_MAX = (*(*info).dispparams).lp_upp, $
				SENSITIVE = ((*(*info).dataparams).nlp GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed, SENSITIVE = ((*(*info).dataswitch).spfile AND (*(*info).pbparams).spmode)
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_slider, SET_SLIDER_MAX = (*(*info).dispparams).lp_range - 1, SET_VALUE = (*(*info).pbparams).lp_step, $
				SENSITIVE = (((*(*info).dispparams).lp_range - 1 NE 1) AND (*(*info).dataparams).nlp GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, SET_BUTTON = (*(*info).pbparams).spmode, SENSITIVE = (((*(*info).dataparams).nlp GT 1) AND ((*(*info).winswitch).showimref EQ 0))
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SET_BUTTON = ((*(*info).ctrlsswitch).lp_ref_lock AND ((*(*info).dataparams).refnlp GT 1)), $
				SENSITIVE = (((*(*info).dataparams).nlp EQ (*(*info).dataparams).refnlp) AND ((*(*info).dataparams).refnlp GT 1))
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SET_VALUE = (*(*info).dataparams).lp_ref, SENSITIVE = ((*(*info).dataswitch).refspfile AND ABS((*(*info).ctrlsswitch).lp_ref_lock-1))
			CRISPEX_PB_BUTTONS_SET, event, bwd_set=((*(*info).pbparams).direction EQ -1), pause_set=((*(*info).pbparams).mode EQ 'PAUSE'), fwd_set=((*(*info).pbparams).direction EQ 1), $
				loop_set=((*(*info).pbparams).lmode EQ 'LOOP'),	cycle_set=((*(*info).pbparams).lmode EQ 'CYCLE'), blink_set=((*(*info).pbparams).lmode EQ 'BLINK')
			; Spatial
			CRISPEX_COORDSLIDERS_SET, ABS((*(*info).curs).lockset-1), ABS((*(*info).curs).lockset-1), event
			WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = (*(*info).curs).lockset
			WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, SET_BUTTON = ABS((*(*info).curs).lockset-1)
			IF ((*(*info).zooming).factor NE 1) THEN BEGIN
				CRISPEX_ZOOM, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, /SENSITIVE
				WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, /SENSITIVE
			ENDIF
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_one, SET_BUTTON = ((*(*info).zooming).factor EQ 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_two, SET_BUTTON = ((*(*info).zooming).factor EQ 2)
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_three, SET_BUTTON = ((*(*info).zooming).factor EQ 3)
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_four, SET_BUTTON = ((*(*info).zooming).factor EQ 4)
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_six, SET_BUTTON = ((*(*info).zooming).factor EQ 6)
			WIDGET_CONTROL, (*(*info).ctrlscp).zoom_eight, SET_BUTTON = ((*(*info).zooming).factor EQ 8)
			WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SET_VALUE = (*(*info).zooming).xpos
			WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SET_VALUE = (*(*info).zooming).ypos
			; Stokes
			stokes_i_available = (WHERE(((*(*info).stokesparams).labels) EQ 'I') GE 0)
			stokes_q_available = (WHERE(((*(*info).stokesparams).labels) EQ 'Q') GE 0)
			stokes_u_available = (WHERE(((*(*info).stokesparams).labels) EQ 'U') GE 0)
			stokes_v_available = (WHERE(((*(*info).stokesparams).labels) EQ 'V') GE 0)
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_xy_i_but, SET_BUTTON = (((*(*info).stokesparams).labels)[(*(*info).dataparams).s] EQ 'I'), SENSITIVE = stokes_i_available
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_xy_q_but, SET_BUTTON = (((*(*info).stokesparams).labels)[(*(*info).dataparams).s] EQ 'Q'), SENSITIVE = stokes_q_available
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_xy_u_but, SET_BUTTON = (((*(*info).stokesparams).labels)[(*(*info).dataparams).s] EQ 'U'), SENSITIVE = stokes_u_available
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_xy_v_but, SET_BUTTON = (((*(*info).stokesparams).labels)[(*(*info).dataparams).s] EQ 'V'), SENSITIVE = stokes_v_available
			pol_sp_i_set = (WHERE(((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))] EQ 'I') GE 0)
			pol_sp_q_set = (WHERE(((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))] EQ 'Q') GE 0)
			pol_sp_u_set = (WHERE(((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))] EQ 'U') GE 0)
			pol_sp_v_set = (WHERE(((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))] EQ 'V') GE 0)
			pol_sp_i_sens = (((TOTAL((*(*info).stokesparams).select_sp) GT 1) AND stokes_i_available) OR ((TOTAL((*(*info).stokesparams).select_sp) EQ 1) AND ABS(pol_sp_i_set-1) AND stokes_i_available))
			pol_sp_q_sens = (((TOTAL((*(*info).stokesparams).select_sp) GT 1) AND stokes_q_available) OR ((TOTAL((*(*info).stokesparams).select_sp) EQ 1) AND ABS(pol_sp_q_set-1) AND stokes_q_available))
			pol_sp_u_sens = (((TOTAL((*(*info).stokesparams).select_sp) GT 1) AND stokes_u_available) OR ((TOTAL((*(*info).stokesparams).select_sp) EQ 1) AND ABS(pol_sp_u_set-1) AND stokes_u_available))
			pol_sp_v_sens = (((TOTAL((*(*info).stokesparams).select_sp) GT 1) AND stokes_v_available) OR ((TOTAL((*(*info).stokesparams).select_sp) EQ 1) AND ABS(pol_sp_v_set-1) AND stokes_v_available))
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_i_but, SENSITIVE = pol_sp_i_sens, SET_BUTTON = pol_sp_i_set
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_q_but, SENSITIVE = pol_sp_q_sens, SET_BUTTON = pol_sp_q_set
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_u_but, SENSITIVE = pol_sp_u_sens, SET_BUTTON = pol_sp_u_set
			WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_v_but, SENSITIVE = pol_sp_v_sens, SET_BUTTON = pol_sp_v_set
			IF ((*(*info).dataparams).nlp EQ 1) THEN WIDGET_CONTROL, (*(*info).ctrlscp).pol_sp_i_but, SET_BUTTON = 0
			; Displays
			WIDGET_CONTROL, (*(*info).ctrlscp).detspect_im_but, SET_BUTTON = ABS((*(*info).ctrlsswitch).imrefdetspect-1)
			WIDGET_CONTROL, (*(*info).ctrlscp).detspect_ref_but, SET_BUTTON = (*(*info).ctrlsswitch).imrefdetspect, SENSITIVE = (*(*info).dataswitch).refspfile
			CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
			WIDGET_CONTROL, (*(*info).ctrlscp).sp_toggle_but, SET_BUTTON = (*(*info).winswitch).showsp
			WIDGET_CONTROL, (*(*info).ctrlscp).phis_toggle_but, SET_BUTTON = (*(*info).winswitch).showphis
			WIDGET_CONTROL, (*(*info).ctrlscp).refsp_toggle_but, SET_BUTTON = (*(*info).winswitch).showrefsp, SENSITIVE = (*(*info).dataswitch).refspfile
			WIDGET_CONTROL, (*(*info).ctrlscp).int_toggle_but, SET_BUTTON = (*(*info).winswitch).showint, SENSITIVE = (*(*info).dataswitch).spfile
			WIDGET_CONTROL, (*(*info).ctrlscp).reference_but, SET_BUTTON = (*(*info).winswitch).showref, SENSITIVE = (*(*info).dataswitch).reffile
			WIDGET_CONTROL, (*(*info).ctrlscp).doppler_but, SET_BUTTON = (*(*info).winswitch).showdop, SENSITIVE = ((*(*info).dataparams).nlp GT 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).param_but, SET_BUTTON = (*(*info).winswitch).showparam
			; Scaling
			WIDGET_CONTROL, (*(*info).ctrlscp).xy_scaling_but, SET_BUTTON = ((*(*info).scaling).imrefscaling EQ 0)
			WIDGET_CONTROL, (*(*info).ctrlscp).ref_scaling_but, SET_BUTTON = ((*(*info).scaling).imrefscaling EQ 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).dop_scaling_but, SET_BUTTON = ((*(*info).scaling).imrefscaling EQ 2)
			CRISPEX_SCALING_SET_BUTTONS, event
			CRISPEX_SCALING_SET_SLIDERS, event
			; Slits
			WIDGET_CONTROL, (*(*info).ctrlscp).phi_slider, SET_VALUE = (*(*info).phiparams).angle, SENSITIVE = (*(*info).winswitch).showphis
			WIDGET_CONTROL, (*(*info).ctrlscp).nphi_slider, SET_VALUE = (*(*info).phiparams).nphi_set, SENSITIVE = (*(*info).winswitch).showphis
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_slit_but, SET_BUTTON = ((*(*info).loopparams).np GE 2), SENSITIVE = ABS((*(*info).meas).spatial_measurement-1)
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_feedb_but, SET_BUTTON = (*(*info).overlayswitch).looppath_feedback
			WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = ((*(*info).loopparams).np GE 3)
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = (ABS((*(*info).winswitch).showloop-1) AND ABS((*(*info).meas).spatial_measurement-1) AND ((*(*info).loopparams).np GE 2))
			; Miscellaneous
			WIDGET_CONTROL, (*(*info).ctrlscp).overlay_but, SET_BUTTON = (*(*info).loopswitch).restore_loops
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_all, SET_BUTTON = (*(*info).overlayswitch).overlalways, SENSITIVE = (*(*info).loopswitch).restore_loops
			WIDGET_CONTROL, (*(*info).ctrlscp).loop_overlay_sav, SET_BUTTON = ABS((*(*info).overlayswitch).overlalways-1), SENSITIVE = (*(*info).loopswitch).restore_loops
			WIDGET_CONTROL, (*(*info).ctrlscp).linestyle_0, SET_BUTTON = ((*(*info).overlayparams).loop_linestyle EQ 0)
			WIDGET_CONTROL, (*(*info).ctrlscp).linestyle_1, SET_BUTTON = ((*(*info).overlayparams).loop_linestyle EQ 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).linestyle_2, SET_BUTTON = ((*(*info).overlayparams).loop_linestyle EQ 2)
			WIDGET_CONTROL, (*(*info).ctrlscp).measure_but, SET_BUTTON = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).apix_label, SENSITIVE = (*(*info).meas).spatial_measurement
	    WIDGET_CONTROL, (*(*info).ctrlscp).apix_unit, SENSITIVE = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).apix_text, SET_VALUE = STRTRIM((*(*info).meas).arcsecpix,2), SENSITIVE = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_lab, SENSITIVE = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SENSITIVE = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_lab, SENSITIVE = (*(*info).meas).spatial_measurement
			WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SENSITIVE = (*(*info).meas).spatial_measurement
			; Mask
			WIDGET_CONTROL, (*(*info).ctrlscp).masks_overlay_ct_cbox, SET_COMBOBOX_SELECT = (*(*info).overlayparams).maskct, SENSITIVE = (*(*info).overlayswitch).mask
			WIDGET_CONTROL, (*(*info).ctrlscp).masks_overlay_col_slid, SET_VALUE = (*(*info).overlayparams).maskcolor, SENSITIVE = (*(*info).overlayswitch).mask
			CRISPEX_MASK_BUTTONS_SET, event
			IF (*(*info).meas).spatial_measurement THEN CRISPEX_MEASURE_CALC, event
			; Open windows for replotting and replot
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=1, feedback_text='Opening windows and refreshing displays...'
			IF (*(*info).winswitch).showparam THEN BEGIN
				(*(*info).winswitch).showparam = 0
				CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showref THEN BEGIN
				(*(*info).winswitch).showref = 0
				CRISPEX_DISPLAYS_REF_TOGGLE, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showdop THEN BEGIN
				(*(*info).winswitch).showdop = 0
				CRISPEX_DISPLAYS_DOPPLER_TOGGLE, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showimref THEN BEGIN
				(*(*info).winswitch).showimref = 0
				CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showsp THEN BEGIN
				(*(*info).winswitch).showsp = 0
				CRISPEX_DISPLAYS_SP_TOGGLE, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showrefsp THEN BEGIN
				(*(*info).winswitch).showrefsp = 0
				CRISPEX_DISPLAYS_REFSP_TOGGLE, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			showimref_ls = [(*(*info).winswitch).showls, (*(*info).winswitch).showrefls]
			FOR i = 0,1 DO BEGIN
				IF showimref_ls[i] THEN BEGIN
					set_imrefdetspect = (*(*info).ctrlsswitch).imrefdetspect
					(*(*info).ctrlsswitch).imrefdetspect = i
					IF (i EQ 0) THEN (*(*info).winswitch).showls = 0 ELSE (*(*info).winswitch).showrefls = 0
					CRISPEX_DISPLAYS_IMREF_LS_TOGGLE, event, /NO_DRAW
					WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
					(*(*info).ctrlsswitch).imrefdetspect = set_imrefdetspect
				ENDIF
			ENDFOR
			IF (*(*info).winswitch).showint THEN BEGIN
				(*(*info).winswitch).showint = 0
				CRISPEX_DISPLAYS_INT_TOGGLE, event, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			old_cfilecount = (*(*info).restoreparams).cfilecount
			CRISPEX_FIND_CSAV, event
			IF ((*(*info).loopswitch).restore_loops AND ((*(*info).restoreparams).cfilecount EQ old_cfilecount)) THEN BEGIN
				CRISPEX_RESTORE_LOOPS_MENU, event, *(*(*info).restoreparams).sel_loops
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
				WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_all, SET_BUTTON = (TOTAL(*(*(*info).restoreparams).sel_loops) EQ N_ELEMENTS(*(*(*info).restoreparams).sel_loops))
				WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_none, SET_BUTTON = (TOTAL(*(*(*info).restoreparams).sel_loops) EQ 0)
				WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
				IF (*(*info).winswitch).showrestloop THEN BEGIN
					FOR i=0,N_ELEMENTS(*(*(*info).restoreparams).disp_loopnr)-1 DO BEGIN
						list_values[(*(*(*info).restoreparams).disp_loopnr)[i]+1] = 'Hide time slice '+STRTRIM((*(*(*info).restoreparams).disp_loopnr)[i],2)
						(*(*info).restoreparams).disp_loopfile = (*(*(*info).restoreparams).cfiles)[(*(*(*info).restoreparams).disp_loopnr)[i]]
						CRISPEX_DISPLAYS_RESTORE_LOOPSLAB, event, /NO_DRAW, INDEX=i
					ENDFOR
				ENDIF ELSE FOR i=1,N_ELEMENTS(list_values)-1 DO list_values[i] = 'Display time slice '+STRTRIM(i-1,2)
				WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, SET_COMBOBOX_SELECT = (*(*(*info).restoreparams).disp_loopnr)[0]+1
			ENDIF ELSE BEGIN			; Add error message
				(*(*info).loopswitch).restore_loops = 0	
				(*(*info).winswitch).showrestloop = 0
			ENDELSE
			old_clfilecount = (*(*info).retrparams).clfilecount
			CRISPEX_FIND_CLSAV, event
			IF ((*(*info).loopswitch).retrieve_loops AND ((*(*info).retrparams).clfilecount EQ old_clfilecount)) THEN BEGIN
				CRISPEX_RETRIEVE_LOOP_MENU, event, *(*(*info).retrparams).sel_loops
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
				WIDGET_CONTROL, (*(*info).ctrlsloop).sel_all, SET_BUTTON = (TOTAL(*(*(*info).retrparams).sel_loops) EQ N_ELEMENTS(*(*(*info).retrparams).sel_loops))
				WIDGET_CONTROL, (*(*info).ctrlsloop).sel_none, SET_BUTTON = (TOTAL(*(*(*info).retrparams).sel_loops) EQ 0)
				WIDGET_CONTROL, (*(*info).ctrlsloop).all_pos, SET_BUTTON = (*(*info).savswitch).all_pos_loops
				WIDGET_CONTROL, (*(*info).ctrlsloop).saved_pos, SET_BUTTON = ABS((*(*info).savswitch).all_pos_loops-1)
				WIDGET_CONTROL, (*(*info).ctrlsloop).del_files, SET_BUTTON = (*(*info).savswitch).delete_clsav
				WIDGET_CONTROL, (*(*info).ctrlsloop).keep_files, SET_BUTTON = ABS((*(*info).savswitch).delete_clsav-1)
				WIDGET_CONTROL, (*(*info).ctrlsloop).save_imonly, SET_BUTTON = ((*(*info).savswitch).imref_only EQ 1)
				WIDGET_CONTROL, (*(*info).ctrlsloop).save_refonly, SET_BUTTON = ((*(*info).savswitch).imref_only EQ 2)
				WIDGET_CONTROL, (*(*info).ctrlsloop).save_imref, SET_BUTTON = ((*(*info).savswitch).imref_only EQ 3)
			ENDIF ELSE (*(*info).loopswitch).retrieve_loops = 0	; Add error message
			detfile = FILE_SEARCH((*(*info).detparams).detfilename, COUNT=detcount)
			IF ((*(*info).loopswitch).retrieve_detfile AND (detcount EQ 1)) THEN BEGIN
				CRISPEX_RETRIEVE_DET_FILE_MENU, event, *(*(*info).detparams).sel_dets, DETFILENAME=(*(*info).detparams).detfilename, /NO_DRAW
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
				WIDGET_CONTROL, (*(*info).ctrlsdet).sel_all, SET_BUTTON = (TOTAL(*(*(*info).detparams).sel_dets) EQ N_ELEMENTS(*(*(*info).detparams).sel_dets))
				WIDGET_CONTROL, (*(*info).ctrlsdet).sel_none, SET_BUTTON = (TOTAL(*(*(*info).detparams).sel_dets) EQ 0)
				WIDGET_CONTROL, (*(*info).ctrlsdet).disp_list, SET_COMBOBOX_SELECT = (*(*info).detparams).idx+1
				WIDGET_CONTROL, (*(*info).ctrlsdet).overlay_all, SET_BUTTON = (*(*info).overlayswitch).det_overlay_all
				WIDGET_CONTROL, (*(*info).ctrlsdet).overlay_sel, SET_BUTTON = ABS((*(*info).overlayswitch).det_overlay_all-1)
				WIDGET_CONTROL, (*(*info).ctrlsdet).width_slider, SET_VALUE = (*(*info).detparams).width
				WIDGET_CONTROL, (*(*info).ctrlsdet).all_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 1)
				WIDGET_CONTROL, (*(*info).ctrlsdet).saved_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 2)
				WIDGET_CONTROL, (*(*info).ctrlsdet).sel_range_pos, SET_BUTTON = ((*(*info).savswitch).pos_dets EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).dlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmin_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).refdlpmax_text, SENSITIVE = ((*(*info).savswitch).pos_dets EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).save_imonly, SET_BUTTON = ((*(*info).savswitch).det_imref_only EQ 1)
				WIDGET_CONTROL, (*(*info).ctrlsdet).save_refonly, SET_BUTTON = ((*(*info).savswitch).det_imref_only EQ 2)
				WIDGET_CONTROL, (*(*info).ctrlsdet).save_imref, SET_BUTTON = ((*(*info).savswitch).det_imref_only EQ 3)
				WIDGET_CONTROL, (*(*info).ctrlsdet).get_dets, SENSITIVE = (TOTAL(*(*(*info).detparams).sel_dets) GE 1)
				IF (*(*info).winswitch).showretrdet THEN CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB, event, /NO_DRAW
				CRISPEX_DISPRANGE_LP_RANGE, event
				CRISPEX_UPDATE_LP, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF ELSE (*(*info).loopswitch).retrieve_detfile = 0	; Add error message
			IF (*(*info).winswitch).showloop THEN BEGIN
				CRISPEX_DISPLAYS_LOOPSLAB, event, /NO_DRAW
				IF (*(*info).winswitch).showrefloop THEN CRISPEX_DISPLAYS_REFLOOPSLAB, event, /NO_DRAW
				CRISPEX_UPDATE_LP, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF
			IF (*(*info).winswitch).showphis THEN BEGIN
				(*(*info).winswitch).showphis = 0
				IF ((*(*info).dataparams).nt EQ 1) THEN CRISPEX_UPDATE_T, event
				CRISPEX_DISPLAYS_PHIS_TOGGLE, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDIF ELSE BEGIN
				CRISPEX_UPDATE_T, event
				CRISPEX_DRAW, event
				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
			ENDELSE
			; Menu
			IF (*(*info).winswitch).dispwids THEN BEGIN
				(*(*info).winswitch).dispwids = 0
				CRISPEX_DISPWIDS, event
			ENDIF
			CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event, /SESSION
			(*(*(*info).sesparams).sessions)[WHERE(*(*(*info).sesparams).sessions EQ 1)] = 0
			WIDGET_CONTROL, (*(*info).winids).restsestlb, /DESTROY
			(*(*info).winids).restsestlb = 0
		ENDIF ELSE BEGIN
			CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event, /SESSION
			CRISPEX_WINDOW_OK, event,'ERROR!','Unable to restore earlier session due to incompatibility','between currently and earlier loaded data.',$
				OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
			(*(*info).winids).errtlb = tlb
		ENDELSE
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event, /SESSION
		IF (N_ELEMENTS(versioninfo) GT 0) THEN message4 = 'Session was saved with CRISPEX v'+versioninfo.version_number+' (rev '+versioninfo.revision_number+').' ELSE message4 = ''
		CRISPEX_WINDOW_OK, event,'ERROR!','Unable to restore earlier session due to incompatibility between','saved and expected session save file format. Running version of',$
			'CRISPEX requires a session saved with CRISPEX v1.6 (rev 542) or later.', message4, OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END

;================================================================================= SAVE LOOPSLICE/SLAB PROCEDURES
PRO CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=infilename, outfilename=outfilename, tlab=tlab, ext=ext, exch_ext=exch_ext, export_id=export_id, import_id=import_id
; Handles the creation of an output filename
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_DETERMINE_FILENAME'
	IF (N_ELEMENTS(INFILENAME) NE 1) THEN infilename = (*(*info).dataparams).imfilename
	infilename = (STRSPLIT(infilename,PATH_SEP(),/EXTRACT))[N_ELEMENTS(STRSPLIT(infilename,PATH_SEP(),/EXTRACT))-1]
	basefstr = STRMID(infilename,0,STRPOS(infilename,'.',/REVERSE_SEARCH))
	IF KEYWORD_SET(EXCH_EXT) THEN outfilename = basefstr ELSE BEGIN
		IF KEYWORD_SET(TLAB) THEN BEGIN
			CRISPEX_SAVE_DETERMINE_SAVEID, event, saveid
			outfilename = basefstr+'_'+saveid
			IF (N_ELEMENTS(EXPORT_ID) EQ 1) THEN export_id = saveid ELSE export_id=''
		ENDIF ELSE IF (N_ELEMENTS(IMPORT_ID) EQ 1) THEN outfilename = basefstr+'_'+import_id ELSE $
			outfilename = basefstr
	ENDELSE
	IF (N_ELEMENTS(EXT) EQ 1) THEN outfilename = outfilename+'.'+ext
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [infilename,outfilename],labels=['base filename', 'output filename']
END

PRO CRISPEX_SAVE_DETERMINE_SAVEID, event, saveid, PREF=pref
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_DETERMINE_SAVEID'
	monthstrarr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	curt = SYSTIME()
	tstr = STRSPLIT(curt,/EXTRACT)
	sstr = STRSPLIT(tstr[3],':',/EXTRACT)
	tstr[2] = STRING(tstr[2],FORMAT='(I02)')								; date
	IF KEYWORD_SET(PREF) THEN defsaveid = (*(*info).prefs).tmp_defsaveid ELSE defsaveid = (*(*info).prefs).defsaveid
	IF (defsaveid GE 2) THEN tstr[1] = STRING(WHERE(monthstrarr EQ tstr[1])+1,FORMAT='(I02)')		; month (numeric)
	IF (defsaveid EQ 0) THEN date_order = [4,1,2] ELSE $							; defsaveid = 0 > YYYYMMMDD_hhmmss, ex: 2011Apr12_114545
		IF (defsaveid EQ 1) THEN date_order = [2,1,4] ELSE $						; defsaveid = 1 > DDMMMYYYY_hhmmss, ex: 12Apr2011_114545
		IF (defsaveid EQ 2) THEN date_order = [4,1,2] ELSE $						; defsaveid = 2 > YYYYMMDD_hhmmss, ex: 20110412_114545
		IF (defsaveid EQ 3) THEN date_order = [2,1,4]							; defsaveid = 3 > DDMMYYYY_hhmmss, ex: 12042011_114545
	saveid = STRJOIN(tstr[date_order])+'_'+STRJOIN(sstr)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [saveid],labels=['save ID']
END

PRO CRISPEX_SAVE_LOOP_PTS, event
; Handles the saving of loop points of a defined path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_LOOP_PTS'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
		WIDGET_CONTROl,/HOURGLASS
		PRINT,'Saving current loop points..'
		x_coords = *(*(*info).loopparams).xp		&	y_coords = *(*(*info).loopparams).yp
		x_loop_pts = *(*(*info).loopparams).xr		&	y_loop_pts = *(*(*info).loopparams).yr
		w_loop_pts = *(*(*info).loopparams).w_lpts	&	spect_pos = (*(*info).dataparams).lp
		t_saved = (*(*info).dataparams).t
		crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
		CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=filename, /tlab, ext='clsav'
		SAVE, crispex_version, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, spect_pos, t_saved, FILENAME = (*(*info).paths).opath+filename
		PRINT,'Saving current loop points done. Saved data to: '+STRTRIM(filename,2)
		WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 1
		WIDGET_CONTROL, (*(*info).ctrlscp).all_saved_loop, SENSITIVE = 1
	ENDIF
END

PRO CRISPEX_SAVE_APPROX_LOOPSLICE, event
; Handles the saving of an approximate (i.e. nearest neighbour interpolated) timeslice along a loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	CRISPEX_SAVE_APPROX_LOOPSLAB, event, /SAVE_SLICE
END

PRO CRISPEX_SAVE_APPROX_LOOPSLAB, event, SAVE_SLICE=save_slice
; Handles the saving of an approximate (i.e. nearest neighbour interpolated) timeslab along a loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_APPROX_LOOPSLAB'
	save_message1 = REPLICATE('Saving current approximate loop ',2)+['slab...','slice...']
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=filename, /tlab, ext='csav', export_id=expid
	CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=filename, outfilename=filename2, ext='clsav', /exch_ext
	save_message2 = REPLICATE('Saving current approximate loop ',2)+['slab','slice']+REPLICATE(' done. Saved data to: '+STRTRIM(filename,2),2)
	IF KEYWORD_SET(SAVE_SLICE) THEN loop_slab = *(*(*info).loopsdata).loopslice ELSE loop_slab = *(*(*info).loopsdata).loopslab
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		PRINT,save_message1[KEYWORD_SET(SAVE_SLICE)]
		type = 1
		average_spectrum = ((*(*info).dataparams).spec)[*,(*(*info).dataparams).s] 
		scaling_factor = (*(*info).dataparams).ms
		vertices = *(*(*info).loopsdata).crossloc
		x_coords = REFORM((*(*(*info).loopparams).xp)[0,*])	&	y_coords = REFORM((*(*(*info).loopparams).yp)[0,*])
		x_loop_pts = *(*(*info).loopparams).xr			&	y_loop_pts = *(*(*info).loopparams).yr
		w_loop_pts = *(*(*info).loopparams).w_lpts		&	spect_pos = (*(*info).dataparams).lp
		t_saved = (*(*info).dataparams).t			&	loop_size= (*(*info).loopsdata).loopsize
		crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
		SAVE, crispex_version, type, average_spectrum, scaling_factor, vertices, x_coords, y_coords, x_loop_pts, y_loop_pts, loop_size, loop_slab, spect_pos, t_saved, FILENAME = (*(*info).paths).opath+filename
		PRINT,save_message2[KEYWORD_SET(SAVE_SLICE)]
		SAVE, crispex_version, x_coords, y_coords, x_loop_pts, y_loop_pts, w_loop_pts, spect_pos, t_saved, FILENAME = (*(*info).paths).opath+filename2
		PRINT,'Saved current loop points for later retrieval to: '+STRTRIM(filename2,2)
		WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 1
		WIDGET_CONTROL, (*(*info).ctrlscp).all_saved_loop, SENSITIVE = 1
	ENDIF
END

PRO CRISPEX_SAVE_EXACT_LOOPSLICE, event
; Handles the saving of an exact (i.e. linearly interpolated) timeslice along a loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_EXACT_LOOPSLICE'
	CRISPEX_SAVE_EXACT_LOOPSLAB, event, /SAVE_SLICE
END

PRO CRISPEX_SAVE_EXACT_LOOPSLAB_CHECK, event, SAVE_SLICE=save_slice
; Opens a time estimate warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_EXACT_LOOPSLAB_CHECK'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
		IF ((*(*info).feedbparams).estimate_run EQ 0) THEN BEGIN
			CRISPEX_ESTIMATE_TIME_WINDOW, event
			CRISPEX_ESTIMATE_TIME_CALCULATION, event
			WIDGET_CONTROL, (*(*info).winids).estimatetlb, /DESTROY
			(*(*info).winids).estimatetlb = 0
			(*(*info).winswitch).estimate_win = 0
		ENDIF
		time = (*(*info).feedbparams).estimate_time * N_ELEMENTS(*(*(*info).loopparams).w_lpts)/FLOAT((*(*info).feedbparams).estimate_lx) * (*(*info).dataparams).nlp * (*(*info).dispparams).t_range
		(*(*info).savswitch).cont = 3
		CRISPEX_ESTIMATE_FULL_TIME, time, denom, units
		CRISPEX_SAVE_WARNING_YESNO, event, 'Saving the exact loop slab may take up to about',STRTRIM(CEIL(time/denom),2)+units+'. Do you wish to continue saving?',$
			OK_EVENT='CRISPEX_SAVE_LOOPSL_CONTINUE', CANCEL_EVENT='CRISPEX_SAVE_LOOPSL_ABORT'
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).cont],labels=['Saving procedure']
	ENDIF
END

PRO CRISPEX_SAVE_EXACT_LOOPSLAB, event, SAVE_SLICE=save_slice
; Handles the saving of an exact (i.e. linearly interpolated) timeslab along a loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_EXACT_LOOPSLAB'
	WIDGET_CONTROL,/HOURGLASS
	save_message1 = ['slab','slice']
	IF KEYWORD_SET(SAVE_SLICE) THEN BEGIN
		CRISPEX_SAVE_CHECK_PATH_WRITE, event
		IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
			PRINT,'Saving current exact loop slice...'
			CRISPEX_LOOP_GET_EXACT_SLICE, event, *(*(*info).data).imagedata, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, $
				*(*(*info).loopparams).w_lpts, *(*(*info).loopsdata).exact_loopslice, *(*(*info).loopsdata).exact_crossloc, loopsize, /im
			(*(*info).loopsdata).exact_loopsize = loopsize
			loop_slab = *(*(*info).loopsdata).exact_loopslice
			lp_low = (*(*info).dataparams).lp	&	lp_upp = lp_low
		ENDIF
	ENDIF ELSE BEGIN
		t_0 = SYSTIME(/SECONDS)
		feedback_text = 'Saving current exact loop slab...'
		PRINT,feedback_text
		(*(*info).savparams).lp_orig = (*(*info).dataparams).lp
		lp_low = 0	&	lp_upp = (*(*info).dataparams).nlp-1
		FOR k=lp_low,lp_upp DO BEGIN
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Saving exact loop slab', var=k, maxvar=(*(*info).dataparams).nlp-1, feedback_text=feedback_text+'            ',/destroy_top
			(*(*info).dataparams).lp = k
			CRISPEX_LOOP_GET_EXACT_SLICE, event, *(*(*info).data).imagedata, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, $
				*(*(*info).loopparams).w_lpts, *(*(*info).loopsdata).exact_loopslice, *(*(*info).loopsdata).exact_crossloc, loopsize, /im
			IF (k EQ 0) THEN BEGIN
				*(*(*info).loopsdata).exact_loopslab = *(*(*info).loopsdata).exact_loopslice
			ENDIF ELSE *(*(*info).loopsdata).exact_loopslab = [[[*(*(*info).loopsdata).exact_loopslab]], [[*(*(*info).loopsdata).exact_loopslice]]]
			t_1 = SYSTIME(/SECONDS)
			CRISPEX_ESTIMATE_FULL_TIME_RUNNING, k+1, (*(*info).dataparams).nlp, t_0, t_1, denom, units, accumsectime, totalsectime
			feedback_text = STRTRIM(k+1,2)+'/'+STRTRIM((*(*info).dataparams).nlp,2)+' slices extracted in '+STRTRIM(STRING(accumsectime/denom,FORMAT='(3(F9.2,x))'),2)+'/'+$
				STRTRIM(STRING(totalsectime/denom,FORMAT='(3(F9.2,x))'),2)+units
		ENDFOR
		(*(*info).dataparams).lp = (*(*info).savparams).lp_orig
		loop_slab= *(*(*info).loopsdata).exact_loopslab
	ENDELSE
	(*(*info).loopsdata).exact_loopsize = loopsize
	type = 0
	average_spectrum = ((*(*info).dataparams).spec)[*,(*(*info).dataparams).s] 
	scaling_factor = (*(*info).dataparams).ms
	vertices = *(*(*info).loopsdata).exact_crossloc
	x_coords = REFORM((*(*(*info).loopparams).xp)[0,*])		&	y_coords = REFORM((*(*(*info).loopparams).yp)[0,*])
	x_loop_pts = *(*(*info).loopparams).xr				&	y_loop_pts = *(*(*info).loopparams).yr
	spect_pos = (*(*info).dataparams).lp				&	loop_size= (*(*info).loopsdata).exact_loopsize
	spect_pos_low = lp_low						&	spect_pos_upp = lp_upp
	t_low = (*(*info).dispparams).t_low				&	t_upp = (*(*info).dispparams).t_upp
	t_saved = (*(*info).dataparams).t
	crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=filename, /tlab, ext='csav'
	IF ~KEYWORD_SET(SAVE_SLICE) THEN CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event
	SAVE, crispex_version, type, average_spectrum, scaling_factor, vertices, x_coords, y_coords, x_loop_pts, y_loop_pts, loop_size, loop_slab, spect_pos, lp_low, lp_upp, t_saved, t_low, t_upp, $
		FILENAME = (*(*info).paths).opath+filename
	PRINT,'Saving current exact loop '+save_message1[KEYWORD_SET(SAVE_SLICE)]+' done. Saved data to: '+STRTRIM(filename,2)
END

PRO CRISPEX_SAVE_LOOPSL_CONTINUE, event
; Handles continue events from the loopslice/slab warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_LOOPSL_CONTINUE'
	IF ((*(*info).savswitch).cont EQ 1) THEN CRISPEX_SAVE_RETRIEVE_LOOPSLAB, event, /SAVE_SLICE		; Save diagram from saved path at single spectral position
	IF ((*(*info).savswitch).cont EQ 2) THEN CRISPEX_SAVE_RETRIEVE_LOOPSLAB, event 				; Save diagram from saved path for all spectral positions
	IF ((*(*info).savswitch).cont EQ 3) THEN CRISPEX_SAVE_EXACT_LOOPSLAB, event 				; Save diagram for all spectral positions
	IF ((*(*info).savswitch).cont EQ 4) THEN CRISPEX_SAVE_RETRIEVE_DET_LOOPSLAB, event, /SAVE_SLICE		; Save diagram from saved detection at single spectral position
	IF ((*(*info).savswitch).cont EQ 5) THEN CRISPEX_SAVE_RETRIEVE_DET_LOOPSLAB, event 			; Save diagram from saved detection at all spectral positions
;	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SAVE_LOOPSL_ABORT, event
; Handles abort event from the loopslice/slab warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_EXACT_LOOPSL_ABORT'
	(*(*info).savswitch).cont = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).cont],labels=['Saving procedure']
	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SAVE_RETRIEVE_LOOPSLAB, event, SAVE_SLICE=save_slice
; Handles the saving of an exact (i.e. linearly interpolated) timeslab along a retrieved loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_RETRIEVE_LOOPSLAB'
	(*(*info).savparams).lp_orig = (*(*info).dataparams).lp
	(*(*info).savparams).lp_ref_orig = (*(*info).dataparams).lp_ref
	nfiles = (*(*info).retrparams).retrieve_filecount * CEIL((*(*info).savswitch).imref_only/2.)
	pass = 0L
	loopdet = ([2,1])[KEYWORD_SET(SAVE_SLICE)]
	IF KEYWORD_SET(SAVE_SLICE) THEN BEGIN
		nlp_pos = 1				&	refnlp_pos = 1
		save_message = 'slice'
	ENDIF ELSE BEGIN
		nlp_pos = (*(*info).dataparams).nlp	&	refnlp_pos = (*(*info).dataparams).refnlp
		save_message = 'slab'
	ENDELSE
	IF ((*(*info).savswitch).imref_only EQ 1) THEN maxpass = nlp_pos * nfiles ELSE IF ((*(*info).savswitch).imref_only EQ 2) THEN maxpass = refnlp_pos * nfiles ELSE $
		maxpass = (nlp_pos + refnlp_pos) * (*(*info).retrparams).retrieve_filecount
	t_0 = SYSTIME(/SECONDS)
	WIDGET_CONTROL,/HOURGLASS
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	FOR i=0,nfiles-1 DO BEGIN
		t_1 = SYSTIME(/SECONDS)
		fstr = (*(*(*info).retrparams).retrieve_files)[(i MOD (*(*info).retrparams).retrieve_filecount)]
		IF ((*(*info).dataparams).reffilename NE '') THEN CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).imref_only, (*(*info).dataparams).imfilename, $
			reffilename=(*(*info).dataparams).reffilename, *(*(*info).data).imagedata, refdata=*(*(*info).data).refdata, (*(*info).paths).opath, i, (*(*info).retrparams).retrieve_filecount, filename, data, nonlp,$
			imref, whichdata, fstr=fstr, loopdet=loopdet ELSE $
			CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).imref_only, (*(*info).dataparams).imfilename, *(*(*info).data).imagedata, (*(*info).paths).opath, i, $
				(*(*info).retrparams).retrieve_filecount, filename, data, nonlp, imref, whichdata, fstr=fstr, loopdet=loopdet
		RESTORE,(*(*(*info).retrparams).retrieve_files)[(i MOD (*(*info).retrparams).retrieve_filecount)]
		*(*(*info).loopparams).xp = x_coords		&	*(*(*info).loopparams).yp = y_coords
		*(*(*info).loopparams).xr = x_loop_pts		&	*(*(*info).loopparams).yr = y_loop_pts
		*(*(*info).loopparams).w_lpts = w_loop_pts
		IF (N_ELEMENTS(T_SAVED) NE 1) THEN t_saved = (*(*info).dataparams).t
		IF KEYWORD_SET(SAVE_SLICE) THEN BEGIN
			lp_low = spect_pos	&	lp_upp = lp_low
		ENDIF ELSE BEGIN
			lp_low = 0
			IF (imref EQ 1) THEN BEGIN
				lp_upp = (*(*info).dataparams).nlp-1 	&	spect_pos = (*(*info).savparams).lp_orig
			ENDIF ELSE BEGIN
				lp_upp = (*(*info).dataparams).refnlp-1	&	spect_pos = (*(*info).savparams).lp_ref_orig
			ENDELSE
		ENDELSE
		FOR k=lp_low,lp_upp DO BEGIN
			pass += 1L
			part = STRTRIM(pass,2)+'/'+STRTRIM(maxpass,2)
			IF (pass EQ 1) THEN feedback_text = ': saving retrieved exact loop '+save_message+'...               '
			CRISPEX_UPDATE_USER_FEEDBACK, event, title='Saving retrieved loop '+save_message+'(s)', var=pass-1, maxvar=maxpass, feedback_text=part+feedback_text, /destroy_top
			IF (imref EQ 1) THEN (*(*info).dataparams).lp = k ELSE (*(*info).dataparams).lp_ref = k
			CRISPEX_LOOP_GET_EXACT_SLICE, event, data, *(*(*info).loopparams).xr, *(*(*info).loopparams).yr, *(*(*info).loopparams).xp, *(*(*info).loopparams).yp, *(*(*info).loopparams).w_lpts, $
				*(*(*info).loopsdata).exact_loopslice, *(*(*info).loopsdata).exact_crossloc, loopsize, no_nlp=nonlp, im=whichdata
			IF ~KEYWORD_SET(SAVE_SLICE) THEN BEGIN
				IF (k EQ lp_low) THEN *(*(*info).loopsdata).exact_loopslab = *(*(*info).loopsdata).exact_loopslice ELSE $
					*(*(*info).loopsdata).exact_loopslab = [[[*(*(*info).loopsdata).exact_loopslab]], [[*(*(*info).loopsdata).exact_loopslice]]]
			ENDIF
			t_1 = SYSTIME(/SECONDS)
			CRISPEX_ESTIMATE_FULL_TIME_RUNNING, pass, maxpass, t_0, t_1, denom, units, accumsectime, totalsectime
			feedback_text = ' slices extracted in '+STRTRIM(STRING(accumsectime/denom,FORMAT='(3(F9.2,x))'),2)+'/'+STRTRIM(STRING(totalsectime/denom,FORMAT='(3(F9.2,x))'),2)+units
		ENDFOR
		IF KEYWORD_SET(SAVE_SLICE) THEN loop_slab = *(*(*info).loopsdata).exact_loopslice ELSE loop_slab= *(*(*info).loopsdata).exact_loopslab
		(*(*info).loopsdata).exact_loopsize = loopsize
		type = 0
		IF (imref EQ 1) THEN BEGIN
			average_spectrum = ((*(*info).dataparams).spec)[*,(*(*info).dataparams).s] 
			scaling_factor = (*(*info).dataparams).ms
		ENDIF ELSE BEGIN
			average_spectrum = ((*(*info).dataparams).refspec)
			scaling_factor = (*(*info).dataparams).refms
		ENDELSE
		vertices = *(*(*info).loopsdata).exact_crossloc
		x_coords = REFORM((*(*(*info).loopparams).xp)[0,*])		&	y_coords = REFORM((*(*(*info).loopparams).yp)[0,*])
		x_loop_pts = *(*(*info).loopparams).xr				&	y_loop_pts = *(*(*info).loopparams).yr
		spect_pos_low = lp_low						&	spect_pos_upp = lp_upp
		loop_size= (*(*info).loopsdata).exact_loopsize
		t_low = (*(*info).dispparams).t_low				&	t_upp = (*(*info).dispparams).t_upp
		crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
		singlefilename = (STRSPLIT(filename, PATH_SEP(), /EXTRACT))[N_ELEMENTS(STRSPLIT(filename, PATH_SEP(), /EXTRACT))-1]
		SAVE, crispex_version, type, average_spectrum, scaling_factor, vertices, x_coords, y_coords, x_loop_pts, y_loop_pts, loop_size, loop_slab, spect_pos, lp_low, lp_upp, t_saved, t_low, t_upp, $
			FILENAME = (*(*info).paths).opath+singlefilename 
		PRINT,'Saving retrieved exact loop '+save_message+' done. Saved data to: '+STRTRIM(singlefilename,2)
		IF (*(*info).savswitch).delete_clsav THEN BEGIN
			IF ((*(*info).savswitch).imref_only LE 2) THEN SPAWN,'rm '+STRTRIM((*(*(*info).retrparams).retrieve_files)[i],2) ELSE BEGIN
				IF (i GE (*(*info).retrparams).retrieve_filecount) THEN SPAWN,'rm '+STRTRIM((*(*(*info).retrparams).retrieve_files)[(i MOD (*(*info).retrparams).retrieve_filecount)],2)
			ENDELSE
		ENDIF
	ENDFOR
	CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event
	(*(*info).dataparams).lp = (*(*info).savparams).lp_orig
	(*(*info).dataparams).lp_ref = (*(*info).savparams).lp_ref_orig
END

PRO CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, imref_only, imfilename, reffilename=reffilename, imdata, refdata=refdata, opath, var, endoflist, outputfilename, outputdata, outputnonlp, outputimref, outputwhichdata, $
	index=index, fstr=fstr, loopdet=loopdet
; Determines the filename of the retrieved loop path or detection timeslice to be saved
	; loopdet = 1 > loop, slice
	; loopdet = 2 > loop, slab
	; loopdet = 3 > det, slice
	; loopdet = 4 > det, slab
	; outputimref = 1 > main
	; outputimref = 2 > reference
	IF (imref_only EQ 1) THEN BEGIN			; Only save from main cube
		IF (loopdet LE 2) THEN CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=fstr, outfilename=outputfilename, ext='csav', /exch_ext ELSE $
			CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=fstr, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
		outputnonlp = (loopdet EQ 1)
		outputimref = 1
		outputdata = imdata
		outputwhichdata = 1
	ENDIF ELSE IF (imref_only EQ 2) THEN BEGIN	; Only save from reference cube
		IF (loopdet LE 2) THEN BEGIN
			clstr = STRSPLIT(fstr,'_',/EXTRACT)
			datestamp = clstr[N_ELEMENTS(clstr)-2]
			timestamp = (STRSPLIT(clstr[N_ELEMENTS(clstr)-1],'.',/EXTRACT))[0]
			CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=reffilename, outfilename=outputfilename, ext='csav', import_id = datestamp+'_'+timestamp
		ENDIF ELSE CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=reffilename, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
		outputnonlp = (loopdet EQ 3)
		outputimref = 2
		outputdata = refdata
		outputwhichdata = 0
	ENDIF ELSE BEGIN				; Save from both main and reference cube
		IF (var LT endoflist) THEN BEGIN
			IF (loopdet LE 2) THEN CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=fstr, outfilename=outputfilename, ext='csav', /exch_ext ELSE $
				CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=fstr, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
			outputnonlp = (loopdet EQ 1)
			outputimref = 1
			outputdata = imdata
			outputwhichdata = 1
		ENDIF ELSE BEGIN
			IF (loopdet LE 2) THEN BEGIN
				clstr = STRSPLIT(fstr,'_',/EXTRACT)
				datestamp = clstr[N_ELEMENTS(clstr)-2]
				timestamp = (STRSPLIT(clstr[N_ELEMENTS(clstr)-1],'.',/EXTRACT))[0]
				CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=reffilename, outfilename=outputfilename, ext='csav', import_id = datestamp+'_'+timestamp
			ENDIF ELSE CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=reffilename, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
			outputnonlp = (loopdet EQ 3)
			outputimref = 2
			outputdata = refdata
			outputwhichdata = 0
		ENDELSE
	ENDELSE
END

PRO CRISPEX_SAVE_RETRIEVE_DET_LOOPSLAB, event, SAVE_SLICE=save_slice
; Handles the saving of an exact (i.e. linearly interpolated) timeslab along a retrieved detection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_RETRIEVE_DET_LOOPSLAB'
	nfiles = (*(*info).detparams).nr_sel_loops * CEIL((*(*info).savswitch).det_imref_only/2.)
	pass = 0
	loopdet = ([4,3])[KEYWORD_SET(SAVE_SLICE)]
	IF KEYWORD_SET(SAVE_SLICE) THEN BEGIN
		lp_dn = (*(*info).dataparams).lp	&	lp_up = lp_dn
		IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
			lp_ref_dn = (*(*info).dataparams).lp_ref	&	lp_ref_up = lp_ref_dn
			(*(*info).detparams).lp_ref_dn = lp_ref_dn	&	(*(*info).detparams).lp_ref_up = lp_ref_up
		ENDIF
		save_message = 'slice'
	ENDIF ELSE BEGIN
		IF ((*(*info).savswitch).pos_dets EQ 1) THEN BEGIN
			lp_dn = (*(*info).dispparams).lp_first	&	lp_up = (*(*info).dispparams).lp_last
			(*(*info).detparams).lp_dn = lp_dn	&	(*(*info).detparams).lp_up = lp_up
			IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
				lp_ref_dn = (*(*info).dispparams).lp_ref_first	&	 lp_ref_up = (*(*info).dispparams).lp_ref_last
				(*(*info).detparams).lp_ref_dn = lp_ref_dn	&	(*(*info).detparams).lp_ref_up = lp_ref_up
			ENDIF
		ENDIF ELSE BEGIN
			lp_dn = (*(*info).detparams).lp_dn	&	lp_up = (*(*info).detparams).lp_up
			IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
				lp_ref_dn = (*(*info).detparams).lp_ref_dn	&	 lp_ref_up = (*(*info).detparams).lp_ref_up
			ENDIF
		ENDELSE
		save_message = 'slab'
	ENDELSE
	IF ((*(*info).savswitch).det_imref_only EQ 1) THEN maxpass = (lp_up-lp_dn+1) * (*(*info).detparams).width * nfiles ELSE $
		IF ((*(*info).savswitch).det_imref_only EQ 2) THEN maxpass = (lp_ref_up-lp_ref_dn+1) * (*(*info).detparams).width * nfiles ELSE $
			maxpass = ((lp_up-lp_dn+1)  + (lp_ref_up-lp_ref_dn+1)) * (*(*info).detparams).width * (*(*info).detparams).nr_sel_loops
	lower_t = (*(*info).dispparams).t_low
	upper_t = (*(*info).dispparams).t_upp
	(*(*info).savparams).lp_orig = (*(*info).dataparams).lp
	IF ((*(*info).savparams).lp_orig LT (*(*info).detparams).lp_dn) THEN lp_im_saved = (*(*info).detparams).lp_dn ELSE $
		IF ((*(*info).savparams).lp_orig GT (*(*info).detparams).lp_up) THEN lp_im_saved = (*(*info).detparams).lp_up ELSE lp_im_saved = (*(*info).savparams).lp_orig
	IF ((*(*info).savswitch).det_imref_only GE 2) THEN BEGIN
		(*(*info).savparams).lp_ref_orig = (*(*info).dataparams).lp_ref
		IF ((*(*info).savparams).lp_ref_orig LT (*(*info).detparams).lp_ref_dn) THEN lp_ref_saved = (*(*info).detparams).lp_ref_dn ELSE $
			IF ((*(*info).savparams).lp_ref_orig GT (*(*info).detparams).lp_ref_up) THEN lp_ref_saved = (*(*info).detparams).lp_ref_up ELSE lp_ref_saved = (*(*info).savparams).lp_ref_orig
	ENDIF
	t_0 = SYSTIME(/SECONDS)
	FOR i=0,nfiles-1 DO BEGIN
		t_1 = SYSTIME(/SECONDS)
		WIDGET_CONTROL,/HOURGLASS
		part = STRTRIM(i+1,2)+'/'+STRTRIM(nfiles,2)
		idx = (*(*(*info).detparams).sel_loops)[(i MOD (*(*info).detparams).nr_sel_loops)]
		IF ((*(*info).dataparams).reffilename NE '') THEN CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).det_imref_only, (*(*info).dataparams).imfilename, $
			reffilename=(*(*info).dataparams).reffilename, *(*(*info).data).imagedata, refdata=*(*(*info).data).refdata, (*(*info).paths).opath, i, (*(*info).detparams).nr_sel_loops, filename, data, nonlp, $
			detimref, detwhichdata, fstr=fstr, loopdet=loopdet, index=idx ELSE $
			CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).det_imref_only, (*(*info).dataparams).imfilename, *(*(*info).data).imagedata, (*(*info).paths).opath, i, $
				(*(*info).detparams).nr_sel_loops, filename, data, nonlp, detimref, detwhichdata, fstr=fstr, loopdet=loopdet, index=idx
		w_lpts = INDGEN((*(*(*info).detparams).lrsizes)[idx+1] - (*(*(*info).detparams).lrsizes)[idx])
		tmp_loopslab = 0
		t_det = (*(*(*info).detparams).t_restored)[idx]
		(*(*info).dispparams).t_low = t_det - (*(*info).detparams).delta_t_dn > (*(*info).dispparams).t_first
		(*(*info).dispparams).t_upp = t_det + (*(*info).detparams).delta_t_up < (*(*info).dispparams).t_last
		delta = FLOOR((*(*info).detparams).width/2.)
		IF (detimref EQ 1) THEN lp_saved = lp_im_saved ELSE BEGIN
			lp_saved = lp_ref_saved
			lp_dn = (*(*info).detparams).lp_ref_dn 	&	lp_up = (*(*info).detparams).lp_ref_up
		ENDELSE
		FOR k=((*(*info).detparams).mid-delta),((*(*info).detparams).mid+delta) DO BEGIN
			xlp_det = (*(*(*info).detparams).xlp_restored)[(0+(*(*(*info).detparams).lpsizes)[idx]):((*(*(*info).detparams).lpsizes)[idx+1]-1)]
			ylp_det = (*(*(*info).detparams).ylp_restored)[(0+(*(*(*info).detparams).lpsizes)[idx]):((*(*(*info).detparams).lpsizes)[idx+1]-1)]
			xlr_det = (*(*(*info).detparams).xlr_restored)[(0+(*(*(*info).detparams).lrsizes)[idx]):((*(*(*info).detparams).lrsizes)[idx+1]-1), k]
			ylr_det = (*(*(*info).detparams).ylr_restored)[(0+(*(*(*info).detparams).lrsizes)[idx]):((*(*(*info).detparams).lrsizes)[idx+1]-1), k]
			no_nlp = (lp_dn EQ lp_up)
			FOR m=lp_dn,lp_up DO BEGIN
				pass += 1L
				part = 'Position '+STRTRIM(k,2)+': '+STRTRIM(m-lp_dn+1,2)+'/'+STRTRIM(lp_up-lp_dn+1,2)+'('+STRTRIM(pass,2)+'/'+STRTRIM(maxpass,2)+')'
				IF (pass EQ 1) THEN feedback_text = ': retrieving detection '+STRTRIM(idx,2)+' and saving exact loop '+save_message+'...'
				CRISPEX_UPDATE_USER_FEEDBACK, event, title='Saving retrieved detection loop '+save_message+'(s)', var=pass-1, maxvar=maxpass, feedback_text=part+feedback_text, /destroy_top
				IF (detimref EQ 1) THEN (*(*info).dataparams).lp = m ELSE (*(*info).dataparams).lp_ref = m
				CRISPEX_LOOP_GET_EXACT_SLICE, event, data, xlr_det, ylr_det, xlp_det, ylp_det, w_lpts, exact_loopslice, exact_crossloc, exact_loopsize, im=detwhichdata, no_nlp=no_nlp
				IF (m EQ lp_dn) THEN exact_loopslab = exact_loopslice ELSE exact_loopslab = [[[exact_loopslab]], [[exact_loopslice]]]
				t_1 = SYSTIME(/SECONDS)
				CRISPEX_ESTIMATE_FULL_TIME_RUNNING, pass, maxpass, t_0, t_1, denom, units, accumsectime, totalsectime
				feedback_text = ' slices extracted in '+STRTRIM(STRING(accumsectime/denom,FORMAT='(3(F9.2,x))'),2)+'/'+STRTRIM(STRING(totalsectime/denom,FORMAT='(3(F9.2,x))'),2)+units
			ENDFOR
			tmp_loopslab += exact_loopslab 
			IF (k EQ (*(*info).detparams).mid) THEN BEGIN
				crossloc = exact_crossloc
				loopsize = exact_loopsize
			ENDIF
		ENDFOR
		IF KEYWORD_SET(SAVE_SLICE) THEN BEGIN
			*(*(*info).loopsdata).exact_loopslice = tmp_loopslab/(*(*info).detparams).width	&	loop_slab = *(*(*info).loopsdata).exact_loopslice
		ENDIF ELSE BEGIN
			*(*(*info).loopsdata).exact_loopslab = tmp_loopslab/(*(*info).detparams).width	&	loop_slab = *(*(*info).loopsdata).exact_loopslab
		ENDELSE
		spect_pos = lp_saved	&	spect_pos_low = lp_dn	&	spect_pos_upp = lp_up
		*(*(*info).loopsdata).exact_crossloc = crossloc
		(*(*info).loopsdata).exact_loopsize = loopsize
		type = 0
		IF (detimref EQ 1) THEN BEGIN
			average_spectrum = ((*(*info).dataparams).spec)[*,(*(*info).dataparams).s] 
			scaling_factor = (*(*info).dataparams).ms
		ENDIF ELSE BEGIN
			average_spectrum = ((*(*info).dataparams).refspec)
			scaling_factor = (*(*info).dataparams).refms
		ENDELSE
		vertices = *(*(*info).loopsdata).exact_crossloc	&	loop_size= (*(*info).loopsdata).exact_loopsize
		x_coords = xlp_det				&	y_coords = ylp_det
		x_loop_pts = xlr_det				&	y_loop_pts = ylr_det
		t_saved = t_det					&	t_low = (*(*info).dispparams).t_low	&	t_upp = (*(*info).dispparams).t_upp
		crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
		singlefilename = (STRSPLIT(filename, PATH_SEP(), /EXTRACT))[N_ELEMENTS(STRSPLIT(filename, PATH_SEP(), /EXTRACT))-1]
		SAVE, crispex_version, type, average_spectrum, scaling_factor, vertices, x_coords, y_coords, x_loop_pts, y_loop_pts, loop_size, loop_slab, spect_pos, spect_pos_low, spect_pos_upp, t_saved, t_low, t_upp, $
			FILENAME = (*(*info).paths).opath+singlefilename 
		PRINT,'Saving retrieved exact loop '+save_message+' done. Saved data to: '+STRTRIM(singlefilename,2)
	ENDFOR
	(*(*info).dispparams).t_low = lower_t
	(*(*info).dispparams).t_upp = upper_t
	(*(*info).dataparams).lp = (*(*info).savparams).lp_orig
	IF ((*(*info).savswitch).det_imref_only GE 2) THEN (*(*info).dataparams).lp_ref = (*(*info).savparams).lp_ref_orig
	CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event
END

;================================================================================= OUTPUT SAVE SETTINGS AND FILENAME CHECK
PRO CRISPEX_SAVE_CHECK_PATH_WRITE, event
; Checks whether the user has write permissions in the current opath
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_CHECK_PATH_WRITE'
	(*(*info).paths).opath_write = FILE_TEST((*(*info).paths).opath, /WRITE)
	IF ((*(*info).paths).opath_write EQ 0) THEN BEGIN
		CRISPEX_WINDOW_OK, event, 'ERROR!','You appear not to have write permissions to the current','output directory ('+(*(*info).paths).opath+').',$
			'Please change the path before continuing saving.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb, /BLOCK
		(*(*info).winids).errtlb = tlb
	END
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).opath_write],labels=['Output path writeable']
END

PRO CRISPEX_SAVE_SET_IPATH, event 
; Sets the output path for all saving procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_SET_IPATH'
	thispath = DIALOG_PICKFILE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Set input path', /DIRECTORY, PATH = (*(*info).paths).ipath)
	IF (thispath EQ '') THEN RETURN ELSE (*(*info).paths).ipath = thispath
	IF ((*(*info).winids).savewintlb GT 0) THEN WIDGET_CONTROL, (*(*info).ctrlssav).path_textlab, SET_VALUE = STRTRIM((*(*info).paths).ipath,2)
	IF (event.TOP EQ (*(*info).winids).preftlb) THEN (*(*info).prefs).tmp_prefipath = (*(*info).paths).ipath
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).ipath],labels=['Input path']
END

PRO CRISPEX_SAVE_SET_OPATH, event 
; Sets the output path for all saving procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_SET_OPATH'
	thispath = DIALOG_PICKFILE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Set output path', /DIRECTORY, PATH = (*(*info).paths).opath)
	IF (thispath EQ '') THEN RETURN ELSE (*(*info).paths).opath = thispath
	IF ((*(*info).winids).savewintlb GT 0) THEN WIDGET_CONTROL, (*(*info).ctrlssav).path_textlab, SET_VALUE = STRTRIM((*(*info).paths).opath,2)
	IF ((*(*info).winids).saveoptwintlb GT 0) THEN WIDGET_CONTROL, (*(*info).ctrlssav).savopt_path_textlab, SET_VALUE = STRTRIM((*(*info).paths).opath,2)
	IF (event.TOP EQ (*(*info).winids).preftlb) THEN (*(*info).prefs).tmp_prefopath = (*(*info).paths).opath
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).opath],labels=['Output path']
END

PRO CRISPEX_SAVE_GET_FILENAME, event, title, standard_filename, ok_event, session_save=session_save, SNAPSHOT=snapshot
; Gets the filename for save routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_GET_FILENAME'
	IF ((*(*info).winids).savewintlb NE 0) THEN BEGIN
		CRISPEX_WINDOW_OK, event, 'WARNING!','You are currently already saving output. Please finish','saving or discard the current saving procedure first,',$
			'before starting a new saving procedure.', OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).warntlb = tlb
		RETURN
	END
	fulltitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+STRTRIM(title,2)
	base = WIDGET_BASE(TITLE = title, GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	path_base = WIDGET_BASE(message_base, /ROW)
	path_label = WIDGET_LABEL(path_base, VALUE = 'Current path:')
	(*(*info).ctrlssav).path_textlab = WIDGET_LABEL(path_base, VALUE = STRTRIM((*(*info).paths).opath,2), /DYNAMIC_RESIZE)
	filename_base = WIDGET_BASE(message_base, /ROW)
	text_label1 = WIDGET_LABEL(filename_base, VALUE = 'Filename (w/o extension):')
	(*(*info).savparams).filename_text = WIDGET_TEXT(filename_base, VALUE = standard_filename,  /EDITABLE, XSIZE = 40)
	WIDGET_CONTROL, (*(*info).savparams).filename_text, SET_TEXT_SELECT = [0,STRLEN(standard_filename)]
	button_base = WIDGET_BASE(disp,COLUMN=3,/GRID_LAYOUT,/ALIGN_CENTER)
	IF KEYWORD_SET(SESSION_SAVE) THEN change_path_but = WIDGET_BUTTON(button_base, VALUE = 'Change path', EVENT_PRO = 'CRISPEX_SAVE_SET_OPATH') ELSE $
		options_but = WIDGET_BUTTON(button_base, VALUE = 'Options...', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS')
	no_but = WIDGET_BUTTON(button_base, VALUE = 'Cancel', EVENT_PRO = 'CRISPEX_CLOSE_EVENT_WINDOW')
	yes_but = WIDGET_BUTTON(button_base, VALUE = 'OK' , EVENT_PRO = ok_event)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).savewintlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).savewintlb],labels=['savewintlb']
END

PRO CRISPEX_SAVE_CHECK_FILENAME, event, extension, ok_event, midtension=midtension
; Checks the save filename for validity and overwrite problems
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_CHECK_FILENAME'
	csesfiles = FILE_SEARCH((*(*info).paths).opath+"*"+extension, COUNT = csesfilecount)
	WIDGET_CONTROL, (*(*info).savparams).filename_text, GET_VALUE = session_filename
	IF (N_ELEMENTS(midtension) GT 0) THEN midtension = midtension ELSE midtension = ''
	full_session_filename = session_filename+midtension+'.'+extension
	compressedfilename = STRCOMPRESS(full_session_filename, /REMOVE_ALL)
	existing = WHERE(STRPOS(csesfiles,(*(*info).paths).opath+full_session_filename) EQ 0)
	IF (existing EQ -1) AND (session_filename NE '') AND (compressedfilename EQ full_session_filename)  AND ((*(*info).paths).opath_write EQ 1) THEN CRISPEX_SAVE_CONTINUE, event, session_filename $
	ELSE IF (existing EQ -1) AND (session_filename EQ '') OR (compressedfilename NE full_session_filename) THEN BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','Invalid filename. Please enter a filename of','at least one character and without any white spaces.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDIF ELSE IF ((*(*info).paths).opath_write EQ 0) THEN CRISPEX_SAVE_CHECK_PATH_WRITE, event $
	ELSE CRISPEX_SAVE_WARNING_YESNO, event,'The file already exists. Do you wish to','continue and overwrite this file?', OK_EVENT=ok_event, CANCEL_EVENT='CRISPEX_CLOSE_EVENT_WINDOW'
END

PRO CRISPEX_SAVE_CONTINUE, event, session_filename
; Selects the continuation of saving session, JPEG or MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_CONTINUE'
	IF ((*(*info).savparams).savpro EQ 'SESSION') THEN CRISPEX_SESSION_SAVE, event, session_filename
	IF ((STRPOS((*(*info).savparams).savpro,'FRAMES') NE -1) OR STRCMP((*(*info).savparams).savpro,'MPEG')) THEN CRISPEX_SAVE_FRAMES_SAVE, event, session_filename ELSE $
		IF (STRPOS((*(*info).savparams).savpro,'LINESCAN') NE -1) THEN CRISPEX_SAVE_LINESCAN_SAVE, event, session_filename
END

PRO CRISPEX_SAVE_WARNING_YESNO, event, warningmessage1, warningmessage2, warningmessage3, ok_event=ok_event, cancel_event=cancel_event
; Creates the loopslice/slab warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_WARNING_YESNO'
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': WARNING!', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	IF (N_ELEMENTS(warningmessage1) GT 0) THEN text_label1 = WIDGET_LABEL(message_base, VALUE = warningmessage1)
	IF (N_ELEMENTS(warningmessage2) GT 0) THEN text_label2 = WIDGET_LABEL(message_base, VALUE = warningmessage2)
	IF (N_ELEMENTS(warningmessage3) GT 0) THEN text_label3 = WIDGET_LABEL(message_base, VALUE = warningmessage3)
	button_base = WIDGET_BASE(disp,/ROW,/ALIGN_CENTER)
	no_but = WIDGET_BUTTON(button_base, VALUE = 'No, cancel', EVENT_PRO = cancel_event)
	yes_but = WIDGET_BUTTON(button_base, VALUE = 'Yes, continue' , EVENT_PRO = ok_event)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
END

;================================================================================= OUTPUT SAVE PROCEDURES
PRO CRISPEX_SAVE_JPEG_SNAPSHOT, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_JPEG_SNAPSHOT'
	(*(*info).savparams).savpro = 'JPEG_FRAMES'
	(*(*info).savparams).snapshot = 1
	CRISPEX_SAVE_FRAMES, event
END	

PRO CRISPEX_SAVE_JPEG_ALL_FRAMES, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_JPEG_ALL_FRAMES'
	(*(*info).savparams).savpro = 'JPEG_FRAMES'
	(*(*info).savparams).snapshot = 0
	CRISPEX_SAVE_FRAMES, event
END

PRO CRISPEX_SAVE_JPEG_LINESCAN, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_JPEG_LINESCAN'
	(*(*info).savparams).savpro = 'JPEG_LINESCAN'
	CRISPEX_SAVE_LINESCAN, event
END	

PRO CRISPEX_SAVE_PNG_SNAPSHOT, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_PNG_SNAPSHOT'
	(*(*info).savparams).savpro = 'PNG_FRAMES'
	(*(*info).savparams).snapshot = 1
	CRISPEX_SAVE_FRAMES, event
END	

PRO CRISPEX_SAVE_PNG_ALL_FRAMES, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_PNG_ALL_FRAMES'
	(*(*info).savparams).savpro = 'PNG_FRAMES'
	(*(*info).savparams).snapshot = 0
	CRISPEX_SAVE_FRAMES, event
END

PRO CRISPEX_SAVE_PNG_LINESCAN, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_PNG_LINESCAN'
	(*(*info).savparams).savpro = 'PNG_LINESCAN'
	CRISPEX_SAVE_LINESCAN, event
END	

PRO CRISPEX_SAVE_CHECK, event
; Checks the jpeg series save filename for validity and overwrite problems
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_CHECK'
	IF ((*(*info).savparams).savpro EQ 'MPEG') THEN BEGIN
		extension = 'mpg'	&	midtension = ''
	ENDIF ELSE IF (STRPOS((*(*info).savparams).savpro,'FRAMES') NE -1) THEN BEGIN		; Check whether saving procedure contains 'FRAMES' (test = 1) or 'LINESCAN' (test = 0)
		IF (*(*info).savparams).snapshot THEN midtension = '' ELSE BEGIN
			nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
			lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
			midtension = '_'+lp_id
		ENDELSE
		IF ((*(*info).savparams).savpro EQ 'JPEG_FRAMES') THEN extension = 'jpg'
		IF ((*(*info).savparams).savpro EQ 'PNG_FRAMES') THEN extension = 'png'
	ENDIF ELSE BEGIN
		ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
		nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
		t_id = 't'+STRING((*(*info).dataparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
		lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
		midtension = '_'+t_id+'_'+lp_id
		IF ((*(*info).savparams).savpro EQ 'JPEG_LINESCAN') THEN extension = 'jpg'
		IF ((*(*info).savparams).savpro EQ 'PNG_LINESCAN') THEN extension = 'png'
	ENDELSE
	CRISPEX_SAVE_CHECK_FILENAME, event, extension, 'CRISPEX_SAVE_OVER', midtension = midtension
END

PRO CRISPEX_SAVE_OVER, event
; Handles the overwriting and activates the subsequent saving of the jpeg frames
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OVER'
	WIDGET_CONTROL, (*(*info).savparams).filename_text, GET_VALUE = filename
	IF ((*(*info).savparams).savpro EQ 'MPEG') THEN CRISPEX_SAVE_MPEG_SAVE, event, filename ELSE $
		IF (STRPOS((*(*info).savparams).savpro,'FRAMES') NE -1) THEN CRISPEX_SAVE_FRAMES_SAVE, event, filename ELSE CRISPEX_SAVE_LINESCAN_SAVE, event, filename
	IF ((*(*info).winids).saveoptwintlb GT 0) THEN BEGIN
		WIDGET_CONTROL, (*(*info).winids).saveoptwintlb,/DESTROY
		(*(*info).winids).saveoptwintlb = 0
	ENDIF
	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SAVE_FRAMES, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_FRAMES'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	label = ['series','snapshot']
	IF (*(*info).savparams).snapshot THEN BEGIN
		ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
		nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
		t_id = 't'+STRING((*(*info).dataparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
		lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
		CRISPEX_SAVE_CHECK_PATH_WRITE, event
		CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename, import_id=lp_id+'_'+t_id
	ENDIF ELSE CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename
	IF ((*(*info).savparams).savpro EQ 'JPEG_FRAMES') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save JPEG '+label[KEYWORD_SET(SNAPSHOT)], standardfilename, 'CRISPEX_SAVE_CHECK'
	IF ((*(*info).savparams).savpro EQ 'PNG_FRAMES') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save PNG '+label[KEYWORD_SET(SNAPSHOT)], standardfilename, 'CRISPEX_SAVE_CHECK'
END

PRO CRISPEX_SAVE_FRAMES_SAVE, event, supplied_filename
; Handles the saving of a series (between temporal boundaries) of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_FRAMES_SAVE'
	ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
	nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
	t_before = (*(*info).dataparams).t
	IF (*(*info).savparams).snapshot THEN BEGIN
		t_low = (*(*info).dataparams).t		&	t_upp = t_low
	ENDIF ELSE BEGIN
		t_low = (*(*info).dispparams).t_low	&	t_upp = (*(*info).dispparams).t_upp 
		lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
	ENDELSE
	WIDGET_CONTROL, /HOURGLASS
	IF STRCMP((*(*info).savparams).savpro,'MPEG') THEN mpegid = MPEG_OPEN([(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny], QUALITY = 100, BITRATE = 1E8)
	FOR i = t_low, t_upp DO BEGIN
		(*(*info).dataparams).t = i
		CRISPEX_UPDATE_T, event
		IF (*(*info).savparams).overlays_incl THEN BEGIN
			CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
				NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
			image=TVRD()
		ENDIF ELSE CRISPEX_DRAW_XYSLICE_SCALING, event, image
		TVLCT,r,g,b,/GET
		IF (*(*info).savparams).snapshot THEN midtension = '' ELSE BEGIN
			t_id = 't'+STRING((*(*info).dataparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
			midtension = '_'+lp_id+'_'+t_id
		ENDELSE
		IF (STRCMP((*(*info).savparams).savpro,'JPEG_FRAMES') OR STRCMP((*(*info).savparams).savpro,'MPEG')) THEN BEGIN
			s = SIZE(image)
			finalimage = BYTARR(3,s[1],s[2])
			finalimage(0,*,*) = r(image)
			finalimage(1,*,*) = g(image)
			finalimage(2,*,*) = b(image)
			IF STRCMP((*(*info).savparams).savpro,'MPEG') THEN BEGIN
				imdisp = REVERSE(finalimage,3)
				MPEG_PUT,mpegid, IMAGE = imdisp, FRAME = i-(*(*info).dispparams).t_low, /COLOR
			ENDIF ELSE BEGIN
				filename = (*(*info).paths).opath+supplied_filename+midtension+'.jpg'
				WRITE_JPEG, filename,finalimage,TRUE=1,QUALITY = 75				
				PRINT, 'Written: '+filename
			ENDELSE
		ENDIF
		IF STRCMP((*(*info).savparams).savpro,'PNG_FRAMES') THEN BEGIN
			filename = (*(*info).paths).opath+supplied_filename+midtension+'.png'
			WRITE_PNG,filename,image,r,g,b
			PRINT, 'Written: '+filename
		ENDIF
	ENDFOR
	IF STRCMP((*(*info).savparams).savpro,'MPEG') THEN BEGIN
		MPEG_SAVE, mpegid, FILENAME = (*(*info).paths).opath+supplied_filename+'.mpg'
		MPEG_CLOSE, mpegid
		PRINT, 'Written: '+(*(*info).paths).opath+supplied_filename+'.mpg'
	ENDIF
	(*(*info).dataparams).t = t_before
	CRISPEX_UPDATE_T, event
	IF (*(*info).savparams).overlays_incl THEN CRISPEX_DRAW_XY, event
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SAVE_LINESCAN, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_LINESCAN'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename
	IF ((*(*info).savparams).savpro EQ 'JPEG_LINESCAN') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save JPEG line scan', standardfilename, 'CRISPEX_SAVE_CHECK'
	IF ((*(*info).savparams).savpro EQ 'PNG_LINESCAN') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save PNG line scan', standardfilename, 'CRISPEX_SAVE_CHECK'
END

PRO CRISPEX_SAVE_LINESCAN_SAVE, event, supplied_filename
; Handles the saving of a series (between temporal boundaries) of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_LINESCAN_SAVE'
	nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
	ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
	lp_before = (*(*info).dataparams).lp
	WIDGET_CONTROL, /HOURGLASS
	t_id = 't'+STRING((*(*info).dataparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
	FOR i = (*(*info).dispparams).lp_low, (*(*info).dispparams).lp_upp DO BEGIN
		(*(*info).dataparams).lp = i
		CRISPEX_UPDATE_T, event
		IF (*(*info).savparams).overlays_incl THEN BEGIN
			CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
				NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
			image=TVRD()
		ENDIF ELSE CRISPEX_DRAW_XYSLICE_SCALING, event, image
		TVLCT,r,g,b,/GET
		IF (*(*info).savparams).linescan_ls THEN BEGIN
			CRISPEX_UPDATE_LP, event
			CRISPEX_DRAW_LS, event
			lsimage = TVRD()
		ENDIF
		lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
		IF ((*(*info).savparams).savpro EQ 'JPEG_LINESCAN') THEN BEGIN
			s = SIZE(image)
			finalimage = BYTARR(3,s[1],s[2])
			finalimage(0,*,*) = r(image)
			finalimage(1,*,*) = g(image)
			finalimage(2,*,*) = b(image)
			filename = (*(*info).paths).opath+supplied_filename+'_'+t_id+'_'+lp_id+'.jpg'
			WRITE_JPEG, filename,finalimage,TRUE=1,QUALITY = 75
			PRINT, 'Written: '+filename
			IF (*(*info).savparams).linescan_ls THEN BEGIN
				lsfilename = (*(*info).paths).opath+supplied_filename+'_detspect_'+t_id+'_'+lp_id+'.jpg'
				WRITE_JPEG, lsfilename,lsimage,QUALITY = 75
				PRINT, 'Written: '+lsfilename
			ENDIF
		ENDIF
		IF ((*(*info).savparams).savpro EQ 'PNG_LINESCAN') THEN BEGIN
			filename = (*(*info).paths).opath+supplied_filename+'_'+t_id+'_'+lp_id+'.png'
			WRITE_PNG,filename,image,r,g,b
			PRINT, 'Written: '+filename
			IF (*(*info).savparams).linescan_ls THEN BEGIN
				lsfilename = (*(*info).paths).opath+supplied_filename+'_detspect_'+t_id+'_'+lp_id+'.png'
				WRITE_PNG,lsfilename,lsimage
				PRINT, 'Written: '+lsfilename
			ENDIF
		ENDIF
	ENDFOR
	(*(*info).dataparams).lp = lp_before
	CRISPEX_UPDATE_T, event
	IF (*(*info).savparams).overlays_incl THEN CRISPEX_DRAW_XY, event
	IF (*(*info).savparams).linescan_ls THEN BEGIN
		CRISPEX_UPDATE_LP, event
		CRISPEX_DRAW_LS, event
	ENDIF
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SAVE_MPEG, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_MPEG'
	(*(*info).savparams).savpro = 'MPEG'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename, import_id='lp'+STRTRIM(LONG((*(*info).dataparams).lp),2)
	CRISPEX_SAVE_GET_FILENAME, event, 'Save MPEG movie', standardfilename, 'CRISPEX_SAVE_CHECK'
END

PRO CRISPEX_SAVE_OPTIONS, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS'
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Saving options', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	path_base = WIDGET_BASE(disp, /COLUMN, /FRAME)
	path_label = WIDGET_LABEL(path_base, VALUE = 'Path options', /ALIGN_LEFT)
	curpath_base = WIDGET_BASE(path_base, /ROW)
	curpath_label = WIDGET_LABEL(curpath_base, VALUE = 'Current path:')
	(*(*info).ctrlssav).savopt_path_textlab = WIDGET_LABEL(curpath_base, VALUE = STRTRIM((*(*info).paths).opath,2), /DYNAMIC_RESIZE)
	pathbut_base = WIDGET_BASE(path_base, COLUMN=1,/GRID_LAYOUT,/ALIGN_LEFT)
	change_path_but = WIDGET_BUTTON(pathbut_base, VALUE = 'Change path', EVENT_PRO = 'CRISPEX_SAVE_SET_OPATH')
	overlays_base = WIDGET_BASE(disp, /COLUMN, /FRAME)
	overlays_label = WIDGET_LABEL(overlays_base, VALUE = 'Overlays options', /ALIGN_LEFT)
	overlays_but_base = WIDGET_BASE(overlays_base, /ROW)
	overlays_buts_base = WIDGET_BASE(overlays_but_base, /COLUMN, /NONEXCLUSIVE)
	overlays_incl_but = WIDGET_BUTTON(overlays_buts_base, VALUE = 'Include overlays', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE')
	overlays_subbuts_base = WIDGET_BASE(overlays_but_base, /COLUMN, /NONEXCLUSIVE)
	WIDGET_CONTROL, overlays_incl_but, SET_BUTTON = (*(*info).savparams).overlays_incl
	(*(*info).ctrlssav).overlays_num_but = WIDGET_BUTTON(overlays_subbuts_base, VALUE = 'Number overlays', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_NUMBER', SENSITIVE = (*(*info).savparams).overlays_incl)
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_num_but, SET_BUTTON = (*(*info).savparams).overlays_num
	(*(*info).ctrlssav).overlays_curs_but = WIDGET_BUTTON(overlays_subbuts_base, VALUE = 'Include cursor', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_CURSOR', SENSITIVE = (*(*info).savparams).overlays_incl)
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_curs_but, SET_BUTTON = (*(*info).savparams).overlays_curs
	(*(*info).ctrlssav).overlays_pts_but = WIDGET_BUTTON(overlays_subbuts_base, VALUE = 'Include endpoints', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ENDPOINTS', SENSITIVE = (*(*info).savparams).overlays_incl)
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_pts_but, SET_BUTTON = (*(*info).savparams).overlays_pts
	(*(*info).ctrlssav).overlays_asecbar_but = WIDGET_BUTTON(overlays_subbuts_base, VALUE = 'Add arcseconds bar', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR', SENSITIVE = (*(*info).savparams).overlays_incl)
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_asecbar_but, SET_BUTTON = (*(*info).savparams).overlays_asecbar
	overlays_combo_base = WIDGET_BASE(overlays_base, /COLUMN)
	overlays_thick_base = WIDGET_BASE(overlays_combo_base, /ROW)
	(*(*info).ctrlssav).overlays_thick_slider = WIDGET_SLIDER(overlays_thick_base, TITLE = 'Overlay thickness', MIN = 1, MAX = 8, VALUE = (*(*info).savparams).overlays_thick, $
		EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_THICK', SENSITIVE = (*(*info).savparams).overlays_incl, /DRAG)
	(*(*info).ctrlssav).overlays_symsize_slider = WIDGET_SLIDER(overlays_thick_base, TITLE = 'Overlay symbol size', MIN = 1, MAX = 8, VALUE = (*(*info).savparams).overlays_symsize, $
		EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_SYMSIZE', SENSITIVE = (*(*info).savparams).overlays_incl, /DRAG)
	(*(*info).ctrlssav).overlays_asecbar_slider = WIDGET_SLIDER(overlays_base, TITLE = 'Arcseconds bar length', MIN = 1, MAX = 8, VALUE = (*(*info).savparams).overlays_asecbar_length, $
		EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR_LENGTH', SENSITIVE = ((*(*info).savparams).overlays_incl AND (*(*info).savparams).overlays_asecbar), /DRAG)
	linescan_base = WIDGET_BASE(disp,/COLUMN, /FRAME)
	linescan_label = WIDGET_LABEL(linescan_base, VALUE = 'Linescan options', /ALIGN_LEFT)
	linescan_but_base = WIDGET_BASE(linescan_base, /NONEXCLUSIVE)
	linescan_incl_ls_but = WIDGET_BUTTON(linescan_but_base, VALUE = 'Save detailed spectrum', EVENT_PRO = 'CRISPEX_SAVE_OPTIONS_INCLUDE_DETSPECT', $
		SENSITIVE = (((*(*info).dataparams).nlp GT 1) AND (((*(*info).savparams).savpro EQ 'JPEG_LINESCAN') OR ((*(*info).savparams).savpro EQ 'PNG_LINESCAN'))))
	WIDGET_CONTROL, linescan_incl_ls_but, SET_BUTTON = (*(*info).savparams).linescan_ls
	button_base = WIDGET_BASE(disp,COLUMN=1,/GRID_LAYOUT,/ALIGN_CENTER)
	ok_but = WIDGET_BUTTON(button_base, VALUE = 'OK', EVENT_PRO = 'CRISPEX_CLOSE_EVENT_WINDOW')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).spxoffset, TLB_SET_YOFFSET = 0, /TLB_KILL_REQUEST_EVENTS
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).saveoptwintlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).saveoptwintlb],labels=['saveoptwintlb']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE'
	(*(*info).savparams).overlays_incl = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_num_but, SENSITIVE = (*(*info).savparams).overlays_incl
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_curs_but, SENSITIVE = (*(*info).savparams).overlays_incl
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_pts_but, SENSITIVE = (*(*info).savparams).overlays_incl
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_asecbar_but, SENSITIVE = (*(*info).savparams).overlays_incl
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_thick_slider, SENSITIVE = (*(*info).savparams).overlays_incl
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_symsize_slider, SENSITIVE = (*(*info).savparams).overlays_incl
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_incl],labels=['Include overlays']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY'
	IF (*(*info).savparams).overlays_incl THEN BEGIN
		CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
			NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
	ENDIF ELSE CRISPEX_DRAW_XY, event
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_NUMBER, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_NUMBER'
	(*(*info).savparams).overlays_num = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_num],labels=['Number overlays']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_CURSOR, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_CURSOR'
	(*(*info).savparams).overlays_curs = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_curs],labels=['Include cursor']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_THICK, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_THICK'
	(*(*info).savparams).overlays_thick = event.VALUE
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_thick],labels=['Overlay thickness']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_SYMSIZE, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_SYMSIZE'
	(*(*info).savparams).overlays_symsize = event.VALUE
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_symsize],labels=['Overlay symbol size']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ENDPOINTS, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ENDPOINTS'
	(*(*info).savparams).overlays_pts = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_pts],labels=['Include endpoints']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR'
	(*(*info).savparams).overlays_asecbar = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_asecbar_slider, SENSITIVE = (*(*info).savparams).overlays_asecbar
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_asecbar],labels=['Add arcseconds bar']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR_LENGTH, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR_LENGTH'
	(*(*info).savparams).overlays_asecbar_length = event.VALUE
	(*(*info).savparams).overlays_asecbar_pix = (*(*info).savparams).overlays_asecbar_length / FLOAT((*(*info).meas).arcsecpix) * (*(*info).zooming).factor
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_asecbar_length, (*(*info).savparams).overlays_asecbar_pix],$
		labels=['Arcseconds bar length','Arcseconds bar length in pixels']
END

PRO CRISPEX_SAVE_OPTIONS_INCLUDE_DETSPECT, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SAVE_OPTIONS_INCLUDE_DETSPECT'
	(*(*info).savparams).linescan_ls = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).linescan_ls],labels=['Save detailed spectrum']
END

;================================================================================= SLIDER CONTROL PROCEDURES
PRO CRISPEX_SLIDER_NPHI, event
; Handles the change in spectral slit length slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_NPHI'
	(*(*info).phiparams).nphi_set = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).phiparams).nphi_set],labels=['Phi-slit length']
	CRISPEX_UPDATE_PHIS, event
END

PRO CRISPEX_SLIDER_PHI, event
; Handles the change in spectral slit angle slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_PHI'
	(*(*info).phiparams).angle = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).phiparams).angle],labels=['Phi-slit angle']
	CRISPEX_PHISLIT_DIRECTION, event
	CRISPEX_UPDATE_PHIS, event
END

PRO CRISPEX_SLIDER_LP, event
; Handles the change in spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP'
	(*(*info).dataparams).lp = event.VALUE
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_DECR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP_DECR'
	(*(*info).dataparams).lp -= (*(*info).pbparams).lp_step 
	IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN (*(*info).dataparams).lp = (*(*info).dispparams).lp_upp
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_INCR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP_INCR'
	(*(*info).dataparams).lp = (((*(*info).dataparams).lp - (*(*info).dispparams).lp_low + (*(*info).pbparams).lp_step) MOD ((*(*info).dispparams).lp_upp - (*(*info).dispparams).lp_low + 1)) + (*(*info).dispparams).lp_low
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_REF, event
; Handles change in reference spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP_REF'
	(*(*info).dataparams).lp_ref = event.VALUE
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_REF_LOCK, event
; Handles (un)locking the reference reference to (from) the main spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP_REF_LOCK'
	(*(*info).ctrlsswitch).lp_ref_lock = event.SELECT
	(*(*info).dataparams).lp_ref = (*(*info).dataparams).lp
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SENSITIVE = ABS((*(*info).ctrlsswitch).lp_ref_lock-1), SET_VALUE = (*(*info).dataparams).lp_ref
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).lp, (*(*info).dataparams).lp_ref], labels=['lp','lp_ref']
	CRISPEX_UPDATE_T, event
	CRISPEX_UPDATE_LP, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_LP_UPDATE, event
; Handles the the update after change in the spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_LP_UPDATE'
	IF (*(*info).ctrlsswitch).lp_ref_lock THEN (*(*info).dataparams).lp_ref = (*(*info).dataparams).lp
	IF (*(*info).ctrlsswitch).lp_ref_lock THEN WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SET_VALUE = (*(*info).dataparams).lp_ref
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).lp, (*(*info).dataparams).lp_ref], labels=['lp','lp_ref']
	CRISPEX_UPDATE_T, event
	CRISPEX_UPDATE_LP, event
	IF (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling] THEN BEGIN
		IF ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 2) THEN CRISPEX_SCALING_MAN_FIRST, event
		IF ((*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] EQ 3) THEN CRISPEX_SCALING_MAN_CURR, event
	ENDIF 
	CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_SPECTSTEP, event
; Handles the change in spectral step (for spectral blink) slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_SPECTSTEP'
	(*(*info).pbparams).lp_step = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).lp_step], labels=['lp_step']
END

PRO CRISPEX_SLIDER_SPEED, event
; Handles the change in playback speed slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_SPEED', /IGNORE_LAST
	(*(*info).pbparams).t_speed = event.VALUE
	IF (*(*info).pbparams).spmode THEN WIDGET_CONTROL, (*(*info).ctrlscp).t_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed ELSE $
		WIDGET_CONTROL, (*(*info).ctrlscp).lp_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).t_speed], labels=['Playback (blink) speed']
END

PRO CRISPEX_SLIDER_STEP, event
; Handles the change in temporal step slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_STEP'
	(*(*info).pbparams).t_step = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).t_step], labels=['t_step']
END

PRO CRISPEX_SLIDER_T, event
; Handles the change in temporal slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_T'
	(*(*info).dataparams).t = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).t], labels=['t']
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW, event
	IF (*(*info).winswitch).showphis THEN BEGIN 
		IF (*(*info).dispparams).phislice_update THEN CRISPEX_UPDATE_SLICES, event ELSE BEGIN		
			IF (*(*info).dataswitch).onecube THEN WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1, SET_VALUE = 'Update spectral windows' ELSE WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1
		ENDELSE
	ENDIF
END

PRO CRISPEX_SLIDER_X, event
; Handles the change in cursor x-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_X'
	(*(*info).dataparams).x = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x], labels=['x']
	CRISPEX_UPDATE_SX, event
	IF (*(*info).winswitch).showphis THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
		CRISPEX_UPDATE_PHIS, event
	ENDIF ELSE CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_XPOS, event
; Handles change in x-slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_XPOS'
	(*(*info).zooming).xpos = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).xpos], labels=['xpos']
	CRISPEX_UPDATE_SX, event
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_Y, event
; Handles the change in cursor y-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_Y'
	(*(*info).dataparams).y = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).y], labels=['y']
	CRISPEX_UPDATE_SY, event
	IF (*(*info).winswitch).showphis THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
		CRISPEX_UPDATE_PHIS, event
	ENDIF ELSE CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_YPOS, event
; Handles change in y-slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_SLIDER_YPOS'
	(*(*info).zooming).ypos = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).ypos], labels=['ypos']
	CRISPEX_UPDATE_SY, event
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW, event
END

;================================================================================= UPDATE SLICES AND PARAMETERS PROCEDURES
PRO CRISPEX_UPDATE_SLICES, event, NO_DRAW=no_draw
; Gets the new spectral phi slit scan for update of the spectral phi slit slice after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_SLICES'
	IF ((*(*info).winswitch).showphis OR ((*(*info).dataswitch).onecube AND (*(*info).winswitch).showls)) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		IF ((*(*info).dataparams).nt GT 1) THEN *(*(*info).data).sspscan = (*(*(*info).data).scan)[(*(*info).dataparams).t] ELSE *(*(*info).data).sspscan = (*(*(*info).data).scan)
		*(*(*info).data).phiscan = (*(*(*info).data).sspscan)[*,*,((*(*info).dataparams).s * (*(*info).dataparams).nlp):(((*(*info).dataparams).s+1)*(*(*info).dataparams).nlp-1)] 
		CRISPEX_UPDATE_PHIS, event, NO_DRAW=no_draw
		WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 0
	ENDIF
END

PRO CRISPEX_UPDATE_PHIS, event, NO_DRAW=no_draw
; Handles the actual update of the spectral phi slit slice after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_PHIS'
	CRISPEX_PHISLIT_GET_SLICE, event
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_UPDATE_LP, event
; Handles the update of displayed data after change in spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_LP'
	IF (*(*info).overlayswitch).loopslit AND ((*(*info).loopparams).np GT 0) THEN BEGIN
		IF (*(*info).dispswitch).exts THEN *(*(*info).loopsdata).loopslice = REFORM((*(*(*info).loopsdata).loopslab)[*,*,(*(*info).dataparams).lp-(*(*info).dispparams).lp_low]) ELSE $
			*(*(*info).loopsdata).loopslice = REFORM((*(*(*info).loopsdata).loopslab)[*,*,(*(*info).dataparams).lp])
		IF (*(*info).winswitch).showrefloop THEN BEGIN
			IF (*(*info).dispswitch).refexts THEN *(*(*info).loopsdata).refloopslice = REFORM((*(*(*info).loopsdata).refloopslab)[*,*,(*(*info).dataparams).lp_ref-(*(*info).dispparams).lp_ref_low]) ELSE $
				*(*(*info).loopsdata).refloopslice = REFORM((*(*(*info).loopsdata).refloopslab)[*,*,(*(*info).dataparams).lp_ref])
		ENDIF
	ENDIF
	IF (*(*info).winswitch).showrestloop THEN FOR k=0,N_ELEMENTS(*(*(*info).restoreparams).disp_loopnr)-1 DO BEGIN
		IF (SIZE(*(*(*(*info).loopsdata).rest_loopslab[k]),/N_DIMENSIONS) EQ 3) THEN BEGIN
			IF (*(*(*info).restoreparams).disp_imref)[k] THEN *(*(*(*info).loopsdata).rest_loopslice[k]) = REFORM((*(*(*(*info).loopsdata).rest_loopslab[k]))[*,*,(*(*info).dataparams).lp_ref-$
				(*(*info).dispparams).lp_ref_low]) ELSE *(*(*(*info).loopsdata).rest_loopslice[k]) = REFORM((*(*(*(*info).loopsdata).rest_loopslab[k]))[*,*,(*(*info).dataparams).lp-(*(*info).dispparams).lp_low])
		ENDIF ELSE BEGIN
			IF (*(*(*info).restoreparams).disp_imref)[k] THEN *(*(*(*info).loopsdata).rest_loopslice[k]) = *(*(*(*info).loopsdata).rest_loopslab[k]) ELSE $
				*(*(*(*info).loopsdata).rest_loopslice[k]) = *(*(*(*info).loopsdata).rest_loopslab[k])
		ENDELSE
	ENDFOR
	IF (*(*info).winswitch).showretrdet THEN *(*(*info).loopsdata).det_loopslice = REFORM((*(*(*info).loopsdata).det_loopslab)[*,*,(*(*info).dataparams).lp-(*(*info).dispparams).lp_low])
END

PRO CRISPEX_UPDATE_T, event
; Handles the updated of displayed data after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_T'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		x_low = (*(*info).zooming).xpos
		x_upp = (*(*info).zooming).xpos + (*(*info).dataparams).d_nx
		y_low = (*(*info).zooming).ypos
		y_upp = (*(*info).zooming).ypos + (*(*info).dataparams).d_ny
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [x_low,x_upp,(*(*info).dataparams).d_nx,(*(*info).dataparams).nx,y_low,y_upp,(*(*info).dataparams).d_ny,(*(*info).dataparams).ny],$
			labels=['x_low','x_upp','d_nx','nx','y_low','y_upp','d_ny','ny']
	ENDIF
	IF (*(*info).winswitch).showdop THEN BEGIN
		(*(*info).dataparams).lp_dop = 2*(*(*info).dataparams).lc - (*(*info).dataparams).lp
		(*(*info).dispswitch).drawdop = (((*(*info).dataparams).lp_dop GE (*(*info).dispparams).lp_low) AND ((*(*info).dataparams).lp_dop LE (*(*info).dispparams).lp_upp) AND $
			((*(*info).dataparams).lp_dop NE (*(*info).dataparams).lc)) 
	ENDIF
	IF ((*(*info).dataswitch).spfile EQ 1) OR (*(*info).dataswitch).onecube THEN BEGIN
		IF((*(*info).zooming).factor NE 1) THEN BEGIN
			*(*(*info).data).xyslice = CONGRID( REFORM( ((*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
			(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) 
			IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN $
				temp_xyslice = CONGRID( REFORM( ((*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
				(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp_dop])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny)
		ENDIF ELSE BEGIN
			*(*(*info).data).xyslice = REFORM((*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + (*(*info).dataparams).s * (*(*info).dataparams).nlp + $
			(*(*info).dataparams).lp])
			IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN temp_xyslice = REFORM((*(*(*info).data).imagedata)[(*(*info).dataparams).t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
				(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp_dop])
		ENDELSE
	ENDIF ELSE BEGIN
		IF ((*(*info).zooming).factor NE 1) THEN BEGIN
			*(*(*info).data).xyslice = CONGRID( REFORM( ((*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + $
			(*(*info).dataparams).lp])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx,(*(*info).winsizes).xywiny ) 
			IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN temp_xyslice = CONGRID( REFORM( ((*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + $
				(*(*info).dataparams).lp_dop])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx,(*(*info).winsizes).xywiny )
		ENDIF ELSE BEGIN
			*(*(*info).data).xyslice = REFORM((*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp] ) 
			IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN $
				temp_xyslice = REFORM((*(*(*info).data).imagedata)[(*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp_dop] ) 
		ENDELSE
	ENDELSE
	IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN BEGIN
		IF ((*(*info).dataparams).lp_dop GT (*(*info).dataparams).lc) THEN *(*(*info).data).dopslice = temp_xyslice - *(*(*info).data).xyslice ELSE *(*(*info).data).dopslice = *(*(*info).data).xyslice - temp_xyslice 
	ENDIF
	IF ((*(*info).winswitch).showref OR (*(*info).winswitch).showimref) THEN BEGIN
		IF (*(*info).dataswitch).refspfile THEN BEGIN
			IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( REFORM( ((*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + $
				(*(*info).dataparams).lp_ref])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE $
				*(*(*info).data).refslice = REFORM( (*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref])
		ENDIF ELSE BEGIN
			IF ((*(*info).dataparams).refnt EQ 0) THEN BEGIN
				IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( (*(*(*info).data).refdata)[x_low:x_upp, y_low:y_upp], (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) $
				ELSE *(*(*info).data).refslice = *(*(*info).data).refdata 
			ENDIF ELSE IF ((*(*info).dataparams).refnt EQ 1) THEN BEGIN
				IF ((*(*info).dataparams).refnlp NE 1) THEN BEGIN
					IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( REFORM( ((*(*(*info).data).refdata)[(*(*info).dataparams).lp_ref])[x_low:x_upp, y_low:y_upp]), $
						(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE *(*(*info).data).refslice = REFORM( (*(*(*info).data).refdata)[(*(*info).dataparams).lp_ref])
				ENDIF ELSE BEGIN
					IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( REFORM( ((*(*(*info).data).refdata)[0])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, $
						(*(*info).winsizes).xywiny) ELSE *(*(*info).data).refslice = REFORM( (*(*(*info).data).refdata)[0]) 
				ENDELSE
			ENDIF ELSE IF ((*(*info).dataparams).refnt EQ (*(*info).dataparams).nt) THEN BEGIN
				IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( REFORM( ((*(*(*info).data).refdata)[(*(*info).dataparams).t])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, $
				(*(*info).winsizes).xywiny) ELSE *(*(*info).data).refslice = REFORM( (*(*(*info).data).refdata)[(*(*info).dataparams).t])
			ENDIF ELSE IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).refslice = CONGRID( REFORM( ((*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + $
				(*(*info).dataparams).lp_ref])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) ELSE $
				*(*(*info).data).refslice = REFORM( (*(*(*info).data).refdata)[(*(*info).dataparams).t * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref])
		ENDELSE
	ENDIF
	IF (*(*info).dataswitch).maskfile THEN BEGIN
		IF ((*(*info).dataparams).masknt GT 1) THEN BEGIN
			IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).maskslice = CONGRID( REFORM( ((*(*(*info).data).maskdata)[(*(*info).dataparams).t])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, $
				(*(*info).winsizes).xywiny) ELSE *(*(*info).data).maskslice = REFORM( (*(*(*info).data).maskdata)[(*(*info).dataparams).t])
		ENDIF ELSE BEGIN
			IF ((*(*info).zooming).factor NE 1) THEN *(*(*info).data).maskslice = CONGRID( REFORM( ((*(*(*info).data).maskdata)[0])[x_low:x_upp, y_low:y_upp]), (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) $
				ELSE *(*(*info).data).maskslice = (*(*(*info).data).maskdata)[0]
		ENDELSE
	ENDIF
END

PRO CRISPEX_UPDATE_SX, event
; Handles the change in xy- and reference image x-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_SX'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		sx = ((*(*info).dataparams).x - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		sxp = (*(*(*info).loopparams).xp - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
	ENDIF ELSE BEGIN
		sx = (*(*info).dataparams).x * (*(*info).winsizes).xywinx / FLOAT((*(*info).dataparams).nx)
		sxp = *(*(*info).loopparams).xp * (*(*info).winsizes).xywinx / FLOAT((*(*info).dataparams).nx)
		sxr = *(*(*info).loopparams).xr * (*(*info).winsizes).xywinx / FLOAT((*(*info).dataparams).nx)
	ENDELSE
	(*(*info).curs).sxlock = sx 
	(*(*info).curs).sx = (*(*info).curs).sxlock
	*(*(*info).overlayparams).sxp = sxp 
	*(*(*info).overlayparams).sxr = sxr 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x, (*(*info).curs).sxlock, (*(*info).curs).sx],labels=['x','sxlock','sx']
END

PRO CRISPEX_UPDATE_SY, event
; Handles the change in xy- and reference image y-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_SY'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		sy = ((*(*info).dataparams).y - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
		syp = (*(*(*info).loopparams).yp - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
		syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
	ENDIF ELSE BEGIN
		sy = (*(*info).dataparams).y * (*(*info).winsizes).xywiny / FLOAT((*(*info).dataparams).ny)
		syp = *(*(*info).loopparams).yp * (*(*info).winsizes).xywiny / FLOAT((*(*info).dataparams).ny)
		syr = *(*(*info).loopparams).yr * (*(*info).winsizes).xywiny / FLOAT((*(*info).dataparams).ny)
	ENDELSE
	(*(*info).curs).sylock = sy 
	(*(*info).curs).sy = (*(*info).curs).sylock
	*(*(*info).overlayparams).syp = syp 
	*(*(*info).overlayparams).syr = syr 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).y, (*(*info).curs).sylock, (*(*info).curs).sy],labels=['y','sylock','sy']
END

PRO CRISPEX_UPDATE_USER_FEEDBACK, event, title=title, var=var, minvar=minvar, maxvar=maxvar, feedback_text=feedback_text, destroy_top=destroy_top, close_button=close_button, session=session
; Handles the update of user feedback while saving timeslices
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_UPDATE_USER_FEEDBACK'
	IF ~KEYWORD_SET(MINVAR) THEN minvar = 0
	IF (var EQ minvar) THEN BEGIN
		CRISPEX_WINDOW_USER_FEEDBACK, event, title, feedback_text+'     ', close_button=close_button, session=session
		IF (N_ELEMENTS(destroy_top) EQ 1) THEN destroy_top = destroy_top ELSE destroy_top = 0
		IF destroy_top THEN BEGIN
			WIDGET_CONTROL, (*(*info).winids).feedbacktlb, SET_UVALUE = info
			WIDGET_CONTROL, event.TOP, /DESTROY
			event.TOP = (*(*info).winids).feedbacktlb
		ENDIF
	ENDIF ELSE BEGIN
		WIDGET_CONTROL,(*(*info).ctrlsfeedb).feedback_text, SET_VALUE = feedback_text
		IF KEYWORD_SET(CLOSE_BUTTON) THEN WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, SENSITIVE = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [feedback_text], labels=['Feedback text']
END

PRO CRISPEX_UPDATE_STARTUP_FEEDBACK, bgim, xout, yout, feedback_text
	LOADCT,3,/SILENT
	TVSCL,bgim
	LOADCT,0,/SILENT
	FOR i=0,N_ELEMENTS(feedback_text)-1 DO XYOUTS, xout[i], yout[i], feedback_text[i], COLOR=255, /DEVICE, CHARSIZE=1.125
END

;================================================================================= VERBOSE PROCEDURES
PRO CRISPEX_VERBOSE_GET, event, variables, labels=labels
; Displays feedback of variables structuredly 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	nvars = N_ELEMENTS(variables)
	IF ((*(*info).feedbparams).last_routine_count NE 0) THEN PRINT,''
	IF (N_ELEMENTS(labels) EQ nvars) THEN BEGIN
		maxchar = MAX(STRLEN(labels))
		FOR i=0,nvars-1 DO BEGIN
			IF (STRLEN(labels[i]) NE maxchar) THEN whitespace = STRJOIN(REPLICATE(' ',(maxchar-STRLEN(labels[i])))) ELSE whitespace = ''
			PRINT,STRJOIN(REPLICATE('  ',SCOPE_LEVEL()-2))+'> '+labels[i]+whitespace+': '+STRTRIM(variables[i],2)
		ENDFOR
	ENDIF ELSE BEGIN
		FOR i=0,nvars-1 DO HELP, variables[i]
	ENDELSE
END

PRO CRISPEX_VERBOSE_GET_ROUTINE, event, rname, IGNORE_LAST=ignore_last
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	prespace = STRJOIN(REPLICATE('  ',SCOPE_LEVEL()-2))
	IF KEYWORD_SET(IGNORE_LAST) THEN (*(*info).feedbparams).last_routine = ''
	IF ((rname NE (*(*info).feedbparams).last_routine) AND ((*(*info).feedbparams).last_routine_count GT 0)) THEN PRINT,''
	IF (rname EQ (*(*info).feedbparams).last_routine) THEN (*(*info).feedbparams).last_routine_count += 1 ELSE (*(*info).feedbparams).last_routine_count = 0
	IF ((*(*info).feedbparams).last_routine_count GT 0) THEN rcount = ' x '+STRTRIM((*(*info).feedbparams).last_routine_count,2)+'.' ELSE rcount = '.'
	IF (rname NE (*(*info).feedbparams).last_routine) THEN PRINT,prespace+'CRISPEX RUN: Called '+rname+'.' ELSE $
		WRITEU,-1,STRING(FORMAT='(%"\r'+prespace+'CRISPEX RUN: Called ",a'+STRTRIM(STRLEN(rname),2)+',a'+STRTRIM(STRLEN(rcount),2)+')',rname,rcount) 
	(*(*info).feedbparams).last_routine = rname
END

PRO CRISPEX_VERBOSE_SET, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_VERBOSE_SET'
	WIDGET_CONTROL, event.ID, GET_UVALUE = add_verbosity
	tmp_verbosity = (*(*info).feedbparams).verbosity
	IF (add_verbosity EQ -1) THEN (*(*info).feedbparams).verbosity = INTARR(N_ELEMENTS((*(*info).feedbparams).verbosity)) ELSE BEGIN
		tmp_verbosity[add_verbosity] = ABS((tmp_verbosity[add_verbosity] EQ 1) - 1)
		(*(*info).feedbparams).verbosity = tmp_verbosity
	ENDELSE
	CRISPEX_VERBOSE_SET_BUTTONS, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [add_verbosity], labels=['Verbosity level']
END

PRO CRISPEX_VERBOSE_SET_BUTTONS, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_VERBOSE_SET'
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[0], SET_BUTTON = (TOTAL((*(*info).feedbparams).verbosity) EQ 0)
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[1], SET_BUTTON = ((*(*info).feedbparams).verbosity)[2]
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[2], SET_BUTTON = ((*(*info).feedbparams).verbosity)[3]
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[3], SET_BUTTON = ((*(*info).feedbparams).verbosity)[4]
END

;================================================================================= GENERAL WINDOW PROCEDURES
PRO CRISPEX_WINDOW, xsize, ysize, leader, title, base, wid, xoffset, yoffset, DRAWID = drawid, DRAWBASE =disp, SCROLL = xscrollsize, YSCROLL = yscrollsize, SCROLL = scroll, RESIZING = resizing, $
	RES_HANDLER = res_handler
; Sets up the display windows
	IF (N_ELEMENTS(RESIZING) EQ 0) THEN resizing = 0
	IF (N_ELEMENTS(LEADER) EQ 0) THEN base = WIDGET_BASE(TITLE = title, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS, TLB_SIZE_EVENTS = resizing) ELSE $
		base = WIDGET_BASE(TITLE = STRTRIM(title), GROUP_LEADER = leader, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS, TLB_SIZE_EVENTS = resizing)
	disp = WIDGET_BASE(base, /COLUMN)
	drawid = WIDGET_DRAW(disp, XSIZE = xsize, YSIZE = ysize, RETAIN = 2)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET=xoffset, TLB_SET_YOFFSET=yoffset
	IF (N_ELEMENTS(RES_HANDLER) GT 0) THEN XMANAGER, 'CRISPEX', base, EVENT_HANDLER = res_handler, /NO_BLOCK
	WIDGET_CONTROL, drawid, GET_VALUE = wid
END

PRO CRISPEX_WINDOW_OK, event, title, message1, message2, message3, message4, message5, OK_EVENT=ok_event, CANCEL_EVENT=cancel_event, CANCEL_LABEL=cancel_label, BASE=base, BLOCK=block
; Sets up the message windows with only an OK-button
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_WINDOW_OK'
	fulltitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+STRTRIM(title,2)
	base = WIDGET_BASE(TITLE = fulltitle, GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	text_label1 = WIDGET_LABEL(message_base, VALUE = message1)
	IF (N_ELEMENTS(message2) GT 0) THEN text_label2 = WIDGET_LABEL(message_base, VALUE = message2)
	IF (N_ELEMENTS(message3) GT 0) THEN text_label3 = WIDGET_LABEL(message_base, VALUE = message3)
	IF (N_ELEMENTS(message4) GT 0) THEN text_label4 = WIDGET_LABEL(message_base, VALUE = message4)
	IF (N_ELEMENTS(message5) GT 0) THEN text_label5 = WIDGET_LABEL(message_base, VALUE = message5)
	IF (N_ELEMENTS(CANCEL_EVENT) GT 0) THEN BEGIN
		IF (N_ELEMENTS(CANCEL_LABEL) NE 1) THEN cancel_label = 'Cancel'
		button_base = WIDGET_BASE(disp,COLUMN=2,/GRID_LAYOUT,/ALIGN_CENTER) 
		cancel_but = WIDGET_BUTTON(button_base, VALUE = cancel_label, EVENT_PRO = cancel_event)
	ENDIF ELSE button_base = WIDGET_BASE(disp,/ROW,/ALIGN_CENTER)
	ok_but = WIDGET_BUTTON(button_base, VALUE = 'OK' , EVENT_PRO = ok_event)
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	IF (N_ELEMENTS(BLOCK) NE 1) THEN block = 0
	XMANAGER, 'CRISPEX', base, NO_BLOCK=ABS(block-1)
END

PRO CRISPEX_WINDOW_USER_FEEDBACK, event, title, initial_feedback, close_button=close_button, session=session
; Sets up the message windows for user feedback 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_WINDOW_USER_FEEDBACK'
	fulltitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+STRTRIM(title,2)
	base = WIDGET_BASE(TITLE = fulltitle, GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
	message_base = WIDGET_BASE(disp, /COLUMN)
	(*(*info).ctrlsfeedb).feedback_text = WIDGET_LABEL(message_base, VALUE = initial_feedback, /ALIGN_LEFT, /DYNAMIC_RESIZE)
	IF KEYWORD_SET(close_button) THEN BEGIN
		button_base = WIDGET_BASE(disp,/ROW,/ALIGN_CENTER)
		(*(*info).ctrlsfeedb).close_button = WIDGET_BUTTON(button_base, VALUE = 'Close', EVENT_PRO = 'CRISPEX_WINDOW_USER_FEEDBACK_CLOSE', SENSITIVE = 0)
	ENDIF ELSE empty_text = WIDGET_LABEL(message_base, VALUE = ' ')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = 500, TLB_SET_YOFFSET = 500
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	IF KEYWORD_SET(SESSION) THEN (*(*info).winids).restsesfeedbtlb = base ELSE (*(*info).winids).feedbacktlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [base], labels=['feedbacktlb']
END

PRO CRISPEX_WINDOW_USER_FEEDBACK_CLOSE, event, session=session
; Handles the closure of the user feedback window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_WINDOW_USER_FEEDBACK_CLOSE'
	IF KEYWORD_SET(SESSION) THEN BEGIN
		tlb = (*(*info).winids).restsesfeedbtlb 
		WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /DESTROY
		(*(*info).winids).restsesfeedbtlb = 0
	ENDIF ELSE BEGIN
		tlb = (*(*info).winids).feedbacktlb
		WIDGET_CONTROL, (*(*info).winids).feedbacktlb, /DESTROY
		(*(*info).winids).feedbacktlb = 0
	ENDELSE
	(*(*info).ctrlsfeedb).feedback_text = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [tlb], labels=['feedbacktlb was']
END

;================================================================================= ZOOM PROCEDURES
PRO CRISPEX_ZOOM, event, NO_DRAW=no_draw
; Handles the zoom event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOM'
	CRISPEX_UPDATE_SX, event
	CRISPEX_UPDATE_SY, event
	IF (*(*info).overlayswitch).loopslit THEN CRISPEX_ZOOM_LOOP, event
	IF ((*(*info).meas).np GE 1) THEN CRISPEX_ZOOM_MEAS, event
	xposconstr 	= ((*(*info).dataparams).nx-1) - (*(*info).dataparams).d_nx
	yposconstr	= ((*(*info).dataparams).ny-1) - (*(*info).dataparams).d_ny
	WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = xposconstr, SET_VALUE = (*(*info).zooming).xpos 
	WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = yposconstr, SET_VALUE = (*(*info).zooming).ypos 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).factor,xposconstr,yposconstr], labels=['Zoomfactor','Maximum xpos','Maximum ypos']
	CRISPEX_UPDATE_T, event
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
; Handles cursor setting as a result of zoom
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOM_CURSORPOS'
	IF (*(*info).curs).lockset THEN BEGIN
		cursor_x = (*(*info).curs).xlock
		cursor_y = (*(*info).curs).ylock
	ENDIF ELSE BEGIN
		cursor_x = (*(*info).dataparams).x
		cursor_y = (*(*info).dataparams).y
	END
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [cursor_x, cursor_y], labels=['cursor_x','cursor_y']
END	

PRO CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
; Handles the update of xpos and ypos sliders when changing zoomfactor
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOM_UPDATE_SLIDERS'
	(*(*info).dataparams).d_nx = (*(*info).dataparams).nx / (*(*info).zooming).factor
	(*(*info).dataparams).d_ny = (*(*info).dataparams).ny / (*(*info).zooming).factor
	(*(*info).phiparams).d_nphi_set = (*(*info).phiparams).nphi_set / (*(*info).zooming).factor
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).d_nx, (*(*info).dataparams).d_ny], labels=['d_nx','d_ny']
	(*(*info).zooming).xpos = (cursor_x - (*(*info).dataparams).d_nx / 2.) 
	(*(*info).zooming).ypos = (cursor_y - (*(*info).dataparams).d_ny / 2.) 
	(*(*info).zooming).xpos = (*(*info).zooming).xpos > 0
	(*(*info).zooming).ypos = (*(*info).zooming).ypos > 0
	IF (((*(*info).zooming).xpos+(*(*info).dataparams).d_nx) GE (*(*info).dataparams).nx) THEN (*(*info).zooming).xpos = ((*(*info).dataparams).nx-1) - (*(*info).dataparams).d_nx
	IF (((*(*info).zooming).ypos+(*(*info).dataparams).d_ny) GE (*(*info).dataparams).ny) THEN (*(*info).zooming).ypos = ((*(*info).dataparams).ny-1) - (*(*info).dataparams).d_ny
	(*(*info).zooming).xpos = FIX((*(*info).zooming).xpos)
	(*(*info).zooming).ypos = FIX((*(*info).zooming).ypos)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).xpos, (*(*info).zooming).ypos], labels=['xpos','ypos']
	WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SENSITIVE = 1
	WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SENSITIVE = 1

END


PRO CRISPEX_ZOOMFAC_ONE, event
; Sets the zoomfactor to 1 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_ONE'
	(*(*info).zooming).factor = 1	&	(*(*info).zooming).xpos = 0	&	(*(*info).zooming).ypos = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SENSITIVE = 0, SET_VALUE = (*(*info).zooming).xpos
	WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SENSITIVE = 0, SET_VALUE = (*(*info).zooming).ypos
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_one, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOMFAC_TWO, event
; Sets the zoomfactor to 2 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_TWO'
	CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
	(*(*info).zooming).factor = 2.
	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_two, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOMFAC_THREE, event
; Sets the zoomfactor to 3 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_THREE'
	CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
	(*(*info).zooming).factor = 3.
	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_three, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOMFAC_FOUR, event
; Sets the zoomfactor to 4 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_FOUR'
	CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
	(*(*info).zooming).factor = 4.
	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_four, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOMFAC_SIX, event
; Sets the zoomfactor to 6 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_SIX'
	CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
	(*(*info).zooming).factor = 6.
	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_six, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOMFAC_EIGHT, event
; Sets the zoomfactor to 8 and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOMFAC_EIGHT'
	CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
	(*(*info).zooming).factor = 8.
	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y
	WIDGET_CONTROL, (*(*info).ctrlscp).zoom_eight, /SET_BUTTON
	CRISPEX_ZOOM, event
END

PRO CRISPEX_ZOOM_MEAS, event
; Handles the change in measurement position after a zoom event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOM_MEAS'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		(*(*(*info).meas).sxp) = ((*(*(*info).meas).xp) - (*(*info).zooming).xpos ) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		(*(*(*info).meas).syp) = ((*(*(*info).meas).yp) - (*(*info).zooming).ypos ) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
	ENDIF ELSE BEGIN
		(*(*(*info).meas).sxp) = (*(*(*info).meas).xp) * (*(*info).winsizes).xywinx / FLOAT((*(*info).dataparams).nx)
		(*(*(*info).meas).syp) = (*(*(*info).meas).yp) * (*(*info).winsizes).xywiny / FLOAT((*(*info).dataparams).ny)
	ENDELSE
END	

PRO CRISPEX_ZOOM_LOOP, event
; Handles the change in loop positions after a zoom event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event, 'CRISPEX_ZOOM_LOOP'
	IF ((*(*info).zooming).factor NE 1) THEN BEGIN
		*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
	ENDIF ELSE BEGIN
		*(*(*info).overlayparams).sxr = *(*(*info).loopparams).xr * (*(*info).winsizes).xywinx / (*(*info).dataparams).nx
		*(*(*info).overlayparams).syr = *(*(*info).loopparams).yr * (*(*info).winsizes).xywiny / (*(*info).dataparams).ny
	ENDELSE
END

;===================================================================================================
;================================== MAIN PROGRAM CODE ==============================================
;===================================================================================================
PRO CRISPEX, imcube,$										; call program / filename of image cube
	spcube, REFCUBE=refcube, SPECTFILE=spectfile, LINE_CENTER=line_center, $		; spectral & reference cube filename, spectral save file, line centre and/or wavelength information
	DT=dt, EXTS=exts, MNSPEC=mnspec, SINGLE_CUBE=single_cube, $				; time step in seconds, exact timeslices keyword, mean spectrum over selected scans, single full cube call
	SCALE_STOKES=scale_stokes, VALS_REF=vals_ref, VALS_IMG=vals_img, $			; Stokes cube call, scale Stokes spectra internally, Doppler map reference values
	NO_WARP=no_warp, SCALE_CUBES=scale_cubes, XTITLE=xtitle, YTITLE=ytitle, $		; Impede warping of temporal spectrum, custom detailed spectrum xtitle and ytitle
	MASKCUBE=maskcube, WINDOW_LARGE=window_large, VERBOSE=verbose				; Mask cube filename, large window for small cubes, program verbosity

;================================================================================= PROGRAM-INFO ON CALL W/O PARAMS
	IF N_PARAMS() LT 1 THEN BEGIN
		PRINT,'CRISPEX, imcube, spcube, REFCUBE=refcube, LINE_CENTER=line_center, DT=dt, EXTS=exts, $'
		PRINT,'	MNSPEC=mnspec, SINGLE_CUBE=single_cube, SCALE_STOKES=scale_stokes, VALS_REF=vals_ref, $'
		PRINT,'	VALS_IMG=vals_img, NO_WARP=no_warp, SCALE_CUBES=scale_cubes, XTITLE=xtitle, YTITLE=ytitle, $'
		PRINT,'	MASKCUBE=maskcube, WINDOW_LARGE=window_large, VERBOSE=verbose'
		RETURN
	ENDIF

;================================================================================= VERSION AND REVISION NUMBER
	version_number = '1.6.3'
	revision_number = '575'

;================================================================================= PROGRAM VERBOSITY CHECK
	IF (N_ELEMENTS(VERBOSE) NE 1) THEN BEGIN			
		IF (N_ELEMENTS(VERBOSE) GT 1) THEN PRINT,'ERROR: The VERBOSE keyword may only be set to a single integer number. Reverting to default verbosity level 0.'
		verbose = 0
		verbosity = [0,0,0,0,0]
	ENDIF ELSE BEGIN
		verbose >= 0	&	verbose <= 26
	ENDELSE
	verbosity = CRISPEX_DEC2BIN(verbose)

;================================================================================= CRISPEX DIRECTORY CHECK
	file_crispex		= (ROUTINE_INFO('CRISPEX',/SOURCE)).PATH
	dir_crispex 		= FILE_DIRNAME(file_crispex,/MARK_DIRECTORY)
	dir_aux			= dir_crispex+'aux'+PATH_SEP()
	dir_resources		= dir_crispex+'resources'+PATH_SEP()
	dir_buttons		= dir_resources+'buttons'+PATH_SEP()
	dir_settings		= dir_crispex+'settings'+PATH_SEP()
	dir_cpft		= dir_settings+'cpft'+PATH_SEP()
	dir_inst		= dir_settings+'inst'+PATH_SEP()
	dir_cpft_write		= FILE_TEST(dir_cpft, /WRITE)
	dir_inst_write		= FILE_TEST(dir_inst, /WRITE)
	IF (verbosity[1] EQ 1) THEN PRINT,'CRISPEX SETUP: CRISPEX has been compiled from: '+file_crispex

;================================================================================= LOAD PREFERENCES
	default_startupwin = 1		&	default_interpspslice = 1
	default_autoplay = 0		&	default_defsaveid = 0								; 0 = yyyymmdd, 1 = ddmmyyyy
	default_defipath = 0		&	default_defopath = 0								; 0 = local working directory, 1 = saved directory
	default_bgplotcol = 255		&	default_plotcol = 0
	default_phislice_update = 0	&	default_slices_imscale = 0
	cpreffiles = FILE_SEARCH(dir_settings+'crispex.cpref', COUNT = cpreffilecount)
	IF (cpreffilecount GE 1) THEN BEGIN
		RESTORE, cpreffiles[0] 
		IF (verbosity[1] EQ 1) THEN PRINT,'CRISPEX SETUP: Preferences restored from: '+dir_settings+'crispex.cpref'
		resave_preferences = ((N_ELEMENTS(phislice_update) NE 1) OR (N_ELEMENTS(slices_imscale) NE 1))
		IF (N_ELEMENTS(phislice_update) NE 1) THEN phislice_update = default_phislice_update				; Failsafe against older preference files
		IF (N_ELEMENTS(slices_imscale) NE 1) THEN slices_imscale = default_slices_imscale				; Failsafe against older preference files
	ENDIF ELSE BEGIN													; If no preference file is present, set defaults
		startupwin = default_startupwin 		&	interpspslice = default_interpspslice
		autoplay = default_autoplay			&	defsaveid = default_defsaveid
		defipath = default_defipath			&	defopath = default_defopath	
		bgplotcol = default_bgplotcol			&	plotcol = default_plotcol
		phislice_update = default_phislice_update	&	slices_imscale = default_slices_imscale
		resave_preferences = 0
	ENDELSE
;================================================================================= START-UP WINDOW
	screeninfo	= OBJ_NEW('IDLsysMonitorInfo')
	monitors	= screeninfo -> GetNumberOfMonitors()
	screensizes	= screeninfo -> GetRectangles()
	IF (monitors GT 1) THEN BEGIN
		monitor_order = (INDGEN(monitors))[SORT(screensizes[0,0:1])]
		monitor_xsize_order = (INDGEN(monitors))[SORT(screensizes[2,0:1])]
		monitor_ysize_order = (INDGEN(monitors))[SORT(screensizes[3,0:1])]
	ENDIF ELSE BEGIN
		monitor_order = 0
		monitor_ysize_order = 0
	ENDELSE
	x_screen_mid	= screensizes[2,monitor_order[0]]/2.
	y_screen_mid	= screensizes[3,monitor_order[0]]/2.
	startup_im 	= REBIN(REFORM(TOTAL((CRISPEX_READ_BMP_BUTTONS('crispex_startup.bmp',dir_resources))[*,*,1:2],3)),400,300)
	startup_nx 	= (SIZE(startup_im))[1]
	startup_ny 	= (SIZE(startup_im))[2]
	startup_xpos 	= FIX(x_screen_mid-startup_nx/2.)
	startup_ypos 	= FIX(y_screen_mid-startup_ny/2.)
	xout 		= REPLICATE(24,9)
	yout 		= REPLICATE(FIX(startup_ny/2.5)+10,9)-INDGEN(9)*15
	IF startupwin THEN BEGIN
		CRISPEX_WINDOW, startup_nx, startup_ny, 0, 'CRISPEX', startuptlb, startupwid, startup_xpos, startup_ypos, DRAWID = startupdrawid, $
			DRAWBASE = drawbase
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Initializing... '
	ENDIF

;========================== READ-IN AND INITIALISATION OF FILES
  IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Reading input files... '
  ; Check whether machine is big or little
  IF ((BYTE(1L,0,1))[0] EQ 1) THEN endian = 'l' ELSE endian = 'b'
  
  IF N_ELEMENTS(SPCUBE) EQ 1 THEN BEGIN 
		spext = STRMID(spcube,STRPOS(spcube,'.',/REVERSE_SEARCH)+1,STRLEN(spcube))
		IF STRMATCH(spext,'fits',/FOLD_CASE) THEN BEGIN   ; Check whether dealing with fits cube
      CRISPEX_READ_FITSHEADER, spcube, datatype=sptype, nlp=nlp, nt=nt, spnt=spnt, lps=lps, $
                               offset=spoffset, lptitle=xtitle, ttitle=spytitle, dt=dt, $
                               lpunit=lpunit, exten_no=0
      ms = 1.0
			swapvalue = 1
    ENDIF ELSE BEGIN
		  CRISPEX_READ_HEADER, spcube, datatype=sptype, dims=imdims, nx=nlp, ny=nt, nt=spnt, $
                           endian=endian_file, stokes=spstokes, ns=spns, diagnostics=spdiagnostics
      spoffset = 512
  		swapvalue = ((sptype GT 1) AND (endian NE endian_file))
      lpunit = ''
    ENDELSE
    ; Actual read-in of the spectral cube
    OPENR, lur, spcube, /get_lun, SWAP_ENDIAN = swapvalue
    ; Read data from associated file, skip first 512 (header)bytes
    IF (sptype EQ 1) THEN spectra = ASSOC(lur,BYTARR(nlp,nt),spoffset) $
      ELSE IF (sptype EQ 2) THEN spectra = ASSOC(lur,INTARR(nlp,nt),spoffset) $
      ELSE IF (sptype EQ 4) THEN spectra = ASSOC(lur,FLTARR(nlp,nt),spoffset)
    spfile = 1	
    IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: Read main spectral cube: '+spcube+'. Dimensions: (nlp,nt,nx*ny*ns) = ('+STRTRIM(nlp,2)+','+STRTRIM(nt,2)+','+STRTRIM(spnt,2)+').'
	ENDIF ELSE BEGIN
		spfile = 0
		nt = 1
		spcube = ''
		lur = 0
		IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: No spectral cube supplied.'
	ENDELSE

	imext = STRMID(imcube,STRPOS(imcube,'.',/REVERSE_SEARCH)+1,STRLEN(imcube))
	IF STRMATCH(imext,'fits',/FOLD_CASE) THEN BEGIN
    CRISPEX_READ_FITSHEADER, imcube, datatype=imtype, nx=nx, ny=ny, imnt=imnt, ns=ns, $
                             offset=imoffset, inttitle=ytitle, bunit=bunit, dx=dx, dy=dy, $
                             xunit=xunit, yunit=yunit, exten_no=0
    dx_fixed = 1
		swapvalue = 1
	ENDIF ELSE BEGIN
    CRISPEX_READ_HEADER, imcube, datatype=imtype, dims=imdims, nx=nx, ny=ny, nt=imnt, $
                         endian=endian_file, stokes=imstokes, ns=imns, diagnostics=imdiagnostics
		imoffset = 512
    bunit = 'counts'
    arcsecpix = 0.0592
    dx_fixed = 0
    dx = arcsecpix
    dy = arcsecpix
    xunit = 'arcsec'
    yunit = 'arcsec'
    ns = imns
		swapvalue = ((imtype GT 1) AND (endian NE endian_file))									
  ENDELSE
	onecube = 0
	stokesfile = (ns GE 2)
	IF spfile THEN BEGIN
		IF stokesfile THEN BEGIN
			IF ((spstokes NE imstokes) OR (spns NE imns)) THEN BEGIN
				PRINT,'ERROR: IMCUBE and SPCUBE have incompatible Stokes dimensions and seem to belong to different datasets.'
				PRINT,'       Please check whether the input is correct (you provided IMCUBE='+STRTRIM(imcube,2)
				PRINT,'       and SPCUBE='+STRTRIM(spcube,2)+').'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF ELSE ns = imns
		ENDIF ELSE ns = 1
		IF ((spcube EQ imcube) OR (nlp EQ nx) AND (nt EQ ny) AND (spnt EQ imnt)) THEN BEGIN
			PRINT,'ERROR: IMCUBE and SPCUBE must be different. Please check input (you seem to have provided the same file twice).'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF
		IF ((nx*ny*ns NE spnt) OR (nt*nlp*ns NE imnt)) THEN BEGIN							; Failsafe against incompatible IMCUBE and SPCUBE
			PRINT,'ERROR: IMCUBE and SPCUBE have incompatible dimensions and seem to belong to different datasets.'
			PRINT,'       Please check whether the input is correct (you provided IMCUBE='+STRTRIM(imcube,2)
			PRINT,'       and SPCUBE='+STRTRIM(spcube,2)+').'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF
		IF (N_ELEMENTS(SINGLE_CUBE) GT 0) THEN BEGIN
			PRINT,'WARNING: Calling CRISPEX with SINGLE_CUBE, while SPCUBE is provided, is not allowed. SINGLE_CUBE keyword will be ignored.'
			onecube = 0
		ENDIF
		single_cube = 0
	ENDIF ELSE BEGIN 
		IF (N_ELEMENTS(SINGLE_CUBE) GT 0) THEN BEGIN
			IF (N_ELEMENTS(SINGLE_CUBE) EQ 1) THEN BEGIN
				nlp = LONG(SINGLE_CUBE)
				onecube = 1
				nt = imnt / nlp / ns
			ENDIF ELSE BEGIN
				PRINT,'ERROR: SINGLE_CUBE must be provided with a single integer.'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDELSE
		ENDIF ELSE IF (imnt GT 500) THEN BEGIN
			PRINT,'ERROR: Third dimension of IMCUBE is too large ('+STRTRIM(imnt,2)+'). In single argument mode, please be sure to provide either a '
			PRINT,'       single scan. A full timeseries is only allowed with a call to SINGLE_CUBE. Alternatively, switch to double argument mode '
			PRINT,'	      by providing SPCUBE.'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF ELSE BEGIN
			nlp = imnt / ns
			single_cube = 0
		ENDELSE
	ENDELSE
	IF stokesfile THEN BEGIN
		stokes_comp = STRSPLIT(STRSPLIT(STRJOIN(STRSPLIT(imstokes,',',/EXTRACT)),']',/EXTRACT),'[',/EXTRACT)
		IF (STRLEN(stokes_comp) NE imns) THEN BEGIN
			PRINT,'ERROR: The number of Stokes components ('+STRTRIM(imns,2)+') does not correspond to the number of Stokes labels ('+STRTRIM(STRLEN(stokes_comp),2)+').'
			PRINT,'       Please check whether the Stokes cube production has proceded correctly.'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF ELSE BEGIN
			ns = imns
			stokes_labels = STRARR(ns)
			stokes_select_sp = INTARR(ns)
			stokes_labels = STRMID(stokes_comp,INDGEN(ns),1)
			IF ((WHERE(stokes_labels EQ 'I') GE 0) AND (WHERE(stokes_labels EQ 'I') LE imns-1)) THEN BEGIN
				stokes_i_enabled = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'I')] = 1
			ENDIF ELSE stokes_i_enabled = 0
			IF ((WHERE(stokes_labels EQ 'Q') GE 0) AND (WHERE(stokes_labels EQ 'Q') LE imns-1)) THEN BEGIN
				stokes_q_enabled = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'Q')] = 1
			ENDIF ELSE stokes_q_enabled = 0
			IF ((WHERE(stokes_labels EQ 'U') GE 0) AND (WHERE(stokes_labels EQ 'U') LE imns-1)) THEN BEGIN
				stokes_u_enabled = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'U')] = 1
			ENDIF ELSE stokes_u_enabled = 0
			IF ((WHERE(stokes_labels EQ 'V') GE 0) AND (WHERE(stokes_labels EQ 'V') LE imns-1)) THEN BEGIN
				stokes_v_enabled = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'V')] = 1
			ENDIF ELSE stokes_v_enabled = 0
		ENDELSE
	ENDIF ELSE BEGIN
		stokes_i_enabled = 0
		stokes_q_enabled = 0
		stokes_u_enabled = 0
		stokes_v_enabled = 0
		stokes_labels = ['I']
		stokes_select_sp = 1
	ENDELSE
	scalestokes_max = (stokes_q_enabled OR stokes_u_enabled OR stokes_v_enabled)
	diagnostics = STRARR(nlp)
	IF (N_ELEMENTS(imdiagnostics) GT 0) THEN BEGIN
		diagsplit = STRSPLIT(imdiagnostics,',',/EXTRACT)
		ndiag = N_ELEMENTS(diagsplit)
		diagsplit[0] = STRMID(diagsplit[0],STRPOS(diagsplit[0],'[')+1,STRLEN(diagsplit[0]))
		diagsplit[ndiag-1] = STRMID(diagsplit[ndiag-1],0,STRLEN(diagsplit[ndiag-1])-(STRPOS(diagsplit[ndiag-1],']') GT 0))
		IF (ndiag GE nlp) THEN diagnostics = diagsplit[0:(nlp-1)] ELSE IF (ndiag LT nlp) THEN BEGIN
			diagnostics[0:(ndiag-1)] = diagsplit[0:(ndiag-1)]
			diagnostics[ndiag:(nlp-1)] = REPLICATE('Undefined',(nlp-ndiag))
		ENDIF 
	ENDIF ELSE diagnostics = REPLICATE('SST ',nlp)+STRTRIM(INDGEN(nlp),2)
	sel_diagnostics = REPLICATE(1,nlp)
	lines_diagnostics = (INDGEN(nlp) MOD 6)
	linlab_diagnostics = ['Solid', 'Dotted', 'Dashed', 'Dash Dot', 'Dash Dot Dot', 'Long Dashes'] 
	colors_diagnostics = [0,200,135,120,100,90,230,40]
	collab_diagnostics = ['Black', 'Red', 'Pink', 'Purple', 'Blue', 'Turquoise', 'Grey', 'Green']
	FOR i=0,FLOOR(nlp/6.) DO BEGIN
		IF (i EQ 0) THEN selcol_diagnostics = REPLICATE(i,6) ELSE $
			IF (i EQ (FLOOR(nlp/6.))) THEN BEGIN
				IF (nlp/6. NE FLOOR(nlp/6.)) THEN selcol_diagnostics = [selcol_diagnostics, REPLICATE(i,ROUND((nlp/6.-FLOOR(nlp/6.))*6))] 
			ENDIF ELSE selcol_diagnostics = [selcol_diagnostics, REPLICATE(i,6)]
	ENDFOR

	nx = nx * 1L														; Convert x-dimension to LONG
	ny = ny * 1L														; Convert y-dimension to LONG
	OPENR, lun, imcube, /get_lun, SWAP_ENDIAN = swapvalue									; Actual read-in of the image cube
	IF (imtype EQ 1) THEN BEGIN												; Read data from associated file,
		imagefile = ASSOC(lun,BYTARR(nx,ny),imoffset)									; skip imoffset header bytes
		scanfile  = ASSOC(lun,BYTARR(nx,ny,nlp*ns),imoffset)			; Re-read in of the image cube for slices
	ENDIF ELSE IF (imtype EQ 2) THEN BEGIN
		imagefile = ASSOC(lun,INTARR(nx,ny),imoffset)
		scanfile  = ASSOC(lun,INTARR(nx,ny,nlp*ns),imoffset)							
	ENDIF ELSE IF (imtype EQ 4) THEN BEGIN
		imagefile = ASSOC(lun,FLTARR(nx,ny),imoffset)
		scanfile  = ASSOC(lun,FLTARR(nx,ny,nlp*ns),imoffset)								
	ENDIF
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN
		IF stokesfile THEN PRINT,'CRISPEX SETUP: Read Stokes image cube: '+imcube+'. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(imnt,2)+').' ELSE $
			PRINT,'CRISPEX SETUP: Read image cube: '+imcube+'. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(imnt,2)+').'
		IF (verbosity[1] EQ 1) THEN PRINT,'CRISPEX SETUP: Main cubes dimensions: (nx,ny,nt,nlp,ns) = ('+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+','+STRTRIM(nlp,2)+','+STRTRIM(ns,2)+')'
	ENDIF
	IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Stokes parameters: '+STRJOIN(stokes_labels,' ')

	IF (N_ELEMENTS(REFCUBE) EQ 1) THEN BEGIN										; Assuming a single reference cube
  	refimext = STRMID(refcube,STRPOS(refcube,'.',/REVERSE_SEARCH)+1,STRLEN(refcube))
  	IF STRMATCH(refimext,'fits',/FOLD_CASE) THEN BEGIN
      CRISPEX_READ_FITSHEADER, refcube, datatype=refimtype, nx=refnx, ny=refny, imnt=refnt, $
                             offset=refimoffset, bunit=refbunit, exten_no=0
  		swapvalue = 1
  	ENDIF ELSE BEGIN
  		CRISPEX_READ_HEADER, refcube, datatype=refimtype, dims=refdims, nx=refnx, ny=refny, nt=refnt,$
                           endian=endian_file	; Calling LP_HEADER.PRO to read the header of the reference cube
  		refimoffset = 512
      refbunit = 'counts'
  		swapvalue = ((refimtype GT 1) AND (endian NE endian_file))									
    ENDELSE
		IF (refnx EQ nx) AND (refny EQ ny) AND (refnt EQ nt) OR (refnt EQ 1) OR (refnt EQ imnt) OR (refnt EQ nlp) OR (nt EQ 1) THEN BEGIN
;			swapvalue = ((refimtype GT 1) AND (endian NE endian_file))						; and determine whether correction is needed
			OPENR, luf, refcube, /get_lun, SWAP_ENDIAN = swapvalue							; Actual read-in of the reference cube
			IF (refimtype EQ 1) THEN referencefile = ASSOC(luf,BYTARR(nx,ny),refimoffset) $				; Read data from associated file, skip first 512 (header)bytes
			ELSE IF (refimtype EQ 2) THEN referencefile = ASSOC(luf,INTARR(nx,ny),refimoffset) $
			ELSE IF (refimtype EQ 4) THEN referencefile = ASSOC(luf,FLTARR(nx,ny),refimoffset)
			showref = 1
			refspfile = 0
			refspcube = ''
      reflpunit = ''
			lufs = 0
			refns = 1
			IF (((refnt GT nt) AND (nt NE 1)) OR (refnt EQ nlp)) THEN refnlp = nlp ELSE BEGIN
				IF (refnt EQ nt) THEN refnlp = 1 ELSE refnlp = refnt
			ENDELSE
			IF (refnt EQ nlp) OR (nt EQ 1) THEN BEGIN
				refnt = 1
				IF (refimtype EQ 1) THEN refscanfile = ASSOC(luf,BYTARR(nx,ny,refnlp*refns),refimoffset) $				; Read data from associated file, skip first 512 (header)bytes
				ELSE IF (refimtype EQ 2) THEN refscanfile = ASSOC(luf,INTARR(nx,ny,refnlp*refns),refimoffset) $
				ELSE IF (refimtype EQ 4) THEN refscanfile = ASSOC(luf,FLTARR(nx,ny,refnlp*refns),refimoffset)
			ENDIF
		ENDIF ELSE BEGIN
			IF (refnt NE nt) AND (refnt NE 1) THEN BEGIN
				PRINT,'ERROR: Dimensions of the reference cube (['+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+']) are not compatible with those of the image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(imnt,2)+'])!'
				PRINT,'       Number of timesteps times the number of spectral positions must be equal to that of the image cube ('+STRTRIM(imnt,2)+').'
			ENDIF ELSE BEGIN
				PRINT,'ERROR: Dimensions of the reference cube (['+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+']) are not compatible with those of the image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
				PRINT,'       Number of timesteps must be either equal to that of the image ('+STRTRIM(nt,2)+') or to 1.'
			ENDELSE
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDELSE
	ENDIF ELSE IF ((N_ELEMENTS(REFCUBE) EQ 2) AND (SIZE(REFCUBE,/TYPE) EQ 7)) THEN BEGIN					; Assuming full dual cube mode
  	refspext = STRMID(refcube[1],STRPOS(refcube[1],'.',/REVERSE_SEARCH)+1,STRLEN(refcube[1]))
  	refimext = STRMID(refcube[0],STRPOS(refcube[0],'.',/REVERSE_SEARCH)+1,STRLEN(refcube[0]))
  	IF STRMATCH(refspext,'fits',/FOLD_CASE) THEN BEGIN
      CRISPEX_READ_FITSHEADER, refcube[1], datatype=refsptype, nlp=refnlp, nt=refnt, spnt=refspnt, $
                             offset=refspoffset, lptitle=refxtitle, lpunit=reflpunit, exten_no=0
      refms = 1.0
  		refspswapvalue = 1
  	ENDIF ELSE BEGIN
  		CRISPEX_READ_HEADER, refcube[1], datatype=refimtype, dims=refdims, nx=refnx, ny=refny, nt=refnt,$
                           endian=endian_file	; Calling LP_HEADER.PRO to read the header of the reference cube
  		refspoffset = 512
  		refspswapvalue = ((refimtype GT 1) AND (endian NE endian_file))						
      reflpunit = ''
    ENDELSE
  	IF STRMATCH(refimext,'fits',/FOLD_CASE) THEN BEGIN
      CRISPEX_READ_FITSHEADER, refcube[0], datatype=refimtype, nx=refnx, ny=refny, imnt=refimnt, $
                             offset=refimoffset, inttitle=refytitle, bunit=refbunit, exten_no=0
  		refimswapvalue = 1
  	ENDIF ELSE BEGIN
  		CRISPEX_READ_HEADER, refcube[0], datatype=refimtype, dims=refdims, nx=refnx, ny=refny, nt=refnt,$
                           endian=endian_file	; Calling LP_HEADER.PRO to read the header of the reference cube
  		refimoffset = 512
      refbunit = 'counts'
  		refimswapvalue = ((refimtype GT 1) AND (endian NE endian_file))									
    ENDELSE
;		CRISPEX_READ_HEADER, refcube[1], datatype=refsptype, dims=refspdims, nx=refnlp, ny=refnt, nt=refspnt, endian=endian_file	; Calling LP_HEADER.PRO to read the header of the reference cube
;		CRISPEX_READ_HEADER, refcube[0], datatype=refimtype, dims=refimdims, nx=refnx, ny=refny, nt=refimnt, endian=endian_file	; Calling LP_HEADER.PRO to read the header of the reference cube
		; First check whether reference cubes are compatible
		IF ((refcube[1] EQ refcube[0]) OR (refnlp EQ refnx) AND (refnt EQ refny) AND (refspnt EQ refimnt)) THEN BEGIN
			PRINT,'ERROR: The reference image and spectral cubes must be different. Please check input (you seem to have provided the same file twice).'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF
		IF ((refnx*refny NE refspnt) OR (refnt*refnlp NE refimnt)) THEN BEGIN							; Failsafe against incompatible IMCUBE and SPCUBE
			PRINT,'ERROR: The reference image and spectral cubes have incompatible dimensions and seem to belong to different datasets.'
			PRINT,'       Please check whether the input is correct (you provided REFCUBE[0]='+STRTRIM(refcube[0],2)
			PRINT,'       and REFCUBE[1]='+STRTRIM(refcube[1],2)+').'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDIF
		; Then check whether reference cubes are compatible with main cubes
		IF ((refnx EQ nx) AND (refny EQ ny) AND ((refnt EQ nt) OR (refnt EQ 1) OR (refnt EQ imnt) OR (refnt EQ nlp))) THEN BEGIN
;			swapvalue = ((refsptype GT 1) AND (endian NE endian_file))							; Determine whether correction for endian is needed
			OPENR, lufs, refcube[1], /get_lun, SWAP_ENDIAN = refspswapvalue								; Actual read-in of the spectral cube
			IF (refsptype EQ 1) THEN referencespectra = ASSOC(lufs,BYTARR(refnlp,refnt),refspoffset) $						; Read data from associated file, skip first 512 (header)bytes
			ELSE IF (refsptype EQ 2) THEN referencespectra = ASSOC(lufs,INTARR(refnlp,refnt),refspoffset) $
			ELSE IF (refsptype EQ 4) THEN referencespectra = ASSOC(lufs,FLTARR(refnlp,refnt),refspoffset)
			refspfile = 1	
;			swapvalue = ((refimtype GT 1) AND (endian NE endian_file)) 						; and determine whether correction is needed
			OPENR, luf, refcube[0], /get_lun, SWAP_ENDIAN = refimswapvalue							; Actual read-in of the reference cube
			IF (refimtype EQ 1) THEN referencefile = ASSOC(luf,BYTARR(nx,ny),refimoffset) $				; Read data from associated file, skip first 512 (header)bytes
			ELSE IF (refimtype EQ 2) THEN referencefile = ASSOC(luf,INTARR(nx,ny),refimoffset) $
			ELSE IF (refimtype EQ 4) THEN referencefile = ASSOC(luf,FLTARR(nx,ny),refimoffset)
			showref = 1
			IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: Read reference spectral cube: '+refcube[1]+$
				'. Dimensions: (nlp,nt,nx*ny) = ('+STRTRIM(refnlp,2)+','+STRTRIM(refnt,2)+','+STRTRIM(refspnt,2)+').'
			refspcube = refcube[1]
			refcube = refcube[0]
			refns = 1
		ENDIF ELSE BEGIN
			IF ((refnx NE nx) OR (refny NE ny)) THEN BEGIN
				PRINT,'ERROR: Dimensions of the reference cube (['+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+']) are not compatible with those of the main image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
				PRINT,'       Number of pixels in the x- and y-dimension must be equal to those of the main image cube.'
			ENDIF ELSE IF (refnt NE nt) AND (refnt NE 1) THEN BEGIN
				PRINT,'ERROR: Dimensions of the reference cube (['+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+']) are not compatible with those of the main image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
				PRINT,'       Number of timesteps must be either equal to that of the main image cube ('+STRTRIM(nt,2)+') or to 1.'
			ENDIF
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDELSE	
	ENDIF ELSE IF (N_ELEMENTS(REFCUBE) GT 1) THEN BEGIN									; Assuming a single image supplied as an array
		refnx = (SIZE(refcube))[1]
		refny = (SIZE(refcube))[2]
    refbunit = 'counts'
    reflpunit = ''
		IF (SIZE(REFCUBE,/N_DIMENSIONS) EQ 3) THEN refnt = (SIZE(refcube))[3] ELSE refnt = 0
		IF (refnx EQ nx) AND (refny EQ ny) THEN BEGIN
			referencefile = refcube
			showref = 1
			luf = 0
			lufs = 0
			refnlp = 0
			refns = 1
			refspfile = 0
			refspcube = ''
      reflpunit = ''
      refbunit = 'counts'
		ENDIF ELSE BEGIN
			PRINT,'ERROR: Dimensions of the reference cube (['+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+']) are not compatible with those of the main image cube (['+STRTRIM(nx,2)+','+$
				STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
			PRINT,'       Number of pixels in the x- and y-direction must be equal.'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDELSE
	ENDIF ELSE BEGIN 
		showref = 0
		luf = 0
		refnt = 0
		lufs = 0
		refnlp = 0
		refns = 0
		refspfile = 0
		refspcube = ''
    reflpunit = ''
    refbunit = ''
	ENDELSE
	IF (refnlp NE nlp) THEN BEGIN
		eqnlps = 0 
		refslid_sens = (showref AND (refnlp GT 1)) 
	ENDIF ELSE BEGIN
		eqnlps = 1
		refslid_sens = 0
	ENDELSE
	showrefls = (refspfile OR (refnlp GT 1))

	IF (N_ELEMENTS(MASKCUBE) EQ 1) THEN BEGIN
		CRISPEX_READ_HEADER, maskcube, datatype=masktype, dims=maskdims, nx=masknx, ny=maskny, nt=masknt, endian=endian_file;, $	; Calling LP_HEADER.PRO to read the header of the spectral cube
		IF ((masknx EQ nx) AND (maskny EQ ny) AND ((masknt EQ nt) OR (masknt EQ 1))) THEN BEGIN
			swapvalue = ((masktype GT 1) AND (endian NE endian_file))								; Determine whether correction for endian is needed
			OPENR, lum, maskcube, /get_lun, SWAP_ENDIAN = swapvalue								; Actual read-in of the spectral cube
			IF (masktype EQ 1) THEN mask = ASSOC(lum,BYTARR(nx,ny),hoffset) $						; Read data from associated file, skip first 512 (header)bytes
			ELSE IF (masktype EQ 2) THEN mask = ASSOC(lum,INTARR(nx,ny),hoffset) $
			ELSE IF (masktype EQ 4) THEN mask = ASSOC(lum,FLTARR(nx,ny),hoffset)
			maskfile = 1	
			IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: Read mask cube: '+maskcube+'. Dimensions: (nx,ny,nt) = ('+STRTRIM(masknx,2)+','+STRTRIM(maskny,2)+','+STRTRIM(masknt,2)+').'
		ENDIF ELSE BEGIN
			IF ((masknx NE nx) OR (maskny NE ny)) THEN BEGIN
				PRINT,'ERROR: Dimensions of the mask cube (['+STRTRIM(masknx,2)+','+STRTRIM(maskny,2)+','+STRTRIM(masknt,2)+']) are not compatible with those of the main image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
				PRINT,'       Number of pixels in the x- and y-dimension must be equal to those of the main image cube.'
			ENDIF ELSE IF ((masknt NE nt) AND (masknt NE 1)) THEN BEGIN
				PRINT,'ERROR: Dimensions of the mask cube (['+STRTRIM(masknx,2)+','+STRTRIM(maskny,2)+','+STRTRIM(masknt,2)+']) are not compatible with those of the main image '+$
					'cube (['+STRTRIM(nx,2)+','+STRTRIM(ny,2)+','+STRTRIM(nt,2)+'])!'
				PRINT,'       Number of timesteps must be either equal to that of the main image cube ('+STRTRIM(nt,2)+') or to 1.'
			ENDIF
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDELSE
	ENDIF ELSE BEGIN
		maskfile = 0
		maskcube = ''
		lum = 0
		mask = 0
		masknt = 0
		IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: No spectral cube supplied.'
	ENDELSE

	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN
		IF (showref EQ 1) THEN BEGIN
			IF (refspfile NE 1) THEN PRINT,'CRISPEX SETUP: Read reference image cube: '+refcube+'. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt*refnlp,2)+').'  ELSE $
				IF (refspfile EQ 1) THEN PRINT,'CRISPEX SETUP: Read reference image cube: '+refcube[0]+'. Dimensions: (nx,ny,nt*nlp) = ('+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+').'
			IF (verbosity[1] EQ 1) THEN PRINT,'CRISPEX SETUP: Reference cubes dimensions: (nx,ny,nt,nlp,ns) = ('+STRTRIM(refnx,2)+','+STRTRIM(refny,2)+','+STRTRIM(refnt,2)+','+STRTRIM(refnlp,2)+','+$
				STRTRIM(refns,2)+')'
		ENDIF ELSE IF (verbosity[1] EQ 1) THEN PRINT,'CRISPEX SETUP: No reference image cube supplied.'
	ENDIF
	IF (showref AND startupwin) THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Reading input files... done!'


	DEVICE, DECOMPOSE = 0	

	SPAWN, 'echo $HOSTNAME', hostname
	cpftfile = FILE_SEARCH(dir_cpft+'crispex.'+hostname+'.cpft', COUNT = cpftfilecount)
	IF cpftfilecount THEN BEGIN
		RESTORE, cpftfile[0] 
		IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Restored '+cpftfile[0]+'.'
	ENDIF ELSE BEGIN
		IF (verbosity[1] EQ 1) THEN BEGIN
			PRINT, 'CRISPEX SETUP: No CRISPEX performance test file (crispex.'+hostname+'.cpft) found to restore '
			PRINT, '               in '+dir_cpft
		ENDIF
		estimate_lx = 0
		estimate_time = 0.
		estimate_run = 0
	ENDELSE 

	instfilename = 'crispex.'+hostname+'.inst'
	IF dir_inst_write THEN BEGIN
		instfile = FILE_SEARCH(dir_inst+instfilename, COUNT = instfilecount)
		IF instfilecount THEN BEGIN
			IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Opening existing instance tracking file: '+instfilename+'.'
			nlines = FILE_LINES(instfile)
			datarr = STRARR(1,nlines)
			OPENR,unit1,instfile,/GET_LUN
			READF,unit1,datarr
			FREE_LUN,unit1
			routine_name = STRARR(nlines)
			instance_id = LONARR(nlines)
			FOR i=1,nlines[0]-1 DO BEGIN
				splitline = STRSPLIT(datarr[i],'	',/EXTRACT)
				routine_name[i] = splitline[0]
				instance_id[i] = splitline[3]
			ENDFOR
			where_crispex = WHERE(routine_name EQ 'CRISPEX')
			OPENU, unit2, dir_inst+instfilename, WIDTH = 360, /GET_LUN, /APPEND
		ENDIF ELSE BEGIN
			IF (verbosity[1] EQ 1) THEN BEGIN
				PRINT, 'CRISPEX SETUP: No CRISPEX instance tracking file ('+instfilename+') found in '+dir_inst
				PRINT, '               Creating file.'
			ENDIF
			where_crispex = -1
			OPENW, unit2, dir_inst+instfilename, WIDTH = 360, /GET_LUN
			PRINTF, unit2, '# routine_name	version		revision	ID'
		ENDELSE
		IF (where_crispex[0] NE -1) THEN set_instance_id = STRTRIM((instance_id[where_crispex])[WHERE(instance_id[where_crispex] EQ MAX(instance_id[where_crispex]))] + 1,2) ELSE set_instance_id = STRTRIM(0,2)
		PRINTF, unit2, 'CRISPEX	'+version_number+'	'+revision_number+'	'+set_instance_id
		FREE_LUN, unit2
		IF (set_instance_id GE 1) THEN instance_label = '-'+set_instance_id ELSE instance_label = ''
		IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Written instance ID ('+set_instance_id+') to '+instfilename+'.'
	ENDIF ELSE BEGIN
		set_instance_id = ''
		instance_label = ''
		instfilecount = 0
		IF (verbosity[1] EQ 1) THEN BEGIN
			PRINT, 'CRISPEX SETUP ERROR: Could not write CRISPEX instance tracking file '+instfilename 
			PRINT, '                     to '+dir_inst+'. Permission denied.'
		ENDIF
	ENDELSE
	
;================================================================================= SETTING START-UP OPTIONS 
;--------------------------------------------------------------------------------- PARAMETERS FROM MEAN SPEC
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (parameters from/for mean spectrum)...",a5)','     ') 
	feedback_text = ['Setting start-up options... ','> Parameters from/for mean spectrum... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	detspect_scale_enable = (nlp GT 1)
	detspect_scale = (nlp GT 1)
	ref_detspect_scale = (refnlp GT 1)

	IF (N_ELEMENTS(SCALE_STOKES) EQ 1) THEN scalestokes = 1 ELSE BEGIN
		IF (N_ELEMENTS(SCALE_STOKES) GT 1) THEN PRINT, 'ERROR: The SCALE_STOKES keyword must either be set or supplied with a single number! Reverting to standard scaling.'
		scalestokes = 0
	ENDELSE

	IF (N_ELEMENTS(XTITLE) EQ 1) THEN BEGIN
    IF (N_ELEMENTS(REFXTITLE) EQ 1) THEN $
      xtitle = [xtitle[0],refxtitle[0]] $
    ELSE $
      xtitle = [xtitle[0],''] 
  ENDIF ELSE IF (N_ELEMENTS(XTITLE) EQ 2) THEN $
    xtitle = xtitle $
  ELSE IF (N_ELEMENTS(REFXTITLE) EQ 1) THEN $
    xtitle = ['',refxtitle[0]] $
  ELSE $
    xtitle = ['','']
	IF (N_ELEMENTS(YTITLE) EQ 1) THEN BEGIN
    IF (N_ELEMENTS(REFYTITLE) EQ 1) THEN $
      ytitle = [ytitle[0],refytitle[0]] $
    ELSE $
      ytitle = [ytitle[0],''] 
  ENDIF ELSE IF (N_ELEMENTS(YTITLE) EQ 2) THEN $
    ytitle = ytitle $
  ELSE IF (N_ELEMENTS(REFYTITLE) EQ 1) THEN $
    ytitle = ['',refytitle[0]] $
  ELSE $
    ytitle = ['','']

	IF ((N_ELEMENTS(SPECTFILE) EQ 1) OR (N_ELEMENTS(SPECTFILE) EQ 2)) THEN BEGIN										; If SPECTFILE is specified, use that to
		IF (N_ELEMENTS(MNSPEC) GT 0) THEN PRINT,'WARNING: Calling CRISPEX with MNSPEC, while SPECTFILE is provided, is not allowed. Using SPECTFILE for mean spectrum determination.'
		IF (spectfile[0] NE '') THEN BEGIN
			RESTORE, spectfile[0]												; determine:
			IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Restored main spectral file: '+spectfile[0]+'.'
			IF (N_ELEMENTS(norm_factor) GE 1) THEN BEGIN										; Failsafe against old spectfiles
				IF (N_ELEMENTS(spect_pos[*,0]) NE nlp) THEN BEGIN										; (after failsave against incompatible SPECTFILE
					PRINT,'ERROR: Number of spectral positions in SPECTFILE (nlp='+STRTRIM(N_ELEMENTS(spect_pos),2)+') is incompatible with that in the datacubes (nlp='+STRTRIM(nlp,2)+').'
					PRINT,'       Please load the correct spectral file.'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				mainspec = norm_spect													; - the normalised spectrum
				spectrum = norm_spect[*,0]
				ms	= norm_factor												; - roughly twice the maximum of "continuum"
				lps	= spect_pos													; - the values of the line positions
				IF (N_ELEMENTS(xtitle_label) EQ 1) THEN BEGIN
					IF ((STRCOMPRESS(xtitle_label) NE '') AND (STRCOMPRESS(xtitle[0]) EQ '')) THEN xtitle[0] = xtitle_label
				ENDIF
				IF (N_ELEMENTS(ytitle_label) EQ 1) THEN BEGIN
					IF ((STRCOMPRESS(ytitle_label) NE '') AND (STRCOMPRESS(ytitle[0]) EQ '')) THEN ytitle[0] = ytitle_label
				ENDIF
			ENDIF ELSE BEGIN
				IF (N_ELEMENTS(ll) NE nlp) THEN BEGIN										; (after failsave against incompatible SPECTFILE
					PRINT,'ERROR: Number of spectral positions in SPECTFILE (nlp='+STRTRIM(N_ELEMENTS(ll),2)+') is incompatible with that in the datacubes (nlp='+STRTRIM(nlp,2)+').'
					PRINT,'       Please load the correct spectral file.'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				mainspec = spec
				spectrum = spec
				ms 	= 1 / mn
				lps 	= ll
				PRINT,'WARNING: The restored main spectral file ('+spectfile[0]+') has a format pre-dating CRISPEX v1.6.'
				PRINT,'         If possible, please update the main spectral file to the newest format. Reading spectral file in compatibility mode.'
			ENDELSE
			spectfile_set = 1
			refspec = 0	
		ENDIF
		IF (N_ELEMENTS(SPECTFILE) EQ 2) THEN BEGIN
			norm_spect = 0	&	norm_factor = [0,0]	&	spect_pos = 0
			RESTORE, spectfile[1]												; determine:
			IF (N_ELEMENTS(norm_factor) EQ 1) THEN BEGIN										; Failsafe against old spectfiles
				IF (verbosity[1] EQ 1) THEN PRINT, 'CRISPEX SETUP: Restored reference spectral file: '+spectfile[1]+'.'
				IF (N_ELEMENTS(spect_pos) NE refnlp) THEN BEGIN										; (after failsave against incompatible SPECTFILE
					PRINT,'ERROR: Number of spectral positions in SPECTFILE (nlp='+STRTRIM(N_ELEMENTS(spect_pos),2)+') is incompatible with that in the reference datacubes (nlp='+STRTRIM(refnlp,2)+').'
					PRINT,'       Please load the correct reference spectral file.'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				refspec	= norm_spect												; - the normalised spectrum
				refms	= norm_factor												; - roughly twice the maximum of "continuum"
				reflps	= spect_pos													; - the values of the line positions
				IF (N_ELEMENTS(xtitle_label) EQ 1) THEN BEGIN
					IF ((STRCOMPRESS(xtitle_label) NE '') AND (STRCOMPRESS(xtitle[1]) EQ '')) THEN xtitle[1] = xtitle_label
				ENDIF
				IF (N_ELEMENTS(ytitle_label) EQ 1) THEN BEGIN
					IF ((STRCOMPRESS(ytitle_label) NE '') AND (STRCOMPRESS(ytitle[1]) EQ '')) THEN ytitle[1] = ytitle_label
				ENDIF
			ENDIF ELSE BEGIN
				IF (N_ELEMENTS(ll) NE refnlp) THEN BEGIN										; (after failsave against incompatible SPECTFILE
					PRINT,'ERROR: Number of spectral positions in SPECTFILE (nlp='+STRTRIM(N_ELEMENTS(ll),2)+') is incompatible with that in the reference datacubes (nlp='+STRTRIM(refnlp,2)+').'
					PRINT,'       Please load the correct reference spectral file.'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				refspec = spec
				refms	= 1 / mn
				reflps	= ll
				PRINT,'WARNING: The restored reference spectral file ('+spectfile[1]+') has a format pre-dating CRISPEX v1.6.'
				PRINT,'         If possible, please update the reference spectral file to the newest format. Reading spectral file in compatibility mode.'
			ENDELSE
			refspectfile_set = 1
		ENDIF ELSE BEGIN
			refspectfile_set = 0
			IF showrefls THEN BEGIN
				refspectrum = DBLARR(refnlp) 
				IF (N_ELEMENTS(MNSPEC) EQ 1) THEN BEGIN										; - the spectrum of the entire image at a given t
					IF (mnspec[0] LT 0) OR (mnspec[0] GE nt) THEN BEGIN						; Failsafe against incorrect MNSPEC values
						PRINT,'ERROR: MNSPEC value ('+STRTRIM(mnspec[0],2)+') must fall within allowed range [0,'+STRTRIM(LONG(nt-1),2)+']!'
						WIDGET_CONTROL, startuptlb, /DESTROY
						RETURN
					ENDIF	
					FOR k=0,refnlp-1 DO refspectrum[k] = MEAN(referencefile[mnspec[0]*refnlp + k])
				ENDIF ELSE IF (N_ELEMENTS(MNSPEC) EQ 2) THEN BEGIN								; - the spectrum of the entire image at a range of t
					IF (mnspec[0] GT mnspec[1]) OR (mnspec[0] LT 0) OR (mnspec[1] GE nt) THEN BEGIN			; Failsafe against incorrect MNSPEC values
						PRINT,'ERROR: MNSPEC values (['+STRTRIM(mnspec[0],2)+','+STRTRIM(mnspec[1],2)+']) must fall within allowed range [0,'+$
							STRTRIM(nt-1,2)+'] and be ordered from lower to higher value!'
						WIDGET_CONTROL, startuptlb, /DESTROY
						RETURN
					ENDIF	
					FOR k=0,refnlp-1 DO BEGIN
						refspect = 0.
						FOR i=mnspec[0],mnspec[1] DO refspect += MEAN(referencefile[i*refnlp + k])
						refspectrum[k] = refspect/FLOAT(mnspec[1]-mnspec[0]+1)
					ENDFOR
				ENDIF ELSE FOR k=0,refnlp-1 DO refspectrum[k] = MEAN(referencefile[k])
				refms = MAX(refspectrum)
				refspec = refspectrum / refms
				reflps	= FINDGEN(refnlp)
			ENDIF ELSE BEGIN
				refms = 0
				refspec = 0
				reflps = 0
				reflc = 0
				refspxtitle = 0
			ENDELSE 
		ENDELSE
	ENDIF ELSE BEGIN													; Else, determine from the datacubes:
		spectrum = DBLARR(nlp,ns) 
		IF showrefls THEN refspectrum = DBLARR(refnlp) ELSE refspectrum = 0
		IF (N_ELEMENTS(MNSPEC) EQ 1) THEN BEGIN										; - the spectrum of the entire image at a given t
			IF (mnspec[0] LT 0) OR (mnspec[0] GE nt) THEN BEGIN						; Failsafe against incorrect MNSPEC values
				PRINT,'ERROR: MNSPEC value ('+STRTRIM(mnspec[0],2)+') must fall within allowed range [0,'+STRTRIM(LONG(nt-1),2)+']!'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF	
			FOR k=0,nlp-1 DO FOR j=0,ns-1 DO spectrum[k,j] = MEAN(imagefile[mnspec[0]*nlp*ns + j*nlp + k])
			IF showrefls THEN FOR k=0,refnlp-1 DO refspectrum[k] = MEAN(referencefile[mnspec[0]*refnlp + k])
		ENDIF ELSE IF (N_ELEMENTS(MNSPEC) EQ 2) THEN BEGIN								; - the spectrum of the entire image at a range of t
			IF (mnspec[0] GT mnspec[1]) OR (mnspec[0] LT 0) OR (mnspec[1] GE nt) THEN BEGIN			; Failsafe against incorrect MNSPEC values
				PRINT,'ERROR: MNSPEC values (['+STRTRIM(mnspec[0],2)+','+STRTRIM(mnspec[1],2)+']) must fall within allowed range [0,'+$
					STRTRIM(nt-1,2)+'] and be ordered from lower to higher value!'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF	
			FOR k=0,nlp-1 DO BEGIN
				FOR j=0,ns-1 DO BEGIN
					spect = 0.
					FOR i=mnspec[0],mnspec[1] DO spect += MEAN(imagefile[i*nlp*ns + j*nlp + k])
					spectrum[k,j] = spect/FLOAT(mnspec[1]-mnspec[0]+1)
				ENDFOR
			ENDFOR
			IF showrefls THEN BEGIN
				FOR k=0,refnlp-1 DO BEGIN
					refspect = 0.
					FOR i=mnspec[0],mnspec[1] DO refspect += MEAN(referencefile[i*refnlp + k])
					refspectrum[k] = refspect/FLOAT(mnspec[1]-mnspec[0]+1)
				ENDFOR
			END
		ENDIF ELSE BEGIN												; - or the spectrum of the entire image at t=0
			FOR k=0,nlp-1 DO FOR j=0,ns-1 DO spectrum[k,j] = MEAN(imagefile[j*nlp + k])
			IF showrefls THEN FOR k=0,refnlp-1 DO refspectrum[k] = MEAN(referencefile[k])
		ENDELSE

		IF scalestokes THEN BEGIN
			ms = scale_stokes 
			detspect_scale_enable = ((ms EQ 1.) AND (nlp GT 1)) 
		ENDIF ELSE ms = DBLARR(ns) 
		mainspec = DBLARR(nlp,ns)
		FOR j=0,ns-1 DO BEGIN
			IF scalestokes THEN mainspec[*,j] = spectrum[*,j]/ms ELSE BEGIN					; - the normalised spectrum
				min_spectrum = MIN(spectrum[*,j], MAX=max_spectrum)
				IF (ABS(max_spectrum) GE ABS(min_spectrum)) THEN extremum = ABS(max_spectrum) ELSE extremum = ABS(min_spectrum)
				ms[j]		= extremum									; - the maximum of "continuum"
				mainspec[*,j]	= spectrum[*,j]/ms[j]								; - the normalised spectrum
			ENDELSE 
		ENDFOR
		spectrum = spectrum[*,0]
		IF (N_ELEMENTS(lps) NE nlp) THEN lps = FINDGEN(nlp)												; - the values of the line positions
		IF showrefls THEN BEGIN
			refms = MAX(refspectrum)
			refspec = refspectrum / refms
			reflps	= FINDGEN(refnlp)
		ENDIF ELSE BEGIN
			refms = 0
			refspec = 0
			reflps = 0
			reflc = 0
			refspxtitle = 0
		ENDELSE
		spectfile_set = 0
		refspectfile_set = 0
	ENDELSE

	c_speed	= 2.99792458D5													; Speed of light indeed in km/s
	v_dop_set = 0														; Standard setting: No Doppler velocity labelling
	v_dop_set_ref = 0													; Standard setting: No reference Doppler velocity labelling
	dlambda_set = 0
	dlambda_set_ref = 0
	IF (N_ELEMENTS(LINE_CENTER) EQ 0) THEN BEGIN										; If the LINE_CENTER keyword is not set:
		lc	= ( WHERE( spectrum EQ MIN(spectrum) ) )[0]								; autodetermine linecentre value from spectrum
		IF (spectfile_set EQ 1) THEN BEGIN
			v_dop_set = 1 
			spxtitle = 'Wavelength'
		ENDIF ELSE BEGIN
			lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
			spxtitle= 'Spectral position'
		ENDELSE
		IF showrefls THEN BEGIN												; do the same for reference spectrum if refspfile is given
			reflc	= ( WHERE( refspec EQ MIN(refspec) ) )[0]							
			IF (spectfile_set EQ 1) THEN BEGIN
				v_dop_set_ref = 1 
				refspxtitle = 'Wavelength'
			ENDIF ELSE BEGIN
				reflps	= ( reflps - reflps[LONG(reflc)] )
				refspxtitle = 'Spectral position'
			ENDELSE
		ENDIF
	ENDIF ELSE IF ((SIZE(LINE_CENTER))[0] LE 1) THEN BEGIN									; if only set for main
		IF (N_ELEMENTS(LINE_CENTER) EQ 1) THEN BEGIN									; else, if the position is supplied
			lc 	= line_center[0]											; use that given value
			IF ((lc GE nlp) OR (lc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
				PRINT,'ERROR: Linecentre index value '+STRTRIM(LONG(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF
			IF (spectfile_set EQ 1) THEN BEGIN
				v_dop_set = 1 
				spxtitle = 'Wavelength'
			ENDIF ELSE BEGIN
				lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
				spxtitle= 'Spectral position'
			ENDELSE
		ENDIF ELSE IF (N_ELEMENTS(LINE_CENTER) EQ 2) THEN BEGIN									; else, if also the wavelength is supplied
			IF (spectfile_set EQ 1) THEN BEGIN
				PRINT,'ERROR: Wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF	
			lc	 = ( WHERE( spectrum EQ MIN(spectrum) ) )[0]								; autodetermine linecentre value from spectrum
			lambda_c= line_center[0]											; get the linecentre wavelength
			dlambda	= line_center[1]											; and get delta lambda per lineposition
			lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
			lps	= lps*dlambda + lambda_c											; reset scale to wavelength scale
			spxtitle= 'Wavelength'												; and adapt the xtitle accordingly
			v_dop_set = 1
			dlambda_set = 1
		ENDIF ELSE BEGIN													; else, if linecentre and wavelength supplied
			IF (spectfile_set EQ 1) THEN BEGIN
				PRINT,'ERROR: Wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF	
			lc 	= line_center[0]											; get the linecentre position
			IF ((lc GE nlp) OR (lc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
				PRINT,'ERROR: Linecentre index value '+STRTRIM(LONG(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
				WIDGET_CONTROL, startuptlb, /DESTROY
				RETURN
			ENDIF
			lambda_c= line_center[1]											; and the linecentre wavelength
			dlambda	= line_center[2]											; and get delta(wavelength) per lineposition
			lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
			lps	= lps*dlambda + lambda_c											; reset scale to wavelength scale
			spxtitle= 'Wavelength'												; and adapt the xtitle accordingly
			v_dop_set = 1
			dlambda_set = 1
		ENDELSE
		IF showrefls THEN BEGIN												; do the same for reference spectrum if refspfile is given
			reflc	= ( WHERE( refspec EQ MIN(refspec) ) )[0]							
			reflps	= ( reflps - reflps[LONG(reflc)] )
			refspxtitle = 'Spectral position'
		ENDIF
	ENDIF ELSE IF ((SIZE(LINE_CENTER))[0] EQ 2) THEN BEGIN										; if set for main and for reference
		IF showrefls THEN BEGIN
			IF ((SIZE(LINE_CENTER))[1] EQ 1) THEN BEGIN									; else, if the position is supplied
				lc 	= line_center[0,0]											; use that given value
				IF ((lc GE nlp) OR (lc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
					PRINT,'ERROR: Linecentre index value '+STRTRIM(LONG(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				IF (spectfile_set EQ 1) THEN BEGIN
					v_dop_set = 1 
					spxtitle = 'Wavelength'
				ENDIF ELSE BEGIN
					lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
					spxtitle= 'Spectral position'
				ENDELSE
				reflc	= line_center[0,1]
				IF ((reflc GE refnlp) OR (reflc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
					PRINT,'ERROR: Reference linecentre index value '+STRTRIM(LONG(reflc),2)+' falls outside of allowed range [0,'+STRTRIM(refnlp-1,2)+']!'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				IF (refspectfile_set EQ 1) THEN BEGIN
					v_dop_set_ref = 1 
					refspxtitle = 'Wavelength'
				ENDIF ELSE BEGIN
					reflps	= ( reflps - reflps[LONG(reflc)] ) 											; reset scale to have lps=0 at lc
					refspxtitle= 'Spectral position'
				ENDELSE
			ENDIF ELSE IF ((SIZE(LINE_CENTER))[1] EQ 2) THEN BEGIN									; else, if also the wavelength is supplied
				IF (spectfile_set EQ 1) THEN BEGIN
					PRINT,'ERROR: Wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF	
				lc	 = ( WHERE( spectrum EQ MIN(spectrum) ) )[0]								; autodetermine linecentre value from spectrum
				lambda_c= line_center[0,0]											; get the linecentre wavelength
				dlambda	= line_center[1,0]											; and get delta lambda per lineposition
				lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
				lps	= lps*dlambda + lambda_c											; reset scale to wavelength scale
				spxtitle= 'Wavelength'												; and adapt the xtitle accordingly
				v_dop_set = 1
				dlambda_set = 1
				IF (refspectfile_set EQ 1) THEN BEGIN
					PRINT,'ERROR: Reference wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF	
				reflc	 	= ( WHERE( refspec EQ MIN(refspec) ) )[0]								; autodetermine linecentre value from spectrum
				ref_lambda_c	= line_center[0,1]											; get the linecentre wavelength
				ref_dlambda	= line_center[1,1]											; and get delta lambda per lineposition
				reflps		= ( reflps - reflps[LONG(reflc)] ) 											; reset scale to have lps=0 at lc
				reflps		= reflps*ref_dlambda + ref_lambda_c											; reset scale to wavelength scale
				refspxtitle	= 'Wavelength'												; and adapt the xtitle accordingly
				v_dop_set_ref 	= 1
				dlambda_set_ref	= 1
			ENDIF ELSE IF ((SIZE(LINE_CENTER))[1] EQ 3) THEN BEGIN													; else, if linecentre and wavelength supplied
				IF (spectfile_set EQ 1) THEN BEGIN
					PRINT,'ERROR: Wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF	
				lc 	= line_center[0,0]											; get the linecentre position
				IF ((lc GE nlp) OR (lc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
					PRINT,'ERROR: Linecentre index value '+STRTRIM(LONG(lc),2)+' falls outside of allowed range [0,'+STRTRIM(nlp-1,2)+']!'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				lambda_c= line_center[1,0]											; and the linecentre wavelength
				dlambda	= line_center[2,0]											; and get delta(wavelength) per lineposition
				lps	= ( lps - lps[LONG(lc)] ) 											; reset scale to have lps=0 at lc
				lps	= lps*dlambda + lambda_c											; reset scale to wavelength scale
				spxtitle= 'Wavelength'												; and adapt the xtitle accordingly
				v_dop_set = 1
				dlambda_set = 1
				IF (refspectfile_set EQ 1) THEN BEGIN
					PRINT,'ERROR: Reference wavelength information may not be specified in both SPECTFILE and LINE_CENTER! '
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF	
				reflc 	= line_center[0,1]											; get the linecentre position
				IF ((reflc GE refnlp) OR (reflc LT 0)) THEN BEGIN										; Check whether 0 LE lc LT nlp
					PRINT,'ERROR: Reference linecentre index value '+STRTRIM(LONG(reflc),2)+' falls outside of allowed range [0,'+STRTRIM(refnlp-1,2)+']!'
					WIDGET_CONTROL, startuptlb, /DESTROY
					RETURN
				ENDIF
				ref_lambda_c	= line_center[1,1]											; and the linecentre wavelength
				ref_dlambda	= line_center[2,1]											; and get delta(wavelength) per lineposition
				reflps		= ( reflps - reflps[LONG(reflc)] ) 											; reset scale to have lps=0 at lc
				reflps		= reflps*ref_dlambda + ref_lambda_c											; reset scale to wavelength scale
				refspxtitle	= 'Wavelength'												; and adapt the xtitle accordingly
				v_dop_set_ref 	= 1
				dlambda_set_ref = 1
			ENDIF
		ENDIF ELSE BEGIN
			PRINT,'ERROR: LINE_CENTER keyword contains too many elements for the number of cubes provided!'
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		ENDELSE 
	ENDIF

	IF (lps[LONG(lc)] EQ 0) THEN v_dop_set = 0
	IF (reflps[LONG(reflc)] EQ 0) THEN v_dop_set_ref = 0

	IF (v_dop_set EQ 1) THEN v_dop = c_speed*(lps/lps[LONG(lc)]-1) ELSE v_dop = FLTARR(nlp)					; array with Doppler velocities in km/s
	IF showrefls THEN BEGIN
		IF (v_dop_set_ref EQ 1) THEN v_dop_ref = c_speed*(reflps/reflps[LONG(reflc)]-1) ELSE v_dop_ref = FLTARR(refnlp)					; array with Doppler velocities in km/s
	ENDIF ELSE v_dop_ref = 0

	IF (nlp GT 1) THEN BEGIN
		IF dlambda_set THEN ndecimals = ABS(FLOOR(ALOG10(ABS(dlambda)))) ELSE ndecimals = 2
		equidist = STRING((SHIFT(FLOAT(lps),-1) - FLOAT(lps))[0:nlp-2],FORMAT='(F8.'+STRTRIM(ndecimals,2)+')')
		IF (((WHERE(equidist NE equidist[0]))[0] NE -1) AND (KEYWORD_SET(NO_WARP) EQ 0)) THEN BEGIN				; Check for non-equidistant spectral positions and allowed consequential warping
			warpspslice = 1													; Temporal spectrum is warped to correct non-equidistant spectral positions
			min_lps = MIN(lps)
			xo = FINDGEN(nlp) # REPLICATE(1,nt)
			xi = ((lps-min_lps) / FLOAT(MAX(lps-min_lps)) * nlp) # REPLICATE(1,nt)
			yo = REPLICATE(1,nlp) # FINDGEN(nt)
			yi = yo
		ENDIF ELSE BEGIN
			warpspslice = 0													; Temporal spectrum is not warped to correct non-equidistant spectral positions
			xi = 0	&	yi = 0
			xo = 0	&	yo = 0
		ENDELSE
	ENDIF ELSE BEGIN
		warpspslice = 0													; Temporal spectrum is not warped to correct non-equidistant spectral positions
		xi = 0	&	yi = 0
		xo = 0	&	yo = 0
	ENDELSE
	IF (refspfile AND (refnlp GT 1)) THEN BEGIN
		IF dlambda_set_ref THEN ndecimals = ABS(FLOOR(ALOG10(ABS(ref_dlambda)))) ELSE ndecimals = 2
		refequidist = STRING((SHIFT(FLOAT(reflps),-1) - FLOAT(reflps))[0:refnlp-2],FORMAT='(F8.'+STRTRIM(ndecimals,2)+')')
		IF (((WHERE(refequidist NE refequidist[0]))[0] NE -1) AND (KEYWORD_SET(NO_WARP) EQ 0)) THEN BEGIN				; Check for non-equidistant spectral positions and allowed consequential warping
			warprefspslice = 1													; Temporal spectrum is warped to correct non-equidistant spectral positions
			min_reflps = MIN(reflps)
			xo_ref = FINDGEN(refnlp) # REPLICATE(1,nt)
			xi_ref = ((reflps-min_reflps) / FLOAT(MAX(reflps-min_reflps)) * refnlp) # REPLICATE(1,nt)
			yo_ref = REPLICATE(1,refnlp) # FINDGEN(nt)
			yi_ref = yo_ref
		ENDIF ELSE BEGIN
			warprefspslice = 0													; Temporal spectrum is not warped to correct non-equidistant spectral positions
			xi_ref = 0	&	yi_ref = 0
			xo_ref = 0	&	yo_ref = 0
		ENDELSE
	ENDIF ELSE BEGIN
		warprefspslice = 0													; Temporal spectrum is not warped to correct non-equidistant spectral positions
		xi_ref = 0	&	yi_ref = 0
		xo_ref = 0	&	yo_ref = 0
	ENDELSE
	
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (parameters from/for mean spectrum)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Parameters from/for mean spectrum... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SLIT PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial slit parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initial slit parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	nphi	= LONG(CEIL(SQRT( FLOAT(nx)^2 + FLOAT(ny)^2 )))									; Determine maximum number of slitpositions
	angle = 45														; Set initial angle of the slit
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial slit parameters)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial slit parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL PLAYBACK PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial playback parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initial playback parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	t_first		= 0													; Set number of first frame		
	IF (nt EQ 1) THEN BEGIN
		t_last = 2
		t_last_tmp = 0
		t_slid_sens = 0
	ENDIF ELSE BEGIN
		t_last = nt-1													; Set number of last frame
		t_last_tmp = t_last
		t_slid_sens = 1
	ENDELSE
	t_step		= 1													; Set initial timestep
	t_speed 	= 10													; Set initial animation speed
	direction 	= 1													; Set initial animation direction
	nt		= FLOAT(nt)												; Convert the number of timesteps to float
	t_start = t_first

	IF (N_ELEMENTS(DT) EQ 1) THEN BEGIN
		IF (spfile OR onecube) THEN BEGIN
			dt_set = 1
			IF (N_ELEMENTS(SPYTITLE) NE 1) THEN spytitle = 'Time (s)'
		ENDIF ELSE BEGIN
			PRINT,'WARNING: Calling CRISPEX with DT has no influence when no SPCUBE is supplied. Setting seconds per timestep to default value.'
			dt_set = 0
			dt = 1.
			spytitle = 'Frame number'
		ENDELSE
	ENDIF ELSE BEGIN
		dt_set = 0
		dt = 1.
		spytitle = 'Frame number'
	ENDELSE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial playback parameters)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial playback parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SPECTRAL PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial spectral parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initial spectral parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	lp_first 	= 0													; Set number of first lineposition
	lp_last		= nlp-1													; Set number of last lineposition
	sp		= nt * nlp												; Set spectral dimension
	lp_start 	= lc
	lp_ref_first = lp_first
	IF showrefls THEN BEGIN
		lp_ref_last = refnlp - 1
		lp_ref_start = reflc
	ENDIF ELSE IF (refnlp GT 1) THEN BEGIN
		lp_ref_last = lp_last
		lp_ref_start = lc
	ENDIF ELSE BEGIN
		lp_ref_last = 1
		lp_ref_start = 0
	ENDELSE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial spectral parameters)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial spectral parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- WINDOW SIZES (CHANGE ONLY
;--------------------------------------------------------------------------------- NUMERICAL VALUES!)
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (window sizes)...",a5)','     ') 
	feedback_text = [feedback_text,'> Window sizes... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	
	heightset 	= 0
	refheightset 	= 0
	IF ((N_ELEMENTS(XTITLE) GE 1) AND (N_ELEMENTS(XTITLE) LE 2)) THEN BEGIN
		IF (STRCOMPRESS(xtitle[0]) NE '') THEN BEGIN
			heightset = (STRCMP(xtitle[0],'Height',6,/FOLD_CASE) OR $
                   STRCMP(xtitle[0],'z',1,/FOLD_CASE))
			v_dop_set = 0
			spxtitle = xtitle[0]
		ENDIF
		IF (N_ELEMENTS(XTITLE) EQ 2) THEN BEGIN
			IF (STRCOMPRESS(xtitle[1]) NE '') THEN BEGIN
				refheightset = (STRCMP(xtitle[1],'Height',6,/FOLD_CASE) OR $
                        STRCMP(xtitle[1],'z',1,/FOLD_CASE))
				v_dop_set_ref = 0
				refspxtitle = xtitle[1]
			ENDIF
		ENDIF
	ENDIF
	IF heightset THEN detspect_scale = 0
	IF refheightset THEN ref_detspect_scale = 0
;  wav_h_unit = STRARR(1+showref)
;  IF ((STRPOS(spxtitle,'['))[0] NE -1) THEN $
;    wav_h_unit[0] = STRMID(spxtitle,STRPOS(spxtitle,'['),STRLEN(spxtitle))
;  IF showref THEN BEGIN
;    IF ((STRPOS(spxtitle,'['))[0] NE -1) THEN $
;      wav_h_unit[1] = STRMID(spxtitle,STRPOS(spxtitle,'['),STRLEN(spxtitle))
;  ENDIF

	sp_h 	= ['Spectral','Height']
	wav_h 	= ['Wavelength','Height']
;	lp_h		= ['lp','h']
	lp_h_capital	= ['S','H']
	but_tooltip = ['Spectrum','Height distribution']

	lsytitle	= 'Intensity'
	reflsytitle	= 'Intensity'
	IF ((N_ELEMENTS(YTITLE) GE 1) AND (N_ELEMENTS(YTITLE) LE 2)) THEN BEGIN
		IF (STRCOMPRESS(ytitle[0]) NE '') THEN lsytitle = ytitle[0]
		IF (N_ELEMENTS(YTITLE) EQ 2) THEN BEGIN
			IF (STRCOMPRESS(ytitle[1]) NE '') THEN reflsytitle = ytitle[1]
		ENDIF
	ENDIF
	IF (lsytitle NE 'Intensity') THEN detspect_scale = 0
	IF (reflsytitle NE 'Intensity') THEN ref_detspect_scale = 0

	xdelta		= 20													; Extra xoffset for positioning of windows
	ydelta		= 40													; Extra yoffset for positioning of windows
	refxoffset	= 0
	refyoffset	= ydelta

	ratio 		= FLOAT(nx) / FLOAT(ny)											; Determine x/y-ratio
	x_scr_size = screensizes[2,monitor_order[0]]
	y_scr_size = screensizes[3,monitor_order[0]]
	IF ((x_scr_size GT (nx+2*xdelta+0.45*x_scr_size)) AND (y_scr_size GT (ny+ydelta+90))) THEN BEGIN		; If xsize > nx+space for spectral windows AND ysize > ny+space for params window, then:
		IF ((nx LT 0.48 * x_scr_size) AND (ny LT (0.48 * x_scr_size / ratio)) AND KEYWORD_SET(WINDOW_LARGE)) THEN BEGIN				; If xsize is small, then still go to old settings procedures
			imwinx 	= 0.48 * x_scr_size											; Set maximum x-extent of image window
			imwiny 	= imwinx / ratio											; Set maximum y-extent of image window
			IF (verbosity[1] EQ 1) THEN BEGIN
				PRINT,''
				PRINT,'CRISPEX SETUP: User screen resolution allows 1:1 image window sizing, but dimensions are small. Reverting. Image window is '+STRTRIM(imwinx,2)+'x'+STRTRIM(imwiny,2)+'.'
			ENDIF
		ENDIF ELSE BEGIN
			imwinx	= nx													; - use actual nx as imwinx
			imwiny	= ny													; - use actual ny as imwiny
			IF (verbosity[1] EQ 1) THEN BEGIN
				PRINT,''
				PRINT,'CRISPEX SETUP: User screen resolution allows 1:1 image window sizing. Image window is '+STRTRIM(imwinx,2)+'x'+STRTRIM(imwiny,2)+'.'
			ENDIF
		ENDELSE
	ENDIF ELSE BEGIN													; Else use the old procedures to determine imwinx and imwiny
		imwinx 	= 0.48 * x_scr_size											; Set maximum x-extent of image window
		imwiny 	= imwinx / ratio											; Set maximum y-extent of image window
		IF (imwiny GT y_scr_size) THEN BEGIN										; Failsafe to avoid a window larger than the 
			imwiny = 0.85 * y_scr_size											; screensize
			imwinx = imwiny * ratio
		ENDIF
		IF (verbosity[1] EQ 1) THEN BEGIN 
			PRINT,'' 
			PRINT,'CRISPEX SETUP: User screen resolution does not allow for 1:1 image window sizing. Image window is '+STRTRIM(imwinx,2)+'x'+STRTRIM(imwiny,2)+'.' 
		ENDIF
	ENDELSE

	windowx		= 0.2 * x_scr_size											; Set maximum x-extent of spectral win
	IF (nt GE 50) THEN windowy = imwiny ELSE windowy = imwiny/2.
	lswinx 		= 0.25 * x_scr_size											; Set maximum x-extent of loc spec win
	
	spxoffset = imwinx + xdelta
	lsxoffset = spxoffset + windowx + xdelta

	xswinx		= windowx												; Set maximum x-extent of x-slice window
	xswiny		= windowy												; Set maximum y-extent of x-slice window

	lswintitle	= ['Detailed spectrum','Height distribution']
	lsmargin 	= 0.1
	lswall 		= 0.03
	ticklen 	= 0.01
	xsize 		= 1.*lswinx

	IF (ns LE 2) THEN BEGIN
		npanels = ns	&	cols = ns	&	rowarr = REPLICATE(0,ns)
	ENDIF ELSE BEGIN
		npanels = 4	&	cols = 2	&	rowarr = [1,1,0,0]
	ENDELSE
	rows = CEIL(npanels / FLOAT(cols))
	lsx0 = FLTARR(npanels)
	lsx1 = FLTARR(npanels)
	lsy0 = FLTARR(npanels)
	lsy1 = FLTARR(npanels)
	lswidth = (xsize/lswinx - (cols*lsmargin + lswall))/FLOAT(cols)
	lsheight = lswidth * 2D / (1 + SQRT(5))
	IF (v_dop_set EQ 1) THEN lswiny = (lsmargin + rows*lsheight + rows*lsmargin + (rows-1)*lswall) * lswinx ELSE lswiny = (lswall + rows*lsheight + rows*lsmargin) * lswinx
	lsx0 		= lsmargin * lswinx/lswinx + (INDGEN(npanels) MOD cols) * (lswidth + lsmargin) * lswinx/lswinx
	lsx1 		= lsx0 + lswidth * lswinx/lswinx
	lsy0 		= lsmargin * lswinx/lswiny + rowarr * (lsheight + lsmargin + v_dop_set*lswall) * lswinx/lswiny
	lsy1 		= lsy0 + lsheight * lswinx/lswiny
	lsxmargin_init	= lsmargin * lswinx
	lsxwall_init	= lswall * lswinx
	lsxticklen 	= ticklen / lsheight
	lsyticklen 	= ticklen / lswidth

	reflswintitle	= ['Reference detailed spectrum','Reference height distribution']
	reflswidth 	= (xsize/lswinx - (lsmargin + lswall))
	reflsheight 	= reflswidth * 2D / (1 + SQRT(5))
	reflswinx	= lswinx
	IF (v_dop_set_ref EQ 1) THEN reflswiny = (lsmargin + reflsheight + lsmargin) * lswinx ELSE reflswiny = (lsmargin + reflsheight + lswall) * lswinx
	reflsx0 	= lsmargin * reflswinx/reflswinx 
	reflsx1 	= reflsx0 + reflswidth * reflswinx/reflswinx
	reflsy0 	= lsmargin * reflswinx/reflswiny
	reflsy1 	= reflsy0 + reflsheight * reflswinx/reflswiny
	reflsxmargin_init= lsmargin * reflswinx
	reflsxwall_init	= lswall * reflswinx
	reflsxticklen 	= ticklen / reflsheight
	reflsyticklen 	= ticklen / reflswidth

	intwidth 	= (xsize/lswinx - (lsmargin + lswall))
	intheight 	= intwidth * 2D / (1 + SQRT(5))
	intwinx		= lswinx
	intwiny 	= (lsmargin + intheight + lswall) * intwinx
	intx0 		= lsmargin * intwinx/intwinx 
	intx1 		= intx0 + intwidth * intwinx/intwinx
	inty0 		= lsmargin * intwinx/intwiny
	inty1	 	= inty0 + intheight * intwinx/intwiny
	intxmargin_init	= lsmargin * intwinx
	intxwall_init	= lswall * intwinx
	intxticklen 	= ticklen / intheight
	intyticklen 	= ticklen / intwidth

	spwintitle	= ['Spectral T-slice','Height T-slice']
	spxmargin_init	= lsmargin * lswinx * 1.3
	spxwall_init	= lswall * lswinx

	spwinx 		= windowx
	spmargin 	= spxmargin_init/spwinx
	spwall 		= spxwall_init/spwinx
	xsize 		= 1.*spwinx

	phiswinx	= spwinx
	phiswiny	= imwiny

	spwidth 	= (1. - (spmargin + spwall))
	IF spfile THEN spwiny = windowy ELSE spwiny = imwiny
	IF ((v_dop_set EQ 1) OR (ns GT 1)) THEN BEGIN
		spheight = (1. - (spmargin * 2.) * spwinx/spwiny)
		phisheight = (1. - (spmargin * 2.) * phiswinx/phiswiny)
	ENDIF ELSE BEGIN
		spheight = (1. - (spmargin + spwall) * spwinx/spwiny)
		phisheight = (1. - (spmargin + spwall) * phiswinx/phiswiny)
	ENDELSE
	IF (v_dop_set_ref EQ 1) THEN refspheight = (1. - (spmargin * 2.) * spwinx/spwiny) ELSE refspheight = (1. - (spmargin + spwall) * spwinx/spwiny)
	refspwiny	= windowy

	spx0 		= spmargin * spwinx/spwinx 
	spx1 		= spx0 + spwidth * spwinx/spwinx
	spy0 		= spmargin * spwinx/spwiny
	spy1 		= spy0 + spheight; * spwinx/spwiny
	xplspw		= spx1 - spx0												; x-extent of the plot
	yplspw		= spy1 - spy0												; y-extent of the plot

	refspwintitle	= ['Reference spectral T-slice','Reference height T-slice']
	refspy0 	= spmargin * spwinx/refspwiny
	refspy1 	= refspy0 + refspheight; * spwinx/refspwiny
	refyplspw	= refspy1 - refspy0												; y-extent of the plot
	
	IF ((spfile EQ 1) OR (SINGLE_CUBE GE 1)) THEN ntreb = yplspw * spwiny ELSE ntreb = 0						; actual nt rebinning factor
	refntreb	= refyplspw * refspwiny
	nlpreb		= xplspw * spwinx							 				; actual nlp rebinning factor
	spxticklen 	= -1. * ticklen / spheight
	spyticklen 	= -1. * ticklen / spwidth

	phiswintitle	= ['Spectral Phi-slice','Height Phi-slice']
	phisx0 		= spmargin * phiswinx/phiswinx 
	phisx1 		= phisx0 + spwidth * phiswinx/phiswinx
	phisy0 		= spmargin * phiswinx/phiswiny
	phisy1 		= phisy0 + phisheight; * spwinx/spwiny
	phisxplspw	= phisx1 - phisx0												; x-extent of the plot
	phisyplspw	= phisy1 - phisy0												; y-extent of the plot
	nphireb		= phisyplspw * phiswiny
	phisxticklen 	= -1. * ticklen / phisheight
	phisyticklen 	= -1. * ticklen / spwidth
	
	loopheight	= (1. - (spmargin + spwall) * spwinx/windowy)
	loopy1		= spy0 + loopheight
	loopyplspw	= loopy1 - spy0
	loopntreb	= loopyplspw * windowy

	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (window sizes)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Window sizes... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SPATIAL PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial spatial parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initial spatial parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	x_first		= 0													; Set number of first x-coordinate
	x_last		= nx-1													; Set number of last x-coordinate
	y_first		= 0													; Set number of first y-coordinate
	y_last		= ny-1													; Set number of last y-coordinate
	x_start		= FLOAT(FLOOR(nx/2))												; Determine the middle x-coordinate
	sx_start	= x_start * imwinx / FLOAT(nx)										; Convert that to device
	y_start		= FLOAT(FLOOR(ny/2))												; Determine the middle y-coordinate
	sy_start	= y_start * imwiny / FLOAT(ny)										; Convert that to device
;	arcsecpix	= 0.0592
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (initial spatial parameters)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial spatial parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- SCALING PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (image scaling parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initial scaling parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	imagescale = PTR_NEW([0,0,0])												; Image scaling based on first image
	relative_scaling = PTR_NEW([0,0,0])
	immin = DBLARR(nlp,ns)
	immax = DBLARR(nlp,ns)
	immean = DBLARR(nlp,ns)
	imsdev = DBLARR(nlp,ns)
	dopplermin = DBLARR(nlp,ns)
	dopplermax = DBLARR(nlp,ns)
	ls_low_y = FLTARR(ns)
	ls_upp_y = FLTARR(ns)
	ls_yrange = FLTARR(ns)
	int_low_y = FLTARR(ns)
	int_upp_y = FLTARR(ns)
	FOR j=0,ns-1 DO BEGIN
		FOR k=0,lp_last DO BEGIN
			temp_image = imagefile[j*nlp + k]
			immin[k,j] = MIN(temp_image, MAX=max_val)
			immax[k,j] = max_val
			immean[k,j] = MEAN(temp_image)
			imsdev[k,j] = STDDEV(temp_image)
			temp_k = 2*lc - k
			IF ((temp_k EQ lc) OR (temp_k LT 0) OR (temp_k GT lp_last)) THEN BEGIN
				dopplermin[k,j] = 0
				dopplermax[k,j] = 0
			ENDIF ELSE BEGIN
				mirror_temp_image = imagefile[j*nlp + temp_k]
				IF (temp_k GT lc) THEN BEGIN
					dopplermin[k,j] = MIN(temp_image - mirror_temp_image, MAX=max_val)
					dopplermax[k,j] = max_val
				ENDIF ELSE BEGIN
					dopplermin[k,j] = MIN(mirror_temp_image - temp_image, MAX=max_val)
					dopplermax[k,j] = max_val
				ENDELSE
			ENDELSE
		ENDFOR
		IF scalestokes THEN BEGIN
			ls_low_y[j] = MIN(immin[*,j])/ms
			ls_upp_y[j] = MAX(immax[*,j])/ms
			ls_low_y[j] = MIN((immean[*,j]-3.*imsdev[*,j])/ms)
			ls_upp_y[j] = MAX((immean[*,j]+3.*imsdev[*,j])/ms)
		ENDIF ELSE BEGIN
			ls_low_y[j] = MIN(immin[*,j])/ms[j]
			ls_upp_y[j] = MAX(immax[*,j])/ms[j]
			ls_low_y[j] = MIN((immean[*,j]-3.*imsdev[*,j])/ms[j])
			ls_upp_y[j] = MAX((immean[*,j]+3.*imsdev[*,j])/ms[j])
		ENDELSE
		ls_yrange[j] = ls_upp_y[j] - ls_low_y[j]
		max_imsdev = MAX(imsdev[*,j])
		int_low_y[j] = (MEAN(immean[*,j])-3.*max_imsdev)/ABS(MEAN(immean[*,j]))
		int_upp_y[j] = (MEAN(immean[*,j])+3.*max_imsdev)/ABS(MEAN(immean[*,j]))
	ENDFOR
	ls_low_y_init = ls_low_y[0]
	ls_upp_y_init = ls_upp_y[0]
	ls_low_y = PTR_NEW(ls_low_y,/NO_COPY)
	ls_upp_y = PTR_NEW(ls_upp_y,/NO_COPY)
	ls_yrange = PTR_NEW(ls_yrange,/NO_COPY)
	int_low_y = PTR_NEW(int_low_y,/NO_COPY)
	int_upp_y = PTR_NEW(int_upp_y,/NO_COPY)

	ls_low_y_ref = 0
	ls_upp_y_ref = 0
	IF showref THEN BEGIN
		IF (refnt EQ 0) THEN BEGIN
			refmin = MIN(referencefile, MAX=refmax)
			refmean = MEAN(referencefile)
			refdev = STDDEV(referencefile)
		ENDIF ELSE IF (((refnt EQ 1) OR (refnt EQ nt)) AND (refnlp EQ 1)) THEN BEGIN
			refmin = MIN(referencefile[0], MAX=refmax)
			refmean = MEAN(referencefile[0])
			refdev = STDDEV(referencefile[0])
		ENDIF ELSE BEGIN
			refmin = FLTARR(refnlp)
			refmax = FLTARR(refnlp)
			refmean = FLTARR(refnlp)
			refdev = FLTARR(refnlp)
			FOR k=0,refnlp-1 DO BEGIN
				temp_referencefile = referencefile[k]
				refmin[k] = MIN(temp_referencefile, MAX=max_val)
				refmax[k] = max_val
				refmean[k] = MEAN(temp_referencefile)
				refdev[k] = STDDEV(temp_referencefile)
			ENDFOR
		ENDELSE
		IF showrefls THEN BEGIN
			ls_low_y_ref = MIN((refmean-3.*refdev)/refms)
			ls_upp_y_ref = MAX((refmean+3.*refdev)/refms)
		ENDIF
	ENDIF ELSE BEGIN
		refmin = 0
		refmax = 0
	ENDELSE
	ls_yrange_ref = ls_upp_y_ref - ls_low_y_ref 

	scale_range = [0.,0.]
	scale_minimum = [0.,0.]
	scale_max_val = [255.,255.]
	scale_min_val = [0., 0.]
	rel_scale_max_val = [100., 100.]
	rel_scale_min_val = [0.,0.]
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (image scaling parameters)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial scaling parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- OTHER PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (other parameters)...",a5)','     ') 
	feedback_text = [feedback_text,'> Other parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	CD, CURRENT=curpath
	default_prefipath = STRTRIM(curpath+PATH_SEP(),2)
	default_prefopath = default_prefipath
	IF (defipath EQ 0) THEN BEGIN
		prefipath = default_prefipath
		ipath = default_prefipath
	ENDIF ELSE ipath = prefipath
	IF (defopath EQ 0) THEN BEGIN
		prefopath = default_prefopath
		opath = default_prefopath
	ENDIF ELSE opath = prefopath
	opath_write = FILE_TEST(opath, /WRITE)
	IF KEYWORD_SET(EXTS) THEN exts_set = 1 ELSE exts_set = 0
	IF (SINGLE_CUBE GE 1) THEN BEGIN
		IF (nt GT 1) THEN BEGIN
			exts_set = 1						; Automatically set EXTS when providing a (single lineposition) 3D cube
			PRINT,'WARNING: The exact timeslice (EXTS) keyword has been automatically set to enable the '
			PRINT,'		drawing of loop paths and extraction of timeslices!'
		ENDIF ELSE exts_set=0
	ENDIF
	refexts_set = (refspfile NE 1)
	lp_slid_sens = (nlp GE 2)
	lp_blink_vals_sens = (nlp GT 2)
	lp_last_slid = (nlp-1) > 1
	lp_last_blink = (nlp-1) > 2
	lp_last_vals = nlp-1

	scale_cubes_vals = [1.,1.]
	IF (N_ELEMENTS(SCALE_CUBES) EQ 1) THEN scale_cubes_vals[0] = scale_cubes ELSE IF (N_ELEMENTS(SCALE_CUBES) EQ 2) THEN scale_cubes_vals = scale_cubes
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting start-up options (other parameters)...",a5)','done!') 
		IF (verbosity[1] EQ 1) THEN BEGIN
			PRINT,'' 
			PRINT,'CRISPEX SETUP: Input path set to: '+STRTRIM(ipath,2)
			PRINT,'CRISPEX SETUP: Output path set to: '+STRTRIM(opath,2)
		ENDIF
	ENDIF
	feedback_text = [feedback_text[0]+'done!',feedback_text[1:N_ELEMENTS(feedback_text)-2],'> Other parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	WAIT,0.1

;================================================================================= SETTING UP WIDGET
;--------------------------------------------------------------------------------- INITIALISE CONTROL PANEL
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (loading BMP buttons)...",a5)','     ') 
	feedback_text = ['Setting up widget... ','> Loading BMP buttons... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	bmpbut_search = FILE_SEARCH(dir_buttons, '*.bmp', COUNT=bmpbut_count)
	IF (bmpbut_count EQ 16) THEN BEGIN
		bmpbut_fbwd_idle 	= CRISPEX_READ_BMP_BUTTONS('fbwd_idle.bmp',dir_buttons)
		bmpbut_fbwd_pressed 	= CRISPEX_READ_BMP_BUTTONS('fbwd_pressed.bmp',dir_buttons)
		bmpbut_bwd_idle 	= CRISPEX_READ_BMP_BUTTONS('bwd_idle.bmp',dir_buttons)
		bmpbut_bwd_pressed 	= CRISPEX_READ_BMP_BUTTONS('bwd_pressed.bmp',dir_buttons)
		bmpbut_pause_idle	= CRISPEX_READ_BMP_BUTTONS('pause_idle.bmp',dir_buttons)
		bmpbut_pause_pressed 	= CRISPEX_READ_BMP_BUTTONS('pause_pressed.bmp',dir_buttons)
		bmpbut_fwd_idle 	= CRISPEX_READ_BMP_BUTTONS('fwd_idle.bmp',dir_buttons)
		bmpbut_fwd_pressed 	= CRISPEX_READ_BMP_BUTTONS('fwd_pressed.bmp',dir_buttons)
		bmpbut_ffwd_idle 	= CRISPEX_READ_BMP_BUTTONS('ffwd_idle.bmp',dir_buttons)
		bmpbut_ffwd_pressed 	= CRISPEX_READ_BMP_BUTTONS('ffwd_pressed.bmp',dir_buttons)
		bmpbut_loop_idle	= CRISPEX_READ_BMP_BUTTONS('loop_idle.bmp',dir_buttons)
		bmpbut_loop_pressed	= CRISPEX_READ_BMP_BUTTONS('loop_pressed.bmp',dir_buttons)
		bmpbut_cycle_idle	= CRISPEX_READ_BMP_BUTTONS('cycle_idle.bmp',dir_buttons)
		bmpbut_cycle_pressed	= CRISPEX_READ_BMP_BUTTONS('cycle_pressed.bmp',dir_buttons)
		bmpbut_blink_idle	= CRISPEX_READ_BMP_BUTTONS('blink_idle.bmp',dir_buttons)
		bmpbut_blink_pressed	= CRISPEX_READ_BMP_BUTTONS('blink_pressed.bmp',dir_buttons)
		donetext = 'done!'
	ENDIF ELSE BEGIN
		bmpbut_fbwd_idle 	= '<<'		&	bmpbut_fbwd_pressed 	= bmpbut_fbwd_idle
		bmpbut_bwd_idle 	= '<'		&	bmpbut_bwd_pressed 	= bmpbut_bwd_idle
		bmpbut_pause_idle	= '||'		&	bmpbut_pause_pressed 	= bmpbut_pause_idle
		bmpbut_fwd_idle 	= '>'		&	bmpbut_fwd_pressed 	= bmpbut_fwd_idle
		bmpbut_ffwd_idle 	= '>>'		&	bmpbut_ffwd_pressed 	= bmpbut_ffwd_idle
		bmpbut_loop_idle	= 'Loop'	&	bmpbut_loop_pressed	= bmpbut_loop_idle
		bmpbut_cycle_idle	= 'Cycle'	&	bmpbut_cycle_pressed	= bmpbut_cycle_idle
		bmpbut_blink_idle	= 'Blink'	&	bmpbut_blink_pressed	= bmpbut_blink_idle
		donetext = 'failed.'
	ENDELSE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (loading BMP buttons)...",a5)',donetext) 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Loading BMP buttons... '+donetext]
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIALISE CONTROL PANEL
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (initializing control panel)...",a5)','     ') 
	feedback_text = [feedback_text,'> Initializing control panel... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	control_panel		= WIDGET_BASE(TITLE = 'CRISPEX'+instance_label+': Control Panel', TLB_FRAME_ATTR = 1, /COLUMN, /FRAME, KILL_NOTIFY = 'CRISPEX_CLOSE_CLEANUP', MBAR = menubar, TLB_SIZE_EVENTS = 1)
	filemenu		= WIDGET_BUTTON(menubar, VALUE = 'File', /MENU, UVALUE = 'file')
	about			= WIDGET_BUTTON(filemenu, VALUE = 'About', EVENT_PRO = 'CRISPEX_ABOUT_WINDOW', ACCELERATOR = 'Ctrl+A')
	preferences		= WIDGET_BUTTON(filemenu, VALUE = 'Preferences', EVENT_PRO = 'CRISPEX_PREFERENCES_WINDOW', ACCELERATOR = 'Ctrl+P')
	dispwid			= WIDGET_BUTTON(filemenu, VALUE = 'Display window IDs', EVENT_PRO = 'CRISPEX_DISPWIDS', ACCELERATOR = 'Ctrl+I', /CHECKED_MENU)
	pathmenu		= WIDGET_BUTTON(filemenu, VALUE = 'Set path...', /SEPARATOR, /MENU)
	setipath		= WIDGET_BUTTON(pathmenu, VALUE = 'Set input path', EVENT_PRO = 'CRISPEX_SAVE_SET_IPATH')
	setopath		= WIDGET_BUTTON(pathmenu, VALUE = 'Set output path', EVENT_PRO = 'CRISPEX_SAVE_SET_OPATH')
	savemenu		= WIDGET_BUTTON(filemenu, VALUE = 'Save current...', /MENU)
	timeslicemenu		= WIDGET_BUTTON(savemenu, VALUE = 'Time slice(s)...', /MENU, SENSITIVE = 0)
	save_loop_pts		= WIDGET_BUTTON(savemenu, VALUE = 'Loop for later retrieval', EVENT_PRO = 'CRISPEX_SAVE_LOOP_PTS', SENSITIVE = 0)
	approxmenu		= WIDGET_BUTTON(timeslicemenu, VALUE = 'Approximated loop...', /MENU)
	save_app_slab_but	= WIDGET_BUTTON(approxmenu, VALUE = 'All '+STRLOWCASE(sp_h[heightset])+' positions', EVENT_PRO = 'CRISPEX_SAVE_APPROX_LOOPSLAB')
	save_app_slice_but	= WIDGET_BUTTON(approxmenu, VALUE = 'Current '+STRLOWCASE(sp_h[heightset])+' position', EVENT_PRO = 'CRISPEX_SAVE_APPROX_LOOPSLICE')
	interpolmenu		= WIDGET_BUTTON(timeslicemenu, VALUE = 'Interpolated loop...', /MENU)
	save_ex_slab_but	= WIDGET_BUTTON(interpolmenu, VALUE = 'All '+STRLOWCASE(sp_h[heightset])+' positions', EVENT_PRO = 'CRISPEX_SAVE_EXACT_LOOPSLAB_CHECK')
	save_ex_slice_but	= WIDGET_BUTTON(interpolmenu, VALUE = 'Current '+STRLOWCASE(sp_h[heightset])+' position', EVENT_PRO = 'CRISPEX_SAVE_EXACT_LOOPSLICE')
	retrievemenu		= WIDGET_BUTTON(filemenu, VALUE = 'Retrieve and save...', /MENU)
	sel_saved_loop		= WIDGET_BUTTON(retrievemenu, VALUE = 'From selected saved loop(s)', SENSITIVE = 0, EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_MENU')
	all_saved_loop		= WIDGET_BUTTON(retrievemenu, VALUE = 'From all saved loops...', /MENU, SENSITIVE = 0)
	all_saved_all_pos	= WIDGET_BUTTON(all_saved_loop, VALUE = 'At all '+STRLOWCASE(sp_h[heightset])+' positions', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLAB')
	all_saved_sel_pos	= WIDGET_BUTTON(all_saved_loop, VALUE = 'At saved '+STRLOWCASE(sp_h[heightset])+' position', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLICE')
	det_file_loop		= WIDGET_BUTTON(retrievemenu, VALUE = 'From detection file', EVENT_PRO = 'CRISPEX_RETRIEVE_DET_FILE_MENU')
	save_as_menu		= WIDGET_BUTTON(filemenu, VALUE = 'Save as...', /MENU)
	save_as_png_menu	= WIDGET_BUTTON(save_as_menu, VALUE = 'PNG...', /MENU)
	save_as_png_sns		= WIDGET_BUTTON(save_as_png_menu, VALUE = 'Snapshot', EVENT_PRO = 'CRISPEX_SAVE_PNG_SNAPSHOT')
	save_as_png_all		= WIDGET_BUTTON(save_as_png_menu, VALUE = 'All frames', EVENT_PRO = 'CRISPEX_SAVE_PNG_ALL_FRAMES', SENSITIVE = (nt GT 1))
	save_as_png_linescan	= WIDGET_BUTTON(save_as_png_menu, VALUE = 'Line scan', EVENT_PRO = 'CRISPEX_SAVE_PNG_LINESCAN', SENSITIVE = (nlp GT 1))
	save_as_jpg_menu	= WIDGET_BUTTON(save_as_menu, VALUE = 'JPEG...', /MENU)
	save_as_jpg_sns 	= WIDGET_BUTTON(save_as_jpg_menu, VALUE = 'Snapshot', EVENT_PRO = 'CRISPEX_SAVE_JPEG_SNAPSHOT')
	save_as_jpg_all		= WIDGET_BUTTON(save_as_jpg_menu, VALUE = 'All frames', EVENT_PRO = 'CRISPEX_SAVE_JPEG_ALL_FRAMES', SENSITIVE = (nt GT 1))
	save_as_jpg_linescan 	= WIDGET_BUTTON(save_as_jpg_menu, VALUE = 'Linescan', EVENT_PRO = 'CRISPEX_SAVE_JPEG_LINESCAN', SENSITIVE = (nlp GT 1))
	save_as_mpeg		= WIDGET_BUTTON(save_as_menu, VALUE = 'MPEG', EVENT_PRO = 'CRISPEX_SAVE_MPEG', SENSITIVE = (nt GT 1))
	session			= WIDGET_BUTTON(filemenu, VALUE = 'Session...', /MENU)
	save_session		= WIDGET_BUTTON(session, VALUE = 'Save current', EVENT_PRO = 'CRISPEX_SESSION_SAVE_WINDOW')
	restore_session		= WIDGET_BUTTON(session, VALUE = 'Restore other', EVENT_PRO = 'CRISPEX_SESSION_RESTORE_WINDOW')
	clear_menu		= WIDGET_BUTTON(filemenu, VALUE = 'Clear...', /MENU)
	clear_current_estimate	= WIDGET_BUTTON(clear_menu, VALUE = 'Current time estimate', EVENT_PRO = 'CRISPEX_CLEAR_CURRENT_ESTIMATE', SENSITIVE = estimate_run)
	clear_current_cpft	= WIDGET_BUTTON(clear_menu, VALUE = 'CPFT file for current machine', EVENT_PRO = 'CRISPEX_CLEAR_CURRENT_CPFT', SENSITIVE = (cpftfilecount EQ 1))
	clear_current_inst	= WIDGET_BUTTON(clear_menu, VALUE = 'Instance file for current machine', EVENT_PRO = 'CRISPEX_CLEAR_CURRENT_INST', SENSITIVE = (instfilecount EQ 1))
	focus_session_windows	= WIDGET_BUTTON(filemenu, VALUE = 'Bring all to front', EVENT_PRO = 'CRISPEX_DISPLAYS_ALL_TO_FRONT', ACCELERATOR = 'Ctrl+F')
	helpmenu		= WIDGET_BUTTON(filemenu, VALUE = 'Help', EVENT_PRO = 'CRISPEX_HELP', /SEPARATOR, ACCELERATOR = 'Ctrl+H')
	mailmenu		= WIDGET_BUTTON(filemenu, VALUE = 'E-mail...', /MENU)
	mailbugmenu		= WIDGET_BUTTON(mailmenu, VALUE = 'A bug report', EVENT_PRO = 'CRISPEX_MAIL_BUG')
	mailsugmenu		= WIDGET_BUTTON(mailmenu, VALUE = 'A suggestion', EVENT_PRO = 'CRISPEX_MAIL_SUGGESTION')
	exitmenu		= WIDGET_BUTTON(filemenu, VALUE = 'Quit', EVENT_PRO = 'CRISPEX_CLOSE', /SEPARATOR, ACCELERATOR = 'Ctrl+Q')

	shortcutmenu		= WIDGET_BUTTON(menubar, VALUE = 'Control shortcuts', /MENU, UVALUE = 'shortcut')
	sh_playbackmenu 	= WIDGET_BUTTON(shortcutmenu, VALUE = 'Playback options', /MENU, UVALUE = 'playback', SENSITIVE = (nt GT 1))
	sh_fbwd_button		= WIDGET_BUTTON(sh_playbackmenu, VALUE = 'One frame backward', EVENT_PRO = 'CRISPEX_PB_FASTBACKWARD', ACCELERATOR = 'Shift+B')
	sh_backward_button	= WIDGET_BUTTON(sh_playbackmenu, VALUE = 'Play backwards', EVENT_PRO = 'CRISPEX_PB_BACKWARD', ACCELERATOR = 'Shift+Backspace')
	sh_pause_button		= WIDGET_BUTTON(sh_playbackmenu, VALUE = 'Pause', EVENT_PRO = 'CRISPEX_PB_PAUSE', ACCELERATOR = 'Shift+Space')
	sh_forward_button	= WIDGET_BUTTON(sh_playbackmenu, VALUE = 'Play forwards', EVENT_PRO = 'CRISPEX_PB_FORWARD', ACCELERATOR = 'Shift+Tab')
	sh_ffwd_button		= WIDGET_BUTTON(sh_playbackmenu, VALUE = 'One frame forward', EVENT_PRO = 'CRISPEX_PB_FASTFORWARD', ACCELERATOR = 'Shift+F')
	sh_spectralmenu 	= WIDGET_BUTTON(shortcutmenu, VALUE = sp_h[heightset]+' options', /MENU, UVALUE = 'spectral', SENSITIVE = (nlp GT 1))
	sh_lp_incr_button 	= WIDGET_BUTTON(sh_spectralmenu, VALUE = sp_h[heightset]+' position +', EVENT_PRO = 'CRISPEX_SLIDER_LP_INCR', ACCELERATOR = 'Shift+S')
	sh_lp_decr_button 	= WIDGET_BUTTON(sh_spectralmenu, VALUE = sp_h[heightset]+' position -', EVENT_PRO = 'CRISPEX_SLIDER_LP_DECR', ACCELERATOR = 'Shift+A')
	sh_zoommenu		= WIDGET_BUTTON(shortcutmenu, VALUE = 'Zoom options',/MENU)
	sh_zoom_one		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 1x', EVENT_PRO = 'CRISPEX_ZOOMFAC_ONE',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+1')
	sh_zoom_two		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 2x', EVENT_PRO = 'CRISPEX_ZOOMFAC_TWO',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+2')
	sh_zoom_three		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 3x', EVENT_PRO = 'CRISPEX_ZOOMFAC_THREE',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+3')
	sh_zoom_four		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 4x', EVENT_PRO = 'CRISPEX_ZOOMFAC_FOUR',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+4')
	sh_zoom_six		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 6x', EVENT_PRO = 'CRISPEX_ZOOMFAC_SIX',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+6')
	sh_zoom_eight		= WIDGET_BUTTON(sh_zoommenu, VALUE = 'Zoom 8x', EVENT_PRO = 'CRISPEX_ZOOMFAC_EIGHT',/NO_RELEASE, ACCELERATOR = 'Ctrl+Shift+8')
	sh_runtime_menu		= WIDGET_BUTTON(shortcutmenu, VALUE = 'Runtime options', /MENU, UVALUE = 'runtime')
	sh_runtime_interrupt	= WIDGET_BUTTON(sh_runtime_menu, VALUE = 'Interrupt', EVENT_PRO = 'CRISPEX_INTERRUPT', ACCELERATOR = 'Ctrl+Shift+C')
	sh_runtime_verb_menu	= WIDGET_BUTTON(sh_runtime_menu, VALUE = 'Verbosity', /MENU, UVALUE = 'verbosity')
	sh_verb_0		= WIDGET_BUTTON(sh_runtime_verb_menu, VALUE = 'No verbosity', EVENT_PRO = 'CRISPEX_VERBOSE_SET', UVALUE = -1, /NO_RELEASE, /CHECKED_MENU, ACCELERATOR = 'Shift+0')
	sh_verb_4		= WIDGET_BUTTON(sh_runtime_verb_menu, VALUE = 'Basic runtime verbosity', EVENT_PRO = 'CRISPEX_VERBOSE_SET', UVALUE = 2, /NO_RELEASE, /CHECKED_MENU, ACCELERATOR = 'Shift+4')
	sh_verb_8		= WIDGET_BUTTON(sh_runtime_verb_menu, VALUE = 'Extended runtime verbosity', EVENT_PRO = 'CRISPEX_VERBOSE_SET', UVALUE = 3, /NO_RELEASE, /CHECKED_MENU, ACCELERATOR = 'Shift+8')
	sh_verb_16		= WIDGET_BUTTON(sh_runtime_verb_menu, VALUE = 'Enable playback statistics', EVENT_PRO = 'CRISPEX_VERBOSE_SET', UVALUE = 4, /NO_RELEASE, /CHECKED_MENU, ACCELERATOR = 'Shift+P')

	tab_tlb			= WIDGET_TAB(control_panel,LOCATION=0, MULTILINE=5)

	playback_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Temporal', /COLUMN)
	playback_contr		= WIDGET_BASE(control_panel, /ROW)
	playback_field_basic	= WIDGET_BASE(playback_contr, /FRAME, /GRID_LAYOUT, COLUMN=5)
	playback_field_add	= WIDGET_BASE(playback_contr, /FRAME, GRID_LAYOUT = (bmpbut_count EQ 16), COLUMN=3, EXCLUSIVE = (bmpbut_count NE 16))
	fbwd_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_fbwd_idle, EVENT_PRO = 'CRISPEX_PB_FASTBACKWARD', TOOLTIP = 'Move one frame backward')
	backward_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_bwd_idle, EVENT_PRO = 'CRISPEX_PB_BACKWARD', TOOLTIP = 'Play backward')
	pause_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_pause_pressed, EVENT_PRO = 'CRISPEX_PB_PAUSE', TOOLTIP = 'Pause')
	forward_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_fwd_idle, EVENT_PRO = 'CRISPEX_PB_FORWARD', TOOLTIP = 'Play forward')
	ffwd_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_ffwd_idle, EVENT_PRO = 'CRISPEX_PB_FASTFORWARD', TOOLTIP = 'Move one frame forward')
	loop_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_loop_pressed, EVENT_PRO = 'CRISPEX_PB_LOOP', TOOLTIP = 'Loop')
	WIDGET_CONTROL, loop_button, SET_BUTTON = (bmpbut_count NE 16)
	cycle_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_cycle_idle, EVENT_PRO = 'CRISPEX_PB_CYCLE', TOOLTIP = 'Cycle')
	blink_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_blink_idle, EVENT_PRO = 'CRISPEX_PB_BLINK', TOOLTIP = 'Blink')

	t_slid			= WIDGET_SLIDER(playback_tab, TITLE = 'Frame number', MIN = t_first, MAX = t_last, VALUE = t_start, EVENT_PRO = 'CRISPEX_SLIDER_T', /DRAG, SENSITIVE = t_slid_sens)
	t_ranges		= WIDGET_BASE(playback_tab, /COLUMN, /FRAME)
	t_range_field		= WIDGET_BASE(t_ranges, /ROW)
	lower_t_label		= WIDGET_LABEL(t_range_field, VALUE = 'Lower index:', /ALIGN_LEFT)
	lower_t_text		= WIDGET_TEXT(t_range_field, VALUE = STRTRIM(t_first,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_T_LOW', SENSITIVE = t_slid_sens)
	upper_t_label		= WIDGET_LABEL(t_range_field, VALUE = 'Upper index:', /ALIGN_LEFT)
	upper_t_text		= WIDGET_TEXT(t_range_field, VALUE = STRTRIM(t_last_tmp,2),  /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_T_UPP', SENSITIVE = t_slid_sens)
	reset_trange_but	= WIDGET_BUTTON(t_ranges, VALUE = 'Reset temporal boundaries', EVENT_PRO = 'CRISPEX_DISPRANGE_T_RESET', SENSITIVE = 0)
	slice_update_but	= WIDGET_BUTTON(playback_tab, VALUE='Update '+STRLOWCASE(sp_h[heightset])+' slice', EVENT_PRO = 'CRISPEX_UPDATE_SLICES', SENSITIVE = 0)
	t_speed_slid		= WIDGET_SLIDER(playback_tab, TITLE = 'Animation speed [frame/s]', MIN = 1, MAX = 100, VALUE = t_speed, EVENT_PRO = 'CRISPEX_SLIDER_SPEED', /DRAG)
	t_step_slid		= WIDGET_SLIDER(playback_tab, TITLE = 'Frame increment', MIN = 1, MAX = t_last, VALUE = t_step, EVENT_PRO = 'CRISPEX_SLIDER_STEP', SENSITIVE = t_slid_sens)
	imref_blink_field	= WIDGET_BASE(playback_tab, /ROW,/NONEXCLUSIVE)
	imref_blink_but		= WIDGET_BUTTON(imref_blink_field, VALUE = 'Blink between main and reference image', EVENT_PRO = 'CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE', SENSITIVE = showref)

	spectral_tab		= WIDGET_BASE(tab_tlb, TITLE = sp_h[heightset], /COLUMN)
	lp_ranges		= WIDGET_BASE(spectral_tab, /COLUMN, /FRAME)
	lp_range_field		= WIDGET_BASE(lp_ranges, /ROW)
;	lower_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Lower '+lp_h[heightset]+'-value:', /ALIGN_LEFT)
	lower_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Lower index:', /ALIGN_LEFT)
	lower_lp_text		= WIDGET_TEXT(lp_range_field, VALUE = STRTRIM(lp_first,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LP_LOW', SENSITIVE = lp_blink_vals_sens)
;	upper_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Upper '+lp_h[heightset]+'-value:', /ALIGN_LEFT)
	upper_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Upper index:', /ALIGN_LEFT)
	upper_lp_text		= WIDGET_TEXT(lp_range_field, VALUE = STRTRIM(lp_last_vals,2),  /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LP_UPP', SENSITIVE = lp_blink_vals_sens)
	reset_lprange_but	= WIDGET_BUTTON(lp_ranges, VALUE = 'Reset '+STRLOWCASE(sp_h[heightset])+' boundaries', EVENT_PRO = 'CRISPEX_DISPRANGE_LP_RESET', SENSITIVE = 0)
	lp_slid			= WIDGET_SLIDER(control_panel, TITLE = 'Main '+STRLOWCASE(sp_h[heightset])+' position', MIN = lp_first, MAX = lp_last_slid, VALUE = lp_start, EVENT_PRO = 'CRISPEX_SLIDER_LP', $
					/DRAG, SENSITIVE = lp_slid_sens)
	lp_speed_slid		= WIDGET_SLIDER(spectral_tab, TITLE = 'Animation speed [blink/s]', MIN = 1, MAX = 100, VALUE = t_speed, EVENT_PRO = 'CRISPEX_SLIDER_SPEED', /DRAG, SENSITIVE = 0)
	lp_blink_slid		= WIDGET_SLIDER(spectral_tab, TITLE = sp_h[heightset]+' increment', MIN = 1, MAX = lp_last_blink, EVENT_PRO = 'CRISPEX_SLIDER_SPECTSTEP',/DRAG, SENSITIVE = lp_blink_vals_sens)
	lp_blink_field		= WIDGET_BASE(spectral_tab, /ROW,/NONEXCLUSIVE)
	lp_blink_but		= WIDGET_BUTTON(lp_blink_field, VALUE = 'Blink between '+STRLOWCASE(sp_h[heightset])+' positions', EVENT_PRO = 'CRISPEX_PB_SPECTBLINK', SENSITIVE = lp_slid_sens)
	lp_ref_but_field	= WIDGET_BASE(spectral_tab, /ROW, /NONEXCLUSIVE)
	IF (heightset NE refheightset) THEN reflab = STRLOWCASE(sp_h[refheightset])+' ' ELSE reflab = ''
	lp_ref_but		= WIDGET_BUTTON(lp_ref_but_field, VALUE = 'Lock reference '+reflab+'to main '+STRLOWCASE(sp_h[heightset])+' position', EVENT_PRO = 'CRISPEX_SLIDER_LP_REF_LOCK', $
					SENSITIVE = (eqnlps AND (refnlp GT 1)))
	WIDGET_CONTROL, lp_ref_but, SET_BUTTON = (eqnlps AND (refnlp GT 1))
	lp_ref_slid = WIDGET_SLIDER(spectral_tab, TITLE = 'Reference '+STRLOWCASE(sp_h[refheightset])+' position', MIN = lp_ref_first, MAX = lp_ref_last, VALUE = lp_ref_start, EVENT_PRO = 'CRISPEX_SLIDER_LP_REF', $
					/DRAG, SENSITIVE = (refslid_sens AND ABS(eqnlps-1)))

	spatial_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Spatial', /COLUMN)
	cursor_frame		= WIDGET_BASE(spatial_tab, /FRAME, /COLUMN)
	x_slid			= WIDGET_SLIDER(cursor_frame, TITLE = 'X position of the cursor [pixel]', MIN = x_first, MAX = x_last, VALUE = x_start, EVENT_PRO = 'CRISPEX_SLIDER_X', /DRAG)
	y_slid			= WIDGET_SLIDER(cursor_frame, TITLE = 'Y position of the cursor [pixel]', MIN = y_first, MAX = y_last, VALUE = y_start, EVENT_PRO = 'CRISPEX_SLIDER_Y', /DRAG)

	lock_field		= WIDGET_BASE(control_panel, /FRAME, /ROW, /EXCLUSIVE)
	lockbut			= WIDGET_BUTTON(lock_field, VALUE = 'Lock to position', EVENT_PRO = 'CRISPEX_CURSOR_LOCK', TOOLTIP = 'Lock cursor to current position')
	WIDGET_CONTROL, lockbut, SET_BUTTON = 0
	unlockbut		= WIDGET_BUTTON(lock_field, VALUE = 'Unlock from position', TOOLTIP = 'Unlock cursor from current position')
	WIDGET_CONTROL, unlockbut, SET_BUTTON = 1

	zoom_sup_frame		= WIDGET_BASE(cursor_frame, /FRAME, /COLUMN)
	zoom_frame		= WIDGET_BASE(zoom_sup_frame, /ROW)
	zoom_label		= WIDGET_LABEL(zoom_frame, VALUE = 'Zoom:', /ALIGN_LEFT)
	zoom_but_field		= WIDGET_BASE(zoom_frame, /ROW, /EXCLUSIVE)
	zoom_one		= WIDGET_BUTTON(zoom_but_field, VALUE = '1x', EVENT_PRO = 'CRISPEX_ZOOMFAC_ONE',/NO_RELEASE)
	zoom_two		= WIDGET_BUTTON(zoom_but_field, VALUE = '2x', EVENT_PRO = 'CRISPEX_ZOOMFAC_TWO',/NO_RELEASE)
	zoom_three		= WIDGET_BUTTON(zoom_but_field, VALUE = '3x', EVENT_PRO = 'CRISPEX_ZOOMFAC_THREE',/NO_RELEASE)
	zoom_four		= WIDGET_BUTTON(zoom_but_field, VALUE = '4x', EVENT_PRO = 'CRISPEX_ZOOMFAC_FOUR',/NO_RELEASE)
	zoom_six		= WIDGET_BUTTON(zoom_but_field, VALUE = '6x', EVENT_PRO = 'CRISPEX_ZOOMFAC_SIX',/NO_RELEASE)
	zoom_eight		= WIDGET_BUTTON(zoom_but_field, VALUE = '8x', EVENT_PRO = 'CRISPEX_ZOOMFAC_EIGHT',/NO_RELEASE)
	WIDGET_CONTROL, zoom_one, SET_BUTTON = 1
	xpos_slider		= WIDGET_SLIDER(zoom_sup_frame, TITLE = 'X position of image (lower left corner) [pixel]', MIN = x_first, MAX = x_last, VALUE = x_first, EVENT_PRO = 'CRISPEX_SLIDER_XPOS',$
					/DRAG, SENSITIVE = 0)
	ypos_slider		= WIDGET_SLIDER(zoom_sup_frame, TITLE = 'Y position of image (lower left corner) [pixel]', MIN = y_first, MAX = y_last, VALUE = y_first, EVENT_PRO = 'CRISPEX_SLIDER_YPOS',$
					/DRAG, SENSITIVE = 0)

	stokes_tab		= WIDGET_BASE(tab_tlb, TITLE='Stokes', /COLUMN)
	pol_frame		= WIDGET_BASE(stokes_tab, /FRAME, /COLUMN)
	pol_main_label		= WIDGET_LABEL(pol_frame, VALUE = 'Stokes component:                                   ', /ALIGN_LEFT)
	pol_xy			= WIDGET_BASE(pol_frame, /ROW)
	pol_xy_label		= WIDGET_LABEL(pol_xy, VALUE = 'Main image:',/ALIGN_LEFT)
	pol_xy_buts		= WIDGET_BASE(pol_xy, /ROW, /EXCLUSIVE)
	pol_xy_i_but		= WIDGET_BUTTON(pol_xy_buts, VALUE = 'I', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_I', SENSITIVE = stokes_i_enabled, /NO_RELEASE)
	WIDGET_CONTROL, pol_xy_i_but, /SET_BUTTON
	pol_xy_q_but		= WIDGET_BUTTON(pol_xy_buts, VALUE = 'Q', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_Q', SENSITIVE = stokes_q_enabled, /NO_RELEASE)
	pol_xy_u_but		= WIDGET_BUTTON(pol_xy_buts, VALUE = 'U', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_U', SENSITIVE = stokes_u_enabled, /NO_RELEASE)
	pol_xy_v_but		= WIDGET_BUTTON(pol_xy_buts, VALUE = 'V', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_XY_V', SENSITIVE = stokes_v_enabled, /NO_RELEASE)
	pol_sp			= WIDGET_BASE(pol_frame, /ROW)
	pol_sp_label		= WIDGET_LABEL(pol_sp, VALUE = 'Detailed spectra:',/ALIGN_LEFT)
	pol_sp_buts		= WIDGET_BASE(pol_sp, /ROW, /NONEXCLUSIVE)
	spconstraint		= (nlp GT 1)
	pol_sp_i_but		= WIDGET_BUTTON(pol_sp_buts, VALUE = 'I', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_I', SENSITIVE = (spconstraint AND stokes_i_enabled))
	pol_sp_q_but		= WIDGET_BUTTON(pol_sp_buts, VALUE = 'Q', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_Q', SENSITIVE = (spconstraint AND stokes_q_enabled))
	pol_sp_u_but		= WIDGET_BUTTON(pol_sp_buts, VALUE = 'U', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_U', SENSITIVE = (spconstraint AND stokes_u_enabled))
	pol_sp_v_but		= WIDGET_BUTTON(pol_sp_buts, VALUE = 'V', EVENT_PRO = 'CRISPEX_DISPLAYS_STOKES_SELECT_SP_V', SENSITIVE = (spconstraint AND stokes_v_enabled))
	IF stokesfile THEN WIDGET_CONTROL, pol_sp_i_but, SET_BUTTON=(spconstraint AND stokes_i_enabled) ELSE WIDGET_CONTROL, pol_sp_i_but, SET_BUTTON=spconstraint
	WIDGET_CONTROL, pol_sp_q_but, SET_BUTTON=(spconstraint AND stokes_q_enabled)
	WIDGET_CONTROL, pol_sp_u_but, SET_BUTTON=(spconstraint AND stokes_u_enabled)
	WIDGET_CONTROL, pol_sp_v_but, SET_BUTTON=(spconstraint AND stokes_v_enabled)

	display_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Displays', /COLUMN)
	detspect_frame		= WIDGET_BASE(display_tab, /FRAME, /COLUMN)
	detspect_label_imref	= WIDGET_BASE(detspect_frame, /ROW)
	detspect_label		= WIDGET_LABEL(detspect_label_imref, VALUE = lswintitle[heightset]+':',/ALIGN_LEFT, /DYNAMIC_RESIZE)
	detspect_imref		= WIDGET_BASE(detspect_label_imref, /ROW, /EXCLUSIVE)
	detspect_im_but		= WIDGET_BUTTON(detspect_imref, VALUE = 'Main', EVENT_PRO = 'CRISPEX_DISPLAYS_DETSPECT_IM_SELECT', /NO_RELEASE, SENSITIVE = (nlp GT 1), $
					TOOLTIP = 'Main '+STRLOWCASE(lswintitle[heightset])+' display options')
	WIDGET_CONTROL, detspect_im_but, SET_BUTTON = (nlp GT 1)
	detspect_ref_but	= WIDGET_BUTTON(detspect_imref, VALUE = 'Reference', EVENT_PRO = 'CRISPEX_DISPLAYS_DETSPECT_REF_SELECT', /NO_RELEASE, SENSITIVE = showrefls, $
					TOOLTIP = 'Reference '+STRLOWCASE(lswintitle[refheightset])+' display options')
	detspect_buts		= WIDGET_BASE(detspect_frame, /ROW, /NONEXCLUSIVE)
	ls_toggle_but		= WIDGET_BUTTON(detspect_buts, VALUE = 'Display '+STRLOWCASE(lswintitle[heightset]), EVENT_PRO = 'CRISPEX_DISPLAYS_IMREF_LS_TOGGLE', /DYNAMIC_RESIZE)
	subtract_but		= WIDGET_BUTTON(detspect_buts, VALUE = 'Subtract average', EVENT_PRO = 'CRISPEX_DISPRANGE_LS_SUBTRACT', TOOLTIP = 'Subtract detailed spectrum from average spectrum')
	detspect_range		= WIDGET_BASE(detspect_frame, /ROW)
	lower_y_label		= WIDGET_LABEL(detspect_range, VALUE = 'Lower y-value:', /ALIGN_LEFT)
	lower_y_text		= WIDGET_TEXT(detspect_range, VALUE = STRTRIM(ls_low_y_init,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LS_LOW')
	upper_y_label		= WIDGET_LABEL(detspect_range, VALUE = 'Upper y-value:', /ALIGN_LEFT)
	upper_y_text		= WIDGET_TEXT(detspect_range, VALUE = STRTRIM(ls_upp_y_init,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LS_UPP')
	scale_detspect_buts	= WIDGET_BASE(detspect_frame, /ROW, /NONEXCLUSIVE)
	scale_detspect_but	= WIDGET_BUTTON(scale_detspect_buts, VALUE = 'Scale '+STRLOWCASE(lswintitle[heightset])+' to maximum of average', EVENT_PRO = 'CRISPEX_DISPRANGE_LS_SCALE_SELECT', $
					SENSITIVE = detspect_scale_enable, /DYNAMIC_RESIZE)
	WIDGET_CONTROL, scale_detspect_but, SET_BUTTON = detspect_scale
	all_other_disp		= WIDGET_BASE(display_tab, /FRAME, /COLUMN)
	other_label		= WIDGET_LABEL(all_other_disp, VALUE = 'Other displays:', /ALIGN_LEFT)
	other_disp		= WIDGET_BASE(all_other_disp, /ROW)
	slices_label		= WIDGET_LABEL(other_disp, VALUE = 'Slices:')
	other_buts		= WIDGET_BASE(other_disp, /ROW, /NONEXCLUSIVE)
	sp_toggle_but		= WIDGET_BUTTON(other_buts, VALUE = lp_h_capital[heightset]+'-t', EVENT_PRO = 'CRISPEX_DISPLAYS_SP_TOGGLE', TOOLTIP = 'Toggle display temporal '+STRLOWCASE(but_tooltip[heightset]))
	phis_toggle_but		= WIDGET_BUTTON(other_buts, VALUE = lp_h_capital[heightset]+'-Phi', EVENT_PRO = 'CRISPEX_DISPLAYS_PHIS_TOGGLE', TOOLTIP = 'Toggle display '+STRLOWCASE(but_tooltip[heightset])+' along a slit')
	refsp_toggle_but	= WIDGET_BUTTON(other_buts, VALUE = 'Reference '+lp_h_capital[refheightset]+'-t', EVENT_PRO = 'CRISPEX_DISPLAYS_REFSP_TOGGLE', $
					TOOLTIP = 'Toggle display reference temporal '+STRLOWCASE(but_tooltip[refheightset]))
	int_toggle_but		= WIDGET_BUTTON(other_buts, VALUE = 'I-t', EVENT_PRO = 'CRISPEX_DISPLAYS_INT_TOGGLE', SENSITIVE=(nt GT 1), TOOLTIP = 'Toggle display intensity versus time plot')
	images_disp		= WIDGET_BASE(all_other_disp, /ROW)
	images_label		= WIDGET_LABEL(images_disp, VALUE = 'Images:')
	images_buts		= WIDGET_BASE(images_disp, /ROW, /NONEXCLUSIVE)
	reference_but		= WIDGET_BUTTON(images_buts, VALUE = 'Reference image', EVENT_PRO = 'CRISPEX_DISPLAYS_REF_TOGGLE', SENSITIVE = showref, TOOLTIP = 'Toggle display reference image')
	WIDGET_CONTROL, reference_but, SET_BUTTON = showref
	doppler_but		= WIDGET_BUTTON(images_buts, VALUE = 'Doppler image', EVENT_PRO = 'CRISPEX_DISPLAYS_DOPPLER_TOGGLE', SENSITIVE = (nlp GT 1), TOOLTIP = 'Toggle display Doppler image')
	param_but_field		= WIDGET_BASE(all_other_disp, /ROW, /NONEXCLUSIVE)
	param_but		= WIDGET_BUTTON(param_but_field, VALUE = 'Parameters overview', EVENT_PRO = 'CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE', TOOLTIP = 'Toggle display parameters overview window')
	WIDGET_CONTROL, param_but, /SET_BUTTON

	scaling_tab 		= WIDGET_BASE(tab_tlb, TITLE = 'Scaling', /COLUMN)
	xy_scaling		= WIDGET_BASE(scaling_tab, /FRAME, /COLUMN)
	xy_ref_scaling		= WIDGET_BASE(xy_scaling, /ROW, /EXCLUSIVE)
	xy_scaling_but		= WIDGET_BUTTON(xy_ref_scaling, VALUE = 'Main image', EVENT_PRO = 'CRISPEX_SCALING_XY_SELECT', /NO_RELEASE)
	WIDGET_CONTROL, xy_scaling_but, /SET_BUTTON
	ref_scaling_but		 = WIDGET_BUTTON(xy_ref_scaling, VALUE = 'Reference image', EVENT_PRO = 'CRISPEX_SCALING_REF_SELECT', /NO_RELEASE, SENSITIVE = showref)
	dop_scaling_but		 = WIDGET_BUTTON(xy_ref_scaling, VALUE = 'Doppler image', EVENT_PRO = 'CRISPEX_SCALING_DOP_SELECT', /NO_RELEASE, SENSITIVE = 0)
	xy_scale_opts		= WIDGET_BASE(xy_scaling, /COLUMN)
	xy_scale_buts		= WIDGET_BASE(xy_scale_opts, /COLUMN, /EXCLUSIVE)
	auto_const		= WIDGET_BUTTON(xy_scale_buts, VALUE = 'Automatically, based on first slice', EVENT_PRO = 'CRISPEX_SCALING_AUTO_CONST', /NO_RELEASE, $
					TOOLTIP = 'Scale each image automatically using extrema from full FOV of first time step')
	WIDGET_CONTROL, auto_const, SET_BUTTON = 1
	auto_sing		= WIDGET_BUTTON(xy_scale_buts, VALUE = 'Automatically, per slice', EVENT_PRO = 'CRISPEX_SCALING_AUTO_SING', /NO_RELEASE, $
					TOOLTIP = 'Scale each image automatically using extrema from full FOV of same image')
	man_first		= WIDGET_BUTTON(xy_scale_buts, VALUE = 'Manually, based on first slice', EVENT_PRO = 'CRISPEX_SCALING_MAN_FIRST', /NO_RELEASE, $
					TOOLTIP = 'Scale each image manually using extrema from full FOV of first time step')
	man_curr		= WIDGET_BUTTON(xy_scale_buts, VALUE = 'Manually, based on current slice', EVENT_PRO = 'CRISPEX_SCALING_MAN_CURR', /NO_RELEASE, $
					TOOLTIP = 'Scale each image manully using extrema from FOV of current time step')
	abs_rel_buts 		= WIDGET_BASE(xy_scale_opts, /ROW, /EXCLUSIVE)
	abs_scale_but		= WIDGET_BUTTON(abs_rel_buts, VALUE = 'Absolute (0-255)', SENSITIVE = 0)
	rel_scale_but		= WIDGET_BUTTON(abs_rel_buts, VALUE = 'Relative (0-100%)', EVENT_PRO = 'CRISPEX_SCALING_REL', SENSITIVE = 0)
	WIDGET_CONTROL, abs_scale_but, /SET_BUTTON
	min_scale_slid		= WIDGET_SLIDER(xy_scaling, TITLE = 'Image minimum', MIN = 0, MAX = 254, VALUE = 0, EVENT_PRO = 'CRISPEX_SCALING_MIN_SLIDER', SENSITIVE = 0, /DRAG)
	max_scale_slid		= WIDGET_SLIDER(xy_scaling, TITLE = 'Image maximum', MIN = 1, MAX = 255, VALUE = 255, EVENT_PRO = 'CRISPEX_SCALING_MAX_SLIDER', SENSITIVE = 0, /DRAG)
	
	slit_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Slits',/COLUMN)
	slit_frame		= WIDGET_BASE(slit_tab, /FRAME, /COLUMN)
	slit_label		= WIDGET_LABEL(slit_frame, VALUE = 'Slit controls:', /ALIGN_LEFT)
	phi_slid		= WIDGET_SLIDER(slit_frame, TITLE = 'Slit angle [degrees]', MIN = 0, MAX = 179, VALUE = angle, EVENT_PRO = 'CRISPEX_SLIDER_PHI', SENSITIVE = 0, /DRAG)
	nphi_slid		= WIDGET_SLIDER(slit_frame, TITLE = 'Slit length [pixel]', MIN = 2, MAX = nphi, VALUE = nphi, EVENT_PRO = 'CRISPEX_SLIDER_NPHI', SENSITIVE = 0, /DRAG)
	slit_move_field		= WIDGET_BASE(slit_frame,/ROW)
	bwd_move_slit		= WIDGET_BUTTON(slit_move_field, VALUE = '< Move slit backwards', EVENT_PRO = 'CRISPEX_PHISLIT_MOVE_BWD', SENSITIVE= 0, TOOLTIP = 'Move slit backward along slit direction')
	fwd_move_slit		= WIDGET_BUTTON(slit_move_field, VALUE = 'Move slit forwards >', EVENT_PRO = 'CRISPEX_PHISLIT_MOVE_FWD', SENSITIVE = 0, TOOLTIP = 'Move slit forward along slit direction')
	loop_frame		= WIDGET_BASE(slit_tab, /FRAME, /COLUMN)
	loop_label		= WIDGET_LABEL(loop_frame, VALUE = 'Time slice along a loop:', /ALIGN_LEFT)
	loop_but_frame		= WIDGET_BASE(loop_frame, /ROW, /NONEXCLUSIVE)
	loop_slit_but		= WIDGET_BUTTON(loop_but_frame, VALUE = 'Draw loop path  ', EVENT_PRO = 'CRISPEX_LOOP_DEFINE') 
	loop_feedb_but		= WIDGET_BUTTON(loop_but_frame, VALUE = 'Path feedback', EVENT_PRO = 'CRISPEX_LOOP_FEEDBACK')
	WIDGET_CONTROL, loop_feedb_but, /SET_BUTTON
	loop_buts_frame		= WIDGET_BASE(loop_frame, /ROW)
	rem_loop_pt_but		= WIDGET_BUTTON(loop_buts_frame, VALUE = 'Remove last loop point', EVENT_PRO = 'CRISPEX_LOOP_REMOVE_POINT', SENSITIVE = 0)
	loop_slice_but		= WIDGET_BUTTON(loop_buts_frame, VALUE = 'Time slice along loop', EVENT_PRO = 'CRISPEX_DISPLAYS_LOOPSLAB_GET', SENSITIVE = 0)

	masks_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Mask', /COLUMN)
	masks			= WIDGET_BASE(masks_tab, /FRAME, /COLUMN)
	masks_overlay		= WIDGET_BASE(masks, /ROW)
	masks_overlay_label	= WIDGET_LABEL(masks_overlay, VALUE = 'Overlay on:',/ALIGN_LEFT)
	masks_overlay_buts	= CW_BGROUP(masks_overlay, ['Main','Reference','Doppler'],BUTTON_UVALUE=INDGEN(3),IDS=mask_button_ids,/NONEXCLUSIVE, /ROW, EVENT_FUNC = 'CRISPEX_MASK_OVERLAY')
	masks_overlay_color	= WIDGET_BASE(masks,/COLUMN)
	LOADCT,GET_NAMES=ctnames,/SILENT
	masks_overlay_ct_cbox	= WIDGET_COMBOBOX(masks_overlay_color, VALUE = STRTRIM(INDGEN(N_ELEMENTS(ctnames)),2)+REPLICATE(': ',N_ELEMENTS(ctnames))+ctnames, EVENT_PRO = 'CRISPEX_MASK_OVERLAY_SELECT_COLOR_TABLE', $
					SENSITIVE = maskfile)
	maskct = 13
	WIDGET_CONTROL, masks_overlay_ct_cbox, SET_COMBOBOX_SELECT = maskct
	masks_overlay_col_slid	= WIDGET_SLIDER(masks_overlay_color, MIN = 0, MAX = 255, VALUE = 255, TITLE = 'Color index', EVENT_PRO = 'CRISPEX_MASK_OVERLAY_COLOR_SLIDER', /DRAG, SENSITIVE = maskfile)

	overlays_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Miscellaneous', /COLUMN)
	overlays		= WIDGET_BASE(overlays_tab, /FRAME, /COLUMN)
	overlay_label		= WIDGET_LABEL(overlays, VALUE = 'Overlays:', /ALIGN_LEFT)
	overlay_buts		= WIDGET_BASE(overlays, /ROW)
	overlay_onebut		= WIDGET_BASE(overlay_buts, /NONEXCLUSIVE)
	overlay_but 		= WIDGET_BUTTON(overlay_onebut, VALUE = 'Saved loops:', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_MAIN')
	overlay_actbuts	 = WIDGET_BASE(overlay_buts, /ROW, /EXCLUSIVE)
	loop_overlay_al		= WIDGET_BUTTON(overlay_actbuts, VALUE = 'Always', EVENT_PRO = 'CRISPEX_RESTORE_LOOPS_ALWAYS', SENSITIVE = 0)
	WIDGET_CONTROL, loop_overlay_al, SET_BUTTON = 1
	loop_overlay_sav	= WIDGET_BUTTON(overlay_actbuts, VALUE = 'At saved '+STRLOWCASE(wav_h[heightset]), SENSITIVE = 0)
	linestyle_base		= WIDGET_BASE(overlays, /ROW)
	linestyle_label		= WIDGET_LABEL(linestyle_base, VALUE = 'Loops linestyle:', /ALIGN_LEFT)
	linestyle_buts		= WIDGET_BASE(linestyle_base, /ROW, /EXCLUSIVE)
	linestyle_0		= WIDGET_BUTTON(linestyle_buts, VALUE = 'solid', EVENT_PRO = 'CRISPEX_DRAW_LOOP_LINESTYLE_0')
	WIDGET_CONTROL, linestyle_0, SET_BUTTON = 1
	linestyle_1		= WIDGET_BUTTON(linestyle_buts, VALUE = 'dotted', EVENT_PRO = 'CRISPEX_DRAW_LOOP_LINESTYLE_1')
	linestyle_2		= WIDGET_BUTTON(linestyle_buts, VALUE = 'dashed', EVENT_PRO = 'CRISPEX_DRAW_LOOP_LINESTYLE_2')

	measuretool		= WIDGET_BASE(overlays_tab, /FRAME, /COLUMN)
	measure_label		= WIDGET_LABEL(measuretool, VALUE = 'Spatial measurement tool:', /ALIGN_LEFT)
	measure_buts		= WIDGET_BASE(measuretool, /ROW, /NONEXCLUSIVE)
	measure_but		= WIDGET_BUTTON(measure_buts, VALUE = 'Start measurement', EVENT_PRO = 'CRISPEX_MEASURE_ENABLE')
	apix_base		= WIDGET_BASE(measuretool, /ROW)
	apix_label		= WIDGET_LABEL(apix_base, VALUE = 'Pixel size:', /ALIGN_LEFT, SENSITIVE = 0)
  IF dx_fixed THEN $
  	apix_text		= WIDGET_LABEL(apix_base, VALUE = STRTRIM(dx,2), /ALIGN_LEFT, SENSITIVE = 0) $
  ELSE $
  	apix_text		= WIDGET_TEXT(apix_base, VALUE = STRTRIM(dx,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_MEASURE_ARCSEC', SENSITIVE = 0)
	apix_unit		= WIDGET_LABEL(apix_base, VALUE = '['+xunit+']', /ALIGN_LEFT, SENSITIVE = 0)
	measure_asec		= WIDGET_BASE(measuretool, /ROW)
	measure_asec_lab	= WIDGET_LABEL(measure_asec, VALUE = 'Distance [arcsec]:', /ALIGN_LEFT, SENSITIVE = 0)
	measure_asec_text	= WIDGET_LABEL(measure_asec, VALUE = '0.00', /DYNAMIC_RESIZE, SENSITIVE = 0)
	measure_km		= WIDGET_BASE(measuretool, /ROW)
	measure_km_lab		= WIDGET_LABEL(measure_km, VALUE = 'Distance [km]:', /ALIGN_LEFT, SENSITIVE = 0)
	measure_km_text		= WIDGET_LABEL(measure_km, VALUE = '0.00', /DYNAMIC_RESIZE, SENSITIVE = 0)

	bg = WIDGET_BASE(control_panel, EVENT_PRO = 'CRISPEX_PB_BG')
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (initializing control panel)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initializing control panel... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- SETTING UP DATA POINTERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (setting up data pointers)...",a5)','     ') 
	feedback_text = [feedback_text,'> Setting up data pointers... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	index = INTARR(nx,ny,2)
	xarr = INDGEN(nx)
	yarr = INDGEN(ny)
	FOR i=0,nx-1 DO BEGIN
		index[i,*,1] = yarr
	ENDFOR
	FOR j=0,ny-1 DO BEGIN
		index[*,j,0] = xarr
	ENDFOR
	indexmap= PTR_NEW(index, /NO_COPY)
	indices = PTR_NEW(INTARR(nx,ny,2))

	imdata	= PTR_NEW(imagefile, /NO_COPY)
	xyslice	= PTR_NEW(BYTARR(nx,ny))
	dopslice= PTR_NEW(BYTARR(nx,ny))
	emptydopslice= PTR_NEW(BYTARR(nx,ny))
	maskslice= PTR_NEW(BYTARR(nx,ny))
	maskdata = PTR_NEW(mask, /NO_COPY)

	IF (spfile EQ 1) THEN spdata = PTR_NEW(spectra, /NO_COPY) ELSE spdata = 0
	IF (refspfile EQ 1) THEN refspdata = PTR_NEW(referencespectra, /NO_COPY) ELSE refspdata = 0
	scan	= PTR_NEW(scanfile, /NO_COPY)
	IF (imtype EQ 1) THEN phiscan = PTR_NEW(BYTARR(nx,ny,nlp)) ELSE IF (imtype EQ 2) THEN phiscan = PTR_NEW(INTARR(nx,ny,nlp)) ELSE IF (imtype EQ 4) THEN phiscan = PTR_NEW(FLTARR(nx,ny,nlp))
	IF (imtype EQ 1) THEN sspscan = PTR_NEW(BYTARR(nx,ny,nlp*ns)) ELSE IF (imtype EQ 2) THEN sspscan = PTR_NEW(INTARR(nx,ny,nlp*ns)) ELSE IF (imtype EQ 4) THEN sspscan = PTR_NEW(FLTARR(nx,ny,nlp*ns))
	phislice= PTR_NEW(BYTARR(nlp,nphi))
	IF ((spfile EQ 1) OR (SINGLE_CUBE GE 1)) THEN BEGIN
		loopslab= PTR_NEW(FLTARR(nlp,nt,nphi))
		loopslice = PTR_NEW(BYTARR(nphi,nt))
		refloopslab= PTR_NEW(FLTARR(nlp,nt,nphi))
		refloopslice = PTR_NEW(BYTARR(nphi,nt))
		crossloc= PTR_NEW(INTARR(nphi))
		exact_loopslab= PTR_NEW(FLTARR(nlp,nt,nphi))
		exact_loopslice = PTR_NEW(BYTARR(nphi,nt))
		exact_crossloc= PTR_NEW(INTARR(nphi))
		rest_loopslab = PTRARR(nphi,/ALLOCATE_HEAP)
		rest_loopslice = PTRARR(nphi,/ALLOCATE_HEAP)
		rest_crossloc = PTRARR(nphi,/ALLOCATE_HEAP)
		FOR i=0,nphi-1 DO BEGIN
			*rest_loopslab[i]= PTR_NEW(0)
			*rest_loopslice[i] = PTR_NEW(0)
			*rest_crossloc[i]= PTR_NEW(0)
		ENDFOR
		det_loopslab= PTR_NEW(FLTARR(nlp,nt,nphi))
		det_loopslice = PTR_NEW(BYTARR(nphi,nt))
		det_crossloc= PTR_NEW(INTARR(nphi))
	ENDIF ELSE BEGIN
		loopslab = 0		&	loopslice = 0		&	crossloc = 0
		rest_loopslab = 0	&	rest_loopslice = 0	& 	rest_crossloc = 0
		det_loopslab = 0	&	det_loopslice = 0	& 	det_crossloc = 0
		exact_loopslice = 0	& 	exact_loopslab = 0	&	exact_crossloc = 0
		refloopslab = 0		&	refloopslice = 0
		WIDGET_CONTROL, loop_slit_but, SENSITIVE = 0
		WIDGET_CONTROL, loop_feedb_but, SENSITIVE = 0
	ENDELSE
	xp = PTR_NEW(FLTARR(1))
	yp = PTR_NEW(FLTARR(1))
	sxp = PTR_NEW(FLTARR(1))
	syp = PTR_NEW(FLTARR(1))
	
	refscan = 0
	refsspscan = 0
	IF (showref EQ 1) THEN BEGIN
		refdata	= PTR_NEW(referencefile, /NO_COPY)
		refslice= PTR_NEW(BYTARR(nx,ny))
		IF ((refnlp GT 1) AND (refspfile EQ 0)) THEN BEGIN
			refscan = PTR_NEW(refscanfile, /NO_COPY)
			IF (refimtype EQ 1) THEN refsspscan = PTR_NEW(BYTARR(nx,ny,refnlp*refns)) ELSE IF (refimtype EQ 2) THEN refsspscan = PTR_NEW(INTARR(nx,ny,refnlp*refns)) ELSE $
				IF (refimtype EQ 4) THEN refsspscan = PTR_NEW(FLTARR(nx,ny,refnlp*refns))
		ENDIF
	ENDIF ELSE BEGIN
		refdata = 0	&	refslice = 0	&	refcube = ''	
	ENDELSE

	imwintitle = 'CRISPEX'+instance_label+': Main image'
	CRISPEX_WINDOW, imwinx, imwiny, control_panel, imwintitle, imwin, xywid, DRAWID = xydrawid, $
		DRAWBASE = drawbase, 0, 0;, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_XYREF_RESIZE'
	WIDGET_CONTROL, xydrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,$
		/DRAW_BUTTON_EVENTS
	
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (setting up data pointers)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	IF startupwin THEN BEGIN
		WSET, startupwid
		WSHOW, startupwid
		feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Setting up data pointers... done!']
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	ENDIF

;--------------------------------------------------------------------------------- DEFINE INFO POINTER
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (defining info pointer)...",a5)','     ') 
	feedback_text = [feedback_text,'> Defining info pointer... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- CONTROL PANEL BUTTON/TEXT REFERENCES
	ctrlscp = { $
		save_loop_pts:save_loop_pts, timeslicemenu:timeslicemenu, clear_current_estimate:clear_current_estimate, $			
		sel_saved_loop:sel_saved_loop, all_saved_loop:all_saved_loop, det_file_but:det_file_loop, $
		fbwd_button:fbwd_button, bwd_button:backward_button, pause_button:pause_button, $
		fwd_button: forward_button, ffwd_button:ffwd_button, $
		loop_button:loop_button, blink_button:blink_button, cycle_button:cycle_button, $
		t_slider:t_slid, lower_t_text:lower_t_text, upper_t_text:upper_t_text, $		
		reset_trange_but:reset_trange_but, slice_button:slice_update_but, $			
		t_speed_slider:t_speed_slid, t_step_slider:t_step_slid, imref_blink_but:imref_blink_but, $					
		lower_lp_text:lower_lp_text, upper_lp_text:upper_lp_text, $	
		reset_lprange_but:reset_lprange_but, lp_slider:lp_slid, $				
		lp_speed_slider:lp_speed_slid, lp_blink_slider:lp_blink_slid, $	
		lp_blink_but:lp_blink_but, lp_ref_but:lp_ref_but, lp_ref_slider:lp_ref_slid, $
		x_slider:x_slid, y_slider:y_slid, lock_button:lockbut, unlock_button:unlockbut, $
		zoom_one:zoom_one, zoom_two:zoom_two, zoom_three:zoom_three, zoom_four:zoom_four,$ 		
		zoom_six:zoom_six, zoom_eight:zoom_eight, xpos_slider:xpos_slider, ypos_slider:ypos_slider, $			
		pol_xy_i_but:pol_xy_i_but, pol_xy_q_but:pol_xy_q_but, pol_xy_u_but:pol_xy_u_but,$		
		pol_sp_i_but:pol_sp_i_but, pol_sp_q_but:pol_sp_q_but, pol_sp_u_but:pol_sp_u_but,$		
		pol_sp_v_but:pol_sp_v_but, pol_xy_v_but:pol_xy_v_but, $						
		detspect_label:detspect_label, scale_detspect_but:scale_detspect_but, $
		detspect_im_but:detspect_im_but, detspect_ref_but:detspect_ref_but, $
		ls_toggle_but:ls_toggle_but, subtract_but:subtract_but, $		
		lower_y_text:lower_y_text, upper_y_text:upper_y_text, $	
		sp_toggle_but:sp_toggle_but, refsp_toggle_but:refsp_toggle_but, int_toggle_but:int_toggle_but, $
		phis_toggle_but:phis_toggle_but, param_but:param_but, reference_but:reference_but, doppler_but:doppler_but, $						
		xy_scaling_but:xy_scaling_but, ref_scaling_but:ref_scaling_but, dop_scaling_but:dop_scaling_but, $				
		auto_const:auto_const, auto_sing:auto_sing, man_first:man_first, man_curr:man_curr, $		
		abs_scale_but:abs_scale_but, rel_scale_but:rel_scale_but, $
		max_scale_slider:max_scale_slid, min_scale_slider:min_scale_slid, $				
		phi_slider:phi_slid, nphi_slider:nphi_slid, loop_feedb_but:loop_feedb_but, $	
		bwd_move_slit:bwd_move_slit, fwd_move_slit:fwd_move_slit, $
		loop_slit_but:loop_slit_but, rem_loop_pt_but:rem_loop_pt_but, loop_slice_but:loop_slice_but, $	
		overlay_but:overlay_but, loop_overlay_all:loop_overlay_al, loop_overlay_sav:loop_overlay_sav, $		
		linestyle_0:linestyle_0, linestyle_1:linestyle_1, linestyle_2:linestyle_2, $			
		measure_but:measure_but, apix_label:apix_label, apix_unit:apix_unit, apix_text:apix_text, $
		measure_asec_lab:measure_asec_lab, measure_asec_text:measure_asec_text, $
		measure_km_lab:measure_km_lab, measure_km_text:measure_km_text, $					
		mask_button_ids:mask_button_ids, masks_overlay_ct_cbox:masks_overlay_ct_cbox, $
		masks_overlay_col_slid:masks_overlay_col_slid, dispwid:dispwid, clear_current_inst:clear_current_inst, $
		verbose_set:PTR_NEW([sh_verb_0,sh_verb_4,sh_verb_8,sh_verb_16]) $
	}
;--------------------------------------------------------------------------------- RETRIEVE LOOPS CONTROLS 
	ctrlsdet = { $
		sel_all:0, sel_none:0, disp_list:0, overlay_all:0, overlay_sel:0, all_pos:0, saved_pos:0, $
		sel_range_pos:0, save_imonly:0, save_refonly:0, save_imref:0, $
		dtmin_text:0, dtmax_text:0, dlpmin_text:0, dlpmax_text:0, refdlpmin_text:0, refdlpmax_text:0, $	
		width_slider:0, get_dets:0 $
	}
;--------------------------------------------------------------------------------- FEEDBACK CONTROLS
	ctrlsfeedb = { $
		estimate_label:0, feedback_text:0, close_button:0 $
	}
;--------------------------------------------------------------------------------- INT CONTROLS 
	ctrlsint = { $
		int_sel_all:0, int_sel_none:0, int_sel_save:0, lower_y_int_text:0, upper_y_int_text:0, $
		lower_t_int_text:0, upper_t_int_text:0, reset_trange_but:0 $
	}
;--------------------------------------------------------------------------------- RETRIEVE DETECTIONS CONTROLS 
	ctrlsloop = { $
		get_loops:0, sel_all:0, sel_none:0, all_pos:0, saved_pos:0, del_files:0, keep_files:0, $
		save_imonly:0, save_refonly:0, save_imref:0 $
	}
;--------------------------------------------------------------------------------- BMP BUTTON IMAGES
	ctrlspbbut = { $
		fbwd_pressed:bmpbut_fbwd_pressed, fbwd_idle:bmpbut_fbwd_idle, $
		bwd_pressed:bmpbut_bwd_pressed, bwd_idle:bmpbut_bwd_idle, $
		pause_pressed:bmpbut_pause_pressed, pause_idle:bmpbut_pause_idle, $
		fwd_pressed:bmpbut_fwd_pressed, fwd_idle:bmpbut_fwd_idle, $
		ffwd_pressed:bmpbut_ffwd_pressed, ffwd_idle:bmpbut_ffwd_idle, $
		loop_pressed:bmpbut_loop_pressed, loop_idle:bmpbut_loop_idle, $
		cycle_pressed:bmpbut_cycle_pressed, cycle_idle:bmpbut_cycle_idle, $
		blink_pressed:bmpbut_blink_pressed, blink_idle:bmpbut_blink_idle $
	}
;--------------------------------------------------------------------------------- PARAMETER WINDOW REFERENCES
	ctrlsparam = { $
		x_coord_val:0, y_coord_val:0, t_coord_val:0, lp_coord_val:0, $					
		zoom_val:0, act_t_val:0, act_lp_val:0, v_dop_val:0, $						
		lp_ref_coord_val:0, v_dop_ref_val:0, act_lp_ref_val:0, $
		img_val:0, imgsc_val:0, ref_val:0, refsc_val:0, stokes_val:0 $			
	}
;--------------------------------------------------------------------------------- PREFERENCE BUTTONS
	ctrlspref = { $
		startup_autopl:0, startup_win:0, displays_bgcols:0, displays_plcols:0, displays_interp:0, $
		displays_phislice:0, displays_slices:0, $								
		paths_i_def_but:0, paths_i_sav_but:0, paths_ipath_text:0, $
		paths_o_def_but:0, paths_o_sav_but:0, paths_opath_text:0, $
		paths_iopath:0, save_defsaveid:0, save_defsaveid_sample:0, $
		set_defaults:0 $
	}
;--------------------------------------------------------------------------------- RESTORE LOOPS CONTROLS 
	ctrlsrestore = { $
		disp_list:0, sel_all:0, sel_none:0, open_tanat:0 $	
	}
;--------------------------------------------------------------------------------- SAVING CONTROLS 
	ctrlssav = { $
		path_textlab:0, savopt_path_textlab:0, overlays_num_but:0, overlays_curs_but:0, overlays_thick_slider:0, $
		overlays_pts_but:0, overlays_symsize_slider:0, overlays_asecbar_but:0, overlays_asecbar_slider:0 $
	}
;--------------------------------------------------------------------------------- CONTROL SWITCHES
	ctrlsswitch = { $
		imrefdetspect:0, lp_ref_lock:eqnlps, bwd_insensitive:1, fwd_insensitive:1 $
	}
;--------------------------------------------------------------------------------- CURSOR 
	curs = { $
		sx:sx_start, sy:sy_start, sxlock:sx_start, sylock:sx_start, $		
		xlock:x_start, ylock:y_start, lockset:0 $				
	}
;--------------------------------------------------------------------------------- DATA 
	data = { $
		imagedata:imdata, xyslice:xyslice, refdata:refdata, refslice:refslice, maskdata:maskdata, maskslice:maskslice, $
		dopslice:dopslice, spdata:spdata, sspscan:sspscan, refspdata:refspdata, refscan:refscan, refsspscan:refsspscan, $						
		emptydopslice:emptydopslice, scan:scan, phiscan:phiscan, phislice:phislice, $				
		indexmap:indexmap, indices:indices, $	
		lur:lur, lun:lun, luf:luf, lufs:lufs, lum:lum $
	}
;--------------------------------------------------------------------------------- DATA PARAMETERS
	dataparams = { $
		imfilename:imcube, spfilename:spcube, reffilename:refcube, refspfilename:refspcube, maskfilename:maskcube, $	
		x:x_start, y:y_start, d_nx:nx, d_ny:ny, nx:nx, ny:ny, $								
		lc:lc, lp:lp_start, lp_ref:lp_ref_start, lp_dop:lp_start, nlp:nlp, refnlp:refnlp, ns:ns, s:0, $					
		lps:lps, ms:ms, spec:mainspec, $
		reflps:reflps, refms:refms, refspec:refspec, $
		t:t_start, nt:nt, refnt:refnt,masknt:masknt, $		
    dx:dx, dy:dy, $
    bunit:[bunit,refbunit], lpunit:[lpunit,reflpunit], xunit:xunit, yunit:yunit $
	}
;--------------------------------------------------------------------------------- DATA SWITCH
	dataswitch = { $
		onecube:onecube, reffile:showref, refspfile:refspfile, spfile:spfile, maskfile:maskfile $							
	}
;--------------------------------------------------------------------------------- DET PARAMS
	detparams = { $
		nr_dets:0, sel_dets:PTR_NEW(INTARR(nphi*3)), t_restored:PTR_NEW(INTARR(nphi)), $
		detfilename:'', idx:0, nr_sel_loops:0, sel_loops:PTR_NEW(INTARR(nphi)), $
		xlr_restored:PTR_NEW(BYTARR(nphi,nphi)), ylr_restored:PTR_NEW(BYTARR(nphi,nphi)), $	
		xlp_restored:PTR_NEW(LONARR(nphi,nphi)), ylp_restored:PTR_NEW(LONARR(nphi,nphi)), $	
		lpsizes:PTR_NEW(LONARR(nphi)), lrsizes:PTR_NEW(LONARR(nphi)), $				
		width:0, mid:0, delta_t_dn:10, delta_t_up:10, lp_dn:lp_first, lp_up:lp_last, $					
		lp_ref_dn:lp_ref_first, lp_ref_up:lp_ref_last $
	}
;--------------------------------------------------------------------------------- DATA DISPLAY PARAMETERS
	dispparams = { $
		t_first:t_first, t_last:t_last, t_range:nt, t_low:t_first, t_upp:t_last, $					
		x_first:x_first, x_last:x_last, y_first:y_first, y_last:y_last, $				
		lp_first:lp_first, lp_last:lp_last, lp_range:nlp, lp_low:lp_first, lp_upp:lp_last, $
		lp_ref_first:lp_ref_first, lp_ref_last:lp_ref_last, lp_ref_low:lp_ref_first, lp_ref_upp:(refnlp-1), lp_ref_range:refnlp, $
		nlpreb:nlpreb, ntreb:ntreb, refnlpreb:nlpreb, refntreb:refntreb, refloopnlxreb:nlpreb, refloopntreb:ntreb, $
		loopnlxreb:nlpreb, loopntreb:loopntreb, restloopnlxreb:nlpreb, restloopntreb:loopntreb, $
		retrdetnlxreb:nlpreb, retrdetntreb:loopntreb, phisnlpreb:nlpreb, nphireb:nphireb, $					
		xi:xi, yi:yi, xo:xo, yo:yo, xi_ref:xi_ref, yi_ref:yi_ref, xo_ref:xo_ref, yo_ref:yo_ref, $
		interpspslice:interpspslice, phislice_update:phislice_update, slices_imscale:slices_imscale $				
	}
;--------------------------------------------------------------------------------- DATA DISPLAY SWITCHES
	dispswitch = { $
		restricted_t_range:PTR_NEW(0), restricted_lp_range:PTR_NEW(0), restricted_lp_ref_range:PTR_NEW(0), $
		exts:exts_set, refexts:refexts_set, warpspslice:warpspslice, warprefspslice:warprefspslice, $
		detspect_scale:detspect_scale, ref_detspect_scale:ref_detspect_scale, drawdop:0 $	
	}
;--------------------------------------------------------------------------------- FEEDBACK PARAMS
	feedbparams = { $
		estimate_lx:estimate_lx, estimate_time:estimate_time, estimate_run:estimate_run, $		
		startup_im:startup_im, xout:xout, yout:yout, verbosity:verbosity, last_routine:'', last_routine_count:0, $
		pbstats:DOUBLE(SYSTIME(/SECONDS)), sum_pbstats:PTR_NEW(DBLARR(10)), av_pbstats:0D, count_pbstats:0 $	
	}
;--------------------------------------------------------------------------------- INT PARAMS
	intparams = { $
		sel_diagnostics:PTR_NEW(sel_diagnostics), lines_diagnostics:PTR_NEW(lines_diagnostics), $
		linlab_diagnostics:linlab_diagnostics, colors_diagnostics:colors_diagnostics, $
		collab_diagnostics:collab_diagnostics, selcol_diagnostics:PTR_NEW(selcol_diagnostics), $
		diagnostics:diagnostics, lock_t:1 $ 
	}
;--------------------------------------------------------------------------------- LOOP PARAMS
	loopparams = { $
		xp:xp, yp:yp, xr:PTR_NEW(FLTARR(nphi)), yr:PTR_NEW(FLTARR(nphi)), np:0, $			
		w_lpts:PTR_NEW(BYTARR(nphi)), nw_lpts:0 $						
	}
;--------------------------------------------------------------------------------- LOOPS DATA
	loopsdata = { $
		loopsize:0, rest_loopsize:PTR_NEW(INTARR(10)), exact_loopsize:0, det_loopsize:0, $				
		crossloc:crossloc, rest_crossloc:rest_crossloc, exact_crossloc:exact_crossloc, $		
		det_crossloc:det_crossloc, $
		loopslab:loopslab, loopslice:loopslice, $			
		refloopslab:refloopslab, refloopslice:refloopslice, $			
		rest_loopslab:rest_loopslab, rest_loopslice:rest_loopslice, $
		exact_loopslab:exact_loopslab, exact_loopslice:exact_loopslice, $	
		det_loopslab:det_loopslab, det_loopslice:det_loopslice $	
	}
;--------------------------------------------------------------------------------- LOOPS SWITCHES
	loopswitch = { $
		retrieve_loops:0, restore_loops:0, was_restore_loops:0, retrieve_detfile:0 $		
	}
;--------------------------------------------------------------------------------- MEASUREMENT
	meas = { $
		arcsecpix:dx, spatial_measurement:0, np:0, $					
		xp:PTR_NEW(FLTARR(1)), yp:PTR_NEW(FLTARR(1)), $					
		sxp:PTR_NEW(FLTARR(1)), syp:PTR_NEW(FLTARR(1)) $					
	}
;--------------------------------------------------------------------------------- OVERLAY PARAMS
	overlayparams = { $
		sxp:sxp, syp:syp, sxr:PTR_NEW(FLTARR(nphi)), syr:PTR_NEW(FLTARR(nphi)), $
		loop_linestyle:0, maskcolor:255, maskct:maskct $				
	}
;--------------------------------------------------------------------------------- OVERLAY SWITCHES
	overlayswitch = { $
		det_overlay_all:1, loopslit:0, overlalways:1, looppath_feedback:1, mask:maskfile, maskim:[maskfile,showref,0] $		
	}
;--------------------------------------------------------------------------------- PARAMETER WINDOW CONTROLS 
	paramparams = { $
		wav_h:wav_h, scale_cubes:scale_cubes_vals, img_get_val:0., $
    imgsc_get_val:0., ref_get_val:0., refsc_get_val:0. $ 
	}
;--------------------------------------------------------------------------------- PARAM SWITCHES
	paramswitch = { $
		img_get:KEYWORD_SET(vals_img), ref_get:KEYWORD_SET(vals_ref), dt_set:dt_set $ 
	}
;--------------------------------------------------------------------------------- PATHS AND DIRECTORIES
	paths = { $
		ipath:ipath, opath:opath, opath_write:opath_write, hostname:hostname, $
		dir_cpft:dir_cpft, dir_cpft_write:dir_cpft_write, dir_settings:dir_settings, $
		dir_aux:dir_aux, dir_tanat:dir_aux, tanat_repointed:0, dir_inst_write:dir_inst_write, $
		dir_inst:dir_inst $
	}
;--------------------------------------------------------------------------------- PLAYBACK
	pbparams = { $
		t_step:t_step, t_speed:t_speed, direction:direction, bg:bg, mode:'PAUSE', lmode:'LOOP', $				
		imrefmode:0, lp_step:1, spmode:0, spdirection:1 $						
	}
;--------------------------------------------------------------------------------- PHI SLIT VARIABLES 
	phiparams = { $
		d_nphi_set:nphi, nphi:nphi, angle:angle, nphi_set:nphi-1, sphi:0, $			
		nw_cur:nphi, curindex:0, maxindex:0 $				
	}
;--------------------------------------------------------------------------------- PLOTAXES
	plotaxes = { $
		ticklen:ticklen, lsxticklen:lsxticklen, lsyticklen:lsyticklen, $
		reflsxticklen:reflsxticklen, reflsyticklen:reflsyticklen, $
		spxticklen:spxticklen, spyticklen:spyticklen, $
		refspxticklen:spxticklen, refspyticklen:spyticklen, $
		phisxticklen:phisxticklen, phisyticklen:phisyticklen, $
		loopxticklen:spxticklen, loopyticklen:spyticklen, $
		refloopxticklen:spxticklen, refloopyticklen:spyticklen, $
		restloopxticklen:spxticklen, restloopyticklen:spyticklen, $
		retrdetxticklen:spxticklen, retrdetyticklen:spyticklen, $
		intxticklen:intxticklen, intyticklen:intyticklen, $
		ls_low_y:ls_low_y, ls_upp_y:ls_upp_y, ls_yrange:ls_yrange, $		
		ls_low_y_ref:ls_low_y_ref, ls_upp_y_ref:ls_upp_y_ref, ls_yrange_ref:ls_yrange_ref, $
		int_low_y:int_low_y, int_upp_y:int_upp_y, int_low_t:t_first, int_upp_t:t_last, $
		dt:dt, v_dop:v_dop, v_dop_ref:v_dop_ref $		
	}
;--------------------------------------------------------------------------------- PLOT PARAMETERS
	plotparams = { $
		bgplotcol:bgplotcol, plotcol:plotcol $
	}
;--------------------------------------------------------------------------------- PLOT POSITION
	plotpos = { $
		lsx0:lsx0, lsy0:lsy0, lsx1:lsx1, lsy1:lsy1, $
		reflsx0:reflsx0, reflsy0:reflsy0, reflsx1:reflsx1, reflsy1:reflsy1, $
		spx0:spx0, spy0:spy0, spx1:spx1, spy1:spy1, $							
		refspx0:spx0, refspy0:refspy0, refspx1:spx1, refspy1:refspy1, $
		phisx0:phisx0, phisx1:phisx1, phisy0:phisy0, phisy1:phisy1, $
		loopx0:spx0, loopx1:spx1, loopy0:spy0, loopy1:loopy1, $
		refloopx0:spx0, refloopx1:spx1, refloopy0:spy0, refloopy1:loopy1, $
		restloopx0:spx0, restloopx1:spx1, restloopy0:spy0, restloopy1:loopy1, $
		retrdetx0:spx0, retrdetx1:spx1, retrdety0:spy0, retrdety1:loopy1, $
		intx0:intx0, intx1:intx1, inty0:inty0, inty1:inty1, $
		xplspw:xplspw, yplspw:yplspw, refxplspw:xplspw, refyplspw:refyplspw, $		
		phisxplspw:phisxplspw, phisyplspw:phisyplspw, loopxplspw:xplspw, loopyplspw:loopyplspw, $
		refloopxplspw:xplspw, refloopyplspw:loopyplspw, restloopxplspw:xplspw, restloopyplspw:loopyplspw, $
		retrdetxplspw:xplspw, retrdetyplspw:loopyplspw, $
		lsxmargin_init:lsxmargin_init, lsxwall_init:lsxwall_init, $
		reflsxmargin_init:lsxmargin_init, reflsxwall_init:lsxwall_init, $
		spxmargin_init:spxmargin_init, spxwall_init:spxwall_init, $
		refspxmargin_init:spxmargin_init, refspxwall_init:spxwall_init, $
		phisxmargin_init:spxmargin_init, phisxwall_init:spxwall_init, $
		loopxmargin_init:spxmargin_init, loopxwall_init:spxwall_init, $
		refloopxmargin_init:spxmargin_init, refloopxwall_init:spxwall_init, $
		restloopxmargin_init:spxmargin_init, restloopxwall_init:spxwall_init, $
		retrdetxmargin_init:spxmargin_init, retrdetxwall_init:spxwall_init, $
		intxmargin_init:intxmargin_init, intxwall_init:intxwall_init $
	}
;--------------------------------------------------------------------------------- PLOT SWITCHES
	plotswitch = { $
		heightset:heightset, refheightset:refheightset, stokesfile:stokesfile, scalestokes:scalestokes, $
		v_dop_set:v_dop_set, v_dop_set_ref:v_dop_set_ref, subtract:0, ref_subtract:0 $						
	}
;--------------------------------------------------------------------------------- PLOT TITLES
	plottitles = { $
		spxtitle:spxtitle, spytitle:spytitle, spwintitle:spwintitle, $
		refspxtitle:refspxtitle, refspwintitle:refspwintitle, $
		lsytitle:lsytitle, reflsytitle:reflsytitle, $
		lswintitle:lswintitle, reflswintitle:reflswintitle, phiswintitle:phiswintitle $
	}
;--------------------------------------------------------------------------------- PREFERENCE SETTINGS
	prefs = { $
		autoplay:autoplay, startupwin:startupwin, defsaveid:defsaveid,  $
		defipath:defipath, prefipath:prefipath, defopath:defopath, prefopath:prefopath, $
		bgplotcol_old:bgplotcol, plotcol_old:plotcol, interpspslice_old:interpspslice, $
		slices_imscale_old:slices_imscale, $
		tmp_autoplay:autoplay, tmp_startupwin:startupwin, tmp_defsaveid:defsaveid, $
		tmp_bgplotcol:bgplotcol, tmp_plotcol:plotcol, tmp_defipath:defipath, tmp_prefipath:prefipath, $
		tmp_defopath:defopath, tmp_prefopath:prefopath, tmp_interpspslice:interpspslice, $
		tmp_phislice_update:phislice_update, tmp_slices_imscale:slices_imscale, $								
		default_autoplay:default_autoplay, default_startupwin:default_startupwin, $
		default_bgplotcol:default_bgplotcol, default_plotcol:default_plotcol, $
		default_defipath:default_defipath, default_prefipath:default_prefipath, $
		default_defopath:default_defopath, default_prefopath:default_prefopath, $
		default_defsaveid:default_defsaveid, default_interpspslice:default_interpspslice, $
		default_phislice_update:default_phislice_update, default_slices_imscale:default_slices_imscale, $
		preview:0 $					
	}
;--------------------------------------------------------------------------------- RESTORED LOOPS PARAMS
	restoreparams = { $
		cfilecount:0, cfiles:PTR_NEW(STRARR(nphi)), sel_loops:PTR_NEW(INTARR(nphi)), $		
		xr_restored:PTR_NEW(BYTARR(nphi,nphi)), yr_restored:PTR_NEW(BYTARR(nphi,nphi)), $		
		xp_restored:PTR_NEW(LONARR(nphi,nphi)), yp_restored:PTR_NEW(LONARR(nphi,nphi)), $		
		psizes:PTR_NEW(LONARR(nphi)), rsizes:PTR_NEW(LONARR(nphi)), lp_restored:PTR_NEW(INTARR(nphi)), $
		disp_loopfile:'0', disp_loopnr:PTR_NEW(-1), disp_imref:PTR_NEW(-1), disp_slices:PTR_NEW(0), $
		disp_ref_slices:PTR_NEW(0)$
	}
;--------------------------------------------------------------------------------- RETRIEVE LOOP PARAMS
	retrparams = { $
		clfilecount:0, clfiles:PTR_NEW(STRARR(nphi)), sel_loops:PTR_NEW(INTARR(nphi)), $	
		retrieve_files:PTR_NEW(STRARR(nphi)), retrieve_filecount:0, $			
		lpsizes:PTR_NEW(LONARR(nphi)), lrsizes:PTR_NEW(LONARR(nphi)), $					
		xlr_restored:PTR_NEW(BYTARR(nphi,nphi)), ylr_restored:PTR_NEW(BYTARR(nphi,nphi)), $		
		xlp_restored:PTR_NEW(LONARR(nphi,nphi)), ylp_restored:PTR_NEW(LONARR(nphi,nphi)) $		
	}
;--------------------------------------------------------------------------------- SAVING PARAMS
	savparams = { $
		savpro:'', snapshot:0, lp_orig:lp_start, lp_ref_orig:lp_ref_start, filename_text:0, overlays_incl:0, overlays_thick:1, $
		overlays_curs:1, overlays_num:1, overlays_pts:1, linescan_ls:0, overlays_symsize:1, $
		overlays_asecbar:0, overlays_asecbar_length:1, overlays_asecbar_pix:1. $
	}
;--------------------------------------------------------------------------------- SAVING SWITCHES
	savswitch = { $
		cont:0, delete_clsav:0, imref_only:1, det_imref_only:1, $				
		all_pos_loops:1, pos_dets:1 $							
	}
;--------------------------------------------------------------------------------- SESSION PARAMS
	sesparams = { $
		rest_sessions:0, sessions:PTR_NEW(INTARR(nphi)), csesfiles:PTR_NEW(STRARR(nphi)), $
		curr_instance_id:LONG(set_instance_id), instance_label:instance_label $		
	}
;--------------------------------------------------------------------------------- SCALING
	scaling = { $
		imagescale:imagescale, imagemin:immin, imagemax:immax, scale_range:PTR_NEW([0.,0.,0.]), scale_minimum:PTR_NEW([0.,0.,0.]), $	
		scale_max_val:PTR_NEW([255.,255.,255.]), scale_min_val:PTR_NEW([0.,0.,0.]), $				
		rel_scale_max_val:PTR_NEW([100.,100.,100.]), rel_scale_min_val:PTR_NEW([0.,0.,0.]), $			
		refmin:refmin, refmax:refmax, imrefscaling:0, relative:relative_scaling, $	
		scalestokes_max:scalestokes_max, dopplermin:dopplermin, dopplermax:dopplermax $
	}
;--------------------------------------------------------------------------------- STOKES PARAMS
	stokesparams = { $
		labels:stokes_labels, select_sp:stokes_select_sp, $	
		prev_select_sp:stokes_select_sp $
	}
;--------------------------------------------------------------------------------- VERSION INFO
	versioninfo = { $
		version_number:version_number, revision_number:revision_number $				
	}
;--------------------------------------------------------------------------------- WINDOW IDs
	winids = { $
		root:control_panel, imtlb:imwin, imwid:xywid, $		
		reftlb:0, refwid:0, refdrawid:0, refdrawbase:0, $				
		doptlb:0, dopwid:0, dopdrawid:0, dopdrawbase:0, $
		imrefdisp:0, imreftlb:0, imrefwid:0, imrefdrawid:0, imrefdrawbase:0, $					
		lstlb:0, lswid:0, lsdrawid:0, reflstlb:0, reflswid:0, reflsdrawid:0, $
		sptlb:0, spwid:0, spdrawid:0, refsptlb:0, refspwid:0, refspdrawid:0, $
		phistlb:0, phiswid:0, phisdrawid:0, $				
		looptlb:0, loopwid:0, loopdrawid:0, $						
		reflooptlb:0, refloopwid:0, refloopdrawid:0, $						
		restlooptlb:PTR_NEW(0), restloopwid:PTR_NEW(0), restloopdrawid:PTR_NEW(0), $					
		retrdettlb:0, retrdetwid:0, retrdetdrawid:0, $					
		intwid:0, inttlb:0, intdrawid:0, intmenutlb:0, $
		savetlb:0, detsavetlb:0, restoretlb:0, preftlb:0, $							
		estimatetlb:0, savewintlb:0, saveoptwintlb:0, restsestlb:0, paramtlb:0, $					
		feedbacktlb:0, abouttlb:0, errtlb:0, warntlb:0, restsesfeedbtlb:0, $
		imwintitle:imwintitle, spwintitle:'',lswintitle:'',refwintitle:'',refspwintitle:'',reflswintitle:'', $
		imrefwintitle:'',dopwintitle:'',phiswintitle:'',restloopwintitle:PTR_NEW(''),retrdetwintitle:'',$
		loopwintitle:'',refloopwintitle:'',intwintitle:'' $
	}
;--------------------------------------------------------------------------------- WINDOW (RE)SIZES 
	winsizes = { $
		aboutwinx:startup_nx, aboutwiny:startup_ny, xywinx:imwinx, xywiny:imwiny, $					
		lswinx:lswinx, lswiny:lswiny, lsxres:lswinx, lsyres:lswiny, $
		spwinx:spwinx, spwiny:spwiny, spxres:spwinx, spyres:spwiny, $
		reflsxres:reflswinx, reflsyres:reflswiny, refspxres:spwinx, refspyres:refspwiny, $
		phisxres:phiswinx, phisyres:phiswiny, refloopxres:spwinx, refloopyres:spwiny, $				
		loopxres:spwinx, loopyres:spwiny, restloopxres:spwinx, restloopyres:spwiny, $			
		retrdetxres:spwinx, retrdetyres:spwiny, intxres:lswinx, intyres:intwiny, $							
		xdelta:xdelta, ydelta:ydelta, spxoffset:spxoffset, lsxoffset:lsxoffset, $
		aboutxoffset:startup_xpos, aboutyoffset:startup_ypos $
		}
;--------------------------------------------------------------------------------- WINDOW SWITCHES
	winswitch = { $
		showls:0, showsp:0, showrefsp:refspfile, showrefls:showrefls, showimref:0, $
		estimate_win:0, showphis:0, showref:showref, showparam:0, showint:0, $
		showloop:0, showrefloop:0, showrestloop:0, showretrdet:0, showdop:0, dispwids:0 $
	}
;--------------------------------------------------------------------------------- ZOOMING 
	zooming = { $
		factor:1, xpos:0., ypos:0. $								
	}
;--------------------------------------------------------------------------------- DEFINE INFO POINTER
	info = { $
		ctrlscp:PTR_NEW(ctrlscp, /NO_COPY), $
		ctrlsdet:PTR_NEW(ctrlsdet, /NO_COPY), $
		ctrlsfeedb:PTR_NEW(ctrlsfeedb, /NO_COPY), $
		ctrlsint:PTR_NEW(ctrlsint, /NO_COPY), $
		ctrlsloop:PTR_NEW(ctrlsloop, /NO_COPY), $
		ctrlspbbut:PTR_NEW(ctrlspbbut, /NO_COPY), $
		ctrlsparam:PTR_NEW(ctrlsparam, /NO_COPY), $
		ctrlspref:PTR_NEW(ctrlspref, /NO_COPY), $
		ctrlsrestore:PTR_NEW(ctrlsrestore, /NO_COPY), $
		ctrlssav:PTR_NEW(ctrlssav, /NO_COPY), $
		ctrlsswitch:PTR_NEW(ctrlsswitch, /NO_COPY), $
		curs:PTR_NEW(curs, /NO_COPY), $
		data:PTR_NEW(data, /NO_COPY), $
		dataparams:PTR_NEW(dataparams, /NO_COPY), $
		dataswitch:PTR_NEW(dataswitch, /NO_COPY), $
		detparams:PTR_NEW(detparams, /NO_COPY), $
		dispparams:PTR_NEW(dispparams, /NO_COPY), $
		dispswitch:PTR_NEW(dispswitch, /NO_COPY), $
		feedbparams:PTR_NEW(feedbparams, /NO_COPY), $
		intparams:PTR_NEW(intparams, /NO_COPY), $
		loopsdata:PTR_NEW(loopsdata, /NO_COPY), $
		loopparams:PTR_NEW(loopparams, /NO_COPY), $
		loopswitch:PTR_NEW(loopswitch, /NO_COPY), $
		meas:PTR_NEW(meas, /NO_COPY), $
		overlayparams:PTR_NEW(overlayparams, /NO_COPY), $
		overlayswitch:PTR_NEW(overlayswitch, /NO_COPY), $
		paramparams:PTR_NEW(paramparams, /NO_COPY), $
		paramswitch:PTR_NEW(paramswitch, /NO_COPY), $
		paths:PTR_NEW(paths, /NO_COPY), $
		pbparams:PTR_NEW(pbparams, /NO_COPY), $
		phiparams:PTR_NEW(phiparams, /NO_COPY), $
		plotaxes:PTR_NEW(plotaxes, /NO_COPY), $
		plotparams:PTR_NEW(plotparams, /NO_COPY), $
		plotpos:PTR_NEW(plotpos, /NO_COPY), $
		plotswitch:PTR_NEW(plotswitch, /NO_COPY), $
		plottitles:PTR_NEW(plottitles, /NO_COPY), $
		prefs:PTR_NEW(prefs, /NO_COPY), $
		restoreparams:PTR_NEW(restoreparams, /NO_COPY), $
		retrparams:PTR_NEW(retrparams, /NO_COPY), $
		savparams:PTR_NEW(savparams, /NO_COPY), $
		savswitch:PTR_NEW(savswitch, /NO_COPY), $
		sesparams:PTR_NEW(sesparams, /NO_COPY), $
		scaling:PTR_NEW(scaling, /NO_COPY), $
		stokesparams:PTR_NEW(stokesparams, /NO_COPY), $
		versioninfo:PTR_NEW(versioninfo, /NO_COPY), $
		winids:PTR_NEW(winids, /NO_COPY), $
		winsizes:PTR_NEW(winsizes, /NO_COPY), $
		winswitch:PTR_NEW(winswitch, /NO_COPY), $
		zooming:PTR_NEW(zooming, /NO_COPY) $
	}
	info = PTR_NEW(info, /NO_COPY)
	WIDGET_CONTROL, control_panel, SET_UVALUE = info

	WIDGET_CONTROL, imwin, SET_UVALUE = info

	pseudoevent = { WIDGET_BUTTON, id:control_panel, top:control_panel, handler:0L, select:1 }
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (defining info pointer)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Defining info pointer... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- DETERMINE DISPLAY OF PLOTS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (determining display of plots)...",a5)','     ')
	feedback_text = [feedback_text,'> Determining display of plots... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	IF ((*(*info).dataswitch).spfile EQ 1) THEN BEGIN
		WIDGET_CONTROL, sp_toggle_but, SET_BUTTON = 1
		(*(*info).winswitch).showsp = 1
	ENDIF ELSE BEGIN
		IF (SINGLE_CUBE EQ 1) THEN BEGIN
			(*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_low
			(*(*info).winswitch).showphis = 0
			(*(*info).winswitch).showls = 0
			WIDGET_CONTROL, approxmenu, SENSITIVE = 0
			WIDGET_CONTROL, save_ex_slab_but, SENSITIVE = 0
			WIDGET_CONTROL, ls_toggle_but, SENSITIVE = 0
			WIDGET_CONTROL, subtract_but, SENSITIVE = 0
			WIDGET_CONTROl, lower_y_text, SENSITIVE = 0
			WIDGET_CONTROL, upper_y_text, SENSITIVE = 0
			WIDGET_CONTROL, phis_toggle_but, SENSITIVE = 0
		ENDIF ELSE BEGIN
			IF (SINGLE_CUBE GE 1) THEN WIDGET_CONTROL, approxmenu, SENSITIVE = 0
			WIDGET_CONTROL, phis_toggle_but, SET_BUTTON = 1
			WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 0
			(*(*info).winswitch).showphis = 1
			WIDGET_CONTROL, (*(*info).ctrlscp).phi_slider, SENSITIVE = 1
			WIDGET_CONTROL, (*(*info).ctrlscp).nphi_slider, SENSITIVE = 1
			(*(*info).ctrlsswitch).bwd_insensitive = 0
			(*(*info).ctrlsswitch).fwd_insensitive = 0
			WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = (*(*info).ctrlsswitch).bwd_insensitive+1
			WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = (*(*info).ctrlsswitch).fwd_insensitive+1
		ENDELSE
		WIDGET_CONTROL, sp_toggle_but, SENSITIVE=0
		IF ((onecube EQ 0) OR (nt EQ 1)) THEN BEGIN
			WIDGET_CONTROL, fbwd_button, SENSITIVE = 0, SET_VALUE = bmpbut_fbwd_idle
			WIDGET_CONTROL, backward_button, SENSITIVE = 0, SET_VALUE = bmpbut_bwd_idle
			WIDGET_CONTROL, pause_button, SENSITIVE = 0, SET_VALUE = bmpbut_pause_idle
			WIDGET_CONTROL, forward_button, SENSITIVE = 0, SET_VALUE = bmpbut_fwd_idle
			WIDGET_CONTROL, ffwd_button, SENSITIVE = 0, SET_VALUE = bmpbut_ffwd_idle
			WIDGET_CONTROL, loop_button, SENSITIVE = 0, SET_VALUE = bmpbut_loop_idle
			WIDGET_CONTROL, cycle_button, SENSITIVE = 0, SET_VALUE = bmpbut_cycle_idle
			WIDGET_CONTROL, blink_button, SENSITIVE = 0, SET_VALUE = bmpbut_blink_idle
			WIDGET_CONTROL, t_slid, SENSITIVE = 0
			WIDGET_CONTROL, t_step_slid, SENSITIVE = 0
			WIDGET_CONTROL, upper_t_text, SET_VALUE = '0', SENSITIVE = 0
			WIDGET_CONTROL, lower_t_text, SET_VALUE = '0', SENSITIVE = 0
			WIDGET_CONTROL, save_as_jpg_all, SENSITIVE = 0
			WIDGET_CONTROL, save_as_mpeg, SENSITIVE = 0
			WIDGET_CONTROL, auto_sing, SENSITIVE = 0
			WIDGET_CONTROL, man_curr, SENSITIVE = 0
			WIDGET_CONTROL, retrievemenu, SENSITIVE = 0
			*(*(*info).data).scan = (*(*(*info).data).scan)[0]
			*(*(*info).data).phiscan = *(*(*info).data).scan
		ENDIF 
		WIDGET_CONTROL, loop_slit_but, SENSITIVE = exts_set
		WIDGET_CONTROL, loop_feedb_but, SENSITIVE = exts_set
		WIDGET_CONTROL, savemenu, SENSITIVE = exts_set
		WIDGET_CONTROL, retrievemenu, SENSITIVE = exts_set
	ENDELSE
	IF (N_ELEMENTS(SINGLE_CUBE) EQ 1) THEN BEGIN
		WIDGET_CONTROL, ls_toggle_but, SET_BUTTON = (SINGLE_CUBE NE 1) 
		(*(*info).winswitch).showls = (SINGLE_CUBE NE 1) 
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, ls_toggle_but, SET_BUTTON = 1
		(*(*info).winswitch).showls = 1
	ENDELSE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (determining display of plots)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Determining display of plots... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
;--------------------------------------------------------------------------------- OPENING WINDOWS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (realizing widget)...",a5)','     ')
	feedback_text = [feedback_text,'> Realizing widget... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	CRISPEX_MASK_BUTTONS_SET, pseudoevent
	IF (*(*info).winswitch).showsp THEN BEGIN
		(*(*info).winswitch).showsp = 0
		CRISPEX_DISPLAYS_SP_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF
	IF (*(*info).winswitch).showls THEN BEGIN
		(*(*info).winswitch).showls = 0
		(*(*info).ctrlsswitch).imrefdetspect = 0
		CRISPEX_DISPLAYS_IMREF_LS_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF
	IF (*(*info).winswitch).showref THEN BEGIN
		(*(*info).winswitch).showref = 0
		CRISPEX_DISPLAYS_REF_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF
	IF (*(*info).winswitch).showrefsp THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).refsp_toggle_but, SET_BUTTON = 1
		(*(*info).winswitch).showrefsp = 0
		CRISPEX_DISPLAYS_REFSP_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF ELSE WIDGET_CONTROL, (*(*info).ctrlscp).refsp_toggle_but, SENSITIVE = 0
	IF (*(*info).winswitch).showrefls THEN BEGIN
		IF ((*(*info).dataswitch).refspfile EQ 0) THEN *(*(*info).data).refsspscan = (*(*(*info).data).refscan)[0]
		(*(*info).winswitch).showrefls = 0
		(*(*info).ctrlsswitch).imrefdetspect = 1
		CRISPEX_DISPLAYS_IMREF_LS_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF
	IF (*(*info).winswitch).showphis THEN BEGIN
		(*(*info).winswitch).showphis = 0
		CRISPEX_DISPLAYS_PHIS_TOGGLE, pseudoevent
		IF startupwin THEN WSHOW, startupwid
	ENDIF

	CRISPEX_FIND_CLSAV, pseudoevent
	IF ((*(*info).retrparams).clfilecount NE 0) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 1
		WIDGET_CONTROL, (*(*info).ctrlscp).all_saved_loop, SENSITIVE = 1
	ENDIF

	CRISPEX_UPDATE_T, pseudoevent
	CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE, pseudoevent
	IF ((*(*info).dataswitch).spfile EQ 1) OR (*(*info).dataswitch).onecube THEN CRISPEX_UPDATE_SLICES, pseudoevent
	IF showrefls THEN BEGIN
		IF (ref_detspect_scale EQ 0) THEN BEGIN
			CRISPEX_DISPRANGE_LS_SCALE_REF, pseudoevent
			CRISPEX_DISPRANGE_LS_RANGE, pseudoevent
		ENDIF
		CRISPEX_DISPLAYS_DETSPECT_IM_SELECT, pseudoevent
	ENDIF
	IF (detspect_scale EQ 0) THEN BEGIN
		CRISPEX_DISPRANGE_LS_SCALE_MAIN, pseudoevent
		CRISPEX_DISPRANGE_LS_RANGE, pseudoevent
	ENDIF
	CRISPEX_DRAW, pseudoevent
	
	IF ((*(*info).winswitch).showsp OR (*(*info).winswitch).showphis OR (*(*info).winswitch).showrefsp) THEN spwset = 1 ELSE spwset = 0
	IF (*(*info).winswitch).showls THEN lsoffset = (*(*info).winswitch).showls * ((*(*info).winsizes).lswiny + (*(*info).winsizes).ydelta) + (*(*info).dataswitch).refspfile * (*(*info).winsizes).ydelta ELSE $
		lsoffset = (*(*info).winswitch).showrefls * (reflswiny + (*(*info).winsizes).ydelta)
	WIDGET_CONTROL, control_panel, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).xywinx + (*(*info).winsizes).xdelta + spwset * ((*(*info).winsizes).xdelta + (*(*info).winsizes).spwinx), $
		TLB_SET_YOFFSET = lsoffset
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (realizing widget)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	IF startupwin THEN BEGIN
		WSET, startupwid
		feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Realizing widget... done!']
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	ENDIF

;--------------------------------------------------------------------------------- START MANAGING
	IF (TOTAL(verbosity[0:1]) GE 1) THEN WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (start managing)...",a5)','     ')
	feedback_text = [feedback_text,'> Start managing... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	XMANAGER, 'CRISPEX', control_panel, /NO_BLOCK
	XMANAGER, 'CRISPEX', imwin, /NO_BLOCK
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN 
		WRITEU,-1,STRING(FORMAT='(%"\rCRISPEX SETUP: Setting up widget (start managing)...",a5)','done!') 
		PRINT,'' 
	ENDIF
	feedback_text = [feedback_text[0]+'done!',feedback_text[1:N_ELEMENTS(feedback_text)-2],'> Start managing... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	WAIT,0.1
	IF (*(*info).prefs).autoplay THEN CRISPEX_PB_FORWARD, pseudoevent
	IF resave_preferences THEN CRISPEX_PREFERENCES_SAVE_SETTINGS, pseudoevent, /RESAVE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN PRINT,'CRISPEX SETUP: Set-up done!'
	IF startupwin THEN BEGIN
		WSET, startupwid
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Set-up done!'
		WAIT,0.25
		WIDGET_CONTROL, startuptlb, /DESTROY
	ENDIF
	CRISPEX_VERBOSE_SET_BUTTONS, pseudoevent
END
