;+
; NAME:
;   CRISPEX: CRIsp SPectral EXplorer
;
; PURPOSE:
;   This procedure allows browsing and analysis of multidimensional data cubes. The procedure was
;   initially written for analysis of temporal spectropolarimetric data obtained with the CRISP 
;   instrument at the Swedish Solar Telescope (SST) at La Palma (Spain), however any data formatted 
;   in a particular way (cf. the online reference pages) may be fed to CRISPEX for browsing and/on
;	  analysis purposes
;
; CATEGORY:
;   Data browsing and analysis
;
; CALLING SEQUENCE:
;   CRISPEX, Imcube
;
; INPUTS:
;	  Imcube        - 3D image datacube (dimensions [nx,ny,nt*nlp*ns]) or, if SPCUBE is not provided,
;                   a scan (dimensions [nx,ny,nlp]). If SINGLE_CUBE is specified a 3D datacube may 
;                   be provided even if SPCUBE is not. Allowed data types: BYTE, INTEGER, FLOAT.
;
; OPTIONAL INPUTS:
;	  Spcube        - 3D spectral datacube (dimensions [nlp,nt,nx*ny*ns]). Required input if 
;                   SINGLE_CUBE is not supplied with the number of linepositions and Imcube is a 
;                   time series. Allowed data types: BYTE, INTEGER, FLOAT.
;
; KEYWORDS:
;	  REFCUBE	      - Reference data of same spatial dimensions as main data. REFCUBE may be supplied
;                   with either:
;			                * 2D data array: dimensions [nx,ny].
;			                * 3D data array: dimensions [nx,ny,refnt], where refnt must be equal to nt.
;			                * Scalar string: pointing to a reference image cube of dimensions [nx,nx,
;                       refnt*refnlp] where refnt must be equal to 1 or to nt, but refnlp need not 
;                       be equal to nlp.
;		 	                * 2-element string array: first element pointing to a reference image cube 
;                       (dimensions [nx,ny,nt*refnlp]), the second element pointing to the 
;                       corresponding reference spectral cube (dimensions [refnlp,refnt,nx*ny]). 
;                       Again, refnt must be equal to 1 or nt, while any (positive integer) value of 
;                       refnlp is allowed. 
;                   Allowed data types: BYTE, INTEGER, FLOAT.
;	  MASKCUBE	    - Mask data of same spatial dimensions as main data. Must be supplied with a scalar
;                   string pointing to a mask image cube of dimensions [nx,nx,masknt] where masknt 
;                   must be equal to nt or 1. Allowed data types: BYTE, INTEGER, FLOAT.
;	  SPECTFILE	    - File containing the normalised spectrum as function of the linepositions or 
;                   wavelength. Contains the required variables:
;			                * NORM_SPECT (normalised spectrum)
;                     * NORM_FACTOR (normalisation factor used to produce NORM_SPECT) 
;                     * SPECT_POS (spectral wavelength positions in any desired unit). 
;                       Optional variables: XTITLE and YTITLE (x- and y-title label for plots, 
;                       respectively). If not set, the mean spectrum is determined from the scan(s) 
;                       determined by the MNSPEC keyword (i.e. from all x,y pixels at all 
;                       linepositions at t given by the setting of MNSPEC).
;                   SPECTFILE may be supplied with either:
;			                * Scalar string: spectral file corresponding to main data.
;			                * 2-element string array: first element is the spectral file corresponding to 
;                       the main data, second element corresponds to the reference data.
;			              Set first element to '' if you only want to set the reference spectral file, e.g. 
;                   SPECTFILE=['','reference.spectfile'].
;                   Setting of SPECTFILE will be ignored if FITS cubes are provided to Imcube, Spcube
;                   and/or REFCUBE.
;	  LINE_CENTER	  - Specifies line centre or wavelength information. LINE_CENTER may be supplied with:
;			                * Integer scalar: linecentre is set to position specified by LINE_CENTER.
;			                * 1D 2-element array (format: [WAVELENGTH, DELTA_LAMBDA]): linecentre is 
;                       determined from the data and set to WAVELENGTH. The distance in wavelength 
;                       between the linepositions is specified by DELTA_LAMBDA.
;			                * 1D 3-element array (format: [Integer scalar, WAVELENGTH, DELTA_LAMBDA]): 
;                       combination of the two above.
;			                * 2D 1-element integer array (format: [[Integer scalar], [Integer scalar]]): 
;                       first element sets the linecentre position for the main data, the second 
;                       element that of the reference data.
;			                * 2D 2-element array (format: [[MAIN_WAVELENGTH, MAIN_DELTA_LAMBDA],
;                       [REF_WAVELENGTH, REF_DELTA_LAMBDA]]): the elements from the first subarray 
;                       set the linecentre wavelength and wavelength spacing for the main data, 
;                       while those of the second subarray set those values for the reference data.
;			                * 2D 3-element array (format: [[Integer scalar, MAIN_WAVELENGTH, 
;                       MAIN_DELTA_LAMBDA],[Integer scalar, REF_WAVELENGTH, REF_DELTA_LAMBDA]]): 
;                       combination of the two above.
;                   If not set: linecentre is determined from the data or from SPECTFILE.
;                   Setting of LINE_CENTER will be ignored if FITS cubes are provided to Imcube,
;                   Spcube and/or REFCUBE.
;	  MNSPEC        - Determines the calculation of the average spectrum for in-program display. 
;                   MNSPEC may be supplied with:
;			                * Integer scalar: mean spectrum is determined from the t=MNSPEC scan.
;			                * 2-element integer array: mean spectrum is determined from the t=MNSPEC[0] 
;                       through t=MNSPEC[1] scans.
;                   If not set: mean spectrum is determined from the t=0 scan.
;	  DT            - Specifies the elapsed time in seconds per time step. Defaults to not defined, 
;                   showing frame number instead of time on the vertical axes.
;	  EXTS          -	If set, the time slices/slabs displayed in the program will be exact timeslices, 
;			              obtained through linear interpolation, rather than approximated timeslices, 
;                   obtained through nearest-neighbour interpolation (which is the default setting). 
;                   Note that setting this keyword may slow down the browsing of the spectral range 
;                   (i.e. movement of the spectral slider) considerably, because of the
;                   computationally more expensive linear interpolation.
;	  SINGLE_CUBE	  - Single integer value specifying the number of spectral positions of the 
;                   datacube. Only to be used when IMCUBE is provided with a 3D spectrotemporal 
;                   datacube and SPCUBE is not specified.
;	  SCALE_STOKES  - If set, the detailed and average spectra of Stokes Q, U and/or V will be scaled
;                   to the maximum of Stokes I (i.e. I/I, Q/I, U/I and/or V/I). If not set, each 
;                   Stokes component will be scaled to its respective maximum.
;	  SCALE_CUBES	  - Specifies the value that the data should be multiplied with. SCALE_CUBES may be
;                   supplied with:
;			                * Integer scalar: value used for the main data.
;			                * 2-element integer array: first element is used for the main data, the second
;                       for the reference data.
;                   Defaults to 1.
;	  VALS_IMG	    - If set, the value of the pixel of the image image under the cursor will be 
;                   returned in the parameters overview window.
;	  VALS_REF	    - If set, the value of the pixel of the reference image under the cursor will be 
;                   returned in the parameters overview window.
;	  NO_WARP		    - Prevents the warping of the temporal spectrum when the wavelength spacing is 
;                   non-equidistant. Applies equally to both sets of data if reference data is 
;                   supplied. Defaults to not set.
;	  XTITLE		    - Sets the x-title of the temporal spectrum and the detailed spectrum. XTITLE may
;                   be supplied with:
;			                * Scalar string: x-title labels for the main temporal spectrum and detailed 
;                       spectrum.
;			                * 2-element string array: first element sets the x-title labels for the main 
;                       temporal spectrum and detailed spectrum, while the second element sets the 
;                       labels for the reference data.
;			              Set first element to '' if you only want to set the reference x-title, e.g. 
;                   XTITLE=['','Height [km]'].
;   YTITLE		    - Sets the y-title of the detailed spectrum. YTITLE may be supplied with:
;			                * Scalar string: y-title label for the detailed spectrum.
;			                * 2-element string array: first element sets the y-title label for the main 
;                       detailed spectrum, while the second element sets the label for the 
;                       reference data.
;			              Set first element to '' if you only want to set the reference y-title, e.g. 
;                   YTITLE=['','Velocity [km/s]'].
;	  WINDOW_LARGE	- Override the "1:1 window scaling whenever possible" setting. Useful for data with
;                   small nx and/or ny, where 1:1 image display is possible, but would yield small 
;                   image display windows. Defaults to not set.
;	  VERBOSE		    - Verbosity setting for program setup and running. Mainly for maintenance purposes. 
;                   Verbosity levels can be set by supplying the keyword with the following values 
;                   (add bitwise):
;				              * 0:  No verbosity.
;				              * 1:  Basic setup verbosity.
;				              * 2:  Extended setup verbosity.
;				              * 4:  Basic runtime verbosity.
;				              * 8:  Extended runtime verbosity.
;				              * 16: Playback statistics verbosity.
;			              In practice the values 3 and 7 are not useful and 1 has become obsolete for all 
;                   practical purposes. All values larger than 26 are reduced to 26, all values 
;                   smaller than 0 are set to 0. Note that verbosity levels may also be set in-program.
;
; OUTPUTS:
;	  Window outputs:
;     Depending on the call. In all cases at least three windows will appear, one of which being the
;		  control panel, one the main image window and one the detailed spectrum window. If called with:
;			  - one datacube, then an additional window will open, showing the spectrum along a slit;
;			  - two datacubes, then an additional window will open, showing the temporal spectrum;
;			  - three datacubes, then two additional windows will open, one showing the temporal
;			    spectrum and one showing the reference image from the third cube;
;			  - four datacubes, then four additional windows will open, one showing the main temporal spectrum,
;			    one showing the reference temporal spectrum, one showing the reference detailed spectrum and
;			    one showing the reference image.
;		  Additional data display windows may be accessed through the tabs in the control panel in all 
;		  cases, although not all data display windows may be available, depending on the number of data 
;		  cubes with which the program is called.
;	  Saveable outputs:
;		  - loop path points (*.CLSAV file);
;		  - loop path or detection space-time diagram (*.CSAV file). Note that only that part of the
;       space-time diagram between the current lower and upper t-values will be saved, i.e., if you
;       wish to save it for the full range in time, be sure to reset the temporal boundaries.
;	  	- intensity versus time data for specific linepositions (*.CINT file);
;		  - (selected) timeseries as MPEG movie;
;		  - selected frame as JPEG snapshot;
;	  	- (selected) timeseries as JPEG files;
;		  - current session (*.CSES file).
;
; RESTRICTIONS:
;   Requires the following procedures and functions:
;     Functions: HISTO_OPT()                  [general]
;                FITSPOINTER(), READFITS()    [if reading FITS cubes]
;
; PROCEDURE:
;   In default setting, four windows are opened, one control panel and three subsidiary windows 
;   containing the main image, the temporal spectrum and the detailed spectrum, respectively. 
;   Additional data browsing and analysis options/windows may be obtained through control panel 
;   options.
;
;	  Some example calling sequences are discussed below. The most common calling sequence is:
;
;		  CRISPEX, 'main.imcube', 'main.spcube'
;
;	  Spectral information may be supplied either through the use of a spectral save file or the use 
;   of the LINE_CENTER keyword, e.g.:
;
;		  CRISPEX, 'main.imcube', 'main.spcube', SPECTFILE='main.spectfile'
;	  or
;		  CRISPEX, 'main.imcube', 'main.spcube', LINE_CENTER=[6562,0.1]
;
;	  Reference data may be viewed by supplying such data to the REFCUBE keyword, either in the 
;   simple reference mode:
;
;		  CRISPEX, 'main.imcube', 'main.spcube', REFCUBE='ref.imcube'
;
;	  or in the dual cube mode:
;
;		  CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube']
;
;	  In addition, when running in dual cube mode, one may provide arrays to certain keywords in 
;   order to set options for both the main and the reference data, e.g.:
;	
;		  CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			  SPECTFILE=['main.spectfile','ref.specftile']
;	  or
;		  CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			  LINE_CENTER=[[6562,0.1],[8542,0.1]]
;	  etc.
;
;	  One may also set some of the options for the reference cube only (while retaining the default 
;   options for the main data) by setting the element corresponding to the main data to an empty 
;   scalar string, e.g.:
;
;		  CRISPEX, 'main.imcube', 'main.spcube', REFCUBE=['ref.imcube', 'ref.spcube'], 
;			  YTITLE=['','Velocity [km/s]']
;
;	  Shortcut controls (requires focus on the control panel to work):
;		  - File menu:
;			  - 'About' (Ctrl+A);
;			  - 'Preferences (Ctrl+P);
;			  - 'Help' (Ctrl+H, redirecting to online reference pages);
;			  - 'Quit' (Ctrl+Q).
;		  - Control menu:
;			  - 'Playback options':
;				  - 'one frame backward' (Shift+B);
;				  - 'play backwards' (Shift+Backspace);
;				  - 'pause' (Shift+Spacebar);
;				  - 'play forwards' (Shift+Tab);
;				  - 'one frame forwards' (Shift+F).
;			  - 'Spectral options':
;				  - 'increase spectral position' (Shift+S); 
;				  - 'decrease spectral position' (Shift+A).
;			  - 'Zoom options':
;		  		- 'In' (Ctrl+Shift+I); 
;				  - 'Out' (Ctrl+Shift+O).
;			  - 'Runtime options:
;				  - 'Interrupt' (Ctrl+Shift+C);
;         - 'Verbosity':
;           - 'No verbosity' (Shift+0);
;           - 'Basic runtime verbosity' (Shift+4);
;           - 'Extended runtime verbosity' (Shift+8);
;           - 'Enable playback statistics' (Shift+P).
;
;	  Further information on the calling sequence and specific options can be found in the online 
;   reference pages, which can be reached through the in-program help function or by going directly 
;   to http://folk.uio.no/~gregal/crispex
;
; MODIFICATION HISTORY:
;	  2009 Feb 11 GV: start buildup of program, based on Øystein Langangen's August 2008 
;			              version of CRISP_SPECTRAL_EXPLORER.PRO. Further (structural) inspiration from 
;                   XIMOVIE.PRO by Oivind Wikstol and XSLICE.PRO by Alfred de Wijn
;	  2009 Feb 21 GV: incorporation of LP_HEADER.PRO functionality
;	  2009 Feb 24 GV: release of beta version (v0.9)
;	  2009 Mar 02 GV: implementation of x- and y-slice display (v0.9.1)
;	  2009 Mar 04 GV: implementation of display of spectral slice along a slit, including controls to
;			              set slit angle and length (v0.9.3)
;	  2009 Mar 09 GV: implementation of movement along the slit direction	(v0.9.4)
;	  2009 Mar 10 GV: implementation of option to study spectral scan and show a reference image or
;			              cube (v0.9.5)
;	  2009 Mar 15 GV: implementation of zoom option (v0.9.6)
;	  2009 Mar 23 GV: implementation of time slice along a segmented line and	option to save the
;			              resulting data(v0.9.7)
;	  2009 Apr 14 GV: corrected extraction and saving of time slice, removed x- and y-slice display
;			              due to redundancy from spectral slice along a slit (v0.9.8)
;	  2009 May 15 GV: implementation of extended save and retrieve options, extended image scaling	
;			              options, save and restore session options, adjustable time and spectral range, 
;			              different loop linestyles, selection menu for loop overlays, parameter overview 
;			              window and display of saved timeslices (v0.9.9)
;	  2009 Aug 24 GV: implementation of extra zoomfactors, save as MPEG and JPEG options, save from	
;			              detection file, spectral and temporal range choice in saving timeslabs, exact 
;			              timeslice display in-program, resizable display windows, read-in of full 
;			              reference cubes and an option to calculate mean spectrum over a range in 
;			              timesteps (v1.0)
;	  2009 Nov 05 GV: implementation of loop path feedback, shortcut controls through keyboard,	
;			              single full cube call and also fixed a number of bugs (v1.1)
;	  2010 Mar 10 GV: enabled visualisation of Stokes cube data, display of reference and image cube	
;			              data values, drawing of loop paths for 3D temporal image cube, retreival and
;			              saving timeslices from reference cube, extended scaling options for reference
;			              cube image, implemented spatial measurement tool and help function, moved user
;			              feedback to pop-up windows, disposed of obsolete keywords and fixed a number
;			              of bugs (v1.5)
;	  2011 Jul 11 GV: enabled dual cube mode, in-program Doppler images, extended (Stokes) spectral	
;			              options, extended plot options, blinking while playing, setting of preferences, 
;			              extraction of intensity-time plots, and made aesthetic improvements (bitmap
;			              play buttons, better cursor visibility and startup screen) (v1.6)
;	  2011 Aug 22 GV: extended save as options, implemented save as PNG, saving of color MPEG, an	
;			              option to open a restored loop in TANAT and fixed a number of bugs(v1.6.1)
;	  2012 Mar 23 GV: extended reference and Stokes data input options, implemented mask overlays,	
;			              enabled display of multiple space-time diagrams and fixed a number of bugs (v1.6.2)
;	  2012 Dec 04 GV: extended in-program analysis options through display of reference space-time 	
;			              diagram, extended image and space-time diagram scaling options and fixed a
;			              number of bugs (v1.6.3)
;
; ACKNOWLEDGEMENTS:
;	  This code would not be present in its current state and with the current functionalities without
;	  the relentless practical testing and the valuable input and ideas of Luc Rouppe van der Voort,
;	  Sven Wedemeyer-Böhm, Mats Carlsson, Patrick Antolin, Jorrit Leenaarts, Bart de Pontieu, 
;	  Eamon Scullion and Jaime de la Cruz Rodriguez. 
;
; AUTHOR:
;	  Gregal Vissers (g.j.m.vissers@astro.uio.no)
;	  @ Institute of Theoretical Astrophysics, University of Oslo
;   $Id$
;-

;========================= CRISPEX FUNCTIONS
;------------------------- APPLICATION USER DIRECTORY FUNCTION
FUNCTION CRISPEX_CONFIG_DIR
; Handles creation and update of ~/.idl/gvissers/crispex/ contents 
  ; Author details
  author_dirname = 'gvissers'
  author_desc = 'Gregal Vissers (g.j.m.vissers@astro.uio.no), '+$
    'Institute of Theoretical Astrophsyics, '+ $
    'University of Oslo'
  author_readme_version = 1

  ; Application details
  app_dirname = 'crispex'
  app_desc = 'CRISPEX: CRIsp SPectral EXplorer'
  app_readme_txt = [$
    'This is the user configuration directory for CRISPEX, ',$
    'written by Gregal Vissers (g.j.m.vissers@astro.uio.no), ',$
    'Institute of Theoretical Astrophysics, University of Oslo.',$
    'It is used to store user preferences and enable instance ',$
    'tracking. Both the directory and its contents can be safely ',$
    'deleted, however, any stored settings will then be lost.']
  app_readme_version = 1

  RETURN, APP_USER_DIR(author_dirname, author_desc, $
    app_dirname, app_desc, app_readme_txt, app_readme_version, $
    AUTHOR_README_VERSION=author_readme_version)

END

;------------------------- BINARY CONVERSION FUNCTION
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

;------------------------- BUTTON GROUP FUNCTIONS
FUNCTION CRISPEX_BGROUP_DIAGNOSTICS_SELECT, event
; Handles selection of diagnostics and calls necessary replot procedures
  WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (event.VALUE EQ 0) THEN BEGIN
    (*(*info).intparams).disp_diagnostics = REPLICATE(1,(*(*info).intparams).ndiagnostics)   
    (*(*info).intparams).ndisp_diagnostics = TOTAL((*(*info).intparams).disp_diagnostics)
    FOR i=0,(*(*info).intparams).ndiagnostics-1 DO $
      WIDGET_CONTROL, (*(*info).ctrlscp).specwin_button_ids[i+1], /SET_BUTTON, /SENSITIVE
  ENDIF ELSE BEGIN
    sel_idx = event.VALUE-1
    (*(*info).intparams).disp_diagnostics[sel_idx] = event.SELECT 
    (*(*info).intparams).ndisp_diagnostics = TOTAL((*(*info).intparams).disp_diagnostics)
    IF ((*(*info).intparams).ndisp_diagnostics EQ 1) THEN $
      WIDGET_CONTROL, (*(*info).ctrlscp).specwin_button_ids[$
        WHERE((*(*info).intparams).disp_diagnostics EQ 1)+1], SENSITIVE=0 $
    ELSE $
      FOR i=0,(*(*info).intparams).ndiagnostics-1 DO $
        WIDGET_CONTROL, (*(*info).ctrlscp).specwin_button_ids[i+1], /SENSITIVE
  ENDELSE
  WIDGET_CONTROL, (*(*info).ctrlscp).specwin_button_ids[0], $
    SENSITIVE=((*(*info).intparams).ndisp_diagnostics NE (*(*info).intparams).ndiagnostics), $
    SET_BUTTON=((*(*info).intparams).ndisp_diagnostics EQ (*(*info).intparams).ndiagnostics)
  ; Adjust slider settings based on available lp-range
  low_sel = (WHERE((*(*info).intparams).disp_diagnostics EQ 1, nwhereq1))[0]
  upp_sel = (WHERE((*(*info).intparams).disp_diagnostics EQ 1))[nwhereq1-1]
  lp_low = (*(*info).intparams).diag_start[low_sel]
  lp_upp = (*(*info).intparams).diag_start[upp_sel] + (*(*info).intparams).diag_width[upp_sel]-1
  IF (lp_low NE (*(*info).dispparams).lp_low) THEN $
    CRISPEX_DISPRANGE_LP_LOW, event, LP_SET=lp_low, /NO_DRAW
  IF (lp_upp NE (*(*info).dispparams).lp_upp) THEN $
    CRISPEX_DISPRANGE_LP_UPP, event, LP_SET=lp_upp, /NO_DRAW
  IF ((*(*info).intparams).ndisp_diagnostics NE (*(*info).intparams).ndisp_refdiagnostics) THEN BEGIN
    IF (*(*info).ctrlsswitch).lp_ref_lock THEN $
      CRISPEX_SLIDER_LP_REF_LOCK, event, /UNLOCK, /NO_DRAW
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SENSITIVE=0
  ENDIF ELSE $
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, /SENSITIVE
  CRISPEX_DRAW_GET_SPECTRAL_AXES, event, /MAIN
  CRISPEX_UPDATE_T, event
  CRISPEX_SCALING_APPLY_SELECTED, event
  IF (*(*info).winswitch).showsp THEN BEGIN
    CRISPEX_UPDATE_SPSLICE, event
    CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
  ENDIF
  IF (*(*info).winswitch).showphis THEN BEGIN
    CRISPEX_UPDATE_PHISLICE, event, /NO_DRAW
    CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event
  ENDIF
  CRISPEX_DRAW, event, /NO_REF
END

FUNCTION CRISPEX_BGROUP_REFDIAGNOSTICS_SELECT, event
; Handles selection of reference diagnostics and calls necessary replot procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info	
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (event.VALUE EQ 0) THEN BEGIN
    (*(*info).intparams).disp_refdiagnostics = REPLICATE(1,(*(*info).intparams).nrefdiagnostics)   
    (*(*info).intparams).ndisp_refdiagnostics = TOTAL((*(*info).intparams).disp_refdiagnostics)
    FOR i=0,(*(*info).intparams).nrefdiagnostics-1 DO $
      WIDGET_CONTROL, (*(*info).ctrlscp).refspecwin_button_ids[i+1], /SET_BUTTON, /SENSITIVE
  ENDIF ELSE BEGIN
    sel_idx = event.VALUE-1
    (*(*info).intparams).disp_refdiagnostics[sel_idx] = event.SELECT
    (*(*info).intparams).ndisp_refdiagnostics = TOTAL((*(*info).intparams).disp_refdiagnostics)
    IF ((*(*info).intparams).ndisp_refdiagnostics EQ 1) THEN $
      WIDGET_CONTROL, (*(*info).ctrlscp).refspecwin_button_ids[$
        WHERE((*(*info).intparams).disp_refdiagnostics EQ 1)+1], SENSITIVE=0 $
    ELSE $
      FOR i=0,(*(*info).intparams).nrefdiagnostics-1 DO $
        WIDGET_CONTROL, (*(*info).ctrlscp).refspecwin_button_ids[i+1], /SENSITIVE
  ENDELSE
  WIDGET_CONTROL, (*(*info).ctrlscp).refspecwin_button_ids[0], $
    SENSITIVE=((*(*info).intparams).ndisp_refdiagnostics NE (*(*info).intparams).nrefdiagnostics), $
    SET_BUTTON=((*(*info).intparams).ndisp_refdiagnostics EQ (*(*info).intparams).nrefdiagnostics)
  ; Adjust slider settings based on available lp-range
  low_sel = (WHERE((*(*info).intparams).disp_refdiagnostics EQ 1, nwhereq1))[0]
  upp_sel = (WHERE((*(*info).intparams).disp_refdiagnostics EQ 1))[nwhereq1-1]
  lp_low = (*(*info).intparams).refdiag_start[low_sel]
  lp_upp = (*(*info).intparams).refdiag_start[upp_sel] + (*(*info).intparams).refdiag_width[upp_sel]-1
  IF (lp_low NE (*(*info).dispparams).lp_ref_low) THEN $
    CRISPEX_DISPRANGE_LP_REF_RANGE, event, LP_LOW=lp_low, /NO_DRAW
  IF (lp_upp NE (*(*info).dispparams).lp_ref_upp) THEN $
    CRISPEX_DISPRANGE_LP_REF_RANGE, event, LP_UPP=lp_upp, /NO_DRAW
  IF ((*(*info).intparams).ndisp_diagnostics NE (*(*info).intparams).ndisp_refdiagnostics) THEN BEGIN
    IF (*(*info).ctrlsswitch).lp_ref_lock THEN $
      CRISPEX_SLIDER_LP_REF_LOCK, event, /UNLOCK, /NO_DRAW
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SENSITIVE=0
  ENDIF ELSE $
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, /SENSITIVE
  CRISPEX_DRAW_GET_SPECTRAL_AXES, event, /REFERENCE
  CRISPEX_UPDATE_T, event
  CRISPEX_SCALING_APPLY_SELECTED, event
  IF (*(*info).winswitch).showrefsp THEN BEGIN
    CRISPEX_UPDATE_REFSPSLICE, event
    CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
  ENDIF
  CRISPEX_DRAW, event, /NO_MAIN
END

FUNCTION CRISPEX_BGROUP_MASTER_TIME, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Get old time 
  CASE (*(*info).dispparams).master_time OF
    0:  t_old = (*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t_main]
    1:  t_old = (*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t_ref]
    2:  t_old = (*(*(*info).dispparams).tarr_sji)[(*(*info).dispparams).t_sji]
  ENDCASE
  ; Get new timing master
  (*(*info).dispparams).master_time = event.VALUE
  ; Reset timing offset to defaults
  (*(*info).dispparams).toffset_main = (*(*info).dataparams).default_toffset_main
  (*(*info).dispparams).toffset_ref = (*(*info).dataparams).default_toffset_ref
  CRISPEX_COORDS_TRANSFORM_T, event, T_OLD=t_old
  CRISPEX_UPDATE_T, event
  IF (*(*info).winswitch).showsp THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
  IF (*(*info).winswitch).showrefsp THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
	CRISPEX_DRAW, event
  WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_SLIDER_MIN=(*(*info).dispparams).t_low, $
    SET_SLIDER_MAX=((*(*info).dispparams).t_upp<(*(*info).dataparams).nt)
END

FUNCTION CRISPEX_BGROUP_STOKES_SELECT_SP, event, NO_DRAW=no_draw
  WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).stokesparams).prev_select_sp = (*(*info).stokesparams).select_sp
		(*(*info).stokesparams).select_sp[WHERE((*(*info).stokesparams).labels EQ $
                                      (*(*info).stokesparams).button_labels[event.VALUE])] = $
                                      event.SELECT
  sens = (TOTAL((*(*info).stokesparams).select_sp) NE 1)
	FOR i=0,TOTAL((*(*info).stokesparams).select_sp)-1 DO BEGIN
    FOR k=0,N_ELEMENTS((*(*info).stokesparams).button_labels)-1 DO BEGIN
		  IF (((*(*info).stokesparams).labels)[(WHERE((*(*info).stokesparams).select_sp EQ 1))[i]] EQ $
           (*(*info).stokesparams).button_labels[k]) THEN $
        WIDGET_CONTROL, (*(*info).ctrlscp).stokes_spbutton_ids[k], SENSITIVE = sens, /SET_BUTTON 
    ENDFOR
	ENDFOR
	CRISPEX_DISPLAYS_LS_RESIZE, event, /STOKES_SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [STRJOIN(((*(*info).stokesparams).labels)[$
      WHERE((*(*info).stokesparams).select_sp EQ 1)],', ')],labels=['Stokes detspect selected']
END

FUNCTION CRISPEX_BGROUP_STOKES_SELECT_XY, event, NO_DRAW=no_draw;, SET_STOKES=set_stokes
  WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
;  IF (N_ELEMENTS(SET_STOKES) NE 1) THEN idx = event.VALUE ELSE idx = set_stokes
  (*(*info).dataparams).s = $
    WHERE((*(*info).stokesparams).labels EQ (*(*info).stokesparams).button_labels[event.VALUE])
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [((*(*info).stokesparams).labels)[(*(*info).dataparams).s]],$
                         labels=['Stokes image selected']
  IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
;  	WIDGET_CONTROL, (*(*info).ctrlsparam).stokes_val, $
;      SET_VALUE = STRTRIM(((*(*info).stokesparams).labels)[(*(*info).dataparams).s],2)
    CRISPEX_SCALING_APPLY_SELECTED, event
  	CRISPEX_DISPLAYS_STOKES_SELECT_XY_RECOVER_YRANGE, event
  	CRISPEX_UPDATE_T, event
  	CRISPEX_UPDATE_SLICES, event
  	IF ((*(*info).winids).sptlb GT 0) THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
	  CRISPEX_DRAW, event
  ENDIF
END

FUNCTION CRISPEX_BGROUP_ZOOMFAC_SET, event, NO_DRAW=no_draw, NO_UPDATE_SLIDERS=no_update_sliders, $
                              SET_FACTOR_IDX=set_factor_idx, UNSET_FACTOR_IDX=unset_factor_idx
; Sets the zoomfactor and changes options and paramters accordingly
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(SET_FACTOR_IDX) EQ 1) THEN BEGIN
    idx = set_factor_idx   & select = 1B
    IF (N_ELEMENTS(UNSET_FACTOR_IDX) EQ 1) THEN (*(*info).zooming).factorswitch[unset_factor_idx]= 0B
  ENDIF ELSE BEGIN
    idx = event.VALUE & select = BYTE(event.SELECT)
  ENDELSE
  (*(*info).zooming).factorswitch[idx] = select
  (*(*info).zooming).factor = ((*(*info).zooming).factors)[idx]
;  IF ((*(*info).zooming).handle_extreme NE 2) THEN BEGIN
    IF ((*(*info).zooming).factor EQ 1) THEN BEGIN
      (*(*info).zooming).xpos = 0	&	(*(*info).zooming).ypos = 0
    ENDIF
;    sensitive = [0,0]
;  ENDIF ELSE $
;    sensitive = [(((*(*info).zooming).handle_extreme EQ 2) AND ((*(*info).data).ratio GT 1) OR $
;                  ((*(*info).zooming).factor NE 1)),$
;                 (((*(*info).zooming).handle_extreme EQ 2) AND ((*(*info).data).ratio LT 1) OR $
;                  ((*(*info).zooming).factor NE 1))]
  sensitive = [(((*(*info).dataparams).nx NE 1) AND ((*(*info).zooming).factor NE 1)),$
    ((*(*info).zooming).factor NE 1)]
  CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
  IF ~KEYWORD_SET(NO_UPDATE_SLIDERS) THEN $
  	CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y, SENSITIVE=sensitive
  FOR i=0,N_ELEMENTS((*(*info).zooming).factors)-1 DO $
    WIDGET_CONTROL, ((*(*info).ctrlscp).zoom_button_ids)[i], $
                    SET_BUTTON = ((*(*info).zooming).factorswitch)[i]
	CRISPEX_ZOOM, event, NO_DRAW=no_draw
END

FUNCTION CRISPEX_BGROUP_MASK_OVERLAY, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	maskim = (*(*info).overlayswitch).maskim
	maskim[event.value] = event.select
	(*(*info).overlayswitch).maskim = maskim
	CRISPEX_DRAW, event
END

;------------------------- READ BMP BUTTONS FUNCTION
FUNCTION CRISPEX_READ_BMP_BUTTONS, filename, srcdir
; Handles the reading of (button) BMP files
	button_dummy = READ_BMP(srcdir+filename)  
	button_dummy = TRANSPOSE(button_dummy, [1,2,0])
	RETURN, button_dummy
END

;------------------------- SCALING FUNCTIONS
FUNCTION CRISPEX_SCALING_CONTRAST, minimum_init, maximum_init, $
  minimum_perc, maximum_perc
  minimum_init = DOUBLE(minimum_init)
  maximum_init = DOUBLE(maximum_init)
  range = maximum_init-minimum_init
  minmax = [minimum_init+range*FLOAT(minimum_perc)/100.,$
            minimum_init+range*FLOAT(maximum_perc)/100.]
;  level = 0.5*(maximum-minimum)+minimum
;  width = (1.-contrast_value/101.)*(maximum-minimum)
;  minmax = [level-width/2.,level+width/2.]
  RETURN, minmax
END

FUNCTION CRISPEX_SCALING_SLICES, dispim, gamma_val, histo_opt_val, $
  default_min, default_max, FORCE_HISTO=force_histo
  IF (gamma_val NE 1) THEN BEGIN
    whereneg = WHERE(dispim LT 0, nwhereneg)
    IF (nwhereneg GT 0) THEN BEGIN
      dispim = (TEMPORARY(ABS(dispim)))^gamma_val
      dispim[whereneg] *= -1.
    ENDIF ELSE $
      dispim = (TEMPORARY((dispim)))^gamma_val
  ENDIF
  IF ((histo_opt_val NE 0) OR KEYWORD_SET(FORCE_HISTO)) THEN BEGIN
    IF (MIN(dispim, MAX=dispmax, /NAN) NE dispmax) THEN BEGIN
;      ; Copied over MISSING handling from modified HISTO_OPT() since it gives  "Floating illegal
;      ; operand" errors otherwise
;      finitvals = FINITE(dispim)
;      dispim[WHERE(finitvals eq 0)]=-32768
;      dispim=dispim[WHERE(dispim ne -32768)]
      dispim = HISTO_OPT(TEMPORARY(dispim), histo_opt_val, MISSING=-32768)
    ENDIF
  ENDIF
  minimum = MIN(dispim,MAX=maximum, /NAN)
  IF ((N_ELEMENTS(default_min) EQ 1) AND (N_ELEMENTS(default_max) EQ 1)) THEN $
    minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
      default_min, default_max) $
  ELSE $
    minmax = [minimum,maximum]
  RETURN, minmax
END

;========================= ABOUT WINDOW PROCEDURES
PRO CRISPEX_ABOUT_WINDOW, event 							
; Creates an about-window displaying code name, version and revision number
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

;================================================================================= CLEAR ESTIMATE PROCEDURES
PRO CRISPEX_CLEAR_CURRENT_ESTIMATE, event								
; Clears current saving time estimate
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	pftfiles = FILE_SEARCH((*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'cpft', COUNT = pftfilecount)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, ['crispex.'+(*(*info).paths).hostname+'cpft',pftfilecount], labels=['File to be deleted','Filecount']
	IF pftfilecount THEN BEGIN
		SPAWN,'rm '+(*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'cpft'
		(*(*info).feedbparams).estimate_lx = 0
		(*(*info).feedbparams).estimate_time = 0.
		(*(*info).feedbparams).estimate_run = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_estimate, SENSITIVE = 0
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','Could not delete crispex.'+$
      ((*(*info).paths).hostname)[0]+'cpft','from '+(*(*info).paths).dir_settings+'.','File does not exist.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END

PRO CRISPEX_CLEAR_CURRENT_INST, event								
; Clears current saving time estimate
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	instfiles = FILE_SEARCH((*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'inst', COUNT = instfilecount)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, ['crispex.'+(*(*info).paths).hostname+'inst',instfilecount], labels=['File to be deleted','Filecount']
	IF instfilecount THEN BEGIN
		SPAWN,'rm '+(*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'inst'
		(*(*info).feedbparams).estimate_lx = 0
		(*(*info).feedbparams).estimate_time = 0.
		(*(*info).feedbparams).estimate_run = 0
		WIDGET_CONTROL, (*(*info).ctrlscp).clear_current_inst, SENSITIVE = 0
	ENDIF ELSE BEGIN
		CRISPEX_WINDOW_OK, event,'ERROR!','Could not delete crispex.'+((*(*info).paths).hostname)[0]+$
      'inst','from '+(*(*info).paths).dir_settings+'.','File does not exist.',$
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
		(*(*info).winids).errtlb = tlb
	ENDELSE
END
;================================================================================= PROGRAM EXIT PROCEDURES
PRO CRISPEX_CLOSE, event								
; Called upon closing program, checks for existence of performance test file; if not present it is written
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).paths).dir_settings_write EQ 1) THEN BEGIN
		pftfiles = FILE_SEARCH((*(*info).paths).dir_settings+'crispex.'+$
      (*(*info).paths).hostname+'cpft', COUNT = pftfilecount)
		IF (pftfilecount EQ 0) AND ((*(*info).feedbparams).estimate_run EQ 1) THEN BEGIN
			estimate_lx = (*(*info).feedbparams).estimate_lx
			estimate_time = (*(*info).feedbparams).estimate_time
			estimate_run = (*(*info).feedbparams).estimate_run
			SAVE, estimate_lx, estimate_time, estimate_run, FILENAME = (*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'cpft'
			PRINT,'Written: '+(*(*info).paths).dir_settings+'crispex.'+(*(*info).paths).hostname+'cpft'
		ENDIF
	ENDIF ELSE BEGIN
		PRINT, 'ERROR: Could not write performance file crispex.'+(*(*info).paths).hostname+'cpft '
		PRINT, '       to '+(*(*info).paths).dir_settings+'. Permission denied.'
	ENDELSE
	FREE_LUN, (*(*info).data).lunim
	IF (*(*info).dataswitch).spfile THEN FREE_LUN, (*(*info).data).lunsp
	IF ((*(*info).winswitch).showref AND ((*(*info).data).lunrefim GT 0)) THEN FREE_LUN, (*(*info).data).lunrefim
	IF ((*(*info).dataswitch).refspfile AND ((*(*info).data).lunrefsp GT 0)) THEN FREE_LUN, (*(*info).data).lunrefsp
	IF ((*(*info).dataswitch).maskfile AND ((*(*info).data).lunmask GT 0)) THEN FREE_LUN, (*(*info).data).lunmask
	WIDGET_CONTROL, (*(*info).winids).root, /DESTROY
	PTR_FREE, info
END

PRO CRISPEX_CLOSE_CLEANUP, base								
; Clean-up upon closing program
	WIDGET_CONTROL, base, GET_UVALUE = info
	FREE_LUN, (*(*info).data).lunim
	IF (*(*info).dataswitch).spfile THEN FREE_LUN, (*(*info).data).lunsp
	IF ((*(*info).winswitch).showref AND ((*(*info).data).lunrefim GT 0)) THEN FREE_LUN, (*(*info).data).lunrefim
	IF ((*(*info).dataswitch).refspfile AND ((*(*info).data).lunrefsp GT 0)) THEN FREE_LUN, (*(*info).data).lunrefsp
	IF ((*(*info).dataswitch).maskfile AND ((*(*info).data).lunmask GT 0)) THEN FREE_LUN, (*(*info).data).lunmask
	CRISPEX_CLOSE_CLEAN_INSTANCE_FILE, (*(*info).paths).dir_settings_write, $
  (*(*info).paths).dir_settings, (*(*info).paths).hostname, ((*(*info).sesparams).curr_instance_id)[0]
	PTR_FREE, info
END

PRO CRISPEX_CLOSE_CLEAN_INSTANCE_FILE, dir_inst_write, dir_inst, hostname, curr_instance_id
; Called upon closing program, checks for existence of performance test file; if not present it is written
	IF (dir_inst_write EQ 1) THEN BEGIN
		instfile = FILE_SEARCH(dir_inst+'crispex.'+hostname+'inst', COUNT = instfilecount)
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  IF (event.TOP EQ (*(*info).winids).shorttlb) THEN (*(*info).winids).shorttlb = 0
  IF (event.TOP EQ (*(*info).winids).headertlb) THEN (*(*info).winids).headertlb = 0
	WIDGET_CONTROL, event.TOP, /DESTROY
END

;================================================================================= CURSOR PROCEDURES
PRO CRISPEX_CURSOR, event								
; Cursor handling procedure, tracks and handles events from the cursor on the main image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_TRACKING' THEN BEGIN
		IF event.ENTER THEN BEGIN
			WIDGET_CONTROL, event.HANDLER, get_value = wid
			WSET, wid
			ci = UINTARR(16) & cim = ci & cim[7] = 1
			DEVICE, CURSOR_IMAGE = ci, CURSOR_MASK = cim, CURSOR_XY = [8,8]
		ENDIF ELSE BEGIN
			IF (((*(*info).loopparams).np GE 1) AND (*(*info).overlayswitch).looppath_feedback AND ((*(*info).curs).lockset GT 0)) THEN BEGIN
				*(*(*info).loopparams).xp = (*(*(*info).loopparams).xp)[0:(*(*info).loopparams).np-1]
				*(*(*info).loopparams).yp = (*(*(*info).loopparams).yp)[0:(*(*info).loopparams).np-1]
				*(*(*info).overlayparams).sxp = (*(*(*info).overlayparams).sxp)[0:(*(*info).loopparams).np-1]
				*(*(*info).overlayparams).syp = (*(*(*info).overlayparams).syp)[0:(*(*info).loopparams).np-1]
				IF ((*(*info).loopparams).np GE 2) THEN CRISPEX_LOOP_GET_PATH, event ELSE BEGIN
					*(*(*info).loopparams).xr = *(*(*info).loopparams).xp
					*(*(*info).loopparams).yr = *(*(*info).loopparams).yp
				ENDELSE
				(*(*info).dataparams).x = (*(*(*info).loopparams).xp)[(*(*info).loopparams).np-1]
				(*(*info).dataparams).y = (*(*(*info).loopparams).yp)[(*(*info).loopparams).np-1]
				CRISPEX_COORDSLIDERS_SET, 0, 0, event
				IF (*(*info).winswitch).showphis THEN BEGIN
					(*(*info).curs).sx = (*(*info).curs).sxlock
					(*(*info).curs).sy = (*(*info).curs).sylock
					CRISPEX_PHISLIT_DIRECTION, event
          CRISPEX_UPDATE_PHISLIT_COORDS, event
					CRISPEX_UPDATE_PHISLICE, event
				ENDIF ELSE CRISPEX_DRAW, event
			ENDIF		
			IF (!D.WINDOW NE -1) THEN DEVICE, /CURSOR_CROSSHAIR 
		ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [event.ENTER], labels=['WIDGET_TRACKING: event.Enter']
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
					(*(*info).curs).xlock = FLOAT((*(*info).curs).sxlock * ((*(*info).dataparams).d_nx+1)) $
                                  / (*(*info).winsizes).xywinx + (*(*info).zooming).xpos
					(*(*info).curs).ylock = FLOAT((*(*info).curs).sylock * ((*(*info).dataparams).d_ny+1)) $
                                  / (*(*info).winsizes).xywiny + (*(*info).zooming).ypos
					(*(*info).dataparams).x = (*(*info).curs).xlock
					(*(*info).dataparams).y = (*(*info).curs).ylock
					IF (*(*info).overlayswitch).loopslit THEN BEGIN
						(*(*info).loopparams).np += 1
						IF ((*(*info).loopparams).np EQ 2) THEN $
              WIDGET_CONTROL, (*(*info).ctrlscp).loop_slit_but, SET_VALUE = 'Erase loop path'
						IF ((*(*info).loopparams).np GE 2) THEN BEGIN
							*(*(*info).loopparams).xp = [*(*(*info).loopparams).xp,(*(*info).curs).xlock]
							*(*(*info).loopparams).yp = [*(*(*info).loopparams).yp,(*(*info).curs).ylock]
							*(*(*info).overlayparams).sxp = [*(*(*info).overlayparams).sxp,(*(*info).curs).sxlock]
							*(*(*info).overlayparams).syp = [*(*(*info).overlayparams).syp,(*(*info).curs).sylock]
							IF ((*(*info).winids).looptlb EQ 0) THEN $
                WIDGET_CONTROL, (*(*info).ctrlscp).loop_slice_but, SENSITIVE = 1
							IF ((*(*info).loopparams).np EQ 3) THEN $
                WIDGET_CONTROL, (*(*info).ctrlscp).rem_loop_pt_but, SENSITIVE = 1
							CRISPEX_LOOP_GET, event
							CRISPEX_UPDATE_LP, event
							*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - $
                (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
							*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - $
                (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
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
				ENDIF ELSE IF ((*(*info).overlayswitch).loopslit AND $
                       (*(*info).overlayswitch).looppath_feedback AND $
                       ((*(*info).loopparams).np GE 1)) THEN BEGIN
					(*(*info).curs).sx = event.X
					(*(*info).curs).sy = event.Y
					CRISPEX_CURSOR_GET_XY, event
					*(*(*info).loopparams).xp = [(*(*(*info).loopparams).xp)[$
                                        0:(*(*info).loopparams).np-1],(*(*info).dataparams).x]
					*(*(*info).loopparams).yp = [(*(*(*info).loopparams).yp)[$
                                        0:(*(*info).loopparams).np-1],(*(*info).dataparams).y]
					*(*(*info).overlayparams).sxp = [(*(*(*info).overlayparams).sxp)[$
                                            0:(*(*info).loopparams).np-1],(*(*info).curs).sx]
					*(*(*info).overlayparams).syp = [(*(*(*info).overlayparams).syp)[$
                                            0:(*(*info).loopparams).np-1],(*(*info).curs).sy]
					CRISPEX_LOOP_GET_PATH, event
					IF ((*(*info).zooming).factor NE 1) THEN BEGIN
						*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) *$
              (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
						*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) *$
              (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
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
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, $
      [event.TYPE,event.PRESS,(*(*info).dataparams).x,(*(*info).dataparams).y,(*(*info).curs).sx,$
      (*(*info).curs).sy], labels=['WIDGET_DRAW: event.TYPE','WIDGET_DRAW: event.PRESS','x','y',$
      'sx','sy']
    IF (*(*info).winswitch).showsp THEN CRISPEX_UPDATE_SPSLICE, event
    IF (*(*info).winswitch).showrefsp THEN CRISPEX_UPDATE_REFSPSLICE, event
		IF (*(*info).winswitch).showphis THEN BEGIN
			CRISPEX_PHISLIT_DIRECTION, event
      CRISPEX_UPDATE_PHISLIT_COORDS, event
			CRISPEX_UPDATE_PHISLICE, event
		ENDIF ELSE CRISPEX_DRAW, event
	ENDIF
END

PRO CRISPEX_CURSOR_GET_XY, event
; Converts the window x and y coordinates to data x and y coordinates
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).x = (*(*info).curs).sx * ((*(*info).dataparams).d_nx+1) / $
                            (*(*info).winsizes).xywinx + (*(*info).zooming).xpos
	(*(*info).dataparams).y = (*(*info).curs).sy * ((*(*info).dataparams).d_ny+1) / $
                            (*(*info).winsizes).xywiny + (*(*info).zooming).ypos
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, $
    [(*(*info).dataparams).x,(*(*info).dataparams).y], labels=['x','y']
END

PRO CRISPEX_CURSOR_LOCK, event								
; Called upon locking/unlocking cursor with 'lock cursor' or 'unlock cursor' button, handles cursor (un)locking
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).curs).lockset = event.SELECT
	IF (*(*info).curs).lockset THEN BEGIN
		(*(*info).curs).xlock = (*(*info).dataparams).x	&	(*(*info).curs).ylock = (*(*info).dataparams).y
		(*(*info).curs).sxlock = (*(*info).curs).sx	    &	(*(*info).curs).sylock = (*(*info).curs).sy
		CRISPEX_COORDSLIDERS_SET, 0, 0, event
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
			CRISPEX_VERBOSE_GET, event, [(*(*info).curs).xlock,(*(*info).curs).ylock,(*(*info).curs).sxlock,(*(*info).curs).sylock], labels=['xlock','ylock','sxlock','sylock']
	ENDIF ELSE CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_COORDSLIDERS_SET, xsensitive, ysensitive, event				
; Adjusts sliders according to change in cursor position or locked/unlocked state
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).x_slider, SET_VALUE = (*(*info).dataparams).x, $
    SENSITIVE = (xsensitive AND ((*(*info).dispparams).x_first NE (*(*info).dispparams).x_last))
	WIDGET_CONTROL, (*(*info).ctrlscp).y_slider, SET_VALUE = (*(*info).dataparams).y, $
    SENSITIVE = (ysensitive AND ((*(*info).dispparams).y_first NE (*(*info).dispparams).y_last))
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x,(*(*info).dataparams).y,xsensitive,ysensitive], labels=['x','y','xsensitive','ysensitive']
END

;========================= COORDINATE TRANSFORMATIONS
PRO CRISPEX_COORDS_TRANSFORM_XY, event, MAIN2SJI=main2sji
; Handles transformation of x/y-coordinates in case of unequal spatial dimensions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF KEYWORD_SET(MAIN2SJI) THEN BEGIN
    (*(*info).dataparams).xsji = $
      (*(*info).dispparams).xyrastersji[(*(*info).dataparams).x,0] + $
      ((*(*info).dataparams).x * (*(*info).dataparams).dx / (*(*info).dataparams).sjidx - $
      ((*(*info).dispparams).xyrastersji[(*(*info).dataparams).x,0] - $
      (*(*info).dispparams).xyrastersji[0,0])) 
    (*(*info).dataparams).ysji = $
      (*(*info).dispparams).xyrastersji[(*(*info).dataparams).x,1] + $
      ((*(*info).dataparams).y * (*(*info).dataparams).dy / (*(*info).dataparams).sjidy - $
      ((*(*info).dispparams).xyrastersji[(*(*info).dataparams).x,1] - $
      (*(*info).dispparams).xyrastersji[0,1])) 
		sxsji = (*(*info).dataparams).xsji * (*(*info).winsizes).sjiwinx / $
          ((*(*info).dataparams).sjinx)
	  (*(*info).curs).sxsji = sxsji
		sysji = (*(*info).dataparams).ysji * (*(*info).winsizes).sjiwiny / $
          ((*(*info).dataparams).sjiny)
	  (*(*info).curs).sysji = sysji
  ENDIF
END

PRO CRISPEX_COORDS_TRANSFORM_T, event, T_OLD=t_old
; Handles transformation of t-coordinates in case of unequal temporal dimensions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
    ; Select temporal array for main and reference depending on dimensions and offset
    IF (SIZE((*(*info).dataparams).tarr_raster_main,/N_DIMENSIONS) EQ 2) THEN $
      tarr_main = REFORM((*(*info).dataparams).tarr_raster_main[$
        (*(*info).dispparams).toffset_main,*]) $
    ELSE $
      tarr_main = (*(*info).dataparams).tarr_raster_main
    IF (SIZE((*(*info).dataparams).tarr_raster_ref,/N_DIMENSIONS) EQ 2) THEN $
      tarr_ref = REFORM((*(*info).dataparams).tarr_raster_ref[$
        (*(*info).dispparams).toffset_ref,*]) $
    ELSE $
      tarr_ref = (*(*info).dataparams).tarr_raster_ref
    ; Set master array and temporal dimension depending on chosen master
    CASE (*(*info).dispparams).master_time OF
      0:  BEGIN
            tarr_master = tarr_main
            (*(*info).dataparams).nt = (*(*info).dataparams).mainnt
            offset_value = (*(*info).dispparams).toffset_main
          END
      1:  BEGIN
            tarr_master = tarr_ref
            (*(*info).dataparams).nt = (*(*info).dataparams).refnt
            offset_value = (*(*info).dispparams).toffset_ref
          END
      2:  BEGIN
            tarr_master = (*(*info).dataparams).tarr_sji
            (*(*info).dataparams).nt = (*(*info).dataparams).sjint
            offset_value = 0
          END
    ENDCASE
    ; Adjust time offset slider according to choices
    WIDGET_CONTROL, (*(*info).ctrlscp).time_offset_slider, $
      SENSITIVE=((*(*info).dispparams).master_time LT 2), SET_VALUE=offset_value
    ; Initialise variables
    tsel_main = LONARR((*(*info).dataparams).nt)
    IF ((*(*info).dataparams).refnt GT 1) THEN tsel_ref = LONARR((*(*info).dataparams).nt)
    IF ((*(*info).dataparams).sjint GT 1) THEN tsel_sji = LONARR((*(*info).dataparams).nt)
    ; Determine frame closest in time to master array
    FOR tt=0,(*(*info).dataparams).nt-1 DO BEGIN
      tdiff_main = ABS(tarr_main - tarr_master[tt])
      tsel_main[tt] = (WHERE(tdiff_main EQ MIN(tdiff_main, /NAN)))[0]
      IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN 
        tdiff_ref = ABS(tarr_ref - tarr_master[tt])
        tsel_ref[tt] = (WHERE(tdiff_ref EQ MIN(tdiff_ref, /NAN)))[0]
      ENDIF
      IF ((*(*info).dataparams).sjint GT 1) THEN BEGIN 
        tdiff_sji = ABS((*(*info).dataparams).tarr_sji - tarr_master[tt])
        tsel_sji[tt] = (WHERE(tdiff_sji EQ MIN(tdiff_sji, /NAN)))[0]
      ENDIF
    ENDFOR
    ; Populate variables with results
    *(*(*info).dispparams).tsel_main = tsel_main
    *(*(*info).dispparams).tarr_main = tarr_main[tsel_main]
    IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN 
      *(*(*info).dispparams).tsel_ref = tsel_ref
      *(*(*info).dispparams).tarr_ref = tarr_ref[tsel_ref]
    ENDIF
    IF ((*(*info).dataparams).sjint GT 1) THEN BEGIN 
      *(*(*info).dispparams).tsel_sji = tsel_sji
      *(*(*info).dispparams).tarr_sji = (*(*info).dataparams).tarr_sji[tsel_sji]
    ENDIF
    ; Reset temporal boundaries and get T_SET
    (*(*info).dispparams).t_last = (*(*info).dataparams).nt-1
    IF (N_ELEMENTS(T_OLD) EQ 1) THEN BEGIN
      tdiff = ABS(tarr_master - t_old)
      t_set = (WHERE(tdiff EQ MIN(tdiff, /NAN)))[0]
    ENDIF
    CRISPEX_DISPRANGE_T_RESET, event, /NO_DRAW, T_SET=t_set
END

;================================================================================= DISPLAYS PROCEDURES
PRO CRISPEX_DISPLAYS_ALL_TO_FRONT, event
; Brings all opened session windows to front
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info	
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	; Data windows
	IF ((*(*info).winids).sjitlb NE 0) THEN WSHOW, (*(*info).winids).sjiwid
;	WSHOW, (*(*info).winids).imwid
	IF ((*(*info).winids).sptlb NE 0) THEN WSHOW, (*(*info).winids).spwid
	IF ((*(*info).winids).lstlb NE 0) THEN WSHOW, (*(*info).winids).lswid
	IF ((*(*info).winids).phistlb NE 0) THEN WSHOW, (*(*info).winids).phiswid
	IF ((*(*info).winids).reftlb NE 0) THEN WSHOW, (*(*info).winids).refwid
	IF ((*(*info).winids).doptlb NE 0) THEN WSHOW, (*(*info).winids).dopwid
	IF ((*(*info).winids).imreftlb NE 0) THEN WSHOW, (*(*info).winids).imrefwid
	IF (TOTAL(*(*(*info).winids).restlooptlb) NE 0) THEN $
    FOR i=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO WSHOW, (*(*(*info).winids).restloopwid)[i]
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
  ; Control panel
  WIDGET_CONTROL, (*(*info).winids).root, /SHOW
END

PRO CRISPEX_DISPWIDS, event
; Brings all opened session windows to front
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info	
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).ctrlsswitch).imrefdetspect = 0
	CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
END

PRO CRISPEX_DISPLAYS_DETSPECT_REF_SELECT, event
; Handles the selection of detspect options for the reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).ctrlsswitch).imrefdetspect = 1
	CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
END

PRO CRISPEX_DISPLAYS_DETSPECT_SET_BUTTONS, event
; Handles the setting of scaling buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
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
;  IF ((*(*info).scaling).imrefscaling EQ 2) THEN BEGIN
;    CRISPEX_SCALING_SET_BUTTONS, event
;    CRISPEX_SCALING_SET_SLIDERS, event
;  ENDIF
	IF (*(*info).overlayswitch).mask THEN CRISPEX_MASK_BUTTONS_SET, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).doptlb,(*(*info).winids).dopwid,(*(*info).winids).dopdrawid], labels=['doptlb','dopwid','dopdrawid']
END

PRO CRISPEX_DISPLAYS_HEADER, event
; Pops up window with header information
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  files = ['Main: ','Reference: ','Slit-jaw image: ']
  hdrlen = [STRLEN((*(*(*(*info).dataparams).hdrs[0])[0])[0]), $
            STRLEN((*(*(*(*info).dataparams).hdrs[1])[0])[0]), $
            STRLEN((*(*(*(*info).dataparams).hdrs[2])[0])[0])]
  wherefileset = WHERE((hdrlen GT 0) EQ 1)
  FOR i=0,N_ELEMENTS(wherefileset)-1 DO BEGIN
    exte_idx = INDGEN((*(*info).dataparams).next[wherefileset[i]])
    tmp_vals = REPLICATE(files[wherefileset[i]],(*(*info).dataparams).next[wherefileset[i]])+$
      REPLICATE('Extension ',(*(*info).dataparams).next[wherefileset[i]])+$
      STRTRIM(exte_idx,2)
    tmp_uvals = [[REPLICATE(wherefileset[i],(*(*info).dataparams).next[wherefileset[i]])],[exte_idx]]
    IF (i EQ 0) THEN BEGIN
      vals = tmp_vals 
      uvals = tmp_uvals
    ENDIF ELSE BEGIN
      vals = [vals,tmp_vals]
      uvals = [uvals,tmp_uvals]
    ENDELSE
  ENDFOR
	base = WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+$
    ': File headers', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, $
    /TLB_KILL_REQUEST_EVENTS)
	disp = WIDGET_BASE(base, /COLUMN)
  (*(*info).ctrlshdr).header_select = WIDGET_COMBOBOX(disp, VALUE=vals, UVALUE=uvals, $
    EVENT_PRO='CRISPEX_DISPLAYS_HEADER_SELECT')
  text_base = WIDGET_BASE(disp, /COLUMN)
  (*(*info).ctrlshdr).header_txt = WIDGET_TEXT(text_base, $
    VALUE=(*(*(*(*info).dataparams).hdrs[0])[0]), XSIZE=CEIL(MAX(hdrlen)*1.1), $
    YSIZE=CEIL(MAX(hdrlen)/2.), /SCROLL, /WRAP)
  close_base = WIDGET_BASE(disp, /ALIGN_CENTER)
  close_button = WIDGET_BUTTON(close_base, VALUE='Close', EVENT_PRO='CRISPEX_CLOSE_EVENT_WINDOW')
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).lsxoffset, $
    TLB_SET_YOFFSET = (*(*info).winsizes).lswiny+1.5*(*(*info).winsizes).ydelta
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
  (*(*info).winids).headertlb = base
END

PRO CRISPEX_DISPLAYS_HEADER_SELECT, event
; Handles changing of header display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  WIDGET_CONTROL, (*(*info).ctrlshdr).header_select, GET_UVALUE=uvals
  WIDGET_CONTROL, (*(*info).ctrlshdr).header_txt, $
    SET_VALUE=(*(*(*(*info).dataparams).hdrs[uvals[event.INDEX,0]])[uvals[event.INDEX,1]])
END

PRO CRISPEX_DISPLAYS_INT_MENU, event, set_but_array
; Sets up the intensity-time plot options menu
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET = (*(*info).winsizes).spxoffset, $
    TLB_SET_YOFFSET = 0 
	WIDGET_CONTROL, base, SET_UVALUE = info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).intmenutlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).intmenutlb], labels=['intmenutlb']
END

PRO CRISPEX_DISPLAYS_INT_MENU_EVENT, event
; Handles the selection of diagnostics to be shown in the intensity-time plot
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	IF ( (*(*(*info).intparams).sel_diagnostics)[eventval] EQ 1) THEN (*(*(*info).intparams).selcol_diagnostics)[eventval] = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).intparams).selcol_diagnostics)[eventval]], labels=['Diagnostic ID','Color index selected']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_SEL_LINE, event
; Handles selection of linestyle of intensity versus time plot
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	IF ( (*(*(*info).intparams).sel_diagnostics)[eventval] EQ 1) THEN (*(*(*info).intparams).lines_diagnostics)[eventval] = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).intparams).lines_diagnostics)[eventval]], labels=['Diagnostic ID','Linestyle selected']
	CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DISPLAYS_INT_SEL_NONE, event
; Handles selection of none intensity versus time plots
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN
		CRISPEX_DISPLAYS_LOOPSLAB, event,/NO_DRAW
		CRISPEX_DISPLAYS_REFLOOPSLAB, event
	ENDIF ELSE CRISPEX_DISPLAYS_LOOPSLAB, event
END

PRO CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event					
; Updates loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).loopwid
;  t_low_y = (*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t_low]
;  tarr_main_sel = (*(*(*info).dispparams).tarr_main)[$
;    (*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
;  t_upp_y = tarr_main_sel[(WHERE(tarr_main_sel NE 0))[-1]]    ; Temp fix for tarr[-1]=0
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), *(*(*info).dispparams).tarr_main, $
    /NODATA, YR=[(*(*info).dispparams).t_low_main,(*(*info).dispparams).t_upp_main], /YS, $
    POS = [(*(*info).plotpos).loopx0,(*(*info).plotpos).loopy0,$
    (*(*info).plotpos).loopx1,(*(*info).plotpos).loopy1], YTICKLEN=(*(*info).plotaxes).loopyticklen,$
    XTICKLEN=(*(*info).plotaxes).loopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, $
    XTITLE='Pixel along loop', BACKGROUND = (*(*info).plotparams).bgplotcol, $
    COLOR=(*(*info).plotparams).plotcol
  IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).loopwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_LOOPSLAB_RESIZE, event
; Loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).winswitch).showloop = 1
	WIDGET_CONTROL,/HOURGLASS
	CRISPEX_LOOP_GET_PATH, event
	CRISPEX_LOOP_GET_SLAB, event
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': T-slice along loop'
	CRISPEX_WINDOW, (*(*info).winsizes).loopxres, (*(*info).winsizes).loopyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta, $
		((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = loopdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_LOOPSLAB_RESIZE'
;  t_low_y = (*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t_low]
;  tarr_main_sel = (*(*(*info).dispparams).tarr_main)[$
;    (*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
;  t_upp_y = tarr_main_sel[(WHERE(tarr_main_sel NE 0))[-1]]    ; Temp fix for tarr[-1]=0
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), *(*(*info).dispparams).tarr_main, $
    /NODATA, YR=[(*(*info).dispparams).t_low_main,(*(*info).dispparams).t_upp_main], /YS, $
    POS=[(*(*info).plotpos).loopx0,(*(*info).plotpos).loopy0,$
    (*(*info).plotpos).loopx1,(*(*info).plotpos).loopy1],$
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).refloopwid
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), *(*(*info).dispparams).tarr_ref, $
    /NODATA, /YS, POS=[(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopy0, $
    (*(*info).plotpos).refloopx1,(*(*info).plotpos).refloopy1], $
		YTICKLEN = (*(*info).plotaxes).refloopyticklen, XTICKLEN = (*(*info).plotaxes).refloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refloopwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_REFLOOPSLAB_RESIZE, event
; Reference loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).winswitch).showrefloop = 1
	WIDGET_CONTROL,/HOURGLASS
	CRISPEX_LOOP_GET_REFSLAB, event		
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Reference T-slice along loop'
	CRISPEX_WINDOW, (*(*info).winsizes).refloopxres, (*(*info).winsizes).refloopyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).xywinx+(*(*info).winsizes).xdelta, $
		((*(*info).winswitch).showsp + (*(*info).winswitch).showphis) * (*(*info).winsizes).ydelta, DRAWID = refloopdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFLOOPSLAB_RESIZE'
	PLOT, FINDGEN((*(*info).loopsdata).loopsize), *(*(*info).dispparams).tarr_ref, $
    /NODATA, YR=[(*(*info).dispparams).t_low_ref, (*(*info).dispparams).t_upp_ref], $
		/YS, POS=[(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopy0,$
    (*(*info).plotpos).refloopx1,(*(*info).plotpos).refloopy1],$
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

PRO CRISPEX_DISPLAYS_PLOT_RESIZE, event, new_xres_tmp, new_yres_tmp, init_xres, init_yres, init_xmargin, init_xwall, new_xres, new_yres, new_width, new_height, $
	x0, x1, y0, y1, v_dop_set, INX0=inx0, INX1=inx1, INY0=iny0, INY1=iny1, ERROR=error, GOLDEN=golden, ACTUAL_RESIZE=actual_resize, DETSPECT=detspect, STOKES_SELECT=stokes_select
; Handles the display plot resizing
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_DRAW_SPECTRAL_MAIN, event, /LS_ONLY
END

PRO CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE, event
; Sets the playback mode to blink of main and reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	IF (*(*info).ctrlsswitch).imrefdetspect THEN BEGIN	; For reference detailed spectrum window
		IF ((*(*info).winswitch).showrefls EQ 0) THEN BEGIN
			title = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+((*(*info).plottitles).reflswintitle)[(*(*info).plotswitch).refheightset]
			CRISPEX_WINDOW, (*(*info).winsizes).reflsxres, (*(*info).winsizes).reflsyres, (*(*info).winids).root, title, tlb, wid, (*(*info).winsizes).lsxoffset, $
				(*(*info).winswitch).showsp * (*(*info).winsizes).ydelta, DRAWID = lsdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFLS_RESIZE'
			(*(*info).winids).reflstlb = tlb		&	(*(*info).winids).reflswid = wid	&	(*(*info).winswitch).showrefls = 1
			(*(*info).winids).reflsdrawid = lsdrawid	&	(*(*info).winids).reflswintitle = title
			WIDGET_CONTROL, (*(*info).winids).reflstlb, SET_UVALUE = info
			IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_SPECTRAL_REF, event, /LS_ONLY
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
			IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_SPECTRAL_MAIN, event, /LS_ONLY
		ENDIF ELSE BEGIN
			WIDGET_CONTROL, (*(*info).winids).lstlb, /DESTROY
			(*(*info).winids).lstlb = 0
			(*(*info).winswitch).showls = 0
		ENDELSE
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).lstlb,(*(*info).winids).lswid,(*(*info).winids).lsdrawid], labels=['lstlb','lswid','lsdrawid']
	ENDELSE
END

PRO CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event, NO_AXES=no_axes
; Updates temporal spectrum display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).phiswid
	IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
	IF (*(*info).plotswitch).multichannel THEN $
    title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle $
  ELSE $
    title = ''
  ytitle = 'Position along slit [pixel]'
  ; Set axes parameters depending on # of diagnostics
  IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
    xticklen = 1E-9   & xtitle = ''
    xtickname = REPLICATE(' ',60)
    IF KEYWORD_SET(NO_AXES) THEN BEGIN
      yticklen = 1E-9 & ytitle = ''
      ytickname = REPLICATE(' ',60)
    ENDIF ELSE BEGIN
      yticklen = (*(*info).plotaxes).phisyticklen
      ytickname = ''
    ENDELSE
  ENDIF ELSE BEGIN
    xticklen = (*(*info).plotaxes).phisxticklen
    xtitle = (*(*info).plottitles).spxtitle
    yticklen = (*(*info).plotaxes).phisyticklen
    xtickname = ''  & ytickname = ''
  ENDELSE
  topxtitle = title
  IF (*(*info).plotswitch).v_dop_set THEN topxtitle += 'Doppler velocity [km/s]'
  ; Determine plot ranges
    (*(*info).plotaxes).phis_yrange = $
              [-(((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + $
                 (*(*info).phiparams).sphi)-0.5, (*(*info).phiparams).nw_cur - $
               (((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + $
                 (*(*info).phiparams).sphi)-0.5]
    phis_xrange = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low],$
          (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]]
  ; Plot basic axes box
	PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).phiparams).nw_cur), /NODATA, $
    YRANGE = (*(*info).plotaxes).phis_yrange, /YS, XRANGE=phis_xrange, $
		XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1,  $
    YTICKLEN = yticklen, YTITLE=ytitle, YTICKNAME=ytickname, $
    XTICKLEN = xticklen, XTITLE=xtitle, XTICKNAME=xtickname, $
		POS = [(*(*info).plotpos).phisx0,(*(*info).plotpos).phisy0,$
             (*(*info).plotpos).phisx1,(*(*info).plotpos).phisy1], $
    BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol;, $
    NOERASE=((*(*info).intparams).ndiagnostics GT 1)
  IF ~KEYWORD_SET(NO_AXES) THEN BEGIN
    ; Plot xtitle(s)
    IF ((*(*info).intparams).ndiagnostics GT 1) THEN $
        XYOUTS,(*(*info).plotpos).phisxplspw/2.+(*(*info).plotpos).phisx0,$
          (*(*info).plotpos).phisy0/3.,(*(*info).plottitles).spxtitle,ALIGNMENT=0.5, $
          COLOR = (*(*info).plotparams).plotcol,/NORMAL
    ; Loop over all diagnostics for plotting of detailed spectrum
    FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
      disp_idx = (WHERE((*(*info).intparams).disp_diagnostics EQ 1))[d]
      ; Determine xrange to display
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        phis_xrange = (*(*info).dataparams).lps[[(*(*info).intparams).diag_start[disp_idx],$
          ((*(*info).intparams).diag_start[disp_idx]+(*(*info).intparams).diag_width[disp_idx]-1)]]
        vdop_xrange = (*(*(*info).plotaxes).v_dop[disp_idx])[$
          [0,(*(*(*info).intparams).diag_widths)[d]-1]]
      ENDIF ELSE $
        vdop_xrange = (*(*(*info).plotaxes).v_dop[0])[$
          [(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp]]
      IF (d EQ 0) THEN offset = 0 ELSE offset = TOTAL((*(*(*info).plotaxes).diag_range_phis)[0:(d-1)])
      ; Determine lower left corner position of plot
      phisx0 = (*(*info).plotpos).phisx0 + offset
      phisx1 = (*(*(*info).plotaxes).diag_range_phis)[d] + phisx0
  	  PLOT, (*(*info).dataparams).lps, FINDGEN((*(*info).phiparams).nw_cur), /NODATA, $
        YRANGE=(*(*info).plotaxes).phis_yrange, /YS, XRANGE=phis_xrange, $
        XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, $
        POS = [phisx0,(*(*info).plotpos).phisy0,phisx1,(*(*info).plotpos).phisy1], $
        YTICKLEN=1E-9, XTICKLEN=(*(*info).plotaxes).phisxticklen, $
        YTICKNAME=REPLICATE(' ',60), $;XTICKNAME=xtickname, $
        XTICKINTERVAL=(*(*info).plotaxes).xtickinterval, $
        BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol,$
        /NOERASE
      ; Display Doppler top axis if Doppler set, else regular top axis
  		IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN BEGIN
        IF (d EQ 0) THEN $
          XYOUTS,((*(*info).plotpos).phisx1-(*(*info).plotpos).phisx0)/2.+$
            (*(*info).plotpos).phisx0,$
            (*(*info).plotpos).phisy0/5.*3+(*(*info).plotpos).phisy1, topxtitle, $
            ALIGNMENT=0.5, COLOR = (*(*info).plotparams).plotcol,/NORMAL
        topxtitle = ''
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, $
          XRANGE=vdop_xrange, XSTYLE=1, $
          COLOR = (*(*info).plotparams).plotcol,$
          XTICKINT=(*(*info).plotaxes).xdoptickinterval
      ENDIF ELSE $
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).phisxticklen, XRANGE = xrange, XSTYLE=1, $
          XTITLE = topxtitle, COLOR = (*(*info).plotparams).plotcol, $
          XTICKNAME=xtickname
    ENDFOR
  ENDIF
END

PRO CRISPEX_DISPLAYS_PHIS_RESIZE, event							
; Spectral phi slice window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_DRAW_GET_SPECTRAL_AXES, event, /MAIN
  CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event
	CRISPEX_DRAW_PHIS, event
END

PRO CRISPEX_DISPLAYS_PHIS_TOGGLE, event, NO_DRAW=no_draw
; Spectral phi slice window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	IF ((*(*info).winswitch).showphis EQ 0) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		wintitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+$
      ((*(*info).plottitles).phiswintitle)[(*(*info).plotswitch).heightset]
		CRISPEX_WINDOW, (*(*info).winsizes).phisxres, (*(*info).winsizes).phisyres, $
      (*(*info).winids).root, wintitle, tlb, wid, (*(*info).winsizes).spxoffset, $
      (*(*info).winswitch).showsp * (*(*info).winsizes).ydelta, DRAWID = phisdrawid, $
      RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_PHIS_RESIZE'
		(*(*info).winids).phistlb = tlb			&	(*(*info).winids).phiswid = wid
		(*(*info).winids).phisdrawid = phisdrawid	&	(*(*info).winids).phiswintitle = wintitle 
		WIDGET_CONTROL, (*(*info).winids).phistlb, SET_UVALUE = info
		(*(*info).ctrlsswitch).bwd_insensitive = 0	
		(*(*info).ctrlsswitch).fwd_insensitive = 0
		(*(*info).winswitch).showphis = 1
    CRISPEX_DISPRANGE_LP_RANGE, event, /NO_DRAW
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).phistlb, /DESTROY
		(*(*info).winids).phistlb = 0
		(*(*info).winswitch).showphis = 0
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).phi_slider, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).nphi_slider, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).bwd_move_slit, SENSITIVE = (*(*info).winswitch).showphis
	WIDGET_CONTROL, (*(*info).ctrlscp).fwd_move_slit, SENSITIVE = (*(*info).winswitch).showphis
  IF (*(*info).winswitch).showphis THEN BEGIN
  	CRISPEX_PHISLIT_DIRECTION, event
    CRISPEX_UPDATE_PHISLIT_COORDS, event
    CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event
    CRISPEX_UPDATE_SLICES, event, NO_DRAW=no_draw
  ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).phistlb,(*(*info).winids).phiswid,$
      (*(*info).winids).phisdrawid], labels=['phistlb','phiswid','phisdrawid']
END

PRO CRISPEX_DISPLAYS_REF_TOGGLE, event, NO_DRAW=no_draw
; Reference image window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	(*(*info).winswitch).showref = event.SELECT
	IF (*(*info).winswitch).showref THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Reference image'
		CRISPEX_WINDOW, (*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny, (*(*info).winids).root,$
                    title, reftlb, refwid, (*(*info).winsizes).refxoffset,$
                    (*(*info).winsizes).refyoffset, DRAWID = refdrawid, DRAWBASE = refdrawbase
		(*(*info).winids).reftlb = reftlb		&	(*(*info).winids).refwid = refwid	
    (*(*info).winids).refdrawid = refdrawid
		(*(*info).winids).refdrawbase = refdrawbase	&	(*(*info).winids).refwintitle = title
		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
			CRISPEX_UPDATE_T, event
			CRISPEX_DRAW_REF, event
		ENDIF
		WIDGET_CONTROL, refdrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, $
                    /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS
		WIDGET_CONTROL, reftlb, SET_UVALUE = info
		XMANAGER, 'CRISPEX', reftlb, /NO_BLOCK
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).reftlb, /DESTROY
		(*(*info).winids).reftlb = 0
	ENDELSE
;  IF ((*(*info).scaling).imrefscaling EQ 1) THEN BEGIN
;    CRISPEX_SCALING_SET_BUTTONS, event
;    CRISPEX_SCALING_SET_SLIDERS, event
;  ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SENSITIVE = $
    (((*(*info).dataparams).nlp EQ (*(*info).dataparams).refnlp) AND ((*(*info).dataparams).refnlp GT 1))
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, $
    SENSITIVE = (((*(*info).dataparams).refnlp GT 1) AND ABS((*(*info).ctrlsswitch).lp_ref_lock-1))
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).reftlb,(*(*info).winids).refwid,$
                         (*(*info).winids).refdrawid], labels=['reftlb','refwid','refdrawid']
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	FOR i=0,N_ELEMENTS(*(*(*info).winids).restloopwid)-1 DO BEGIN
		WSET, (*(*(*info).winids).restloopwid)[i]
		PLOT, FINDGEN((*(*(*info).loopsdata).rest_loopsize)[i]), *(*(*info).dispparams).tarr_main, $
      /NODATA, YR=[(*(*info).dispparams).t_low_main, (*(*info).dispparams).t_upp_main], $
			/YS, POS=[(*(*info).plotpos).restloopx0,(*(*info).plotpos).restloopy0,	$
      (*(*info).plotpos).restloopx1,(*(*info).plotpos).restloopy1], $
			YTICKLEN = (*(*info).plotaxes).restloopyticklen, XTICKLEN = (*(*info).plotaxes).restloopxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
			BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).winids).restloopwid)[i]], labels=['Window ID for replot']
	ENDFOR
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_RESIZE, event				
; Restored loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
			refbase = FILE_BASENAME(STRMID((*(*info).dataparams).refimfilename,0,$
                  STRPOS((*(*info).dataparams).refimfilename,'.',/REVERSE_SEARCH)))
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
END

PRO CRISPEX_DISPLAYS_RESTORE_LOOPSLAB, event, NO_DRAW=no_draw, INDEX=index
; Restored loopslab display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
		update_t_range = 0
		update_lp_range = 0
		restricted_t_range = 0
		restricted_lp_range = 0
		restricted_lp_ref_range = 0
		ref_slice_only = 0
		WIDGET_CONTROL,/HOURGLASS
		RESTORE, (*(*info).restoreparams).disp_loopfile
		IF (N_ELEMENTS(INDEX) EQ 1) THEN $
      idx = index $
    ELSE $
      idx = N_ELEMENTS(*(*(*info).restoreparams).disp_loopnr)-1
		*(*(*(*info).loopsdata).rest_crossloc[idx]) = vertices
		IF (N_ELEMENTS(loop_slab) GT 0) THEN $
      loopslab = loop_slab $
    ELSE $
      loopslab = loop_slice
		*(*(*(*info).loopsdata).rest_loopslab[idx]) = loopslab
		slice_only = (SIZE(*(*(*(*info).loopsdata).rest_loopslab[idx]),/N_DIMENSIONS) LT 3)
		IF slice_only THEN BEGIN			; Only a slice
			*(*(*(*info).loopsdata).rest_loopslice[idx]) = *(*(*(*info).loopsdata).rest_loopslab[idx])
			IF ((SIZE(*(*(*(*info).loopsdata).rest_loopslice[idx])))[2] NE $
        (*(*info).dataparams).nt) THEN BEGIN
				IF (N_ELEMENTS(t_low) EQ 0) THEN BEGIN
					CRISPEX_WINDOW_OK, event, 'WARNING!', $
						'The slice to be loaded has a reduced temporal range,',$
            'however the format in which it was saved does not',$
						'allow for correct slice restoration. If display',$
            'is required, please save the slice again.', $
            OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
					(*(*info).winids).warntlb = tlb
					RETURN
				ENDIF ELSE BEGIN
					IF ((t_low GE (*(*info).dispparams).t_upp) OR $
              (t_upp LT (*(*info).dispparams).t_low)) THEN BEGIN
						CRISPEX_WINDOW_OK, event, 'WARNING!', $
							'The temporal range of the loaded slice falls outside',$
              'the range set by the currently loaded slices. If',$
							'display is required, please close all currently',$
              'loaded slices before proceeding.', $
              OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
						(*(*info).winids).warntlb = tlb
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
						caseidx = (*(*(*info).restoreparams).disp_loopnr)[idx]
						list_values[caseidx+1] = 'Display time slice '+STRTRIM(caseidx,2)
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, $
              SET_COMBOBOX_SELECT = caseidx+1
						*(*(*info).restoreparams).disp_loopnr = (*(*(*info).restoreparams).disp_loopnr)[0:idx-1]
						*(*(*info).restoreparams).disp_imref = (*(*(*info).restoreparams).disp_imref)[0:idx-1]
						RETURN
					ENDIF ELSE BEGIN
						update_t_range = 1
						restricted_t_range = 1
						(*(*info).dispparams).t_low = t_low
						(*(*info).dispparams).t_upp = t_upp
						IF (N_ELEMENTS(t_saved) EQ 0) THEN $
              (*(*info).dispparams).t = (*(*info).dispparams).t_low $
            ELSE $
              (*(*info).dispparams).t = t_saved
						WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, $
              SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, $
              SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = 0
					ENDELSE
				ENDELSE
			ENDIF
			IF (*(*(*info).restoreparams).disp_imref)[idx] THEN BEGIN
				restricted_lp_ref_range = 1
				(*(*info).dataparams).lp_ref = spect_pos > (*(*info).dispparams).lp_ref_low < (*(*info).dispparams).lp_ref_upp 
				WIDGET_CONTROL,(*(*info).ctrlscp).lp_ref_slider, SENSITIVE = 0, $
          SET_VALUE = (*(*info).dataparams).lp_ref
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
				IF (N_ELEMENTS(t_saved) EQ 0) THEN $
          (*(*info).dispparams).t = (*(*info).dispparams).t_low $
        ELSE $
          (*(*info).dispparams).t = t_saved
				IF (N_ELEMENTS(t_low) EQ 0) THEN BEGIN
					CRISPEX_WINDOW_OK, event, 'WARNING!', $
						'The slice to be loaded has a reduced temporal range,',$
            'however the format in which it was saved does not',$
						'allow for correct slice restoration. If display',$
            'is required, please save the slice again.', $
            OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
					(*(*info).winids).warntlb = tlb
					RETURN
				ENDIF ELSE BEGIN
					IF ((t_low GE (*(*info).dispparams).t_upp) OR $
            (t_upp LT (*(*info).dispparams).t_low)) THEN BEGIN
						CRISPEX_WINDOW_OK, event, 'WARNING!', $
							'The temporal range of the loaded slice falls outside',$
              'the range set by the currently loaded slices. If',$
							'display is required, please close all currently',$
              'loaded slices before proceeding.', $
              OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', BASE=tlb
						(*(*info).winids).warntlb = tlb
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, GET_VALUE = list_values
						caseidx = (*(*(*info).restoreparams).disp_loopnr)[idx]
						list_values[caseidx+1] = 'Display time slice '+STRTRIM(caseidx,2)
						WIDGET_CONTROL, (*(*info).ctrlsrestore).disp_list, SET_VALUE = list_values, $
              SET_COMBOBOX_SELECT = caseidx+1
						*(*(*info).restoreparams).disp_loopnr = (*(*(*info).restoreparams).disp_loopnr)[0:idx-1]
						*(*(*info).restoreparams).disp_imref = (*(*(*info).restoreparams).disp_imref)[0:idx-1]
						RETURN
					ENDIF ELSE BEGIN
						update_t_range = 1
						restricted_t_range = 1
						(*(*info).dispparams).t_low = t_low
						(*(*info).dispparams).t_upp = t_upp
						IF (N_ELEMENTS(t) EQ 0) THEN $
              (*(*info).dispparams).t = (*(*info).dispparams).t_low $
            ELSE $
              (*(*info).dispparams).t = t
						WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, $
              SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE = 0
						WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, $
              SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE = 0
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
		PLOT, FINDGEN((*(*(*info).loopsdata).rest_loopsize)[idx]), *(*(*info).dispparams).tarr_main, $
      /NODATA, YR=[(*(*info).dispparams).t_low_main, (*(*info).dispparams).t_upp_main], $
			/YS, POS=[(*(*info).plotpos).restloopx0,(*(*info).plotpos).restloopy0,$
      (*(*info).plotpos).restloopx1,(*(*info).plotpos).restloopy1], $
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).retrdetwid
	PLOT, FINDGEN((*(*info).loopsdata).det_loopsize), *(*(*info).dispparams).tarr_main, $
    /NODATA, YR=[(*(*info).dispparams).t_low_main, (*(*info).dispparams).t_upp_main], $
    /YS, POS=[(*(*info).plotpos).retrdetx0,(*(*info).plotpos).retrdety0,$
    (*(*info).plotpos).retrdetx1,(*(*info).plotpos).retrdety1], $
		YTICKLEN = (*(*info).plotaxes).retrdetyticklen, XTICKLEN = (*(*info).plotaxes).retrdetxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
		BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).retrdetwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_RESIZE, event	
; Retrieved detection loopslab window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
		PLOT, FINDGEN((*(*info).loopsdata).det_loopsize), *(*(*info).dispparams).tarr_main, $
      /NODATA, YR=[(*(*info).dispparams).t_low_main, (*(*info).dispparams).t_upp_main], $
      /YS, POS=[(*(*info).plotpos).retrdetx0,(*(*info).plotpos).retrdety0,$
      (*(*info).plotpos).retrdetx1,(*(*info).plotpos).retrdety1], $
			YTICKLEN = (*(*info).plotaxes).retrdetyticklen, XTICKLEN = (*(*info).plotaxes).retrdetxticklen, /XS, YTITLE = (*(*info).plottitles).spytitle, XTITLE = 'Pixel along loop', $
			BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol
		(*(*info).winids).retrdettlb = tlb	&	(*(*info).winids).retrdetwid = wid	&	(*(*info).winids).retrdetdrawid = disp_retr_detdrawid
		(*(*info).winids).retrdetwintitle = title
		WIDGET_CONTROL, (*(*info).winids).retrdettlb, SET_UVALUE = info
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DISPRANGE_T_RANGE, event
		WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, $
      SET_VALUE=STRTRIM((*(*info).dispparams).t_low,2), SENSITIVE=0
		WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, $
      SET_VALUE=STRTRIM((*(*info).dispparams).t_upp,2), SENSITIVE=0 
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_DRAW_SPECTRAL_REF, event, /LS_ONLY
END

PRO CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event, NO_AXES=no_axes
; Updates reference temporal spectrum display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).refspwid
  ; Set axes parameters depending on # of diagnostics
  IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
    xticklen = 1E-9   & xtitle = ''
    xtickname = REPLICATE(' ',60)
    IF KEYWORD_SET(NO_AXES) THEN BEGIN
      yticklen = 1E-9 & ytitle = ''
      ytickname = REPLICATE(' ',60)
    ENDIF ELSE BEGIN
      yticklen = (*(*info).plotaxes).refspyticklen
      ytitle = (*(*info).plottitles).spytitle
      ytickname = ''
    ENDELSE
  ENDIF ELSE BEGIN
    xticklen = (*(*info).plotaxes).refspxticklen
    xtitle = (*(*info).plottitles).refspxtitle
    yticklen = (*(*info).plotaxes).refspyticklen
    ytitle = (*(*info).plottitles).spytitle
    xtickname = ''  & ytickname = ''
  ENDELSE
  t_low_y = (*(*info).dispparams).t_low_ref
  ;(*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t_low]
  t_upp_y = (*(*info).dispparams).t_upp_ref
  ;(*(*(*info).dispparams).tarr_ref)[$
  ;            (WHERE(*(*(*info).dispparams).tarr_ref NE 0))[-1]]  ; Temp fix for tarr[-1]=0
  xrange = [(*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low], $
            (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_upp]]
  correct_axes = (FLOOR(ALOG10(ABS(t_upp_y))) GE 3)
  IF ((t_low_y NE 0.) AND ~KEYWORD_SET(correct_axes)) THEN $
    correct_axes = (FLOOR(ALOG10(ABS(t_low_y))) LE -2) 
  IF correct_axes THEN BEGIN
 	  order_corr = FLOOR(ALOG10(ABS(t_upp_y)))
 	  IF ~KEYWORD_SET(NO_AXES) THEN ytitle += ' (x10!U'+STRTRIM(order_corr,2)+'!N)'
    t_low_y /= (10.^(order_corr))
    t_upp_y /= (10.^(order_corr))
  ENDIF 
  IF (*(*info).plotswitch).v_dop_set_ref THEN topxtitle = 'Doppler velocity [km/s]'
  ; Plot basic axes box
  PLOT, (*(*info).dataparams).lps, *(*(*info).dispparams).tarr_ref, $
    YR = [t_low_y,t_upp_y], /YS, XR=xrange, $
    XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, $
  	YTICKLEN = yticklen, YTITLE = ytitle, YTICKNAME = ytickname, $
    XTICKLEN = xticklen, XTITLE = xtitle, XTICKNAME = xtickname, $
    POS = [(*(*info).plotpos).refspx0,(*(*info).plotpos).refspy0,(*(*info).plotpos).refspx1,$
    (*(*info).plotpos).refspy1], BACKGROUND=(*(*info).plotparams).bgplotcol, $
    COLOR=(*(*info).plotparams).plotcol, /NODATA, NOERASE=KEYWORD_SET(NO_AXES)
  ; In case of multiple diagnostics and regular replotting of axes
  ; Determine proportional spectral window sizes
  diag_widths = (*(*info).intparams).refdiag_width[$
    WHERE((*(*info).intparams).disp_refdiagnostics EQ 1)]
  diag_ratio = diag_widths / FLOAT(TOTAL(diag_widths))
  diag_range = diag_ratio * (*(*info).plotpos).refxplspw
  IF ~KEYWORD_SET(NO_AXES) THEN BEGIN
    ; Plot xtitle(s)
    IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN $
      XYOUTS,(*(*info).plotpos).refxplspw/2.+(*(*info).plotpos).refspx0,$
        (*(*info).plotpos).refspy0/3.,(*(*info).plottitles).refspxtitle,ALIGNMENT=0.5,$
        COLOR = (*(*info).plotparams).plotcol,/NORMAL
    ; Loop over all diagnostics for plotting of detailed spectrum
    FOR d=0,(*(*info).intparams).ndisp_refdiagnostics-1 DO BEGIN
      disp_idx = (WHERE((*(*info).intparams).disp_refdiagnostics EQ 1))[d]
      ; Determine xrange to display
      IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN $
        xrange = (*(*info).dataparams).reflps[[(*(*info).intparams).refdiag_start[disp_idx],$
          ((*(*info).intparams).refdiag_start[disp_idx]+(*(*info).intparams).refdiag_width[disp_idx]-1)]]
      IF (d EQ 0) THEN offset = 0 ELSE offset = TOTAL(diag_range[0:(d-1)])
      ; Determine lower left corner position of plot
      refspx0 = (*(*info).plotpos).refspx0 + offset
      refspx1 = diag_range[d] + refspx0
      PLOT, (*(*info).dataparams).lps, *(*(*info).dispparams).tarr_ref, $
        /NODATA, YR=[t_low_y,t_upp_y], /YS, $
  			XRANGE=xrange, XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, $
        POS = [refspx0,(*(*info).plotpos).refspy0,refspx1,(*(*info).plotpos).refspy1], $
        YTICKLEN=(*(*info).plotaxes).refspyticklen, XTICKLEN=(*(*info).plotaxes).refspxticklen, $
        YTICKNAME=REPLICATE(' ',60), $;XTICKNAME=xtickname, $
        XTICKINTERVAL=(*(*info).plotaxes).xreftickinterval, $
        BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol,$
        /NOERASE
      ; Display Doppler top axis if Doppler set, else regular top axis
  		IF ((*(*info).plotswitch).v_dop_set_ref EQ 1) THEN BEGIN
        IF (d EQ 0) THEN $
          XYOUTS,((*(*info).plotpos).refspx1-(*(*info).plotpos).refspx0)/2.+$
            (*(*info).plotpos).refspx0,$
            (*(*info).plotpos).refspy0/5.*3+(*(*info).plotpos).refspy1, topxtitle, $
            ALIGNMENT=0.5, COLOR = (*(*info).plotparams).plotcol,/NORMAL
        topxtitle = ''
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).refspxticklen, $
          XRANGE = [(*(*(*info).plotaxes).v_dop_ref[disp_idx])[0], $
          (*(*(*info).plotaxes).v_dop_ref[disp_idx])[diag_widths[d]-1]], XSTYLE=1, $
          COLOR = (*(*info).plotparams).plotcol, XTITLE=topxtitle, $
          XTICKINT=(*(*info).plotaxes).xrefdoptickinterval 
      ENDIF ELSE $
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).refspxticklen, XRANGE = xrange, XSTYLE=1, $
          XTITLE = topxtitle, COLOR = (*(*info).plotparams).plotcol, $
          XTICKNAME=xtickname
    ENDFOR
  ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refspwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_REFSP_RESIZE, event
; Reference temporal spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_DRAW_SPECTRAL_REF, event, /SP_ONLY
END

PRO CRISPEX_DISPLAYS_REFSP_TOGGLE, event, NO_DRAW=no_draw
; Reference temporal spectrum display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
  ; If window doesn't exist, create it
	IF ((*(*info).winswitch).showrefsp EQ 0) THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+$
      ((*(*info).plottitles).refspwintitle)[(*(*info).plotswitch).refheightset]
    ; Create window
		CRISPEX_WINDOW, (*(*info).winsizes).refspxres, (*(*info).winsizes).refspyres, $
      (*(*info).winids).root, title, refsptlb, refspwid, (*(*info).winsizes).spxoffset, $
      (*(*info).winsizes).spyoffset+((*(*info).winswitch).showsp * (*(*info).winsizes).ydelta),$
      DRAWID = refspdrawid, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_REFSP_RESIZE'
    ; Save window variables
		(*(*info).winids).refsptlb = refsptlb		&	(*(*info).winids).refspwid = refspwid	&	(*(*info).winswitch).showrefsp = 1
		(*(*info).winids).refspdrawid = refspdrawid	&	(*(*info).winids).refspwintitle = title 
		WIDGET_CONTROL, (*(*info).winids).refsptlb, SET_UVALUE = info
    ; Fill window
    CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
    CRISPEX_UPDATE_REFSPSLICE, event
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_SPECTRAL_REF, event, /SP_ONLY
  ; If window exists, destroy it
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).refsptlb, /DESTROY
		(*(*info).winids).refsptlb = 0
		(*(*info).winswitch).showrefsp = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refsptlb,(*(*info).winids).refspwid,$
      (*(*info).winids).refspdrawid], labels=['refsptlb','refspwid','refspdrawid']
END

PRO CRISPEX_DISPLAYS_SJI_TOGGLE, event, NO_DRAW=no_draw
; Reference image window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	(*(*info).winswitch).showsji = event.SELECT
	IF (*(*info).winswitch).showsji THEN BEGIN
		title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Slit-jaw image'
		CRISPEX_WINDOW, (*(*info).winsizes).sjiwinx, (*(*info).winsizes).sjiwiny, $
      (*(*info).winids).root, title, sjitlb, sjiwid, $
      (*(*info).winsizes).lsxoffset+(*(*info).winsizes).lswinx+(*(*info).winsizes).xdelta, 0, $
      DRAWID = sjidrawid, DRAWBASE = sjidrawbase
		(*(*info).winids).sjitlb = sjitlb		&	(*(*info).winids).sjiwid = sjiwid	
    (*(*info).winids).sjidrawid = sjidrawid
		(*(*info).winids).sjidrawbase = sjidrawbase	&	(*(*info).winids).sjiwintitle = title
		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
			CRISPEX_UPDATE_T, event
			CRISPEX_DRAW_SJI, event
		ENDIF
;		WIDGET_CONTROL, dopdrawid, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, /TRACKING_EVENTS,/DRAW_BUTTON_EVENTS
;		WIDGET_CONTROL, doptlb, SET_UVALUE = info
;		XMANAGER, 'CRISPEX', doptlb, /NO_BLOCK
	ENDIF ELSE BEGIN
;		(*(*info).dispswitch).drawsji = 0
;		IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
;			CRISPEX_DRAW_SPECTRAL, event
;			CRISPEX_DRAW_TIMESLICES, event
;		ENDIF
		WIDGET_CONTROL, (*(*info).winids).sjitlb, /DESTROY
		(*(*info).winids).sjitlb = 0
	ENDELSE
;	IF (*(*info).overlayswitch).mask THEN CRISPEX_MASK_BUTTONS_SET, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).sjitlb,(*(*info).winids).sjiwid,$
      (*(*info).winids).sjidrawid], labels=['sjitlb','sjiwid','sjidrawid']
END

PRO CRISPEX_DISPLAYS_SP_REPLOT_AXES, event, NO_AXES=no_axes
; Updates temporal spectrum display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).spwid
	IF (*(*info).plotswitch).v_dop_set THEN extratitle = '!C' ELSE extratitle = ''
	IF (*(*info).plotswitch).multichannel THEN $
    title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]+extratitle $
  ELSE $
    title = ''
  ; Set axes parameters depending on # of diagnostics
  IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
    xticklen = 1E-9   & xtitle = ''
    xtickname = REPLICATE(' ',60)
    IF KEYWORD_SET(NO_AXES) THEN BEGIN
      yticklen = 1E-9 & ytitle = ''
      ytickname = REPLICATE(' ',60)
    ENDIF ELSE BEGIN
      yticklen = (*(*info).plotaxes).spyticklen
      ytitle = (*(*info).plottitles).spytitle
      ytickname = ''
    ENDELSE
  ENDIF ELSE BEGIN
    xticklen = (*(*info).plotaxes).spxticklen
    xtitle = (*(*info).plottitles).spxtitle
    yticklen = (*(*info).plotaxes).spyticklen
    ytitle = (*(*info).plottitles).spytitle
    xtickname = ''  & ytickname = ''
  ENDELSE
  t_low_y = (*(*info).dispparams).t_low_main
  t_upp_y = (*(*info).dispparams).t_upp_main
  xrange = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], $
            (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]]
  correct_axes = (FLOOR(ALOG10(ABS(t_upp_y))) GE 3)
  IF ((t_low_y NE 0.) AND ~KEYWORD_SET(correct_axes)) THEN $
    correct_axes = (FLOOR(ALOG10(ABS(t_low_y))) LE -2) 
  IF correct_axes THEN BEGIN
 	  order_corr = FLOOR(ALOG10(ABS(t_upp_y)))
 	  IF ~KEYWORD_SET(NO_AXES) THEN ytitle += ' (x10!U'+STRTRIM(order_corr,2)+'!N)'
    t_low_y /= (10.^(order_corr))
    t_upp_y /= (10.^(order_corr))
  ENDIF 
  topxtitle = title
  IF (*(*info).plotswitch).v_dop_set THEN topxtitle += 'Doppler velocity [km/s]'
  ; Plot basic axes box
  PLOT, (*(*info).dataparams).lps, *(*(*info).dispparams).tarr_main, $
    YR = [t_low_y,t_upp_y], /YS, XR=xrange, $
    XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, $
  	YTICKLEN = yticklen, YTITLE = ytitle, YTICKNAME = ytickname, $
    XTICKLEN = xticklen, XTITLE = xtitle, XTICKNAME = xtickname, $
    POS = [(*(*info).plotpos).spx0,(*(*info).plotpos).spy0,$
           (*(*info).plotpos).spx1,(*(*info).plotpos).spy1], $
    BACKGROUND=(*(*info).plotparams).bgplotcol, COLOR=(*(*info).plotparams).plotcol, $
    /NODATA, NOERASE=KEYWORD_SET(NO_AXES)
  ; In case of multiple diagnostics and regular replotting of axes
  IF ~KEYWORD_SET(NO_AXES) THEN BEGIN
    ; Determine proportional spectral window sizes
    diag_widths = (*(*info).intparams).diag_width[$
      WHERE((*(*info).intparams).disp_diagnostics EQ 1)]
    diag_ratio = diag_widths / FLOAT(TOTAL(diag_widths))
    diag_range = diag_ratio * (*(*info).plotpos).xplspw 
    whererangemin = WHERE(diag_widths EQ MIN(diag_widths, /NAN))
    ; Plot xtitle(s)
    IF ((*(*info).intparams).ndiagnostics GT 1) THEN $
      XYOUTS,(*(*info).plotpos).xplspw/2.+(*(*info).plotpos).spx0,$
        (*(*info).plotpos).spy0/3.,(*(*info).plottitles).spxtitle,ALIGNMENT=0.5, $
        COLOR = (*(*info).plotparams).plotcol,/NORMAL
    ; Loop over all diagnostics for plotting of detailed spectrum
    FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
      disp_idx = (WHERE((*(*info).intparams).disp_diagnostics EQ 1))[d]
      ; Determine xrange to display
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        xrange = (*(*info).dataparams).lps[[(*(*info).intparams).diag_start[disp_idx],$
          ((*(*info).intparams).diag_start[disp_idx]+(*(*info).intparams).diag_width[disp_idx]-1)]]
        vdop_xrange = (*(*(*info).plotaxes).v_dop[disp_idx])[[0,diag_widths[d]-1]]
      ENDIF ELSE $
        vdop_xrange = (*(*(*info).plotaxes).v_dop[0])[$
          [(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp]]
      IF (d EQ 0) THEN offset = 0 ELSE offset = TOTAL(diag_range[0:(d-1)])
      ; Determine lower left corner position of plot
      spx0 = (*(*info).plotpos).spx0 + offset
      spx1 = diag_range[d] + spx0
  		PLOT, (*(*info).dataparams).lps, *(*(*info).dispparams).tarr_main, $
        /NODATA, YR=[t_low_y,t_upp_y], $
  			/YS, XRANGE=xrange, XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, $
        POS = [spx0,(*(*info).plotpos).spy0,spx1,(*(*info).plotpos).spy1], $
        YTICKLEN=(*(*info).plotaxes).spyticklen, XTICKLEN=(*(*info).plotaxes).spxticklen, $
        YTICKNAME=REPLICATE(' ',60), $;XTICKNAME=xtickname, $
        XTICKINTERVAL=(*(*info).plotaxes).xtickinterval, $;,XMINOR=, $
        BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol,$
        /NOERASE
      ; Display Doppler top axis if Doppler set, else regular top axis
  		IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN BEGIN
        IF (d EQ 0) THEN $
          XYOUTS,((*(*info).plotpos).spx1-(*(*info).plotpos).spx0)/2.+$
            (*(*info).plotpos).spx0,$
            (*(*info).plotpos).spy0/5.*3+(*(*info).plotpos).spy1, topxtitle, $
            ALIGNMENT=0.5, COLOR = (*(*info).plotparams).plotcol,/NORMAL
        topxtitle = ''
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, $
          XRANGE = vdop_xrange, XSTYLE=1, $
          COLOR=(*(*info).plotparams).plotcol, XTITLE=topxtitle, $
          XTICKINTERVAL=(*(*info).plotaxes).xdoptickinterval
      ENDIF ELSE $
  			AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).spxticklen, XRANGE = xrange, XSTYLE=1, $
          XTITLE = topxtitle, COLOR = (*(*info).plotparams).plotcol, $
          XTICKNAME=xtickname
    ENDFOR
  ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).spwid], labels=['Window ID for replot']
END

PRO CRISPEX_DISPLAYS_SP_RESIZE, event
; Temporal spectrum window resize handler, gets new window dimensions and calls (re)display routines
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_DRAW_GET_SPECTRAL_AXES, event, /MAIN
	CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
  CRISPEX_DRAW_SPECTRAL_MAIN, event,/SP_ONLY
END

PRO CRISPEX_DISPLAYS_SP_TOGGLE, event, NO_DRAW=no_draw
; Temporal spectrum display window creation procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
  ; If window doesn't exist, create it
	IF ((*(*info).winswitch).showsp EQ 0) THEN BEGIN
		wintitle = 'CRISPEX'+(*(*info).sesparams).instance_label+': '+$
      ((*(*info).plottitles).spwintitle)[(*(*info).plotswitch).heightset]
    ; Create window
		CRISPEX_WINDOW, (*(*info).winsizes).spxres, (*(*info).winsizes).spyres, $
      (*(*info).winids).root, wintitle, tlb, wid, (*(*info).winsizes).spxoffset, $
      (*(*info).winsizes).spyoffset, DRAWID = spdrawid, RESIZING = 1, $
      RES_HANDLER = 'CRISPEX_DISPLAYS_SP_RESIZE'
    ; Save window variables
		(*(*info).winids).sptlb = tlb		&	(*(*info).winids).spwid = wid	&	(*(*info).winswitch).showsp = 1
		(*(*info).winids).spdrawid = spdrawid	&	(*(*info).winids).spwintitle = wintitle
		WIDGET_CONTROL, (*(*info).winids).sptlb, SET_UVALUE = info
    ; Fill window
    CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
    CRISPEX_UPDATE_SPSLICE, event
		IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW_SPECTRAL_MAIN, event,/SP_ONLY
  ; If window does exist, destroy it
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, (*(*info).winids).sptlb, /DESTROY
		(*(*info).winids).sptlb = 0
		(*(*info).winswitch).showsp = 0
	ENDELSE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).sptlb,(*(*info).winids).spwid,$
      (*(*info).winids).spdrawid], labels=['sptlb','spwid','spdrawid']
END

PRO CRISPEX_DISPLAYS_STOKES_SELECT_XY_RECOVER_YRANGE, event
; Restores lower/upper y-values of specific Stokes component detailed spectrum plot for input
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_y_text, $
    SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],2)
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_y_text, $
    SET_VALUE = STRTRIM((*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],2)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],$
                                 (*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s]], $
                            		labels=['Lower detspect y-value','Upper detspect y-value']
	IF (*(*info).winswitch).showint THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlsint).lower_y_int_text, $
      SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],2)
		WIDGET_CONTROL, (*(*info).ctrlsint).upper_y_int_text, $
      SET_VALUE = STRTRIM((*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s],2)
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s],$
                                   (*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]], $
                            			labels=['Lower int y-value','Upper int y-value']
	ENDIF
END

;================================================================================= DISPLAY RANGE PROCEDURES
PRO CRISPEX_DISPRANGE_INT_LOW, event
; Handles change in lower y-value of intensity versus time display window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).intparams).lock_t = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).intparams).lock_t], labels=['Lock main to int temporal range']
	CRISPEX_DISPRANGE_INT_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_INT_T_RANGE, event
; Locks main temporal range to intensity versus time display temporal range
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

PRO CRISPEX_DISPRANGE_LS_RANGE, event, NO_DRAW=no_draw
; Determines range from change in lower or upper y-value of detailed spectrum display window and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).plotaxes).ls_yrange_ref = (*(*info).plotaxes).ls_upp_y_ref - (*(*info).plotaxes).ls_low_y_ref ELSE $
		(*(*(*info).plotaxes).ls_yrange)[(*(*info).dataparams).s] = (*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s] - (*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s] 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*(*info).plotaxes).ls_low_y)[(*(*info).dataparams).s],(*(*(*info).plotaxes).ls_upp_y)[(*(*info).dataparams).s],$
		(*(*info).plotaxes).ls_low_y_ref,(*(*info).plotaxes).ls_upp_y_ref], labels=['Lower main y-value','Upper main y-value','Lower ref y-value','Upper ref y-value']
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_DISPRANGE_LS_SCALE_SELECT, event
; Handles the selection of scaling (or not) of the detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).dispswitch).ref_detspect_scale = event.SELECT ELSE (*(*info).dispswitch).detspect_scale = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dispswitch).detspect_scale,(*(*info).dispswitch).ref_detspect_scale], labels=['Scale main detspect','Scale ref detspect']
	CRISPEX_DISPRANGE_LS_SCALE, event
END


PRO CRISPEX_DISPRANGE_LS_SCALE, event
; Handles the scaling (or not) of the detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).ctrlsswitch).imrefdetspect THEN CRISPEX_DISPRANGE_LS_SCALE_REF, event ELSE CRISPEX_DISPRANGE_LS_SCALE_MAIN, event
	CRISPEX_DISPRANGE_LS_RANGE, event
END

PRO CRISPEX_DISPRANGE_LS_SCALE_MAIN, event
; Handles the scaling (or not) of the main detailed spectrum to the maximum of the average spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE=info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).ctrlsswitch).imrefdetspect THEN (*(*info).plotswitch).ref_subtract = event.SELECT ELSE (*(*info).plotswitch).subtract = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).plotswitch).subtract,(*(*info).plotswitch).ref_subtract], labels=['Subtract from main','Subtract from ref']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DISPRANGE_T_LOW, event
; Handles change in lower t-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, GET_VALUE = textvalue
	(*(*info).dispparams).t_low = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).t_low GE (*(*info).dispparams).t_upp) THEN (*(*info).dispparams).t_low = (*(*info).dispparams).t_upp - 1
	IF ((*(*info).dispparams).t_low LT (*(*info).dispparams).t_first) THEN (*(*info).dispparams).t_low = (*(*info).dispparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	IF ((*(*info).winswitch).showint AND (*(*info).intparams).lock_t) THEN BEGIN
		(*(*info).plotaxes).int_low_t = (*(*info).dispparams).t_low
		CRISPEX_DISPRANGE_INT_T_RANGE, event
		IF (*(*info).winswitch).showint THEN WIDGET_CONTROL, (*(*info).ctrlsint).lower_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_low,2)
	ENDIF ELSE CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_T_UPP, event
; Handles change in upper t-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, GET_VALUE = textvalue
	(*(*info).dispparams).t_upp = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).t_upp LE (*(*info).dispparams).t_low) THEN (*(*info).dispparams).t_upp = (*(*info).dispparams).t_low + 1
	IF ((*(*info).dispparams).t_upp GT (*(*info).dispparams).t_last) THEN (*(*info).dispparams).t_upp = (*(*info).dispparams).t_last
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	IF ((*(*info).winswitch).showint AND (*(*info).intparams).lock_t) THEN BEGIN
		(*(*info).plotaxes).int_upp_t = (*(*info).dispparams).t_upp
		CRISPEX_DISPRANGE_INT_T_RANGE, event
		IF (*(*info).winswitch).showint THEN WIDGET_CONTROL, (*(*info).ctrlsint).upper_t_int_text, SET_VALUE = STRTRIM((*(*info).dispparams).t_upp,2)
	ENDIF ELSE CRISPEX_DISPRANGE_T_RANGE, event
END

PRO CRISPEX_DISPRANGE_T_RANGE, event, NO_DRAW=no_draw, T_SET=t_set, RESET=reset
; Determines range from change in lower or upper t-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).winswitch).showretrdet THEN BEGIN
		(*(*info).dispparams).t_low = (*(*(*info).detparams).t_restored)[(*(*info).detparams).idx] - $
      (*(*info).detparams).delta_t_dn > (*(*info).dispparams).t_first
		(*(*info).dispparams).t_upp = (*(*(*info).detparams).t_restored)[(*(*info).detparams).idx] + $
      (*(*info).detparams).delta_t_up < (*(*info).dispparams).t_last
	ENDIF
	(*(*info).dispparams).t_range = (*(*info).dispparams).t_upp - (*(*info).dispparams).t_low + 1
	IF ((*(*info).dispparams).t_range NE (*(*info).dataparams).nt) THEN $
    WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, /SENSITIVE $
  ELSE $
    WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
  IF KEYWORD_SET(T_SET) THEN $
    (*(*info).dispparams).t = t_set $
  ELSE IF (~KEYWORD_SET(RESET) AND ((*(*info).winswitch).showretrdet EQ 0)) THEN BEGIN
    IF ((*(*info).dispparams).t LT (*(*info).dispparams).t_low) THEN $
      (*(*info).dispparams).t = (*(*info).dispparams).t_low $
    ELSE IF ((*(*info).dispparams).t GT (*(*info).dispparams).t_upp) THEN $
      (*(*info).dispparams).t = (*(*info).dispparams).t_upp
  ENDIF
	WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_SLIDER_MIN=(*(*info).dispparams).t_low, $
    SET_SLIDER_MAX=(*(*info).dispparams).t_upp, SET_VALUE=(*(*info).dispparams).t
	IF ((*(*info).dispparams).t_range - 1 EQ 1) THEN BEGIN
		t_step = 1
		t_sens = 0 
	ENDIF ELSE BEGIN
		t_step = (*(*info).pbparams).t_step
		t_sens = 1
	ENDELSE
  ; Set real lower/upper time values for main
  (*(*info).dispparams).t_low_main = (*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t_low]
  tarr_main_sel = (*(*(*info).dispparams).tarr_main)[$
    (*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
  (*(*info).dispparams).t_upp_main = tarr_main_sel[(WHERE(tarr_main_sel NE 0))[-1]] ;Fix tarr[-1]=0
  IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN
    ; Set real lower/upper time values for reference
    (*(*info).dispparams).t_low_ref = (*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t_low]
    tarr_ref = (*(*(*info).dispparams).tarr_ref)[$
      (*(*info).dispparams).t_low:(*(*info).dispparams).t_upp]
    (*(*info).dispparams).t_upp_ref = tarr_ref[(WHERE(tarr_ref NE 0))[-1]] ;Fix tarr[-1]=0
  ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp], $
      labels=['Lower t-value','Upper t-value']
	WIDGET_CONTROL, (*(*info).ctrlscp).t_step_slider, SET_SLIDER_MAX=(*(*info).dispparams).t_range-1,$
    SET_VALUE=t_step, SENSITIVE=t_sens
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		CRISPEX_UPDATE_T, event
		IF (*(*info).winswitch).showsp THEN BEGIN
      CRISPEX_UPDATE_SPSLICE, event
      CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
    ENDIF
		IF (*(*info).winswitch).showrefsp THEN BEGIN
      CRISPEX_UPDATE_REFSPSLICE, event
      CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
    ENDIF
    IF (*(*info).winswitch).showphis THEN CRISPEX_UPDATE_SLICES, event, /NO_DRAW
		IF (*(*info).winswitch).showloop THEN CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event
		IF (*(*info).winswitch).showrestloop THEN CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event
		IF (*(*info).winswitch).showretrdet THEN CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
		CRISPEX_DRAW, event
	ENDIF
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF (*(*info).dataswitch).onecube THEN $
      WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SET_VALUE='Update spectral windows', $
        SENSITIVE=1 ;$
;    ELSE $
;      WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE=1
	ENDIF
END

PRO CRISPEX_DISPRANGE_T_RESET, event, NO_DRAW=no_draw, T_SET=t_set
; Handles reset of temporal boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dispparams).t_upp = (*(*info).dispparams).t_last
	(*(*info).dispparams).t_low = (*(*info).dispparams).t_first
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_t_text, $
    SET_VALUE=STRTRIM((*(*info).dispparams).t_upp,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_t_text, $
    SET_VALUE=STRTRIM((*(*info).dispparams).t_low,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
	IF ((*(*info).winswitch).showint AND (*(*info).intparams).lock_t) THEN $
    CRISPEX_DISPRANGE_INT_T_RESET, event $
  ELSE $
    CRISPEX_DISPRANGE_T_RANGE, event, NO_DRAW=no_draw, T_SET=t_set, /RESET
  IF (N_ELEMENTS(T_SET) EQ 1) THEN print,t_set
END

PRO CRISPEX_DISPRANGE_GET_WARP, event, PHIS=phis
; Handles determining warping tie points
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  lps_sel = ((*(*info).dataparams).lps)[(*(*info).dispparams).lp_low:(*(*info).dispparams).lp_upp]
  min_lps = MIN(lps_sel, /NAN)
  IF KEYWORD_SET(PHIS) THEN $
    nt = (*(*info).phiparams).nphi $
  ELSE $
    nt = (*(*info).dataparams).nt
    xi = FINDGEN((*(*info).dispparams).lp_range) # REPLICATE(1,nt)
    xo = ((lps_sel-min_lps) / FLOAT(MAX(lps_sel-min_lps, /NAN)) * $
          (*(*info).dispparams).lp_range) # REPLICATE(1,nt)
    yi = REPLICATE(1,(*(*info).dispparams).lp_range) # FINDGEN(nt)
  IF KEYWORD_SET(PHIS) THEN BEGIN
    *(*(*info).dispparams).phisxi = xi
    *(*(*info).dispparams).phisyi = yi
    *(*(*info).dispparams).phisxo = xo
    *(*(*info).dispparams).phisyo = *(*(*info).dispparams).phisyi
  ENDIF ELSE BEGIN
    (*(*info).dispparams).xi = xi
    (*(*info).dispparams).yi = yi
    (*(*info).dispparams).xo = xo
    (*(*info).dispparams).yo = (*(*info).dispparams).yi
  ENDELSE
END

PRO CRISPEX_DISPRANGE_GET_WARP_TRIANGULATE, event, xo, yo, xi, yi, slice, PHIS=phis, NX=nx, NY=ny
  WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  WIDGET_CONTROL, /HOURGLASS
  ; Functionality taken from IDL's WARP_TRI
  IF (N_ELEMENTS(NX) NE 1) THEN nx = (SIZE(slice))[1]
  IF (N_ELEMENTS(NY) NE 1) THEN ny = (SIZE(slice))[2]
  gs = [1,1]				        ;Grid spacing
  b = [0,0, nx-1, ny-1]			;Bounds
  ; Triangulate given the input and output tie points
  TRIANGULATE, xo, yo, tr, bounds
  xtri_new = TRIGRID(xo, yo, xi, tr, gs, b)
  ytri_new = TRIGRID(xo, yo, yi, tr, gs, b)
  IF KEYWORD_SET(PHIS) THEN BEGIN
    *(*(*info).dispparams).phisxtri = xtri_new
    *(*(*info).dispparams).phisytri = ytri_new
  ENDIF ELSE BEGIN
    *(*(*info).dispparams).xtri = xtri_new
    *(*(*info).dispparams).ytri = ytri_new
  ENDELSE
END

PRO CRISPEX_DISPRANGE_LP_LOW, event, LP_SET=lp_set, NO_DRAW=no_draw
; Handles change in lower lp-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(LP_SET) EQ 1) THEN $
    textvalue = lp_set $
  ELSE $
	  WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, GET_VALUE = textvalue
	(*(*info).dispparams).lp_low = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).lp_low GE (*(*info).dispparams).lp_upp) THEN $
    (*(*info).dispparams).lp_low = (*(*info).dispparams).lp_upp - 1
	IF ((*(*info).dispparams).lp_low LT (*(*info).dispparams).lp_first) THEN $
    (*(*info).dispparams).lp_low = (*(*info).dispparams).lp_first
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2)
	CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
END

PRO CRISPEX_DISPRANGE_LP_UPP, event, LP_SET=lp_set, NO_DRAW=no_draw
; Handles change in upper lp-value of accessed data cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(LP_SET) EQ 1) THEN $
    textvalue = lp_set $
  ELSE $
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, GET_VALUE = textvalue
	(*(*info).dispparams).lp_upp = FLOAT(textvalue[0])
	IF ((*(*info).dispparams).lp_upp LE (*(*info).dispparams).lp_low) THEN $
    (*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_low + 1
	IF ((*(*info).dispparams).lp_upp GT (*(*info).dispparams).lp_last) THEN $
    (*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_last
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2)
	CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
END

PRO CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
; Determines range from change in lower or upper s-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dispparams).lp_range = (*(*info).dispparams).lp_upp - (*(*info).dispparams).lp_low + 1
	IF (((*(*info).intparams).ndiagnostics EQ 1) AND $
    ((*(*info).dispparams).lp_range NE (*(*info).dataparams).nlp)) THEN $
    WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, /SENSITIVE $
  ELSE $
    WIDGET_CONTROL, (*(*info).ctrlscp).reset_lprange_but, SENSITIVE = 0
	IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN $
    (*(*info).dataparams).lp = (*(*info).dispparams).lp_low $
  ELSE IF ((*(*info).dataparams).lp GT (*(*info).dispparams).lp_upp) THEN $
    (*(*info).dataparams).lp = (*(*info).dispparams).lp_upp $
  ELSE $
		(*(*info).dataparams).lp = (*(*info).dataparams).lp
  CRISPEX_SLIDER_LP_UPDATE, event, /NO_DRAW
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_SLIDER_MIN = (*(*info).dispparams).lp_low, $
                  SET_SLIDER_MAX = (*(*info).dispparams).lp_upp, SET_VALUE = (*(*info).dataparams).lp
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_slider, SET_SLIDER_MIN=(*(*info).dispparams).lp_low, $
    SET_SLIDER_MAX=(*(*info).dispparams).lp_upp, SET_VALUE=(*(*info).dataparams).lp, $
    SENSITIVE=(((*(*info).dispparams).lp_range-1) NE 1)
  IF ((*(*info).dispswitch).warpspslice AND (*(*info).winswitch).showphis) THEN BEGIN
    WIDGET_CONTROL, /HOURGLASS
    CRISPEX_DISPRANGE_GET_WARP, event, /PHIS
    CRISPEX_DISPRANGE_GET_WARP_TRIANGULATE, event, $
      *(*(*info).dispparams).phisxo,*(*(*info).dispparams).phisyo,$
      *(*(*info).dispparams).phisxi,*(*(*info).dispparams).phisyi, $
      *(*(*info).data).phislice, /PHIS, NX=(*(*info).dispparams).lp_range, NY=(*(*info).phiparams).nphi
  ENDIF
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		CRISPEX_UPDATE_T, event
		IF (*(*info).winswitch).showsp THEN BEGIN
      CRISPEX_UPDATE_SPSLICE, event
      CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
    ENDIF
		IF (*(*info).winswitch).showphis THEN BEGIN
      CRISPEX_UPDATE_PHISLICE, event, /NO_DRAW
      CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event
    ENDIF
		CRISPEX_DRAW, event
	ENDIF
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF (*(*info).dataswitch).onecube THEN $
      WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SET_VALUE = 'Update spectral windows', $
                      SENSITIVE = 1  $
    ELSE $
      WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp], $
                        labels=['Lower lp-value','Upper lp-value']
END

PRO CRISPEX_DISPRANGE_LP_RESET, event, NO_DRAW=no_draw
; Handles reset of spectral boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dispparams).lp_upp = (*(*info).dispparams).lp_last
	(*(*info).dispparams).lp_low = (*(*info).dispparams).lp_first
	WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, $
                  SET_VALUE = STRTRIM((*(*info).dispparams).lp_upp,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, $
                  SET_VALUE = STRTRIM((*(*info).dispparams).lp_low,2), /SENSITIVE
	WIDGET_CONTROL, (*(*info).ctrlscp).reset_trange_but, SENSITIVE = 0
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, /SENSITIVE
	CRISPEX_DISPRANGE_LP_RANGE, event, NO_DRAW=no_draw
END

PRO CRISPEX_DISPRANGE_LP_REF_RANGE, event, NO_DRAW=no_draw, LP_LOW_SET=lp_low_set, $
  LP_UPP_SET=lp_upp_set
; Determines range from change in lower or upper s-value and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(LP_LOW_SET) EQ 1) THEN (*(*info).dispparams).lp_ref_low = lp_low_set
  IF (N_ELEMENTS(LP_UPP_SET) EQ 1) THEN (*(*info).dispparams).lp_ref_upp = lp_upp_set
	(*(*info).dispparams).lp_ref_range = (*(*info).dispparams).lp_ref_upp - (*(*info).dispparams).lp_ref_low + 1
	IF ((*(*info).dataparams).lp_ref LT (*(*info).dispparams).lp_ref_low) THEN $
    (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_low $
  ELSE IF ((*(*info).dataparams).lp_ref GT (*(*info).dispparams).lp_ref_upp) THEN $
    (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_upp 
  CRISPEX_SLIDER_LP_UPDATE, event, /NO_DRAW
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, $
    SET_SLIDER_MIN = (*(*info).dispparams).lp_ref_low, $
    SET_SLIDER_MAX = (*(*info).dispparams).lp_ref_upp, SET_VALUE = (*(*info).dataparams).lp_ref
	IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
		IF (*(*info).winswitch).showrefsp THEN BEGIN
      CRISPEX_UPDATE_REFSPSLICE, event
      CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
    ENDIF
		CRISPEX_UPDATE_T, event
		CRISPEX_DRAW, event
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).lp_ref_low,$
    (*(*info).dispparams).lp_ref_upp], labels=['Lower lp-value','Upper lp-value']
END

PRO CRISPEX_DISPRANGE_LP_REF_RESET, event, NO_DRAW=no_draw
; Handles reset of reference spectral boundaries and calls (re)display routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dispparams).lp_ref_upp = (*(*info).dispparams).lp_ref_last
	(*(*info).dispparams).lp_ref_low = (*(*info).dispparams).lp_ref_first
	CRISPEX_DISPRANGE_LP_REF_RANGE, event, NO_DRAW=no_draw
END
;================================================================================= DISPLAY DRAW PROCEDURES
PRO CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, no_cursor=no_cursor, no_number=no_number,$
      thick=thick, no_endpoints=no_endpoints, symsize=symsize, draw_mask=draw_mask, SJI=sji
  ; Handles overplotting of the cursor, slits and loop paths
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).phiparams).d_nphi_set = (*(*info).zooming).factor * (*(*info).phiparams).nphi_set
  IF KEYWORD_SET(SJI) THEN BEGIN
    sx_loc = (*(*info).curs).sxsji
    sy_loc = (*(*info).curs).sysji
    IF (*(*info).overlayswitch).sjiraster THEN CRISPEX_DRAW_RASTER_OVERLAYS, event
  ENDIF ELSE BEGIN
    sx_loc = (*(*info).curs).sx
    sy_loc = (*(*info).curs).sy
  ENDELSE
	IF (*(*info).winswitch).showphis THEN BEGIN
		IF ((*(*info).zooming).factor EQ 1) THEN $
			PLOTS, ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywinx/$
              (2.*(*(*info).dataparams).nx))*COS((*(*info).phiparams).angle*!DTOR) + $
              sx_loc, $
				      ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywiny/$
              (2.*(*(*info).dataparams).ny))*SIN((*(*info).phiparams).angle*!DTOR) + $
              sy_loc, /DEVICE, COLOR = !P.COLOR $
		ELSE $
      PLOTS, ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywinx/$
              (2.*(*(*info).dataparams).nx))*COS((*(*info).phiparams).angle*!DTOR) + $
              sx_loc, $
			        ([-1,1]*(*(*info).phiparams).d_nphi_set*(*(*info).winsizes).xywiny/$
              (2.*(*(*info).dataparams).ny))*SIN((*(*info).phiparams).angle*!DTOR) + $
              sy_loc, /DEVICE, COLOR = !P.COLOR
	ENDIF
	IF (*(*info).overlayswitch).loopslit AND ((*(*info).loopparams).np GT 0) THEN BEGIN
		CRISPEX_ZOOM_LOOP, event
		IF ~KEYWORD_SET(NO_ENDPOINTS) THEN $
      PLOTS, *(*(*info).overlayparams).sxp, *(*(*info).overlayparams).syp, /DEVICE, $
              COLOR=!P.COLOR, PSYM=1, THICK=thick, SYMSIZE=symsize
		IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN $
      PLOTS,*(*(*info).overlayparams).sxr,*(*(*info).overlayparams).syr, /DEVICE, $
            COLOR=!P.COLOR, PSYM=3, THICK=thick $
		ELSE $
      PLOTS,*(*(*info).overlayparams).sxr,*(*(*info).overlayparams).syr,/DEVICE, $
            COLOR=!P.COLOR, LINESTYLE=(*(*info).overlayparams).loop_linestyle, THICK=thick
	ENDIF ELSE IF ((*(*info).meas).np GE 1) THEN BEGIN
		CRISPEX_ZOOM_MEAS, event
		IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
			PLOTS, *(*(*info).meas).sxp,*(*(*info).meas).syp, /DEVICE, COLOR=!P.COLOR, $
        PSYM=1, THICK=thick, SYMSIZE=symsize
			PLOTS, *(*(*info).meas).sxp,*(*(*info).meas).syp, /DEVICE, COLOR=!P.COLOR, $
        LINESTYLE=0, THICK=thick, SYMSIZE=symsize
		ENDIF
	ENDIF ELSE IF ~KEYWORD_SET(NO_CURSOR) THEN $
    PLOTS, sx_loc,sy_loc, /DEVICE, PSYM=1, COLOR=curscolor, $
      THICK=thick, SYMSIZE=symsize
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).curs).sx,(*(*info).curs).sy,(*(*info).winsizes).xywinx,$
                                  (*(*info).winsizes).xywiny], labels=['sx','sy','xywinx','xywiny']
	CRISPEX_DRAW_LOOP_OVERLAYS, event, NO_NUMBER=no_number, THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize
	IF draw_mask THEN CRISPEX_DRAW_MASK_OVERLAYS, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_0, event
; Handles setting of loop path linestyle to solid
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayparams).loop_linestyle = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_1, event
; Handles setting of loop path linestyle to dotted (i.e. actual loop points, not dots between the vertices)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayparams).loop_linestyle = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_LINESTYLE_2, event
; Handles setting of loop path linestyle to dashed
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayparams).loop_linestyle = 2
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayparams).loop_linestyle], labels=['Loop linestyle']
	CRISPEX_DRAW, event
END

PRO CRISPEX_DRAW_LOOP_OVERLAYS, event, no_number=no_number, thick=thick, $
  no_endpoints=no_endpoints, symsize=symsize
; Handles overplotting of loop paths from the restored and retrieved loops as well as from the retrieved detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (((*(*info).loopswitch).restore_loops EQ 1) AND $
      ((*(*info).restoreparams).cfilecount GT 0)) THEN BEGIN
		drawrestore = 0
		FOR i=0,(*(*info).restoreparams).cfilecount-1 DO BEGIN
			IF (((*(*info).overlayswitch).overlalways EQ 1) OR $
          ((*(*(*info).restoreparams).lp_restored)[i] EQ (*(*info).dataparams).lp) AND $
          ((*(*(*info).restoreparams).sel_loops)[i] EQ 1)) THEN BEGIN
				drawrestore += 1
				xp_orig = (*(*(*info).restoreparams).xp_restored)[$
          (0+(*(*(*info).restoreparams).psizes)[i]):((*(*(*info).restoreparams).psizes)[i+1]-1)]
				yp_orig = (*(*(*info).restoreparams).yp_restored)[$
          (0+(*(*(*info).restoreparams).psizes)[i]):((*(*(*info).restoreparams).psizes)[i+1]-1)]
				xr_orig = (*(*(*info).restoreparams).xr_restored)[$
          (0+(*(*(*info).restoreparams).rsizes)[i]):((*(*(*info).restoreparams).rsizes)[i+1]-1)]
				yr_orig = (*(*(*info).restoreparams).yr_restored)[$
          (0+(*(*(*info).restoreparams).rsizes)[i]):((*(*(*info).restoreparams).rsizes)[i+1]-1)]
				sxp = (xp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
          ((*(*info).dataparams).d_nx+1)
				syp = (yp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
          ((*(*info).dataparams).d_ny+1)
				sxr = (xr_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
          ((*(*info).dataparams).d_nx+1)
				syr = (yr_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
          ((*(*info).dataparams).d_ny+1)
				sxp_last = sxp[(SIZE(sxp))[1]-1]-1.5*(*(*info).zooming).factor
				syp_last = syp[(SIZE(sxp))[1]-1]+1.5*(*(*info).zooming).factor
				IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
					PLOTS, sxp, syp, PSYM = 1, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
					PLOTS, sxp, syp, PSYM = 4, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
				ENDIF
				IF ~KEYWORD_SET(NO_NUMBER) THEN XYOUTS,sxp_last,syp_last,STRTRIM(i,2), /DEVICE
				IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN $
          PLOTS, sxr, syr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
				ELSE $
          PLOTS, sxr, syr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR,$
            /DEVICE, THICK=thick
			ENDIF
		ENDFOR
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, [(*(*info).restoreparams).cfilecount,drawrestore], $
                           labels=['Loops restored','Loops drawn']
	ENDIF 
	IF (((*(*info).retrparams).clfilecount GT 0) AND $
      ((*(*info).loopswitch).retrieve_loops EQ 1)) THEN BEGIN
		drawretr = 0
		FOR k=0,(*(*info).retrparams).clfilecount-1 DO BEGIN
			IF ((*(*(*info).retrparams).sel_loops)[k] EQ 1) THEN BEGIN
				drawretr += 1
				xlp_orig = (*(*(*info).retrparams).xlp_restored)[$
          (0+(*(*(*info).retrparams).lpsizes)[k]):((*(*(*info).retrparams).lpsizes)[k+1]-1)]
				ylp_orig = (*(*(*info).retrparams).ylp_restored)[$
          (0+(*(*(*info).retrparams).lpsizes)[k]):((*(*(*info).retrparams).lpsizes)[k+1]-1)]
				xlr_orig = (*(*(*info).retrparams).xlr_restored)[$
          (0+(*(*(*info).retrparams).lrsizes)[k]):((*(*(*info).retrparams).lrsizes)[k+1]-1)]
				ylr_orig = (*(*(*info).retrparams).ylr_restored)[$
          (0+(*(*(*info).retrparams).lrsizes)[k]):((*(*(*info).retrparams).lrsizes)[k+1]-1)]
				sxlp = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
          ((*(*info).dataparams).d_nx+1)
				sylp = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
          ((*(*info).dataparams).d_ny+1)
				sxlr = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
          ((*(*info).dataparams).d_nx+1)
				sylr = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
          ((*(*info).dataparams).d_ny+1)
				sxlp_last = sxlp[(SIZE(sxlp))[1]-1]-1.5*(*(*info).zooming).factor
				sylp_last = sylp[(SIZE(sxlp))[1]-1]+1.5*(*(*info).zooming).factor
				IF ~KEYWORD_SET(NO_ENDPOINTS) THEN BEGIN
					PLOTS, sxlp, sylp, PSYM = 1, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
					PLOTS, sxlp, sylp, PSYM = 4, COLOR = !P.COLOR, /DEVICE, THICK=thick, SYMSIZE=symsize
				ENDIF
				IF ~KEYWORD_SET(NO_NUMBER) THEN XYOUTS,sxlp_last,sylp_last,'L'+STRTRIM(k,2), /DEVICE
				IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN $
          PLOTS, sxlr, sylr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
				ELSE $
          PLOTS, sxlr, sylr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR,$
            /DEVICE, THICK=thick
			ENDIF
		ENDFOR
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, [(*(*info).retrparams).clfilecount,drawretr], $
                           labels=['Loops retrieved','Retrieved loops drawn']
	ENDIF 
	IF ((*(*info).loopswitch).retrieve_detfile EQ 1) THEN BEGIN
		low = (*(*info).detparams).mid-FLOOR((*(*info).detparams).width/2.)
		upp = (*(*info).detparams).mid+FLOOR((*(*info).detparams).width/2.)
		condition = WHERE(*(*(*info).detparams).sel_dets EQ 1, conditioncount)
		IF (*(*info).overlayswitch).det_overlay_all THEN BEGIN
			FOR j=0,(*(*info).detparams).nr_dets-1 DO $
        CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, j, low, upp, NO_NUMBER=no_number, $
                                                  THICK=thick, NO_ENDPOINTS=no_endpoints, $
                                                  SYMSIZE=symsize
			detdrawn = (*(*info).detparams).nr_dets
		ENDIF ELSE IF (((*(*info).overlayswitch).det_overlay_all EQ 0) AND $
                   (N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) THEN BEGIN
				indices = WHERE((*(*(*info).detparams).sel_dets) EQ 1)
				FOR m=0,N_ELEMENTS(condition)-1 DO $
          CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, indices[m], low, upp, $
                                                    NO_NUMBER=no_number, THICK=thick, $
                                                    NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize
				detdrawn = N_ELEMENTS(indices)
		ENDIF ELSE detdrawn = 0
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, [(*(*info).detparams).nr_dets,conditioncount,detdrawn], $
                           labels=['Detections retrieved','Detections selected','Detections drawn']
	ENDIF ELSE RETURN
END

PRO CRISPEX_DRAW_RETRIEVED_DET_LOOP_OVERLAYS, event, j, low, upp, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize
; Handles the actual overplotting of loop paths from the retrieved detections
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	xlp_orig = (*(*(*info).detparams).xlp_restored)[(0+(*(*(*info).detparams).lpsizes)[j]):((*(*(*info).detparams).lpsizes)[j+1]-1)]
	ylp_orig = (*(*(*info).detparams).ylp_restored)[(0+(*(*(*info).detparams).lpsizes)[j]):((*(*(*info).detparams).lpsizes)[j+1]-1)]
	sxlp = (xlp_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
	sylp = (ylp_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
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
		sxlr = (xlr_orig - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
		sylr = (ylr_orig - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
		IF ((*(*info).overlayparams).loop_linestyle EQ 1) THEN PLOTS, sxlr, sylr, PSYM = 3, COLOR = !P.COLOR, /DEVICE, THICK=thick $
		ELSE PLOTS, sxlr, sylr, LINESTYLE = (*(*info).overlayparams).loop_linestyle, COLOR = !P.COLOR,/DEVICE, THICK=thick
	ENDFOR
END

PRO CRISPEX_DRAW_MASK_OVERLAYS, event
; Handles the overlay of a mask
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	x_low = (*(*info).zooming).xpos
	x_upp = (*(*info).zooming).xpos + (*(*info).dataparams).d_nx
	y_low = (*(*info).zooming).ypos
	y_upp = (*(*info).zooming).ypos + (*(*info).dataparams).d_ny
	LOADCT, (*(*info).overlayparams).maskct, /SILENT
	CONTOUR,*(*(*info).data).maskslice,COLOR=(*(*info).overlayparams).maskcolor, LEVELS = 1, /ISOTROPIC, XS=13,YS=13,POSITION=[0,0,1,1],/NORMAL, /NOERASE
;	CONTOUR,*(*(*info).data).maskslice,COLOR=(*(*info).overlayparams).maskcolor-100, LEVELS = 2, /ISOTROPIC, XS=13,YS=13,POSITION=[0,0,1,1],/NORMAL, /NOERASE
	LOADCT, 0, /SILENT
END

PRO CRISPEX_DRAW_RASTER_OVERLAYS, event
; Handles the overlay of raster contours on SJI
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	LOADCT, (*(*info).overlayparams).maskct, /SILENT
  FOR i=0,N_ELEMENTS((*(*info).dispparams).rastercont[*,0])-1 DO BEGIN
    xlow = (*(*info).dispparams).rastercont[i,0] * (*(*info).winsizes).sjiwinx / $
      ((*(*info).dataparams).sjinx)
    xupp = (*(*info).dispparams).rastercont[i,2] * (*(*info).winsizes).sjiwinx / $
      ((*(*info).dataparams).sjinx)
    ylow = (*(*info).dispparams).rastercont[i,1] * (*(*info).winsizes).sjiwiny / $
      ((*(*info).dataparams).sjiny)
    yupp = (*(*info).dispparams).rastercont[i,3] * (*(*info).winsizes).sjiwiny / $
      ((*(*info).dataparams).sjiny)
    ; plot lower vertical boundary
    PLOTS, REPLICATE(xlow,2), [ylow,yupp], COLOR=(*(*info).overlayparams).maskcolor, /DEVICE
    ; plot upper horizontal boundary
    PLOTS, [xlow,xupp], REPLICATE(yupp,2), COLOR=(*(*info).overlayparams).maskcolor, /DEVICE
    ; plot upper vertical boundary
    PLOTS, REPLICATE(xupp,2), [ylow,yupp], COLOR=(*(*info).overlayparams).maskcolor, /DEVICE
    ; plot lower horizontal boundary
    PLOTS, [xlow,xupp], REPLICATE(ylow,2), COLOR=(*(*info).overlayparams).maskcolor, /DEVICE
  ENDFOR
	LOADCT, 0, /SILENT
END

PRO CRISPEX_DRAW_GET_SPECTRAL_AXES, event, MAIN=main, REFERENCE=reference
; Handles determination of TICKINTERVAL of spectral axes
  ; should be called whenever changing the number of displayed diagnostics
  ; AND on startup
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Write current variables to dummy variables depending on case
  IF KEYWORD_SET(MAIN) THEN BEGIN
    disp_diagnostics = (*(*info).intparams).disp_diagnostics
    ndisp_diagnostics = (*(*info).intparams).ndisp_diagnostics
    diag_width = (*(*info).intparams).diag_width
    diag_start = (*(*info).intparams).diag_start
    lps = (*(*info).dataparams).lps
    v_dop = (*(*info).plotaxes).v_dop
  ENDIF ELSE IF KEYWORD_SET(REFERENCE) THEN BEGIN
    disp_diagnostics = (*(*info).intparams).disp_refdiagnostics
    ndisp_diagnostics = (*(*info).intparams).ndisp_refdiagnostics
    diag_width = (*(*info).intparams).refdiag_width
    diag_start = (*(*info).intparams).refdiag_start
    lps = (*(*info).dataparams).reflps
    v_dop = (*(*info).plotaxes).v_dop_ref
  ENDIF
  ; Determine settings from displayed diagnostics
  wheredisp = WHERE(disp_diagnostics EQ 1)
  diag_widths = diag_width[wheredisp]
  diag_starts = diag_start[wheredisp]
  diag_ratio = diag_widths / FLOAT(TOTAL(diag_widths))
  IF (ndisp_diagnostics GT 1) THEN BEGIN
    lambda_widths = FLTARR(ndisp_diagnostics)
    dop_widths = FLTARR(ndisp_diagnostics)
    FOR d=0,ndisp_diagnostics-1 DO BEGIN
      lambda_widths[d] = lps[diag_starts[d]+diag_widths[d]-1] - lps[diag_starts[d]]
      dop_widths[d] = (*v_dop[wheredisp[d]])[diag_widths[d]-1] - $
        (*v_dop[wheredisp[d]])[0]
    ENDFOR
    order = FLOOR(ALOG10(lambda_widths))
    int_widths = FLOOR(lambda_widths/10^FLOAT(order))*10^FLOAT(order)
    wheremax = WHERE(int_widths EQ MAX(int_widths, /NAN))
    ; Have at least two tickmarks per range for the biggest range
    xtickinterval = int_widths[wheremax]/2.
    doporder = FLOOR(ALOG10(dop_widths))
    int_dopwidths = FLOOR(dop_widths/10^FLOAT(doporder))*10^FLOAT(doporder)
    wheredopmax = WHERE(int_dopwidths EQ MAX(int_dopwidths, /NAN))
    ; Have at least two tickmarks per range for the biggest range
    xdoptickinterval = int_dopwidths[wheredopmax]/2.
  ENDIF ELSE BEGIN
    xtickinterval = 0
    xdoptickinterval = 0
  ENDELSE
  ; Save results to appropriate variables
  IF KEYWORD_SET(MAIN) THEN BEGIN
    (*(*info).plotaxes).xtickinterval = xtickinterval[0]
    (*(*info).plotaxes).xdoptickinterval = xdoptickinterval[0]
    ; Determine proportional spectral window sizes
    *(*(*info).intparams).wheredispdiag = wheredisp
    *(*(*info).intparams).diag_widths = diag_widths
    *(*(*info).intparams).diag_starts = diag_starts 
    *(*(*info).plotaxes).diag_ratio = diag_ratio 
    IF (*(*info).winswitch).showsp THEN $
      *(*(*info).plotaxes).diag_range_sp = *(*(*info).plotaxes).diag_ratio * $
        (*(*info).plotpos).xplspw 
    IF (*(*info).winswitch).showphis THEN $
      *(*(*info).plotaxes).diag_range_phis = *(*(*info).plotaxes).diag_ratio * $
        (*(*info).plotpos).phisxplspw
  ENDIF ELSE IF KEYWORD_SET(REFERENCE) THEN BEGIN
    (*(*info).plotaxes).xreftickinterval = xtickinterval[0]
    (*(*info).plotaxes).xrefdoptickinterval = xdoptickinterval[0]
    ; Determine proportional spectral window sizes
    *(*(*info).intparams).wheredisprefdiag = wheredisp 
    *(*(*info).intparams).refdiag_widths = diag_widths
    *(*(*info).intparams).refdiag_starts = diag_starts 
    *(*(*info).plotaxes).refdiag_ratio = diag_ratio 
    IF (*(*info).winswitch).showrefsp THEN $
      *(*(*info).plotaxes).refdiag_range_sp = *(*(*info).plotaxes).refdiag_ratio * $
        (*(*info).plotpos).refxplspw 
  ENDIF
END

PRO CRISPEX_DRAW, event, NO_MAIN=no_main, NO_REF=no_ref, NO_PHIS=no_phis
; Handles the actual drawing of the data into the respective open display windows
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).curs).lockset AND ((*(*info).overlayswitch).loopslit NE 1)) THEN BEGIN
		(*(*info).curs).sx = (*(*info).curs).sxlock	
    (*(*info).curs).sy = (*(*info).curs).sylock
	ENDIF 
  CRISPEX_DRAW_FEEDBPARAMS, event
	CRISPEX_DRAW_IMREF, event, NO_MAIN=no_main, NO_REF=no_ref
	IF ((*(*info).winswitch).showls OR (*(*info).winswitch).showsp OR $
    (*(*info).winswitch).showrefls OR (*(*info).winswitch).showrefsp) THEN $
      CRISPEX_DRAW_SPECTRAL, event, NO_MAIN=no_main, NO_REF=no_ref, NO_PHIS=no_phis
	IF ((*(*info).winswitch).showloop OR (*(*info).winswitch).showrefloop OR $
    (*(*info).winswitch).showrestloop OR (*(*info).winswitch).showretrdet) THEN $
      CRISPEX_DRAW_TIMESLICES, event
	IF (*(*info).winswitch).showint THEN CRISPEX_DRAW_INT, event
END

PRO CRISPEX_DRAW_FEEDBPARAMS, event
; Prints all feedback parameters to appropriate fields
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Position parameters
  ; Main
  xycoord_txt = '('+STRING(LONG((*(*info).dataparams).x),$
    FORMAT=(*(*info).paramparams).xcoord_format)+$
    ','+STRING(LONG((*(*info).dataparams).y),$
    FORMAT=(*(*info).paramparams).ycoord_format)+')'
  WIDGET_CONTROL, (*(*info).ctrlsparam).xycoord_val, $
    SET_VALUE=xycoord_txt
  xycoord_real_txt = '('+STRING(FLOAT((*(*info).dataparams).x*(*(*info).dataparams).dx),$
    FORMAT=(*(*info).paramparams).xcoord_real_format)+$
    ','+STRING(FLOAT((*(*info).dataparams).y*(*(*info).dataparams).dy),$
    FORMAT=(*(*info).paramparams).ycoord_real_format)+')'
  WIDGET_CONTROL, (*(*info).ctrlsparam).xycoord_real_val, $
    SET_VALUE=xycoord_real_txt
  ; Reference
  IF (*(*info).winswitch).showref THEN BEGIN
    refxycoord_txt = '('+STRING(LONG((*(*info).dataparams).x),$
      FORMAT=(*(*info).paramparams).refxcoord_format)+$
      ','+STRING(LONG((*(*info).dataparams).y),$
      FORMAT=(*(*info).paramparams).refycoord_format)+')'
    WIDGET_CONTROL, (*(*info).ctrlsparam).refxycoord_val, $
      SET_VALUE=refxycoord_txt
    refxycoord_real_txt = '('+STRING(FLOAT($
      (*(*info).dataparams).x*(*(*info).dataparams).dx),$
      FORMAT=(*(*info).paramparams).refxcoord_real_format)+$
      ','+STRING(FLOAT((*(*info).dataparams).y*(*(*info).dataparams).dy),$
      FORMAT=(*(*info).paramparams).refycoord_real_format)+')'
    WIDGET_CONTROL, (*(*info).ctrlsparam).refxycoord_real_val, $
      SET_VALUE=refxycoord_real_txt
  ENDIF
  ; SJI
  IF (*(*info).winswitch).showsji THEN BEGIN
    sjixycoord_txt = '('+STRING(LONG((*(*info).dataparams).xsji),$
      FORMAT=(*(*info).paramparams).sjixcoord_format)+$
      ','+STRING(LONG((*(*info).dataparams).ysji),$
      FORMAT=(*(*info).paramparams).sjiycoord_format)+')'
    WIDGET_CONTROL, (*(*info).ctrlsparam).sjixycoord_val, $
      SET_VALUE=sjixycoord_txt
    sjixycoord_real_txt = '('+STRING(FLOAT($
      (*(*info).dataparams).xsji*(*(*info).dataparams).sjidx),$
      FORMAT=(*(*info).paramparams).sjixcoord_real_format)+$
      ','+STRING(FLOAT((*(*info).dataparams).ysji*(*(*info).dataparams).sjidy),$
      FORMAT=(*(*info).paramparams).sjiycoord_real_format)+')'
    WIDGET_CONTROL, (*(*info).ctrlsparam).sjixycoord_real_val, $
      SET_VALUE=sjixycoord_real_txt
  ENDIF

  ; Spectral parameters
  ; Main
  lp_idx_txt = STRING((*(*info).dataparams).lp, FORMAT=(*(*info).paramparams).lp_idx_format)
  WIDGET_CONTROL, (*(*info).ctrlsparam).lp_idx_val, SET_VALUE=lp_idx_txt
  IF (*(*info).plotswitch).v_dop_set THEN BEGIN
    lp_real_txt = STRING((*(*info).dataparams).lps[(*(*info).dataparams).lp], $
      FORMAT=(*(*info).paramparams).lp_real_format)
    lp_vdop_txt = STRING((*(*(*info).plotaxes).v_dop[(*(*info).intparams).lp_diag_all])[$
      (*(*info).dataparams).lp-(*(*info).intparams).diag_start[$
      (*(*info).intparams).lp_diag_all]], FORMAT=(*(*info).paramparams).lp_vdop_format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).lp_real_val, SET_VALUE=lp_real_txt
    WIDGET_CONTROL, (*(*info).ctrlsparam).lp_vdop_val, SET_VALUE=lp_vdop_txt
  ENDIF
  ; Reference
  IF ((*(*info).dataparams).refnlp GT 1) THEN BEGIN
    lp_ref_idx_txt = STRING((*(*info).dataparams).lp_ref, $
      FORMAT=(*(*info).paramparams).lp_ref_idx_format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).lp_ref_idx_val, SET_VALUE=lp_ref_idx_txt
    IF (*(*info).plotswitch).v_dop_set_ref THEN BEGIN
      lp_ref_real_txt = STRING((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref], $
        FORMAT=(*(*info).paramparams).lp_ref_real_format)
      lp_ref_vdop_txt = STRING((*(*(*info).plotaxes).v_dop_ref[$
        (*(*info).intparams).lp_ref_diag_all])[(*(*info).dataparams).lp_ref-$
        (*(*info).intparams).refdiag_start[(*(*info).intparams).lp_ref_diag_all]], $
        FORMAT=(*(*info).paramparams).lp_ref_vdop_format)
      WIDGET_CONTROL, (*(*info).ctrlsparam).lp_ref_real_val, SET_VALUE=lp_ref_real_txt
      WIDGET_CONTROL, (*(*info).ctrlsparam).lp_ref_vdop_val, SET_VALUE=lp_ref_vdop_txt
    ENDIF
  ENDIF

  ; Time parameters
  ; Main
  t_idx_txt = STRING(LONG((*(*info).dispparams).t_main),$
    FORMAT=(*(*info).paramparams).t_idx_format)
  WIDGET_CONTROL, (*(*info).ctrlsparam).t_idx_val, SET_VALUE=t_idx_txt
  IF (*(*info).paramswitch).dt_set THEN BEGIN
    t_real_txt = STRING((*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t],$
      FORMAT=(*(*info).paramparams).t_real_format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).t_real_val, SET_VALUE=t_real_txt
  ENDIF
  IF (*(*info).paramswitch).t_raster THEN BEGIN
    t_raster_real_txt = STRING((*(*info).dataparams).tarr_raster_main[(*(*info).dataparams).x,$
      (*(*info).dispparams).t_main], FORMAT=(*(*info).paramparams).t_raster_real_format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).t_raster_real_val, SET_VALUE=t_raster_real_txt
  ENDIF
  ; Reference
  IF (*(*info).winswitch).showref THEN BEGIN
    IF ((*(*info).dataparams).refnt GT 1) THEN BEGIN
      t_ref_idx_txt = STRING(LONG((*(*info).dispparams).t_ref),$
        FORMAT=(*(*info).paramparams).t_ref_idx_format)
      WIDGET_CONTROL, (*(*info).ctrlsparam).t_ref_idx_val, SET_VALUE=t_ref_idx_txt
      IF (*(*info).paramswitch).dt_set THEN BEGIN
        t_ref_real_txt = STRING((*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t],$
          FORMAT=(*(*info).paramparams).t_ref_real_format)
        WIDGET_CONTROL, (*(*info).ctrlsparam).t_ref_real_val, SET_VALUE=t_ref_real_txt
      ENDIF
      IF (*(*info).paramswitch).t_raster THEN BEGIN
        t_raster_ref_real_txt = STRING((*(*info).dataparams).tarr_raster_ref[$
          (*(*info).dataparams).x,(*(*info).dispparams).t_ref], $
          FORMAT=(*(*info).paramparams).t_raster_ref_real_format)
        WIDGET_CONTROL, (*(*info).ctrlsparam).t_raster_ref_real_val, $
          SET_VALUE=t_raster_ref_real_txt
      ENDIF
    ENDIF
  ENDIF
  ; SJI
  IF (*(*info).winswitch).showsji THEN BEGIN
    IF ((*(*info).dataparams).sjint GT 1) THEN BEGIN
      t_sji_idx_txt = STRING(LONG((*(*info).dispparams).t_sji), $
        FORMAT=(*(*info).paramparams).t_sji_idx_format)
      WIDGET_CONTROL, (*(*info).ctrlsparam).t_sji_idx_val, SET_VALUE=t_sji_idx_txt
      IF (*(*info).paramswitch).dt_set THEN BEGIN
        t_sji_real_txt = STRING((*(*(*info).dispparams).tarr_sji)[(*(*info).dispparams).t],$
          FORMAT=(*(*info).paramparams).t_sji_real_format)
        WIDGET_CONTROL, (*(*info).ctrlsparam).t_sji_real_val, SET_VALUE=t_sji_real_txt
      ENDIF
    ENDIF
  ENDIF
  
  ; Data values parameters
  ; Main
  datadims = SIZE(*(*(*info).data).xyslice,/N_DIMENSIONS)
  IF (datadims EQ 2) THEN $
    act_dataval = (*(*(*info).data).xyslice)[LONG((*(*info).dataparams).x), $
      LONG((*(*info).dataparams).y)] $
  ELSE $  ; Failsafe for IRIS sit-and-stare
    act_dataval = (*(*(*info).data).xyslice)[LONG((*(*info).dataparams).y)]
  IF ((FINITE(act_dataval) EQ 1) AND (act_dataval NE 0)) THEN $
    order = FLOOR(ALOG10(ABS(act_dataval))) $
  ELSE $
    order = 0
  IF ((order LE -2) OR (order GE 3)) THEN format = '(E10.4)' ELSE format = '(F9.2)'
  dataval_real_txt = STRING(act_dataval, FORMAT=format) 
  WIDGET_CONTROL, (*(*info).ctrlsparam).dataval_real_val, SET_VALUE=dataval_real_txt
  ; Reference
  IF (*(*info).winswitch).showref THEN BEGIN
    datadims = SIZE(*(*(*info).data).refslice,/N_DIMENSIONS)
    IF (datadims EQ 2) THEN $
      act_ref_dataval = (*(*(*info).data).refslice)[LONG((*(*info).dataparams).x), $
        LONG((*(*info).dataparams).y)] $
    ELSE $  ; Failsafe for IRIS sit-and-stare
      act_ref_dataval = (*(*(*info).data).refslice)[LONG((*(*info).dataparams).y)]
    IF ((FINITE(act_ref_dataval) EQ 1) AND (act_ref_dataval NE 0)) THEN $
      order = FLOOR(ALOG10(ABS(act_ref_dataval))) $
    ELSE $
      order = 0
    IF ((order LE -2) OR (order GE 3)) THEN format = '(E10.4)' ELSE format = '(F9.2)'
    dataval_ref_real_txt = STRING(act_ref_dataval, FORMAT=format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).dataval_ref_real_val, $
      SET_VALUE=dataval_ref_real_txt
  ENDIF
  ; SJI
  IF (*(*info).winswitch).showsji THEN BEGIN
    act_sji_dataval = (*(*(*info).data).sjislice)[LONG((*(*info).dataparams).xsji), $
      LONG((*(*info).dataparams).ysji)]
    IF ((FINITE(act_sji_dataval) EQ 1) AND (act_sji_dataval NE 0)) THEN $
      order = FLOOR(ALOG10(ABS(act_sji_dataval))) $
    ELSE $
      order = 0
    IF ((order LE -2) OR (order GE 3)) THEN format = '(E10.4)' ELSE format = '(F9.2)'
    dataval_sji_real_txt = STRING(act_sji_dataval, FORMAT=format)
    WIDGET_CONTROL, (*(*info).ctrlsparam).dataval_sji_real_val, $
      SET_VALUE=dataval_sji_real_txt
  ENDIF

  ; Zoom value
	WIDGET_CONTROL, (*(*info).ctrlsparam).zoom_val, $
    SET_VALUE = STRING((*(*info).zooming).factor*100.,FORMAT='(I4)')+'%'
END

PRO CRISPEX_DRAW_IMREF, event, NO_MAIN=no_main, NO_REF=no_ref
; (Re)draw main and reference image window procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ~KEYWORD_SET(NO_MAIN) THEN CRISPEX_DRAW_XY, event
	IF (~KEYWORD_SET(NO_REF) AND (*(*info).winswitch).showref) THEN CRISPEX_DRAW_REF, event
	IF (*(*info).winswitch).showimref THEN CRISPEX_DRAW_IMREF_BLINK, event
	IF (*(*info).winswitch).showdop THEN CRISPEX_DRAW_DOPPLER, event
  IF (*(*info).winswitch).showsji THEN CRISPEX_DRAW_SJI, event
END

PRO CRISPEX_DRAW_SPECTRAL, event, NO_MAIN=no_main, NO_REF=no_ref, NO_PHIS=no_phis
; (Re)draw spectral windows procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (~KEYWORD_SET(NO_MAIN) AND $
    ((*(*info).winswitch).showls OR (*(*info).winswitch).showsp)) THEN $
      CRISPEX_DRAW_SPECTRAL_MAIN, event
	IF (~KEYWORD_SET(NO_REF) AND $
    ((*(*info).winswitch).showrefls OR (*(*info).winswitch).showrefsp)) THEN $
      CRISPEX_DRAW_SPECTRAL_REF, event
	IF (~KEYWORD_SET(NO_PHIS) AND $
      ((*(*info).winswitch).showphis AND (((*(*info).pbparams).mode EQ 'PAUSE') OR $
      (*(*info).dispparams).phislice_update))) THEN CRISPEX_DRAW_PHIS, event		
END

PRO CRISPEX_DRAW_TIMESLICES, event
; (Re)draw timeslices procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).winswitch).showloop THEN CRISPEX_DRAW_LOOPSLAB, event 
	IF (*(*info).winswitch).showrefloop THEN CRISPEX_DRAW_REFLOOPSLAB, event 
	IF (*(*info).winswitch).showrestloop THEN CRISPEX_DRAW_REST_LOOP, event
	IF (*(*info).winswitch).showretrdet THEN CRISPEX_DRAW_RETR_DET, event
END

PRO CRISPEX_DRAW_SUBCOLOR, event, imref, subcolor, minimum, maximum, XYRANGE=xyrange,$
      MAIN2SJI=main2sji
; Determines the color beneath the cursor to get the cursor color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (N_ELEMENTS(XYRANGE) EQ 0) THEN BEGIN
    IF KEYWORD_SET(MAIN2SJI) THEN BEGIN
      CRISPEX_COORDS_TRANSFORM_XY, event, /MAIN2SJI
  		x_upp = LONG((*(*info).dataparams).xsji) + 5 < LONG((*(*info).dispparams).xsji_last)
  		x_low = LONG((*(*info).dataparams).xsji) - 5 > LONG((*(*info).dispparams).xsji_first)
  		y_upp = LONG((*(*info).dataparams).ysji) + 5 < LONG((*(*info).dispparams).ysji_last)
  		y_low = LONG((*(*info).dataparams).ysji) - 5 > LONG((*(*info).dispparams).ysji_first)
    ENDIF
		x_upp = LONG((*(*info).dataparams).x) + 5 < LONG((*(*info).dispparams).x_last)
		x_low = LONG((*(*info).dataparams).x) - 5 > LONG((*(*info).dispparams).x_first)
		y_upp = LONG((*(*info).dataparams).y) + 5 < LONG((*(*info).dispparams).y_last)
		y_low = LONG((*(*info).dataparams).y) - 5 > LONG((*(*info).dispparams).y_first)
	ENDIF ELSE BEGIN
		x_low = xyrange[0]	&	x_upp = xyrange[1]
		y_low = xyrange[2]	&	y_upp = xyrange[3]
	ENDELSE
  scale_idx = ((imref EQ 0) OR (imref EQ 2)) * (*(*info).intparams).lp_diag_all + $
    ((imref GT 0)+(imref GT 2)) * (*(*info).intparams).ndiagnostics + $
    (imref EQ 1) * (*(*info).intparams).lp_ref_diag_all + $
    (imref GT 1) * (*(*info).intparams).nrefdiagnostics ;+ $
;    (imref GT 2)
	IF (imref EQ 0) THEN $
		selected_data = (*(*(*info).data).imagedata)[$
      (*(*info).dispparams).t_main * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
      (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp]$
	ELSE IF (imref EQ 1) THEN BEGIN
		IF ((*(*info).dataparams).refnt GT 1) THEN $
      selected_data = (*(*(*info).data).refdata)[(*(*info).dispparams).t_ref *  $
        (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref]$
    ELSE $
			selected_data = (*(*(*info).data).refdata)[(*(*info).dataparams).lp_ref]
	ENDIF ELSE IF (imref EQ 2) THEN	$
    selected_data = *(*(*info).data).dopslice $
  ELSE IF (imref EQ 3) THEN $
    selected_data = *(*(*info).data).sjislice 
  minimum = MIN(selected_data, MAX=maximum, /NAN)
  minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
    (*(*info).scaling).minimum[scale_idx],(*(*info).scaling).maximum[scale_idx])
	scaled_data = BYTSCL(selected_data[x_low:x_upp, y_low:y_upp], MIN=minmax[0], MAX=minmax[1], /NAN)
	subcolor = MEAN(scaled_data, /NAN)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [x_low,x_upp,y_low,y_upp,subcolor], $
      labels=['x_low','x_upp','y_low','y_upp','subcolor']
END

PRO CRISPEX_DRAW_SCALING, event, finalimage, minimum, maximum, $
  MAIN=main, DOPPLER=doppler, REFERENCE=reference, SJI=sji, $
  SELECTED_DATA=selected_data
; Determines the minimum and maximum value for the image scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	x_low = (*(*info).zooming).xpos
	y_low = (*(*info).zooming).ypos
  d_nx = (*(*info).dataparams).d_nx
  d_ny = (*(*info).dataparams).d_ny
  x_upp = (*(*info).zooming).xpos + d_nx
  y_upp = (*(*info).zooming).ypos + d_ny
  ; sel = 0 -> main
  ; sel = 1 -> reference
  ; sel = 2 -> doppler
  ; sel = 3 -> sji
  ; imagescale = 0 -> based on first
  ; imagescale = 1 -> based on current
  ; imagescale = 2 -> per time step
  IF KEYWORD_SET(MAIN) THEN sel = 0 $
    ELSE IF KEYWORD_SET(REFERENCE) THEN sel = 1 $
    ELSE IF KEYWORD_SET(DOPPLER) THEN sel = 2 $
    ELSE IF KEYWORD_SET(SJI) THEN sel = 3
  scale_idx = ((sel EQ 0) OR (sel EQ 2)) * (*(*info).intparams).lp_diag_all + $
    ((sel GT 0) + (sel GT 2)) * (*(*info).intparams).ndiagnostics + $
    (sel EQ 1) * (*(*info).intparams).lp_ref_diag_all + $
    (sel GT 1) * (*(*info).intparams).nrefdiagnostics ;+ $
;    (sel GT 2)
  IF KEYWORD_SET(MAIN) THEN BEGIN
   ; sel = 0
    datadims = SIZE(*(*(*info).data).xyslice,/N_DIMENSIONS)
    IF (N_ELEMENTS(SELECTED_DATA) LT 1) THEN BEGIN
      IF (datadims EQ 2) THEN $
        selected_data = (*(*(*info).data).xyslice)[x_low:x_upp,y_low:y_upp] $
      ELSE $
        selected_data = (*(*(*info).data).xyslice)[y_low:y_upp]
    ENDIF
		IF ((*(*(*info).scaling).imagescale)[sel] EQ 0) THEN BEGIN
      minimum = (*(*info).scaling).imagemin
      maximum = (*(*info).scaling).imagemax
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[sel] EQ 1) THEN BEGIN
      minimum = (*(*info).scaling).imagemin_curr
      maximum = (*(*info).scaling).imagemax_curr
    ENDIF
  ENDIF ELSE IF KEYWORD_SET(REFERENCE) THEN BEGIN
    ; sel = 1
    datadims = SIZE(*(*(*info).data).refslice,/N_DIMENSIONS)
    IF (N_ELEMENTS(SELECTED_DATA) LT 1) THEN BEGIN
      IF (datadims EQ 2) THEN $
        selected_data = (*(*(*info).data).refslice)[x_low:x_upp,y_low:y_upp] $
      ELSE $
        selected_data = (*(*(*info).data).refslice)[y_low:y_upp]
    ENDIF
  	IF ((*(*(*info).scaling).imagescale)[sel] EQ 0) THEN BEGIN
      minimum = (*(*info).scaling).refmin
      maximum = (*(*info).scaling).refmax
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[sel] EQ 1) THEN BEGIN
      minimum = (*(*info).scaling).refmin_curr
      maximum = (*(*info).scaling).refmax_curr
    ENDIF 
  ENDIF ELSE IF KEYWORD_SET(DOPPLER) THEN BEGIN
    ; sel = 2
    datadims = SIZE(*(*(*info).data).dopslice,/N_DIMENSIONS)
    IF (N_ELEMENTS(SELECTED_DATA) LT 1) THEN BEGIN
      IF (datadims EQ 2) THEN $
        selected_data = (*(*(*info).data).dopslice)[x_low:x_upp,y_low:y_upp] $
      ELSE $
        selected_data = (*(*(*info).data).dopslice)[y_low:y_upp]
    ENDIF
		IF ((*(*(*info).scaling).imagescale)[sel] EQ 0) THEN BEGIN
      minimum = (*(*info).scaling).dopmin
      maximum = (*(*info).scaling).dopmax
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[sel] EQ 1) THEN BEGIN
      minimum = (*(*info).scaling).dopmin_curr
      maximum = (*(*info).scaling).dopmax_curr
    ENDIF
  ENDIF ELSE IF KEYWORD_SET(SJI) THEN BEGIN
   ; sel = 3
    IF (N_ELEMENTS(SELECTED_DATA) LT 1) THEN selected_data = *(*(*info).data).sjislice 
    IF ((*(*(*info).scaling).imagescale)[sel] EQ 0) THEN BEGIN
      minimum = (*(*info).scaling).sjimin
      maximum = (*(*info).scaling).sjimax
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[sel] EQ 1) THEN BEGIN
      minimum = (*(*info).scaling).sjimin_curr
      maximum = (*(*info).scaling).sjimax_curr
    ENDIF
  ENDIF
  IF ((*(*info).scaling).gamma[scale_idx] NE 1.) THEN BEGIN
    wherelt0 = WHERE(selected_data LT 0)
    IF (wherelt0[0] EQ -1) THEN $
      selected_data = (TEMPORARY(selected_data))^(*(*info).scaling).gamma[scale_idx] $
    ELSE BEGIN
      selected_data = (TEMPORARY(ABS(selected_data)))^(*(*info).scaling).gamma[scale_idx]
      selected_data[wherelt0] *= -1
    ENDELSE
  ENDIF
  IF ((*(*(*info).scaling).imagescale)[sel] EQ 2) THEN BEGIN
;    ; Copied over MISSING handling from modified HISTO_OPT() since it gives  "Floating illegal
;    ; operand" errors otherwise
;    finitvals = FINITE(selected_data)
;    selected_data[WHERE(finitvals eq 0)]=-32768
;    selected_data = selected_data[WHERE(selected_data ne -32768)]
    selected_data = HISTO_OPT(TEMPORARY(selected_data), $
      (*(*info).scaling).histo_opt_val[scale_idx], MISSING=-32768)
    minimum = MIN(selected_data, MAX=maximum, /NAN)
  ENDIF
  minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
    (*(*info).scaling).minimum[scale_idx],(*(*info).scaling).maximum[scale_idx])
	finalimage = BYTSCL(selected_data, MIN = minmax[0], MAX = minmax[1], /NAN) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [minimum,maximum], labels=['minimum','maximum']
END

PRO CRISPEX_DRAW_XY, event, no_cursor=no_cursor, no_number=no_number, thick=thick, no_endpoints=no_endpoints, symsize=symsize, asecbar=asecbar
; (Re)draw main image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).imwid
	CRISPEX_DRAW_SCALING, event, imdisp, minimum, maximum, /MAIN
  TV, CONGRID(imdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny)
	CRISPEX_DRAW_SUBCOLOR, event, 0, subcolor, minimum, maximum
	IF (subcolor GE 122) THEN curscolor = 0 ELSE curscolor = 255
	CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, NO_CURSOR=no_cursor, NO_NUMBER=no_number, $
    THICK=thick, NO_ENDPOINTS=no_endpoints, SYMSIZE=symsize, $
		DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[0])
	IF KEYWORD_SET(ASECBAR) THEN BEGIN
		xlow = 25							& 	ylow = 25 
		xupp = xlow + (*(*info).savparams).overlays_asecbar_pix		
    yupp = ylow + 15 + (*(*info).savparams).overlays_thick
		asecbarcol = 255
		PLOTS, [xlow,xupp],[ylow,ylow], THICK=(*(*info).savparams).overlays_thick, COLOR=asecbarcol, $
      /DEVICE
		XYOUTS, (xupp-xlow)/2.+xlow, ylow+5, $
      STRTRIM((*(*info).savparams).overlays_asecbar_length,2)+'"', /DEVICE, COLOR=asecbarcol, $
      ALIGN=0.5,	CHARSIZE=(*(*info).savparams).overlays_symsize, $
      CHARTHICK=(*(*info).savparams).overlays_thick
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).imwid,curscolor], $
      labels=['Window ID for draw','Main curscolor']
END

PRO CRISPEX_DRAW_DOPPLER, event
; (Re)draw Doppler-image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).dopwid
	IF (*(*info).dispswitch).drawdop THEN BEGIN
		CRISPEX_DRAW_SCALING, event, dopdisp, dopminimum, dopmaximum, /DOPPLER
    TV, CONGRID(dopdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny)
		CRISPEX_DRAW_SUBCOLOR, event, 2, dopsubcolor, dopminimum, dopmaximum
		IF (dopsubcolor GE 122) THEN $
      dopcurscolor = 0 $
    ELSE $
      dopcurscolor = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, dopcurscolor, $
      DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[2])
	ENDIF ELSE BEGIN
		TV, CONGRID(*(*(*info).data).emptydopslice,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny)
		IF ((*(*info).dataparams).lp GT (*(*info).dataparams).lp_dop) THEN BEGIN
			lp_blue = (*(*info).dataparams).lp_dop		&	lp_red = (*(*info).dataparams).lp
		ENDIF ELSE BEGIN
			lp_blue = (*(*info).dataparams).lp		&	lp_red = (*(*info).dataparams).lp_dop
		ENDELSE
		IF (lp_red EQ lp_blue) THEN $
      extramessage = 'Same spectral position' $
    ELSE $
      extramessage = 'Spectral position outside set spectral range'
		XYOUTS,(*(*info).winsizes).xywinx/2.,(*(*info).winsizes).xywiny/2.,$
      'Could not create Doppler image for selected spectral positions.!C'+extramessage+': (lp_blue,lp_red)=('+$
			STRTRIM(lp_blue,2)+','+STRTRIM(lp_red,2)+')', CHARSIZE = 1.2, COLOR = 255, ALIGNMENT = 0.5, /DEVICE
		dopcurscolor = 255
	ENDELSE
;  IF ((*(*info).scaling).imrefscaling EQ 2) THEN BEGIN
;    CRISPEX_SCALING_SET_BUTTONS, event
;    CRISPEX_SCALING_SET_SLIDERS, event
;  ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).dopwid,(*(*info).dispswitch).drawdop,dopcurscolor], $
		labels=['Window ID for draw','Drawing Doppler image','Doppler curscolor']
END

PRO CRISPEX_DRAW_REF, event
; (Re)draw reference image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).refwid
	CRISPEX_DRAW_SCALING, event, refdisp, refmin, refmax, /REFERENCE
  TV, CONGRID(refdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) 
	CRISPEX_DRAW_SUBCOLOR, event, 1, subcolor_ref, refmin, refmax
	IF (subcolor_ref GE 122) THEN curscolor_ref = 0 ELSE curscolor_ref = 255
	CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor_ref, DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[1])
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refwid,curscolor_ref], labels=['Window ID for draw','Reference curscolor']
END

PRO CRISPEX_DRAW_IMREF_BLINK, event
; Handles the blinking of the main and reference image
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET,(*(*info).winids).imrefwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).imrefwid], labels=['Window ID for draw']
	IF (*(*info).winids).imrefdisp THEN BEGIN
		CRISPEX_DRAW_SCALING, event, refdisp, refmin, refmax, /REFERENCE
    TV, CONGRID(refdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) 
		CRISPEX_DRAW_SUBCOLOR, event, 1, subcolor_ref, refmin, refmax
		IF (subcolor_ref GE 122) THEN curscolor_ref = 0 ELSE curscolor_ref = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor_ref, $
      DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[1])
		CRISPEX_DRAW_SUBCOLOR, event, 1, color_reftxt, refmin, refmax;, xyrange=[10,100,5,20]
		IF (color_reftxt GE 122) THEN reftxtcol = 0 ELSE reftxtcol = 255
    label = 'Reference'
    time_val = (*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t]
	ENDIF ELSE BEGIN
		CRISPEX_DRAW_SCALING, event, imdisp, minimum, maximum, /MAIN
    TV, CONGRID(imdisp,(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny) 
		CRISPEX_DRAW_SUBCOLOR, event, 0, subcolor, minimum, maximum
		IF (subcolor GE 122) THEN curscolor = 0 ELSE curscolor = 255
		CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor, $
      DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[0])
		CRISPEX_DRAW_SUBCOLOR, event, 0, color_txt, minimum, maximum;, xyrange=[10,70,5,20]
		IF (color_txt GE 122) THEN txtcol = 0 ELSE txtcol = 255
    label = 'Main'
    time = (*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t]
	ENDELSE
	IF (SIZE((*(*info).plotaxes).dt,/TYPE) NE 2) THEN $ 
    time = STRING(time_val, FORMAT='(F'+STRTRIM(FLOOR(ALOG10(time_val))+4,2)+'.2)')+' s' $
  ELSE $
    time = STRTRIM(time_val, 2)
  XYOUTS, 10, 10, label+' image, t='+time,/DEVICE, COLOR=txtcol
END

PRO CRISPEX_DRAW_SJI, event
; (Re)draw reference image procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).sjiwid
;	sjimin = MIN( *(*(*info).data).sjislice, MAX=sjimax) 
	CRISPEX_DRAW_SCALING, event, sjidisp, minimum, maximum, /SJI
  TV, CONGRID(sjidisp,(*(*info).winsizes).sjiwinx, (*(*info).winsizes).sjiwiny)
  ; Determine cursor colour and overplot cursors and masks
	CRISPEX_DRAW_SUBCOLOR, event, 3, subcolor_sji, sjimin, sjimax, /MAIN2SJI
	IF (subcolor_sji GE 122) THEN curscolor_sji = 0 ELSE curscolor_sji = 255
	CRISPEX_DRAW_CURSCROSS_PLOT, event, curscolor_sji, $
    DRAW_MASK=((*(*info).overlayswitch).mask AND ((*(*info).overlayswitch).maskim)[1]), /SJI
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).sjiwid,curscolor_sji], $
                          labels=['Window ID for draw','Slit-jaw image curscolor']
END

PRO CRISPEX_DRAW_INT, event
; (Re)draw intensity versus time plot procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).intwid
;	title = 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s]
	int_low_y = (*(*(*info).plotaxes).int_low_y)[(*(*info).dataparams).s] 
	int_upp_y = (*(*(*info).plotaxes).int_upp_y)[(*(*info).dataparams).s]
	condition = WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)
	PLOT, FINDGEN((*(*info).dataparams).nt)*(*(*info).plotaxes).dt, $
    FINDGEN((*(*info).dataparams).nt), /NODATA, /NORM, CHARSIZE=1, YR=[int_low_y, int_upp_y], $
    XR = [(*(*info).plotaxes).int_low_t*(*(*info).plotaxes).dt, $
    (*(*info).plotaxes).int_upp_t*(*(*info).plotaxes).dt], /YS, /XS, $
    XTITLE = (*(*info).plottitles).spytitle, YTITLE = 'Counts/Mean Counts', $
    BACKGROUND = (*(*info).plotparams).bgplotcol, COLOR = (*(*info).plotparams).plotcol, $
    LINESTYLE = 0, $
    POSITION = [(*(*info).plotpos).intx0,(*(*info).plotpos).inty0,$
                (*(*info).plotpos).intx1,(*(*info).plotpos).inty1], $
    XTICKLEN = (*(*info).plotaxes).intxticklen, YTICKLEN = (*(*info).plotaxes).intyticklen
	IF (condition[0] NE -1) THEN BEGIN
		selcol = (*(*(*info).intparams).selcol_diagnostics)[condition]
		FOR i=0,N_ELEMENTS(condition)-1 DO BEGIN
			IF (*(*info).dataswitch).spfile THEN BEGIN
				ssp = REFORM( ( ( *(*(*info).data).spdata)[ $
          FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + $
          FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + $
					(*(*info).dataparams).s ] )[condition[i],*] )
			ENDIF ELSE BEGIN
				ssp = FLTARR((*(*info).dataparams).nt)
				FOR t=0,(*(*info).dataparams).nt-1 DO BEGIN
					ssp[t] = ((*(*(*info).data).imagedata)[$
            t * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
            (*(*info).dataparams).s * (*(*info).dataparams).nlp + $
						(*(*info).dataparams).lp])[(*(*info).dataparams).x,(*(*info).dataparams).y]
				ENDFOR
			ENDELSE
			avgdint = ssp / ABS(MEAN(ssp, /NAN))
			LOADCT, 12, /SILENT
			plotcol = ((*(*info).intparams).colors_diagnostics)[selcol[i]]
			OPLOT, FINDGEN((*(*info).dataparams).nt)*(*(*info).plotaxes).dt, avgdint, COLOR=plotcol, $
        LINESTYLE = (*(*(*info).intparams).lines_diagnostics)[condition[i]]
			LOADCT, 0, /SILENT
		ENDFOR
	ENDIF
	t_range = (*(*info).plotaxes).int_upp_t - (*(*info).plotaxes).int_low_t + 1
  IF (*(*info).plotswitch).multichannel THEN $
	  XYOUTS,(t_range*(*(*info).plotaxes).dt)*0.1+(*(*info).plotaxes).int_low_t*(*(*info).plotaxes).dt, $
    (int_upp_y-int_low_y)*0.9+int_low_y, 'Stokes '+((*(*info).stokesparams).labels)[(*(*info).dataparams).s], $
    COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).dispparams).t GE (*(*info).plotaxes).int_low_t) AND $
      ((*(*info).dispparams).t LE (*(*info).plotaxes).int_upp_t)) THEN $
    PLOTS, [1,1] * (*(*info).dispparams).t * (*(*info).plotaxes).dt, [int_upp_y, int_low_y], $
      COLOR = (*(*info).plotparams).plotcol
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).intwid,int_low_y,int_upp_y], $
                          labels=['Window ID for draw','Lower y-value','Upper y-value']
END

PRO CRISPEX_DRAW_SPECTRAL_MAIN, event, LS_ONLY=ls_only, SP_ONLY=sp_only
; (Re)draw detailed spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Pre-draw procedures
  IF (~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showsp) THEN BEGIN
  	IF (*(*info).dispparams).slices_imscale THEN BEGIN
      CRISPEX_DRAW_SCALING,event,imdisp,minimum,maximum, /MAIN
      minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
        (*(*info).scaling).minimum[(*(*info).intparams).lp_diag_all],$
        (*(*info).scaling).maximum[(*(*info).intparams).lp_diag_all])
    ENDIF 
  ENDIF
  pass = 0L
  IF KEYWORD_SET(SP_ONLY) THEN $
    ns = 1 $
  ELSE $
    ns = TOTAL((*(*info).stokesparams).select_sp)
  ; Loop over all Stokes parameters for detailed spectrum plots
	FOR i=0,ns-1 DO BEGIN
    IF (~KEYWORD_SET(SP_ONLY) AND (*(*info).winswitch).showls) THEN BEGIN
   		s = (WHERE((*(*info).stokesparams).select_sp EQ 1))[i]
   		spec = ((*(*info).dataparams).spec)[*,s]
   		IF (*(*info).plotswitch).scalestokes THEN $
         ms = (*(*info).dataparams).ms $
       ELSE $
         ms = ((*(*info).dataparams).ms)[s]
   		ls_low_y = (*(*(*info).plotaxes).ls_low_y)[s] 
   		ls_upp_y = (*(*(*info).plotaxes).ls_upp_y)[s]
   		order_corr=0.
   		IF (*(*info).dispswitch).detspect_scale THEN $
         lsytitle = 'Scaled '+STRLOWCASE((*(*info).plottitles).lsytitle) $
      ELSE BEGIN
   			IF ((FLOOR(ALOG10(ABS(ls_low_y))) LE -2) OR $
          (FLOOR(ALOG10(ABS(ls_upp_y))) GE 3)) THEN BEGIN
   				order_corr = FLOOR(ALOG10(ABS(ls_upp_y)))
   				lsytitle = (*(*info).plottitles).lsytitle+' (x10!U'+STRTRIM(order_corr,2)+'!N)'
   			ENDIF ELSE lsytitle = (*(*info).plottitles).lsytitle
   			ls_low_y /= (10.^(order_corr))
   			ls_upp_y /= (10.^(order_corr))
   		ENDELSE
   		IF ((*(*info).dispswitch).detspect_scale EQ 0) THEN BEGIN
   			spec *= ms * ((*(*info).paramparams).scale_cubes)[0] / (10.^(order_corr))
   			ms = 1
   		ENDIF
      ; Set x-axes parameters depending on # of diagnostics
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        xticklen = 1E-9
        xtitle = ''
        xtickname = REPLICATE(' ',60)
      ENDIF ELSE BEGIN
        xticklen = (*(*info).plotaxes).lsxticklen
        xtitle = (*(*info).plottitles).spxtitle
        xtickname = ''
        yticklen = (*(*info).plotaxes).lsyticklen
      ENDELSE
      IF (*(*info).plotswitch).v_dop_set THEN topxtitle = 'Doppler velocity [km/s]'
      ; Get data for display
  		IF ((*(*info).dataswitch).spfile EQ 1) THEN BEGIN
        spidx = FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * $
                    (*(*info).dataparams).ns + $
  				      FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + s  
  			ssp = ( ( *(*(*info).data).spdata)[spidx] )[*,(*(*info).dispparams).t_main]/ms
  		ENDIF ELSE BEGIN    ; If no spectral cube supplied, determine from image cube
  			IF (*(*info).dataswitch).onecube THEN $
				  ssp = (*(*(*info).data).sspscan)[$
            FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y),$
            (s * (*(*info).dataparams).nlp):((s+1) * (*(*info).dataparams).nlp - 1)]/ms  $
        ELSE $
				  ssp = (*(*(*info).data).scan)[$
            FIX((*(*info).dataparams).x),FIX((*(*info).dataparams).y),$
            (s * (*(*info).dataparams).nlp):((s+1) * (*(*info).dataparams).nlp -1)]/ms
  		ENDELSE			
    ENDIF
    ; Determine proportional spectral window size for LS
    diag_range_ls = *(*(*info).plotaxes).diag_ratio * $
      (((*(*info).plotpos).lsx1)[i] - ((*(*info).plotpos).lsx0)[i])
    ; Determine current diagnostic window LP is in
    lp_diag = TOTAL((*(*info).dataparams).lp GE *(*(*info).intparams).diag_starts)-1
    FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
      IF (~KEYWORD_SET(SP_ONLY) AND (*(*info).winswitch).showls) THEN BEGIN
        WSET, (*(*info).winids).lswid
        ; Determine xrange to display
        IF ((*(*info).intparams).ndiagnostics GT 1) THEN $
          xrange = (*(*info).dataparams).lps[[$
            (*(*info).intparams).diag_start[(*(*(*info).intparams).wheredispdiag)[d]],$
            ((*(*info).intparams).diag_start[(*(*(*info).intparams).wheredispdiag)[d]]+$
            (*(*info).intparams).diag_width[(*(*(*info).intparams).wheredispdiag)[d]]-1)]] $
        ELSE $
          xrange = [(*(*info).dataparams).lps[(*(*info).dispparams).lp_low], $
            (*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]]
        ; Set y-axes parameters based on diagnostic plot
        IF (d EQ 0) THEN BEGIN
          offset = 0 
          ytickname = '' 
          ytitle = lsytitle
        ENDIF ELSE BEGIN
          offset = TOTAL(diag_range_ls[0:(d-1)])
          ytickname = REPLICATE(' ',60)
          ytitle = ''
        ENDELSE
        ; Determine lower left corner position of plot
        lsx0 = ((*(*info).plotpos).lsx0)[i] + offset
        lsx1 = diag_range_ls[d] + lsx0
        ; Plot basic window with average spectrum
     		PLOT, (*(*info).dataparams).lps, spec*(*(*info).scaling).mult_val[$
          (*(*(*info).intparams).wheredispdiag)[d]], $
          /NORM, CHARSIZE=1, YS=1, YR=[ls_low_y,ls_upp_y], $
               XR=xrange, YTICKNAME=ytickname, XTICKINTERVAL=(*(*info).plotaxes).xtickinterval, $
               XSTYLE = (*(*info).plotswitch).v_dop_set * 8 + 1, $
               BACKGROUND = (*(*info).plotparams).bgplotcol, $
               XTITLE = xtitle, YTITLE=ytitle, $
               POSITION = [lsx0,((*(*info).plotpos).lsy0)[i],lsx1,((*(*info).plotpos).lsy1)[i]], $
               XTICKLEN = (*(*info).plotaxes).lsxticklen, YTICKLEN = (*(*info).plotaxes).lsyticklen, $
               COLOR = (*(*info).plotparams).plotcol, LINE=3, $
               NOERASE=((((*(*info).intparams).ndiagnostics GT 1) AND (d GT 0)) OR (i GT 0))
        IF ((*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]] NE 1) THEN BEGIN
          XYOUTS,0.1*(xrange[1]-xrange[0])+xrange[0], 0.9*(ls_upp_y-ls_low_y)+ls_low_y, $
            STRING((*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]], FORMAT='(F'+$
            STRTRIM(FLOOR(ALOG10(ABS($
            (*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]])))+3+$
            ((*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]] LT 0),2)+'.1)')+'x',$
            COLOR=(*(*info).plotparams).plotcol, /DATA
        ENDIF
        ; In case of multiple diagnostics
        IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
          ; Draw x-title centered on plot box on first pass
          IF (d EQ 0) THEN $
            XYOUTS,((*(*info).plotpos).lsx1-(*(*info).plotpos).lsx0)/2.+(*(*info).plotpos).lsx0,$
              (*(*info).plotpos).lsy0/3., (*(*info).plottitles).spxtitle,ALIGNMENT=0.5, $
              COLOR = (*(*info).plotparams).plotcol,/NORMAL
          ; Set range for Doppler axis
          vdop_xrange = (*(*(*info).plotaxes).v_dop[(*(*(*info).intparams).wheredispdiag)[d]])[$
            [0,(*(*(*info).intparams).diag_widths)[d]-1]]
        ENDIF ELSE $
          vdop_xrange = (*(*(*info).plotaxes).v_dop[0])[$
            [(*(*info).dispparams).lp_low,(*(*info).dispparams).lp_upp]]
        ; Display Stokes label if Stokes data
        IF (*(*info).plotswitch).multichannel THEN $
    		  XYOUTS, ((*(*info).dataparams).lps[(*(*info).dispparams).lp_upp]-$
                   (*(*info).dataparams).lps[(*(*info).dispparams).lp_low])*0.1+$
                   (*(*info).dataparams).lps[(*(*info).dispparams).lp_low], $
    			        (ls_upp_y-ls_low_y)*0.9+ls_low_y, 'Stokes '+((*(*info).stokesparams).labels)[s], $
                  COLOR = (*(*info).plotparams).plotcol
        ; Display Doppler velocities on top axis if info available
    	  IF ((*(*info).plotswitch).v_dop_set EQ 1) THEN BEGIN
          IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
            XYOUTS,((*(*info).plotpos).lsx1-(*(*info).plotpos).lsx0)/2.+(*(*info).plotpos).lsx0,$
              (*(*info).plotpos).lsy0/5.*3+(*(*info).plotpos).lsy1, topxtitle, ALIGNMENT=0.5, $
              COLOR = (*(*info).plotparams).plotcol,/NORMAL
            topxtitle = ''
          ENDIF 
          ; Draw top axis
    	  	AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).lsxticklen, $
;              XRANGE = [(*(*(*info).plotaxes).v_dop[(*(*(*info).intparams).wheredispdiag)[d]])[0], $
;              (*(*(*info).plotaxes).v_dop[(*(*(*info).intparams).wheredispdiag)[d]])[$
;                (*(*(*info).intparams).diag_widths)[d]-1]], XSTYLE=1, $
              XRANGE=vdop_xrange, XSTYLE=1, $
    	  		XTITLE = topxtitle, COLOR = (*(*info).plotparams).plotcol,$
            XTICKINTERVAL=(*(*info).plotaxes).xdoptickinterval
        ENDIF 
        ; Overplot detailed spectrum
    		IF ((*(*info).dispswitch).detspect_scale EQ 0) THEN $
          ssp *= ((*(*info).paramparams).scale_cubes)[0] / (10.^(order_corr))
    		OPLOT, (*(*info).dataparams).lps, $
          ssp*(*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]], $
          LINE=0, COLOR = (*(*info).plotparams).plotcol
        ; Overplot average minus detailed spectrum
    		IF (*(*info).plotswitch).subtract THEN $
    			OPLOT, (*(*info).dataparams).lps, $
            (spec-ssp)*(*(*info).scaling).mult_val[(*(*(*info).intparams).wheredispdiag)[d]], $
            COLOR = (*(*info).plotparams).plotcol, LINE=2
        ; Draw line through y=0
    		IF ((ls_low_y LT 0.) AND (ls_upp_y GT 0.)) THEN $
          PLOTS, xrange, [0.,0.], COLOR = (*(*info).plotparams).plotcol
        IF (d EQ lp_diag) THEN BEGIN
          ; Overplot spectral indicator
    		  PLOTS, [1,1] * (*(*info).dataparams).lps[(*(*info).dataparams).lp],[ls_low_y,ls_upp_y], $
                COLOR = (*(*info).plotparams).plotcol
          ; Overplot 2nd spectral indicator in case of Doppler image display
      		IF (*(*info).dispswitch).drawdop THEN $
            PLOTS, [1,1] * (*(*info).dataparams).lps[(*(*info).dataparams).lp_dop],[ls_low_y,ls_upp_y], $
                    COLOR = (*(*info).plotparams).plotcol
        ENDIF
        ; Overplot reference spectral indicator if same number of spectral positions
    		IF ((*(*info).winswitch).showref AND ((*(*info).ctrlsswitch).lp_ref_lock EQ 0) AND $
           ((*(*info).dataswitch).refspfile EQ 0) AND ((*(*info).dataparams).refnlp GT 1)) THEN BEGIN
    			PLOTS, [1,1] * (*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref],$
                [ls_low_y,ls_upp_y], COLOR = (*(*info).plotparams).plotcol
    		ENDIF
      ENDIF
      IF ((pass EQ 0) AND ~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showsp) THEN BEGIN
        WSET, (*(*info).winids).spwid
        ; Get spectrum-time diagram by diagnostic
        IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
          lp_lower = TOTAL((*(*(*info).intparams).diag_widths)[0:d]) - $
            (*(*(*info).intparams).diag_widths)[d]
          lp_upper = lp_lower+(*(*(*info).intparams).diag_widths)[d]-1
        ENDIF ELSE BEGIN
          lp_lower = 0
          lp_upper = (*(*info).dispparams).lp_upp-(*(*info).dispparams).lp_low
        ENDELSE
        tmp_disp = (*(*(*info).data).spslice)[lp_lower:lp_upper,*]
        ; Display spectrum-time diagram by diagnostic
        TV,(CONGRID( BYTSCL(tmp_disp, $
          MIN=(*(*info).scaling).spslice_min[(*(*(*info).intparams).wheredispdiag)[d]],$ 
          MAX=(*(*info).scaling).spslice_max[(*(*(*info).intparams).wheredispdiag)[d]],$
          /NAN), (*(*info).dispparams).nlpreb*(*(*(*info).plotaxes).diag_ratio)[d], $
          (*(*info).dispparams).ntreb, INTERP = (*(*info).dispparams).interpspslice, $
          /CENTER) ), (d GE 1)*TOTAL((*(*(*info).plotaxes).diag_range_sp)[0:(d-1)])+$
          (*(*info).plotpos).spx0, (*(*info).plotpos).spy0, /NORM
      ENDIF
    ENDFOR
    IF (~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showsp) THEN BEGIN
      !X.WINDOW = [(*(*info).plotpos).spx0,(*(*info).plotpos).spx1]
      !Y.WINDOW = [(*(*info).plotpos).spy0,(*(*info).plotpos).spy1]
      ; In case of multiple diagnostics, overplot separator axes
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        ; Loop over all but the first diagnostic for plotting separators
        FOR d=1,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
          AXIS,(d GE 1)*TOTAL((*(*(*info).plotaxes).diag_range_sp)[0:(d-1)])+$
            (*(*info).plotpos).spx0, YAXIS=0, $
            YTICKLEN=(*(*info).plotaxes).spyticklen, YTICKNAME = REPLICATE(' ',60), $
            COLOR=100,/NORMAL, $
            YRANGE=[(*(*info).dispparams).t_low_main, $
            (*(*info).dispparams).t_upp_main], YSTYLE=1, /NOERASE
          AXIS,(d GE 1)*TOTAL((*(*(*info).plotaxes).diag_range_sp)[0:(d-1)])+$
            (*(*info).plotpos).spx0, YAXIS=1, $
            YTICKLEN=(*(*info).plotaxes).spyticklen, YTICKNAME=REPLICATE(' ',60), $
            COLOR=100,/NORMAL,$
            YRANGE = [(*(*info).dispparams).t_low_main, (*(*info).dispparams).t_upp_main], $
            YSTYLE=1, /NOERASE
        ENDFOR
      ENDIF
      ; Overplot time indicator
    	PLOTS, [(*(*info).plotpos).spx0,(*(*info).plotpos).spx1], $
              [1,1]*( ((*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t]-$
              (*(*info).dispparams).t_low_main) / $
              FLOAT((*(*info).dispparams).t_upp_main-(*(*info).dispparams).t_low_main) * $
              (*(*info).plotpos).yplspw + (*(*info).plotpos).spy0), /NORMAL, COLOR = 100
      ; Overplot lp indicator
;      lp_diag = TOTAL((*(*info).dataparams).lp GE *(*(*info).intparams).diag_starts)-1
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        ; Determine in which diagnostic the current LP lies
        IF (lp_diag EQ 0) THEN $
          offset = 0 $
        ELSE $
          offset = TOTAL((*(*(*info).plotaxes).diag_range_sp)[0:(lp_diag-1)])
        spx0 = (*(*info).plotpos).spx0 + offset
        xplspw = (*(*(*info).plotaxes).diag_range_sp)[lp_diag]
        ; Get the corresponding lp_disprange
        lp_disprange = (*(*info).dataparams).lps[$
                        (*(*(*info).intparams).diag_starts)[lp_diag]+$
                        (*(*(*info).intparams).diag_widths)[lp_diag]-1] - $
                       (*(*info).dataparams).lps[$
                        (*(*(*info).intparams).diag_starts)[lp_diag]]
        lp_lower = (*(*info).dataparams).lps[(*(*(*info).intparams).diag_starts)[lp_diag]]
      ENDIF ELSE BEGIN
        lp_disprange = ((*(*info).dataparams).lps[(*(*info).dispparams).lp_upp] - $
                        (*(*info).dataparams).lps[(*(*info).dispparams).lp_low])
        spx0 = (*(*info).plotpos).spx0
        xplspw = (*(*info).plotpos).xplspw 
        lp_lower = (*(*info).dataparams).lps[(*(*info).dispparams).lp_low]
      ENDELSE
    	PLOTS, [1,1] * ( ((*(*info).dataparams).lps[(*(*info).dataparams).lp] - $
                        lp_lower) / lp_disprange * xplspw + spx0 ), $
    		[(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
      ; If drawing Doppler, overplot lp_dop indicator
    	IF (*(*info).dispswitch).drawdop THEN $         
        PLOTS, [1,1] * ( ((*(*info).dataparams).lps[(*(*info).dataparams).lp_dop] - $
                        (*(*info).dataparams).lps[$
                          (*(*(*info).intparams).diag_starts)[lp_diag]]) / $
                        lp_disprange * xplspw + spx0 ), $
    		  [(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
      ; If drawing reference, and refnlp=nlp but lp_ref != lp, then overplot lp_ref indicator
    	IF ((*(*info).winswitch).showref AND ((*(*info).ctrlsswitch).lp_ref_lock EQ 0) AND $
         ((*(*info).dataswitch).refspfile EQ 0) AND ((*(*info).dataparams).refnlp GT 1)) THEN $
    		PLOTS, [1,1] * ( ((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref] - $
                          (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low]) / $
                         ((*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_upp] - $
                          (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low]) * $
          (*(*info).plotpos).xplspw + (*(*info).plotpos).spx0 ), $
    			[(*(*info).plotpos).spy0, (*(*info).plotpos).spy1], /NORMAL, COLOR = 100
    ENDIF
    pass += 1L
  ENDFOR
END

PRO CRISPEX_DRAW_SPECTRAL_REF, event, LS_ONLY=ls_only, SP_ONLY=sp_only
; (Re)draw detailed spectrum procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Pre-draw procedures: spectrum-time diagram (SP)
  IF (~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showrefsp) THEN BEGIN
  	IF (*(*info).dispparams).slices_imscale THEN BEGIN
      sel_idx = (*(*info).intparams).lp_ref_diag_all+(*(*info).intparams).ndiagnostics
      CRISPEX_DRAW_SCALING,event,refdisp,minimum,maximum, /REF
      refminmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
        (*(*info).scaling).minimum[sel_idx], (*(*info).scaling).maximum[sel_idx])
    ENDIF 
  ENDIF
  ; Pre-draw procedures: detailed spectrum (LS)
	refspec = ((*(*info).dataparams).refspec)
	refms = (*(*info).dataparams).refms
	ls_low_y = (*(*info).plotaxes).ls_low_y_ref
	ls_upp_y = (*(*info).plotaxes).ls_upp_y_ref
	order_corr=0.
	IF (*(*info).dispswitch).ref_detspect_scale THEN $
    reflsytitle = 'Scaled '+STRLOWCASE((*(*info).plottitles).reflsytitle) $
  ELSE BEGIN
		IF ((FLOOR(ALOG10(ABS(ls_low_y))) LE -2) OR (FLOOR(ALOG10(ABS(ls_upp_y))) GE 3)) THEN BEGIN
			order_corr = FLOOR(ALOG10(ABS(ls_upp_y)))
			reflsytitle = (*(*info).plottitles).reflsytitle+' (x10!U'+STRTRIM(order_corr,2)+'!N)'
		ENDIF ELSE reflsytitle = (*(*info).plottitles).reflsytitle
		ls_low_y /= (10.^(order_corr))
		ls_upp_y /= (10.^(order_corr))
	ENDELSE
	IF ((*(*info).dispswitch).ref_detspect_scale EQ 0) THEN BEGIN
		refspec *= refms * ((*(*info).paramparams).scale_cubes)[1] / (10.^(order_corr))
		refms = 1
	ENDIF
  ; Set x-axes parameters depending on # of diagnostics
  IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
    xticklen = 1E-9
    xtitle = ''
    xtickname = REPLICATE(' ',60)
  ENDIF ELSE BEGIN
    xticklen = (*(*info).plotaxes).reflsxticklen
    xtitle = (*(*info).plottitles).refspxtitle
    xtickname = ''
    yticklen = (*(*info).plotaxes).reflsyticklen
  ENDELSE
  IF (*(*info).plotswitch).v_dop_set_ref THEN topxtitle = 'Doppler velocity [km/s]'
  ; Get data for display
  IF ((*(*info).dataswitch).refspfile EQ 1) THEN BEGIN
    sspidx = FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx + FIX((*(*info).dataparams).x)
    refssp = ( ( *(*(*info).data).refspdata)[sspidx] ) [*,(*(*info).dispparams).t_ref]/refms 
  ENDIF ELSE $
    refssp = (*(*(*info).data).refsspscan)[FIX((*(*info).dataparams).x),$
      FIX((*(*info).dataparams).y),*]/refms
  pass = 0L
  ; Determine proportional spectral window size for LS
  diag_range_ls = *(*(*info).plotaxes).refdiag_ratio * $
    ((*(*info).plotpos).reflsx1 - (*(*info).plotpos).reflsx0)
  lp_ref_diag = TOTAL((*(*info).dataparams).lp_ref GE *(*(*info).intparams).refdiag_starts)-1
  FOR d=0,(*(*info).intparams).ndisp_refdiagnostics-1 DO BEGIN
    mult_idx = (*(*info).intparams).ndiagnostics+(*(*(*info).intparams).wheredisprefdiag)[d]
    ; Draw REFLS window
    IF (~KEYWORD_SET(SP_ONLY) AND (*(*info).winswitch).showrefls) THEN BEGIN
	    WSET, (*(*info).winids).reflswid
      ; Determine xrange to display
      IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN $
        xrange = (*(*info).dataparams).reflps[[$
          (*(*info).intparams).refdiag_start[(*(*(*info).intparams).wheredisprefdiag)[d]],$
          ((*(*info).intparams).refdiag_start[(*(*(*info).intparams).wheredisprefdiag)[d]]+$
          (*(*info).intparams).refdiag_width[(*(*(*info).intparams).wheredisprefdiag)[d]]-1)]] $
      ELSE $
        xrange = [(*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low], $
          (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_upp]]
      ; Set y-axes parameters based on diagnostic plot
      IF (d EQ 0) THEN BEGIN
        offset = 0 
        ytickname = '' 
        ytitle = reflsytitle
      ENDIF ELSE BEGIN
        offset = TOTAL(diag_range_ls[0:(d-1)])
        ytickname = REPLICATE(' ',60)
        ytitle = ''
      ENDELSE
      ; Determine lower left corner position of plot
      reflsx0 = (*(*info).plotpos).reflsx0 + offset
      reflsx1 = diag_range_ls[d] + reflsx0
      ; Plot basic window with average spectrum
     	PLOT, (*(*info).dataparams).reflps, refspec*(*(*info).scaling).mult_val[mult_idx], $
        /NORM, CHARSIZE=1, YS=1, YR=[ls_low_y,ls_upp_y], $
             XR=xrange, YTICKNAME=ytickname, XTICKINTERVAL=(*(*info).plotaxes).xreftickinterval,$
             XSTYLE = (*(*info).plotswitch).v_dop_set_ref * 8 + 1, $
             XTICK_GET=xtickvals, BACKGROUND = (*(*info).plotparams).bgplotcol, $
             XTITLE = xtitle, YTITLE=ytitle, $
             POSITION = [reflsx0,((*(*info).plotpos).reflsy0),reflsx1,((*(*info).plotpos).reflsy1)], $
             XTICKLEN = (*(*info).plotaxes).reflsxticklen, YTICKLEN = (*(*info).plotaxes).reflsyticklen, $
             COLOR = (*(*info).plotparams).plotcol, LINE=3, $
             NOERASE=(((*(*info).intparams).nrefdiagnostics GT 1) AND (d GT 0))
        IF ((*(*info).scaling).mult_val[mult_idx] NE 1) THEN BEGIN
          XYOUTS,0.1*(xrange[1]-xrange[0])+xrange[0], 0.9*(ls_upp_y-ls_low_y)+ls_low_y, $
            STRING((*(*info).scaling).mult_val[mult_idx], $
            FORMAT='(F'+STRTRIM(FLOOR(ALOG10(ABS((*(*info).scaling).mult_val[mult_idx])))+3+$
            ((*(*info).scaling).mult_val[mult_idx] LT 0),2)+'.1)')+'x',$
            COLOR=(*(*info).plotparams).plotcol, /DATA
        ENDIF
      ; In case of multiple diagnostics
      IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
        ; Draw x-title centered on plot box on first pass
        IF (d EQ 0) THEN $
          XYOUTS,((*(*info).plotpos).reflsx1-(*(*info).plotpos).reflsx0)/2.+(*(*info).plotpos).reflsx0,$
            (*(*info).plotpos).reflsy0/3.,(*(*info).plottitles).refspxtitle,ALIGNMENT=0.5, $
            COLOR=(*(*info).plotparams).plotcol, /NORMAL
      ENDIF
      ; Display Doppler velocities on top axis if info available
    	IF ((*(*info).plotswitch).v_dop_set_ref EQ 1) THEN BEGIN
        IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
          XYOUTS,((*(*info).plotpos).reflsx1-(*(*info).plotpos).reflsx0)/2.+(*(*info).plotpos).reflsx0,$
            (*(*info).plotpos).reflsy0/5.*3+(*(*info).plotpos).reflsy1, topxtitle, ALIGNMENT=0.5, $
            COLOR = (*(*info).plotparams).plotcol,/NORMAL
          topxtitle = ''
        ENDIF 
        ; Draw top axis
    		AXIS, XAXIS=1, XTICKLEN = (*(*info).plotaxes).reflsxticklen, $
              XRANGE=[(*(*(*info).plotaxes).v_dop_ref[(*(*(*info).intparams).wheredisprefdiag)[d]])[0], $
              (*(*(*info).plotaxes).v_dop_ref[(*(*(*info).intparams).wheredisprefdiag)[d]])[$
                (*(*(*info).intparams).refdiag_widths)[d]-1]], XSTYLE=1, $
    			XTITLE = topxtitle, COLOR = (*(*info).plotparams).plotcol,$
          XTICKINTERVAL=(*(*info).plotaxes).xrefdoptickinterval
      ENDIF 
      ; Overplot detailed spectrum
  	  IF ((*(*info).dispswitch).ref_detspect_scale EQ 0) THEN $
        refssp *= ((*(*info).paramparams).scale_cubes)[1] / (10.^(order_corr))
  	  OPLOT, (*(*info).dataparams).reflps, refssp*(*(*info).scaling).mult_val[mult_idx], $
        LINE=0, COLOR = (*(*info).plotparams).plotcol
      ; Overplot average minus detailed spectrum
  	  IF (*(*info).plotswitch).ref_subtract THEN $
  	  	OPLOT, (*(*info).dataparams).reflps, (refspec-refssp)*(*(*info).scaling).mult_val[mult_idx], $
          COLOR=(*(*info).plotparams).plotcol, LINE=2
      ; Draw line through y=0
  	  IF ((ls_low_y LT 0.) AND (ls_upp_y GT 0.)) THEN $
        PLOTS, xrange, [0.,0.], COLOR = (*(*info).plotparams).plotcol
      IF (d EQ lp_ref_diag) THEN $
      ; Overplot spectral indicator
  	  PLOTS, [1,1] * (*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref],[ls_low_y,ls_upp_y],$
        COLOR = (*(*info).plotparams).plotcol
    ENDIF
    ; Draw REFSP window
    IF ((pass EQ 0) AND ~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showrefsp) THEN BEGIN
      WSET, (*(*info).winids).refspwid
      ; Get spectrum-time diagram by diagnostic
      IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
        lp_lower = TOTAL((*(*(*info).intparams).refdiag_widths)[0:d]) - $
          (*(*(*info).intparams).refdiag_widths)[d]
        lp_upper = lp_lower+(*(*(*info).intparams).refdiag_widths)[d]-1
      ENDIF ELSE BEGIN
        lp_lower = 0
        lp_upper = (*(*info).dispparams).lp_ref_upp-(*(*info).dispparams).lp_ref_low
      ENDELSE
      tmp_disp = (*(*(*info).data).refspslice)[lp_lower:lp_upper,*]
      ; Display spectrum-time diagram by diagnostic
      TV,(CONGRID( BYTSCL(tmp_disp, $
        MIN=(*(*info).scaling).spslice_min[mult_idx], $
        MAX=(*(*info).scaling).spslice_max[mult_idx], $
        /NAN), (*(*info).dispparams).refnlpreb*(*(*(*info).plotaxes).refdiag_ratio)[d], $
        (*(*info).dispparams).refntreb, INTERP = (*(*info).dispparams).interpspslice, $
        /CENTER) ), (d GE 1)*TOTAL((*(*(*info).plotaxes).refdiag_range_sp)[0:(d-1)])+$
        (*(*info).plotpos).refspx0, (*(*info).plotpos).refspy0, /NORM
    ENDIF
  ENDFOR
  IF (~KEYWORD_SET(LS_ONLY) AND (*(*info).winswitch).showrefsp) THEN BEGIN
    !X.WINDOW = [(*(*info).plotpos).refspx0,(*(*info).plotpos).refspx1]
    !Y.WINDOW = [(*(*info).plotpos).refspy0,(*(*info).plotpos).refspy1]
    ; In case of multiple diagnostics, overplot separator axes
    IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
      ; Loop over all but the last diagnostic for plotting separators
      FOR d=1,(*(*info).intparams).ndisp_refdiagnostics-1 DO BEGIN
        AXIS,(d GE 1)*TOTAL((*(*(*info).plotaxes).refdiag_range_sp)[0:(d-1)])+$
          (*(*info).plotpos).refspx0, YAXIS=0, $
          YTICKLEN=(*(*info).plotaxes).refspyticklen, YTICKNAME = REPLICATE(' ',60), $
          COLOR=100,/NORMAL, $
          YRANGE=[(*(*info).dispparams).t_low_ref, $
          (*(*info).dispparams).t_upp_ref], YSTYLE=1, /NOERASE
        AXIS,(d GE 1)*TOTAL((*(*(*info).plotaxes).refdiag_range_sp)[0:(d-1)])+$
          (*(*info).plotpos).refspx0, YAXIS=1, $
          YTICKLEN=(*(*info).plotaxes).refspyticklen, YTICKNAME=REPLICATE(' ',60), $
          COLOR=100,/NORMAL,$
          YRANGE = [(*(*info).dispparams).t_low_ref, (*(*info).dispparams).t_upp_ref], $
          YSTYLE=1, /NOERASE
      ENDFOR
    ENDIF
    ; Overplot time indicator
  	PLOTS, [(*(*info).plotpos).refspx0,(*(*info).plotpos).refspx1], $
            [1,1]*( ((*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t]-$
            (*(*info).dispparams).t_low_ref) / $
            FLOAT((*(*info).dispparams).t_upp_ref-(*(*info).dispparams).t_low_ref) * $
            (*(*info).plotpos).refyplspw + (*(*info).plotpos).refspy0), /NORMAL, COLOR = 100
    ; Overplot lp indicator
    IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN BEGIN
      ; Determine in which diagnostic the current LP lies
      IF (lp_ref_diag EQ 0) THEN $
        offset = 0 $
      ELSE $
        offset = TOTAL((*(*(*info).plotaxes).refdiag_range_sp)[0:(lp_ref_diag-1)])
      refspx0 = (*(*info).plotpos).refspx0 + offset
      refxplspw = (*(*(*info).plotaxes).refdiag_range_sp)[lp_ref_diag]
      ; Get the corresponding lp_disprange
      lp_ref_disprange = (*(*info).dataparams).reflps[$
                      (*(*(*info).intparams).refdiag_starts)[lp_ref_diag]+$
                      (*(*(*info).intparams).refdiag_widths)[lp_ref_diag]-1] - $
                     (*(*info).dataparams).reflps[$
                      (*(*(*info).intparams).refdiag_starts)[lp_ref_diag]]
    ENDIF ELSE BEGIN
      lp_ref_disprange = ((*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_upp] - $
                      (*(*info).dataparams).reflps[(*(*info).dispparams).lp_ref_low])
      refspx0 = (*(*info).plotpos).refspx0
      refxplspw = (*(*info).plotpos).refxplspw 
    ENDELSE
  	PLOTS, [1,1] * ( ((*(*info).dataparams).reflps[(*(*info).dataparams).lp_ref] - $
                      (*(*info).dataparams).reflps[$
                        (*(*(*info).intparams).refdiag_starts)[lp_ref_diag]]) / $
                      lp_ref_disprange * refxplspw + refspx0 ), $
  		[(*(*info).plotpos).refspy0, (*(*info).plotpos).refspy1], /NORMAL, COLOR = 100
  ENDIF
END

PRO CRISPEX_DRAW_PHIS, event
; (Re)draw spectral slice along a slit procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).phiswid
  ; If 2D, then display
  IF ((SIZE(*(*(*info).data).phislice))[0] EQ 2) THEN BEGIN
    FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
      ; Get spectrum-slit diagram by diagnostics
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
        lp_lower = TOTAL((*(*(*info).intparams).diag_widths)[0:d]) - $
          (*(*(*info).intparams).diag_widths)[d]
        lp_upper = lp_lower+(*(*(*info).intparams).diag_widths)[d]-1
      ENDIF ELSE BEGIN
        lp_lower = 0
        lp_upper = (*(*info).dispparams).lp_upp-(*(*info).dispparams).lp_low
      ENDELSE
      tmp_disp = (*(*(*info).data).phislice)[lp_lower:lp_upper,*]
      ; Display spectrum-slit diagram by diagnostics
      TV,(CONGRID( BYTSCL(tmp_disp, $
        MIN=(*(*info).scaling).phislice_min[(*(*(*info).intparams).wheredispdiag)[d]],$ 
        MAX=(*(*info).scaling).phislice_max[(*(*(*info).intparams).wheredispdiag)[d]],$
        /NAN), (*(*info).dispparams).phisnlpreb*(*(*(*info).plotaxes).diag_ratio)[d], $
        (*(*info).dispparams).nphireb, INTERP = (*(*info).dispparams).interpspslice, /CENTER) ), $
        (d GE 1)*TOTAL((*(*(*info).plotaxes).diag_range_phis)[0:(d-1)])+$
        (*(*info).plotpos).phisx0, (*(*info).plotpos).phisy0, /NORMAL
    ENDFOR
    !X.WINDOW = [(*(*info).plotpos).phisx0,(*(*info).plotpos).phisx1]
    !Y.WINDOW = [(*(*info).plotpos).phisy0,(*(*info).plotpos).phisy1]
    IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
      ; Loop over all but the last diagnostic for plotting separators
      FOR d=1,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
        AXIS,TOTAL((*(*(*info).plotaxes).diag_range_phis)[0:(d-1)])+(*(*info).plotpos).phisx0, YAXIS=0, $
          YTICKLEN=(*(*info).plotaxes).phisyticklen, YTICKNAME = REPLICATE(' ',60), COLOR=100,/NORMAL,$
          YRANGE = (*(*info).plotaxes).phis_yrange, /YSTYLE, /NOERASE
        AXIS,TOTAL((*(*(*info).plotaxes).diag_range_phis)[0:(d-1)])+(*(*info).plotpos).phisx0, YAXIS=1, $
          YTICKLEN=(*(*info).plotaxes).phisyticklen, YTICKNAME = REPLICATE(' ',60), COLOR=100,/NORMAL,$
          YRANGE = (*(*info).plotaxes).phis_yrange, /YSTYLE, /NOERASE
      ENDFOR
    ENDIF
    ; Overplot slit center indicator
  	PLOTS, [(*(*info).plotpos).phisx0, (*(*info).plotpos).phisx1], $
      [1,1] * ( (((*(*info).phiparams).nw_cur - (*(*info).phiparams).nphi)/2. + $
                  (*(*info).phiparams).sphi+0.5)/FLOAT((*(*info).phiparams).nw_cur)*$
                  (*(*info).plotpos).phisyplspw + (*(*info).plotpos).phisy0), /NORMAL, COLOR = 100
    ; Overplot lp indicator
    lp_diag = TOTAL((*(*info).dataparams).lp GE *(*(*info).intparams).diag_starts)-1
    IF ((*(*info).intparams).ndiagnostics GT 1) THEN BEGIN
      ; Determine in which diagnostic the current LP lies
      IF (lp_diag EQ 0) THEN $
        offset = 0 $
      ELSE $
        offset = TOTAL((*(*(*info).plotaxes).diag_range_phis)[0:(lp_diag-1)])
      phisx0 = (*(*info).plotpos).phisx0 + offset
      phisxplspw = (*(*(*info).plotaxes).diag_range_phis)[lp_diag]
      ; Get the corresponding lp_disprange
      lp_disprange = (*(*info).dataparams).lps[$
                      (*(*(*info).intparams).diag_starts)[lp_diag]+$
                      (*(*(*info).intparams).diag_widths)[lp_diag]-1] - $
                      (*(*info).dataparams).lps[$
                      (*(*(*info).intparams).diag_starts)[lp_diag]]
      lp_lower = (*(*info).dataparams).lps[(*(*(*info).intparams).diag_starts)[lp_diag]]
    ENDIF ELSE BEGIN
      lp_disprange = ((*(*info).dataparams).lps[(*(*info).dispparams).lp_upp] - $
                      (*(*info).dataparams).lps[(*(*info).dispparams).lp_low])
      phisx0 = (*(*info).plotpos).phisx0
      phisxplspw = (*(*info).plotpos).phisxplspw 
      lp_lower = (*(*info).dataparams).lps[(*(*info).dispparams).lp_low]
    ENDELSE
  	PLOTS, [1,1] * ( ((*(*info).dataparams).lps[(*(*info).dataparams).lp] - $
                      lp_lower) / lp_disprange * phisxplspw + phisx0 ), $
  		[(*(*info).plotpos).phisy0, (*(*info).plotpos).phisy1], /NORMAL, COLOR = 100
  ; If not 2D, then show black screen and "error" message
  ENDIF ELSE BEGIN
    TV, CONGRID(*(*(*info).data).emptydopslice,(*(*info).dispparams).phisnlpreb, $
      (*(*info).dispparams).nphireb), (*(*info).plotpos).phisx0, (*(*info).plotpos).phisy0, /NORMAL
		XYOUTS, (*(*info).plotpos).phisxplspw/2.+(*(*info).plotpos).phisx0,$
            (*(*info).plotpos).phisyplspw/2.+(*(*info).plotpos).phisy0,$
            'No Phi-slice availabe for !C selected pixel position', COLOR = 255, ALIGNMENT = 0.5, $
            CHARSIZE=1.2, /NORMAL
  ENDELSE
END

PRO CRISPEX_DRAW_LOOPSLAB, event
; (Re)draw loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((SIZE(*(*(*info).loopsdata).crossloc))[1] NE (*(*info).loopparams).np) THEN BEGIN
		CRISPEX_LOOP_GET, event
		CRISPEX_UPDATE_LP, event
		IF ((SIZE(*(*(*info).loopsdata).crossloc))[1] NE (*(*info).loopparams).np) THEN RETURN
	ENDIF
	WSET, (*(*info).winids).loopwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).loopwid], labels=['Window ID for draw']
	dispslice = (*(*(*info).loopsdata).loopslice)[*,$
    (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_low]:$
    (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_upp]]
	IF (*(*info).dispparams).slices_imscale THEN BEGIN
    CRISPEX_DRAW_SCALING,event,imdisp,minimum,maximum,/MAIN 
    minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
      (*(*info).scaling).minimum[(*(*info).intparams).lp_diag_all],$
      (*(*info).scaling).maximum[(*(*info).intparams).lp_diag_all])
  ENDIF ELSE BEGIN
;    minimum = MIN(dispslice,MAX=maximum)
    minmax = CRISPEX_SCALING_SLICES(dispslice, $
      (*(*info).scaling).gamma[(*(*info).intparams).lp_diag_all], $
      (*(*info).scaling).histo_opt_val[(*(*info).intparams).lp_diag_all], $
      (*(*info).scaling).minimum[(*(*info).intparams).lp_diag_all],$
      (*(*info).scaling).maximum[(*(*info).intparams).lp_diag_all])
  ENDELSE
	TV, CONGRID( BYTSCL(dispslice, MIN=minmax[0], MAX=minmax[1], /NAN), (*(*info).dispparams).loopnlxreb, $
    (*(*info).dispparams).loopntreb, /INTERP), $
		(*(*info).plotpos).loopx0, (*(*info).plotpos).loopy0,/NORM
  ; Overplot time indicator
;	PLOTS, [(*(*info).plotpos).loopx0, (*(*info).plotpos).loopx1], $
;    [1,1] * ( ((*(*info).dispparams).t-(*(*info).dispparams).t_low) / $
;    FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).loopyplspw + $
	PLOTS, [(*(*info).plotpos).loopx0,(*(*info).plotpos).loopx1], $
          [1,1]*( ((*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t]-$
          (*(*info).dispparams).t_low_main) / $
          FLOAT((*(*info).dispparams).t_upp_main-(*(*info).dispparams).t_low_main) * $
	        (*(*info).plotpos).loopyplspw + (*(*info).plotpos).loopy0), /NORMAL, COLOR = 100
	FOR i=0,(*(*info).loopparams).np-1 DO BEGIN
		PLOTS, [1,1] * ( FLOAT((*(*(*info).loopsdata).crossloc)[i]) / $
      FLOAT((*(*info).loopsdata).loopsize-1) * (*(*info).plotpos).loopxplspw + $
      (*(*info).plotpos).loopx0 ), $
			[(*(*info).plotpos).loopy0, (*(*info).plotpos).loopy1], /NORMAL, COLOR = 100
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, (*(*(*info).loopsdata).crossloc)[i], $
        labels='crossloc['+STRTRIM(i,2)+']'
	ENDFOR
END

PRO CRISPEX_DRAW_REFLOOPSLAB, event
; (Re)draw reference loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).refloopwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).refloopwid], labels=['Window ID for draw']
	dispslice = (*(*(*info).loopsdata).refloopslice)[*,$
    (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).t_low]:$
    (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).t_upp]]
  sel_idx = (*(*info).intparams).lp_ref_diag_all+(*(*info).intparams).ndiagnostics
	IF (*(*info).dispparams).slices_imscale THEN BEGIN
    CRISPEX_DRAW_SCALING,event,refdisp,minimum,maximum,/REFERENCE 
    minmax = CRISPEX_SCALING_CONTRAST(minimum,maximum,$
      (*(*info).scaling).minimum[sel_idx], (*(*info).scaling).maximum[sel_idx])
  ENDIF ELSE BEGIN
    minmax = CRISPEX_SCALING_SLICES(dispslice, $
      (*(*info).scaling).gamma[sel_idx],(*(*info).scaling).histo_opt_val[sel_idx], $
      (*(*info).scaling).minimum[sel_idx],(*(*info).scaling).maximum[sel_idx])
  ENDELSE
	TV, CONGRID( BYTSCL(dispslice, MIN=minmax[0], MAX=minmax[1], /NAN), (*(*info).dispparams).refloopnlxreb, $
    (*(*info).dispparams).refloopntreb, /INTERP), $
		(*(*info).plotpos).refloopx0, (*(*info).plotpos).refloopy0,/NORM
  ; Overplot time indicator
;	PLOTS, [(*(*info).plotpos).refloopx0, (*(*info).plotpos).refloopx1], $
;    [1,1] * ( ((*(*info).dispparams).t-(*(*info).dispparams).t_low) / $
;    FLOAT((*(*info).dispparams).t_range-1) * (*(*info).plotpos).refloopyplspw + $
;		(*(*info).plotpos).refloopy0), /NORMAL, COLOR = 100
	PLOTS, [(*(*info).plotpos).refloopx0,(*(*info).plotpos).refloopx1], $
          [1,1]*( ((*(*(*info).dispparams).tarr_ref)[(*(*info).dispparams).t]-$
          (*(*info).dispparams).t_low_ref) / $
          FLOAT((*(*info).dispparams).t_upp_ref-(*(*info).dispparams).t_low_ref) * $
	        (*(*info).plotpos).refloopyplspw + (*(*info).plotpos).refloopy0), /NORMAL, COLOR = 100
	FOR i=0,(*(*info).loopparams).np-1 DO BEGIN
		PLOTS, [1,1] * ( FLOAT((*(*(*info).loopsdata).crossloc)[i]) / $
      FLOAT((*(*info).loopsdata).loopsize-1) * (*(*info).plotpos).refloopxplspw + $
      (*(*info).plotpos).refloopx0 ), $
			[(*(*info).plotpos).refloopy0, (*(*info).plotpos).refloopy1], /NORMAL, COLOR = 100
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, (*(*(*info).loopsdata).crossloc)[i], $
        labels='crossloc['+STRTRIM(i,2)+']'
	ENDFOR
END

PRO CRISPEX_DRAW_REST_LOOP, event
; (Re)draw restored loop timeslice procedure
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	FOR k=0,N_ELEMENTS(*(*(*info).winids).restlooptlb)-1 DO BEGIN
		WSET, (*(*(*info).winids).restloopwid)[k]
		IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
      CRISPEX_VERBOSE_GET, event, [(*(*(*info).winids).restloopwid)[k]], labels=['Window ID for draw']
    IF (*(*(*info).restoreparams).disp_imref)[k] THEN BEGIN
      tsel = *(*(*info).dispparams).tsel_ref 
      tarr = *(*(*info).dispparams).tarr_ref 
    ENDIF ELSE BEGIN
      tsel = *(*(*info).dispparams).tsel_main
      tarr = *(*(*info).dispparams).tarr_main
    ENDELSE
		IF (*(*(*info).dispswitch).restricted_t_range)[k] THEN BEGIN
			lower_t = 0
			upper_t = tsel[(*(*info).dispparams).t_upp] - tsel[(*(*info).dispparams).t_low]
		ENDIF ELSE BEGIN
			lower_t = tsel[(*(*info).dispparams).t_low]
			upper_t = tsel[(*(*info).dispparams).t_upp]
		ENDELSE
    lower_t_val = tarr[(*(*info).dispparams).t_low]
    upper_t_val = tarr[(*(*info).dispparams).t_upp]
    ; Get dispslice
		dispslice = (*(*(*(*info).loopsdata).rest_loopslice[k]))[*,lower_t:upper_t]
		IF (*(*(*info).restoreparams).disp_imref)[k] THEN BEGIN
			IF (*(*info).dispparams).slices_imscale THEN $
        CRISPEX_DRAW_SCALING,event,refdisp,minimum,maximum, /REFERENCE $
      ELSE $
        minimum = MIN(dispslice,MAX=maximum, /NAN)
		ENDIF ELSE BEGIN
			IF (*(*info).dispparams).slices_imscale THEN $
        CRISPEX_DRAW_SCALING,event,imdisp,minimum,maximum,/MAIN $
      ELSE $
        minimum = MIN(dispslice,MAX=maximum, /NAN)
		ENDELSE
    ; Display slice
		TV, CONGRID( BYTSCL(dispslice, MIN=minimum, MAX=maximum, /NAN), (*(*info).dispparams).restloopnlxreb,$
      (*(*info).dispparams).restloopntreb, /INTERP), $
			(*(*info).plotpos).restloopx0, (*(*info).plotpos).restloopy0, /NORM
    ; Overplot time indicator
		PLOTS, [(*(*info).plotpos).restloopx0, (*(*info).plotpos).restloopx1], $
      ;[1,1] * ( ((*(*info).dispparams).t-(*(*info).dispparams).t_low) / $
      ;FLOAT((*(*info).dispparams).t_range-1) * $
      [1,1]*( (tarr[(*(*info).dispparams).t]-lower_t_val) / FLOAT(upper_t_val-lower_t_val) * $
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WSET, (*(*info).winids).retrdetwid
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).retrdetwid], labels=['Window ID for draw']
	TVSCL, CONGRID((*(*(*info).loopsdata).det_loopslice)[*,$
    (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_low]:$
    (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_upp]], $
    (*(*info).dispparams).retrdetnlxreb, (*(*info).dispparams).retrdetntreb, /INTERP), $
		(*(*info).plotpos).retrdetx0, (*(*info).plotpos).retrdety0, /NORM
  ; Overplot time indicator
	PLOTS, [(*(*info).plotpos).retrdetx0, (*(*info).plotpos).retrdetx1], $
    ;[1,1] * ( ((*(*info).dispparams).t-(*(*info).dispparams).t_low) / $
    ;FLOAT((*(*info).dispparams).t_range-1) * $
    [1,1]*( ((*(*(*info).dispparams).tarr_main)[(*(*info).dispparams).t]-$
    (*(*info).dispparams).t_low_main) / $
    FLOAT((*(*info).dispparams).t_upp_main-(*(*info).dispparams).t_low_main) * $
    (*(*info).plotpos).retrdetyplspw + (*(*info).plotpos).retrdety0), /NORMAL, COLOR = 100
END

;================================================================================= ESTIMATE SAVING TIME PROCEDURES
PRO CRISPEX_ESTIMATE_TIME_WINDOW, event
; Opens the calculating saving time estimate window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
END

;================================================================================= FIND CRISPEX OUTPUT FILE PROCEDURES
PRO CRISPEX_FIND_CLSAV, event
; Finds CLSAV output files (i.e. saved loop points files)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	imagefilename = STRMID((*(*info).dataparams).imfilename,STRPOS((*(*info).dataparams).imfilename,PATH_SEP(),/REVERSE_SEARCH)+1,STRLEN((*(*info).dataparams).imfilename))
	firstsplit = STRMID(imagefilename,0,STRPOS(imagefilename,'.',/REVERSE_SEARCH))
	fstr = STRSPLIT(firstsplit[0],'_',/EXTRACT)
	IF (N_ELEMENTS(fstr) GE 2) THEN filename = fstr[0]+'_'+fstr[1] ELSE filename = fstr[0]
	cfiles = FILE_SEARCH((*(*info).paths).ipath+filename+"*csav", COUNT = cfilecount)
	IF (*(*info).dataswitch).reffile THEN BEGIN
		refimfilename = STRMID((*(*info).dataparams).refimfilename,$
                      STRPOS((*(*info).dataparams).refimfilename,PATH_SEP(),/REVERSE_SEARCH)+1,$
                      STRLEN((*(*info).dataparams).refimfilename))
		reffirstsplit = STRMID(refimfilename,0,STRPOS(refimfilename,'.',/REVERSE_SEARCH))
		reffstr = STRSPLIT(reffirstsplit[0],'_',/EXTRACT)
		IF (N_ELEMENTS(reffstr) GE 2) THEN refimfilename = reffstr[0]+'_'+reffstr[1] ELSE refimfilename = reffstr[0]
		refcfiles = FILE_SEARCH((*(*info).paths).ipath+refimfilename+"*csav", COUNT = refcfilecount)
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	tempfile = FILEPATH('temp_crispex_redirect.html', /TMP)
	OPENW, lun, tempfile, /GET_LUN
	PRINTF, lun, '<HTML><HEADER><TITLE>CRISPEX help pages</TITLE><META HTTP-EQUIV="REFRESH" CONTENT="0;URL=http://folk.uio.no/gregal/crispex">'+$
		'</HEADER><BODY></BODY></HTML>'
	FREE_LUN, lun
	ONLINE_HELP, BOOK=tempfile
	WAIT, 10.0
	FILE_DELETE, tempfile
END

PRO CRISPEX_HELP_MAIL_BUG, event
; Opens a new message for bug reporting
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

PRO CRISPEX_HELP_MAIL_SUGGESTION, event
; Opens a new message for suggestion reporting
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

PRO CRISPEX_HELP_SHORTCUTS, event
; Opens a window with an overview over the shortcuts
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Populate shortcut elements
  kb_shortcuts = [{sh:'Ctrl+Shift+I', label:'Zoom in'}, $
                  {sh:'Ctrl+Shift+O', label:'Zoom out'}, $
                  {sh:'Shift+B', label:'Step to previous frame'}, $
                  {sh:'Shift+Backspace', label:'Play backwards'}, $
                  {sh:'Shift+Spacebar', label:'Pause'}, $
                  {sh:'Shift+Tab', label:'Play forwards'}, $
                  {sh:'Shift+F', label:'Step to next frame'}, $
                  {sh:'Shift+A', label:'Decrease main '+$
                    STRLOWCASE((*(*info).paramparams).sp_h[(*(*info).plotswitch).heightset])+$
                    ' position'}, $
                  {sh:'Shift+S', label:'Increase main '+$
                    STRLOWCASE((*(*info).paramparams).sp_h[(*(*info).plotswitch).heightset])+$
                    ' position'}, $
                  {sh:'Ctrl+A', label:'Decrease reference '+$
                    STRLOWCASE((*(*info).paramparams).sp_h[(*(*info).plotswitch).refheightset])+$
                    ' position'}, $
                  {sh:'Ctrl+S', label:'Increase reference '+$
                    STRLOWCASE((*(*info).paramparams).sp_h[(*(*info).plotswitch).refheightset])+$
                    ' position'} ]
  ms_shortcuts = [{sh:'Left click', label:'Lock cursor to current position /'}, $
                  {sh:' ', label:'Add current position to path'}, $
                  {sh:'Middle click', label:'Fix cursor position for measurement'}, $
                  {sh:'Right click', label:'Unlock cursor'} ]
	title = 'CRISPEX'+(*(*info).sesparams).instance_label+': Shortcuts'
  base  = WIDGET_BASE(TITLE=title, GROUP_LEADER=(*(*info).winids).root, $
;            XSIZE=(*(*info).winsizes).aboutwinx*1.5, YSIZE=(*(*info).winsizes).aboutwiny, $
            /TLB_FRAME_ATTR, /TLB_KILL_REQUEST_EVENTS)
  disp  = WIDGET_BASE(base, /COLUMN)
  cols  = WIDGET_BASE(disp, /GRID_LAYOUT, COLUMN=2)
  col1  = WIDGET_BASE(cols, /COLUMN);, /FRAME)
  kb_lab= WIDGET_LABEL(col1, VALUE='Keyboard shortcuts:', /ALIGN_LEFT)
  bases1 = WIDGET_BASE(col1, /ROW)
  sh_base1 = WIDGET_BASE(bases1, /COLUMN)
  lb_base1 = WIDGET_BASE(bases1, /COLUMN)
  FOR i=0,N_ELEMENTS(kb_shortcuts.sh)-1 DO BEGIN
    kb_label= WIDGET_LABEL(sh_base1, VALUE=(kb_shortcuts.sh)[i], /ALIGN_LEFT)
    lb_label= WIDGET_LABEL(lb_base1, VALUE=(kb_shortcuts.label)[i], /ALIGN_LEFT)
  ENDFOR
  col2  = WIDGET_BASE(cols, /COLUMN);, /FRAME)
  ms_lab= WIDGET_LABEL(col2, VALUE='Mouse shortcuts:', /ALIGN_LEFT)
  bases2 = WIDGET_BASE(col2, /ROW)
  sh_base2 = WIDGET_BASE(bases2, /COLUMN)
  lb_base2 = WIDGET_BASE(bases2, /COLUMN)
  FOR i=0,N_ELEMENTS(ms_shortcuts.sh)-1 DO BEGIN
    ms_label= WIDGET_LABEL(sh_base2, VALUE=(ms_shortcuts.sh)[i], /ALIGN_LEFT)
    lb_label= WIDGET_LABEL(lb_base2, VALUE=(ms_shortcuts.label)[i], /ALIGN_LEFT)
  ENDFOR
  close_base = WIDGET_BASE(disp, /ALIGN_CENTER)
  close_button = WIDGET_BUTTON(close_base, VALUE='Close', EVENT_PRO='CRISPEX_CLOSE_EVENT_WINDOW')
  WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET=(*(*info).winsizes).aboutxoffset, $
    TLB_SET_YOFFSET=(*(*info).winsizes).aboutyoffset
  (*(*info).winids).shorttlb = base
  WIDGET_CONTROL, base, SET_UVALUE=info
  XMANAGER, 'CRISPEX', base, /NO_BLOCK
END

;================================================================================= INTENSITY-TIME SAVE PROCEDURES
PRO CRISPEX_INT_SAVE, event
; Handles the actual saving of the intensity-time plots
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=intfilename, /tlab, ext='cint'
	condition = WHERE(*(*(*info).intparams).sel_diagnostics EQ 1)
	intensities = FLTARR((*(*info).dispparams).t_range,N_ELEMENTS(condition))
	avg_intensity = FLTARR(N_ELEMENTS(condition))
	FOR i=0,N_ELEMENTS(condition)-1 DO BEGIN
		intensities[0,i] = REFORM( ( ( *(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * (*(*info).dataparams).ns + FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns + $
			(*(*info).dataparams).s ] )[condition[i],(*(*info).dispparams).t_low:(*(*info).dispparams).t_upp] )
		avg_intensity[i] = MEAN(intensities[*,i], /NAN)
	ENDFOR
	diagnostics = ((*(*info).intparams).diagnostics)[condition]
	x = (*(*info).dataparams).x		&	y = (*(*info).dataparams).y
	nt = (*(*info).dispparams).t_range	& 	dt = (*(*info).plotaxes).dt
	t_low = (*(*info).dispparams).t_low	&	t_upp = (*(*info).dispparams).t_upp
	t_saved = (*(*info).dispparams).t	&	crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
	SAVE, crispex_version, intensities, avg_intensity, diagnostics, nt, dt, t_low, t_upp, t_saved, x, y, FILENAME=(*(*info).paths).opath+intfilename
	PRINT, 'Written: '+(*(*info).paths).opath+intfilename+'.cint'
END


;================================================================================= INTERRUPT PROCEDURE
PRO CRISPEX_INTERRUPT, event
; Handles interrupting at runtime
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	PRINT,'Interrupted CRISPEX at runtime. Type [.c] to continue...'
	STOP
END

;========================= INPUT/OUTPUT PROCEDURES
PRO CRISPEX_IO_FAILSAFES_MAIN, imcube, spcube, input_single_cube, $
                               HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                               STARTUPTLB=startuptlb, $
                               IO_FAILSAFE_ERROR=io_failsafe_error
  hdr_out = hdr_in
  io_failsafe_error = 0
  ; If SPCUBE has been supplied, check IMCUBE and SPCUBE compatibility
	IF hdr_out.spfile THEN BEGIN         
    ; Check whether SPCUBE and IMCUBE are actually the same
		IF ((spcube EQ imcube) OR ((hdr_in.nlp EQ hdr_in.nx) AND (hdr_in.mainnt EQ hdr_in.ny) AND $
       (hdr_in.spnt EQ hdr_in.imnt))) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'IMCUBE and SPCUBE must be different. Please check '+$
        'input (you seem to have provided the same file twice).',/ERROR,/NO_ROUTINE
			IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_error = 1
			RETURN
		ENDIF
    ; Check whether condensed third dimensions of SPCUBE and IMCUBE are incompatible
		IF ((hdr_in.nx*hdr_in.ny*hdr_in.ns NE hdr_in.spnt) OR $
       (hdr_in.mainnt*hdr_in.nlp*hdr_in.ns NE hdr_in.imnt)) THEN BEGIN							
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'IMCUBE and SPCUBE have incompatible dimensions and '+$
        'seem to belong to different datasets. Please check whether the input is correct (you '+$
        'provided IMCUBE='+STRTRIM(imcube,2)+' and SPCUBE='+STRTRIM(spcube,2)+').',/ERROR,$
        /NO_ROUTINE
			IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_error = 1
			RETURN
		ENDIF
    ; Check, for multichannel cube, whether channels are incompatible
		IF hdr_out.multichannel THEN BEGIN
			IF ((hdr_in.spstokes NE hdr_in.imstokes) OR (hdr_in.spns NE hdr_in.imns)) THEN BEGIN
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'IMCUBE and SPCUBE have incompatible number of '+$
          'channels and seem to belong to different datasets. Please check whether the input is '+$
          'correct (you provided IMCUBE='+STRTRIM(imcube,2)+' and SPCUBE='+STRTRIM(spcube,2)+').',$
          /ERROR,/NO_ROUTINE
				IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
        io_failsafe_error = 1
				RETURN
			ENDIF ELSE $
        hdr_out.ns = hdr_in.imns  ; If they are compatible, set general number of channels
		ENDIF ELSE hdr_out.ns = 1L       ; If not multiple channels, set general number channels to 1
    ; Check whether SINGLE_CUBE keyword has been set in combination with provided SPCUBE
		IF (N_ELEMENTS(INPUT_SINGLE_CUBE) GT 0) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Calling CRISPEX with SINGLE_CUBE, while SPCUBE is '+$
        'provided, is not allowed. SINGLE_CUBE keyword will be ignored.', /WARNING, /NO_ROUTINE
      io_failsafe_error = 2
			hdr_out.onecube = 0
		ENDIF
		hdr_out.single_cube[0] = 0
	ENDIF ;ELSE BEGIN  ; If no SPCUBE has been supplied, check other settings
;    ; Check setting of SINGLE_CUBE keyword
;		IF (N_ELEMENTS(INPUT_SINGLE_CUBE) EQ 1) THEN BEGIN  ; If SINGLE_CUBE properly set, set nlp and nt
;			hdr_out.nlp = LONG(INPUT_SINGLE_CUBE)
;			hdr_out.onecube = 1
;      hdr_out.single_cube[0] = hdr_out.nlp
;			hdr_out.nt = hdr_in.imnt / hdr_in.nlp / hdr_in.ns
;		ENDIF ELSE BEGIN  ; If no SPCUBE or SINGLE_CUBE are set, we are dealing with a snapshot
;			hdr_out.nlp = hdr_in.imnt / hdr_in.ns
;			hdr_out.single_cube[0] = 0
;		ENDELSE
;	ENDELSE
END

PRO CRISPEX_IO_FAILSAFES_MAIN_REF, HDR=hdr, STARTUPTLB=startuptlb, $
                                   IO_FAILSAFE_ERROR=io_failsafe_error
;  io_failsafe_error = 0
  IF (((hdr.refnx EQ hdr.nx) AND (hdr.refny EQ hdr.ny)) $   ; Require equal dimensions and require:
    AND ((hdr.refnt EQ hdr.mainnt) OR $                         ; either same number of timesteps
         (hdr.refnt EQ 1) OR $                              ; or reference nt being 1
         (hdr.refimnt EQ hdr.imnt) OR $                     ; or third dimensions being equal
         (hdr.refimnt EQ hdr.nlp) OR $                      ; or reference being snapshot
         (hdr.mainnt EQ 1))) THEN BEGIN                         ; or main being snapshot
    io_failsafe_error = 0
  ENDIF ELSE BEGIN                                          ; If not, throw error messages and exit
    IF ((hdr.refnt NE hdr.mainnt) AND (hdr.refnt NE 1)) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the reference cube (['+$
        STRTRIM(hdr.refnx,2)+','+STRTRIM(hdr.refny,2)+','+STRTRIM(hdr.refnt,2)+']) are not '+$
        'compatible with those of the main cube (['+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+','+$
        STRTRIM(hdr.imnt,2)+'])! Number of timesteps times the number of spectral positions '+$
        'must be equal to that of the main cube ('+STRTRIM(hdr.imnt,2)+').', /ERROR, /NO_ROUTINE, $
        /NEWLINE
    ENDIF ELSE BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the reference cube (['+$
        STRTRIM(hdr.refnx,2)+','+STRTRIM(hdr.refny,2)+','+STRTRIM(hdr.refnt,2)+']) are not '+$
        'compatible with those of the main cube (['+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+$
        ','+STRTRIM(hdr.mainnt,2)+'])! Number of timesteps must be either equal to that of the main '+$
        'cube ('+STRTRIM(hdr.mainnt,2)+') or to 1.', /ERROR, /NO_ROUTINE, /NEWLINE
	  ENDELSE
		IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
    io_failsafe_error = 1
		RETURN
  ENDELSE
END

PRO CRISPEX_IO_FAILSAFES_REF, refcube, input_single_cube, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                              STARTUPTLB=startuptlb, $
                              IO_FAILSAFE_ERROR=io_failsafe_error
; Handles failsafes against wrongly supplied reference image and spectral cubes
  hdr_out = hdr_in
  io_failsafe_error = 0
  IF (N_ELEMENTS(REFCUBE) EQ 2) THEN BEGIN
    ; Failsafe against providing the same file as REFIMCUBE and REFSPCUBE
  	IF ((refcube[1] EQ refcube[0]) OR $         ; Check whether input cube names are the same, or
      (hdr_out.refnlp EQ hdr_out.refspnx) AND $         ; the cubes have same first dimensions
      (hdr_out.refnt EQ hdr_out.refspny) AND $          ; the cubes have same second dimenions
      (hdr_out.refspnt EQ hdr_out.refimnt)) THEN BEGIN  ; the cubes have same third dimensions
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'The reference image and spectral cubes must be '+$
        'different. Please check input (you seem to have provided the same file twice).', /ERROR, $
        /NO_ROUTINE, /NEWLINE
  		IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_error = 1
  		RETURN
  	ENDIF
    ; Failsafe against providing incompatible REFIMCUBE and REFSPCUBE
  	IF ((hdr_out.refnx*hdr_out.refny NE hdr_out.refspnt) OR $
        (hdr_out.refnt*hdr_out.refnlp NE hdr_out.refimnt)) THEN BEGIN							
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'The reference image and spectral cubes have '+$
        'incompatible dimensions and seem to belong to different datasets. Please check whether '+$
        'the input is correct (you provided REFCUBE[0]='+STRTRIM(refcube[0],2)+' and REFCUBE[1]='+$
        STRTRIM(refcube[1],2)+').', /ERROR, /NO_ROUTINE, /NEWLINE
  	  IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_error = 1
  		RETURN
  	ENDIF
    ; Check whether SINGLE_CUBE keyword has been set in combination with provided SPCUBE
		IF (N_ELEMENTS(INPUT_SINGLE_CUBE) GT 0) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Calling CRISPEX with SINGLE_CUBE with a 2-element '+$
        'array, while a reference SPCUBE is provided, is not allowed. SINGLE_CUBE keyword will '+$
        'be ignored.', /WARNING, /NO_ROUTINE
      io_failsafe_error = 2
		ENDIF
		hdr_out.single_cube[1] = 0
	ENDIF ;ELSE BEGIN  ; If no REFSPCUBE has been supplied, check other settings
;    ; Check setting of SINGLE_CUBE keyword
;		IF (N_ELEMENTS(INPUT_SINGLE_CUBE) EQ 1) THEN BEGIN  ; If SINGLE_CUBE properly set, set refnlp 
;			hdr_out.refnlp = LONG(INPUT_SINGLE_CUBE)
;			hdr_out.onecube = 1
;      hdr_out.single_cube[1] = hdr_out.refnlp
;		ENDIF 
;  ENDELSE
END

PRO CRISPEX_IO_FAILSAFES_MAIN_SJI, sjicube, HDR=hdr, STARTUPTLB=startuptlb, $
                              IO_FAILSAFE_ERROR=io_failsafe_error
; Handles failsafes against wrongly supplied slit-jaw image cube
  IF ((hdr.sjinx GE hdr.nx) AND (hdr.sjiny GE hdr.ny)) THEN $ ; Require ge xy-dimensions than main
    io_failsafe_error = 0 $
	ELSE BEGIN
	  IF ((hdr.sjinx LT hdr.nx) OR (hdr.sjiny LT hdr.ny)) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the slit-jaw image cube (['+$
        STRTRIM(hdr.sjinx,2)+','+STRTRIM(hdr.sjiny,2)+','+STRTRIM(hdr.sjint,2)+']) are not '+$
        'compatible with those of the main image cube (['+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+$
        ','+STRTRIM(hdr.mainnt,2)+'])! Number of pixels in the x- and y-dimension must greater than '+$
        'or equal to those of the main image cube.', /ERROR, /NO_ROUTINE, /NEWLINE
	  ENDIF 
    IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
    io_failsafe_error = 1
		RETURN
  ENDELSE
END

PRO CRISPEX_IO_FAILSAFES_MASK, maskcube, HDR=hdr, STARTUPTLB=startuptlb, $
                              IO_FAILSAFE_ERROR=io_failsafe_error
; Handles failsafes against wrongly supplied mask cube
  IF ((hdr.masknx EQ hdr.nx) AND (hdr.maskny EQ hdr.ny) AND $ ; Require same xy-dimensions as main
    ((hdr.masknt EQ hdr.mainnt) OR (hdr.masknt EQ 1))) THEN BEGIN ; Require same as main or 1 timestep
    io_failsafe_error = 0
	ENDIF ELSE BEGIN
	  IF ((hdr.masknx NE hdr.nx) OR (hdr.maskny NE hdr.ny)) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the mask cube (['+$
        STRTRIM(hdr.masknx,2)+','+STRTRIM(hdr.maskny,2)+','+STRTRIM(hdr.masknt,2)+']) are not '+$
        'compatible with those of the main image cube (['+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+$
        ','+STRTRIM(hdr.mainnt,2)+'])! Number of pixels in the x- and y-dimension must be equal to '+$
        'those of the main image cube.', /ERROR, /NO_ROUTINE, /NEWLINE
	  ENDIF ELSE IF ((hdr.masknt NE hdr.mainnt) AND (hdr.masknt NE 1)) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the mask cube (['+$
        STRTRIM(hdr.masknx,2)+','+STRTRIM(hdr.maskny,2)+','+STRTRIM(hdr.masknt,2)+']) are not '+$
        'compatible with those of the main image cube (['+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+$
        ','+STRTRIM(hdr.mainnt,2)+'])! Number of timesteps must be either equal to that of the main '+$
        'image cube ('+STRTRIM(hdr.mainnt,2)+') or to 1.', /ERROR, /NO_ROUTINE, /NEWLINE
	  ENDIF
    IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
    io_failsafe_error = 1
		RETURN
  ENDELSE
END

PRO CRISPEX_IO_FAILSAFES_MNSPEC, mnspec, hdr, STARTUPTLB=startuptlb, $
                                 IO_FAILSAFE_ERROR=io_failsafe_error
; Handles failsafes against wrongly supplied MNSPEC values                                 
  extra_text = ''
  cumulative_error = (mnspec[0] LT 0)                               ; Check whether lower value >= 0
  cumulative_error += (mnspec[(N_ELEMENTS(MNSPEC) EQ 2)] GE hdr.mainnt) ; Check whether upper value < nt
  IF (N_ELEMENTS(MNSPEC) EQ 2) THEN BEGIN
    feedback_text = 's (['+STRTRIM(mnspec[0],2)+','+STRTRIM(mnspec[1],2)+'])'
    IF (mnspec[0] GT mnspec[1]) THEN $                     ; Check whether upper value > lower value
      extra_text = ' and be ordered from lower to higher value'
  ENDIF ELSE feedback_text = ' ('+STRTRIM(mnspec[0],2)+')' 
  io_failsafe_error = (TOTAL(cumulative_error) NE 0)
  IF io_failsafe_error THEN BEGIN
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'MNSPEC value'+feedback_text+' must fall within '+$
      'allowed range [0,'+STRTRIM(LONG(hdr.mainnt-1),2)+']'+extra_text+'!',/ERROR,/NO_ROUTINE,/NEWLINE
	  IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
	  RETURN
  ENDIF
END

PRO CRISPEX_IO_FAILSAFES_LINE_CENTER, line_center, hdr, NFILES=nfiles, STARTUPTLB=startuptlb, $
                                      SPECTFILE_SET=spectfile_set, $
                                      REFSPECTFILE_SET=refspectfile_set, $
                                      IO_FAILSAFE_ERROR=io_failsafe_error
; Handles failsafes against wrongly supplied LINE_CENTER values or formatting
  ndims = SIZE(LINE_CENTER,/N_DIMENSIONS) > 1
  nelem = N_ELEMENTS(LINE_CENTER)
  lcase = nelem / ndims
  nlp_select = [hdr.nlp,hdr.refnlp] & feedback_text = ['Main','Reference']
  io_failsafe_error = 0
  ; Check whether conflicting with SPECTFILE
  IF ((spectfile_set AND (ndims GE 1)) OR (refspectfile_set AND (ndims EQ 2))) THEN BEGIN
    idx = (refspectfile_set AND (ndims EQ 2))
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Detailed '+STRLOWCASE(feedback_text[idx])+$
      ' information may not be specified in both SPECTFILE and LINE_CENTER! Settings from '+$
      'LINE_CENTER will be ignored.',/WARNING,/NO_ROUTINE
    io_failsafe_error = 2
		RETURN
  ENDIF
  ; Check proper setting of LINE_CENTER
  IF (((lcase EQ 1) OR (lcase EQ 3)) AND (ndims LE nfiles)) THEN BEGIN
  ; lcase of 1: LINE_CENTER = lc or [[lc],[lc]]
  ; lcase of 3: LINE_CENTER = [lc,cwav,dwav] or [[lc,cwav,dwav],[lc,cwav,dwav]]
    FOR d=0,ndims-1 DO BEGIN
      IF ((line_center[0,d] GE nlp_select[d]) OR (line_center[0,d] LT 0)) THEN BEGIN
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, feedback_text[d]+' linecentre index value '+$
          STRTRIM(line_center[0,d],2)+' falls outside of allowed range [0,'+$
          STRTRIM(nlp_select[d]-1,2)+']!',/ERROR,/NO_ROUTINE,NEWLINE=(io_failsafe_error EQ 2)
        IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
        io_failsafe_error = 1
        RETURN
      ENDIF
    ENDFOR
  ENDIF ELSE IF (lcase NE 2) THEN BEGIN
  ; lcase of 2: LINE_CENTER = [cwav,dwav] or [[cwav,dwav],[cwav,dwav]]
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'LINE_CENTER keyword contains too many elements for '+$
      'the number of cubes provided!', /ERROR, /NO_ROUTINE, NEWLINE=(io_failsafe_error EQ 2)
    IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
    io_failsafe_error = 1
    RETURN
  ENDIF
END

PRO CRISPEX_IO_MULTICHANNEL, CUBE_COMPATIBILITY=cube_compatibility, CHANNELS_LABELS=channels_labels
; Handles setting of parameters in case of multichannel input
  IF KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN
  ENDIF
END

PRO CRISPEX_IO_FEEDBACK, verbosity, hdr, IMCUBE=imcube, SPCUBE=spcube, REFIMCUBE=refimcube, $
                         REFSPCUBE=refspcube, MASKCUBE=maskcube, SJICUBE=sjicube
	multichannel = (hdr.ns GE 2)
	IF (TOTAL(verbosity[0:1]) GE 1) THEN BEGIN
    IF (N_ELEMENTS(IMCUBE) EQ 1) THEN BEGIN
      IF ((SIZE(IMCUBE,/TYPE) NE 0) AND (STRCOMPRESS(IMCUBE) NE '')) THEN BEGIN
		    IF multichannel THEN $
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read Stokes image cube: '+imcube+$
            '. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+','+$
            STRTRIM(hdr.imnt,2)+').' $
        ELSE $
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read image cube: '+imcube+$
            '. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+','+$
            STRTRIM(hdr.imnt,2)+').'
		    IF (verbosity[1] EQ 1) THEN $
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Main cube(s) dimensions: (nx,ny,nt,nlp,ns) = ('+$
            STRTRIM(hdr.nx,2)+','+STRTRIM(hdr.ny,2)+','+STRTRIM(hdr.mainnt,2)+','+STRTRIM(hdr.nlp,2)+$
            ','+STRTRIM(hdr.ns,2)+')'
	    ENDIF ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No main image cube supplied.'
    ENDIF ELSE IF (N_ELEMENTS(SPCUBE) EQ 1) THEN BEGIN
      IF ((SIZE(SPCUBE,/TYPE) NE 0) AND (STRCOMPRESS(SPCUBE) NE '')) THEN $
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read main spectral cube: '+spcube+$
              '. Dimensions: (nlp,nt,nx*ny*ns) = ('+STRTRIM(hdr.nlp,2)+','+STRTRIM(hdr.mainnt,2)+$
              ','+STRTRIM(hdr.spnt,2)+').' $
      ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No main spectral cube supplied.'
    ENDIF ELSE IF (N_ELEMENTS(REFIMCUBE) EQ 1) THEN BEGIN
      IF ((SIZE(REFIMCUBE,/TYPE) NE 0) AND (STRCOMPRESS(REFIMCUBE) NE '')) THEN BEGIN
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read reference image cube: '+refimcube+$
              '. Dimensions: (nx,ny,nt*nlp*ns) = ('+STRTRIM(hdr.refnx,2)+','+STRTRIM(hdr.refny,2)+$
              ','+STRTRIM(hdr.refnt*hdr.refnlp*hdr.refns,2)+').'  
			  IF (verbosity[1] EQ 1) THEN $
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Reference cubes dimensions: (nx,ny,nt,nlp,ns) = ('+$
                STRTRIM(hdr.refnx,2)+','+STRTRIM(hdr.refny,2)+','+STRTRIM(hdr.refnt,2)+','+$
                STRTRIM(hdr.refnlp,2)+','+STRTRIM(hdr.refns,2)+')'
      ENDIF ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No reference image cube supplied.'
    ENDIF ELSE IF (N_ELEMENTS(REFSPCUBE) EQ 1) THEN BEGIN
      IF ((SIZE(REFSPCUBE,/TYPE) NE 0) AND (STRCOMPRESS(REFSPCUBE) NE '')) THEN $
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read reference spectral cube: '+refspcube+$
              '. Dimensions: (nlp,nt,nx*ny) = ('+STRTRIM(hdr.refnlp,2)+','+STRTRIM(hdr.refnt,2)+$
              ','+STRTRIM(hdr.refspnt,2)+').' $
      ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No reference spectral cube supplied.'
    ENDIF ELSE IF (N_ELEMENTS(MASKCUBE) EQ 1) THEN BEGIN
      IF ((SIZE(MASKCUBE,/TYPE) NE 0) AND (STRCOMPRESS(MASKCUBE) NE '')) THEN $
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read mask cube: '+maskcube+$
          '. Dimensions: (nx,ny,nt) = ('+STRTRIM(hdr.masknx,2)+','+STRTRIM(hdr.maskny,2)+','+$
          STRTRIM(hdr.masknt,2)+').' $
      ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No mask cube supplied.'
    ENDIF ELSE IF (N_ELEMENTS(SJICUBE) EQ 1) THEN BEGIN
      IF ((SIZE(SJICUBE,/TYPE) NE 0) AND (STRCOMPRESS(SJICUBE) NE '')) THEN $
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Read slit-jaw image cube: '+sjicube+$
          '. Dimensions: (nx,ny,nt) = ('+STRTRIM(hdr.sjinx,2)+','+STRTRIM(hdr.sjiny,2)+','+$
          STRTRIM(hdr.sjint,2)+').' $
      ELSE CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No slit-jaw image cube supplied.'
    ENDIF
  ENDIF
END

PRO CRISPEX_IO_OPEN_MAINCUBE, IMCUBE=imcube, SPCUBE=spcube, SINGLE_CUBE=single_cube, $
                              HDR_IN=hdr_in, HDR_OUT=hdr_out, STARTUPTLB=startuptlb, $
                              IO_FAILSAFE_MAIN_ERROR=io_failsafe_main_error
  hdr_out = hdr_in
  ipath = hdr_out.ipath
  instance_label = hdr_out.instance_label
  hdr_out.imfilename = imcube
  ; Determine cube compatibility mode for inputfiles (0: running FITS cubes, 1: running old cubes)
	imext = STRMID(hdr_out.imfilename,STRPOS(hdr_out.imfilename,'.',/REVERSE_SEARCH)+1,$
                  STRLEN(hdr_out.imfilename))  ; Process extension
	hdr_out.imcube_compatibility = ABS(STRMATCH(imext,'fits',/FOLD_CASE)-1)             ; Determine comp mode
  IF (N_ELEMENTS(SPCUBE) EQ 1) THEN BEGIN
    hdr_out.spfilename = spcube
    spext = STRMID(hdr_out.spfilename,STRPOS(hdr_out.spfilename,'.',/REVERSE_SEARCH)+1,$
                    STRLEN(hdr_out.spfilename))
    hdr_out.spcube_compatibility = ABS(STRMATCH(spext,'fits',/FOLD_CASE)-1)
    hdr_out.spfile = 1
;  ENDIF 
;  IF N_ELEMENTS(SPCUBE) EQ 1 THEN BEGIN
    CRISPEX_IO_PARSE_HEADER, hdr_out.spfilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                            CUBE_COMPATIBILITY=hdr_out.spcube_compatibility, EXTEN_NO=0, /SPCUBE
	ENDIF ELSE hdr_out.onecube = 1                       ; onecube switch if no SPCUBE has been provided
  ; If single_cube value has been set from single FITS cube, use that
  IF ((hdr_out.imcube_compatibility EQ 0) AND (N_ELEMENTS(SPCUBE) NE 1)) THEN $
    main_single_cube = hdr_out.single_cube[0] $
  ELSE IF (N_ELEMENTS(SINGLE_CUBE) GE 1) THEN $
    main_single_cube = single_cube[0]
  CRISPEX_IO_PARSE_HEADER, hdr_out.imfilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                            CUBE_COMPATIBILITY=hdr_out.imcube_compatibility, EXTEN_NO=0, /IMCUBE, $
                            SINGLE_CUBE=main_single_cube
	hdr_out.multichannel = (hdr_out.ns GE 2)
;  ; If single_cube value has been set from single FITS cube, use that
;  IF ((hdr_out.single_cube[0] NE 0) AND (N_ELEMENTS(SPCUBE) NE 1)) THEN $
;    main_single_cube = hdr_out.single_cube[0] $
;  ELSE IF (N_ELEMENTS(SINGLE_CUBE) GE 1) THEN $
;    main_single_cube = single_cube[0]
  CRISPEX_IO_FAILSAFES_MAIN, hdr_out.imfilename, hdr_out.spfilename, main_single_cube, $
                             HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                             STARTUPTLB=startuptlb, $
                             IO_FAILSAFE_ERROR=io_failsafe_main_error
  IF (io_failsafe_main_error EQ 1) THEN RETURN
  ; Define axes titles based on IMCUBE and SPCUBE headers
  IF (STRCOMPRESS(hdr_out.bunit,/REMOVE_ALL) NE '') THEN $
    ytitle_unit = ' ['+hdr_out.bunit+']' ELSE ytitle_unit = ''
  IF (STRCOMPRESS(hdr_out.lpunit,/REMOVE_ALL) NE '') THEN $
    xtitle_unit = ' ['+hdr_out.lpunit+']' ELSE xtitle_unit = ''
  IF (STRCOMPRESS(hdr_out.tunit,/REMOVE_ALL) NE '') THEN $
    spytitle_unit = ' ['+hdr_out.tunit+']' ELSE spytitle_unit = ''
  hdr_out.ytitle[0] = hdr_out.blabel+ytitle_unit
  hdr_out.xtitle[0] = hdr_out.lplabel+xtitle_unit
  hdr_out.spytitle = hdr_out.tlabel+spytitle_unit
  ; Handle Stokes 
  IF hdr_out.multichannel THEN BEGIN
    stokes_labels = STRSPLIT(STRMID(hdr_out.imstokes,1,STRLEN(hdr_out.imstokes)-2),',',/EXTRACT)
		IF (N_ELEMENTS(stokes_labels) NE hdr_out.imns) THEN BEGIN
			PRINT,'ERROR: The number of Stokes components ('+STRTRIM(hdr_out.imns,2)+') does not '+$
            'correspond to the number of Stokes labels ('+STRTRIM(N_ELEMENTS(stokes_labels),2)+').'
			PRINT,'       Please check whether the Stokes cube production has proceded correctly.'
			WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_main_error = 1
			RETURN
		ENDIF ELSE BEGIN
			stokes_select_sp = INTARR(hdr_out.ns)
			IF ((WHERE(stokes_labels EQ 'I') GE 0) AND $
          (WHERE(stokes_labels EQ 'I') LE hdr_out.imns-1)) THEN BEGIN
				hdr_out.stokes_enabled[0] = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'I')] = 1
			ENDIF 
			IF ((WHERE(stokes_labels EQ 'Q') GE 0) AND $
          (WHERE(stokes_labels EQ 'Q') LE hdr_out.imns-1)) THEN BEGIN
				hdr_out.stokes_enabled[1] = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'Q')] = 1
			ENDIF 
			IF ((WHERE(stokes_labels EQ 'U') GE 0) AND $
          (WHERE(stokes_labels EQ 'U') LE hdr_out.imns-1)) THEN BEGIN
				hdr_out.stokes_enabled[2] = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'U')] = 1
			ENDIF 
			IF ((WHERE(stokes_labels EQ 'V') GE 0) AND $
          (WHERE(stokes_labels EQ 'V') LE hdr_out.imns-1)) THEN BEGIN
				hdr_out.stokes_enabled[3] = 1 
				stokes_select_sp[WHERE(stokes_labels EQ 'V')] = 1
			ENDIF
		ENDELSE
	ENDIF ELSE BEGIN
    stokes_labels = ['I']
    stokes_select_sp = 1
  ENDELSE
	hdr_out.scalestokes_max = (TOTAL(hdr_out.stokes_enabled[1:3]) GE 1)
  hdr_out = CREATE_STRUCT(hdr_out, 'stokes_labels', stokes_labels, 'stokes_select_sp', $
                          stokes_select_sp)
	IF (hdr_out.verbosity[1] EQ 1) THEN $
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Stokes parameters: '+STRJOIN(hdr_out.stokes_labels,' ')
  ; Handle diagnostics
  IF (WHERE(TAG_NAMES(hdr_out) EQ 'DIAG_START') EQ -1) THEN BEGIN
    wstart = 0
    wwidth = hdr_out.nlp
    hdr_out = CREATE_STRUCT(hdr_out, 'diag_start', wstart, 'diag_width', wwidth)
  ENDIF
;	diagnostics = STRARR(hdr_out.nlp)
;	IF (hdr_out.imdiagnostics NE '') THEN BEGIN
;		diagsplit = STRSPLIT(hdr_out.imdiagnostics,',',/EXTRACT)
;		ndiag = N_ELEMENTS(diagsplit)
;		diagsplit[0] = STRMID(diagsplit[0],STRPOS(diagsplit[0],'[')+1,STRLEN(diagsplit[0]))
;		diagsplit[ndiag-1] = STRMID(diagsplit[ndiag-1],0,$
;                                STRLEN(diagsplit[ndiag-1])-(STRPOS(diagsplit[ndiag-1],']') GT 0))
;		IF (ndiag GE hdr_out.nlp) THEN $
;      diagnostics = diagsplit[0:(hdr_out.nlp-1)] $
;    ELSE IF (ndiag LT hdr_out.nlp) THEN BEGIN
;			diagnostics[0:(ndiag-1)] = diagsplit[0:(ndiag-1)]
;			diagnostics[ndiag:(hdr_out.nlp-1)] = REPLICATE('Undefined',(hdr_out.nlp-ndiag))
;		ENDIF 
;	ENDIF ELSE diagnostics = REPLICATE('SST ',hdr_out.nlp)+STRTRIM(INDGEN(hdr_out.nlp),2)
	sel_diagnostics = REPLICATE(1,hdr_out.nlp)
	lines_diagnostics = (INDGEN(hdr_out.nlp) MOD 6)
	FOR i=0,FLOOR(hdr_out.nlp/6.) DO BEGIN
		IF (i EQ 0) THEN selcol_diagnostics = REPLICATE(i,6) ELSE $
			IF (i EQ (FLOOR(hdr_out.nlp/6.))) THEN BEGIN
				IF (hdr_out.nlp/6. NE FLOOR(hdr_out.nlp/6.)) THEN selcol_diagnostics = [selcol_diagnostics,$
        REPLICATE(i,ROUND((hdr_out.nlp/6.-FLOOR(hdr_out.nlp/6.))*6))] 
			ENDIF ELSE selcol_diagnostics = [selcol_diagnostics, REPLICATE(i,6)]
	ENDFOR
  tsel_main = INDGEN(hdr_out.mainnt)
;  hdr_out = CREATE_STRUCT(hdr_out, 'diagnostics', diagnostics, 'sel_diagnostics', sel_diagnostics, $
  hdr_out = CREATE_STRUCT(hdr_out, 'sel_diagnostics', sel_diagnostics, $
    'lines_diagnostics', lines_diagnostics, 'selcol_diagnostics', selcol_diagnostics, $
    'tsel_main', tsel_main)
  CRISPEX_IO_OPEN_MAINCUBE_READ, HDR_IN=hdr_out, HDR_OUT=hdr_out
END

PRO CRISPEX_IO_OPEN_MAINCUBE_READ, HDR_IN=hdr_in, HDR_OUT=hdr_out
  hdr_out = hdr_in
  ; Actual read-in of the data cubes  
  ; Read in main image cube 
	OPENR, lunim, hdr_out.imfilename, /get_lun, $
         SWAP_ENDIAN = ((hdr_out.imtype GT 1) AND (hdr_out.endian NE hdr_out.imendian))									
	imagefile = ASSOC(lunim,MAKE_ARRAY(hdr_out.nx,hdr_out.ny, TYPE=hdr_out.imtype),hdr_out.imoffset)
  ; Re-read in of the image cube for slices
	scanfile  = ASSOC(lunim,MAKE_ARRAY(hdr_out.nx,hdr_out.ny,hdr_out.nlp*hdr_out.ns, $
                      TYPE=hdr_out.imtype),hdr_out.imoffset)			
	hdr_out.imdata	= PTR_NEW(imagefile, /NO_COPY)
	hdr_out.scan	= PTR_NEW(scanfile, /NO_COPY)
  hdr_out.lunim = lunim
  hdr_out.dx_fixed = ABS(KEYWORD_SET(hdr_out.imcube_compatibility)-1)
  CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, IMCUBE=hdr_out.imfilename
  ; Read in main spectral cube 
  IF hdr_out.spfile THEN BEGIN
    OPENR, lunsp, hdr_out.spfilename, /get_lun, $
           SWAP_ENDIAN = ((hdr_out.sptype GT 1) AND (hdr_out.endian NE hdr_out.spendian))
    spectra = ASSOC(lunsp,MAKE_ARRAY(hdr_out.nlp,hdr_out.mainnt, TYPE=hdr_out.sptype),hdr_out.spoffset)
    hdr_out.lunsp = lunsp
	  hdr_out.spdata = PTR_NEW(spectra, /NO_COPY) 
  ENDIF
  CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, SPCUBE=hdr_out.spfilename
END

PRO CRISPEX_IO_OPEN_REFCUBE, REFCUBE=refcube, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                              SINGLE_CUBE=single_cube, $
                              CUBE_COMPATIBILITY=cube_compatibility, $
                              IO_FAILSAFE_REF_ERROR=io_failsafe_ref_error, $
                              IO_FAILSAFE_MAIN_REF_ERROR=io_failsafe_main_ref_error
  io_failsafe_ref_error = 0
  io_failsafe_main_ref_error = 0
;  IF (N_ELEMENTS(event) EQ 1) THEN BEGIN
;    WIDGET_CONTROL, event.TOP, GET_UVALUE = info
;    ipath = (*(*info).paths).ipath
;    instance_label = (*(*info).sesparams).instance_label
;    hdr_out = (*(*info).ioparams).hdr
;  ENDIF ELSE BEGIN
    hdr_out = hdr_in
    ipath = hdr_out.ipath
    instance_label = hdr_out.instance_label
;  ENDELSE
;  IF ((N_ELEMENTS(event) EQ 1) AND (N_ELEMENTS(REFCUBE) LT 1)) THEN BEGIN
;    refcube = STRARR(2)
;	  refcube[0] = DIALOG_PICKFILE(/READ,/MUST_EXIST, PATH = ipath, TITLE='CRISPEX'+instance_label+$
;                                ': Select reference image cube', FILTER=['*cube','*fits'])
;	  refcube[1] = DIALOG_PICKFILE(/READ,/MUST_EXIST, PATH = ipath, TITLE='CRISPEX'+instance_label+$
;                                ': Select reference spectral cube', FILTER=['*cube','*fits'])
;  ENDIF
	IF ((N_ELEMENTS(REFCUBE) GE 1) AND (SIZE(REFCUBE,/TYPE) EQ 7)) THEN BEGIN					
    hdr_out.refimfilename = refcube[0]
  	refimext = STRMID(hdr_out.refimfilename,STRPOS(hdr_out.refimfilename,'.',/REVERSE_SEARCH)+1,$
                      STRLEN(hdr_out.refimfilename))
  	hdr_out.refimcube_compatibility = ABS(STRMATCH(refimext,'fits',/FOLD_CASE)-1)
    ; If single_cube value has been set from single FITS cube, use that
    IF ((hdr_out.refimcube_compatibility EQ 0) AND (N_ELEMENTS(REFCUBE) NE 2)) THEN $
      ref_single_cube = hdr_out.single_cube[1] $
    ELSE IF (N_ELEMENTS(SINGLE_CUBE) EQ 2) THEN $
      ref_single_cube = single_cube[1]
  ; Handle reference image cube first, only after that check for reference spectral cube
    CRISPEX_IO_PARSE_HEADER, hdr_out.refimfilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                           CUBE_COMPATIBILITY=hdr_out.refimcube_compatibility, EXTEN_NO=0, /REFIMCUBE
    IF (N_ELEMENTS(REFCUBE) EQ 2) THEN BEGIN
      hdr_out.refspfilename = refcube[1]
     	refspext = STRMID(hdr_out.refspfilename,STRPOS(hdr_out.refspfilename,'.',/REVERSE_SEARCH)+1,$
                        STRLEN(hdr_out.refspfilename))
    	hdr_out.refspcube_compatibility = ABS(STRMATCH(refspext,'fits',/FOLD_CASE)-1)
      CRISPEX_IO_PARSE_HEADER, hdr_out.refspfilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                             CUBE_COMPATIBILITY=hdr_out.refspcube_compatibility, EXTEN_NO=0, /REFSPCUBE
    ENDIF 
;    ; If single_cube value has been set from single FITS cube, use that
;    IF ((hdr_out.single_cube[1] NE 0) AND (N_ELEMENTS(REFCUBE) NE 2)) THEN $
;      ref_single_cube = hdr_out.single_cube[1] $
;    ELSE IF (N_ELEMENTS(SINGLE_CUBE) EQ 2) THEN $
;      ref_single_cube = single_cube[1]
    CRISPEX_IO_FAILSAFES_REF, refcube, ref_single_cube, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                            STARTUPTLB=startuptlb, IO_FAILSAFE_ERROR=io_failsafe_ref_error
    IF (io_failsafe_ref_error EQ 1) THEN RETURN
    CRISPEX_IO_FAILSAFES_MAIN_REF, HDR=hdr_out, STARTUPTLB=startuptlb, $
                                 IO_FAILSAFE_ERROR=io_failsafe_main_ref_error
    IF (io_failsafe_main_ref_error EQ 1) THEN RETURN
    ; Add reference axes titles based on REFIMCUBE and REFSPCUBE headers
    IF (STRCOMPRESS(hdr_out.refbunit,/REMOVE_ALL) NE '') THEN $
      refytitle_unit = ' ['+hdr_out.refbunit+']' ELSE refytitle_unit = ''
    IF (STRCOMPRESS(hdr_out.reflpunit,/REMOVE_ALL) NE '') THEN $
      refxtitle_unit = ' ['+hdr_out.reflpunit+']' ELSE refxtitle_unit = ''
    hdr_out.ytitle[1] = hdr_out.refblabel+refytitle_unit
    hdr_out.xtitle[1] = hdr_out.reflplabel+refxtitle_unit
  ENDIF
;  IF (N_ELEMENTS(event) EQ 1) THEN BEGIN
;    (*(*info).ioparams).hdr = hdr_out
;    CRISPEX_IO_OPEN_REFCUBE_READ, event, REFCUBE=refcube
;    CRISPEX_IO_HANDLE_HDR, event, hdr_out, /REFERENCE
;  ENDIF ELSE 
  IF (hdr_out.refnt GT 1) THEN BEGIN
    tsel_ref = LONARR(hdr_out.mainnt)
    FOR tt=0,hdr_out.mainnt-1 DO BEGIN
      tdiff = ABS(hdr_out.tarr_ref - hdr_out.tarr_main[tt])
      tsel_ref[tt] = (WHERE(tdiff EQ MIN(tdiff, /NAN)))[0]
    ENDFOR
    hdr_out = CREATE_STRUCT(hdr_out, 'tsel_ref', tsel_ref)
  ENDIF ELSE $
    hdr_out = CREATE_STRUCT(hdr_out, 'tsel_ref', 0)
  CRISPEX_IO_OPEN_REFCUBE_READ, REFCUBE=refcube, HDR_IN=hdr_out, HDR_OUT=hdr_out
END

PRO CRISPEX_IO_OPEN_REFCUBE_READ, event, REFCUBE=refcube, HDR_IN=hdr_in, HDR_OUT=hdr_out
  IF (N_ELEMENTS(event) EQ 1) THEN BEGIN
    WIDGET_CONTROL, event.TOP, GET_UVALUE = info
    hdr_out = (*(*info).ioparams).hdr
  ENDIF ELSE hdr_out = hdr_in
  ; Read in reference image cube
	IF ((N_ELEMENTS(REFCUBE) GE 1) AND (SIZE(REFCUBE,/TYPE) EQ 7)) THEN BEGIN	; REFIMCUBE as filename
		OPENR, lunrefim, refcube[0], /get_lun, $
           SWAP_ENDIAN = ((hdr_out.refimtype GT 1) AND (hdr_out.endian NE hdr_out.refimendian))
    referencefile = ASSOC(lunrefim,MAKE_ARRAY(hdr_out.refnx,hdr_out.refny, TYPE=hdr_out.refimtype),$
                      hdr_out.refimoffset)
		hdr_out.showref = 1
    IF ((hdr_out.refnlp GT 1) AND (N_ELEMENTS(REFCUBE) NE 2)) THEN BEGIN
      refscanfile = ASSOC(lunrefim,MAKE_ARRAY(hdr_out.refnx,hdr_out.refny,$
                      hdr_out.refnlp*hdr_out.refns, TYPE=hdr_out.refimtype),hdr_out.refimoffset)
    ENDIF 
    hdr_out.lunrefim = lunrefim
    CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, REFIMCUBE=hdr_out.refimfilename
    ; Read in reference spectral cube if so provided
    IF (N_ELEMENTS(REFCUBE) EQ 2) THEN BEGIN                                ; REFSPCUBE as filename
			OPENR, lunrefsp, refcube[1], /get_lun, $
             SWAP_ENDIAN = ((hdr_out.refsptype GT 1) AND (hdr_out.endian NE hdr_out.refspendian))						
      referencespectra = ASSOC(lunrefsp,MAKE_ARRAY(hdr_out.refnlp,hdr_out.refnt, $
                            TYPE=hdr_out.refsptype),hdr_out.refspoffset)
			hdr_out.refspfile = 1	
      hdr_out.lunrefsp = lunrefsp
      CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, REFSPCUBE=hdr_out.refspfilename
    ENDIF 
	ENDIF ELSE BEGIN
    IF (N_ELEMENTS(REFCUBE) GT 1) THEN BEGIN								                ; REFCUBE as image array
		  hdr_out.refnx = LONG((SIZE(refcube))[1])
		  hdr_out.refny = LONG((SIZE(refcube))[2])
		  IF (SIZE(REFCUBE,/N_DIMENSIONS) EQ 3) THEN $
        hdr_out.refnt = LONG((SIZE(refcube))[3]) $
      ELSE $
        hdr_out.refnt = 1L
		  IF (hdr_out.refnx EQ hdr_out.nx) AND (hdr_out.refny EQ hdr_out.ny) THEN BEGIN
			  referencefile = refcube
			  hdr_out.showref = 1
			  hdr_out.refnlp = 1
			  hdr_out.refns = 1
		  ENDIF ELSE BEGIN
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Dimensions of the reference cube (['+$
          STRTRIM(hdr_out.refnx,2)+','+STRTRIM(hdr_out.refny,2)+','+STRTRIM(hdr_out.refnt,2)+']) are not '+$
          'compatible with those of the main image cube (['+STRTRIM(hdr_out.nx,2)+','+$
          STRTRIM(hdr_out.ny,2)+','+STRTRIM(hdr_out.mainnt,2)+'])! Number of pixels in the x- and '+$
          'y-direction must be the same for both.', /ERROR, /NO_ROUTINE, /NEWLINE
			WIDGET_CONTROL, startuptlb, /DESTROY
			RETURN
		  ENDELSE
	  ENDIF 
	ENDELSE
  ; Create pointers to data
	IF (hdr_out.showref EQ 1) THEN BEGIN
		hdr_out.refdata	= PTR_NEW(referencefile, /NO_COPY)
		hdr_out.refslice= PTR_NEW(BYTARR(hdr_out.nx,hdr_out.ny))
		IF ((hdr_out.refnlp GT 1) AND (hdr_out.refspfile EQ 0)) THEN BEGIN
			hdr_out.refscan = PTR_NEW(refscanfile, /NO_COPY)
      hdr_out.refsspscan = PTR_NEW(MAKE_ARRAY(hdr_out.nx,hdr_out.ny,hdr_out.refnlp*hdr_out.refns, $
                            TYPE=hdr_out.refimtype))
		ENDIF 
	  IF (hdr_out.refspfile EQ 1) THEN hdr_out.refspdata = PTR_NEW(referencespectra, /NO_COPY) 
	ENDIF 
END

PRO CRISPEX_IO_OPEN_SJICUBE, SJICUBE=sjicube, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                              STARTUPTLB=startuptlb, $
                              IO_FAILSAFE_SJI_ERROR=io_failsafe_sji_error
  io_failsafe_sji_error = 0
  io_failsafe_main_ref_error = 0
  hdr_out = hdr_in
	IF ((N_ELEMENTS(SJICUBE) GE 1) AND (SIZE(SJICUBE,/TYPE) EQ 7)) THEN BEGIN					
    hdr_out.sjifilename = sjicube[0]
  	sjiext = STRMID(hdr_out.sjifilename,STRPOS(hdr_out.sjifilename,'.',/REVERSE_SEARCH)+1,$
                      STRLEN(hdr_out.sjifilename))
  	sjicube_compatibility = ABS(STRMATCH(sjiext,'fits',/FOLD_CASE)-1)
    IF sjicube_compatibility THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'The slit-jaw image cube must be in FITS format.', $
        /ERROR, /NO_ROUTINE, /NEWLINE
  		IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
      io_failsafe_sji_error = 1
    ENDIF ELSE BEGIN
      ; Parse the SJICUBE header
      CRISPEX_IO_PARSE_HEADER, hdr_out.sjifilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                             CUBE_COMPATIBILITY=sjicube_compatibility, EXTEN_NO=0, /SJICUBE
      ; Check whether dimensions are compatible with main cube
      CRISPEX_IO_FAILSAFES_MAIN_SJI, HDR=hdr_out, STARTUPTLB=startuptlb, $
                                   IO_FAILSAFE_ERROR=io_failsafe_sji_error
    ENDELSE
    IF (io_failsafe_sji_error EQ 1) THEN RETURN
    IF (hdr_out.sjint GT 1) THEN BEGIN
      tsel_sji = LONARR(hdr_out.mainnt)
      FOR tt=0,hdr_out.mainnt-1 DO BEGIN
        tdiff = ABS(hdr_out.tarr_sji - hdr_out.tarr_main[tt])
        tsel_sji[tt] = (WHERE(tdiff EQ MIN(tdiff, /NAN)))[0]
      ENDFOR
      hdr_out = CREATE_STRUCT(hdr_out, 'tsel_sji', tsel_sji)
    ENDIF ELSE $
      hdr_out = CREATE_STRUCT(hdr_out, 'tsel_sji', 0)
    CRISPEX_IO_OPEN_SJICUBE_READ, HDR_IN=hdr_out, HDR_OUT=hdr_out
  ENDIF
END

PRO CRISPEX_IO_OPEN_SJICUBE_READ, HDR_IN=hdr_in, HDR_OUT=hdr_out
  hdr_out = hdr_in
		OPENR, lunsji, hdr_out.sjifilename, /get_lun, $
           SWAP_ENDIAN = ((hdr_out.sjitype GT 1) AND (hdr_out.endian NE hdr_out.sjiendian))
    sji = ASSOC(lunsji,MAKE_ARRAY(hdr_out.sjinx,hdr_out.sjiny, TYPE=hdr_out.sjitype),$
             hdr_out.sjioffset)
		hdr_out.sjifile = 1	
	  hdr_out.sjidata = PTR_NEW(sji, /NO_COPY)
    hdr_out.lunsji = lunsji
	  hdr_out.sjislice	= PTR_NEW(BYTARR(hdr_out.sjinx,hdr_out.sjiny))
    nrasters = N_ELEMENTS(hdr_out.xyrastersji[*,0])
    rastercont = FLTARR(nrasters,4)
    raster_width = hdr_out.dx / FLOAT(hdr_out.sjidx)
    raster_height = hdr_out.ny * (hdr_out.dy / FLOAT(hdr_out.sjidy))
    FOR i=0,nrasters-1 DO BEGIN
      rastercont[i,0] = hdr_out.xyrastersji[i,0]
      rastercont[i,1] = hdr_out.xyrastersji[i,1]
      rastercont[i,2] = hdr_out.xyrastersji[i,0]+raster_width
      rastercont[i,3] = hdr_out.xyrastersji[i,1]+raster_height
    ENDFOR
    hdr_out = CREATE_STRUCT(hdr_out, 'rastercont', rastercont)
    CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, SJICUBE=hdr_out.sjifilename
END

PRO CRISPEX_IO_OPEN_MASKCUBE, MASKCUBE=maskcube, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                              STARTUPTLB=startuptlb, IO_FAILSAFE_MASK_ERROR=io_failsafe_mask_error
  io_failsafe_mask_error = 0
  hdr_out = hdr_in
  IF (N_ELEMENTS(MASKCUBE) GE 1) THEN BEGIN
    hdr_out.maskfilename = maskcube[0]
  	maskext = STRMID(hdr_out.maskfilename,STRPOS(hdr_out.maskfilename,'.',/REVERSE_SEARCH)+1,$
                      STRLEN(hdr_out.maskfilename))
  	hdr_out.maskcube_compatibility = ABS(STRMATCH(maskext,'fits',/FOLD_CASE)-1)
    CRISPEX_IO_PARSE_HEADER, hdr_out.maskfilename, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                            CUBE_COMPATIBILITY=hdr_out.maskcube_compatibility, EXTEN_NO=0, /MASKCUBE
    CRISPEX_IO_FAILSAFES_MASK, HDR=hdr_out, STARTUPTLB=startuptlb, $
                                IO_FAILSAFE_ERROR=io_failsafe_mask_error
    IF (io_failsafe_mask_error EQ 1) THEN RETURN
    CRISPEX_IO_OPEN_MASKCUBE_READ, HDR_IN=hdr_out, HDR_OUT=hdr_out
  ENDIF 
END

PRO CRISPEX_IO_OPEN_MASKCUBE_READ, HDR_IN=hdr_in, HDR_OUT=hdr_out
  hdr_out = hdr_in
		OPENR, lunmask, hdr_out.maskfilename, /get_lun, $
           SWAP_ENDIAN = ((hdr_out.masktype GT 1) AND (hdr_out.endian NE hdr_out.maskendian))
    mask = ASSOC(lunmask,MAKE_ARRAY(hdr_out.masknx,hdr_out.maskny, TYPE=hdr_out.masktype),$
             hdr_out.maskoffset)
		hdr_out.maskfile = 1	
	  hdr_out.maskdata = PTR_NEW(mask, /NO_COPY)
    hdr_out.lunmask = lunmask
    CRISPEX_IO_FEEDBACK, hdr_out.verbosity, hdr_out, MASKCUBE=hdr_out.maskfilename
END

PRO CRISPEX_IO_HANDLE_HDR, event, hdr, MAIN=main, REFERENCE=reference, MASK=mask
  WIDGET_CONTROL, event.TOP, GET_UVALUE = info
  IF KEYWORD_SET(REFERENCE) THEN BEGIN
    (*(*info).data).refdata = hdr.refdata
    (*(*info).data).refslice = hdr.refslice
    (*(*info).data).refspdata = hdr.refspdata
    (*(*info).data).refscan = hdr.refscan
    (*(*info).data).refsspscan = hdr.refsspscan
    (*(*info).data).lunrefim = hdr.lunrefim
    (*(*info).data).lunrefsp = hdr.lunrefsp
  
;   (*(*info).dataparams).refimfilename = hdr.refcube
;   (*(*info).dataparams).refspfilename = hdr.refspcube
    (*(*info).dataparams).refnlp = hdr.refnlp
    (*(*info).dataparams).reflps = hdr.reflps
    (*(*info).dataparams).refms = hdr.refms
    (*(*info).dataparams).refspec = hdr.refspect
    (*(*info).dataparams).refnt = hdr.refnt
    (*(*info).dataparams).bunit[1] = hdr.refbunit
    (*(*info).dataparams).lpunit[1] = hdr.reflpunit
    
    (*(*info).dataswitch).reffile = hdr.showref
    (*(*info).dataswitch).refspfile = hdr.refspfile
    
    (*(*info).dispparams).lp_ref_upp = hdr.refnlp-1
    (*(*info).dispparams).lp_ref_range = hdr.refnlp
  
    (*(*info).overlayswitch).maskim[1] = hdr.showref
  
    (*(*info).plotaxes).v_dop_ref = hdr.v_dop_ref
  
    (*(*info).winswitch).showrefsp = hdr.refspfile
    (*(*info).winswitch).showref = hdr.showref
  ENDIF
END

PRO CRISPEX_IO_SETTINGS_SPECTRAL, event, HDR_IN=hdr_in, HDR_OUT=hdr_out, MNSPEC=mnspec, $
                                  SPECTFILE=spectfile, LINE_CENTER=line_center, NO_WARP=no_warp, $
                                  STARTUPTLB=startuptlb, $
                                  IO_FAILSAFE_MNSPEC_ERROR=io_failsafe_mnspec_error, $
                                  IO_FAILSAFE_IMSPECTFILE_ERROR=io_failsafe_imspectfile_error,$
                                  IO_FAILSAFE_REFSPECTFILE_ERROR=io_failsafe_refspectfile_error, $
                                  IO_FAILSAFE_LINE_CENTER_ERROR=io_failsafe_line_center_error
  ; Handles spectral settings on loading cubes
  IF (N_ELEMENTS(event) EQ 1) THEN $
    WIDGET_CONTROL, event.TOP, GET_UVALUE = info $
  ELSE $
    hdr_out = hdr_in

  ; Check for correct MNSPEC setting
  io_failsafe_mnspec_error = 0
  IF (N_ELEMENTS(MNSPEC) NE 0) THEN BEGIN
    CRISPEX_IO_FAILSAFES_MNSPEC, mnspec, hdr_out, STARTUPTLB=startuptlb, $
                                 IO_FAILSAFE_ERROR=io_failsafe_mnspec_error
    IF io_failsafe_mnspec_error THEN RETURN
  ENDIF 

  ; Handle SPECTFILE input; will only be done when cube is read in compatibility mode (i.e., for
  ; non-FITS cubes). Check whether SPECTFILE has been provided with save files
  io_failsafe_imspectfile_error = 0
  io_failsafe_refspectfile_error = 0
  IF ((N_ELEMENTS(SPECTFILE) GE 1) AND (SIZE(SPECTFILE,/TYPE) EQ 7)) THEN $
    spectfile_set = (spectfile[0] NE '') ELSE BEGIN
    spectfile_set = 0
  ENDELSE
  CRISPEX_IO_PARSE_SPECTFILE, spectfile, *hdr_out.imdata, hdr_out.verbosity, HDR_IN=hdr_out, HDR_OUT=hdr_out, $
                              MNSPEC=mnspec, /IMCUBE, CUBE_COMPATIBILITY=hdr_out.imcube_compatibility, $
                              STARTUPTLB=startuptlb, IO_FAILSAFE_ERROR=io_failsafe_imspectfile_error
  IF (io_failsafe_imspectfile_error EQ 1) THEN RETURN 
  refspectfile_set = ((N_ELEMENTS(SPECTFILE) EQ 2) AND (SIZE(SPECTFILE,/TYPE) EQ 7)) 
  IF (STRCOMPRESS(hdr_out.refimfilename) NE '') THEN BEGIN
    IF ((N_ELEMENTS(SPECTFILE) EQ 1) AND (refspectfile_set EQ 0)) THEN spectfile = [spectfile, '']
    CRISPEX_IO_PARSE_SPECTFILE, spectfile, *hdr_out.refdata, hdr_out.verbosity, HDR_IN=hdr_out, $
                                HDR_OUT=hdr_out, $
                                MNSPEC=mnspec, /REFCUBE, STARTUPTLB=startuptlb, $
                                CUBE_COMPATIBILITY=hdr_out.refimcube_compatibility, $
                                IO_FAILSAFE_ERROR=io_failsafe_refspectfile_error
    IF io_failsafe_refspectfile_error THEN RETURN
  ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, 'refspec', 0, 'refms', 0, 'reflps', 0)

  ; Handle LINE_CENTER input
  io_failsafe_line_center_error = 0
  IF (N_ELEMENTS(LINE_CENTER) NE 0) THEN BEGIN
    CRISPEX_IO_FAILSAFES_LINE_CENTER, line_center, hdr_out, NFILES=(hdr_out.showref+1), $
                                      STARTUPTLB=startuptlb, SPECTFILE_SET=spectfile_set, $
                                      REFSPECTFILE_SET=refspectfile_set,$
                                      IO_FAILSAFE_ERROR=io_failsafe_line_center_error
  ENDIF ;ELSE io_failsafe_line_center_error = 0
  cube_compatibility = [hdr_out.imcube_compatibility,hdr_out.refimcube_compatibility]
  IF (io_failsafe_line_center_error EQ 1) THEN RETURN ELSE $
    CRISPEX_IO_PARSE_LINE_CENTER, line_center, NFILES=(hdr_out.refnlp GT 1)+1, HDR_IN=hdr_out, $
                                  HDR_OUT=hdr_out,$
                                  SPECTFILE_SET=spectfile_set, REFSPECTFILE_SET=refspectfile_set,$
                                  CUBE_COMPATIBILITY=cube_compatibility, DLAMBDA_VAL=dlambda, $
                                  DLAMBDA_SET=dlambda_set, V_DOP_SET=v_dop_set, $
                                  IO_FAILSAFE_ERROR=io_failsafe_line_center_error
  hdr_out = CREATE_STRUCT(hdr_out, 'dlambda', dlambda, 'dlambda_set', dlambda_set, $
                          'v_dop_set', v_dop_set)
  
  ; Process settings from LINE_CENTER parsing 
  IF (hdr_out.refnlp GT 1) THEN BEGIN
    IF hdr_out.refimcube_compatibility THEN $
      hdr_out.refspxtitle = (['Spectral position','Wavelength'])[hdr_out.v_dop_set[1]] $
    ELSE $
      hdr_out.refspxtitle = hdr_out.xtitle[1]
  ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, 'reflc', 0L, 'v_dop_ref', 0)
;  ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, 'v_dop_ref', 0)
  IF hdr_out.imcube_compatibility THEN $
    hdr_out.spxtitle = (['Spectral position','Wavelength'])[hdr_out.v_dop_set[0]] $
  ELSE $
    hdr_out.spxtitle = hdr_out.xtitle[0]
  
  ; Process settings based on LPS variable and NO_WARP keyword to (not) warp of spectral slices
  IF (hdr_out.ndiagnostics GT 1) THEN BEGIN
    idx = LINDGEN(hdr_out.nlp)
    FOR d=0,hdr_out.ndiagnostics-1 DO BEGIN
      IF (d EQ 0) THEN $
        idx_select = idx[hdr_out.diag_start[d]:(hdr_out.diag_start[d]+hdr_out.diag_width[d]-2)] $
      ELSE $
        idx_select = [idx_select,idx[hdr_out.diag_start[d]:(hdr_out.diag_start[d]+$
          hdr_out.diag_width[d]-2)]]
    ENDFOR
    idx_select = [idx_select,idx(hdr_out.diag_start[d-1]+hdr_out.diag_width[d-1]-1)]
  ENDIF
  IF ~KEYWORD_SET(NO_WARP) THEN no_warp = (hdr_out.ndiagnostics GT 1)
  CRISPEX_IO_PARSE_WARPSLICE, hdr_out.lps, hdr_out.nlp, hdr_out.mainnt, hdr_out.dlambda[0], $
                              hdr_out.dlambda_set[0], hdr_out.verbosity, NO_WARP=no_warp, $
                              WARPSPSLICE=warpspslice, XO=xo, XI=xi, YO=yo, YI=yi, $
                              IDX_SELECT=idx_select;, DIAG_START=hdr_out.diag_start, $
;                              DIAG_WIDTH=hdr_out.diag_width
  IF (hdr_out.nrefdiagnostics GT 1) THEN BEGIN
    refidx = LINDGEN(hdr_out.refnlp)
    FOR d=0,hdr_out.nrefdiagnostics-1 DO BEGIN
      IF (d EQ 0) THEN $
        refidx_select = refidx[hdr_out.refdiag_start[d]:(hdr_out.refdiag_start[d]+hdr_out.refdiag_width[d]-2)] $
      ELSE $
        refidx_select = [refidx_select,refidx[hdr_out.refdiag_start[d]:(hdr_out.refdiag_start[d]+$
          hdr_out.refdiag_width[d]-2)]]
    ENDFOR
    refidx_select = [refidx_select,refidx(hdr_out.refdiag_start[d-1]+hdr_out.refdiag_width[d-1]-1)]
  ENDIF
  IF ~KEYWORD_SET(NO_WARP) THEN no_warp = (hdr_out.nrefdiagnostics GT 1)
  CRISPEX_IO_PARSE_WARPSLICE, hdr_out.reflps, hdr_out.refnlp, hdr_out.refnt, hdr_out.dlambda[1], $   
                              hdr_out.dlambda_set[1], hdr_out.verbosity, NO_WARP=no_warp, $
                              WARPSPSLICE=warprefspslice, XO=xo_ref, XI=xi_ref, YO=yo_ref, $
                              YI=yi_ref, IDX_SELECT=refidx_select, /REFERENCE;, $
;                              DIAG_START=hdr_out.refdiag_start, DIAG_WIDTH=hdr_out.refdiag_width
  hdr_out = CREATE_STRUCT(hdr_out, 'xo', xo, 'xi', xi, 'yo', yo, 'yi', yi, 'xo_ref', xo_ref, $
                                    'xi_ref', xi_ref, 'yo_ref', yo_ref, 'yi_ref', yi_ref, $
;                                    'phisxo',phisxo,'phisxi',phisxi,'phisyo',phisyo,'phisyi',phisyi,$
                                    'warpspslice', warpspslice, 'warprefspslice', warprefspslice)
END

PRO CRISPEX_IO_PARSE_LINE_CENTER, line_center, NFILES=nfiles, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                                  SPECTFILE_SET=spectfile_set, REFSPECTFILE_SET=refspectfile_set, $
                                  CUBE_COMPATIBILITY=cube_compatibility, DLAMBDA_VAL=dlambda_val, $
                                  DLAMBDA_SET=dlambda_set, V_DOP_SET=v_dop_set, $
                                  IO_FAILSAFE_ERROR=io_failsafe_error
; Handles parsing of LINE_CENTER keyword
  ; lcase of 0: LINE_CENTER = <Undefined>
  ; lcase of 1: LINE_CENTER = lc or [[lc],[lc]]
  ; lcase of 2: LINE_CENTER = [cwav,dwav] or [[cwav,dwav],[cwav,dwav]]
  ; lcase of 3: LINE_CENTER = [lc,cwav,dwav] or [[lc,cwav,dwav],[lc,cwav,dwav]]
  hdr_out = hdr_in	
  IF (N_ELEMENTS(LINE_CENTER) GT 0) THEN BEGIN
    ndims = SIZE(LINE_CENTER,/N_DIMENSIONS) > 1
    nelem = N_ELEMENTS(LINE_CENTER)
    IF (io_failsafe_error NE 2) THEN lcase = nelem / ndims ELSE lcase = 0
  ENDIF ELSE BEGIN
    ndims = 1
    lcase = 0
  ENDELSE
  c_speed	= 2.99792458D5													; Speed of light in km/s
;  lc = LONARR(nfiles)   &  
  dlambda_val = FLTARR(2)  &  lambda_c = FLTARR(2)
	v_dop_set = BYTARR(2) &  dlambda_set = BYTARR(2)
  nlp_select = [hdr_out.nlp,hdr_out.refnlp]
  spectfile_set = [spectfile_set,refspectfile_set]
  
  FOR d=0,nfiles-1 DO BEGIN
    IF (d EQ 0) THEN tag2 = 'lc' ELSE tag2 = 'reflc'
    IF ((cube_compatibility[d] NE 1) AND (N_ELEMENTS(LINE_CENTER) GT 0) AND (d NE ndims-1)) THEN BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Calling CRISPEX with LINE_CENTER, while FITS cubes '+$
        'are provided, is not allowed. Settings from LINE_CENTER will be ignored.', /WARNING, $
        /NO_ROUTINE
        lcase = 0
    ENDIF 
    IF (d EQ 0) THEN BEGIN
      lps_select = hdr_out.lps    &   spec_select = hdr_out.spectrum 
    ENDIF ELSE BEGIN
      IF ((nfiles EQ 2) AND (ndims NE 2)) THEN lcase = 0
      lps_select = hdr_out.reflps &   spec_select = hdr_out.refspec
    ENDELSE
    IF ((lcase EQ 1) OR (lcase EQ 3)) THEN $            ; First read parameters from LINE_CENTER
      lc = line_center[0,d] $
    ELSE BEGIN              ; Or determine from mean spectrum if line centre position is not supplied
      IF (WHERE(TAG_NAMES(hdr_out) EQ STRUPCASE(tag2)) EQ -1) THEN $          ; Check whether LPS exists
        lc = ( WHERE( spec_select EQ MIN(spec_select, /NAN) ) )[0] $
      ELSE BEGIN
        IF (d EQ 0) THEN lc = hdr_out.lc ELSE lc = hdr_out.reflc
      ENDELSE
    ENDELSE
    IF (lcase GE 2) THEN BEGIN
      lambda_c[d]    = line_center[(lcase EQ 3),d]
      dlambda_val[d] = line_center[1+(lcase EQ 3),d]
      dlambda_set[d] = 1
    ENDIF 
    ; Next determine whether LPS needs to be recalculated, only when in compatibility mode
    IF cube_compatibility[d] THEN BEGIN
      v_dop_set[d] = ((spectfile_set[d] EQ 1) OR (lcase GT 2))
      IF (v_dop_set[d] EQ 0) THEN lps_select -= lps_select[LONG(lc)]
      IF (lcase GE 2) THEN lps_select = lps_select*dlambda_val[d] + lambda_c[d]
      ndiagnostics = 1
      diag_width = nlp_select[d]
      diag_start = 0
    ENDIF ELSE BEGIN
      IF (d EQ 0) THEN BEGIN
        ndiagnostics = hdr_out.ndiagnostics 
        twave = hdr_out.twave
        diag_start = hdr_out.diag_start
        diag_width = hdr_out.diag_width
      ENDIF ELSE BEGIN
        ndiagnostics = hdr_out.nrefdiagnostics
        twave = hdr_out.twave_ref
        diag_start = hdr_out.refdiag_start
        diag_width = hdr_out.refdiag_width
      ENDELSE
      v_dop_set[d] = (lps_select[LONG(lc)] NE 0)
    ENDELSE
    v_dop = PTRARR(ndiagnostics)
    FOR dd=0,ndiagnostics-1 DO BEGIN
  	  IF v_dop_set[d] THEN BEGIN
        lps_select_loc = lps_select[diag_start[dd]:(diag_start[dd]+(diag_width[dd]-1))]
        IF cube_compatibility[d] THEN $
          v_dop_loc = c_speed*(lps_select_loc/lps_select_loc[LONG(lc)]-1) $
        ELSE BEGIN
          lc = twave[dd]
          v_dop_loc = c_speed*(lps_select_loc/FLOAT(lc)-1)
        ENDELSE
      ENDIF ELSE BEGIN
        nlp_select_loc = diag_width[dd]
        v_dop_loc = FLTARR(nlp_select_loc)					; array with Doppler velocities in km/s
      ENDELSE
      v_dop[dd] = PTR_NEW(v_dop_loc)
    ENDFOR
    IF (d EQ 0) THEN BEGIN      ; Fill hdr tags with main results
      hdr_out.lps = lps_select
      tag = 'v_dop'
    ENDIF ELSE BEGIN            ; Fill hdr tags with reference results
      hdr_out.reflps = lps_select   
      tag = 'v_dop_ref'
    ENDELSE
    IF (WHERE(TAG_NAMES(hdr_out) EQ STRUPCASE(tag2)) EQ -1) THEN BEGIN
      lc_sel = lc
      hdr_out = CREATE_STRUCT(hdr_out, tag, v_dop, tag2, lc_sel)
    ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, tag, v_dop)
  ENDFOR
END

PRO CRISPEX_IO_PARSE_SINGLE_CUBE, input_single_cube, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
  MAIN=main, REFERENCE=reference
  hdr_out = hdr_in
  ; Check setting of SINGLE_CUBE keyword
  IF KEYWORD_SET(MAIN) THEN BEGIN
  	IF (N_ELEMENTS(INPUT_SINGLE_CUBE) EQ 1) THEN BEGIN  ; If SINGLE_CUBE properly set, set nlp and nt
  		hdr_out.nlp = LONG(INPUT_SINGLE_CUBE)
  		hdr_out.onecube = 1
      hdr_out.single_cube[0] = hdr_out.nlp
  		hdr_out.mainnt = hdr_in.imnt / hdr_in.nlp / hdr_in.ns
  	ENDIF ELSE IF (hdr_out.spfile EQ 0) THEN BEGIN  
      ; If no SPCUBE or SINGLE_CUBE are set, we are dealing with a snapshot
  		hdr_out.nlp = hdr_in.imnt / hdr_in.ns
  		hdr_out.single_cube[0] = 0
    ENDIF
  ENDIF ELSE IF KEYWORD_SET(REFERENCE) THEN BEGIN
		IF (N_ELEMENTS(INPUT_SINGLE_CUBE) EQ 1) THEN BEGIN  ; If SINGLE_CUBE properly set, set refnlp 
			hdr_out.refnlp = LONG(INPUT_SINGLE_CUBE)
			hdr_out.onecube = 1
      hdr_out.single_cube[1] = hdr_out.refnlp
		ENDIF 
  ENDIF
END

PRO CRISPEX_IO_PARSE_SPECTFILE, spectfile, datafile, verbosity, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                                MNSPEC=mnspec, IMCUBE=imcube, REFCUBE=refcube, $
                                CUBE_COMPATIBILITY=cube_compatibility, STARTUPTLB=startuptlb, $
                                IO_FAILSAFE_ERROR=io_failsafe_error
; Handles parsing the spectral save files into the appropriate variables
  hdr_out = hdr_in
  io_failsafe_error = 0 &  calc_from_cubes = 0
  spectfile_set = 0     &  refspectfile_set = 0
  IF KEYWORD_SET(IMCUBE) THEN BEGIN
    nlp_select = hdr_out.nlp
    ns_select = hdr_out.ns
    feedback_text = 'main '
  ENDIF 
  IF KEYWORD_SET(REFCUBE) THEN BEGIN
    nlp_select = hdr_out.refnlp
    ns_select = hdr_out.refns
    feedback_text = 'reference '
  ENDIF
  IF (N_ELEMENTS(SPECTFILE) NE 0) THEN BEGIN
    IF KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN
      spectfile_set = KEYWORD_SET(IMCUBE)
      refspectfile_set = KEYWORD_SET(REFCUBE)
      IF (spectfile[refspectfile_set] NE '') THEN BEGIN ; Check whether SPECTFILE even has been given
  		  IF (N_ELEMENTS(MNSPEC) GT 0) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
          'WARNING: Calling CRISPEX with MNSPEC, while SPECTFILE is provided, is not allowed.'+$
          ' Using SPECTFILE for mean spectrum determination.', /NO_ROUTINE
          io_failsafe_error = 2
        RESTORE, spectfile[refspectfile_set];,VERBOSE=(TOTAL(verbosity[0:1]) GE 1) 
  			IF (verbosity[1] EQ 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Restored '+feedback_text+$
                                           'spectral file: '+spectfile[refspectfile_set]+'.'
        IF (N_ELEMENTS(norm_factor) NE 1) THEN BEGIN   ; Failsafe against old spectfiles; convert vars
          spect_pos = ll
          norm_factor = 1/mn
          norm_spect = spec
  				CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'WARNING: The restored '+feedback_text+$
            'spectral file ('+spectfile[refspectfile_set]+$
            ') has a format pre-dating CRISPEX v1.6. If possible, please update the '+$
            feedback_text+'spectral file to the newest format. Reading spectral file in '+$
            'compatibility mode.', /NO_ROUTINE, NEWLINE=(io_failsafe_error NE 2)
          io_failsafe_error = 2
        ENDIF ELSE BEGIN                       ; Else, dealing with new spectfile, check title labels
  				IF (N_ELEMENTS(xtitle_label) EQ 1) THEN BEGIN                     ; If xtitle_label exists
  					IF ((STRCOMPRESS(xtitle_label) NE '') AND $                     ; and it's non-zero
                (STRCOMPRESS(hdr_out.xtitle[refspectfile_set]) EQ '')) THEN $ ; and hdr.xtitle empty 
              hdr_out.xtitle[refspectfile_set] = xtitle_label
  				ENDIF
  				IF (N_ELEMENTS(ytitle_label) EQ 1) THEN BEGIN                     ; If ytitle_label exists
  					IF ((STRCOMPRESS(ytitle_label) NE '') AND $                     ; and it's non-zero
                (STRCOMPRESS(hdr_out.ytitle[refspectfile_set]) EQ '')) THEN $ ; and hdr.ytitle empty
              hdr_out.ytitle[refspectfile_set] = ytitle_label
  				ENDIF
        ENDELSE
        ; Failsafe against data-incompatible SPECTFILE
        IF (N_ELEMENTS(spect_pos[*,0]) NE nlp_select) THEN BEGIN  
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'ERROR: Number of spectral positions in '+$
            feedback_text+'SPECTFILE (nlp='+STRTRIM(N_ELEMENTS(spect_pos[*,0]),2)+$
            ') is incompatible with that in the '+feedback_text+'datacubes (nlp='+$
            STRTRIM(nlp_select,2)+'). Please load the correct '+feedback_text+'spectral file.', $
            /NO_ROUTINE, NEWLINE=(io_failsafe_error NE 2)
  					IF (N_ELEMENTS(STARTUPTLB) EQ 1) THEN WIDGET_CONTROL, startuptlb, /DESTROY
            io_failsafe_error = 1
  					RETURN
        ENDIF
        IF KEYWORD_SET(IMCUBE) THEN tag = 'lps'
        IF KEYWORD_SET(REFCUBE) THEN tag = 'reflps'
        hdr_out = CREATE_STRUCT(hdr_out, tag, spect_pos)
      ENDIF ELSE calc_from_cubes = 1
    ENDIF ELSE BEGIN
      IF ((KEYWORD_SET(IMCUBE) AND (N_ELEMENTS(SPECTFILE) EQ 1)) OR $
          (KEYWORD_SET(REFCUBE) AND (N_ELEMENTS(SPECTFILE) EQ 2))) THEN $
          CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'WARNING: Calling CRISPEX with SPECTFILE, while'+$
            ' FITS cubes are provided, is not allowed. Settings from SPECTFILE will be ignored.', $
            /NO_ROUTINE, NEWLINE=(io_failsafe_error NE 2)
        io_failsafe_error = 2
      calc_from_cubes = 1
    ENDELSE
  ENDIF ELSE calc_from_cubes = 1
  IF calc_from_cubes THEN BEGIN
    IF (N_ELEMENTS(MNSPEC) GE 1) THEN BEGIN                       ; If MNSPEC set
      t_lower = mnspec[0]                                         ; Set lower t-value to mnspec[0]
      IF (N_ELEMENTS(MNSPEC) EQ 1) THEN t_upper = mnspec[0] ELSE t_upper = mnspec[1]
    ENDIF ELSE BEGIN                                              ; If not MNSPEC set, default to 0
      t_lower = 0 & t_upper = 0
    ENDELSE
    spectrum = DBLARR(nlp_select,ns_select)
    norm_spect = DBLARR(nlp_select,ns_select)
    IF ((ns_select GT 1) AND (N_ELEMENTS(scale_stokes) EQ 1)) THEN $  ; If ns>1 and scale_stokes set
      norm_factor = REPLICATE(scale_stokes,ns_select) $               ; fill norm_factor with those
    ELSE $
      norm_factor = FLTARR(ns_select)
    FOR s=0,ns_select-1 DO BEGIN        ; Loop over all "Stokes" positions
      FOR lp=0,nlp_select-1 DO BEGIN    ; Loop over all "spectral" positions
        spect = 0.                      ; (Re-)define spect variable
        FOR t=t_lower,t_upper DO $      ; Loop over all time steps
          spect += MEAN(datafile[t*nlp_select*ns_select + s*nlp_select + lp], /NAN)
        spectrum[lp,s] = spect / FLOAT(t_upper - t_lower + 1)
      ENDFOR
      IF (ns_select GT 1) THEN BEGIN                                  
        IF (N_ELEMENTS(scale_stokes) NE 1) THEN BEGIN ; If ns>1 & not scale_stokes, find norm_factor
          min_spectrum = MIN(spectrum[*,s], MAX=max_spectrum, /NAN)
          norm_factor[s] = (ABS([min_spectrum,max_spectrum]))[ABS(max_spectrum) GE ABS(min_spectrum)]
        ENDIF
        norm_spect[*,s] = spectrum[*,s] / norm_factor[s]              ; Determine norm_spect
      ENDIF 
    ENDFOR
    IF (ns_select EQ 1) THEN BEGIN      ; If ns=1, norm_factor and norm_spect follow easily
      norm_factor = MAX(spectrum, /NAN)
      norm_spect = spectrum / FLOAT(norm_factor)
    ENDIF
  ENDIF
  IF KEYWORD_SET(IMCUBE) THEN BEGIN     ; Fill the related IMCUBE variables with results
    reform_norm_spect = REFORM(norm_spect[*,0])
    hdr_out = CREATE_STRUCT(hdr_out, 'mainspec', norm_spect, 'spectrum', reform_norm_spect, $
                                     'ms', norm_factor)
    IF (WHERE(TAG_NAMES(hdr_out) EQ 'LPS') NE -1) THEN BEGIN          ; Check whether LPS exists
      IF (TOTAL(hdr_out.lps) EQ 0) THEN hdr_out.lps = FINDGEN(hdr_out.nlp)
    ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, 'lps', FINDGEN(hdr_out.nlp))
  ENDIF
  IF KEYWORD_SET(REFCUBE) THEN BEGIN    ; Fill the related REFCUBE variables with results
    hdr_out = CREATE_STRUCT(hdr_out, 'refspec', norm_spect, 'refms', norm_factor)
    IF (WHERE(TAG_NAMES(hdr_out) EQ 'REFLPS') NE -1) THEN BEGIN       ; Check whether REFLPS exists
      IF (TOTAL(hdr_out.reflps) EQ 0) THEN hdr_out.reflps = FINDGEN(hdr_out.refnlp)
    ENDIF ELSE hdr_out = CREATE_STRUCT(hdr_out, 'reflps', FINDGEN(hdr_out.refnlp))
  ENDIF
END

PRO CRISPEX_IO_PARSE_WARPSLICE, lps, nlp, nt, dlambda, dlambda_set, verbosity, NO_WARP=no_warp, $
                                WARPSPSLICE=warpspslice, XO=xo, XI=xi, YO=yo, YI=yi, $
                                REFERENCE=reference, IDX_SELECT=idx_select;, $
;                                DIAG_START=diag_start, DIAG_WIDTH=diag_width
; Handles warping of the slice, also given keyword NO_WARP
  warpspslice = 0				; Temporal spectrum is not warped to correct non-equidistant spectral positions
  xi = 0  &   yi = 0  &  xo = 0  &  yo = 0
  feedback_text = ['main','reference']  &  warp_text = ['No warp','Warp']
	IF (nlp GT 1) THEN BEGIN
		IF dlambda_set THEN ndecimals = ABS(FLOOR(ALOG10(ABS(dlambda)))) ELSE ndecimals = 2
		equidist = STRING((SHIFT(FLOAT(lps),-1) - FLOAT(lps))[0:nlp-2],$
                        FORMAT='(F8.'+STRTRIM(ndecimals,2)+')')
    ; if ndiagnostics > 1, check for equidistancy diagnostic excluding jumps
;    IF (N_ELEMENTS(DIAG_START) LT 1) THEN BEGIN
;      diag_start = 0
;      diag_width = nlp
;    ENDIF ;ELSE equidist[diag_start-1] = 0
    IF (N_ELEMENTS(IDX_SELECT) GT 0) THEN equidist = equidist[idx_select] 
    ; Check for non-equidistant spectral positions and allowed consequential warping
		warpspslice = (((WHERE(equidist NE equidist[0]))[0] NE -1) AND ~KEYWORD_SET(NO_WARP))
    IF warpspslice THEN BEGIN
		  min_lps = MIN(lps, /NAN)
		  xi = FINDGEN(nlp) # REPLICATE(1,nt)
		  xo = ((lps-min_lps) / FLOAT(MAX(lps-min_lps, /NAN)) * (nlp-1)) # REPLICATE(1,nt)
		  yi = REPLICATE(1,nlp) # FINDGEN(nt)
		  yo = yi
;		  xi = FINDGEN(nlp) # REPLICATE(1,nt)
;		  yi = REPLICATE(1,nlp) # FINDGEN(nt)
;      xo = FLTARR(nlp,nt)
;      lps_rewrite = FLTARR(nlp)
;      lps_rewrite[0] = equidist[0]
;      FOR i=0,N_ELEMENTS(equidist)-1 DO lps_rewrite[i] = lps_rewrite[i-1] + equidist[i]
;      min_lps = MIN(lps_rewrite)
;		  xo = ((lps_rewrite-min_lps) / FLOAT(MAX(lps_rewrite-min_lps)) * nlp) # REPLICATE(1,nt)
;;      FOR i=0,N_ELEMENTS(DIAG_START)-1 DO BEGIN
;;        lps_sel = lps[diag_start[i]:(diag_start[i]+diag_width[i]-1)]
;;        min_lps = MIN(lps_sel)
;;  		  xo[diag_start[i]:(diag_start[i]+diag_width[i]-1),*] = $
;;          ((lps_sel-min_lps) / FLOAT(MAX(lps_sel-min_lps)) * diag_width[i]) # REPLICATE(1,nt)
;;        stop
;;      ENDFOR
;		  yo = yi
		ENDIF 
    IF verbosity[1] THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, warp_text[warpspslice]+$
                          ' applied to '+feedback_text[KEYWORD_SET(REFERENCE)]+$
                          ' spectral slice.', /NO_ROUTINE
  ENDIF 
END

;================================================================================= LOOP PROCEDURES
PRO CRISPEX_LOOP_DEFINE, event
; Handles the start of loop definition procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	(*(*info).overlayswitch).looppath_feedback = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayswitch).looppath_feedback], labels=['Loop path feedback']
END

PRO CRISPEX_LOOP_GET, event
; Gets the loop path and (if the loop display window is open) also the loopslab for display
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_LOOP_GET_PATH, event
	IF (*(*info).winswitch).showloop THEN BEGIN
		CRISPEX_LOOP_GET_SLAB, event
		CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event					
	ENDIF
END

PRO CRISPEX_LOOP_GET_PATH, event
; Gets the actual loop path from spline interpolation
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	np_local = (SIZE(*(*(*info).loopparams).xp))[1] 
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF KEYWORD_SET(IM) THEN BEGIN
    tlow = (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).tlow]
    tupp = (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).tupp]
  ENDIF ELSE BEGIN
    tlow = (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).tlow]
    tupp = (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).tupp]
  ENDELSE
	IF KEYWORD_SET(NO_NLP) THEN BEGIN
;		FOR t=(*(*info).dispparams).t_low,(*(*info).dispparams).t_upp DO BEGIN    
		FOR tt=tlow,tupp DO BEGIN
			IF (tt EQ tlow) THEN $
        tmp = INTERPOLATE( extractdata[tt], (xrs)[w_lpts],(yrs)[w_lpts]) $
      ELSE $
        tmp = [[tmp], [INTERPOLATE( extractdata[tt], (xrs)[w_lpts],(yrs)[w_lpts])]]
		ENDFOR
	ENDIF ELSE BEGIN
		IF KEYWORD_SET(IM) THEN BEGIN		; Get space-time diagram from main cube
			FOR tt=tlow,tupp DO BEGIN
				IF (tt EQ tlow) THEN $
          tmp = INTERPOLATE( extractdata[$
            tt * (*(*info).dataparams).nlp * (*(*info).dataparams).ns+$
            (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp], $
            (xrs)[w_lpts],(yrs)[w_lpts])  $
        ELSE $
					tmp = [[tmp], [INTERPOLATE( extractdata[$
            tt * (*(*info).dataparams).nlp * (*(*info).dataparams).ns + $
            (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp], $
            (xrs)[w_lpts],(yrs)[w_lpts])]] 
			ENDFOR
		ENDIF ELSE BEGIN			; Get space-time diagram from reference cube
			FOR tt=tlow,tupp DO BEGIN
				IF (tt EQ tlow) THEN $
          tmp = INTERPOLATE( extractdata[$
            tt * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref], $
            (xrs)[w_lpts],(yrs)[w_lpts]) $
        ELSE $
					tmp = [[tmp], [INTERPOLATE( extractdata[$
            tt * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref], $
            (xrs)[w_lpts],(yrs)[w_lpts])]]
			ENDFOR
		ENDELSE
	ENDELSE
	ex_loopslice = tmp
	n = 1 + (ABS(((SIZE(xps))[0] EQ 1)-1))
	FOR k=0,(SIZE(xps))[n]-1 DO BEGIN
		IF (k EQ 0) THEN $
      ex_crossloc = [0] $
    ELSE $
      ex_crossloc = [ex_crossloc,WHERE( (xrs EQ (xps)[k]) AND (yrs EQ (yps)[k]))]
	ENDFOR
	ex_loopsize = (SIZE(ex_loopslice))[1]
END

PRO CRISPEX_LOOP_REMOVE_POINT, event
; Handles the removal of a loop point
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).loopparams).np -= 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).loopparams).np], labels=['np']
	*(*(*info).loopparams).xp = (*(*(*info).loopparams).xp)[0:(*(*info).loopparams).np-1]
	*(*(*info).loopparams).yp = (*(*(*info).loopparams).yp)[0:(*(*info).loopparams).np-1]
	*(*(*info).overlayparams).sxp = (*(*(*info).overlayparams).sxp)[0:(*(*info).loopparams).np-1]
	*(*(*info).overlayparams).syp = (*(*(*info).overlayparams).syp)[0:(*(*info).loopparams).np-1]
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
    CRISPEX_UPDATE_PHISLIT_COORDS, event
		CRISPEX_UPDATE_PHISLICE, event
	ENDIF ELSE CRISPEX_DRAW, event
END
;================================================================================= MASK PROCEDURES
PRO CRISPEX_MASK_OVERLAY_COLOR_SLIDER, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayparams).maskcolor = event.VALUE
	CRISPEX_DRAW, event
END

PRO CRISPEX_MASK_OVERLAY_SELECT_COLOR_TABLE, event
; Handles the change mask overlay window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayparams).maskct = event.index
	CRISPEX_DRAW, event
END

PRO CRISPEX_MASK_BUTTONS_SET, event
; Handles the setting of mask buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[0],SET_BUTTON = ((*(*info).overlayswitch).maskim)[0], SENSITIVE = (*(*info).overlayswitch).mask
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[1],SET_BUTTON = ((*(*info).overlayswitch).maskim)[1], SENSITIVE = ((*(*info).overlayswitch).mask AND (*(*info).dataswitch).reffile)
	WIDGET_CONTROL,((*(*info).ctrlscp).mask_button_ids)[2],SET_BUTTON = ((*(*info).overlayswitch).maskim)[2], SENSITIVE = ((*(*info).overlayswitch).mask AND (*(*info).winswitch).showdop)
END

PRO CRISPEX_MASK_OVERLAY_RASTER_TOGGLE, event
; Updates reference loopslab display window plot axes range according to set parameters
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  (*(*info).overlayswitch).sjiraster = event.SELECT
  CRISPEX_DRAW_SJI, event
END

;================================================================================= MEASUREMENT PROCEDURES
PRO CRISPEX_MEASURE_ARCSEC, event
; Handles the change of arcseconds per pixel resolution
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_text, GET_VALUE = textvalue
	(*(*info).meas).arcsecpix = FLOAT(textvalue[0])
	IF ((*(*info).meas).np GT 0) THEN CRISPEX_MEASURE_CALC, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).meas).arcsecpix], labels=['arcsecpix']
END

PRO CRISPEX_MEASURE_ENABLE, event
; Enables the spatial measurement tool
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_label, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_text, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).apix_unit, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_lab, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_asec_text, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_lab, SENSITIVE = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).measure_km_text, SENSITIVE = event.SELECT
	(*(*info).meas).spatial_measurement = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).meas).spatial_measurement], labels=['Enabled measurement']
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [delta_asec,delta_km], labels=['Arcseconds','Kilometers']
END

;================================================================================= PLAYBACK PROCEDURES
PRO CRISPEX_PB_BG, event
; Handles the actual playback, given the mode of playback
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CASE (*(*info).pbparams).mode OF
		'PAUSE'	: BEGIN
				IF (*(*info).pbparams).spmode THEN BEGIN
          IF ((*(*info).dataparams).lp NE (*(*info).pbparams).lp_blink) THEN $
            (*(*info).dataparams).lp = (*(*info).pbparams).lp_blink $
          ELSE $
            (*(*info).dataparams).lp = (*(*info).pbparams).lp_blink_init
			  ENDIF ELSE IF (*(*info).pbparams).imrefmode THEN BEGIN
					(*(*info).winids).imrefdisp = ABS((*(*info).winids).imrefdisp-1)
				ENDIF ELSE RETURN
				WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
			  END
		'PLAY'	: BEGIN
				IF ((*(*info).feedbparams).count_pbstats EQ 0) THEN (*(*info).feedbparams).pbstats = SYSTIME(/SECONDS)
				(*(*info).dispparams).t += (*(*info).pbparams).direction * (*(*info).pbparams).t_step
				IF (*(*info).pbparams).spmode THEN BEGIN
          IF ((*(*info).dataparams).lp NE (*(*info).pbparams).lp_blink) THEN $
            (*(*info).dataparams).lp = (*(*info).pbparams).lp_blink $
          ELSE $
            (*(*info).dataparams).lp = (*(*info).pbparams).lp_blink_init
					WIDGET_CONTROL,(*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
			  	ENDIF ELSE IF (*(*info).pbparams).imrefmode THEN BEGIN
					(*(*info).winids).imrefdisp = ABS((*(*info).winids).imrefdisp-1)
				ENDIF
				CASE (*(*info).pbparams).lmode OF
					'LOOP'	: BEGIN
							IF (*(*info).dispparams).t GT (*(*info).dispparams).t_upp THEN (*(*info).dispparams).t -= (*(*info).dispparams).t_range
							IF (*(*info).dispparams).t LT (*(*info).dispparams).t_low THEN (*(*info).dispparams).t += (*(*info).dispparams).t_range
						  END
					'CYCLE'	: BEGIN
							IF (*(*info).dispparams).t GT (*(*info).dispparams).t_upp THEN BEGIN
								(*(*info).pbparams).direction = -1
								(*(*info).dispparams).t = (*(*info).dispparams).t_upp - ((*(*info).dispparams).t MOD (*(*info).dispparams).t_upp)
							ENDIF ELSE IF (*(*info).dispparams).t LT (*(*info).dispparams).t_low THEN BEGIN
								(*(*info).pbparams).direction = 1
								IF ((*(*info).dispparams).t_low EQ (*(*info).dispparams).t_first) THEN (*(*info).dispparams).t *= -1 ELSE $
									(*(*info).dispparams).t = -1 * ((*(*info).dispparams).t - (*(*info).dispparams).t_low) + (*(*info).dispparams).t_low
							ENDIF
							IF ((*(*info).pbparams).direction GT 0) THEN CRISPEX_PB_BUTTONS_SET, event, /FWD_SET ELSE CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
						  END
					'BLINK'	: BEGIN
							(*(*info).pbparams).direction *= -1
							IF (*(*info).dispparams).t GT (*(*info).dispparams).t_upp THEN (*(*info).dispparams).t -=(*(*info).dispparams).t_range ELSE $
								IF (*(*info).dispparams).t LT (*(*info).dispparams).t_low THEN (*(*info).dispparams).t += (*(*info).dispparams).t_range
							IF ((*(*info).pbparams).direction GT 0) THEN CRISPEX_PB_BUTTONS_SET, event, /FWD_SET ELSE CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
						  END
				ENDCASE
			  END
	ENDCASE
	WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dispparams).t
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 1. / (*(*info).pbparams).t_speed
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2),STRTRIM((*(*info).pbparams).direction,2),$
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
		IF ((*(*info).feedbparams).count_pbstats GE 10) THEN $
      average = STRING(MEAN(*(*(*info).feedbparams).sum_pbstats, /NAN),FORMAT='(F6.4)')+' s' $
    ELSE $
      average = 'N/A'
		CRISPEX_UPDATE_USER_FEEDBACK, event, title='CRISPEX DEBUGGING: Playback statistics', var=((*(*info).feedbparams).count_pbstats+(*(*info).winids).feedbacktlb), minvar=1, /close_button, $
			feedback_text='Time elapsed since last update: '+STRING(timediff,FORMAT='(F6.4)')+' s, average (over last 10): '+average+', (over last '+$
			STRING((*(*info).feedbparams).count_pbstats,FORMAT='(I'+STRTRIM(FLOOR(ALOG10((*(*info).feedbparams).count_pbstats))+1,2)+')')+'): '+STRING((*(*info).feedbparams).av_pbstats,FORMAT='(F6.4)')+' s'
	ENDIF
END

PRO CRISPEX_PB_BUTTONS_SET, event, fbwd_set=fbwd_set, bwd_set=bwd_set, pause_set=pause_set, fwd_set=fwd_set, ffwd_set=ffwd_set, loop_set=loop_set, cycle_set=cycle_set, blink_set=blink_set
; Sets playback buttons according to actions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).pbparams).mode EQ 'PLAY' AND (*(*info).pbparams).direction EQ -1 THEN RETURN
	(*(*info).pbparams).direction = -1			&	(*(*info).pbparams).mode = 'PLAY'
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	CRISPEX_PB_BUTTONS_SET, event, /BWD_SET
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,STRTRIM((*(*info).pbparams).direction,2)], labels=['Play mode','Play direction']
END

PRO CRISPEX_PB_FORWARD, event
; Sets the playback mode to forward play
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).pbparams).mode EQ 'PLAY' AND (*(*info).pbparams).direction EQ 1 THEN RETURN
	(*(*info).pbparams).direction = 1			&	(*(*info).pbparams).mode = 'PLAY'
	WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	CRISPEX_PB_BUTTONS_SET, event, /FWD_SET
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,STRTRIM((*(*info).pbparams).direction,2)], labels=['Play mode','Play direction']
END

PRO CRISPEX_PB_FASTBACKWARD, event
; Sets the playback mode to single step backward
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).pbparams).mode EQ 'PAUSE' OR (*(*info).pbparams).lmode EQ 'BLINK' THEN BEGIN
		CRISPEX_PB_BUTTONS_SET, event, /FBWD_SET
		(*(*info).dispparams).t -= (*(*info).pbparams).t_step
		IF ((*(*info).dispparams).t LT (*(*info).dispparams).t_low) THEN (*(*info).dispparams).t = (*(*info).dispparams).t_upp
		CRISPEX_UPDATE_T, event
		WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dispparams).t
		CRISPEX_UPDATE_SLICES, event
		CRISPEX_DRAW, event
		CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
END

PRO CRISPEX_PB_FASTFORWARD, event
; Sets the playback mode to single step forward
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).pbparams).mode EQ 'PAUSE' OR (*(*info).pbparams).lmode EQ 'BLINK' THEN BEGIN
		CRISPEX_PB_BUTTONS_SET, event, /FFWD_SET
		(*(*info).dispparams).t = (((*(*info).dispparams).t - (*(*info).dispparams).t_low + (*(*info).pbparams).t_step) MOD ((*(*info).dispparams).t_upp - (*(*info).dispparams).t_low + 1)) + (*(*info).dispparams).t_low
		CRISPEX_UPDATE_T, event
		WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dispparams).t
		CRISPEX_UPDATE_SLICES, event
		CRISPEX_DRAW, event
		CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
END

PRO CRISPEX_PB_PAUSE, event
; Sets the playback mode pause
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).pbparams).mode EQ 'PAUSE' THEN RETURN
	(*(*info).pbparams).mode = 'PAUSE'
	IF ((*(*info).winids).feedbacktlb NE 0) THEN BEGIN
		(*(*info).feedbparams).count_pbstats = 0
		WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, /SENSITIVE
	ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /PAUSE_SET
	IF ((*(*info).dispparams).phislice_update NE 1) THEN CRISPEX_UPDATE_SLICES, event		
END

PRO CRISPEX_PB_BLINK, event
; Sets the playback mode to temporal blink
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).lmode = 'BLINK'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /BLINK_SET
END

PRO CRISPEX_PB_CYCLE, event
; Sets the playback mode to cycle
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).lmode = 'CYCLE'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /CYCLE_SET
END

PRO CRISPEX_PB_LOOP, event
; Sets the playback mode to loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).lmode = 'LOOP'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).mode,(*(*info).pbparams).lmode,STRTRIM((*(*info).dispparams).t,2)], labels=['Play mode','Loop mode','t']
	CRISPEX_PB_BUTTONS_SET, event, /LOOP_SET
END

PRO CRISPEX_PB_SPECTBLINK, event
; Sets the playback mode to spectral blink
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).spmode = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlscp).imref_blink_but, SENSITIVE=ABS((*(*info).pbparams).spmode-1)
	IF (*(*info).pbparams).spmode THEN BEGIN
		(*(*info).pbparams).spmode = 1	&	(*(*info).pbparams).spdirection = 1
    (*(*info).pbparams).lp_blink_init = (*(*info).dataparams).lp
		IF ((*(*info).feedbparams).count_pbstats EQ 0) THEN (*(*info).feedbparams).pbstats = SYSTIME(/SECONDS)
		WIDGET_CONTROL, (*(*info).pbparams).bg, TIMER = 0.0
	ENDIF ELSE BEGIN
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, $
      SENSITIVE=((*(*info).pbparams).lp_blink NE (*(*info).dataparams).lp)
		IF ((*(*info).winids).feedbacktlb NE 0) THEN BEGIN
			(*(*info).feedbparams).count_pbstats = 0
			WIDGET_CONTROL, (*(*info).ctrlsfeedb).close_button, /SENSITIVE
		ENDIF
	ENDELSE
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_speed_slider, SENSITIVE = (*(*info).pbparams).spmode
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).spmode,(*(*info).pbparams).spdirection], labels=['Spectral blink mode','Blink direction']
END

;========================= SPECTRAL PHI SLIT PROCEDURES
PRO CRISPEX_PHISLIT_DIRECTION, event
; Determines the direction of the spectral phi slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_UPDATE_PHISLIT_COORDS, event
	CRISPEX_UPDATE_PHISLICE, event 
	CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_PHISLIT_MOVE_FWD, event
; Handles the movement of the spectral phi slit forward along the slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  CRISPEX_UPDATE_PHISLIT_COORDS, event
	CRISPEX_UPDATE_PHISLICE, event 
	CRISPEX_COORDSLIDERS_SET, 1, 1, event
END

PRO CRISPEX_PHISLIT_MOVE_SET_COORDS, event
; Handles the coordinate update after the movement of the spectral phi slit along the slit
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

;========================= PREFERENCES ROUTINES
PRO CRISPEX_PREFERENCES_WINDOW, event
; Opens up the preferences window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	base 		= WIDGET_BASE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+$
    ': Preferences', GROUP_LEADER = (*(*info).winids).root, TLB_FRAME_ATTR = 1, $
    /TLB_KILL_REQUEST_EVENTS)
	main		= WIDGET_BASE(base, /COLUMN)
 
  ; Start-up preferences
  startup_base   = WIDGET_BASE(main, /COLUMN, /FRAME)
	startup_lab 	= WIDGET_LABEL(startup_base, VALUE = 'Start-up and Playback', /ALIGN_LEFT)
	startup_buts 	= WIDGET_BASE(startup_base, /COLUMN, /NONEXCLUSIVE)
	(*(*info).ctrlspref).startup_win = $
                  WIDGET_BUTTON(startup_buts, VALUE = 'Show start-up window', $
                    EVENT_PRO='CRISPEX_PREFERENCES_SET_STARTUPWIN')
	(*(*info).ctrlspref).startup_autopl	= $
                  WIDGET_BUTTON(startup_buts, VALUE='Start playing automatically', $
                  EVENT_PRO='CRISPEX_PREFERENCES_SET_AUTOPLAY')
	(*(*info).ctrlspref).displays_phislice = $
                    WIDGET_BUTTON(startup_buts, $
                      VALUE='Automatically update spectral Phi-slice', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_PHISLICE_UPDATE')
  ; Set buttons according to preferences	
  WIDGET_CONTROL, (*(*info).ctrlspref).startup_win, SET_BUTTON=(*(*info).prefs).startupwin
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_autopl, SET_BUTTON=(*(*info).prefs).autoplay
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_phislice, $
    SET_BUTTON=(*(*info).dispparams).phislice_update

  ; Displays preferences
	layout_base   	= WIDGET_BASE(main,/COLUMN, /FRAME)
	layout_lab 	    = WIDGET_LABEL(layout_base, VALUE='Layout', /ALIGN_LEFT)
	displays_buts 	= WIDGET_BASE(layout_base, /GRID_LAYOUT, COLUMN=2)
	(*(*info).prefs).bgplotcol_old = (*(*info).plotparams).bgplotcol
	(*(*info).ctrlspref).displays_bgcols = $
                    WIDGET_SLIDER(displays_buts, TITLE='Default background plot color', $
                      MIN=0, MAX=255, VALUE=(*(*info).plotparams).bgplotcol, /DRAG, $
		                  EVENT_PRO='CRISPEX_PREFERENCES_SET_BGPLOTCOL')
	(*(*info).ctrlspref).displays_plcols = $
                    WIDGET_SLIDER(displays_buts, TITLE='Default line plot color', $
                      MIN=0, MAX=55, VALUE=(*(*info).plotparams).plotcol, /DRAG, $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_PLOTCOL')
	displays_int_base = WIDGET_BASE(layout_base, /ROW, /NONEXCLUSIVE)
	(*(*info).ctrlspref).displays_interp = $
                    WIDGET_BUTTON(displays_int_base, VALUE='Interpolate spectral slices',$
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_INTERPOLATE')
	displays_preview= WIDGET_BUTTON(displays_int_base, VALUE='Preview', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_PREVIEW')

  scaling_base    = WIDGET_BASE(main, /COLUMN, /FRAME)
	scaling_lab 	  = WIDGET_LABEL(scaling_base, VALUE='Scaling', /ALIGN_LEFT)
  histo_base      = WIDGET_BASE(scaling_base, /ROW)
  histo_opt_lab   = WIDGET_LABEL(histo_base, VALUE='Default histogram optimisation value:', /ALIGN_LEFT)
  (*(*info).ctrlspref).histo_opt_txt = $
                    WIDGET_TEXT(histo_base, VALUE=STRTRIM((*(*info).prefs).histo_opt_val,2), $
                      /EDITABLE, XSIZE=11, EVENT_PRO='CRISPEX_PREFERENCES_SET_SCALING_HISTO_OPT')
  (*(*info).ctrlspref).gamma_label = $
                    WIDGET_LABEL(scaling_base, VALUE=STRING((*(*info).prefs).gamma_val, $
                      FORMAT='(F6.3)'), /ALIGN_CENTER)
  (*(*info).ctrlspref).gamma_slid = $
                    WIDGET_SLIDER(scaling_base, TITLE='Default gamma', MIN=0, MAX=1000, $
                      VALUE=500*(ALOG10((*(*info).prefs).gamma_val)+1), $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_SCALING_GAMMA', /SUPPRESS, /DRAG)
	displays_opts	  = WIDGET_BASE(scaling_base, /COLUMN, /NONEXCLUSIVE)
	(*(*info).ctrlspref).displays_slices = $
                    WIDGET_BUTTON(displays_opts, $
                      VALUE='Scale slices according main/reference image', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_SLICES_IMSCALE')
  ; Set buttons according to settings
	(*(*info).prefs).plotcol_old = (*(*info).plotparams).plotcol
	(*(*info).prefs).interpspslice_old = (*(*info).dispparams).interpspslice
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_interp, $
    SET_BUTTON=(*(*info).dispparams).interpspslice
	WIDGET_CONTROL, displays_preview, SET_BUTTON=(*(*info).prefs).preview
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_slices, $
    SET_BUTTON=(*(*info).dispparams).slices_imscale

  ; IO preferences: inputs
  paths_base      = WIDGET_BASE(main, /COLUMN, /FRAME)
	paths_lab	      = WIDGET_LABEL(paths_base, VALUE = 'Input/Output paths', /ALIGN_LEFT)
	paths_io_base	  = WIDGET_BASE(paths_base, /COLUMN)
	paths_i_labbuts = WIDGET_BASE(paths_io_base, /ROW)
	paths_i_lab	    = WIDGET_LABEL(paths_i_labbuts, VALUE = 'Default input path:')
	paths_i_buts	  = WIDGET_BASE(paths_i_labbuts, /ROW, /EXCLUSIVE)
	(*(*info).ctrlspref).paths_i_def_but = $
                    WIDGET_BUTTON(paths_i_buts, VALUE='Local working directory', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_IPATH_SEL_DEFAULT', /NO_RELEASE)
	(*(*info).ctrlspref).paths_i_sav_but = $
                    WIDGET_BUTTON(paths_i_buts, VALUE='Other directory', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_IPATH_SEL_OTHER', /NO_RELEASE)
	(*(*info).ctrlspref).paths_ipath_text = $
                    WIDGET_TEXT(paths_io_base, VALUE=(*(*info).prefs).prefipath, /EDITABLE,$
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_IPATH_OTHER', $
                      SENSITIVE=(*(*info).prefs).defipath)
  ; Set buttons according to settings
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_def_but, $
    SET_BUTTON=ABS((*(*info).prefs).defipath-1)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_i_sav_but, $
    SET_BUTTON=(*(*info).prefs).defipath
  ; IO preferences: outputs
	paths_o_labbuts = WIDGET_BASE(paths_io_base, /ROW)
	paths_o_lab	    = WIDGET_LABEL(paths_o_labbuts, VALUE = 'Default output path:')
	paths_o_buts	  = WIDGET_BASE(paths_o_labbuts, /ROW, /EXCLUSIVE)
	(*(*info).ctrlspref).paths_o_def_but = $
                    WIDGET_BUTTON(paths_o_buts, VALUE='Local working directory', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_OPATH_SEL_DEFAULT', /NO_RELEASE)
	(*(*info).ctrlspref).paths_o_sav_but = $
                    WIDGET_BUTTON(paths_o_buts, VALUE='Other directory', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_OPATH_SEL_OTHER', /NO_RELEASE)
	(*(*info).ctrlspref).paths_opath_text = $
                    WIDGET_TEXT(paths_io_base, VALUE=(*(*info).prefs).prefopath, /EDITABLE,$
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_OPATH_OTHER', $
                      SENSITIVE=(*(*info).prefs).defopath)
	(*(*info).ctrlspref).paths_iopath = $
                    WIDGET_BUTTON(paths_base, VALUE='Set output path to input path', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_IOPATH', $
                      SENSITIVE=((*(*info).prefs).prefipath NE (*(*info).prefs).prefopath))
  ; Set buttons according to settings
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_def_but, $
    SET_BUTTON=ABS((*(*info).prefs).defopath-1)
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_o_sav_but, $
    SET_BUTTON=(*(*info).prefs).defopath
  ; File-ID
	save_base	      = WIDGET_BASE(main, /COLUMN, /FRAME)
	save_lab	      = WIDGET_LABEL(save_base, VALUE='Filenaming', /ALIGN_LEFT)
	saveid_buts 	  = WIDGET_BASE(save_base, /ROW)
	saveid_lab	    = WIDGET_LABEL(saveid_buts, VALUE='Default unique file ID:', /ALIGN_LEFT)
	saveids	= ['YYYYMMMDD_hhmmss (default)','DDMMMYYYY_hhmmss', $
              'YYYYMMDD_hhmmss','DDMMYYYY_hhmmss']
	(*(*info).ctrlspref).save_defsaveid = $
                    WIDGET_COMBOBOX(saveid_buts, VALUE=saveids, /DYNAMIC_RESIZE, $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_SAVEID')
	CRISPEX_SAVE_DETERMINE_SAVEID, event, defsaveid_sample, /PREF
	saveid_sample	  = WIDGET_BASE(save_base, /ROW)
	saveid_sample_lab = WIDGET_LABEL(saveid_sample, VALUE = 'Example:')
	(*(*info).ctrlspref).save_defsaveid_sample = $
                    WIDGET_LABEL(saveid_sample, VALUE=defsaveid_sample, /DYNAMIC_RESIZE)
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid, $
    SET_COMBOBOX_SELECT=(*(*info).prefs).defsaveid

  ; Defaults / Cancel / Accept settings buttons
	dec_buts 	      = WIDGET_BASE(main, /ALIGN_CENTER, /GRID_LAYOUT, COLUMN=3)
	(*(*info).ctrlspref).set_defaults	= $
                    WIDGET_BUTTON(dec_buts, VALUE='Default settings', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SET_DEFAULTS', SENSITIVE=nondefault)
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
	cancel 		      = WIDGET_BUTTON(dec_buts, VALUE='Cancel', $
                      EVENT_PRO='CRISPEX_PREFERENCES_CANCEL')
	save_settings 	= WIDGET_BUTTON(dec_buts, VALUE='Save settings', $
                      EVENT_PRO='CRISPEX_PREFERENCES_SAVE_SETTINGS')

  ; Realize widget window
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET=(*(*info).winsizes).spxoffset, $
    TLB_SET_YOFFSET=0
	WIDGET_CONTROL, base, SET_UVALUE=info
	XMANAGER, 'CRISPEX', base, /NO_BLOCK
	(*(*info).winids).preftlb = base
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).winids).preftlb], labels=['preftlb']
END

PRO CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
; Handles the checking whether preference buttons and values are in their default position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (((*(*info).prefs).tmp_autoplay NE (*(*info).prefs).default_autoplay) OR $
		((*(*info).prefs).tmp_startupwin NE (*(*info).prefs).default_startupwin) OR $
		((*(*info).prefs).tmp_bgplotcol NE (*(*info).prefs).default_bgplotcol) OR $
		((*(*info).prefs).tmp_plotcol NE (*(*info).prefs).default_plotcol) OR $
		((*(*info).prefs).tmp_interpspslice NE (*(*info).prefs).default_interpspslice) OR $
		((*(*info).prefs).tmp_slices_imscale NE (*(*info).prefs).default_slices_imscale) OR $			
		((*(*info).prefs).tmp_histo_opt_val NE (*(*info).prefs).default_histo_opt_val) OR $			
		((*(*info).prefs).tmp_gamma_val NE (*(*info).prefs).default_gamma_val) OR $			
		((*(*info).prefs).tmp_phislice_update NE (*(*info).prefs).default_phislice_update) OR $			
		((*(*info).prefs).tmp_defipath NE (*(*info).prefs).default_defipath) OR $
		((*(*info).prefs).tmp_prefipath NE (*(*info).prefs).default_prefipath) OR $
		((*(*info).prefs).tmp_defopath NE (*(*info).prefs).default_defopath) OR $
		((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).default_prefopath) OR $
		((*(*info).prefs).tmp_defsaveid NE (*(*info).prefs).default_defsaveid)) THEN $
      nondefault = 1 $
    ELSE $
      nondefault = 0
	WIDGET_CONTROL, (*(*info).ctrlspref).set_defaults, SENSITIVE=nondefault
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [ABS(nondefault-1)], labels=['Buttons on default']
END

PRO CRISPEX_PREFERENCES_SET_STARTUPWIN, event
; Handles the toggle on/off setting of start-up window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_startupwin = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_startupwin], labels=['Startup window']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_AUTOPLAY, event
; Handles the toggle on/off setting of autoplay at start-up
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_autoplay = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_autoplay], labels=['Autoplay']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_BGPLOTCOL, event
; Handles the setting of the background plot color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_phislice_update = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_phislice_update], labels=['Live update spectral Phi-slice']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_PLOTCOL, event
; Handles the setting of the line plot color
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_slices_imscale = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_slices_imscale], $
      labels=['Scale slices with main/reference image']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_SCALING_HISTO_OPT, event
; Handles setting the default histrogram optimisation value
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlspref).histo_opt_txt, GET_VALUE = textvalue
	(*(*info).prefs).tmp_histo_opt_val = textvalue[0]
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_histo_opt_val], $
      labels=['Histogram optimisation value']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_SCALING_GAMMA, event
; Handles setting the default histrogram optimisation value
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  (*(*info).prefs).tmp_gamma_val = 10.^((FLOAT(event.VALUE)/500.) - 1.)
  WIDGET_CONTROL, (*(*info).ctrlspref).gamma_label, $
    SET_VALUE=STRING((*(*info).prefs).tmp_gamma_val,FORMAT='(F6.3)')
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_gamma_val], $
      labels=['Gamma']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_INTERPOLATE, event
; Handles the toggle on/off setting of interpolating the spectral slices
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_ipath_text, GET_VALUE = textvalue
	(*(*info).prefs).tmp_prefipath = textvalue[0]
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_prefipath], labels=['Input path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_OPATH_SEL_DEFAULT, event
; Handles the toggle on/off selection of the default output path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_opath_text, GET_VALUE = textvalue
	(*(*info).prefs).tmp_prefopath = textvalue[0]
	WIDGET_CONTROL, (*(*info).ctrlspref).paths_iopath, SENSITIVE = ((*(*info).prefs).tmp_prefopath NE (*(*info).prefs).tmp_prefipath)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_prefopath], labels=['Output path']
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_IOPATH, event
; Handles the setting of the output path to the input path
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_defsaveid = event.INDEX
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).prefs).tmp_defsaveid], labels=['Default save ID']
	CRISPEX_SAVE_DETERMINE_SAVEID, event, defsaveid_sample, /PREF
	WIDGET_CONTROL, (*(*info).ctrlspref).save_defsaveid_sample, SET_VALUE = defsaveid_sample
	CRISPEX_PREFERENCES_CHECK_DEFAULT_BUTTON, event
END

PRO CRISPEX_PREFERENCES_SET_DEFAULTS, event
; Handles the setting of all defaults
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).prefs).tmp_startupwin = (*(*info).prefs).default_startupwin
  (*(*info).prefs).tmp_interpspslice = (*(*info).prefs).default_interpspslice
	(*(*info).prefs).tmp_autoplay = (*(*info).prefs).default_autoplay
  (*(*info).prefs).tmp_defsaveid = (*(*info).prefs).default_defsaveid
	(*(*info).prefs).tmp_defipath = (*(*info).prefs).default_defipath
  (*(*info).prefs).tmp_prefipath = (*(*info).prefs).default_prefipath
	(*(*info).prefs).tmp_defopath = (*(*info).prefs).default_defopath
  (*(*info).prefs).tmp_prefopath = (*(*info).prefs).default_prefopath
	(*(*info).prefs).tmp_bgplotcol = (*(*info).prefs).default_bgplotcol
  (*(*info).prefs).tmp_plotcol = (*(*info).prefs).default_plotcol
	(*(*info).prefs).tmp_phislice_update = (*(*info).prefs).default_phislice_update
  (*(*info).prefs).tmp_slices_imscale = (*(*info).prefs).default_slices_imscale		
  (*(*info).prefs).tmp_histo_opt_val = (*(*info).prefs).default_histo_opt_val
  (*(*info).prefs).tmp_gamma_val = (*(*info).prefs).default_gamma_val
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_win, SET_BUTTON = (*(*info).prefs).tmp_startupwin
	WIDGET_CONTROL, (*(*info).ctrlspref).startup_autopl, SET_BUTTON = (*(*info).prefs).tmp_autoplay
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_bgcols, SET_VALUE = (*(*info).prefs).tmp_bgplotcol
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_plcols, SET_VALUE = (*(*info).prefs).tmp_plotcol
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_interp, SET_BUTTON = (*(*info).prefs).tmp_interpspslice
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_phislice, SET_BUTTON = (*(*info).prefs).tmp_phislice_update		
	WIDGET_CONTROL, (*(*info).ctrlspref).displays_slices, SET_BUTTON = (*(*info).prefs).tmp_slices_imscale
	WIDGET_CONTROL, (*(*info).ctrlspref).histo_opt_txt, $
    SET_VALUE=STRTRIM((*(*info).prefs).tmp_histo_opt_val,2)
	WIDGET_CONTROL, (*(*info).ctrlspref).gamma_label, $
    SET_VALUE=STRING((*(*info).prefs).tmp_gamma_val, FORMAT='(F6.3)') 
	WIDGET_CONTROL, (*(*info).ctrlspref).gamma_slid, $
    SET_VALUE=(500*(ALOG10((*(*info).prefs).tmp_gamma_val)+1))
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
  histo_opt_val = (*(*info).prefs).tmp_histo_opt_val
  gamma_val = (*(*info).prefs).tmp_gamma_val
	crispex_version = [(*(*info).versioninfo).version_number, (*(*info).versioninfo).revision_number]
	SAVE, crispex_version, startupwin, interpspslice, phislice_update, slices_imscale, histo_opt_val,$
    gamma_val, autoplay, defsaveid, defipath, defopath, bgplotcol, plotcol, prefipath, prefopath, $
    FILENAME = (*(*info).paths).dir_settings+'crispex.cpref'
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).dir_settings+'crispex.cpref'],labels=['Written']
	(*(*info).prefs).startupwin = startupwin		&	(*(*info).dispparams).interpspslice = interpspslice
	(*(*info).prefs).autoplay = autoplay			&	(*(*info).prefs).defsaveid = defsaveid
	(*(*info).prefs).defipath = defipath			&	(*(*info).prefs).defopath = defopath
	(*(*info).plotparams).bgplotcol = bgplotcol		&	(*(*info).plotparams).plotcol = plotcol
	(*(*info).prefs).prefipath = prefipath			&	(*(*info).paths).ipath = prefipath
	(*(*info).prefs).prefopath = prefopath			&	(*(*info).paths).opath = prefopath
	(*(*info).dispparams).phislice_update = phislice_update	&	(*(*info).dispparams).slices_imscale = slices_imscale	
  (*(*info).prefs).histo_opt_val = histo_opt_val
  (*(*info).prefs).gamma_val = gamma_val
	IF ~KEYWORD_SET(RESAVE) THEN BEGIN
		CRISPEX_PREFERENCES_REDRAW, event
		WIDGET_CONTROL, (*(*info).winids).preftlb, /DESTROY
		(*(*info).winids).preftlb = 0
	ENDIF
END

PRO CRISPEX_PREFERENCES_CANCEL, event
; Handles the exiting from the preferences window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).winids).sptlb NE 0) THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
	IF (TOTAL(*(*(*info).winids).restlooptlb) NE 0) THEN CRISPEX_DISPLAYS_RESTORE_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).retrdettlb NE 0) THEN CRISPEX_DISPLAYS_RETRIEVE_DET_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).looptlb NE 0) THEN CRISPEX_DISPLAYS_LOOPSLAB_REPLOT_AXES, event
	IF ((*(*info).winids).refsptlb NE 0) THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
	CRISPEX_DRAW_SPECTRAL, event
	CRISPEX_DRAW_TIMESLICES, event
	IF (*(*info).winswitch).showint THEN CRISPEX_DRAW_INT, event
END

;================================================================================= READ HEADER PROCEDURE
PRO CRISPEX_IO_PARSE_HEADER, filename, HDR_IN=hdr_in, HDR_OUT=hdr_out, $
                         IMCUBE=imcube, SPCUBE=spcube, REFIMCUBE=refimcube, REFSPCUBE=refspcube, $
                         SJICUBE=sjicube, MASKCUBE=maskcube, CUBE_COMPATIBILITY=cube_compatibility,$
                         EXTEN_NO=exten_no, SINGLE_CUBE=single_cube
; Handles read-in of file header, running different parsing depending on CUBE_COMPATIBILITY setting
  ; Start filling data header structure with header info from inputfile
  IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN
    offset = FITSPOINTER(filename, EXTEN_NO=exten_no, header, /SILENT)   ; Get header offset to data
    CRISPEX_READ_FITSHEADER, header, key, filename, $                    ; Parse FITS header into key struct
      IMCUBE=imcube, SPCUBE=spcube, REFIMCUBE=refimcube, REFSPCUBE=refspcube, $
      SJICUBE=sjicube, VERBOSE=hdr_out.verbosity[1]
  ENDIF ELSE BEGIN
    offset = 512                                                ; Set header offset to data
    CRISPEX_READ_HEADER, filename, datatype=datatype, $         ; Parse old header into variables
                         dims=dims, nx=nx, ny=ny, nt=nt, endian=endian, stokes=stokes, ns=ns, $
                         diagnostics=diagnostics
  ENDELSE
  hdr_out = hdr_in                                              ; Set output hdr to output hdr
  IF KEYWORD_SET(IMCUBE) THEN BEGIN                             ; Fill hdr_out parameters for IMCUBE
    hdr_out.imoffset = offset
    IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN              ; In case of FITS cube
      hdr_out.imtype = key.datatype     &  hdr_out.imnt = key.nlp * key.nt * key.ns
      hdr_out.nx = key.nx               &  hdr_out.dx = key.dx
      hdr_out.ny = key.ny               &  hdr_out.dy = key.dy
      hdr_out.imns = key.ns             &  hdr_out.ns = key.ns
      hdr_out.nlp = key.nlp             &  hdr_out.mainnt = key.nt
      hdr_out.blabel = key.btype        &  hdr_out.bunit = key.bunit
      hdr_out.xlabel = key.xlab         &  hdr_out.xunit = key.xunit     
      hdr_out.ylabel = key.ylab         &  hdr_out.yunit = key.yunit
      hdr_out.tlabel = key.tlab         &  hdr_out.tunit = key.tunit
      hdr_out.lplabel = key.lplab       &  hdr_out.lpunit = key.lpunit
      hdr_out.dt = key.dt               &  hdr_out.single_cube = key.nlp
      lc = key.lc
      hdr_out = CREATE_STRUCT(hdr_out, 'lps', key.lam, 'lc', lc)
      ; Handle spectral windows, if present
      hdr_out.ndiagnostics = key.ndiagnostics
      diagnostics = key.diagnostics
      wstart = key.wstart
      wwidth = key.wwidth
      twave = key.twave
      tarr_main = key.tarr_sel
      tarr_raster_main = key.tarr_raster
      toffset_main = key.tini_col
      xyrastersji = key.xyrastersji
      headers = key.headers
    ENDIF ELSE BEGIN                                            ; In case of compatibility mode
      hdr_out.imtype = datatype         &  hdr_out.imendian = endian
      hdr_out.nx = nx                   &  hdr_out.dx = 0.0592
      hdr_out.ny = ny                   &  hdr_out.dy = 0.0592
      hdr_out.imns = ns                 &  hdr_out.ns = ns
      hdr_out.imstokes = stokes         &  hdr_out.imnt = nt
      ndiagnostics = N_ELEMENTS(diagnostics)
      CRISPEX_IO_PARSE_SINGLE_CUBE, single_cube, HDR_IN=hdr_out, HDR_OUT=hdr_out,/MAIN
      IF (hdr_out.dt EQ 0) THEN $
        tarr_main = FINDGEN(hdr_out.mainnt) $
      ELSE $
        tarr_main = FINDGEN(hdr_out.mainnt) * hdr_out.dt
      tarr_raster_main = tarr_main
      toffset_main = 0
      IF (ndiagnostics GT 0) THEN BEGIN
        diagnostics = STRTRIM(STRSPLIT(STRMID(diagnostics,1,strlen(diagnostics)-2),',',$
                                /EXTRACT),2)
        hdr_out.ndiagnostics = ndiagnostics
        wstart = INDGEN(ndiagnostics)
        wwidth = REPLICATE(1,ndiagnostics)
      ENDIF ELSE BEGIN
        wstart = 0
        wwidth = hdr_out.nlp
        hdr_out.ndiagnostics = 1
        diagnostics = 'CRISP'
      ENDELSE
      xyrastersji = 0
      twave = 0
      headers = PTR_NEW('')
    ENDELSE
    hdr_out = CREATE_STRUCT(hdr_out, 'diagnostics', diagnostics, 'diag_start', wstart, $
      'diag_width', wwidth, 'tarr_main', tarr_main, 'tarr_raster_main', tarr_raster_main, $
      'toffset_main', toffset_main, 'xyrastersji', xyrastersji, 'twave', twave, 'hdrs_main', headers)
  ENDIF ELSE IF KEYWORD_SET(SPCUBE) THEN BEGIN                  ; Fill hdr parameters for SPCUBE
    hdr_out.spoffset = offset
    IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN              ; In case of FITS cube
      hdr_out.sptype = key.datatype     &  hdr_out.spnt = key.nx * key.ny * key.ns
      hdr_out.spns = key.ns
    ENDIF ELSE BEGIN                                            ; In case of compatibility mode
      hdr_out.sptype = datatype         &  hdr_out.spendian = endian
      hdr_out.nlp = nx                  &  hdr_out.mainnt = ny
      hdr_out.spns = ns                 &  hdr_out.spstokes = stokes
      hdr_out.spnt = nt
    ENDELSE
  ENDIF ELSE IF KEYWORD_SET(REFIMCUBE) THEN BEGIN               ; Fill hdr parameters for REFIMCUBE
    hdr_out.refimoffset = offset
    IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN              ; In case of FITS cube
      hdr_out.refimtype = key.datatype  &  hdr_out.refimnt = key.nlp * key.nt * key.ns
      hdr_out.refnx = key.nx            &  hdr_out.refny = key.ny
      hdr_out.refnlp = key.nlp          &  hdr_out.refnt = key.nt
      hdr_out.refns = key.ns
      hdr_out.refbunit = key.bunit      &  hdr_out.refblabel = key.btype
      hdr_out.refxlabel = key.xlab      &  hdr_out.refxunit = key.xunit
      hdr_out.refylabel = key.ylab      &  hdr_out.refyunit = key.yunit
      hdr_out.reflplabel = key.lplab    &  hdr_out.reflpunit = key.lpunit
      reflc = key.lc
      hdr_out = CREATE_STRUCT(hdr_out, 'reflps', key.lam, 'reflc', reflc)
      ; Handle spectral windows, if present
      hdr_out.nrefdiagnostics = key.ndiagnostics
      diagnostics = key.diagnostics
      wstart = key.wstart
      wwidth = key.wwidth
      twave = key.twave
      tarr_ref = key.tarr_sel
      tarr_raster_ref = key.tarr_raster
      toffset_ref = key.tini_col
      headers = key.headers
    ENDIF ELSE BEGIN                                            ; In case of compatibility mode
      hdr_out.refimtype = datatype      &  hdr_out.refimendian = endian
      hdr_out.refnx = nx                &  hdr_out.refny = ny
      hdr_out.refimnt = nt              &  hdr_out.refns = 1L
      IF (FLOOR(hdr_out.refimnt/FLOAT(hdr_out.mainnt)) EQ hdr_out.refimnt/FLOAT(hdr_out.mainnt)) THEN $
        hdr_out.refnt = hdr_out.mainnt ELSE hdr_out.refnt = 1L
      hdr_out.refnlp = LONG(hdr_out.refimnt/FLOAT(hdr_out.refnt))
      CRISPEX_IO_PARSE_SINGLE_CUBE, single_cube, HDR_IN=hdr_out, HDR_OUT=hdr_out,/REFERENCE
      IF (hdr_out.dt EQ 0) THEN $
        tarr_ref = FINDGEN(hdr_out.refnt) $
      ELSE $
        tarr_ref = FINDGEN(hdr_out.refnt) * hdr_out.dt
      tarr_raster_ref = tarr_ref
      toffset_ref = 0
      ndiagnostics = N_ELEMENTS(diagnostics)
      IF (ndiagnostics GT 0) THEN BEGIN
        diagnostics = STRTRIM(STRSPLIT(STRMID(diagnostics,1,strlen(diagnostics)-2),',',$
                                /EXTRACT),2)
        hdr_out.nrefdiagnostics = ndiagnostics
        wstart = INDGEN(ndiagnostics)
        wwidth = REPLICATE(1,ndiagnostics)
      ENDIF ELSE BEGIN
        wstart = 0
        wwidth = hdr_out.refnlp
        hdr_out.nrefdiagnostics = 1
        diagnostics = 'CRISP'
      ENDELSE
      twave = 0
      headers = PTR_NEW('')
    ENDELSE
    hdr_out = CREATE_STRUCT(hdr_out, 'refdiagnostics', diagnostics, 'refdiag_start', wstart, $
      'refdiag_width', wwidth, 'tarr_ref', tarr_ref, 'tarr_raster_ref', tarr_raster_ref, $
      'toffset_ref', toffset_ref, 'twave_ref', twave, 'hdrs_ref', headers)
  ENDIF ELSE IF KEYWORD_SET(REFSPCUBE) THEN BEGIN               ; Fill hdr parameters for REFSPCUBE
    hdr_out.refspoffset = offset
    IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN              ; In case of FITS cube
      hdr_out.refsptype = key.datatype  &  hdr_out.refspnt = key.nx * key.ny * key.ns
    ENDIF ELSE BEGIN                                            ; In case of compatibility mode
      hdr_out.refsptype = datatype      &  hdr_out.refspendian = endian
      hdr_out.refnlp = nx               &  hdr_out.refnt = ny
      hdr_out.refspnx = nx              &  hdr_out.refspny = ny
      hdr_out.refspnt = nt
    ENDELSE
  ENDIF ELSE IF KEYWORD_SET(SJICUBE) THEN BEGIN                ; Fill hdr parameters for SJICUBE
      hdr_out.sjioffset = offset
      hdr_out.sjitype = key.datatype    &  hdr_out.sjint = key.nt
      hdr_out.sjinx = key.nx            &  hdr_out.sjiny = key.ny
      hdr_out.sjidx = key.dx            &  hdr_out.sjidy = key.dy
      hdr_out.sjix0 = key.sjix0         &  hdr_out.sjiy0 = key.sjiy0
      hdr_out.sjibunit = key.bunit      &  headers = key.headers
      hdr_out = CREATE_STRUCT(hdr_out, 'tarr_sji', key.tarr_sel, $
        'sjixoff', key.sjixoff, 'sjiyoff', key.sjiyoff, 'hdrs_sji', headers)
  ENDIF ELSE IF KEYWORD_SET(MASKCUBE) THEN BEGIN                ; Fill hdr parameters for MASKCUBE
    hdr_out.maskoffset = offset
    IF ~KEYWORD_SET(CUBE_COMPATIBILITY) THEN BEGIN              ; In case of FITS cube
      hdr_out.masktype = key.datatype   &  hdr_out.masknt = key.nt
      hdr_out.masknx = key.nx           &  hdr_out.maskny = key.ny
    ENDIF ELSE BEGIN                                            ; In case of compatibility mode
      hdr_out.masktype = datatype       &  hdr_out.maskendian = endian
      hdr_out.masknx = nx               &  hdr_out.maskny = ny
      hdr_out.masknt = nt
    ENDELSE
  ENDIF
END

PRO CRISPEX_READ_FITSHEADER, header, key, filename, $
  IMCUBE=imcube, SPCUBE=spcube, REFIMCUBE=refcube, REFSPCUBE=refspcube, SJICUBE=sjicube, $
  VERBOSE=verbose
; Handle parsing of FITS file header
; Based on earlier PARSEHEADER.PRO, modification history:
;   v1.1 21-Sep-2012 Viggo Hansteen - first version
;   v1.2 26-Jan-2013 Mats Carlsson - works for both Bifrost and Iris cubes
;   v1.3 13-Feb-2013 Mats Carlsson 
; Incorporated functionality into CRISPEX on 29-May-2013 and extended subsequently
  hdr1 = ''   ; initialise header placeholders
  hdr2 = ''
  naxis = SXPAR(header,'NAXIS*')
  CASE (strsplit(SXPAR(header,'CTYPE2'),' ',/extract))[0] OF
    'y': sortorder = INDGEN(4)       ; CRISPEX imcube
    'HPLT-TAN': sortorder = INDGEN(3); IRIS SJI-file
    'time': sortorder = [2,3,0,1]    ; CRISPEX refcube
    'SolarY': sortorder = INDGEN(3)  ; Bifrost simcube
    ELSE: BEGIN
      MESSAGE,'Valid object not found! Correct fits file?',/cont
      key = {nx:-1}
      RETURN
    ENDELSE
  ENDCASE
  nx = naxis[sortorder[0]]
  ny = naxis[sortorder[1]]
  nlp = naxis[sortorder[2]]
  IF (N_ELEMENTS(naxis) EQ 4) THEN $
    nt = naxis[sortorder[3]] $
  ELSE $
    nt = 1
  ns = 1 ; number of Stokes posisitons
  cslab = [' ']
  ; Convert FITS datatype to IDL datatype
  CASE SXPAR(header,'BITPIX') OF
          8:      datatype = 1
         16:      datatype = 2  
         32:      datatype = 3  
        -32:      datatype = 4 
        -64:      datatype = 5 
          8:      datatype = 7
         16:      datatype = 12
         32:      datatype = 13
         64:      datatype = 14
          ELSE:   BEGIN
            MESSAGE,'ERROR: Illegal Image Datatype',/CONT
            datatype = -1
          ENDELSE
       endcase
  ; Read in header keywords      
  cdelt = SXPAR(header,'CDELT*')
  crpix = SXPAR(header,'CRPIX*')
  crval = SXPAR(header,'CRVAL*')
  ctype = SXPAR(header,'CTYPE*')
  cunit = SXPAR(header,'CUNIT*')
  btype = STRTRIM(SXPAR(header,'BTYPE'),2)
  bunit = STRTRIM(SXPAR(header,'BUNIT'))
  ; Assign values to variables
  dx = cdelt[sortorder[0]]
  dy = cdelt[sortorder[1]]
  IF (nt GT 1) THEN BEGIN
    dt = cdelt[sortorder[3]]
    tlab = STRTRIM(ctype[sortorder[3]],2)
    tunit = STRTRIM(cunit[sortorder[3]],2)
  ENDIF ELSE BEGIN
    dt = 0
    tlab = 't'
;    tunit = '[s]'
  ENDELSE
  tini_col = 0    ; Default raster timing column
  IF ~KEYWORD_SET(SJICUBE) THEN BEGIN
    ; Get time array (assuming each raster is co-temporal)
    IF (nt GT 1) THEN BEGIN
      tarr = READFITS(filename, hdr2, EXTEN_NO=2, SILENT=~KEYWORD_SET(VERBOSE))
      ntarrdims = SIZE(tarr,/N_DIMENSIONS)
      tarr_raster = tarr
      tval = SXPAR(header, 'CRVAL4')    ; tini_col = toffset_main/ref defaults to CRVAL4
      dum = MIN(ABS(tarr-tval),wheretval, /NAN)
;      wheretval = (WHERE(tarr EQ tval))[0]
      IF (wheretval EQ -1) THEN BEGIN
        nrasterpos = (SIZE(tarr))[ntarrdims-1]
        tini_col = FLOOR(nrasterpos/2.)
      ENDIF ELSE $
        tini_col = (ARRAY_INDICES(tarr,wheretval))[0]
      IF (ntarrdims EQ 2) THEN $
        tarr_sel = REFORM(tarr[tini_col,*]) $
      ELSE $
        tarr_sel = tarr
    ENDIF ELSE BEGIN
      tarr_sel = [0] 
      tarr_raster = 0
    ENDELSE
    sjixoff = 0
    sjiyoff = 0
    sjix0 = 0
    sjiy0 = 0
  ENDIF ELSE BEGIN
    offsetarray = READFITS(filename, hdr1, EXTEN_NO=1, SILENT=~KEYWORD_SET(VERBOSE))
    tarr_sel = REFORM(offsetarray[0,*]) ; TIME
    tarr_raster = tarr_sel
    sjixoff = REFORM(offsetarray[1,*])  ; PZTX
    sjiyoff = REFORM(offsetarray[2,*])  ; PZTY
    sjix0 = SXPAR(header, 'ROWSTAR')
    sjiy0 = SXPAR(header, 'COLSTAR')
    nt = naxis[sortorder[2]]
    tunit = STRTRIM(cunit[sortorder[2]],2)
  ENDELSE
  ; Determine tfactor and set time units to seconds
  common_tunit = [(tunit EQ 's'),(tunit EQ 'ms'),(tunit EQ 'hs')]
  where_common_tunit = WHERE(common_tunit EQ 1)
  IF (where_common_tunit NE -1) THEN $
    tfactor = ([1,0.001,100.])[where_common_tunit] $
  ELSE $
    tfactor = 1.
  tunit = 's'
  tarr_sel *= REPLICATE(tfactor,nt)
  tarr_raster *= REPLICATE(tfactor,(SIZE(tarr_raster))[1],nt)
  ; Determine plot labels
  xlab = STRTRIM(ctype[sortorder[0]],2)
  ylab = STRTRIM(ctype[sortorder[1]],2)
  lplab = STRTRIM(ctype[sortorder[2]],2)
  xunit = STRTRIM(cunit[sortorder[0]],2)
  yunit = STRTRIM(cunit[sortorder[1]],2)
  lpunit = STRTRIM(cunit[sortorder[2]],2)
  ; Determine spectral parameters
  IF (N_PARAMS() EQ 3) THEN BEGIN
    lam = READFITS(filename,hdr1,EXTEN_NO=1,SILENT=~KEYWORD_SET(VERBOSE))
  ENDIF ELSE BEGIN
    lam = (FINDGEN(nlp)+1-crpix[sortorder[2]])*cdelt[sortorder[2]]+crval[sortorder[2]]
  ENDELSE
  lcval = crval[sortorder[2]]
  lc = (WHERE(lam EQ lcval))[0]
  IF (lc EQ -1) THEN lc = 0
  ; Determine number of diagnostics
  ndiagnostics = SXPAR(header,'NWIN')
  IF (ndiagnostics GT 0) THEN BEGIN
    wstart = SXPAR(header,'WSTART*')
    wwidth = SXPAR(header,'WWIDTH*')
    wdesc = SXPAR(header,'WDESC*')
    whereselect = WHERE(wdesc NE '')
    wstart = wstart[whereselect]
    wwidth = wwidth[whereselect]
    diagnostics = wdesc[whereselect]
  ENDIF ELSE BEGIN
    wstart = 0
    wwidth = nlp
    whereselect = 0
    ndiagnostics = 1
    diagnostics = btype
  ENDELSE
  twave = SXPAR(header,'TWAVE*')
  twave  = twave[whereselect]
  ; Get third header too, even if not used in CRISPEX
  IF KEYWORD_SET(SJICUBE) THEN $
    dummy = READFITS(filename, hdr2, EXTEN_NO=2, SILENT=~KEYWORD_SET(VERBOSE))
  ; Initialise headers variable
  headers = [PTR_NEW(header),PTR_NEW(hdr1),PTR_NEW(hdr2)]
  ; Get slit coordinates on SJI image if raster
  IF ((SIZE(SXPAR(header,'SJIFIL*'),/TYPE) EQ 7) AND $
    (~KEYWORD_SET(SPCUBE) AND ~KEYWORD_SET(REFSPCUBE))) THEN BEGIN
    raster_coords = READFITS(filename, hdr3, EXTEN_NO=3, SILENT=~KEYWORD_SET(VERBOSE))
    xyrastersji = FLTARR(nx,2)
    xyrastersji[*,*] = ABS(raster_coords[*,0,*])
    headers = [headers, PTR_NEW(hdr3)]
  ENDIF ELSE xyrastersji = 0
  ;
  key = {nx:nx,ny:ny,nlp:nlp,nt:nt,ns:ns,cslab:cslab, $
       datatype:datatype,dx:dx,dy:dy,dt:dt,lam:lam,lc:lc, $
       tarr_sel:tarr_sel, tarr_raster:tarr_raster, tini_col:tini_col, xyrastersji:xyrastersji, $
       sjixoff:sjixoff, sjiyoff:sjiyoff, sjix0:sjix0, sjiy0:sjiy0, $
       xlab:xlab,ylab:ylab,lplab:lplab,tlab:tlab, $
       btype:btype,bunit:bunit, $ 
       xunit:xunit,yunit:yunit,lpunit:lpunit,tunit:tunit,$
       wstart:wstart, wwidth:wwidth, diagnostics:diagnostics, $
       ndiagnostics:ndiagnostics, twave:twave, headers:headers $
       }
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
	IF pos NE -1 THEN diagnostics = STRMID(header, pos + STRLEN(search), pos1-pos-STRLEN(search)+1)
END
 
;================================================================================= RESTORE LOOPS PROCEDURES
PRO CRISPEX_RESTORE_LOOPS_MAIN, event
; Start the restore loops procedures, opens the menu if CLSAV files are present or otherwise returns an error message
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).restoreparams).sel_loops)[eventval] = ( (*(*(*info).restoreparams).sel_loops)[eventval] EQ 0) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).restoreparams).sel_loops)[eventval]], labels=['Loop ID','Loop selected']
	CRISPEX_RESTORE_LOOPS_BUTTON_CONDITION, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_RESTORE_LOOPS_BUTTON_CONDITION, event
; Handles the update of buttons after selection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	condition = WHERE(*(*(*info).restoreparams).sel_loops EQ 1)
	WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_none, SET_BUTTON = ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1)
	WIDGET_CONTROL, (*(*info).ctrlsrestore).sel_all, SET_BUTTON = (((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).restoreparams).cfilecount))
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) AND (N_ELEMENTS(condition) EQ (*(*info).restoreparams).cfilecount)),$
		ABS(((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1))-1),N_ELEMENTS(condition)-(TOTAL(condition) EQ -1)], labels=['All selected','None selected','Total selected']
END

PRO CRISPEX_RESTORE_LOOPS_MENU_CLOSE, event
; Handles the closing of the restored loops menu and clean-up of display afterwards
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayswitch).overlalways = event.SELECT
	IF ((*(*info).loopswitch).restore_loops EQ 0) THEN CRISPEX_RESTORE_LOOPS_MAIN, event ELSE CRISPEX_DRAW, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).overlayswitch).overlalways], labels=['Overlay loops always']
END

PRO CRISPEX_RESTORE_LOOPS_OPEN_TANAT, event			
; Handles all prior to loading TANAT to analyse the selected loop slice/slab
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).overlayswitch).det_overlay_all = event.SELECT
	CRISPEX_DRAW, event
END

PRO CRISPEX_RETRIEVE_DET_WIDTH, event
; Handles the width of the retrieved detection
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	(*(*info).dispparams).t = t_det
END

PRO CRISPEX_RETRIEVE_DET_CANCEL, event
; Handles the closing of the detection file menu
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).winids).root, SET_UVALUE = info
	WIDGET_CONTROL, event.TOP, /DESTROY
	event.TOP = (*(*info).winids).root
	CRISPEX_RETRIEVE_LOOP_MENU, event
	CRISPEX_DRAW_IMREF, event
END

PRO CRISPEX_RETRIEVE_LOOP_DELETE_CLSAV, event
; Enables or disables the deletion of CLSAV files after saving the respective loopslabs
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).delete_clsav = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).delete_clsav],labels=['Delete loop path file']
END

PRO CRISPEX_RETRIEVE_LOOP_ALL_POS, event
; Enables the extraction of the retrieved loop at all or only the saved spectral positions
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).all_pos_loops = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).all_pos_loops],labels=['Saving spectral position setting']
END

PRO CRISPEX_RETRIEVE_LOOP_IMCUBE_ONLY, event
; Enables the retreival of loop paths from the image cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).imref_only = 1
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_REFCUBE_ONLY, event
; Enables the retreival of loop paths from the reference cube only
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).imref_only = 2
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_IMREF, event
; Enables the retreival of loop paths from both the image and reference cube
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).imref_only = 3
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).imref_only],labels=['Saving from cube setting']
END

PRO CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLAB, event
; Opens the warning windows giving the saving time estimate of the procedure, intermediate step towards saving all retrieved loopslabs (i.e. at all spectral positions)
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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

;========================= IMAGE SCALING PROCEDURES
PRO CRISPEX_SCALING_SELECT_DATA, event
; Handles the selection of scaling 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ;index = {0,1,2,3} -> Main, reference, Doppler, slit-jaw
  (*(*info).scaling).imrefscaling = event.INDEX
  (*(*info).scaling).idx = $
    (((*(*info).scaling).imrefscaling EQ 0) OR ((*(*info).scaling).imrefscaling EQ 2)) * $
    (*(*info).intparams).lp_diag_all + $
    (((*(*info).scaling).imrefscaling GT 0) + ((*(*info).scaling).imrefscaling GT 2)) * $
    (*(*info).intparams).ndiagnostics + $
    ((*(*info).scaling).imrefscaling EQ 1) * (*(*info).intparams).lp_ref_diag_all + $
    ((*(*info).scaling).imrefscaling GT 1) * (*(*info).intparams).nrefdiagnostics ;+ $
;    ((*(*info).scaling).imrefscaling GT 2)
;	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
;    CRISPEX_VERBOSE_GET, event, [((*(*info).scaling).imrefscaling EQ 2)],labels=['Doppler scaling select']
	CRISPEX_SCALING_SET_BOXBUTTONS, event;, SENSITIVITY_OVERRIDE=sensitivity_override
	CRISPEX_SCALING_SET_SLIDERS, event
END

PRO CRISPEX_SCALING_SELECT_TYPE, event
; Handles the selection of scaling 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ;index = {0,1,2} -> Based on first, based on current, per time step
  (*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling] = event.INDEX
;  IF (event.INDEX EQ 0) THEN $
;    CRISPEX_SCALING_SELECT_FIRST, event $
;  ELSE IF (event.INDEX EQ 1) THEN $
;    CRISPEX_SCALING_SELECT_CURRENT, event
  CRISPEX_SCALING_APPLY_SELECTED, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,$
    (*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling]],$
		labels=['Scaling image','Scaling setting']
	CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_APPLY_SELECTED, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Main image
  idx = (*(*info).intparams).lp_diag_all
  IF ((*(*(*info).scaling).imagescale)[0] EQ 0) THEN BEGIN
    minmax_data = (*(*(*info).data).imagedata)[$
      (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp]
    minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
      (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
    (*(*info).scaling).imagemin = minmax[0]
    (*(*info).scaling).imagemax = minmax[1]
  ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[0] EQ 1) THEN BEGIN
    minmax_data = *(*(*info).data).xyslice
    minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
      (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
    (*(*info).scaling).imagemin_curr = minmax[0]
    (*(*info).scaling).imagemax_curr = minmax[1]
  ENDIF
  ; Reference image
  IF (*(*info).winswitch).showref THEN BEGIN
    idx = (*(*info).intparams).ndiagnostics + (*(*info).intparams).lp_ref_diag_all
    IF ((*(*(*info).scaling).imagescale)[1] EQ 0) THEN BEGIN
      minmax_data = (*(*(*info).data).refdata)[(*(*info).dataparams).lp_ref]
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).refmin = minmax[0]
      (*(*info).scaling).refmax = minmax[1]
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[1] EQ 1) THEN BEGIN
      minmax_data = *(*(*info).data).refslice
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).refmin_curr = minmax[0]
      (*(*info).scaling).refmax_curr = minmax[1]
    ENDIF
  ENDIF
  ; Doppler image
  IF (*(*info).winswitch).showdop THEN BEGIN
    idx = (*(*info).intparams).ndiagnostics + $
      (*(*info).intparams).nrefdiagnostics + (*(*info).intparams).lp_diag_all
    IF ((*(*(*info).scaling).imagescale)[2] EQ 0) THEN BEGIN
      minmax_data = (*(*info).data).dopplerscan[*,*,$
        (*(*info).dataparams).s * (*(*info).dataparams).nlp + (*(*info).dataparams).lp]
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).dopmin = minmax[0]
      (*(*info).scaling).dopmax = minmax[1]
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[2] EQ 1) THEN BEGIN
      minmax_data = *(*(*info).data).dopslice
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).dopmin_curr = minmax[0]
      (*(*info).scaling).dopmax_curr = minmax[1]
    ENDIF
  ENDIF
  ; SJI image
  IF (*(*info).winswitch).showsji THEN BEGIN
    idx = 2*(*(*info).intparams).ndiagnostics + (*(*info).intparams).nrefdiagnostics
    IF ((*(*(*info).scaling).imagescale)[3] EQ 0) THEN BEGIN
      minmax_data = (*(*(*info).data).sjidata)[0]
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).sjimin = minmax[0]
      (*(*info).scaling).sjimax = minmax[1]
    ENDIF ELSE IF ((*(*(*info).scaling).imagescale)[3] EQ 1) THEN BEGIN
      minmax_data = *(*(*info).data).sjislice
      minmax = CRISPEX_SCALING_SLICES(minmax_data, (*(*info).scaling).gamma[idx], $
        (*(*info).scaling).histo_opt_val[idx], /FORCE_HISTO)
      (*(*info).scaling).sjimin_curr = minmax[0]
      (*(*info).scaling).sjimax_curr = minmax[1]
    ENDIF
  ENDIF
END

PRO CRISPEX_SCALING_HISTO_OPT_VALUE, event
; Handles the setting of the HISTO_OPT parameter
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).histo_opt_txt, GET_VALUE = textvalue
  textvalue = FLOAT(textvalue[0])
	(*(*info).scaling).histo_opt_val[(*(*info).scaling).idx] = FLOAT(textvalue[0])
  CRISPEX_SCALING_APPLY_SELECTED, event
  CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_SLIDER_MIN, event
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  idx = (*(*info).scaling).idx
  (*(*info).scaling).minimum[idx] = event.VALUE
  IF ((*(*info).scaling).minimum[idx] GE (*(*info).scaling).maximum[idx]) THEN BEGIN
    (*(*info).scaling).maximum[idx] = (*(*info).scaling).minimum[idx] + 1
    CRISPEX_SCALING_SET_SLIDERS, event
  ENDIF
  CRISPEX_SCALING_APPLY_SELECTED, event
  CRISPEX_SCALING_REDRAW, event
  ; CRISPEX_DRAW, event
END

PRO CRISPEX_SCALING_SLIDER_MAX, event
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  idx = (*(*info).scaling).idx
  (*(*info).scaling).maximum[idx] = event.VALUE
  IF ((*(*info).scaling).maximum[idx] LE (*(*info).scaling).minimum[idx]) THEN BEGIN
    (*(*info).scaling).minimum[idx] = (*(*info).scaling).maximum[idx] - 1
    CRISPEX_SCALING_SET_SLIDERS, event
  ENDIF
  CRISPEX_SCALING_APPLY_SELECTED, event
  CRISPEX_SCALING_REDRAW, event
  ; CRISPEX_DRAW, event
END

PRO CRISPEX_SCALING_GAMMA_SLIDER, event
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  (*(*info).scaling).gamma[(*(*info).scaling).idx] = 10.^((FLOAT(event.VALUE)/500.) - 1.)
  WIDGET_CONTROL, (*(*info).ctrlscp).gamma_label, $
    SET_VALUE=STRING((*(*info).scaling).gamma[(*(*info).scaling).idx],FORMAT='(F6.3)')
  CRISPEX_SCALING_APPLY_SELECTED, event
  CRISPEX_SCALING_REDRAW, event
END

PRO CRISPEX_SCALING_RESET_DEFAULTS, event, IDX=idx, NO_DRAW=no_draw
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(IDX) NE 1) THEN idx = (*(*info).scaling).idx
  ; Reset contrast value
  (*(*info).scaling).minimum[idx] = 0
  (*(*info).scaling).maximum[idx] = 100
  (*(*info).scaling).gamma[idx] = (*(*info).prefs).default_gamma_val
  (*(*info).scaling).histo_opt_val[idx] = (*(*info).prefs).default_histo_opt_val
  IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
    CRISPEX_SCALING_APPLY_SELECTED, event
    CRISPEX_SCALING_SET_BOXBUTTONS, event
    CRISPEX_SCALING_SET_SLIDERS, event;, SET_GAMMA=50, /GAMMA_ONLY
    CRISPEX_SCALING_REDRAW, event
  ENDIF
END 

PRO CRISPEX_SCALING_RESET_ALL_DEFAULTS, event
; Handles events from the minimum scaling value slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  FOR i=0,N_ELEMENTS((*(*info).scaling).gamma)-1 DO $
    CRISPEX_SCALING_RESET_DEFAULTS, event, IDX=i, $
    NO_DRAW=(i NE (N_ELEMENTS((*(*info).scaling).gamma)-1))
END 

PRO CRISPEX_SCALING_REDRAW, event
; Handles the redrawing of window contents after adjustment of scaling
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).scaling).imrefscaling EQ 0) THEN BEGIN
		CRISPEX_DRAW_XY, event 
		IF (*(*info).winswitch).showsp THEN BEGIN
      CRISPEX_UPDATE_SPSLICE, event
      CRISPEX_DRAW_SPECTRAL_MAIN, event, /SP_ONLY
    ENDIF
		IF (*(*info).winswitch).showphis THEN BEGIN
      CRISPEX_UPDATE_PHISLICE, event, /NO_DRAW
      CRISPEX_DRAW_PHIS, event 
    ENDIF
	ENDIF ELSE BEGIN
		IF (*(*info).winswitch).showref THEN CRISPEX_DRAW_REF, event
		IF (*(*info).winswitch).showrefsp THEN BEGIN
      CRISPEX_UPDATE_REFSPSLICE, event
      CRISPEX_DRAW_SPECTRAL_REF, event, /SP_ONLY
    ENDIF
		IF (*(*info).winswitch).showdop THEN CRISPEX_DRAW_DOPPLER, event
    IF (*(*info).winswitch).showsji THEN CRISPEX_DRAW_SJI, event
	ENDELSE
;	IF ((*(*info).dispparams).slices_imscale AND (*(*info).winswitch).showrestloop) THEN $
	IF (*(*info).winswitch).showrestloop THEN $
    CRISPEX_DRAW_REST_LOOP, event
	IF (*(*info).winswitch).showloop THEN CRISPEX_DRAW_LOOPSLAB, event 
	IF (*(*info).winswitch).showrefloop THEN CRISPEX_DRAW_REFLOOPSLAB, event 
END

PRO CRISPEX_SCALING_SET_BOXBUTTONS, event, SENSITIVITY_OVERRIDE=sensitivity_override
; Handles the setting of scaling buttons
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
;  IF (N_ELEMENTS(SENSITIVITY_OVERRIDE) NE 1) THEN BEGIN
  showarr = [1,(*(*info).winswitch).showref,(*(*info).dispswitch).drawdop,$
    (*(*info).winswitch).showsji]
  WIDGET_CONTROL, (*(*info).ctrlscp).imagescale_cbox, $
    SET_COMBOBOX_SELECT=(*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling], $
    SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
  WIDGET_CONTROL, (*(*info).ctrlscp).scaling_reset_button, $
      SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
  WIDGET_CONTROL, (*(*info).ctrlscp).diagscale_label, $
    SET_VALUE=(*(*info).scaling).diagscale_label_vals[(*(*info).scaling).idx]
  ; HISTO_OPT controls
  WIDGET_CONTROL, (*(*info).ctrlscp).histo_opt_txt, $
    SET_VALUE=STRTRIM((*(*info).scaling).histo_opt_val[(*(*info).scaling).idx],2), $
    SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
  ; Multiply controls
	WIDGET_CONTROL, (*(*info).ctrlscp).ls_mult_cbox, $
    SET_COMBOBOX_SELECT=(*(*info).scaling).idx, SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
	WIDGET_CONTROL, (*(*info).ctrlscp).ls_mult_txt, $
    SET_VALUE=STRTRIM((*(*info).scaling).mult_val[(*(*info).scaling).mult_diag],2), $
    SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
;	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
;    CRISPEX_VERBOSE_GET, event, [(*(*info).scaling).imrefscaling,$
;    (*(*(*info).scaling).imagescale)[(*(*info).scaling).imrefscaling],$
;		ABS((*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling]-1),$
;    (*(*(*info).scaling).relative)[(*(*info).scaling).imrefscaling]],$
;		labels=['Image scaling','Scaling button set','Absolute scaling set','Relative scaling set']
END

PRO CRISPEX_SCALING_SET_SLIDERS, event, GAMMA_ONLY=gamma_only, SET_GAMMA=set_gamma, $
  SENSITIVITY_OVERRIDE=sensitivity_override
; Handles the setting of scaling sliders
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  showarr = [1,(*(*info).winswitch).showref,(*(*info).dispswitch).drawdop,$
    (*(*info).winswitch).showsji]
  IF ~KEYWORD_SET(GAMMA_ONLY) THEN BEGIN
    WIDGET_CONTROL, (*(*info).ctrlscp).scalemin_slider,$
      SET_VALUE=(*(*info).scaling).minimum[(*(*info).scaling).idx], $
      SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
    WIDGET_CONTROL, (*(*info).ctrlscp).scalemax_slider,$
      SET_VALUE=(*(*info).scaling).maximum[(*(*info).scaling).idx], $
      SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
  ENDIF
  IF (N_ELEMENTS(SET_GAMMA) EQ 1) THEN BEGIN
    gamma_slider_val = set_gamma
    gamma_val = 10.^((FLOAT(gamma_slider_val)/500.) - 1.)
    (*(*info).scaling).gamma[(*(*info).scaling).idx] = gamma_val
  ENDIF ELSE BEGIN
    gamma_val = (*(*info).scaling).gamma[(*(*info).scaling).idx]
    gamma_slider_val = 500 * (ALOG10(gamma_val) + 1)
  ENDELSE
  WIDGET_CONTROL, (*(*info).ctrlscp).gamma_slider, SET_VALUE=gamma_slider_val, $
    SENSITIVE=showarr[(*(*info).scaling).imrefscaling]
  WIDGET_CONTROL, (*(*info).ctrlscp).gamma_label, SET_VALUE=STRING(gamma_val,FORMAT='(F6.3)')
END

PRO CRISPEX_SCALING_MULTIPLY_LS_SELECT, event
; Handles the selection of diagnostic for multiplying detailed spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  (*(*info).scaling).mult_diag = event.INDEX
  WIDGET_CONTROL, (*(*info).ctrlscp).ls_mult_txt, $
    SET_VALUE=STRTRIM((*(*info).scaling).mult_val[(*(*info).scaling).mult_diag],2)
END  

PRO CRISPEX_SCALING_MULTIPLY_LS_VALUE, event
; Handles the selection of diagnostic for multiplying detailed spectrum
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).ctrlscp).ls_mult_txt, GET_VALUE = textvalue
  textvalue = FLOAT(textvalue[0])
  IF (textvalue EQ 0) THEN BEGIN
    textvalue = 1.
	  WIDGET_CONTROL, (*(*info).ctrlscp).ls_mult_txt, SET_VALUE=STRTRIM(textvalue,2)
  ENDIF
	(*(*info).scaling).mult_val[(*(*info).scaling).mult_diag] = FLOAT(textvalue[0])
  CRISPEX_DRAW_SPECTRAL, event
  CRISPEX_DRAW_TIMESLICES, event
END

;================================================================================= SESSION SAVE/RESTORE PROCEDURES
PRO CRISPEX_SESSION_SAVE_WINDOW, event
; Gets the filename for the session save routine
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_GET_FILENAME, event, 'Save session', 'crispex_session',$
    'CRISPEX_SESSION_SAVE_CONTINUE', /SESSION_SAVE
END

PRO CRISPEX_SESSION_SAVE_CONTINUE, event
; Checks the session save filename for validity and overwrite problems
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'SESSION'
	CRISPEX_SAVE_CHECK_FILENAME, event, 'cses', 'CRISPEX_SESSION_SAVE_OVER_CONTINUE'
END

PRO CRISPEX_SESSION_SAVE_OVER_CONTINUE, event
; Handles the overwriting and activates the subsequent saving of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, (*(*info).savparams).filename_text, GET_VALUE = session_filename
	CRISPEX_SESSION_SAVE, event, session_filename
	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SESSION_SAVE, event, sesfilename
; Handles the actual saving of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	WIDGET_CONTROL, /HOURGLASS
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	ctrlsswitch = *(*info).ctrlsswitch	    &	curs = *(*info).curs
	dataparams = *(*info).dataparams	      &	dataswitch = *(*info).dataswitch		
  detparams = *(*info).detparams          & dispparams = *(*info).dispparams	
  dispswitch = *(*info).dispswitch		    &	intparams = *(*info).intparams
  ioparams = *(*info).ioparams            &	loopparams = *(*info).loopparams	
  loopswitch = *(*info).loopswitch		    &	meas = *(*info).meas
	overlayparams = *(*info).overlayparams	&	overlayswitch = *(*info).overlayswitch		
  paramswitch = *(*info).paramswitch      & pbparams = *(*info).pbparams
  phiparams = *(*info).phiparams          & plotaxes = *(*info).plotaxes
  plotparams = *(*info).plotparams        & plotpos = *(*info).plotpos		
  plotswitch = *(*info).plotswitch		    & plottitles = *(*info).plottitles
	restoreparams = *(*info).restoreparams	&	retrparams = *(*info).retrparams		
  savswitch = *(*info).savswitch          & scaling = *(*info).scaling		
  stokesparams = *(*info).stokesparams		& versioninfo = *(*info).versioninfo
	winsizes = *(*info).winsizes		        &	winswitch = *(*info).winswitch			
  zooming = *(*info).zooming
	SAVE, ctrlsswitch, curs, dataparams, dataswitch, detparams, dispparams, dispswitch, intparams, $
        ioparams, loopparams, loopswitch, meas, overlayparams, overlayswitch, paramswitch, $
        pbparams, phiparams, plotaxes, plotparams, plotpos, plotswitch, plottitles, restoreparams,$
        retrparams, savswitch, scaling, stokesparams, versioninfo, winsizes, winswitch, zooming, $
		    FILENAME = (*(*info).paths).opath+sesfilename+'.cses'
	PRINT,'Written: '+(*(*info).paths).opath+sesfilename+'.cses'
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SESSION_RESTORE_WINDOW, event
; Opens the session restore window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, event.ID, GET_UVALUE = eventval
	(*(*(*info).sesparams).sessions)[eventval] = ( (*(*(*info).sesparams).sessions)[eventval] EQ 0) 
	condition = WHERE(*(*(*info).sesparams).sessions EQ 1)
	WIDGET_CONTROL, (*(*info).sesparams).rest_sessions, SENSITIVE = ((N_ELEMENTS(condition) GT 0) AND (TOTAL(condition) NE -1)) 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [eventval,(*(*(*info).retrparams).sel_loops)[eventval]], labels=['Restored session ID','Restored session selected']
END

PRO CRISPEX_SESSION_RESTORE_READ_POINTER, event, currpointer, restpointer, NO_RESTORE=no_restore
; Handles the actual restoration of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	currtags = TAG_NAMES(currpointer)	&	resttags = TAG_NAMES(restpointer)
	ncurr = N_ELEMENTS(currtags)		  &	nrest = N_ELEMENTS(resttags)
	no_rest = N_ELEMENTS(NO_RESTORE)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [ncurr,nrest,no_rest-1],labels=['Current tags','Restored tags',$
                          'Prevent replace tag']
  ; Both the restored and current pointer have the same number of tags
	IF (ncurr EQ nrest) THEN BEGIN											
		IF (no_rest EQ 0) THEN BEGIN										; Pointer to be restored without skipping tags
      ; If all tagnames are equal, then just replace the pointer
			IF (WHERE((resttags EQ currtags) EQ 0) EQ -1) THEN currpointer = restpointer ELSE BEGIN		
        ; Else go through all tags and replace only where tagnames are the same
				FOR i=0,nrest-1 DO (currpointer).(WHERE(currtags EQ resttags[i])) = restpointer.(i)	
			ENDELSE
		ENDIF ELSE BEGIN											; Pointer to be restored while skipping tags
			no_replace = WHERE((STRLOWCASE(resttags) EQ STRLOWCASE(no_restore)) EQ 1)
			FOR i=0,nrest-1 DO BEGIN
        ; If pointer tag is not the one to be skipped, replace its value
				IF (no_replace NE i) THEN (currpointer).(i) = restpointer.(i)				
			ENDFOR
		ENDELSE
	ENDIF ELSE BEGIN												; Unequal number of tags between current and restored pointer
		FOR i=0,nrest-1 DO BEGIN
      ; Go through all tags and replace only where tagnames are the same
			IF (WHERE(currtags EQ resttags[i]) NE -1) THEN $
        (currpointer).(WHERE(currtags EQ resttags[i])) = restpointer.(i)			
		ENDFOR
	ENDELSE
END

PRO CRISPEX_SESSION_RESTORE, event
; Handles the actual restoration of the session
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL, /HOURGLASS
	CRISPEX_UPDATE_USER_FEEDBACK, event, title='Restoring session...', var=0, feedback_text='Restoring session file and checking version...', /SESSION
	restore_session = (*(*(*info).sesparams).csesfiles)[WHERE(*(*(*info).sesparams).sessions EQ 1)]
	RESTORE, restore_session
	; Check revision number
	IF (N_ELEMENTS(versioninfo) GT 0) THEN cont = (versioninfo.revision_number GE '546') ELSE cont = 0
	IF cont THEN BEGIN
		IF (((*(*info).dataparams).imfilename EQ dataparams.imfilename) AND ((*(*info).dataparams).spfilename EQ dataparams.spfilename) AND $
			((*(*info).dataparams).refimfilename EQ dataparams.refimfilename) AND ((*(*info).dataparams).refspfilename EQ dataparams.refspfilename)) THEN BEGIN
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
			CRISPEX_SESSION_RESTORE_READ_POINTER, event, *(*info).ioparams, ioparams
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
			WIDGET_CONTROL, (*(*info).ctrlscp).t_slider, SET_VALUE = (*(*info).dispparams).t, SET_SLIDER_MIN = (*(*info).dispparams).t_low, SET_SLIDER_MAX = (*(*info).dispparams).t_upp, $
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
  	WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_slider, SET_SLIDER_MIN=(*(*info).dispparams).lp_low, $
      SET_SLIDER_MAX=(*(*info).dispparams).lp_upp, SET_VALUE=(*(*info).dataparams).lp, $
      SENSITIVE=(((*(*info).dispparams).lp_range-1) NE 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, SET_BUTTON = (*(*info).pbparams).spmode, SENSITIVE = (((*(*info).dataparams).nlp GT 1) AND ((*(*info).winswitch).showimref EQ 0))
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SET_BUTTON = ((*(*info).ctrlsswitch).lp_ref_lock AND ((*(*info).dataparams).refnlp GT 1)), $
				SENSITIVE = (((*(*info).dataparams).nlp EQ (*(*info).dataparams).refnlp) AND ((*(*info).dataparams).refnlp GT 1))
			WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SET_VALUE = (*(*info).dataparams).lp_ref, SENSITIVE = ((*(*info).dataswitch).refspfile AND ABS((*(*info).ctrlsswitch).lp_ref_lock-1))
			CRISPEX_PB_BUTTONS_SET, event, bwd_set=((*(*info).pbparams).direction EQ -1), pause_set=((*(*info).pbparams).mode EQ 'PAUSE'), fwd_set=((*(*info).pbparams).direction EQ 1), $
				loop_set=((*(*info).pbparams).lmode EQ 'LOOP'),	cycle_set=((*(*info).pbparams).lmode EQ 'CYCLE'), blink_set=((*(*info).pbparams).lmode EQ 'BLINK')
			; ==================== Spatial Tab ====================
      ; Set cursor sliders and lock/unlock buttons
			CRISPEX_COORDSLIDERS_SET, ABS((*(*info).curs).lockset-1), ABS((*(*info).curs).lockset-1), event
			WIDGET_CONTROL, (*(*info).ctrlscp).lock_button, SET_BUTTON = (*(*info).curs).lockset
			WIDGET_CONTROL, (*(*info).ctrlscp).unlock_button, SET_BUTTON = ABS((*(*info).curs).lockset-1)
      ; Set zooming buttons
			IF ((*(*info).zooming).factor NE 1) THEN CRISPEX_ZOOM, event, /NO_DRAW
      set_zoomfac = CRISPEX_BGROUP_ZOOMFAC_SET(event, /NO_DRAW, /NO_UPDATE_SLIDERS, $
        SET_FACTOR=WHERE((*(*info).zooming).factors EQ (*(*info).zooming).factor))
      ; Set scrolling sliders
			WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SET_VALUE = (*(*info).zooming).xpos, $
        SENSITIVE = ((*(*info).zooming).factor NE 1)
			WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SET_VALUE = (*(*info).zooming).ypos, $
        SENSITIVE = ((*(*info).zooming).factor NE 1)
			; ==================== Stokes Tab ====================
      FOR i=0,N_ELEMENTS((*(*info).stokesparams).button_labels)-1 DO BEGIN
        ; Stokes parameter available in data?
        stokes_enabled = (WHERE(((*(*info).stokesparams).labels) EQ $
                                 (*(*info).stokesparams).button_labels[i]) GE 0)
        ; Stokes parameter selected for detailed spectrum plot?
			  stokes_sp_select = (*(*info).stokesparams).labels[$
                           (WHERE((*(*info).stokesparams).select_sp EQ 1))]
        ; Determine setting of buttons accordingly
        stokes_sp_set = (WHERE(stokes_sp_select EQ (*(*info).stokesparams).button_labels[i]) GE 0)
        ; Total number of selected Stokes parameters for detailed spectrum plot?
        total_sp_select = TOTAL((*(*info).stokesparams).select_sp)
        ; Set Stokes image buttons
        WIDGET_CONTROL, (*(*info).ctrlscp).stokes_button_ids[i], SENSITIVE=stokes_enabled, $
          SET_BUTTON = (((*(*info).stokesparams).labels)[(*(*info).dataparams).s] EQ $
                         (*(*info).stokesparams).button_labels[i])
        ; Set Stokes detailed spectrum buttons
        WIDGET_CONTROL, (*(*info).ctrlscp).stokes_spbutton_ids[i], SET_BUTTON=stokes_sp_set, $
          SENSITIVE=(stokes_enabled AND ((total_sp_select GT 1) OR $
                                         (total_sp_select AND ABS(stokes_sp_set-1))))
      ENDFOR
			IF ((*(*info).dataparams).nlp EQ 1) THEN $
        WIDGET_CONTROL, (*(*info).ctrlscp).stokes_spbutton_ids[0], SET_BUTTON = 0
			; ==================== Displays Tab ====================
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
      WIDGET_CONTROL, (*(*info).ctrlscp).scaling_cbox, $
        SET_COMBOBOX_SELECT=(*(*info).scaling).imrefscaling
			CRISPEX_SCALING_SET_BOXBUTTONS, event
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
;			IF (*(*info).winswitch).showparam THEN BEGIN
;				(*(*info).winswitch).showparam = 0
;				CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE, event
;				WIDGET_CONTROL, (*(*info).winids).restsesfeedbtlb, /SHOW
;			ENDIF
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	IF ((*(*info).paths).opath_write EQ 1) THEN BEGIN
		WIDGET_CONTROl,/HOURGLASS
		PRINT,'Saving current loop points..'
		x_coords = *(*(*info).loopparams).xp		&	y_coords = *(*(*info).loopparams).yp
		x_loop_pts = *(*(*info).loopparams).xr		&	y_loop_pts = *(*(*info).loopparams).yr
		w_loop_pts = *(*(*info).loopparams).w_lpts	&	spect_pos = (*(*info).dataparams).lp
		t_saved = (*(*info).dispparams).t
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
		x_coords = *(*(*info).loopparams).xp	&	y_coords = *(*(*info).loopparams).yp
		x_loop_pts = *(*(*info).loopparams).xr			&	y_loop_pts = *(*(*info).loopparams).yr
		w_loop_pts = *(*(*info).loopparams).w_lpts		&	spect_pos = (*(*info).dataparams).lp
		t_saved = (*(*info).dispparams).t			&	loop_size= (*(*info).loopsdata).loopsize
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_EXACT_LOOPSLAB, event, /SAVE_SLICE
END

PRO CRISPEX_SAVE_EXACT_LOOPSLAB_CHECK, event, SAVE_SLICE=save_slice
; Opens a time estimate warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	x_coords = *(*(*info).loopparams).xp  &	y_coords = *(*(*info).loopparams).yp
	x_loop_pts = *(*(*info).loopparams).xr				&	y_loop_pts = *(*(*info).loopparams).yr
	spect_pos = (*(*info).dataparams).lp				&	loop_size= (*(*info).loopsdata).exact_loopsize
	spect_pos_low = lp_low						&	spect_pos_upp = lp_upp
	t_low = (*(*info).dispparams).t_low				&	t_upp = (*(*info).dispparams).t_upp
	t_saved = (*(*info).dispparams).t
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savswitch).cont = 0
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savswitch).cont],labels=['Saving procedure']
	WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO CRISPEX_SAVE_RETRIEVE_LOOPSLAB, event, SAVE_SLICE=save_slice
; Handles the saving of an exact (i.e. linearly interpolated) timeslab along a retrieved loop
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
		IF ((*(*info).dataparams).refimfilename NE '') THEN CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).imref_only, (*(*info).dataparams).imfilename, $
			refimfilename=(*(*info).dataparams).refimfilename, *(*(*info).data).imagedata, refdata=*(*(*info).data).refdata, (*(*info).paths).opath, i, (*(*info).retrparams).retrieve_filecount, filename, data, nonlp,$
			imref, whichdata, fstr=fstr, loopdet=loopdet ELSE $
			CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).imref_only, (*(*info).dataparams).imfilename, *(*(*info).data).imagedata, (*(*info).paths).opath, i, $
				(*(*info).retrparams).retrieve_filecount, filename, data, nonlp, imref, whichdata, fstr=fstr, loopdet=loopdet
		RESTORE,(*(*(*info).retrparams).retrieve_files)[(i MOD (*(*info).retrparams).retrieve_filecount)]
		*(*(*info).loopparams).xp = x_coords		&	*(*(*info).loopparams).yp = y_coords
		*(*(*info).loopparams).xr = x_loop_pts		&	*(*(*info).loopparams).yr = y_loop_pts
		*(*(*info).loopparams).w_lpts = w_loop_pts
		IF (N_ELEMENTS(T_SAVED) NE 1) THEN t_saved = (*(*info).dispparams).t
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
		x_coords = *(*(*info).loopparams).xp		&	y_coords = *(*(*info).loopparams).yp
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

PRO CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, imref_only, imfilename, refimfilename=refimfilename, imdata, refdata=refdata, opath, var, endoflist, outputfilename, outputdata, outputnonlp, outputimref, outputwhichdata, $
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
			CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=refimfilename, outfilename=outputfilename, ext='csav', import_id = datestamp+'_'+timestamp
		ENDIF ELSE CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=refimfilename, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
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
				CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=refimfilename, outfilename=outputfilename, ext='csav', import_id = datestamp+'_'+timestamp
			ENDIF ELSE CRISPEX_SAVE_DETERMINE_FILENAME, event, infilename=refimfilename, outfilename=outputfilename, ext='csav', import_id = 'D'+STRTRIM(index,2)
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
		IF ((*(*info).dataparams).refimfilename NE '') THEN CRISPEX_SAVE_RETRIEVE_DETERMINE_FILENAME, event, (*(*info).savswitch).det_imref_only, (*(*info).dataparams).imfilename, $
			refimfilename=(*(*info).dataparams).refimfilename, *(*(*info).data).imagedata, refdata=*(*(*info).data).refdata, (*(*info).paths).opath, i, (*(*info).detparams).nr_sel_loops, filename, data, nonlp, $
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	thispath = DIALOG_PICKFILE(TITLE = 'CRISPEX'+(*(*info).sesparams).instance_label+': Set input path', /DIRECTORY, PATH = (*(*info).paths).ipath)
	IF (thispath EQ '') THEN RETURN ELSE (*(*info).paths).ipath = thispath
	IF ((*(*info).winids).savewintlb GT 0) THEN WIDGET_CONTROL, (*(*info).ctrlssav).path_textlab, SET_VALUE = STRTRIM((*(*info).paths).ipath,2)
	IF (event.TOP EQ (*(*info).winids).preftlb) THEN (*(*info).prefs).tmp_prefipath = (*(*info).paths).ipath
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).paths).ipath],labels=['Input path']
END

PRO CRISPEX_SAVE_SET_OPATH, event 
; Sets the output path for all saving procedures
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).savparams).savpro EQ 'SESSION') THEN CRISPEX_SESSION_SAVE, event, session_filename
	IF ((STRPOS((*(*info).savparams).savpro,'FRAMES') NE -1) OR STRCMP((*(*info).savparams).savpro,'MPEG')) THEN CRISPEX_SAVE_FRAMES_SAVE, event, session_filename ELSE $
		IF (STRPOS((*(*info).savparams).savpro,'LINESCAN') NE -1) THEN CRISPEX_SAVE_LINESCAN_SAVE, event, session_filename
END

PRO CRISPEX_SAVE_WARNING_YESNO, event, warningmessage1, warningmessage2, warningmessage3, ok_event=ok_event, cancel_event=cancel_event
; Creates the loopslice/slab warning window
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'JPEG_FRAMES'
	(*(*info).savparams).snapshot = 1
	CRISPEX_SAVE_FRAMES, event
END	

PRO CRISPEX_SAVE_JPEG_ALL_FRAMES, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'JPEG_FRAMES'
	(*(*info).savparams).snapshot = 0
	CRISPEX_SAVE_FRAMES, event
END

PRO CRISPEX_SAVE_JPEG_LINESCAN, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'JPEG_LINESCAN'
	CRISPEX_SAVE_LINESCAN, event
END	

PRO CRISPEX_SAVE_PNG_SNAPSHOT, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'PNG_FRAMES'
	(*(*info).savparams).snapshot = 1
	CRISPEX_SAVE_FRAMES, event
END	

PRO CRISPEX_SAVE_PNG_ALL_FRAMES, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'PNG_FRAMES'
	(*(*info).savparams).snapshot = 0
	CRISPEX_SAVE_FRAMES, event
END

PRO CRISPEX_SAVE_PNG_LINESCAN, event
; Handles the saving of a single main image (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'PNG_LINESCAN'
	CRISPEX_SAVE_LINESCAN, event
END	

PRO CRISPEX_SAVE_CHECK, event
; Checks the jpeg series save filename for validity and overwrite problems
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
		t_id = 't'+STRING((*(*info).dispparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	label = ['series','snapshot']
	IF (*(*info).savparams).snapshot THEN BEGIN
		ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
		nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
		t_id = 't'+STRING((*(*info).dispparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
	nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
	t_before = (*(*info).dispparams).t
	IF (*(*info).savparams).snapshot THEN BEGIN
		t_low = (*(*info).dispparams).t		&	t_upp = t_low
	ENDIF ELSE BEGIN
		t_low = (*(*info).dispparams).t_low	&	t_upp = (*(*info).dispparams).t_upp 
		lp_id = 'lp'+STRING((*(*info).dataparams).lp,FORMAT='(I0'+STRTRIM(nlpos,2)+')')
	ENDELSE
	WIDGET_CONTROL, /HOURGLASS
	IF STRCMP((*(*info).savparams).savpro,'MPEG') THEN mpegid = MPEG_OPEN([(*(*info).winsizes).xywinx, (*(*info).winsizes).xywiny], QUALITY = 100, BITRATE = 1E8)
	FOR i = t_low, t_upp DO BEGIN
		(*(*info).dispparams).t = i
		CRISPEX_UPDATE_T, event
		IF (*(*info).savparams).overlays_incl THEN BEGIN
			CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
				NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
			image=TVRD()
		ENDIF ELSE CRISPEX_DRAW_SCALING, event, image, /MAIN
		TVLCT,r,g,b,/GET
		IF (*(*info).savparams).snapshot THEN midtension = '' ELSE BEGIN
			t_id = 't'+STRING((*(*info).dispparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
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
	(*(*info).dispparams).t = t_before
	CRISPEX_UPDATE_T, event
	IF (*(*info).savparams).overlays_incl THEN CRISPEX_DRAW_XY, event
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SAVE_LINESCAN, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename
	IF ((*(*info).savparams).savpro EQ 'JPEG_LINESCAN') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save JPEG line scan', standardfilename, 'CRISPEX_SAVE_CHECK'
	IF ((*(*info).savparams).savpro EQ 'PNG_LINESCAN') THEN CRISPEX_SAVE_GET_FILENAME, event, 'Save PNG line scan', standardfilename, 'CRISPEX_SAVE_CHECK'
END

PRO CRISPEX_SAVE_LINESCAN_SAVE, event, supplied_filename
; Handles the saving of a series (between temporal boundaries) of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	nlpos = CEIL(ALOG10((*(*info).dataparams).nlp))
	ntpos = CEIL(ALOG10((*(*info).dataparams).nt))
	lp_before = (*(*info).dataparams).lp
	WIDGET_CONTROL, /HOURGLASS
	t_id = 't'+STRING((*(*info).dispparams).t,FORMAT='(I0'+STRTRIM(ntpos,2)+')')
	FOR i = (*(*info).dispparams).lp_low, (*(*info).dispparams).lp_upp DO BEGIN
		(*(*info).dataparams).lp = i
		CRISPEX_UPDATE_T, event
		IF (*(*info).savparams).overlays_incl THEN BEGIN
			CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
				NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
			image=TVRD()
		ENDIF ELSE CRISPEX_DRAW_SCALING, event, image, /MAIN
		TVLCT,r,g,b,/GET
		IF (*(*info).savparams).linescan_ls THEN BEGIN
			CRISPEX_UPDATE_LP, event
      CRISPEX_DRAW_SPECTRAL_MAIN, event, /LS_ONLY
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
    CRISPEX_DRAW_SPECTRAL_MAIN, event, /LS_ONLY
	ENDIF
	WIDGET_CONTROL, (*(*info).winids).savewintlb, /DESTROY
	(*(*info).winids).savewintlb = 0
END

PRO CRISPEX_SAVE_MPEG, event
; Handles the saving of a series of main images (as in display) as JPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).savpro = 'MPEG'
	CRISPEX_SAVE_CHECK_PATH_WRITE, event
	CRISPEX_SAVE_DETERMINE_FILENAME, event, outfilename=standardfilename, import_id='lp'+STRTRIM(LONG((*(*info).dataparams).lp),2)
	CRISPEX_SAVE_GET_FILENAME, event, 'Save MPEG movie', standardfilename, 'CRISPEX_SAVE_CHECK'
END

PRO CRISPEX_SAVE_OPTIONS, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).savparams).overlays_incl THEN BEGIN
		CRISPEX_DRAW_XY, event, NO_CURSOR=ABS((*(*info).savparams).overlays_curs-1), NO_NUMBER=ABS((*(*info).savparams).overlays_num-1), THICK=(*(*info).savparams).overlays_thick, $
			NO_ENDPOINTS=ABS((*(*info).savparams).overlays_pts-1), SYMSIZE=(*(*info).savparams).overlays_symsize, ASECBAR=(*(*info).savparams).overlays_asecbar
	ENDIF ELSE CRISPEX_DRAW_XY, event
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_NUMBER, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_num = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_num],labels=['Number overlays']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_CURSOR, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_curs = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_curs],labels=['Include cursor']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_THICK, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_thick = event.VALUE
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_thick],labels=['Overlay thickness']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_SYMSIZE, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_symsize = event.VALUE
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_symsize],labels=['Overlay symbol size']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ENDPOINTS, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_pts = event.SELECT
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_pts],labels=['Include endpoints']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_asecbar = event.SELECT
	WIDGET_CONTROL, (*(*info).ctrlssav).overlays_asecbar_slider, SENSITIVE = (*(*info).savparams).overlays_asecbar
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_asecbar],labels=['Add arcseconds bar']
END

PRO CRISPEX_SAVE_OPTIONS_OVERLAYS_ASECBAR_LENGTH, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).overlays_asecbar_length = event.VALUE
	(*(*info).savparams).overlays_asecbar_pix = (*(*info).savparams).overlays_asecbar_length / FLOAT((*(*info).meas).arcsecpix) * (*(*info).zooming).factor
	CRISPEX_SAVE_OPTIONS_OVERLAYS_INCLUDE_UPDATE_DISPLAY, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).overlays_asecbar_length, (*(*info).savparams).overlays_asecbar_pix],$
		labels=['Arcseconds bar length','Arcseconds bar length in pixels']
END

PRO CRISPEX_SAVE_OPTIONS_INCLUDE_DETSPECT, event
; Handles the extra saving options for save as PNG/JPEG/MPEG
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).savparams).linescan_ls = event.SELECT
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).savparams).linescan_ls],labels=['Save detailed spectrum']
END

;================================================================================= SLIDER CONTROL PROCEDURES
PRO CRISPEX_SLIDER_NPHI, event
; Handles the change in spectral slit length slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).phiparams).nphi_set = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).phiparams).nphi_set],labels=['Phi-slit length']
  CRISPEX_UPDATE_PHISLIT_COORDS, event
	CRISPEX_UPDATE_PHISLICE, event
END

PRO CRISPEX_SLIDER_PHI_ANGLE, event
; Handles the change in spectral slit angle slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).phiparams).angle = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).phiparams).angle],labels=['Phi-slit angle']
  CRISPEX_UPDATE_PHISLIT_COORDS, event
	CRISPEX_PHISLIT_DIRECTION, event
	CRISPEX_UPDATE_PHISLICE, event
END

PRO CRISPEX_SLIDER_LP, event
; Handles the change in spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp = event.VALUE
  WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, $
    SENSITIVE=((*(*info).pbparams).lp_blink NE (*(*info).dataparams).lp)
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_DECR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp -= 1L 
	IF ((*(*info).dataparams).lp LT (*(*info).dispparams).lp_low) THEN $
    (*(*info).dataparams).lp = (*(*info).dispparams).lp_upp
	CRISPEX_SLIDER_LP_UPDATE, event, /OVERRIDE_DIAGNOSTIC
END

PRO CRISPEX_SLIDER_LP_INCR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp += 1L
	IF ((*(*info).dataparams).lp GT (*(*info).dispparams).lp_upp) THEN $
    (*(*info).dataparams).lp = (*(*info).dispparams).lp_low
	CRISPEX_SLIDER_LP_UPDATE, event, /OVERRIDE_DIAGNOSTIC
END

PRO CRISPEX_SLIDER_LP_REF, event
; Handles change in reference spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp_ref = event.VALUE
	CRISPEX_SLIDER_LP_UPDATE, event
END

PRO CRISPEX_SLIDER_LP_REF_DECR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp_ref -= 1L 
	IF ((*(*info).dataparams).lp_ref LT (*(*info).dispparams).lp_ref_low) THEN $
    (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_upp
	CRISPEX_SLIDER_LP_UPDATE, event, /OVERRIDE_DIAGNOSTIC
END

PRO CRISPEX_SLIDER_LP_REF_INCR, event
; Handles increase of spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).lp_ref += 1L
	IF ((*(*info).dataparams).lp_ref GT (*(*info).dispparams).lp_ref_upp) THEN $
    (*(*info).dataparams).lp_ref = (*(*info).dispparams).lp_ref_low
	CRISPEX_SLIDER_LP_UPDATE, event, /OVERRIDE_DIAGNOSTIC
END

PRO CRISPEX_SLIDER_LP_REF_LOCK, event, UNLOCK=unlock, NO_DRAW=no_draw
; Handles (un)locking the reference reference to (from) the main spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF KEYWORD_SET(UNLOCK) THEN BEGIN
	  (*(*info).ctrlsswitch).lp_ref_lock = 0
    WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_but, SET_BUTTON=0
  ENDIF ELSE $
	  (*(*info).ctrlsswitch).lp_ref_lock = event.SELECT
  IF (*(*info).ctrlsswitch).lp_ref_lock THEN $
	  (*(*info).dataparams).lp_ref = (*(*info).dataparams).lp
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, $
    SENSITIVE = ABS((*(*info).ctrlsswitch).lp_ref_lock-1), SET_VALUE = (*(*info).dataparams).lp_ref
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).lp, (*(*info).dataparams).lp_ref], labels=['lp','lp_ref']
  IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
  	CRISPEX_UPDATE_T, event
  	CRISPEX_UPDATE_LP, event
  	CRISPEX_DRAW, event
  ENDIF
END

PRO CRISPEX_SLIDER_LP_UPDATE, event, OVERRIDE_DIAGNOSTIC=override_diagnostic, $
  NO_DRAW=no_draw
; Handles the the update after change in the spectral position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Determine whether lp falls in not displayed diagnostic window, act accordingly
  lp_diag_all = TOTAL((*(*info).dataparams).lp GE (*(*info).intparams).diag_start)-1
  IF ((*(*info).intparams).disp_diagnostics[lp_diag_all] EQ 0) THEN BEGIN
    ; Determine distance to lower and upper boundaries of current diagnostic window
    dist_low = (*(*info).dataparams).lp - (*(*info).intparams).diag_start[lp_diag_all]
    dist_upp = ABS((*(*info).dataparams).lp - ((*(*info).intparams).diag_start[lp_diag_all]+$
      (*(*info).intparams).diag_width[lp_diag_all]-1))
    IF (((dist_low LT dist_upp) AND KEYWORD_SET(OVERRIDE_DIAGNOSTIC)) OR $
        ((dist_low GE dist_upp) AND ~KEYWORD_SET(OVERRIDE_DIAGNOSTIC))) THEN BEGIN
      (*(*info).dataparams).lp = (*(*info).intparams).diag_start[lp_diag_all+1]
      lp_diag_all += 1
    ENDIF ELSE BEGIN
      (*(*info).dataparams).lp = (*(*info).intparams).diag_start[lp_diag_all]-1 
      lp_diag_all -= 1
    ENDELSE
  ENDIF
  (*(*info).intparams).lp_diag_all = lp_diag_all
	IF (*(*info).ctrlsswitch).lp_ref_lock THEN $
    (*(*info).dataparams).lp_ref = (*(*info).dataparams).lp
  ; Determine whether lp_ref falls in not displayed diagnostic window, act accordingly
  lp_ref_diag_all = TOTAL((*(*info).dataparams).lp_ref GE (*(*info).intparams).refdiag_start)-1
  IF ((*(*info).intparams).disp_refdiagnostics[lp_ref_diag_all] EQ 0) THEN BEGIN
    ; Determine distance to lower and upper boundaries of current diagnostic window
    dist_low = (*(*info).dataparams).lp_ref - (*(*info).intparams).refdiag_start[lp_ref_diag_all]
    dist_upp = ABS((*(*info).dataparams).lp_ref - ((*(*info).intparams).refdiag_start[lp_ref_diag_all]+$
      (*(*info).intparams).refdiag_width[lp_ref_diag_all]-1))
    IF (dist_low LT dist_upp) THEN BEGIN
      (*(*info).dataparams).lp_ref = (*(*info).intparams).refdiag_start[lp_ref_diag_all+1]
      lp_ref_diag_all += 1
    ENDIF ELSE BEGIN
      (*(*info).dataparams).lp_ref = (*(*info).intparams).refdiag_start[lp_ref_diag_all]-1 
      lp_ref_diag_all -= 1
    ENDELSE
  ENDIF
  (*(*info).intparams).lp_ref_diag_all = lp_ref_diag_all
  (*(*info).scaling).idx = $
    (((*(*info).scaling).imrefscaling EQ 0) OR ((*(*info).scaling).imrefscaling EQ 2)) * $
    (*(*info).intparams).lp_diag_all + $
    (((*(*info).scaling).imrefscaling GT 0) + ((*(*info).scaling).imrefscaling GT 2)) * $
    (*(*info).intparams).ndiagnostics + $
    ((*(*info).scaling).imrefscaling EQ 1) * (*(*info).intparams).lp_ref_diag_all + $
    ((*(*info).scaling).imrefscaling GT 1) * (*(*info).intparams).nrefdiagnostics 
	WIDGET_CONTROL, (*(*info).ctrlscp).lp_slider, SET_VALUE = (*(*info).dataparams).lp
  WIDGET_CONTROL, (*(*info).ctrlscp).lp_ref_slider, SET_VALUE = (*(*info).dataparams).lp_ref
  IF (((*(*info).intparams).ndiagnostics GT 1) OR $
    ((*(*info).intparams).nrefdiagnostics GT 1)) THEN BEGIN
    CRISPEX_SCALING_SET_BOXBUTTONS, event
    CRISPEX_SCALING_SET_SLIDERS, event
  ENDIF
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).lp, (*(*info).dataparams).lp_ref], $
      labels=['lp','lp_ref']
  IF ~KEYWORD_SET(NO_DRAW) THEN BEGIN
  	CRISPEX_UPDATE_T, event
  	CRISPEX_UPDATE_LP, event
;    IF (*(*info).winswitch).showsp THEN CRISPEX_UPDATE_SPSLICE, event
;    IF (*(*info).winswitch).showrefsp THEN CRISPEX_UPDATE_REFSPSLICE, event
  	CRISPEX_DRAW, event
  ENDIF
END

PRO CRISPEX_SLIDER_SPECTBLINK, event
; Handles the change in spectral step (for spectral blink) slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).lp_blink = event.VALUE
  WIDGET_CONTROL, (*(*info).ctrlscp).lp_blink_but, $
    SENSITIVE=((*(*info).pbparams).lp_blink NE (*(*info).dataparams).lp)
  (*(*info).dataparams).lp = (*(*info).pbparams).lp_blink
  CRISPEX_SLIDER_LP_UPDATE, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).lp_blink], labels=['lp_blink']
END

PRO CRISPEX_SLIDER_SPEED, event
; Handles the change in playback speed slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN $
    CRISPEX_VERBOSE_GET_ROUTINE, event, /IGNORE_LAST
	(*(*info).pbparams).t_speed = event.VALUE
	IF (*(*info).pbparams).spmode THEN WIDGET_CONTROL, (*(*info).ctrlscp).t_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed ELSE $
		WIDGET_CONTROL, (*(*info).ctrlscp).lp_speed_slider, SET_VALUE = (*(*info).pbparams).t_speed
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).t_speed], labels=['Playback (blink) speed']
END

PRO CRISPEX_SLIDER_STEP, event
; Handles the change in temporal step slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).pbparams).t_step = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).pbparams).t_step], labels=['t_step']
END

PRO CRISPEX_SLIDER_T, event
; Handles the change in temporal slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dispparams).t = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).dispparams).t], labels=['t']
	CRISPEX_UPDATE_T, event
	IF (*(*info).winswitch).showphis THEN BEGIN 
		IF ((*(*info).dispparams).phislice_update OR (event.DRAG EQ 0)) THEN $
      CRISPEX_UPDATE_SLICES, event, /NO_DRAW $
    ELSE BEGIN		
			IF (*(*info).dataswitch).onecube THEN $
        WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1, $
          SET_VALUE = 'Update spectral windows' ;$
     ; ELSE WIDGET_CONTROL, (*(*info).ctrlscp).slice_button, SENSITIVE = 1
		ENDELSE
	ENDIF 
  CRISPEX_DRAW, event, NO_PHIS=event.DRAG
END

PRO CRISPEX_SLIDER_TIME_OFFSET, event
; Handles the change in temporal slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  CASE (*(*info).dispparams).master_time OF
    0:  (*(*info).dispparams).toffset_main = event.VALUE
    1:  (*(*info).dispparams).toffset_ref = event.VALUE
  ENDCASE
  CRISPEX_COORDS_TRANSFORM_T, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [event.VALUE], labels=['Raster timing offset']
	CRISPEX_UPDATE_T, event
  IF (*(*info).winswitch).showsp THEN CRISPEX_DISPLAYS_SP_REPLOT_AXES, event
  IF (*(*info).winswitch).showrefsp THEN CRISPEX_DISPLAYS_REFSP_REPLOT_AXES, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_X, event
; Handles the change in cursor x-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).x = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).x], labels=['x']
	CRISPEX_UPDATE_SX, event
  IF (*(*info).winswitch).showsp THEN CRISPEX_UPDATE_SPSLICE, event
  IF (*(*info).winswitch).showrefsp THEN CRISPEX_UPDATE_REFSPSLICE, event
	IF (*(*info).winswitch).showphis THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
    CRISPEX_UPDATE_PHISLIT_COORDS, event
		CRISPEX_UPDATE_PHISLICE, event
	ENDIF ELSE CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_XPOS, event
; Handles change in x-slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).zooming).xpos = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).xpos], labels=['xpos']
	CRISPEX_UPDATE_SX, event
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_Y, event
; Handles the change in cursor y-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).dataparams).y = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN CRISPEX_VERBOSE_GET, event, [(*(*info).dataparams).y], labels=['y']
	CRISPEX_UPDATE_SY, event
  IF (*(*info).winswitch).showsp THEN CRISPEX_UPDATE_SPSLICE, event
  IF (*(*info).winswitch).showrefsp THEN CRISPEX_UPDATE_REFSPSLICE, event
	IF (*(*info).winswitch).showphis THEN BEGIN
		CRISPEX_PHISLIT_DIRECTION, event
    CRISPEX_UPDATE_PHISLIT_COORDS, event
		CRISPEX_UPDATE_PHISLICE, event
	ENDIF ELSE CRISPEX_DRAW, event
END

PRO CRISPEX_SLIDER_YPOS, event
; Handles change in y-slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*info).zooming).ypos = event.VALUE
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).ypos], labels=['ypos']
	CRISPEX_UPDATE_SY, event
	CRISPEX_UPDATE_T, event
	CRISPEX_DRAW, event
END

;========================= UPDATE SLICES AND PARAMETERS PROCEDURES
PRO CRISPEX_UPDATE_SLICES, event, NO_DRAW=no_draw
; Gets the new spectral phi slit scan for update of the spectral phi slit slice after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF ((*(*info).winswitch).showphis OR $
      ((*(*info).dataswitch).onecube AND (*(*info).winswitch).showls)) THEN BEGIN
		WIDGET_CONTROL,/HOURGLASS
		IF ((*(*info).dataparams).nt GT 1) THEN $
      *(*(*info).data).sspscan = (*(*(*info).data).scan)[(*(*info).dispparams).t_main] $
    ELSE $
      *(*(*info).data).sspscan = (*(*(*info).data).scan)
    IF ((*(*info).winswitch).showrefls AND ((*(*info).dataswitch).refspfile EQ 0)) THEN BEGIN
  		IF ((*(*info).dataparams).refnt GT 1) THEN $
        *(*(*info).data).refsspscan = (*(*(*info).data).refscan)[(*(*info).dispparams).t_ref] 
    ENDIF
    IF (*(*info).winswitch).showphis THEN BEGIN
  		*(*(*info).data).phiscan = (*(*(*info).data).sspscan)[*,*,$
         ((*(*info).dataparams).s * (*(*info).dataparams).nlp):$
        (((*(*info).dataparams).s+1)*(*(*info).dataparams).nlp-1)] 
      CRISPEX_UPDATE_PHISLICE, event, NO_DRAW=no_draw
    ENDIF
	ENDIF
END

PRO CRISPEX_UPDATE_SPSLICE, event
; Handles updating spectrum-time diagram
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  t_low = (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_low]
  t_upp = (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t_upp]
	spslice = ((*(*(*info).data).spdata)[ FIX((*(*info).dataparams).y) * (*(*info).dataparams).nx * $
            (*(*info).dataparams).ns + FIX((*(*info).dataparams).x) * (*(*info).dataparams).ns +$
		        (*(*info).dataparams).s ])
	IF (*(*info).dispswitch).warpspslice THEN $                   ; Warp slice if non-equidistant lp
    dispslice = WARP_TRI( (*(*info).dispparams).xo,(*(*info).dispparams).yo,$
                          (*(*info).dispparams).xi,(*(*info).dispparams).yi,$
                          spslice) $
  ELSE $
    dispslice = spslice
  dispslice = dispslice[(*(*info).dispparams).lp_low:(*(*info).dispparams).lp_upp,t_low:t_upp]
  FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
    IF ((*(*info).intparams).ndiagnostics GT 1) THEN $
      tmp_disp = dispslice[((*(*(*info).intparams).diag_starts)[d]-$
                            (*(*(*info).intparams).diag_starts)[0]):$
                            ((*(*(*info).intparams).diag_starts)[d]-$
                            (*(*(*info).intparams).diag_starts)[0]+$
                            (*(*(*info).intparams).diag_widths)[d]-1),*] $
    ELSE $
      tmp_disp = dispslice
  	IF ((*(*info).dispparams).slices_imscale EQ 0) THEN BEGIN
      minmax = CRISPEX_SCALING_SLICES(tmp_disp, $
        (*(*info).scaling).gamma[(*(*(*info).intparams).wheredispdiag)[d]], $
        (*(*info).scaling).histo_opt_val[(*(*(*info).intparams).wheredispdiag)[d]],$
        (*(*info).scaling).minimum[(*(*(*info).intparams).wheredispdiag)[d]],$
        (*(*info).scaling).maximum[(*(*(*info).intparams).wheredispdiag)[d]])
        (*(*info).scaling).spslice_min[(*(*(*info).intparams).wheredispdiag)[d]] = minmax[0]
        (*(*info).scaling).spslice_max[(*(*(*info).intparams).wheredispdiag)[d]] = minmax[1]
    ENDIF
    IF (d EQ 0) THEN $
      final_disp = tmp_disp $
    ELSE $
      final_disp = [final_disp,tmp_disp]
  ENDFOR
  *(*(*info).data).spslice = final_disp
END

PRO CRISPEX_UPDATE_REFSPSLICE, event
; Handles updating reference spectrum-time diagram
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  t_low = (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).t_low]
  t_upp = (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).t_upp]
	refspslice = ((*(*(*info).data).refspdata)[ FIX((*(*info).dataparams).y) * $
                (*(*info).dataparams).nx + FIX((*(*info).dataparams).x) ])
	IF (*(*info).dispswitch).warprefspslice THEN $                ; Warp slice if non-equidistant lp
    dispslice = WARP_TRI( (*(*info).dispparams).xo_ref,(*(*info).dispparams).yo_ref,$
                          (*(*info).dispparams).xi_ref,(*(*info).dispparams).yi_ref,$
		                      refspslice) $
	ELSE $
    dispslice = refspslice
  dispslice = dispslice[(*(*info).dispparams).lp_ref_low:(*(*info).dispparams).lp_ref_upp,t_low:t_upp]
  FOR d=0,(*(*info).intparams).ndisp_refdiagnostics-1 DO BEGIN
    sel_idx = (*(*info).intparams).ndiagnostics+(*(*(*info).intparams).wheredisprefdiag)[d]
    IF ((*(*info).intparams).nrefdiagnostics GT 1) THEN $
      tmp_disp = dispslice[((*(*(*info).intparams).refdiag_starts)[d]-$
                            (*(*(*info).intparams).refdiag_starts)[0]):$
                            ((*(*(*info).intparams).refdiag_starts)[d]-$
                            (*(*(*info).intparams).refdiag_starts)[0]+$
                            (*(*(*info).intparams).refdiag_widths)[d]-1),*] $
    ELSE $
      tmp_disp = dispslice
    IF ((*(*info).dispparams).slices_imscale EQ 0) THEN BEGIN
      refminmax = CRISPEX_SCALING_SLICES(tmp_disp, $
        (*(*info).scaling).gamma[sel_idx], $
        (*(*info).scaling).histo_opt_val[sel_idx], $
        (*(*info).scaling).minimum[sel_idx],(*(*info).scaling).maximum[sel_idx])
        (*(*info).scaling).spslice_min[sel_idx] = refminmax[0]
        (*(*info).scaling).spslice_max[sel_idx] = refminmax[1]
    ENDIF
    IF (d EQ 0) THEN $
      final_disp = tmp_disp $
    ELSE $
      final_disp = [final_disp,tmp_disp]
  ENDFOR
  *(*(*info).data).refspslice = final_disp
END

PRO CRISPEX_UPDATE_PHISLIT_COORDS, event
; Handles the update of the slit coordinates 
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  nw_prev = (*(*info).phiparams).nw_cur
  ; Failsafe necessary because of numerical accuracy
  cosval = COS(!DTOR * (*(*info).phiparams).angle)
  sinval = SIN(!DTOR * (*(*info).phiparams).angle)
  IF (((*(*info).phiparams).angle EQ 0) OR ((*(*info).phiparams).angle EQ 90)) THEN BEGIN
    cosval = ROUND(cosval)
    sinval = ROUND(sinval)
  ENDIF
  x_pts = cosval * (FINDGEN( 2*(*(*info).phiparams).nphi_set/2) - $
              (*(*info).phiparams).nphi_set/2) + LONG((*(*info).dataparams).x)
  y_pts = sinval * (FINDGEN( 2*(*(*info).phiparams).nphi_set/2) - $
              (*(*info).phiparams).nphi_set/2) + LONG((*(*info).dataparams).y)
	w = WHERE((x_pts GE 0) AND (x_pts LE (*(*info).dispparams).x_last) AND $
            (y_pts GE 0) AND (y_pts LE (*(*info).dispparams).y_last), nw)
	(*(*info).phiparams).nw_cur = nw
	*(*(*info).phiparams).x_pts = REBIN(x_pts[w], nw, (*(*info).dataparams).nlp)
	*(*(*info).phiparams).y_pts = REBIN(y_pts[w], nw, (*(*info).dataparams).nlp)
  midphi = WHERE(((*(*(*info).phiparams).x_pts)[*,0] EQ LONG((*(*info).dataparams).x)) AND $
                  ((*(*(*info).phiparams).y_pts)[*,0] EQ LONG((*(*info).dataparams).y)))
	(*(*info).phiparams).sphi = midphi +  ((*(*info).phiparams).nphi - nw)/2
  ; Replot axes only under conditions specified below
  IF (((*(*info).phiparams).nw_cur NE nw_prev) AND $                        ; change in nw_cur 
    (((*(*info).phiparams).nw_cur NE (*(*info).phiparams).nphi_set) OR $  ; while nw_cur lt nphi_set
    (nw_prev NE (*(*info).phiparams).nphi_set))) THEN $                   ; nw_cur becomes nphi_set
      CRISPEX_DISPLAYS_PHIS_REPLOT_AXES, event
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [nw,midphi], labels=['nw','midphi']
END

PRO CRISPEX_UPDATE_PHISLICE, event, NO_DRAW=no_draw
; Handles the actual update of the spectral phi slit slice after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	lp_pts = REBIN(FINDGEN(1,(*(*info).dataparams).nlp), (*(*info).phiparams).nw_cur, $
    (*(*info).dataparams).nlp)
	tmp = INTERPOLATE( *(*(*info).data).phiscan, *(*(*info).phiparams).x_pts, $
    *(*(*info).phiparams).y_pts, lp_pts) ;$
	phislice = (TRANSPOSE(tmp, [1,0]))[(*(*info).dispparams).lp_low:(*(*info).dispparams).lp_upp,*]
  ; If 2D, then display
  IF ((SIZE(phislice))[0] EQ 2) THEN BEGIN
    ; Warp slice if non-equidistant lp
  	IF (*(*info).dispswitch).warpspslice THEN BEGIN
      dispslice = INTERPOLATE(phislice,$
        (*(*(*info).dispparams).phisxtri)[0:((SIZE(phislice))[1]-1),0:((SIZE(phislice))[2]-1)],$
        (*(*(*info).dispparams).phisytri)[0:((SIZE(phislice))[1]-1),0:((SIZE(phislice))[2]-1)])
    ENDIF ELSE $
      dispslice = phislice
    FOR d=0,(*(*info).intparams).ndisp_diagnostics-1 DO BEGIN
      IF ((*(*info).intparams).ndiagnostics GT 1) THEN $
        tmp_disp = dispslice[((*(*(*info).intparams).diag_starts)[d]-$
                              (*(*(*info).intparams).diag_starts)[0]):$
                              ((*(*(*info).intparams).diag_starts)[d]-$
                              (*(*(*info).intparams).diag_starts)[0]+$
                              (*(*(*info).intparams).diag_widths)[d]-1),*] $
      ELSE $
        tmp_disp = dispslice
    	IF ((*(*info).dispparams).slices_imscale EQ 0) THEN BEGIN
        minmax = CRISPEX_SCALING_SLICES(tmp_disp, $
          (*(*info).scaling).gamma[(*(*(*info).intparams).wheredispdiag)[d]], $
          (*(*info).scaling).histo_opt_val[(*(*(*info).intparams).wheredispdiag)[d]],$
          (*(*info).scaling).minimum[(*(*(*info).intparams).wheredispdiag)[d]],$
          (*(*info).scaling).maximum[(*(*(*info).intparams).wheredispdiag)[d]])
        (*(*info).scaling).phislice_min[(*(*(*info).intparams).wheredispdiag)[d]] = minmax[0]
        (*(*info).scaling).phislice_max[(*(*(*info).intparams).wheredispdiag)[d]] = minmax[1]
      ENDIF
      IF (d EQ 0) THEN $
        final_disp = tmp_disp $
      ELSE $
        final_disp = [final_disp,tmp_disp]
    ENDFOR
    *(*(*info).data).phislice = final_disp
  ENDIF ELSE *(*(*info).data).phislice = phislice
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_UPDATE_LP, event
; Handles the update of displayed data after change in spectral position
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).overlayswitch).loopslit AND ((*(*info).loopparams).np GT 0) THEN BEGIN
    CRISPEX_LOOP_GET, event
    IF (*(*info).winswitch).showloop THEN BEGIN
  		IF (*(*info).dispswitch).exts THEN $
        *(*(*info).loopsdata).loopslice = REFORM((*(*(*info).loopsdata).loopslab)[*,*,$
          (*(*info).dataparams).lp-(*(*info).dispparams).lp_low]) $
      ELSE $
  			*(*(*info).loopsdata).loopslice = REFORM((*(*(*info).loopsdata).loopslab)[*,*,$
          (*(*info).dataparams).lp])
    ENDIF
		IF (*(*info).winswitch).showrefloop THEN BEGIN
			IF (*(*info).dispswitch).refexts THEN $
        *(*(*info).loopsdata).refloopslice = REFORM((*(*(*info).loopsdata).refloopslab)[*,*,$
          (*(*info).dataparams).lp_ref-(*(*info).dispparams).lp_ref_low]) $
      ELSE $
				*(*(*info).loopsdata).refloopslice = REFORM((*(*(*info).loopsdata).refloopslab)[*,*,$
          (*(*info).dataparams).lp_ref])
		ENDIF
	ENDIF
	IF (*(*info).winswitch).showrestloop THEN $
    FOR k=0,N_ELEMENTS(*(*(*info).restoreparams).disp_loopnr)-1 DO BEGIN
		IF (SIZE(*(*(*(*info).loopsdata).rest_loopslab[k]),/N_DIMENSIONS) EQ 3) THEN BEGIN
			IF (*(*(*info).restoreparams).disp_imref)[k] THEN $
        *(*(*(*info).loopsdata).rest_loopslice[k]) = $
          REFORM((*(*(*(*info).loopsdata).rest_loopslab[k]))[*,*,(*(*info).dataparams).lp_ref-$
				    (*(*info).dispparams).lp_ref_low]) $
       ELSE $
        *(*(*(*info).loopsdata).rest_loopslice[k]) = $
          REFORM((*(*(*(*info).loopsdata).rest_loopslab[k]))[*,*,(*(*info).dataparams).lp-$
            (*(*info).dispparams).lp_low])
		ENDIF ELSE BEGIN
			IF (*(*(*info).restoreparams).disp_imref)[k] THEN $
        *(*(*(*info).loopsdata).rest_loopslice[k]) = *(*(*(*info).loopsdata).rest_loopslab[k]) $
      ELSE $
				*(*(*(*info).loopsdata).rest_loopslice[k]) = *(*(*(*info).loopsdata).rest_loopslab[k])
		ENDELSE
	ENDFOR
	IF (*(*info).winswitch).showretrdet THEN $
    *(*(*info).loopsdata).det_loopslice = REFORM((*(*(*info).loopsdata).det_loopslab)[*,*,$
      (*(*info).dataparams).lp-(*(*info).dispparams).lp_low])
  CRISPEX_SCALING_APPLY_SELECTED, event
END

PRO CRISPEX_UPDATE_T, event
; Handles the updated of displayed data after change in framenumber
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  ; Update time indices
  (*(*info).dispparams).t_main = (*(*(*info).dispparams).tsel_main)[(*(*info).dispparams).t]
  IF ((*(*info).dataparams).refnt GT 1) THEN $
    (*(*info).dispparams).t_ref = (*(*(*info).dispparams).tsel_ref)[(*(*info).dispparams).t]
  IF ((*(*info).dataparams).sjint GT 1) THEN $
    (*(*info).dispparams).t_sji = (*(*(*info).dispparams).tsel_sji)[(*(*info).dispparams).t]
	IF (*(*info).winswitch).showdop THEN BEGIN
		(*(*info).dataparams).lp_dop = 2*(*(*info).dataparams).lc - (*(*info).dataparams).lp
		(*(*info).dispswitch).drawdop = (((*(*info).dataparams).lp_dop GE (*(*info).dispparams).lp_low)$
                                 AND ((*(*info).dataparams).lp_dop LE (*(*info).dispparams).lp_upp)$
                                 AND ((*(*info).dataparams).lp_dop NE (*(*info).dataparams).lc)) 
	ENDIF
  ; Determine main image, in case the cube has a spectral dimension
	IF ((*(*info).dataswitch).spfile EQ 1) OR (*(*info).dataswitch).onecube THEN BEGIN
    basecubeidx = (*(*info).dispparams).t_main * (*(*info).dataparams).nlp * (*(*info).dataparams).ns+ $
			            (*(*info).dataparams).s * (*(*info).dataparams).nlp
		*(*(*info).data).xyslice = (*(*(*info).data).imagedata)[$
      basecubeidx + (*(*info).dataparams).lp]
		IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN $
			temp_xyslice = REFORM( ((*(*(*info).data).imagedata)[$
        basecubeidx + (*(*info).dataparams).lp_dop]))
  ; Determine main image, in case the cube has no spectral dimension
	ENDIF ELSE BEGIN
    basecubeidx = (*(*info).dataparams).s * (*(*info).dataparams).nlp 
		*(*(*info).data).xyslice = REFORM( ((*(*(*info).data).imagedata)[basecubeidx + $
      (*(*info).dataparams).lp]))
		IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN $
      temp_xyslice = REFORM( ((*(*(*info).data).imagedata)[basecubeidx + $
      (*(*info).dataparams).lp_dop]))
	ENDELSE
  ; Determine Doppler image
	IF ((*(*info).winswitch).showdop AND (*(*info).dispswitch).drawdop) THEN BEGIN
		IF ((*(*info).dataparams).lp_dop GT (*(*info).dataparams).lc) THEN $
      *(*(*info).data).dopslice = temp_xyslice - *(*(*info).data).xyslice $
    ELSE $
      *(*(*info).data).dopslice = *(*(*info).data).xyslice - temp_xyslice 
	ENDIF
  ; Determine reference image
	IF ((*(*info).winswitch).showref OR (*(*info).winswitch).showimref) THEN BEGIN
		IF (((*(*info).dataparams).refnlp GT 1) AND ((*(*info).dataparams).refnt GT 1)) THEN BEGIN
      refidx = (*(*info).dispparams).t_ref * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref
      *(*(*info).data).refslice = REFORM( ((*(*(*info).data).refdata)[refidx]))
		ENDIF ELSE BEGIN
			IF ((*(*info).dataparams).refnt EQ 0) THEN $
        *(*(*info).data).refslice = (*(*(*info).data).refdata) $
			ELSE IF ((*(*info).dataparams).refnt EQ 1) THEN BEGIN
				IF ((*(*info).dataparams).refnlp NE 1) THEN $
          *(*(*info).data).refslice = REFORM( ((*(*(*info).data).refdata)[$
            (*(*info).dataparams).lp_ref])) $
				ELSE $
          *(*(*info).data).refslice = REFORM( ((*(*(*info).data).refdata)[0]))
			ENDIF ELSE IF ((*(*info).dataparams).refnt EQ (*(*info).dataparams).nt) THEN BEGIN
        *(*(*info).data).refslice = REFORM( ((*(*(*info).data).refdata)[$
          (*(*info).dispparams).t_ref]))
			ENDIF ELSE $
        *(*(*info).data).refslice = REFORM( ((*(*(*info).data).refdata)[$
          (*(*info).dispparams).t_ref * (*(*info).dataparams).refnlp + (*(*info).dataparams).lp_ref]))
		ENDELSE
	ENDIF
  ; Determine mask image
	IF (*(*info).dataswitch).maskfile THEN BEGIN
		IF ((*(*info).dataparams).masknt GT 1) THEN BEGIN
      *(*(*info).data).maskslice = REFORM( ((*(*(*info).data).maskdata)[$
        (*(*info).dispparams).t_main]))
		ENDIF ELSE BEGIN
      *(*(*info).data).maskslice = REFORM( ((*(*(*info).data).maskdata)[0]))
		ENDELSE
	ENDIF
  ; Determine sji image
	IF (*(*info).dataswitch).sjifile THEN BEGIN
		IF ((*(*info).dataparams).sjint GT 1) THEN $
        *(*(*info).data).sjislice = ((*(*(*info).data).sjidata)[(*(*info).dispparams).t_sji]) $
    ELSE $
        *(*(*info).data).sjislice = ((*(*(*info).data).sjidata)[0])
	ENDIF
END

PRO CRISPEX_UPDATE_SX, event
; Handles the change in xy- and reference image x-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	sx = ((*(*info).dataparams).x - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
        ((*(*info).dataparams).d_nx+1)
	sxp = (*(*(*info).loopparams).xp - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
        ((*(*info).dataparams).d_nx+1)
	sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * (*(*info).winsizes).xywinx / $
        ((*(*info).dataparams).d_nx+1)
	sxsji = (*(*info).dataparams).xsji * (*(*info).winsizes).sjiwinx / $
        ((*(*info).dataparams).sjinx+1)
	(*(*info).curs).sxlock = sx 
	(*(*info).curs).sx = (*(*info).curs).sxlock
	(*(*info).curs).sxsji = sxsji
	*(*(*info).overlayparams).sxp = sxp 
	*(*(*info).overlayparams).sxr = sxr 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET,event,[(*(*info).dataparams).x,(*(*info).curs).sxlock,(*(*info).curs).sx],$
      labels=['x','sxlock','sx']
END

PRO CRISPEX_UPDATE_SY, event
; Handles the change in xy- and reference image y-position slider
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	sy = ((*(*info).dataparams).y - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
        ((*(*info).dataparams).d_ny+1)
	syp = (*(*(*info).loopparams).yp - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
        ((*(*info).dataparams).d_ny+1)
	syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * (*(*info).winsizes).xywiny / $
        ((*(*info).dataparams).d_ny+1)
	sysji = (*(*info).dataparams).ysji * (*(*info).winsizes).sjiwiny / $
        ((*(*info).dataparams).sjiny+1)
	(*(*info).curs).sylock = sy 
	(*(*info).curs).sy = (*(*info).curs).sylock
	(*(*info).curs).sysji = sysji
	*(*(*info).overlayparams).syp = syp 
	*(*(*info).overlayparams).syr = syr 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET,event,[(*(*info).dataparams).y,(*(*info).curs).sylock,(*(*info).curs).sy],$
      labels=['y','sylock','sy']
END

PRO CRISPEX_UPDATE_USER_FEEDBACK, event, title=title, var=var, minvar=minvar, maxvar=maxvar, feedback_text=feedback_text, destroy_top=destroy_top, close_button=close_button, session=session
; Handles the update of user feedback while saving timeslices
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	FOR i=0,N_ELEMENTS(feedback_text)-1 DO XYOUTS, xout[i], yout[i], feedback_text[i], COLOR=255, $
                                                 /DEVICE, CHARSIZE=1.125
END

PRO CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, msg, OVERWRITEABLE=overwriteable, DONE=done, OPT=opt, $
                                           WIDGET=widget, NEWLINE=newline, FAILED=failed, $
                                           REPEAT_STAGE=repeat_stage, NO_ROUTINE=no_routine, $
                                           WARNING=warning, ERROR=error
  routine = 'CRISPEX SETUP: '
  stages = [KEYWORD_SET(OPT),KEYWORD_SET(WIDGET),KEYWORD_SET(WARNING),KEYWORD_SET(ERROR)]
  stage_set = (WHERE(stages EQ 1))[0]+1
  stage_lab = ['','Setting start-up options ','Setting up widget ','WARNING: ','ERROR: ']
  IF KEYWORD_SET(NO_ROUTINE) THEN routine = STRJOIN(REPLICATE(' ',STRLEN(routine)))
  IF KEYWORD_SET(FAILED) THEN donetext = 'failed.' ELSE donetext = 'done!'
  IF KEYWORD_SET(REPEAT_STAGE) THEN outputtext = routine+stage_lab[stage_set]+msg+'... ' $
    ELSE outputtext = ''
  IF KEYWORD_SET(OVERWRITEABLE) THEN BEGIN
    IF KEYWORD_SET(DONE) THEN BEGIN
	    WRITEU,-1,outputtext+donetext
      PRINT,''
    ENDIF ELSE WRITEU,-1,routine+stage_lab[stage_set]+msg+'... '
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(NEWLINE) THEN PRINT,''
    FOR i=0,N_ELEMENTS(msg)-1 DO BEGIN
      IF (i NE 0) THEN temp_routine = STRJOIN(REPLICATE(' ',STRLEN(routine))) $
        ELSE temp_routine = routine
      PRINT, temp_routine+stage_lab[stage_set]+msg[i]
    ENDFOR
  ENDELSE
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
  rname = (SCOPE_TRACEBACK(/STRUCTURE))[-2].ROUTINE
	IF KEYWORD_SET(IGNORE_LAST) THEN (*(*info).feedbparams).last_routine = ''
	IF ((rname NE (*(*info).feedbparams).last_routine) AND $
    ((*(*info).feedbparams).last_routine_count GT 0)) THEN PRINT,''
	IF (rname EQ (*(*info).feedbparams).last_routine) THEN $
    (*(*info).feedbparams).last_routine_count += 1 $
  ELSE $
    (*(*info).feedbparams).last_routine_count = 0
	IF ((*(*info).feedbparams).last_routine_count GT 0) THEN $
    rcount = ' x '+STRTRIM((*(*info).feedbparams).last_routine_count,2)+'.' $
  ELSE $
    rcount = '.'
	IF (rname NE (*(*info).feedbparams).last_routine) THEN $
    PRINT,prespace+'CRISPEX RUN: Called '+rname+'.' $
  ELSE $
		WRITEU,-1,STRING(FORMAT='(%"\r'+prespace+'CRISPEX RUN: Called ",a'+$
      STRTRIM(STRLEN(rname),2)+',a'+STRTRIM(STRLEN(rcount),2)+')',rname,rcount) 
	(*(*info).feedbparams).last_routine = rname
END

PRO CRISPEX_VERBOSE_SET, event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[0], SET_BUTTON = (TOTAL((*(*info).feedbparams).verbosity) EQ 0)
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[1], SET_BUTTON = ((*(*info).feedbparams).verbosity)[2]
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[2], SET_BUTTON = ((*(*info).feedbparams).verbosity)[3]
	WIDGET_CONTROL,(*(*(*info).ctrlscp).verbose_set)[3], SET_BUTTON = ((*(*info).feedbparams).verbosity)[4]
END

;================================================================================= GENERAL WINDOW PROCEDURES
PRO CRISPEX_WINDOW, xsize, ysize, leader, title, base, wid, xoffset, yoffset, DRAWID = drawid, $
                    DRAWBASE =disp, XSCROLL = xscroll, YSCROLL = yscroll, SCROLL = scroll, $
                    RESIZING = resizing, RES_HANDLER = res_handler
; Sets up the display windows
	IF (N_ELEMENTS(RESIZING) EQ 0) THEN resizing = 0
	IF (N_ELEMENTS(LEADER) EQ 0) THEN $
    base = WIDGET_BASE(TITLE = title, TLB_FRAME_ATTR = 1, /TLB_KILL_REQUEST_EVENTS, $
                        TLB_SIZE_EVENTS = resizing) $
  ELSE $
		base = WIDGET_BASE(TITLE = STRTRIM(title), GROUP_LEADER = leader, TLB_FRAME_ATTR = 1, $
                        /TLB_KILL_REQUEST_EVENTS, TLB_SIZE_EVENTS = resizing)
	disp = WIDGET_BASE(base, /COLUMN)
  IF KEYWORD_SET(SCROLL) THEN BEGIN
    draw_verslid_base = WIDGET_BASE(disp,/ROW)
    draw_horslid_base = WIDGET_BASE(disp,/ROW)
  ENDIF ELSE draw_verslid_base = disp
	drawid = WIDGET_DRAW(draw_verslid_base, XSIZE = xsize, YSIZE = ysize, RETAIN = 2)
  IF KEYWORD_SET(SCROLL) THEN BEGIN
    yscroll = WIDGET_SLIDER(draw_verslid_base,VALUE=0,MIN=0,MAX=1,/SUPPRESS,/DRAG,$
                            EVENT_PRO='CRISPEX_SLIDER_YPOS',/VERTICAL, YSIZE=ysize)
    xscroll = WIDGET_SLIDER(draw_horslid_base,VALUE=0,MIN=0,MAX=1,/SUPPRESS,/DRAG,$
                            EVENT_PRO='CRISPEX_SLIDER_XPOS', XSIZE=xsize)
  ENDIF
	WIDGET_CONTROL, base, /REALIZE, TLB_SET_XOFFSET=xoffset, TLB_SET_YOFFSET=yoffset
	IF (N_ELEMENTS(RES_HANDLER) GT 0) THEN $
    XMANAGER, 'CRISPEX', base, EVENT_HANDLER = res_handler, /NO_BLOCK
	WIDGET_CONTROL, drawid, GET_VALUE = wid
END

PRO CRISPEX_WINDOW_OK, event, title, message1, message2, message3, message4, message5, OK_EVENT=ok_event, CANCEL_EVENT=cancel_event, CANCEL_LABEL=cancel_label, BASE=base, BLOCK=block
; Sets up the message windows with only an OK-button
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
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
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	CRISPEX_UPDATE_SX, event
	CRISPEX_UPDATE_SY, event
	IF (*(*info).overlayswitch).loopslit THEN CRISPEX_ZOOM_LOOP, event
	IF ((*(*info).meas).np GE 1) THEN CRISPEX_ZOOM_MEAS, event
	xposconstr 	= ((*(*info).dataparams).nx-1) - (*(*info).dataparams).d_nx
	yposconstr	= ((*(*info).dataparams).ny-1) - (*(*info).dataparams).d_ny
	WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = xposconstr,$
                  SET_VALUE = (*(*info).zooming).xpos 
	WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SET_SLIDER_MIN = 0, SET_SLIDER_MAX = yposconstr,$
                  SET_VALUE = (*(*info).zooming).ypos 
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).factor,xposconstr,yposconstr], $
                         labels=['Zoomfactor','Maximum xpos','Maximum ypos']
	CRISPEX_UPDATE_T, event
	IF ~KEYWORD_SET(NO_DRAW) THEN CRISPEX_DRAW, event
END

PRO CRISPEX_ZOOM_CURSORPOS, event, cursor_x, cursor_y
; Handles cursor setting as a result of zoom
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	IF (*(*info).curs).lockset THEN BEGIN
		cursor_x = (*(*info).curs).xlock
		cursor_y = (*(*info).curs).ylock
	ENDIF ELSE BEGIN
		cursor_x = (*(*info).dataparams).x
		cursor_y = (*(*info).dataparams).y
	END
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [cursor_x, cursor_y], labels=['cursor_x','cursor_y']
END	

PRO CRISPEX_ZOOM_UPDATE_SLIDERS, event, cursor_x=cursor_x, cursor_y=cursor_y, SENSITIVE=sensitive
; Handles the update of xpos and ypos sliders when changing zoomfactor
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  IF (N_ELEMENTS(SENSITIVE) EQ 1) THEN $
    sensitive = REPLICATE(sensitive,2)  $
  ELSE IF (N_ELEMENTS(SENSITIVE) LT 1) THEN $
    sensitive = [1,1]
  nx_max = (*(*info).dataparams).nx-1
  ny_max = (*(*info).dataparams).ny-1
;  IF ((*(*info).zooming).handle_extreme EQ 2) THEN BEGIN
;    IF ((*(*info).data).ratio GT 1) THEN BEGIN
;      (*(*info).dataparams).d_nx = (*(*info).winsizes).xywinx / (8. * (*(*info).zooming).factor)
;      (*(*info).dataparams).d_ny = ((*(*info).dataparams).ny / (*(*info).zooming).factor) < ny_max
;    ENDIF ELSE BEGIN
;      (*(*info).dataparams).d_nx = ((*(*info).dataparams).nx / (*(*info).zooming).factor) < nx_max
;      (*(*info).dataparams).d_ny = (*(*info).winsizes).xywiny / (8. * (*(*info).zooming).factor)
;    ENDELSE
;  ENDIF ELSE BEGIN
  	(*(*info).dataparams).d_nx = ((*(*info).dataparams).nx / (*(*info).zooming).factor) < nx_max
  	(*(*info).dataparams).d_ny = ((*(*info).dataparams).ny / (*(*info).zooming).factor) < ny_max
;  ENDELSE
	(*(*info).phiparams).d_nphi_set = (*(*info).phiparams).nphi_set / (*(*info).zooming).factor
	(*(*info).zooming).xpos = (cursor_x - (*(*info).dataparams).d_nx / 2.) > 0
	(*(*info).zooming).ypos = (cursor_y - (*(*info).dataparams).d_ny / 2.) > 0
;	(*(*info).zooming).xpos = (*(*info).zooming).xpos > 0
;	(*(*info).zooming).ypos = (*(*info).zooming).ypos > 0
	IF (((*(*info).zooming).xpos+(*(*info).dataparams).d_nx) GE (*(*info).dataparams).nx) THEN $
    (*(*info).zooming).xpos = (((*(*info).dataparams).nx-1) - (*(*info).dataparams).d_nx) > 0
	IF (((*(*info).zooming).ypos+(*(*info).dataparams).d_ny) GE (*(*info).dataparams).ny) THEN $
    (*(*info).zooming).ypos = (((*(*info).dataparams).ny-1) - (*(*info).dataparams).d_ny) > 0
	(*(*info).zooming).xpos = FIX((*(*info).zooming).xpos)
	(*(*info).zooming).ypos = FIX((*(*info).zooming).ypos)
	IF (((*(*info).feedbparams).verbosity)[3] EQ 1) THEN $
    CRISPEX_VERBOSE_GET, event, [(*(*info).zooming).xpos, (*(*info).zooming).ypos, $
      (*(*info).dataparams).d_nx, (*(*info).dataparams).d_ny], labels=['xpos','ypos','d_nx','d_ny']
	WIDGET_CONTROL, (*(*info).ctrlscp).xpos_slider, SENSITIVE = sensitive[0], $
    SET_VALUE = (*(*info).zooming).xpos
	WIDGET_CONTROL, (*(*info).ctrlscp).ypos_slider, SENSITIVE = sensitive[1], $
    SET_VALUE = (*(*info).zooming).ypos
END

PRO CRISPEX_ZOOMFAC_DECR, event
; Decreases the zoomfactor and calls CRISPEX_BGROUP_ZOOMFAC_SET
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  old_factor_idx = WHERE((*(*info).zooming).factorswitch EQ 1)
  set_factor_idx = ( old_factor_idx - 1 ) > 0
  IF (set_factor_idx NE old_factor_idx) THEN $
    set_zoomfac = CRISPEX_BGROUP_ZOOMFAC_SET(event, SET_FACTOR_IDX=set_factor_idx, $
                                      UNSET_FACTOR_IDX=old_factor_idx)
END

PRO CRISPEX_ZOOMFAC_INCR, event
; Increases the zoomfactor and calls CRISPEX_BGROUP_ZOOMFAC_SET
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
  old_factor_idx = WHERE((*(*info).zooming).factorswitch EQ 1)
  set_factor_idx = ( old_factor_idx + 1 ) < (N_ELEMENTS((*(*info).zooming).factorswitch)-1)
  IF (set_factor_idx NE old_factor_idx) THEN $
    set_zoomfac = CRISPEX_BGROUP_ZOOMFAC_SET(event, SET_FACTOR_IDX=set_factor_idx, $
                                      UNSET_FACTOR_IDX=old_factor_idx)
END

PRO CRISPEX_ZOOM_MEAS, event
; Handles the change in measurement position after a zoom event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	(*(*(*info).meas).sxp) = ((*(*(*info).meas).xp) - (*(*info).zooming).xpos ) * $
                            (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
	(*(*(*info).meas).syp) = ((*(*(*info).meas).yp) - (*(*info).zooming).ypos ) * $
                            (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
END	

PRO CRISPEX_ZOOM_LOOP, event
; Handles the change in loop positions after a zoom event
	WIDGET_CONTROL, event.TOP, GET_UVALUE = info
	IF (TOTAL(((*(*info).feedbparams).verbosity)[2:3]) GE 1) THEN CRISPEX_VERBOSE_GET_ROUTINE, event
	*(*(*info).overlayparams).sxr = (*(*(*info).loopparams).xr - (*(*info).zooming).xpos) * $
                                  (*(*info).winsizes).xywinx / ((*(*info).dataparams).d_nx+1)
	*(*(*info).overlayparams).syr = (*(*(*info).loopparams).yr - (*(*info).zooming).ypos) * $
                                  (*(*info).winsizes).xywiny / ((*(*info).dataparams).d_ny+1)
END

;===================================================================================================
;================================== MAIN PROGRAM CODE ==============================================
;===================================================================================================
PRO CRISPEX, imcube, spcube, $                ; filename of main image cube, spectral cube
              REFCUBE=refcube, $              ; filename(s) of reference image (& spectral) cube(s)
              SJICUBE=sjicube, $              ; filename of slit-jaw image cube
              MASKCUBE=maskcube, $            ; filename of mask cube
              SPECTFILE=spectfile, $          ; filename(s) of spectral save file(s)
              LINE_CENTER=line_center, $		  ; line centre and/or wavelength information
	            DT=dt, $                        ; time step in seconds
              EXTS=exts, $                    ; exact timeslices keyword
              MNSPEC=mnspec, $                ; mean spectrum over selected scans
              SINGLE_CUBE=single_cube, $      ; single full cube call
              SCALE_STOKES=scale_stokes, $    ; scale Stokes spectra internally
              VALS_IMG=vals_img, $            ; get main cube values under cursor
              VALS_REF=vals_ref, $            ; get reference cube values under cursor
              NO_WARP=no_warp, $              ; don't warp nonequidistant spectral slices
              SCALE_CUBES=scale_cubes, $      ; scale cubes
              XTITLE=xtitle, YTITLE=ytitle,$; custom detailed spectrum xtitle and ytitle
              WINDOW_LARGE=window_large, $    ; draw large windows for small cubes
              VERBOSE=verbose                 ; program verbosity

;========================= PROGRAM-INFO ON CALL W/O PARAMS
	IF N_PARAMS() LT 1 THEN BEGIN
		MESSAGE,'Syntax: CRISPEX, imcube, spcube, REFCUBE=refcube, MASKCUBE=maskcube, '+$
            'SPECTFILE=spectfile, LINE_CENTER=line_center, DT=dt, EXTS=exts, MNSPEC=mnspec, '+$
            'SINGLE_CUBE=single_cube, SCALE_STOKES=scale_stokes, VALS_IMG=vals_img, '+$
            'VALS_REF=vals_ref, NO_WARP=no_warp, SCALE_CUBES=scale_cubes, XTITLE=xtitle, '+$
            'YTITLE=ytitle, WINDOW_LARGE=window_large, VERBOSE=verbose', /INFO
		RETURN
	ENDIF

;========================= PROGRAM VERBOSITY CHECK
	IF (N_ELEMENTS(VERBOSE) NE 1) THEN BEGIN			
		IF (N_ELEMENTS(VERBOSE) GT 1) THEN $
      MESSAGE,'ERROR: The VERBOSE keyword may only be set to a single integer number. Reverting '+$
              'to default verbosity level 0.'
		verbose = 0
		verbosity = [0,0,0,0,0]
	ENDIF ELSE BEGIN
		verbose >= 0	&	verbose <= 26
	ENDELSE
	verbosity = CRISPEX_DEC2BIN(verbose)   ; Convert verbosity value to binary array

;========================= CRISPEX DIRECTORY TREE CHECK
	file_crispex    = (ROUTINE_INFO('CRISPEX',/SOURCE)).PATH      ; Path to compiled CRISPEX
	dir_crispex     = FILE_DIRNAME(file_crispex,/MARK_DIRECTORY)
	dir_aux         = dir_crispex+'aux'+PATH_SEP()                ; Path to auxiliary routines
	dir_resources   = dir_crispex+'resources'+PATH_SEP()          ; Path to resources container
	dir_buttons     = dir_resources+'buttons'+PATH_SEP()          ; Path to button images
  dir_settings    = CRISPEX_CONFIG_DIR()+PATH_SEP()             ; Path to settings container
	dir_settings_write  = FILE_TEST(dir_settings, /WRITE)                 ; Check for cpft dir writeability
	IF (verbosity[1] EQ 1) THEN $
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK,'CRISPEX has been compiled from: '+file_crispex

;========================= VERSION AND REVISION NUMBER
	version_number = '1.6.3'
  ; Get revision number from CVS $Id
  SPAWN,"grep '$Id' "+file_crispex, id_string
  cvs_idn = (STRSPLIT(id_string[0],' ',/EXTRACT))[3]
  cvs_rev = (STRSPLIT(cvs_idn,'.',/EXTRACT))[1]
  revision_number = STRTRIM(634L+LONG(cvs_rev)-64L,2)   ; rev_nr=634, cvs_rev=64 when implemented
	IF (verbosity[1] EQ 1) THEN $
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK,'Version and revision number: '+version_number+$
      ' ('+revision_number+')'

;========================= LOAD PREFERENCES
  ; Define default preferences
	default_startupwin = 1      &  default_interpspslice = 1  ; Show startup win, interpolate slices
	default_autoplay = 0        &  default_defsaveid = 0      ; 0 = yyyymmdd, 1 = ddmmyyyy
	default_defipath = 0        &  default_defopath = 0			  ; 0 = local working directory, 1 = saved directory
	default_bgplotcol = 255     &  default_plotcol = 0
	default_phislice_update = 0 &  default_slices_imscale = 0
  default_histo_opt_val = 0.0001  & default_gamma_val = 1.0
	cpreffiles = FILE_SEARCH(dir_settings+'crispex.cpref', COUNT = cpreffilecount)
	IF (cpreffilecount GE 1) THEN BEGIN     ; If preference file is present, load preference file
		RESTORE, cpreffiles[0] 
		IF (verbosity[1] EQ 1) THEN $
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Preferences restored from: '+dir_settings+'crispex.cpref'
		resave_preferences = ((N_ELEMENTS(phislice_update) NE 1) OR (N_ELEMENTS(slices_imscale) NE 1))
    ; Failsafe inheritances from older CRISPEX versions
    ; Automatic phislice update
		IF (N_ELEMENTS(phislice_update) NE 1) THEN phislice_update = default_phislice_update
    ; Scale slices with image scaling
		IF (N_ELEMENTS(slices_imscale) NE 1) THEN slices_imscale = default_slices_imscale
    ; HISTO_OPT value
		IF (N_ELEMENTS(histo_opt_val) NE 1) THEN histo_opt_val = default_histo_opt_val
    ; Gamma value
		IF (N_ELEMENTS(gamma_val) NE 1) THEN gamma_val = default_gamma_val
	ENDIF ELSE BEGIN                        ; If no preference file is present, set defaults
		startupwin = default_startupwin           &  interpspslice = default_interpspslice
		autoplay = default_autoplay               &  defsaveid = default_defsaveid
		defipath = default_defipath               &  defopath = default_defopath
		bgplotcol = default_bgplotcol             &  plotcol = default_plotcol
		phislice_update = default_phislice_update &  slices_imscale = default_slices_imscale
    histo_opt_val = default_histo_opt_val     &  gamma_val = default_gamma_val
		resave_preferences = 0
	ENDELSE

;------------------------- SETTINGS FOR PERFORMANCE SAVE FILE
	hostname = GETENV('HOSTNAME')
  IF (STRLEN(STRCOMPRESS(hostname)) NE 0) THEN hostname += '.'
	cpftfile = FILE_SEARCH(dir_settings+'crispex.'+hostname+'cpft', COUNT = cpftfilecount)
	IF cpftfilecount THEN BEGIN   ; If cpft file is present, restore
		RESTORE, cpftfile[0] 
		IF (verbosity[1] EQ 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Restored '+cpftfile[0]+'.'
	ENDIF ELSE BEGIN              ; If not, then initialise variables
		IF (verbosity[1] EQ 1) THEN BEGIN
			CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, ['No CRISPEX performance test file (crispex.'+hostname+$
                                             'cpft) found to restore in ',dir_settings]
		ENDIF
		estimate_lx = 0             ; Size variable for time estimate
		estimate_run = 0            ; Run counter for time estimate
		estimate_time = 0.          ; Updated time estimate per unit run and unit slice size
	ENDELSE 

;------------------------- SETTINGS FOR INSTANCES SAVE FILE
	instfilename = 'crispex.'+hostname+'inst'
	IF dir_settings_write THEN BEGIN    ; If instances directory is writeable, start procedures
		instfile = FILE_SEARCH(dir_settings+instfilename, COUNT = instfilecount)
		IF instfilecount THEN BEGIN   ; If inst file is present for current hostname, add current
			IF (verbosity[1] EQ 1) THEN $
        CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Opening existing instance tracking file: '+$
                                               instfilename+'.'
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
			OPENU, unit2, dir_settings+instfilename, WIDTH = 360, /GET_LUN, /APPEND
		ENDIF ELSE BEGIN              ; If no inst file present for current hostname, make one
			IF (verbosity[1] EQ 1) THEN BEGIN
				CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'No CRISPEX instance tracking file ('+instfilename+$
                                               ') found in '+dir_settings+'. Creating file.'
			ENDIF
			where_crispex = -1
			OPENW, unit2, dir_settings+instfilename, WIDTH = 360, /GET_LUN
			PRINTF, unit2, '# routine_name	version		revision	ID'
		ENDELSE
		IF (where_crispex[0] NE -1) THEN $
      set_instance_id = STRTRIM((instance_id[where_crispex])[WHERE(instance_id[where_crispex] EQ $
                        MAX(instance_id[where_crispex], /NAN))] + 1,2) $
    ELSE $
      set_instance_id = STRTRIM(0,2)
		PRINTF, unit2, 'CRISPEX	'+version_number+'	'+revision_number+'	'+set_instance_id
		FREE_LUN, unit2
		IF (set_instance_id GE 1) THEN instance_label = '-'+set_instance_id ELSE instance_label = ''
		IF (verbosity[1] EQ 1) THEN $
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Written instance ID ('+set_instance_id+') to '+$
                                             instfilename+'.'
	ENDIF ELSE BEGIN
		set_instance_id = ''
		instance_label = ''
		instfilecount = 0
		IF (verbosity[1] EQ 1) THEN BEGIN
 			CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'ERROR: Could not write CRISPEX instance tracking '+$
                                             'file '+instfilename+'to '+dir_settings+'. Permission denied.'
		ENDIF
	ENDELSE
	
;------------------------- INPUT/OUTPUT PATH SETTINGS
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
  IF verbosity[1] THEN BEGIN
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Input path set to: '+STRTRIM(ipath,2), /NEWLINE, $
                                           /NO_ROUTINE
    CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Output path set to: '+STRTRIM(opath,2), /NO_ROUTINE
  ENDIF

;========================= START-UP WINDOW AND SCREEN PARAMETERS
	screeninfo  = OBJ_NEW('IDLsysMonitorInfo')
	nmonitors    = screeninfo -> GetNumberOfMonitors()
	screensizes = screeninfo -> GetRectangles()
	IF (nmonitors GT 1) THEN BEGIN   ; Define monitor order if multiple monitors attached
		monitor_order       = (INDGEN(nmonitors))[SORT(screensizes[0,0:1])]
		monitor_xsize_order = (INDGEN(nmonitors))[SORT(screensizes[2,0:1])]
		monitor_ysize_order = (INDGEN(nmonitors))[SORT(screensizes[3,0:1])]
	ENDIF ELSE BEGIN                ; Default monitor order if only one monitor attached
		monitor_order       = 0
		monitor_ysize_order = 0
	ENDELSE
	x_screen_mid  = screensizes[2,monitor_order[0]]/2.    ; x-coordinate of central screen pixel
	y_screen_mid  = screensizes[3,monitor_order[0]]/2.    ; y-coordinate of central screen pixel
	startup_im    = REBIN(REFORM(TOTAL((CRISPEX_READ_BMP_BUTTONS('crispex_startup.bmp',$
                  dir_resources))[*,*,1:2],3)),400,300)
	startup_nx    = (SIZE(startup_im))[1]                 ; x-size of statup window image
	startup_ny    = (SIZE(startup_im))[2]                 ; y-size of statup window image
	startup_xpos  = FIX(x_screen_mid-startup_nx/2.)       ; x-position of startup window
	startup_ypos  = FIX(y_screen_mid-startup_ny/2.)       ; y-position of startup window
	xout          = REPLICATE(24,9)
	yout          = REPLICATE(FIX(startup_ny/2.5)+10,9)-INDGEN(9)*15
	IF startupwin THEN BEGIN  ; If startup window is to be shown, launch window
		CRISPEX_WINDOW, startup_nx, startup_ny, 0, 'CRISPEX', startuptlb, startupwid, startup_xpos, $
                    startup_ypos, DRAWID = startupdrawid, DRAWBASE = drawbase
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Initializing... '
	ENDIF

;========================= READ-IN AND INITIALISATION OF FILES
  ; N.B.: After CRISPEX v1.6.3 FITS cubes have become the standard. Old read-in procedures 
  ; are retained in compatability for older cubes. Differentiation is performed based on filename
  ; extension, where FITS cubes are assumed to have a *.fits extension (case insensitive).
  IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Reading input files... '
  IF ((BYTE(1L,0,1))[0] EQ 1) THEN endian = 'l' ELSE endian = 'b' ; Check endianness of machine
  IF (N_ELEMENTS(DT) NE 1) THEN dt = 0.
  
  ; Handle input file headers by parsing them into the hdr structure; first initialise hdr
  hdr = {$
            imtype:0, sptype:0, refimtype:0, refsptype:0, sjitype:0, masktype:0, $
            imoffset:0L, spoffset:0L, refimoffset:0L, refspoffset:0L, sjioffset:0L, maskoffset:0L, $
            imendian:'b', spendian:'b', refimendian:'b', refspendian:'b', sjiendian:'b', $
            maskendian:'b',endian:endian, $
            imcube_compatibility:0, spcube_compatibility:0, refimcube_compatibility:0, $
            refspcube_compatibility:0, maskcube_compatibility:0, multichannel:0, $
            nx:0L, ny:0L, nlp:1L, mainnt:1L, ns:1L, imnt:0L, spnt:0L, refspnx:0L, refspny:0L, $
            refnx:0L, refny:0L, refnlp:0L, refnt:0L, refns:0L, refimnt:0L, refspnt:0L, $
            sjinx:0L, sjiny:0L, sjint:1L, sjidx:0., sjidy:0., sjix0:0L, sjiy0:0L, $
            masknx:0L, maskny:0L, masknt:0L, dx:0., dy:0., dt:dt, dx_fixed:0, $
            imns:0L, spns:0L, imstokes:'', spstokes:'', imdiagnostics:'', $
            xunit:'arcsec', yunit:'arcsec', tunit:'', lpunit:'', sunit:'', bunit:'counts', $
            xlabel:'x', ylabel:'y', tlabel:'Frame number', lplabel:'Spectral position', $
            slabel:'', blabel:'Intensity', $
            stokes_enabled:[0,0,0,0], scalestokes_max:0, ndiagnostics:1, nrefdiagnostics:1, $
            refxunit:'arcsec', refyunit:'arcsec', reflpunit:'', refbunit:'counts', $
            refxlabel:'x', refylabel:'y', reftlabel:'', reflplabel:'', refblabel:'Intensity', $
            sjibunit:'counts', $
            xtitle:STRARR(2), ytitle:STRARR(2), refspxtitle:'', spxtitle:'', spytitle:'', $
            ipath:ipath, opath:opath, instance_label:instance_label, $
            lunim:0, lunsp:0, lunrefim:0, lunrefsp:0, lunsji:0, lunmask:0, $
            imfilename:'', spfilename:'', refimfilename:'', refspfilename:'', sjifilename:'', $
            maskfilename:'', $
            spfile:0, onecube:0, single_cube:[0,0], showref:0, refspfile:0, sjifile:0, maskfile:0, $
            imdata:PTR_NEW(0), scan:PTR_NEW(0), spdata:PTR_NEW(0), spectra:PTR_NEW(0), $
            refdata:PTR_NEW(0), refslice:PTR_NEW(0), refspdata:PTR_NEW(0), refscan:PTR_NEW(0), $
            refsspscan:PTR_NEW(0), sjidata:PTR_NEW(0), sjislice:PTR_NEW(0), maskdata:PTR_NEW(0), $
            verbosity:verbosity $
            }

  CRISPEX_IO_OPEN_MAINCUBE, IMCUBE=imcube, SPCUBE=spcube, HDR_IN=hdr, HDR_OUT=hdr, $
                            SINGLE_CUBE=single_cube, STARTUPTLB=startuptlb, $
                            IO_FAILSAFE_MAIN_ERROR=io_failsafe_main_error
  IF (io_failsafe_main_error EQ 1) THEN RETURN

  CRISPEX_IO_OPEN_REFCUBE, REFCUBE=refcube, HDR_IN=hdr, HDR_OUT=hdr, $
                            SINGLE_CUBE=single_cube, $
                            CUBE_COMPATIBILITY=refcube_compatibility, $
                            IO_FAILSAFE_REF_ERROR=io_failsafe_ref_error, $
                            IO_FAILSAFE_MAIN_REF_ERROR=io_failsafe_main_ref_error
  IF ((io_failsafe_ref_error EQ 1) OR (io_failsafe_main_ref_error EQ 1)) THEN RETURN
  IF (N_ELEMENTS(REFCUBE) LT 1) THEN $
    hdr = CREATE_STRUCT(hdr, 'refdiagnostics', 'N/A', 'refdiag_start', 0, 'refdiag_width', 1, $
            'tarr_ref', 0, 'tarr_raster_ref', 0, 'toffset_ref', 0, 'hdrs_ref', PTR_NEW(''))

  CRISPEX_IO_OPEN_SJICUBE, SJICUBE=sjicube, HDR_IN=hdr, HDR_OUT=hdr, $  
                            STARTUPTLB=startuptlb, $
                            IO_FAILSAFE_SJI_ERROR=io_failsafe_sji_error
  IF (io_failsafe_sji_error EQ 1) THEN RETURN
  IF (N_ELEMENTS(sjicube) NE 1) THEN $
    hdr = CREATE_STRUCT(hdr, 'tarr_sji', 0, 'tsel_sji', 0, 'rastercont', 0, 'hdrs_sji', PTR_NEW(''))
  
  CRISPEX_IO_OPEN_MASKCUBE, MASKCUBE=maskcube, HDR_IN=hdr, HDR_OUT=hdr, STARTUPTLB=startuptlb, $
                            IO_FAILSAFE_MASK_ERROR=io_failsafe_mask_error
  IF (io_failsafe_mask_error EQ 1) THEN RETURN

	IF (hdr.refnlp NE hdr.nlp) THEN BEGIN
		eqnlps = 0 
		refslid_sens = (hdr.showref AND (hdr.refnlp GT 1)) 
	ENDIF ELSE BEGIN
		eqnlps = 1
		refslid_sens = 0
	ENDELSE
	showrefls = (hdr.refspfile OR (hdr.refnlp GT 1))
	
  IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, $
                                                      'Reading input files... done!'

	DEVICE, DECOMPOSE = 0	

;========================= SETTING START-UP OPTIONS 
;------------------------- PARAMETERS FROM MEAN SPEC
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(parameters from/for mean spectrum)', /OPT, /OVER
	feedback_text = ['Setting start-up options... ','> Parameters from/for mean spectrum... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	detspect_scale_enable = (hdr.nlp GT 1)
	detspect_scale = (hdr.nlp GT 1)
	ref_detspect_scale = (hdr.refnlp GT 1)

	scalestokes = (N_ELEMENTS(SCALE_STOKES) EQ 1) 
	IF (N_ELEMENTS(SCALE_STOKES) GT 1) THEN $
    PRINT, 'WARNING: The SCALE_STOKES keyword must either be set or supplied with a single '+$
           'number! Reverting to standard scaling.'

  ; Set final XTITLE and YTITLE variables for LS and SP windows
	IF (N_ELEMENTS(XTITLE) EQ 1) THEN BEGIN         ; Check whether XTITLE has 1 element
    IF (N_ELEMENTS(REFXTITLE) EQ 1) THEN $        ; Check whether REFXTITLE has 1 element
      hdr.xtitle = [xtitle[0],refxtitle[0]] $     ; If so, fill xtitle variable with both
    ELSE $                                        ; Else, fill xtitle variable only with main
      hdr.xtitle = [xtitle[0],''] 
  ENDIF ELSE IF (N_ELEMENTS(XTITLE) EQ 2) THEN $  ; Check whether XTITLE has 2 elements
    hdr.xtitle = xtitle $                         ; If so, use that
  ELSE IF (N_ELEMENTS(REFXTITLE) EQ 1) THEN $     ; If not, check whether REFXTITLE has 1 element
    hdr.xtitle = ['',refxtitle[0]]                ; If so, set 2nd xtitle element to that
	IF (N_ELEMENTS(YTITLE) EQ 1) THEN BEGIN         ; Check whether YTITLE has 1 element
    IF (N_ELEMENTS(REFYTITLE) EQ 1) THEN $        ; Check whether REFYTITLE has 1 element
      hdr.ytitle = [ytitle[0],refytitle[0]] $     ; If so, fill ytitle variable with both
    ELSE $                                        ; Else, fill ytitle variable only with main
      hdr.ytitle = [ytitle[0],''] 
  ENDIF ELSE IF (N_ELEMENTS(YTITLE) EQ 2) THEN $  ; Check whether YTITLE has 2 elements
    hdr.ytitle = ytitle $                         ; If so, use that
  ELSE IF (N_ELEMENTS(REFYTITLE) EQ 1) THEN $     ; If not, check whether REFYTITLE has 1 element
    hdr.ytitle = ['',refytitle[0]]                ; If so, set 2nd ytitle element to that
  
  CRISPEX_IO_SETTINGS_SPECTRAL, HDR_IN=hdr, HDR_OUT=hdr, MNSPEC=mnspec, SPECTFILE=spectfile,$
                                LINE_CENTER=line_center, NO_WARP=no_warp, STARTUPTLB=startuptlb, $
                                IO_FAILSAFE_MNSPEC_ERROR=io_failsafe_mnspec_error, $
                                IO_FAILSAFE_IMSPECTFILE_ERROR=io_failsafe_imspectfile_error, $
                                IO_FAILSAFE_REFSPECTFILE_ERROR=io_failsafe_refspectfile_error, $
                                IO_FAILSAFE_LINE_CENTER_ERROR=io_failsafe_line_center_error

  IF ((io_failsafe_mnspec_error EQ 1) OR (io_failsafe_imspectfile_error EQ 1) OR $
      (io_failsafe_refspectfile_error EQ 1) OR (io_failsafe_line_center_error EQ 1)) THEN RETURN
	
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(parameters from/for mean spectrum)', /OPT, /OVER, /DONE, $
                                        REPEAT_STAGE=(verbosity[1] OR $
                                          (io_failsafe_imspectfile_error EQ 2) OR $
                                          (io_failsafe_refspectfile_error EQ 2))
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],$
                   '> Parameters from/for mean spectrum... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SLIT PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial slit parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Initial slit parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

  nphi	= LONG(CEIL(SQRT( FLOAT(hdr.nx)^2 + FLOAT(hdr.ny)^2 )))	; Determine maximum number of slitpositions
	angle = 90														; Set initial angle of the slit

  IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial slit parameters)', /OPT, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial slit parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL PLAYBACK PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial playback parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Initial playback parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	t_first		= 0L													; Set number of first frame		
	IF (hdr.mainnt EQ 1) THEN BEGIN
		t_last = 2L
		t_last_tmp = 0L
		t_slid_sens = 0
	ENDIF ELSE BEGIN
		t_last = hdr.mainnt-1													; Set number of last frame
		t_last_tmp = t_last
		t_slid_sens = 1
	ENDELSE
	t_step		= 1L													; Set initial timestep
	t_speed 	= 10													; Set initial animation speed
	direction 	= 1													; Set initial animation direction
	;nt		= FLOAT(hdr.nt)												; Convert the number of timesteps to float
	t_start = t_first

	IF (hdr.dt EQ 1) THEN BEGIN
		IF (hdr.spfile OR hdr.onecube) THEN BEGIN
			dt_set = 1
			IF (N_ELEMENTS(SPYTITLE) NE 1) THEN spytitle = hdr.spytitle ;'Time (s)'
		ENDIF ELSE BEGIN
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Calling CRISPEX with DT has no influence when '+$
        'no SPCUBE is supplied. Setting seconds per timestep to default value.', /WARNING, $
        /NO_ROUTINE, /NEWLINE
			dt_set = 0
			hdr.dt = 1.
			spytitle = 'Frame number'
		ENDELSE
	ENDIF ELSE BEGIN
		dt_set = 0
		hdr.dt = 1.
		spytitle = 'Frame number'
	ENDELSE

	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial playback parameters)', /OPT, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial playback parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SPECTRAL PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial spectral parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Initial spectral parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	lp_first 	= 0L													; Set number of first lineposition
	lp_last		= hdr.nlp-1L													; Set number of last lineposition
	sp		= hdr.mainnt * hdr.nlp												; Set spectral dimension
	lp_start 	= hdr.lc
	lp_ref_first = lp_first
	IF showrefls THEN BEGIN
		lp_ref_last = hdr.refnlp - 1L
		lp_ref_start = hdr.reflc
	ENDIF ELSE IF (hdr.refnlp GT 1) THEN BEGIN
		lp_ref_last = lp_last
		lp_ref_start = hdr.lc
	ENDIF ELSE BEGIN
		lp_ref_last = 1L
		lp_ref_start = 0L
	ENDELSE
  lp_ref_lock = eqnlps ;(eqnlps AND (lp_ref_start EQ lp_start))
  IF (lp_ref_lock) THEN lp_ref_start = lp_start
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial spectral parameters)', /OPT, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial spectral parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- WINDOW SIZES (CHANGE ONLY
;--------------------------------------------------------------------------------- NUMERICAL VALUES!)
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(window sizes)', /OPT, /OVER
	feedback_text = [feedback_text,'> Window sizes... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	
	heightset 	= 0
	refheightset 	= 0
	IF ((N_ELEMENTS(hdr.xtitle) GE 1) AND (N_ELEMENTS(hdr.xtitle) LE 2)) THEN BEGIN
		IF (STRCOMPRESS(hdr.xtitle[0]) NE '') THEN BEGIN
			heightset = (STRCMP(hdr.xtitle[0],'Height',6,/FOLD_CASE) OR $
                   STRCMP(hdr.xtitle[0],'z',1,/FOLD_CASE))
      IF heightset THEN BEGIN
  			hdr.v_dop_set[0] = 0
  			hdr.spxtitle = hdr.xtitle[0]
      ENDIF
		ENDIF
		IF (N_ELEMENTS(hdr.xtitle) EQ 2) THEN BEGIN
			IF (STRCOMPRESS(hdr.xtitle[1]) NE '') THEN BEGIN
				refheightset = (STRCMP(hdr.xtitle[1],'Height',6,/FOLD_CASE) OR $
                        STRCMP(hdr.xtitle[1],'z',1,/FOLD_CASE))
        IF refheightset THEN BEGIN
				  hdr.v_dop_set[1] = 0
				  hdr.refspxtitle = hdr.xtitle[1]
        ENDIF
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

  ; Load (ref)lsytitles from FITS cubes or SPECTFILE
	lsytitle	= hdr.ytitle[0]
	reflsytitle	= hdr.ytitle[1]
  ; Override by XTITLE and YTITLE keywords, if set
	IF ((N_ELEMENTS(YTITLE) GE 1) AND (N_ELEMENTS(YTITLE) LE 2)) THEN BEGIN
		IF (STRCOMPRESS(ytitle[0]) NE '') THEN lsytitle = ytitle[0]
		IF (N_ELEMENTS(YTITLE) EQ 2) THEN BEGIN
			IF (STRCOMPRESS(ytitle[1]) NE '') THEN reflsytitle = ytitle[1]
		ENDIF
	ENDIF

	xdelta		= 20													; Extra xoffset for positioning of windows
	ydelta		= 40													; Extra yoffset for positioning of windows
  minsize   = 200.
 
  ; Determine image and pixel aspect ratio
  ratio      = FLOAT(ABS(hdr.nx)) / FLOAT(ABS(hdr.ny))                  
  pixelratio = FLOAT(ABS(hdr.dx)) / FLOAT(ABS(hdr.dy))                  

  ; Determine default window sizes (i.e., regardless of monitor size)
  imwinx_default = hdr.nx
  imwiny_default = hdr.ny

  ; Handle pixel aspect ratio for default sizes
  IF (pixelratio GT 1) THEN imwinx_default *= pixelratio ELSE $
    IF (pixelratio LT 1) THEN imwiny_default *= pixelratio
  
  ; Check for extreme aspect ratio and small dimensions
  extreme_aspect = (((ratio GT 5.) AND (imwiny_default LT minsize)) OR $
                    ((ratio LT 0.2) AND (imwinx_default LT minsize))) 

  ; Get main screen sizes
  x_scr_size = screensizes[2,monitor_order[0]]
	y_scr_size = screensizes[3,monitor_order[0]]

  ; Check whether window would in principle fit (only if xsize > nx+space for spectral 
  ; windows AND ysize > ny):
	IF ((x_scr_size GT (imwinx_default+2*xdelta+0.45*x_scr_size)) AND $
      (y_scr_size GT imwiny_default)) THEN BEGIN		
    ; If xsize is small, then still go to old settings procedures
		IF ((imwinx_default LT 0.48 * x_scr_size) AND $
        (imwiny_default LT (0.48 * x_scr_size / ratio)) AND $
        (extreme_aspect OR KEYWORD_SET(WINDOW_LARGE))) THEN BEGIN				
      ; Set maximum x- and y-extent of image window
      IF (extreme_aspect AND (hdr.nx EQ 1)) THEN $
        imwinx = 5*hdr.nx $
      ELSE $
			  imwinx 	= 0.48 * x_scr_size											
			imwiny 	= imwinx / ratio
      ; Failsafe to avoid a window larger than the screensize
		  IF (imwiny GT y_scr_size) THEN BEGIN										
			  imwiny = 0.9 * y_scr_size											
        IF (extreme_aspect AND (hdr.nx EQ 1)) THEN $
          imwinx = 5*hdr.nx $
        ELSE $
			    imwinx = imwiny * ratio
		  ENDIF
			IF (verbosity[1] EQ 1) THEN $
        msg = 'User screen resolution allows 1:1 image window sizing, '+$
          'but dimensions are small. '
    ; Else fit the window with actual data dimensions
		ENDIF ELSE BEGIN
      ; Use actual nx/ny as imwinx/imwiny
			imwinx	= hdr.nx                  
      IF (extreme_aspect AND (hdr.nx EQ 1)) THEN imwinx *= 5
			imwiny	= hdr.ny                 
			IF (verbosity[1] EQ 1) THEN $
        msg = 'User screen resolution allows 1:1 image window sizing. '
		ENDELSE
  ; Else use the old procedures to determine imwinx and imwiny
	ENDIF ELSE BEGIN													
    ; Set maximum x- and y-extent of image window
		imwinx 	= 0.48 * x_scr_size
		imwiny 	= imwinx / ratio	 
    ; Failsafe to avoid a window larger than the screensize
		IF (imwiny GT y_scr_size) THEN BEGIN										
			imwiny = 0.9 * y_scr_size
			imwinx = imwiny * ratio
		ENDIF
		IF (verbosity[1] EQ 1) THEN $
      msg = 'User screen resolution does not allow for 1:1 image window sizing. '
	ENDELSE
  ; Fix window sizes for pixelratio
  IF (pixelratio GT 1) THEN $
    imwinx *= pixelratio $
  ELSE IF (pixelratio LT 1) THEN $
    imwiny *= pixelratio
  IF (verbosity[1] EQ 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, msg+'Image window set to '+$
                                STRTRIM(imwinx,2)+'x'+STRTRIM(imwiny,2)+'.', /NEWLINE, /NO_ROUTINE
                                
  ; Set zoomfactors
  zoomfactors = [1.,2.,3.,4.,6.,8.]
  factorswitch = [1B,BYTARR(N_ELEMENTS(zoomfactors)-1)]

  ; If SJI cube supplied, determine sizes
  IF hdr.sjifile THEN BEGIN
    ; Determine image and pixel aspect ratio
    sjiratio      = FLOAT(ABS(hdr.sjinx)) / FLOAT(ABS(hdr.sjiny))
    sjipixelratio = FLOAT(ABS(hdr.sjidx)) / FLOAT(ABS(hdr.sjidy))
    sjiwinx_default = hdr.sjinx
    sjiwiny_default = hdr.sjiny
    ; Handle pixel aspect ratio
;    IF (sjipixelratio GT 1) THEN sjiwinx_default *= sjipixelratio ELSE $
;      IF (sjipixelratio LT 1) THEN sjiwiny_default *= sjipixelratio
;    print,hdr.sjidx,hdr.sjidy,hdr.dx,hdr.dy
      sjiwinx = sjiwinx_default/imwinx_default*imwinx
      sjiwiny = sjiwiny_default/imwiny_default*imwiny
      IF (sjiwinx GT 0.75*x_scr_size) THEN BEGIN
        sjiwinx = 0.75*x_scr_size
        sjiwiny = sjiwinx / sjiratio
        imwinx = imwinx_default/sjiwinx_default*sjiwinx
        imwiny = imwiny_default/sjiwiny_default*sjiwiny
      ENDIF
  ENDIF ELSE BEGIN
    sjiwinx = 0
    sjiwiny = 0
  ENDELSE

	windowx		= 0.2 * x_scr_size											; Set maximum x-extent of spectral win
	IF (hdr.mainnt GE 50) THEN windowy = imwiny ELSE windowy = imwiny/2. > (y_scr_size/2.)
	lswinx 		= 0.2 * x_scr_size											; Set maximum x-extent of loc spec win

  ;; If reference cube present, check if it would fit next to main image
  ;refxoffset = 0
  ;refyoffset = 0
  ;IF hdr.showref THEN BEGIN
  ;  windows_xextent = imwinx*2 + windowx + lswinx + xdelta*3
  ;  IF (windows_xextent LE x_scr_size) THEN refxoffset = imwinx + 2*xdelta $
  ;    ELSE refyoffset = ydelta
  ;ENDIF ELSE windows_xextent = imwinx + windowx + lswinx + xdelta*2

	;spxoffset = imwinx + refxoffset + (1+(windows_xextent GT x_scr_size))*xdelta
	;lsxoffset = spxoffset + windowx + xdelta

	xswinx		= windowx												; Set maximum x-extent of x-slice window
	xswiny		= windowy												; Set maximum y-extent of x-slice window

	lswintitle	= ['Detailed spectrum','Height distribution']
	lsmargin 	= 0.1
	lswall 		= 0.03
	ticklen 	= 0.01
	xsize 		= 1.*lswinx

	IF (hdr.ns LE 2) THEN BEGIN
		npanels = hdr.ns	&	cols = hdr.ns	&	rowarr = REPLICATE(0,hdr.ns)
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
	IF (hdr.v_dop_set[0] EQ 1) THEN $
    lswiny = (lsmargin + rows*lsheight + rows*lsmargin + (rows-1)*lswall) * lswinx $
  ELSE $
    lswiny = (lswall + rows*lsheight + rows*lsmargin) * lswinx
	lsx0 		= lsmargin * lswinx/lswinx + (INDGEN(npanels) MOD cols) * (lswidth + lsmargin) * lswinx/lswinx
	lsx1 		= lsx0 + lswidth * lswinx/lswinx
	lsy0 		= lsmargin * lswinx/lswiny + rowarr * (lsheight + lsmargin + hdr.v_dop_set[0]*lswall) * lswinx/lswiny
	lsy1 		= lsy0 + lsheight * lswinx/lswiny
	lsxmargin_init	= lsmargin * lswinx
	lsxwall_init	= lswall * lswinx
	lsxticklen 	= ticklen / lsheight
	lsyticklen 	= ticklen / lswidth

	reflswintitle	= ['Reference detailed spectrum','Reference height distribution']
	reflswidth 	= (xsize/lswinx - (lsmargin + lswall))
	reflsheight 	= reflswidth * 2D / (1 + SQRT(5))
	reflswinx	= lswinx
	IF (hdr.v_dop_set[1] EQ 1) THEN $
    reflswiny = (lsmargin + reflsheight + lsmargin) * lswinx $
  ELSE $
    reflswiny = (lsmargin + reflsheight + lswall) * lswinx
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
	IF hdr.spfile THEN spwiny = windowy ELSE spwiny = imwiny
  spwiny = imwiny - lswiny
	IF ((hdr.v_dop_set[0] EQ 1) OR (hdr.ns GT 1)) THEN BEGIN
		spheight = (1. - (spmargin * 2.) * spwinx/spwiny)
		phisheight = (1. - (spmargin * 2.) * phiswinx/phiswiny)
	ENDIF ELSE BEGIN
		spheight = (1. - (spmargin + spwall) * spwinx/spwiny)
		phisheight = (1. - (spmargin + spwall) * phiswinx/phiswiny)
	ENDELSE
	IF (hdr.v_dop_set[1] EQ 1) THEN $
    refspheight = (1. - (spmargin * 2.) * spwinx/spwiny) $
  ELSE $
    refspheight = (1. - (spmargin + spwall) * spwinx/spwiny)
	refspwiny	= imwiny-lswiny ;windowy

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
	
	IF ((hdr.spfile EQ 1) OR (hdr.single_cube[0] GE 1)) THEN ntreb = yplspw * spwiny ELSE ntreb = 0						; actual nt rebinning factor
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

	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(window sizes)', /OPT, /OVER, /DONE, $
                                        REPEAT_STAGE=verbosity[1]
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Window sizes... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIAL SPATIAL PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial spatial parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Initial spatial parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	x_first		= 0L													; Set number of first x-coordinate
	x_last		= hdr.nx-1													; Set number of last x-coordinate
	y_first		= 0L												; Set number of first y-coordinate
	y_last		= hdr.ny-1													; Set number of last y-coordinate
	x_start		= FLOAT(FLOOR(hdr.nx/2))												; Determine the middle x-coordinate
	sx_start	= x_start * imwinx / FLOAT(hdr.nx)										; Convert that to device
	y_start		= FLOAT(FLOOR(hdr.ny/2))												; Determine the middle y-coordinate
	sy_start	= y_start * imwiny / FLOAT(hdr.ny)										; Convert that to device
  arcsecpix	= hdr.dx
  IF hdr.sjifile THEN BEGIN
    xsji_start = $
      hdr.xyrastersji[x_start,0] + (x_start * hdr.dx / hdr.sjidx - $
      (hdr.xyrastersji[x_start,0] - hdr.xyrastersji[0,0]))
    ysji_start = $
      hdr.xyrastersji[x_start,1] + (y_start * hdr.dy / hdr.sjidy - $
      (hdr.xyrastersji[x_start,1] - hdr.xyrastersji[0,1])) 
    sxsji_start = xsji_start * sjiwinx / hdr.sjinx
    sysji_start = ysji_start * sjiwiny / hdr.sjiny
  ENDIF ELSE BEGIN
    xsji_start = 0L
    ysji_start = 0L
    sxsji_start = 0.
    sysji_start = 0.
  ENDELSE

	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial spatial parameters)', /OPT, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial spatial parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- SCALING PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial scaling parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Initial scaling parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	imagescale = PTR_NEW([0,0,0,0])												; Image scaling based on first image
	relative_scaling = PTR_NEW([0,0,0,0])
	immin = DBLARR(hdr.nlp,hdr.ns)
	immax = DBLARR(hdr.nlp,hdr.ns)
	immean = DBLARR(hdr.nlp,hdr.ns)
	imsdev = DBLARR(hdr.nlp,hdr.ns)
	dopplermin = DBLARR(hdr.nlp,hdr.ns)
	dopplermax = DBLARR(hdr.nlp,hdr.ns)
  dopplerscan = FLTARR(hdr.nx,hdr.ny,hdr.nlp*hdr.ns)
	ls_low_y = FLTARR(hdr.ns)
	ls_upp_y = FLTARR(hdr.ns)
	ls_yrange = FLTARR(hdr.ns)
	int_low_y = FLTARR(hdr.ns)
	int_upp_y = FLTARR(hdr.ns)
	FOR j=0,hdr.ns-1 DO BEGIN
		FOR k=0,lp_last DO BEGIN
			temp_image = (*hdr.imdata)[j*hdr.nlp + k]
			immean[k,j] = MEAN(temp_image, /NAN)
			imsdev[k,j] = STDDEV(temp_image, /NAN)
  		immin[k,j] = MIN(temp_image, MAX=max_val, /NAN)
  		immax[k,j] = max_val
			temp_k = 2*hdr.lc - k
			IF ((temp_k EQ hdr.lc) OR (temp_k LT 0) OR (temp_k GT lp_last)) THEN BEGIN
				dopplermin[k,j] = 0
				dopplermax[k,j] = 0
			ENDIF ELSE BEGIN
				mirror_temp_image = (*hdr.imdata)[j*hdr.nlp + temp_k]
  			IF (temp_k LT hdr.lc) THEN $
          dopplerim = temp_image - mirror_temp_image $
  			ELSE $
          dopplerim = mirror_temp_image - temp_image
        dopplerscan[*,*,j*hdr.nlp+k] = dopplerim
  			dopplermin[k,j] = MIN(dopplerim, MAX=max_val, /NAN)
  			dopplermax[k,j] = max_val
      ENDELSE
		ENDFOR
		IF scalestokes THEN BEGIN
;			ls_low_y[j] = MIN(immin[*,j])/hdr.ms
;			ls_upp_y[j] = MAX(immax[*,j])/hdr.ms
			ls_low_y[j] = MIN((immean[*,j]-3.*imsdev[*,j])/hdr.ms, /NAN)
			ls_upp_y[j] = MAX((immean[*,j]+3.*imsdev[*,j])/hdr.ms, /NAN)
		ENDIF ELSE BEGIN
;			ls_low_y[j] = MIN(immin[*,j])/hdr.ms[j]
;			ls_upp_y[j] = MAX(immax[*,j])/hdr.ms[j]
			ls_low_y[j] = MIN((immean[*,j]-3.*imsdev[*,j])/hdr.ms[j], /NAN)
			ls_upp_y[j] = MAX((immean[*,j]+3.*imsdev[*,j])/hdr.ms[j], /NAN)
		ENDELSE
		ls_yrange[j] = ls_upp_y[j] - ls_low_y[j]
		max_imsdev = MAX(imsdev[*,j], /NAN)
		int_low_y[j] = (MEAN(immean[*,j], /NAN)-3.*max_imsdev)/ABS(MEAN(immean[*,j], /NAN))
		int_upp_y[j] = (MEAN(immean[*,j], /NAN)+3.*max_imsdev)/ABS(MEAN(immean[*,j], /NAN))
	ENDFOR
	ls_low_y_init = ls_low_y[0]
	ls_upp_y_init = ls_upp_y[0]
  IF (hdr.ndiagnostics GT 1) THEN BEGIN
    ; Determine multiplicative factors
    ls_range_tmp = FLTARR(hdr.ndiagnostics,hdr.ns)
    FOR j=0,hdr.ns-1 DO BEGIN
      FOR d=0,hdr.ndiagnostics-1 DO BEGIN
        ls_low_tmp=MIN((immin[hdr.diag_start[d]:(hdr.diag_start[d]+(hdr.diag_width[d]-1)),$
          j]-3.*imsdev[hdr.diag_start[d]:(hdr.diag_start[d]+(hdr.diag_width[d]-1)),j])/$
          hdr.ms[j], /NAN)
        ls_upp_tmp=MAX((immax[hdr.diag_start[d]:(hdr.diag_start[d]+(hdr.diag_width[d]-1)),$
          j]+3.*imsdev[hdr.diag_start[d]:(hdr.diag_start[d]+(hdr.diag_width[d]-1)),j])/$
          hdr.ms[j], /NAN)
        ls_range_tmp[d,j] = ls_upp_tmp - ls_low_tmp
      ENDFOR
    ENDFOR
    main_mult_val = REPLICATE(ls_range_tmp[WHERE(ls_range_tmp EQ MAX(ls_range_tmp, /NAN))],$
      hdr.ndiagnostics)/ls_range_tmp
  ENDIF ELSE main_mult_val = REPLICATE(1.,hdr.ns)
	ls_low_y = PTR_NEW(ls_low_y,/NO_COPY)
	ls_upp_y = PTR_NEW(ls_upp_y,/NO_COPY)
	ls_yrange = PTR_NEW(ls_yrange,/NO_COPY)
	int_low_y = PTR_NEW(int_low_y,/NO_COPY)
	int_upp_y = PTR_NEW(int_upp_y,/NO_COPY)

	ls_low_y_ref = 0
	ls_upp_y_ref = 0
	IF hdr.showref THEN BEGIN
		IF (hdr.refnt EQ 0) THEN BEGIN
			refmin = MIN(*hdr.refdata, MAX=refmax, /NAN)
			refmean = MEAN(*hdr.refdata, /NAN)
			refdev = STDDEV(*hdr.refdata, /NAN)
		ENDIF ELSE IF (((hdr.refnt EQ 1) OR (hdr.refnt EQ hdr.mainnt)) AND (hdr.refnlp EQ 1)) THEN BEGIN
			refmin = MIN((*hdr.refdata)[0], MAX=refmax, /NAN)
			refmean = MEAN((*hdr.refdata)[0], /NAN)
			refdev = STDDEV((*hdr.refdata)[0], /NAN)
		ENDIF ELSE BEGIN
			refmin = FLTARR(hdr.refnlp)
			refmax = FLTARR(hdr.refnlp)
			refmean = FLTARR(hdr.refnlp)
			refdev = FLTARR(hdr.refnlp)
			FOR k=0,hdr.refnlp-1 DO BEGIN
				temp_referencefile = (*hdr.refdata)[k]
				refmean[k] = MEAN(temp_referencefile, /NAN)
				refdev[k] = STDDEV(temp_referencefile, /NAN)
				refmin[k] = MIN(temp_referencefile, MAX=max_val, /NAN)
				refmax[k] = max_val
			ENDFOR
		ENDELSE
		IF showrefls THEN BEGIN
			ls_low_y_ref = MIN((refmean-3.*refdev)/hdr.refms, /NAN)
			ls_upp_y_ref = MAX((refmean+3.*refdev)/hdr.refms, /NAN)
		ENDIF
    IF (hdr.nrefdiagnostics GT 1) THEN BEGIN
      ; Determine multiplicative factors
      ls_range_tmp = FLTARR(hdr.nrefdiagnostics,hdr.ns)
      FOR d=0,hdr.nrefdiagnostics-1 DO BEGIN
        ls_low_tmp=MIN((refmin[hdr.refdiag_start[d]:$
          (hdr.refdiag_start[d]+(hdr.refdiag_width[d]-1))]-$
          3.*imsdev[hdr.refdiag_start[d]:(hdr.refdiag_start[d]+(hdr.refdiag_width[d]-1))])/$
          hdr.refms, /NAN)
        ls_upp_tmp=MAX((immax[hdr.refdiag_start[d]:$
          (hdr.refdiag_start[d]+(hdr.refdiag_width[d]-1))]+$
          3.*imsdev[hdr.refdiag_start[d]:(hdr.refdiag_start[d]+(hdr.refdiag_width[d]-1))])/$
          hdr.refms, /NAN)
        ls_range_tmp[d] = ls_upp_tmp - ls_low_tmp
      ENDFOR
      ref_mult_val = REPLICATE(ls_range_tmp[WHERE(ls_range_tmp EQ MAX(ls_range_tmp, /NAN))],$
        hdr.nrefdiagnostics)/ls_range_tmp
    ENDIF ELSE ref_mult_val = 1.
	ENDIF ELSE BEGIN
		refmin = 0
		refmax = 0
    ref_mult_val = 0
	ENDELSE
	ls_yrange_ref = ls_upp_y_ref - ls_low_y_ref 

  sjimin = MIN((*hdr.sjidata)[0], MAX=sjimax_val, /NAN)
  sjimax = sjimax_val
	
  IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initial scaling parameters)', /OPT, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initial scaling parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- OTHER PARAMETERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(other parameters)', /OPT, /OVER
	feedback_text = [feedback_text,'> Other parameters... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	IF KEYWORD_SET(EXTS) THEN exts_set = 1 ELSE exts_set = 0
	IF (hdr.single_cube[0] GE 1) THEN BEGIN
		IF (hdr.mainnt GT 1) THEN BEGIN
			exts_set = 1						; Automatically set EXTS when providing a (single lineposition) 3D cube
      CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'The exact timeslice (EXTS) keyword has been '+$
        'automatically set to enable the drawing of loop paths and extraction of timeslices!', $
        /WARNING, /NO_ROUTINE, /NEWLINE
		ENDIF ELSE exts_set=0
	ENDIF
	refexts_set = (hdr.refspfile NE 1)
	lp_slid_sens = (hdr.nlp GE 2)
	lp_blink_vals_sens = (hdr.nlp GT 2)
	lp_last_slid = (hdr.nlp-1) > 1
	lp_last_blink = (hdr.nlp-1) > 2
	lp_last_vals = hdr.nlp-1

	scale_cubes_vals = [1.,1.]
	IF (N_ELEMENTS(SCALE_CUBES) EQ 1) THEN scale_cubes_vals[0] = scale_cubes ELSE $
    IF (N_ELEMENTS(SCALE_CUBES) EQ 2) THEN scale_cubes_vals = scale_cubes
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(other parameters)', /OPT, /OVER, /DONE,$
                                        REPEAT_STAGE=verbosity[1]
	feedback_text = [feedback_text[0]+'done!',feedback_text[1:N_ELEMENTS(feedback_text)-2],$
    '> Other parameters... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	WAIT,0.1

;================================================================================= SETTING UP WIDGET
;--------------------------------------------------------------------------------- INITIALISE CONTROL PANEL
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(loading BMP buttons)', /WIDGET, /OVER
	feedback_text = ['Setting up widget... ','> Loading BMP buttons... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	bmpbut_search = FILE_SEARCH(dir_buttons, '*.bmp', COUNT=bmpbut_count)
	IF (bmpbut_count EQ 16) THEN BEGIN
		bmpbut_fbwd_idle     = CRISPEX_READ_BMP_BUTTONS('fbwd_idle.bmp',dir_buttons)
		bmpbut_fbwd_pressed  = CRISPEX_READ_BMP_BUTTONS('fbwd_pressed.bmp',dir_buttons)
		bmpbut_bwd_idle      = CRISPEX_READ_BMP_BUTTONS('bwd_idle.bmp',dir_buttons)
		bmpbut_bwd_pressed   = CRISPEX_READ_BMP_BUTTONS('bwd_pressed.bmp',dir_buttons)
		bmpbut_pause_idle    = CRISPEX_READ_BMP_BUTTONS('pause_idle.bmp',dir_buttons)
		bmpbut_pause_pressed = CRISPEX_READ_BMP_BUTTONS('pause_pressed.bmp',dir_buttons)
		bmpbut_fwd_idle      = CRISPEX_READ_BMP_BUTTONS('fwd_idle.bmp',dir_buttons)
		bmpbut_fwd_pressed   = CRISPEX_READ_BMP_BUTTONS('fwd_pressed.bmp',dir_buttons)
		bmpbut_ffwd_idle     = CRISPEX_READ_BMP_BUTTONS('ffwd_idle.bmp',dir_buttons)
		bmpbut_ffwd_pressed  = CRISPEX_READ_BMP_BUTTONS('ffwd_pressed.bmp',dir_buttons)
		bmpbut_loop_idle     = CRISPEX_READ_BMP_BUTTONS('loop_idle.bmp',dir_buttons)
		bmpbut_loop_pressed  = CRISPEX_READ_BMP_BUTTONS('loop_pressed.bmp',dir_buttons)
		bmpbut_cycle_idle    = CRISPEX_READ_BMP_BUTTONS('cycle_idle.bmp',dir_buttons)
		bmpbut_cycle_pressed = CRISPEX_READ_BMP_BUTTONS('cycle_pressed.bmp',dir_buttons)
		bmpbut_blink_idle    = CRISPEX_READ_BMP_BUTTONS('blink_idle.bmp',dir_buttons)
		bmpbut_blink_pressed = CRISPEX_READ_BMP_BUTTONS('blink_pressed.bmp',dir_buttons)
    failed = 0
	ENDIF ELSE BEGIN
		bmpbut_fbwd_idle  = '<<'    & bmpbut_fbwd_pressed  = bmpbut_fbwd_idle
		bmpbut_bwd_idle   = '<'     & bmpbut_bwd_pressed   = bmpbut_bwd_idle
		bmpbut_pause_idle = '||'    & bmpbut_pause_pressed = bmpbut_pause_idle
		bmpbut_fwd_idle   = '>'     & bmpbut_fwd_pressed   = bmpbut_fwd_idle
		bmpbut_ffwd_idle  = '>>'    & bmpbut_ffwd_pressed  = bmpbut_ffwd_idle
		bmpbut_loop_idle  = 'Loop'  & bmpbut_loop_pressed  = bmpbut_loop_idle
		bmpbut_cycle_idle = 'Cycle' & bmpbut_cycle_pressed = bmpbut_cycle_idle
		bmpbut_blink_idle = 'Blink' & bmpbut_blink_pressed = bmpbut_blink_idle
    failed = 1
	ENDELSE
  IF failed THEN donetext = 'failed.' ELSE donetext = 'done!'
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(loading BMP buttons)', /WIDGET, /OVER, /DONE, FAIL=failed
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Loading BMP buttons... '+donetext]
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- INITIALISE CONTROL PANEL
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initialising control panel)', /WIDGET, /OVER
	feedback_text = [feedback_text,'> Initializing control panel... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	cpanel		    = WIDGET_BASE(TITLE = 'CRISPEX'+instance_label+': Control Panel', $
;	cpanel		    = WIDGET_BASE(TITLE = 'CRISPEX: CRIsp SPectral EXplorer'+instance_label, $
    TLB_FRAME_ATTR = 1, /ROW, /FRAME, KILL_NOTIFY = 'CRISPEX_CLOSE_CLEANUP', $
    APP_MBAR = menubar, TLB_SIZE_EVENTS = 1)
  control_panel       = WIDGET_BASE(cpanel, /COLUMN)

  ; CRISPEX
  crispex_menu        = WIDGET_BUTTON(menubar, VALUE='CRISPEX', /MENU, UVALUE='crispex')
	about			          = WIDGET_BUTTON(crispex_menu, VALUE='About CRISPEX', $
                          EVENT_PRO='CRISPEX_ABOUT_WINDOW')
	preferences		      = WIDGET_BUTTON(crispex_menu, VALUE='Preferences', /SEPARATOR, $
                          EVENT_PRO='CRISPEX_PREFERENCES_WINDOW', ACCELERATOR='Ctrl+P')
	exitmenu		        = WIDGET_BUTTON(crispex_menu, VALUE='Quit', /SEPARATOR, $
                          EVENT_PRO='CRISPEX_CLOSE', ACCELERATOR='Ctrl+Q')
  ; Submenus in menu bar
  ; File
	filemenu		        = WIDGET_BUTTON(menubar, VALUE='File', /MENU, UVALUE='file')
;  openmenu            = WIDGET_BUTTON(filemenu, VALUE = 'Open', /MENU)
;  openref             = WIDGET_BUTTON(openmenu, VALUE = 'Reference cube(s)...', $
;                                      EVENT_PRO = 'CRISPEX_IO_OPEN_REFCUBE')
  header_button       = WIDGET_BUTTON(filemenu, VALUE='Show file header(s)', $
                          EVENT_PRO='CRISPEX_DISPLAYS_HEADER', SENSITIVE=($
                            (STRLEN((*hdr.hdrs_main[0])[0]) GT 0) OR $
                            (STRLEN((*hdr.hdrs_ref[0])[0]) GT 0) OR $
                            (STRLEN((*hdr.hdrs_sji[0])[0]) GT 0)))
	restore_session		  = WIDGET_BUTTON(filemenu, VALUE='Load session...', $
                          EVENT_PRO='CRISPEX_SESSION_RESTORE_WINDOW')
	save_session		    = WIDGET_BUTTON(filemenu, VALUE='Save current...', $
                          EVENT_PRO='CRISPEX_SESSION_SAVE_WINDOW')
	save_as_menu		    = WIDGET_BUTTON(filemenu, VALUE='Save as', /SEPARATOR, /MENU)
	save_as_png_menu	  = WIDGET_BUTTON(save_as_menu, VALUE='PNG', /MENU)
	save_as_png_sns		  = WIDGET_BUTTON(save_as_png_menu, VALUE='Snapshot', $
                          EVENT_PRO='CRISPEX_SAVE_PNG_SNAPSHOT')
	save_as_png_all		  = WIDGET_BUTTON(save_as_png_menu, VALUE='All frames', $
                          EVENT_PRO='CRISPEX_SAVE_PNG_ALL_FRAMES', SENSITIVE=(hdr.mainnt GT 1))
	save_as_png_linescan= WIDGET_BUTTON(save_as_png_menu, VALUE='Line scan', $
                          EVENT_PRO='CRISPEX_SAVE_PNG_LINESCAN', SENSITIVE=(hdr.nlp GT 1))
	save_as_jpg_menu	  = WIDGET_BUTTON(save_as_menu, VALUE='JPEG', /MENU)
	save_as_jpg_sns 	  = WIDGET_BUTTON(save_as_jpg_menu, VALUE='Snapshot', $
                          EVENT_PRO='CRISPEX_SAVE_JPEG_SNAPSHOT')
	save_as_jpg_all		  = WIDGET_BUTTON(save_as_jpg_menu, VALUE='All frames', $
                          EVENT_PRO='CRISPEX_SAVE_JPEG_ALL_FRAMES', SENSITIVE=(hdr.mainnt GT 1))
	save_as_jpg_linescan= WIDGET_BUTTON(save_as_jpg_menu, VALUE='Linescan', $
                          EVENT_PRO='CRISPEX_SAVE_JPEG_LINESCAN', SENSITIVE=(hdr.nlp GT 1))
	save_as_mpeg		    = WIDGET_BUTTON(filemenu, VALUE='Save movie...', $
                          EVENT_PRO='CRISPEX_SAVE_MPEG', SENSITIVE=(hdr.mainnt GT 1))
  ; View
  viewmenu            = WIDGET_BUTTON(menubar, VALUE='View', /MENU)
	sh_zoom_in		      = WIDGET_BUTTON(viewmenu, VALUE='Zoom in', EVENT_PRO='CRISPEX_ZOOMFAC_INCR',$
                          ACCELERATOR='Ctrl+Shift+I')
	sh_zoom_out		      = WIDGET_BUTTON(viewmenu, VALUE='Zoom out', EVENT_PRO='CRISPEX_ZOOMFAC_DECR',$
                          ACCELERATOR='Ctrl+Shift+O')
	focus_session_windows= WIDGET_BUTTON(viewmenu, VALUE='Bring all to front', /SEPARATOR, $
                          EVENT_PRO='CRISPEX_DISPLAYS_ALL_TO_FRONT', ACCELERATOR = 'Ctrl+F')
  ; Movie
  moviemenu           = WIDGET_BUTTON(menubar, VALUE='Movie', /MENU)
	sh_fbwd_button		  = WIDGET_BUTTON(moviemenu, VALUE='Step to previous frame', $
                          EVENT_PRO='CRISPEX_PB_FASTBACKWARD', ACCELERATOR='Shift+B')
	sh_backward_button	= WIDGET_BUTTON(moviemenu, VALUE='Play backwards', $
                          EVENT_PRO='CRISPEX_PB_BACKWARD', ACCELERATOR='Shift+Backspace')
	sh_pause_button		  = WIDGET_BUTTON(moviemenu, VALUE='Pause', $
                          EVENT_PRO='CRISPEX_PB_PAUSE', ACCELERATOR='Shift+Space')
	sh_forward_button	  = WIDGET_BUTTON(moviemenu, VALUE='Play forwards', $
                          EVENT_PRO='CRISPEX_PB_FORWARD', ACCELERATOR='Shift+Tab')
	sh_ffwd_button		  = WIDGET_BUTTON(moviemenu, VALUE='Step to next frame', $
                          EVENT_PRO='CRISPEX_PB_FASTFORWARD', ACCELERATOR='Shift+F')
	sh_lp_incr_button 	= WIDGET_BUTTON(moviemenu, VALUE='Main '+$
                          STRLOWCASE(sp_h[heightset])+' position +', $
                          EVENT_PRO='CRISPEX_SLIDER_LP_INCR', ACCELERATOR='Shift+S', /SEPARATOR)
	sh_lp_decr_button 	= WIDGET_BUTTON(moviemenu, VALUE='Main '+$
                          STRLOWCASE(sp_h[heightset])+' position -', $
                          EVENT_PRO='CRISPEX_SLIDER_LP_DECR', ACCELERATOR='Shift+A')
	sh_lp_ref_incr_button= WIDGET_BUTTON(moviemenu, VALUE='Reference '+$
                          STRLOWCASE(sp_h[refheightset])+' position +', $
                          EVENT_PRO='CRISPEX_SLIDER_LP_REF_INCR', ACCELERATOR='Ctrl+S', $
                          SENSITIVE=hdr.showref)
	sh_lp_ref_decr_button= WIDGET_BUTTON(moviemenu, VALUE='Reference '+$
                          STRLOWCASE(sp_h[refheightset])+' position -', $
                          EVENT_PRO='CRISPEX_SLIDER_LP_REF_DECR', ACCELERATOR='Ctrl+A', $
                          SENSITIVE=hdr.showref)
  ; Analysis                       
  analysismenu        = WIDGET_BUTTON(menubar, VALUE='Analysis', /MENU)
  timeslicemenu		    = WIDGET_BUTTON(analysismenu, VALUE = 'Save current space-time diagram', $
                          /MENU, SENSITIVE=0)
	approxmenu		      = WIDGET_BUTTON(timeslicemenu, VALUE = 'Approximated loop', /MENU)
	save_app_slab_but	  = WIDGET_BUTTON(approxmenu, VALUE='All '+STRLOWCASE(sp_h[heightset])+$
                          ' positions', EVENT_PRO='CRISPEX_SAVE_APPROX_LOOPSLAB')
	save_app_slice_but	= WIDGET_BUTTON(approxmenu, VALUE='Current '+STRLOWCASE(sp_h[heightset])+$
                          ' position', EVENT_PRO='CRISPEX_SAVE_APPROX_LOOPSLICE')
	interpolmenu		    = WIDGET_BUTTON(timeslicemenu, VALUE='Interpolated loop', /MENU)
	save_ex_slab_but	  = WIDGET_BUTTON(interpolmenu, VALUE='All '+STRLOWCASE(sp_h[heightset])+$
                          ' positions', EVENT_PRO='CRISPEX_SAVE_EXACT_LOOPSLAB_CHECK')
	save_ex_slice_but	  = WIDGET_BUTTON(interpolmenu, VALUE='Current '+STRLOWCASE(sp_h[heightset])+$
                          ' position', EVENT_PRO='CRISPEX_SAVE_EXACT_LOOPSLICE')
	save_loop_pts		    = WIDGET_BUTTON(analysismenu, VALUE='Save current path for later retrieval', $
                          EVENT_PRO='CRISPEX_SAVE_LOOP_PTS', SENSITIVE=0)
	sel_saved_loop		  = WIDGET_BUTTON(analysismenu, VALUE='Save from selected path(s)', /SEPARATOR, $
                          EVENT_PRO='CRISPEX_RETRIEVE_LOOP_MENU', SENSITIVE=0)
	all_saved_loop		  = WIDGET_BUTTON(analysismenu, VALUE='Save from all paths', /MENU, SENSITIVE=0)
	all_saved_all_pos	  = WIDGET_BUTTON(all_saved_loop, VALUE='At all '+STRLOWCASE(sp_h[heightset])+$
                          ' positions', EVENT_PRO = 'CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLAB')
	all_saved_sel_pos	  = WIDGET_BUTTON(all_saved_loop, VALUE = 'At saved '+$
                          STRLOWCASE(sp_h[heightset])+' position', $
                          EVENT_PRO='CRISPEX_RETRIEVE_LOOP_ALL_LOOPSLICE')
	det_file_loop		    = WIDGET_BUTTON(analysismenu, VALUE='From detection file...', $
                          EVENT_PRO='CRISPEX_RETRIEVE_DET_FILE_MENU')
  ; Help
  helpmenu            = WIDGET_BUTTON(menubar, VALUE='Help', /MENU)
	open_help		        = WIDGET_BUTTON(helpmenu, VALUE='Open online help', $
                          EVENT_PRO='CRISPEX_HELP', ACCELERATOR='Ctrl+H')
	mailbugmenu		      = WIDGET_BUTTON(helpmenu, VALUE='Submit a bug report', $
                          EVENT_PRO='CRISPEX_HELP_MAIL_BUG')
	mailsugmenu		      = WIDGET_BUTTON(helpmenu, VALUE='Submit a suggestion', $
                          EVENT_PRO='CRISPEX_HELP_MAIL_SUGGESTION')
  shortcuts           = WIDGET_BUTTON(helpmenu, VALUE='Show shortcuts', /SEPARATOR, $
                          EVENT_PRO='CRISPEX_HELP_SHORTCUTS')
  ; Help: Developer
	developermenu   		= WIDGET_BUTTON(helpmenu, VALUE = 'Developer', /MENU, /SEPARATOR, $
                          UVALUE='developer')
	sh_runtime_interrupt= WIDGET_BUTTON(developermenu, VALUE='Interrupt', $
                          EVENT_PRO='CRISPEX_INTERRUPT', ACCELERATOR='Ctrl+Shift+C')
	sh_runtime_verb_menu= WIDGET_BUTTON(developermenu, VALUE='Verbosity', /MENU, UVALUE='verbosity')
	sh_verb_0		        = WIDGET_BUTTON(sh_runtime_verb_menu, VALUE='No verbosity', $
                          EVENT_PRO='CRISPEX_VERBOSE_SET', UVALUE=-1, /NO_RELEASE, /CHECKED_MENU, $
                          ACCELERATOR='Shift+0')
	sh_verb_4		        = WIDGET_BUTTON(sh_runtime_verb_menu, VALUE='Basic runtime verbosity', $
                          EVENT_PRO='CRISPEX_VERBOSE_SET', UVALUE=2, /NO_RELEASE, /CHECKED_MENU, $
                          ACCELERATOR='Shift+4')
	sh_verb_8		        = WIDGET_BUTTON(sh_runtime_verb_menu, VALUE='Extended runtime verbosity', $
                          EVENT_PRO='CRISPEX_VERBOSE_SET', UVALUE=3, /NO_RELEASE, /CHECKED_MENU, $
                          ACCELERATOR='Shift+8')
	sh_verb_16		      = WIDGET_BUTTON(sh_runtime_verb_menu, VALUE='Enable playback statistics', $
                          EVENT_PRO='CRISPEX_VERBOSE_SET', UVALUE=4, /NO_RELEASE, /CHECKED_MENU, $
                          ACCELERATOR='Shift+P')
  clear_menu		      = WIDGET_BUTTON(developermenu, VALUE='Clear', /MENU)
	clear_current_estimate= WIDGET_BUTTON(clear_menu, VALUE='Current time estimate', $
                          EVENT_PRO='CRISPEX_CLEAR_CURRENT_ESTIMATE', SENSITIVE=estimate_run)
	clear_current_cpft	= WIDGET_BUTTON(clear_menu, VALUE='CPFT file for current machine', $
                          EVENT_PRO='CRISPEX_CLEAR_CURRENT_CPFT', SENSITIVE=(cpftfilecount EQ 1))
	clear_current_inst	= WIDGET_BUTTON(clear_menu, VALUE='Instance file for current machine', $
                          EVENT_PRO='CRISPEX_CLEAR_CURRENT_INST', SENSITIVE=(instfilecount EQ 1))
	dispwid			        = WIDGET_BUTTON(developermenu, VALUE='Display window IDs', $
                          EVENT_PRO='CRISPEX_DISPWIDS', ACCELERATOR='Ctrl+I', /CHECKED_MENU)

	tab_tlb			= WIDGET_TAB(control_panel, LOCATION=0, MULTILINE=5)

  ; Playback controls
	playback_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Temporal', /COLUMN)
	playback_contr		= WIDGET_BASE(control_panel, /ROW)
	playback_field_basic	= WIDGET_BASE(playback_contr, /FRAME, /GRID_LAYOUT, COLUMN=5)
	playback_field_add	= WIDGET_BASE(playback_contr, /FRAME, $
    GRID_LAYOUT = (bmpbut_count EQ 16), COLUMN=3, EXCLUSIVE = (bmpbut_count NE 16))
	fbwd_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_fbwd_idle, $
    EVENT_PRO = 'CRISPEX_PB_FASTBACKWARD', TOOLTIP = 'Move one frame backward')
	backward_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_bwd_idle, $
    EVENT_PRO = 'CRISPEX_PB_BACKWARD', TOOLTIP = 'Play backward')
	pause_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_pause_pressed, $
    EVENT_PRO = 'CRISPEX_PB_PAUSE', TOOLTIP = 'Pause')
	forward_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_fwd_idle, $
    EVENT_PRO = 'CRISPEX_PB_FORWARD', TOOLTIP = 'Play forward')
	ffwd_button		= WIDGET_BUTTON(playback_field_basic, VALUE = bmpbut_ffwd_idle, $
    EVENT_PRO = 'CRISPEX_PB_FASTFORWARD', TOOLTIP = 'Move one frame forward')
	loop_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_loop_pressed, $
    EVENT_PRO = 'CRISPEX_PB_LOOP', TOOLTIP = 'Loop')
	WIDGET_CONTROL, loop_button, SET_BUTTON = (bmpbut_count NE 16)
	cycle_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_cycle_idle, $
    EVENT_PRO = 'CRISPEX_PB_CYCLE', TOOLTIP = 'Cycle')
	blink_button		= WIDGET_BUTTON(playback_field_add, VALUE = bmpbut_blink_idle, $
    EVENT_PRO = 'CRISPEX_PB_BLINK', TOOLTIP = 'Blink')

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
  t_last += (t_last EQ 1)
	t_step_slid		= WIDGET_SLIDER(playback_tab, TITLE = 'Frame increment', MIN = 1, MAX = t_last, VALUE = t_step, EVENT_PRO = 'CRISPEX_SLIDER_STEP', SENSITIVE = t_slid_sens)
	imref_blink_field	= WIDGET_BASE(playback_tab, /ROW,/NONEXCLUSIVE)
	imref_blink_but		= WIDGET_BUTTON(imref_blink_field, $
    VALUE = 'Blink between main and reference image', EVENT_PRO = 'CRISPEX_DISPLAYS_IMREFBLINK_TOGGLE', SENSITIVE = hdr.showref)
  master_time_base = WIDGET_BASE(playback_tab, /COLUMN, /FRAME)
  master_time_opts = WIDGET_BASE(master_time_base, /ROW)
  master_time_lab   = WIDGET_LABEL(master_time_opts, VALUE='Master time:', /ALIGN_LEFT)
  master_time_labels = ['Main', 'Reference', 'SJI']
  master_time_buts  = CW_BGROUP(master_time_opts, master_time_labels, $
                        BUTTON_UVALUE=INDGEN(N_ELEMENTS(master_time_labels)), IDS=master_time_ids, $
                        /EXCLUSIVE, /ROW, EVENT_FUNC='CRISPEX_BGROUP_MASTER_TIME')
  showdata = [(hdr.showref OR hdr.sjifile), hdr.showref, hdr.sjifile]
  nrasterdims = [SIZE(hdr.tarr_raster_main,/N_DIMENSIONS), $
                  SIZE(hdr.tarr_raster_ref,/N_DIMENSIONS), 2]
  setbutton = [1,0,0]
  FOR i=0,N_ELEMENTS(master_time_labels)-1 DO $
    WIDGET_CONTROL, master_time_ids[i], SENSITIVE=(showdata[i] AND (nrasterdims[i] EQ 2)), $
    SET_BUTTON=setbutton[i]
  IF (nrasterdims[0] EQ 2) THEN $
    toffset_max = N_ELEMENTS(hdr.tarr_raster_main[*,0])-1 $
  ELSE $
    toffset_max = 1
  time_offset_slid  = WIDGET_SLIDER(master_time_base, TITLE='Raster timing offset', $
                        VALUE=hdr.toffset_main, MIN=0, MAX=toffset_max>1,$
                        EVENT_PRO='CRISPEX_SLIDER_TIME_OFFSET', $
                        SENSITIVE=((toffset_max GT 1) AND (TOTAL(showdata) GT 1)), /DRAG)

  ; Spectral controls
	spectral_tab		= WIDGET_BASE(tab_tlb, TITLE = sp_h[heightset], /COLUMN)
	lp_ranges		= WIDGET_BASE(spectral_tab, /COLUMN, /FRAME)
	lp_range_field		= WIDGET_BASE(lp_ranges, /ROW)
	lower_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Lower index:', /ALIGN_LEFT)
	lower_lp_text		= WIDGET_TEXT(lp_range_field, VALUE = STRTRIM(lp_first,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LP_LOW', SENSITIVE = lp_blink_vals_sens)
	upper_lp_label		= WIDGET_LABEL(lp_range_field, VALUE = 'Upper index:', /ALIGN_LEFT)
	upper_lp_text		= WIDGET_TEXT(lp_range_field, VALUE = STRTRIM(lp_last_vals,2),  /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_DISPRANGE_LP_UPP', SENSITIVE = lp_blink_vals_sens)
	reset_lprange_but	= WIDGET_BUTTON(lp_ranges, VALUE = 'Reset '+STRLOWCASE(sp_h[heightset])+' boundaries', EVENT_PRO = 'CRISPEX_DISPRANGE_LP_RESET', SENSITIVE = 0)
	lp_slid			= WIDGET_SLIDER(control_panel, TITLE = 'Main '+STRLOWCASE(sp_h[heightset])+' position', MIN = lp_first, MAX = lp_last_slid, VALUE = lp_start, EVENT_PRO = 'CRISPEX_SLIDER_LP', $
					/DRAG, SENSITIVE = lp_slid_sens)
	lp_speed_slid		= WIDGET_SLIDER(spectral_tab, TITLE = 'Animation speed [blink/s]', MIN = 1, MAX = 100, VALUE = t_speed, EVENT_PRO = 'CRISPEX_SLIDER_SPEED', /DRAG, SENSITIVE = 0)
	lp_blink_slid		= WIDGET_SLIDER(spectral_tab, $
    TITLE=sp_h[heightset]+' position to blink against', MIN=lp_first, MAX=lp_last_slid, $
    VALUE=lp_start, EVENT_PRO='CRISPEX_SLIDER_SPECTBLINK', /DRAG, SENSITIVE = lp_blink_vals_sens)
	lp_blink_field		= WIDGET_BASE(spectral_tab, /ROW,/NONEXCLUSIVE)
	lp_blink_but		= WIDGET_BUTTON(lp_blink_field, VALUE = 'Blink between '+STRLOWCASE(sp_h[heightset])+' positions', EVENT_PRO = 'CRISPEX_PB_SPECTBLINK', SENSITIVE = lp_slid_sens)
	lp_ref_but_field	= WIDGET_BASE(spectral_tab, /ROW, /NONEXCLUSIVE)
	IF (heightset NE refheightset) THEN reflab = STRLOWCASE(sp_h[refheightset])+' ' ELSE reflab = ''
	lp_ref_but		= WIDGET_BUTTON(lp_ref_but_field, VALUE = 'Lock reference '+reflab+'to main '+$
                  STRLOWCASE(sp_h[heightset])+' position', EVENT_PRO = 'CRISPEX_SLIDER_LP_REF_LOCK', $
					        SENSITIVE = lp_ref_lock) ;(eqnlps AND (hdr.refnlp GT 1)))
	WIDGET_CONTROL, lp_ref_but, SET_BUTTON = lp_ref_lock ;(eqnlps AND (hdr.refnlp GT 1) AND $
;                                           (lp_start EQ lp_ref_start))
	lp_ref_slid = WIDGET_SLIDER(spectral_tab, TITLE = 'Reference '+STRLOWCASE(sp_h[refheightset])+' position', MIN = lp_ref_first, MAX = lp_ref_last, VALUE = lp_ref_start, EVENT_PRO = 'CRISPEX_SLIDER_LP_REF', $
					/DRAG, SENSITIVE = (refslid_sens AND ABS(eqnlps-1)))

  ; Spatial controls
	spatial_tab		= WIDGET_BASE(tab_tlb, TITLE='Spatial', /COLUMN)
	cursor_frame	= WIDGET_BASE(spatial_tab, /FRAME, /COLUMN)
	x_slid			  = WIDGET_SLIDER(cursor_frame, TITLE='X position of the cursor [pixel]', $
                    MIN=x_first, MAX=(x_last > 1), VALUE=x_start, $
                    EVENT_PRO='CRISPEX_SLIDER_X', /DRAG, SENSITIVE=(x_last GT x_first))
	y_slid			  = WIDGET_SLIDER(cursor_frame, TITLE='Y position of the cursor [pixel]', $
                    MIN=y_first, MAX=(y_last > 1), VALUE=y_start, $
                    EVENT_PRO = 'CRISPEX_SLIDER_Y', /DRAG, SENSITIVE=(y_last GT y_first))

	lock_field		= WIDGET_BASE(control_panel, /FRAME, /ROW, /EXCLUSIVE)
	lockbut			= WIDGET_BUTTON(lock_field, VALUE = 'Lock to position', EVENT_PRO = 'CRISPEX_CURSOR_LOCK', TOOLTIP = 'Lock cursor to current position')
	WIDGET_CONTROL, lockbut, SET_BUTTON = 0
	unlockbut		= WIDGET_BUTTON(lock_field, VALUE = 'Unlock from position', TOOLTIP = 'Unlock cursor from current position')
	WIDGET_CONTROL, unlockbut, SET_BUTTON = 1

	zoom_frame		  = WIDGET_BASE(cursor_frame, /ROW)
	zoom_label		  = WIDGET_LABEL(zoom_frame, VALUE = 'Zoom:', /ALIGN_LEFT)
	zoom_but_field  = WIDGET_BASE(zoom_frame, /ROW )
	zoom_buts	      = CW_BGROUP(zoom_but_field,STRTRIM(FIX(zoomfactors),2)+REPLICATE('x',$
                      N_ELEMENTS(zoomfactors)), BUTTON_UVALUE=INDGEN(N_ELEMENTS(zoomfactors)), $
                      IDS=zoom_button_ids,/EXCLUSIVE, /ROW, EVENT_FUNC = 'CRISPEX_BGROUP_ZOOMFAC_SET')

  ; Diagnostics tab: Stokes and spectral windows
;  stokes_tab		= WIDGET_BASE(tab_tlb, TITLE='Stokes', /COLUMN)
  diagnostics_tab = WIDGET_BASE(tab_tlb, TITLE='Diagnostics', /COLUMN)
  ; Stokes part
	stokes_frame    = WIDGET_BASE(diagnostics_tab, /FRAME, /COLUMN)
	stokes_disp_label		= WIDGET_LABEL(stokes_frame, $
    VALUE = 'Stokes parameter:                                    ', /ALIGN_LEFT)
	stokes_main			= WIDGET_BASE(stokes_frame, /ROW)
	stokes_main_label		= WIDGET_LABEL(stokes_main, VALUE = 'Main image:',/ALIGN_LEFT)
	stokes_xy_but_field= WIDGET_BASE(stokes_main, /ROW )
  stokes_button_labels = ['I','Q','U','V']
  stokes_xy_buts     = CW_BGROUP(stokes_xy_but_field, stokes_button_labels, $
                      BUTTON_UVALUE=INDGEN(N_ELEMENTS(stokes_button_labels)), IDS=stokes_button_ids,$
                      /EXCLUSIVE, /ROW, EVENT_FUNC = 'CRISPEX_BGROUP_STOKES_SELECT_XY')
	stokes_sp			= WIDGET_BASE(stokes_frame, /ROW)
	stokes_sp_label		= WIDGET_LABEL(stokes_sp, VALUE = 'Detailed spectra:',/ALIGN_LEFT)
  stokes_sp_buts     = CW_BGROUP(stokes_sp, stokes_button_labels, $
                      BUTTON_UVALUE=INDGEN(N_ELEMENTS(stokes_button_labels)), IDS=stokes_spbutton_ids,$
                      /NONEXCLUSIVE, /ROW, EVENT_FUNC = 'CRISPEX_BGROUP_STOKES_SELECT_SP')
	spconstraint		= (hdr.nlp GT 1)
  FOR i=0,N_ELEMENTS(stokes_button_labels)-1 DO BEGIN
    WIDGET_CONTROL, stokes_button_ids[i], SENSITIVE=hdr.stokes_enabled[i], SET_BUTTON=(i EQ 0)
    IF (hdr.multichannel OR (i GT 0)) THEN $
      set_constraint = (spconstraint AND hdr.stokes_enabled[i]) $
    ELSE $
      set_constraint = spconstraint
    WIDGET_CONTROL, stokes_spbutton_ids[i], SENSITIVE=(spconstraint AND hdr.stokes_enabled[i]), $
      SET_BUTTON=set_constraint
  ENDFOR
  ; Spectral window part
  specwin_frame       = WIDGET_BASE(diagnostics_tab, /FRAME, /COLUMN)
  specwin_disp_label  = WIDGET_LABEL(specwin_frame, $
    VALUE = 'Spectral windows:                                    ', /ALIGN_LEFT)
  specwin_sub_frame   = WIDGET_BASE(specwin_frame, /GRID_LAYOUT, COLUMN=2)
  main_select_base    = WIDGET_BASE(specwin_sub_frame,/COLUMN,/FRAME)
	main_specwin_label  = WIDGET_LABEL(main_select_base, VALUE = 'Main:', /ALIGN_LEFT)
  IF (hdr.ndiagnostics GT 1) THEN $
    vals = ['Display all',hdr.diagnostics] $
  ELSE $
    vals = 'N/A        '
  specwin_buts        = CW_BGROUP(main_select_base, vals, $
                          BUTTON_UVALUE=INDGEN(N_ELEMENTS(vals)), IDS=specwin_button_ids, $
                          /NONEXCLUSIVE, /COLUMN, EVENT_FUNC='CRISPEX_BGROUP_DIAGNOSTICS_SELECT')
  FOR i=0,N_ELEMENTS(vals)-1 DO $
    WIDGET_CONTROL, specwin_button_ids[i], SENSITIVE=(i GT 0), /SET_BUTTON
  ref_select_base    = WIDGET_BASE(specwin_sub_frame,/COLUMN,/FRAME)
	ref_specwin_label  = WIDGET_LABEL(ref_select_base, VALUE = 'Reference:', /ALIGN_LEFT)
  IF (hdr.nrefdiagnostics GT 1) THEN $
    vals = ['Display all',hdr.refdiagnostics] $
  ELSE $
    vals = 'N/A        '
  refspecwin_buts   = CW_BGROUP(ref_select_base, vals, $
                        BUTTON_UVALUE=INDGEN(N_ELEMENTS(vals)), IDS=refspecwin_button_ids, $
                        /NONEXCLUSIVE, /COLUMN, EVENT_FUNC='CRISPEX_BGROUP_REFDIAGNOSTICS_SELECT')
  FOR i=0,N_ELEMENTS(vals)-1 DO $
    WIDGET_CONTROL, refspecwin_button_ids[i], SENSITIVE=(i GT 0), /SET_BUTTON

	display_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Displays', /COLUMN)
	detspect_frame		= WIDGET_BASE(display_tab, /FRAME, /COLUMN)
	detspect_label_imref	= WIDGET_BASE(detspect_frame, /ROW)
	detspect_label		= WIDGET_LABEL(detspect_label_imref, VALUE = lswintitle[heightset]+':',/ALIGN_LEFT, /DYNAMIC_RESIZE)
	detspect_imref		= WIDGET_BASE(detspect_label_imref, /ROW, /EXCLUSIVE)
	detspect_im_but		= WIDGET_BUTTON(detspect_imref, VALUE = 'Main', EVENT_PRO = $
  'CRISPEX_DISPLAYS_DETSPECT_IM_SELECT', /NO_RELEASE, SENSITIVE = (hdr.nlp GT 1), $
					TOOLTIP = 'Main '+STRLOWCASE(lswintitle[heightset])+' display options')
	WIDGET_CONTROL, detspect_im_but, SET_BUTTON = (hdr.nlp GT 1)
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
	int_toggle_but		= WIDGET_BUTTON(other_buts, VALUE = 'I-t', EVENT_PRO = $
  'CRISPEX_DISPLAYS_INT_TOGGLE', SENSITIVE=(hdr.mainnt GT 1), TOOLTIP = 'Toggle display intensity versus time plot')
	images_disp		= WIDGET_BASE(all_other_disp, /ROW)
	images_label		= WIDGET_LABEL(images_disp, VALUE = 'Images:')
	images_buts		= WIDGET_BASE(images_disp, /ROW, /NONEXCLUSIVE)
	reference_but		= WIDGET_BUTTON(images_buts, VALUE = 'Reference', $
    EVENT_PRO = 'CRISPEX_DISPLAYS_REF_TOGGLE', SENSITIVE = hdr.showref, TOOLTIP = 'Toggle display reference image')
	WIDGET_CONTROL, reference_but, SET_BUTTON = hdr.showref
	doppler_but		= WIDGET_BUTTON(images_buts, VALUE = 'Doppler', EVENT_PRO = $
  'CRISPEX_DISPLAYS_DOPPLER_TOGGLE', SENSITIVE = (hdr.nlp GT 1), TOOLTIP = 'Toggle display Doppler image')
  sji_but       = WIDGET_BUTTON(images_buts, VALUE = 'Slit-jaw', EVENT_PRO = $
    'CRISPEX_DISPLAYS_SJI_TOGGLE', SENSITIVE = hdr.sjifile, TOOLTIP='Toggle display slit-jaw image')
  WIDGET_CONTROL, sji_but, SET_BUTTON = hdr.sjifile

	scaling_tab 		= WIDGET_BASE(tab_tlb, TITLE = 'Scaling', /COLUMN)
	xy_scaling		  = WIDGET_BASE(scaling_tab, /FRAME, /COLUMN)
  scaling_cbox    = WIDGET_COMBOBOX(xy_scaling, $
                    VALUE=['Main image','Reference image','Doppler image','Slit-jaw image'], $
                    EVENT_PRO='CRISPEX_SCALING_SELECT_DATA')
	xy_scale_opts		= WIDGET_BASE(xy_scaling, /COLUMN)
  imagescale_cbox = WIDGET_COMBOBOX(xy_scaling, $
    VALUE=['Based on first image','Based on current image','Per time step'],$
    EVENT_PRO='CRISPEX_SCALING_SELECT_TYPE')
  diagscale_label_vals = REPLICATE('Spectral window: ',2*hdr.ndiagnostics+hdr.nrefdiagnostics+1)+$
    [hdr.diagnostics,hdr.refdiagnostics,hdr.diagnostics,'N/A']
  diagscale_label = WIDGET_LABEL(xy_scaling, VALUE=diagscale_label_vals[0], /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  histo_base = WIDGET_BASE(xy_scaling, /ROW)
  histo_opt_label = WIDGET_LABEL(histo_base, VALUE='Histogram optimisation', /ALIGN_LEFT)
  histo_opt_txt   = WIDGET_TEXT(histo_base, VALUE=STRTRIM(histo_opt_val,2), /EDITABLE, $
    XSIZE=11, EVENT_PRO='CRISPEX_SCALING_HISTO_OPT_VALUE')
  minmax_sliders = WIDGET_BASE(xy_scaling, /ROW)
  scalemin_slider = WIDGET_SLIDER(minmax_sliders, TITLE='Image minimum [%]', MIN=0, MAX=99, $
    VALUE=0, EVENT_PRO='CRISPEX_SCALING_SLIDER_MIN', /DRAG, XSIZE=131)
  scalemax_slider = WIDGET_SLIDER(minmax_sliders, TITLE='Image maximum [%]', MIN=1, MAX=100, $
    VALUE=100, EVENT_PRO='CRISPEX_SCALING_SLIDER_MAX', /DRAG, XSIZE=131)
  gamma_label = WIDGET_LABEL(xy_scaling, VALUE=STRING(gamma_val, FORMAT='(F6.3)'), $
    /ALIGN_CENTER,XSIZE=250)
  gamma_slider = WIDGET_SLIDER(xy_scaling, TITLE='Gamma', MIN=0, MAX=1000, $
    VALUE=500*(ALOG10(gamma_val)+1), EVENT_PRO='CRISPEX_SCALING_GAMMA_SLIDER', /SUPPRESS, $
    /DRAG,XSIZE=250)
  reset_buts = WIDGET_BASE(xy_scaling, /ROW, /GRID, /ALIGN_CENTER)
  scaling_reset_but = WIDGET_BUTTON(reset_buts, VALUE='Reset current', $
    EVENT_PRO='CRISPEX_SCALING_RESET_DEFAULTS', TOOLTIP='Reset scaling of current diagnostic to '+$
    'defaults')
  scaling_reset_all_but = WIDGET_BUTTON(reset_buts, VALUE='Reset all', $
    EVENT_PRO='CRISPEX_SCALING_RESET_ALL_DEFAULTS', TOOLTIP='Reset scaling of all diagnostics to '+$
    'defaults', SENSITIVE=(hdr.ndiagnostics GT 1))
  ls_scaling    = WIDGET_BASE(scaling_tab, /FRAME,/COLUMN)
  ls_scale_opts = WIDGET_BASE(ls_scaling, /COLUMN)
  ls_scale_label= WIDGET_LABEL(ls_scale_opts, VALUE='Detailed spectrum:', /ALIGN_LEFT)
  ls_mult_opts  = WIDGET_BASE(ls_scale_opts, /ROW)
  ls_mult_label = WIDGET_LABEL(ls_mult_opts, VALUE='Multiply', /ALIGN_LEFT, $
    SENSITIVE=((hdr.ndiagnostics GT 1) OR (hdr.nrefdiagnostics GT 1)))
  IF (hdr.refdiagnostics[0] NE 'N/A') THEN $
    ls_mult_list  = [REPLICATE('Main ',hdr.ndiagnostics)+hdr.diagnostics, $
                      REPLICATE('Reference ',hdr.nrefdiagnostics)+hdr.refdiagnostics] $
  ELSE $
    ls_mult_list  = [REPLICATE('Main ',hdr.ndiagnostics)+hdr.diagnostics]
  ls_mult_cbox  = WIDGET_COMBOBOX(ls_mult_opts, VALUE=ls_mult_list, $
    EVENT_PRO='CRISPEX_SCALING_MULTIPLY_LS_SELECT', /DYNAMIC_RESIZE, $
    SENSITIVE=((hdr.ndiagnostics GT 1) OR (hdr.nrefdiagnostics GT 1)))
  ls_mult_by    = WIDGET_LABEL(ls_mult_opts, VALUE='by', /ALIGN_CENTER, $
    SENSITIVE=((hdr.ndiagnostics GT 1) OR (hdr.nrefdiagnostics GT 1)))
  ls_mult_txt   = WIDGET_TEXT(ls_mult_opts, VALUE=STRTRIM(main_mult_val[0],2), /EDITABLE, $
    XSIZE=5, EVENT_PRO='CRISPEX_SCALING_MULTIPLY_LS_VALUE', $
    SENSITIVE=((hdr.ndiagnostics GT 1) OR (hdr.nrefdiagnostics GT 1)))
	
	slit_tab		= WIDGET_BASE(tab_tlb, TITLE = 'Slits',/COLUMN)
	slit_frame		= WIDGET_BASE(slit_tab, /FRAME, /COLUMN)
	slit_label		= WIDGET_LABEL(slit_frame, VALUE = 'Slit controls:', /ALIGN_LEFT)
	phi_slid		= WIDGET_SLIDER(slit_frame, TITLE = 'Slit angle [degrees]', MIN = 0, MAX = 179, $
    VALUE=angle, EVENT_PRO = 'CRISPEX_SLIDER_PHI_ANGLE', SENSITIVE = 0, /DRAG)
	nphi_slid		= WIDGET_SLIDER(slit_frame, TITLE = 'Slit length [pixel]', MIN = 2, MAX = nphi, $
    VALUE=LONG(hdr.ny/3.), EVENT_PRO = 'CRISPEX_SLIDER_NPHI', SENSITIVE = 0, /DRAG)
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

	masks_tab		            = WIDGET_BASE(tab_tlb, TITLE = 'Mask', /COLUMN)
	masks			              = WIDGET_BASE(masks_tab, /FRAME, /COLUMN)
	masks_overlay		        = WIDGET_BASE(masks, /ROW)
	masks_overlay_label	    = WIDGET_LABEL(masks_overlay, VALUE = 'Overlay on:',/ALIGN_LEFT)
	masks_overlay_buts	    = CW_BGROUP(masks_overlay, ['Main','Reference','Doppler'],$
                              BUTTON_UVALUE=INDGEN(3),IDS=mask_button_ids,/NONEXCLUSIVE, /ROW, $
                              EVENT_FUNC = 'CRISPEX_BGROUP_MASK_OVERLAY')
	masks_overlay_color	    = WIDGET_BASE(masks,/COLUMN)
	LOADCT,GET_NAMES=ctnames,/SILENT
	masks_overlay_ct_cbox	  = WIDGET_COMBOBOX(masks_overlay_color, $
                              VALUE = STRTRIM(INDGEN(N_ELEMENTS(ctnames)),2)+REPLICATE(': ',$
                              N_ELEMENTS(ctnames))+ctnames, $
                              EVENT_PRO = 'CRISPEX_MASK_OVERLAY_SELECT_COLOR_TABLE', $
                              SENSITIVE = maskfile)
	maskct                  = 13
	WIDGET_CONTROL, masks_overlay_ct_cbox, SET_COMBOBOX_SELECT = maskct
	masks_overlay_col_slid  = WIDGET_SLIDER(masks_overlay_color, MIN = 0, MAX = 255, VALUE = 255, $
                              TITLE = 'Color index', EVENT_PRO='CRISPEX_MASK_OVERLAY_COLOR_SLIDER',$
                              /DRAG, SENSITIVE = maskfile)
  raster_base = WIDGET_BASE(masks_tab,/FRAME,/COLUMN)
  raster_overlay_label = WIDGET_LABEL(raster_base, VALUE='Raster:',/ALIGN_LEFT)
  raster_overlay = WIDGET_BASE(raster_base, /ROW, /NONEXCLUSIVE)
  raster_but = WIDGET_BUTTON(raster_overlay, VALUE='Overlay boundaries on slit-jaw image', $
    EVENT_PRO='CRISPEX_MASK_OVERLAY_RASTER_TOGGLE', SENSITIVE=hdr.sjifile)
  WIDGET_CONTROL, raster_but, SET_BUTTON=hdr.sjifile

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
  IF hdr.dx_fixed THEN $
  	apix_text		= WIDGET_LABEL(apix_base, VALUE = STRTRIM(hdr.dx,2), /ALIGN_LEFT, SENSITIVE = 0) $
  ELSE $
  	apix_text		= WIDGET_TEXT(apix_base, VALUE = STRTRIM(hdr.dx,2), /EDITABLE, XSIZE = 5, EVENT_PRO = 'CRISPEX_MEASURE_ARCSEC', SENSITIVE = 0)
	apix_unit		= WIDGET_LABEL(apix_base, VALUE = '['+hdr.xunit+']', /ALIGN_LEFT, SENSITIVE = 0)
	measure_asec		= WIDGET_BASE(measuretool, /ROW)
	measure_asec_lab	= WIDGET_LABEL(measure_asec, VALUE = 'Distance [arcsec]:', /ALIGN_LEFT, SENSITIVE = 0)
	measure_asec_text	= WIDGET_LABEL(measure_asec, VALUE = '0.00', /DYNAMIC_RESIZE, SENSITIVE = 0)
	measure_km		= WIDGET_BASE(measuretool, /ROW)
	measure_km_lab		= WIDGET_LABEL(measure_km, VALUE = 'Distance [km]:', /ALIGN_LEFT, SENSITIVE = 0)
	measure_km_text		= WIDGET_LABEL(measure_km, VALUE = '0.00', /DYNAMIC_RESIZE, SENSITIVE = 0)

  ; Parameters overview
    ; Position parameters
    params_position_base = WIDGET_BASE(control_panel, /ROW,/FRAME, /GRID_LAYOUT)
    verlabel_base = WIDGET_BASE(params_position_base, /COLUMN)
      no_label    = WIDGET_LABEL(verlabel_base, VALUE='Position', /ALIGN_LEFT)
      pixel_label = WIDGET_LABEL(verlabel_base, VALUE='Index [px]', /ALIGN_RIGHT)
      real_label  = WIDGET_LABEL(verlabel_base, VALUE='Value ["]', /ALIGN_RIGHT)
    params_main_base = WIDGET_BASE(params_position_base, /COLUMN)
      main_label  = WIDGET_LABEL(params_main_base, VALUE='Main', /ALIGN_RIGHT)
        xcoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.nx))+1,2)+')'
        ycoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.ny))+1,2)+')'
        coord_txt = '    ('+STRING(LONG(x_start),FORMAT=xcoord_format)+$
          ','+STRING(LONG(y_start),FORMAT=ycoord_format)+')'
		    xycoord_val = WIDGET_LABEL(params_main_base, VALUE=coord_txt, $
          /ALIGN_RIGHT)
        xcoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.nx*hdr.dx))+3+(hdr.nx*hdr.dx LT 1),2)+'.1)'
        ycoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.ny*hdr.dy))+3+(hdr.ny*hdr.dy LT 1),2)+'.1)'
        real_coord_txt = '('+STRING(FLOAT(x_start*hdr.dx),FORMAT=xcoord_real_format)+$
          ','+STRING(FLOAT(y_start*hdr.dy),FORMAT=ycoord_real_format)+')'
		    xycoord_real_val = WIDGET_LABEL(params_main_base, VALUE=real_coord_txt, $
          /ALIGN_RIGHT)
    params_ref_base = WIDGET_BASE(params_position_base, /COLUMN, /ALIGN_RIGHT)
      ref_label   = WIDGET_LABEL(params_ref_base, VALUE='Reference', /ALIGN_RIGHT)
      IF hdr.showref THEN BEGIN
        refxcoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.refnx))+1,2)+')'
        refycoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.refny))+1,2)+')'
        refcoord_txt = '    ('+STRING(LONG(x_start),FORMAT=refxcoord_format)+$
          ','+STRING(LONG(y_start),FORMAT=refycoord_format)+')'
      ENDIF ELSE BEGIN
        refcoord_txt = 'N/A'
        refxcoord_format = ''
        refycoord_format = ''
      ENDELSE
		    refxycoord_val = WIDGET_LABEL(params_ref_base, VALUE=refcoord_txt, $
          /ALIGN_RIGHT)
      IF hdr.showref THEN BEGIN
        refxcoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.refnx*hdr.dx))+3+$
          (hdr.refnx*hdr.dx LT 1),2)+'.1)'
        refycoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.refny*hdr.dy))+3+$
          (hdr.refny*hdr.dy LT 1),2)+'.1)'
        refcoord_real_txt = '('+STRING(FLOAT(x_start*hdr.dx),FORMAT=refxcoord_real_format)+$
          '",'+STRING(FLOAT(y_start*hdr.dy),FORMAT=refycoord_real_format)+'")'
      ENDIF ELSE BEGIN
        refcoord_real_txt = 'N/A'
        refxcoord_real_format = ''
        refycoord_real_format = ''
      ENDELSE
		    refxycoord_real_val = WIDGET_LABEL(params_ref_base, VALUE=refcoord_real_txt, $
          /ALIGN_RIGHT)
    params_sji_base = WIDGET_BASE(params_position_base, /COLUMN, /ALIGN_RIGHT)
      sji_label   = WIDGET_LABEL(params_sji_base, VALUE='Slit-jaw', /ALIGN_RIGHT)
      IF hdr.sjifile THEN BEGIN
        sjixcoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.sjinx))+1,2)+')'
        sjiycoord_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.sjiny))+1,2)+')'
        sjicoord_txt = '    ('+STRING(LONG(xsji_start),FORMAT=sjixcoord_format)+$
          ','+STRING(LONG(ysji_start),FORMAT=sjiycoord_format)+')'
      ENDIF ELSE BEGIN
        sjicoord_txt = 'N/A'
        sjixcoord_format = ''
        sjiycoord_format = ''
      ENDELSE
		    sjixycoord_val = WIDGET_LABEL(params_sji_base, VALUE=sjicoord_txt, $
          /ALIGN_RIGHT)
      IF hdr.sjifile THEN BEGIN
        sjixcoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.sjinx*hdr.sjidx))+3,2)+'.1)'
        sjiycoord_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.sjiny*hdr.sjidy))+3,2)+'.1)'
        sjicoord_real_txt = '('+STRING(FLOAT(xsji_start*hdr.dx),$
          FORMAT=sjixcoord_real_format)+'",'+STRING(FLOAT(ysji_start*hdr.dy),$
          FORMAT=sjiycoord_real_format)+'")'
      ENDIF ELSE BEGIN
        sjicoord_real_txt = 'N/A'
        sjixcoord_real_format = ''
        sjiycoord_real_format = ''
      ENDELSE
		    sjixycoord_real_val = WIDGET_LABEL(params_sji_base, VALUE=sjicoord_real_txt, $
          /ALIGN_RIGHT)

    ; Spectral parameters
;      divider_label = WIDGET_LABEL(verlabel_base, VALUE=' ')
      no_label    = WIDGET_LABEL(verlabel_base, VALUE=wav_h[heightset[0]], /ALIGN_LEFT)
      pixel_label = WIDGET_LABEL(verlabel_base, VALUE='Index [px]', /ALIGN_RIGHT)
		  IF ((TOTAL(heightset) GE 1) OR (TOTAL(hdr.v_dop_set) GE 1)) THEN $
        real_label  = WIDGET_LABEL(verlabel_base, VALUE='Value', /ALIGN_RIGHT)
      IF ((heightset[0] EQ 0) AND (TOTAL(hdr.v_dop_set) GE 1)) THEN $
        vdop_label  = WIDGET_LABEL(verlabel_base, VALUE='Doppler [km/s]', /ALIGN_RIGHT)
      main_label  = WIDGET_LABEL(params_main_base, VALUE=' ', /ALIGN_RIGHT)
      lp_idx_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.nlp))+1,2)+')'
      IF hdr.v_dop_set[0] THEN BEGIN
        lp_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.lps[lp_last]))+3,2)+'.1)'
        lp_vdop_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.lps[lp_last]))+4,2)+'.2)'
        lp_real_txt = STRING(hdr.lps[lp_start], FORMAT=lp_real_format)
        lp_vdop_txt = STRING((*hdr.v_dop[0])[lp_start-hdr.diag_start[0]],$
          FORMAT=lp_vdop_format)
      ENDIF ELSE BEGIN
        lp_real_format = ''
        lp_vdop_format = ''
        lp_real_val = 0
        lp_vdop_val = 0
        lp_real_txt = 'N/A'
        lp_vdop_txt = 'N/A'
      ENDELSE
		  lp_idx_val = WIDGET_LABEL(params_main_base, VALUE=STRING(LONG(lp_start),$
        FORMAT=lp_idx_format), /ALIGN_RIGHT)
		  IF ((TOTAL(heightset) GE 1) OR (TOTAL(hdr.v_dop_set) GE 1)) THEN $
		    lp_real_val = WIDGET_LABEL(params_main_base, VALUE=lp_real_txt, /ALIGN_RIGHT)
		  IF (TOTAL(hdr.v_dop_set) GE 1) THEN $
        lp_vdop_val = WIDGET_LABEL(params_main_base, VALUE=lp_vdop_txt, /ALIGN_RIGHT)
      ref_label   = WIDGET_LABEL(params_ref_base, VALUE=' ', /ALIGN_RIGHT)
      IF (hdr.refnlp GT 1) THEN BEGIN
        lp_ref_idx_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.refnlp))+1,2)+')'
        lp_ref_idx_txt = STRING(LONG(lp_ref_start), FORMAT=lp_ref_idx_format)
      ENDIF ELSE BEGIN
        lp_ref_idx_format = ''
        lp_ref_idx_txt = 'N/A'
      ENDELSE
      IF ((hdr.refnlp GT 1) AND hdr.v_dop_set[1]) THEN BEGIN
        lp_ref_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.reflps[lp_ref_last]))+3,2)+'.1)'
        lp_ref_vdop_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.reflps[lp_ref_last]))+4,2)+'.2)'
        lp_ref_real_txt = STRING(hdr.reflps[lp_ref_start], FORMAT=lp_ref_real_format)
        lp_ref_vdop_txt = STRING((*hdr.v_dop_ref[0])[lp_ref_start-hdr.diag_start[0]],$
          FORMAT=lp_ref_vdop_format)
      ENDIF ELSE BEGIN
        lp_ref_real_format = ''
        lp_ref_vdop_format = ''
        lp_ref_real_val = 0
        lp_ref_vdop_val = 0
        lp_ref_real_txt = 'N/A'
        lp_ref_vdop_txt = 'N/A'
      ENDELSE
		  lp_ref_idx_val = WIDGET_LABEL(params_ref_base, VALUE=lp_ref_idx_txt, /ALIGN_RIGHT)
		  IF ((TOTAL(heightset) GE 1) OR (TOTAL(hdr.v_dop_set) GE 1)) THEN $
        lp_ref_real_val = WIDGET_LABEL(params_ref_base, VALUE=lp_ref_real_txt, $
          /ALIGN_RIGHT)
		  IF (TOTAL(hdr.v_dop_set) GE 1) THEN $
        lp_ref_vdop_val = WIDGET_LABEL(params_ref_base, VALUE=lp_ref_vdop_txt, $
          /ALIGN_RIGHT)
      ; SJI placeholders
      sji_label   = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
		  lp_sji_idx_val = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
		  IF ((TOTAL(heightset) GE 1) OR (TOTAL(hdr.v_dop_set) GE 1)) THEN $
        lp_sji_real_val = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
		  IF (TOTAL(hdr.v_dop_set) GE 1) THEN $
        lp_sji_vdop_val = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
    
    ; Time parameters
;    params_time_base = WIDGET_BASE(control_panel, /ROW,/FRAME)
;    verlabel_base = WIDGET_BASE(params_time_base, /COLUMN)
      no_label    = WIDGET_LABEL(verlabel_base, VALUE='Time', /ALIGN_LEFT)
      pixel_label = WIDGET_LABEL(verlabel_base, VALUE='Index [px]', /ALIGN_RIGHT)
      IF dt_set THEN $
        real_label  = WIDGET_LABEL(verlabel_base, VALUE='Value [s]', /ALIGN_RIGHT)
      raster_time_fb = (N_ELEMENTS(hdr.tarr_raster_main) NE N_ELEMENTS(hdr.tarr_main)) 
      refraster_time_fb = (N_ELEMENTS(hdr.tarr_raster_ref) NE N_ELEMENTS(hdr.tarr_ref))
      IF (raster_time_fb OR refraster_time_fb) THEN $
        raster_label = WIDGET_LABEL(verlabel_base, VALUE='Raster [s]', /ALIGN_RIGHT)
;    params_main_base = WIDGET_BASE(params_time_base, /COLUMN)
      main_label  = WIDGET_LABEL(params_main_base, VALUE=' ', /ALIGN_RIGHT)
      t_idx_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.mainnt))+1,2)+')'
		  t_idx_val = WIDGET_LABEL(params_main_base, VALUE=STRING(LONG(t_start),$
        FORMAT=t_idx_format), /ALIGN_RIGHT)
      IF dt_set THEN BEGIN
        t_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.tarr_main[$
          (WHERE(hdr.tarr_main GT 0))[-1]]))+3,2)+'.1)'
        t_real_txt = STRING(hdr.tarr_main[t_start], FORMAT=t_real_format)
		    t_real_val = WIDGET_LABEL(params_main_base, VALUE=t_real_txt, /ALIGN_RIGHT)
      ENDIF ELSE BEGIN
        t_real_format = ''
        t_real_txt = 'N/A'
        t_real_val = 0
      ENDELSE
      IF raster_time_fb THEN BEGIN
        t_raster_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.tarr_raster_main[x_start,$
          (WHERE(hdr.tarr_raster_main[x_start,*] GT 0))[-1]]))+3,2)+'.1)'
        t_raster_real_txt = STRING(hdr.tarr_raster_main[x_start,t_start], $
          FORMAT=t_raster_real_format)
      ENDIF ELSE BEGIN
        t_raster_real_format = ''
        t_raster_real_txt = 'N/A'
      ENDELSE
      IF (raster_time_fb OR refraster_time_fb) THEN $
        t_raster_real_val = WIDGET_LABEL(params_main_base, VALUE=t_raster_real_txt, $
          /ALIGN_RIGHT) $
      ELSE $
        t_raster_real_val = 0
      ref_label   = WIDGET_LABEL(params_ref_base, VALUE=' ', /ALIGN_RIGHT)
      IF (hdr.refnt GT 1) THEN BEGIN
        t_ref_idx_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.refnt))+1,2)+')'
        t_ref_idx_txt = STRING(LONG(t_start), FORMAT=t_ref_idx_format)
      ENDIF ELSE BEGIN
        t_ref_idx_format = ''
        t_ref_idx_txt = 'N/A'
;		    t_ref_idx_val = 0
      ENDELSE
		  t_ref_idx_val = WIDGET_LABEL(params_ref_base, VALUE=t_ref_idx_txt, /ALIGN_RIGHT)
      IF ((hdr.refnt GT 1) AND dt_set) THEN BEGIN
        t_ref_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.tarr_ref[$
          (WHERE(hdr.tarr_ref GT 0))[-1]]))+3,2)+'.1)'
        t_ref_real_txt = STRING(hdr.tarr_ref[t_start], FORMAT=t_ref_real_format)
      ENDIF ELSE BEGIN
        t_ref_real_format = ''
        t_ref_real_txt = 'N/A'
        t_ref_real_val = 0
      ENDELSE
      IF dt_set THEN $
        t_ref_real_val = WIDGET_LABEL(params_ref_base, VALUE=t_ref_real_txt, $
          /ALIGN_RIGHT)
      IF ((hdr.refnt GT 1) AND refraster_time_fb) THEN BEGIN
        t_raster_ref_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.tarr_raster_ref[$
          x_start,(WHERE(hdr.tarr_raster_ref[x_start,*] GT 0))[-1]]))+3,2)+'.1)'
        t_raster_ref_real_txt = STRING(hdr.tarr_raster_ref[x_start,t_start], $
          FORMAT=t_raster_ref_real_format)
      ENDIF ELSE BEGIN
        t_raster_ref_real_format = ''
        t_raster_ref_real_txt = 'N/A'
      ENDELSE
      IF (raster_time_fb OR refraster_time_fb) THEN $
        t_raster_ref_real_val = WIDGET_LABEL(params_ref_base, VALUE=t_raster_ref_real_txt, $
          /ALIGN_RIGHT) $
      ELSE $
        t_raster_ref_real_val = 0
      sji_label   = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
      IF (hdr.sjint GT 1) THEN BEGIN
        t_sji_idx_format = '(I'+STRTRIM(FLOOR(ALOG10(hdr.sjint))+1,2)+')'
        t_sji_idx_txt = STRING(LONG(t_start), FORMAT=t_sji_idx_format)
      ENDIF ELSE BEGIN
        t_sji_idx_format = ''
        t_sji_idx_txt = 'N/A'
      ENDELSE
      IF ((hdr.sjint GT 1) AND dt_set) THEN BEGIN
        t_sji_real_format = '(F'+STRTRIM(FLOOR(ALOG10(hdr.tarr_sji[$
          (WHERE(hdr.tarr_sji GT 0))[-1]]))+3,2)+'.1)'
        t_sji_real_txt = STRING(hdr.tarr_sji[hdr.tsel_sji[t_start]], $
          FORMAT=t_sji_real_format)
      ENDIF ELSE BEGIN
        t_sji_real_format = ''
        t_sji_real_val = 0
        t_sji_real_txt = 'N/A'
      ENDELSE
		  t_sji_idx_val = WIDGET_LABEL(params_sji_base, VALUE=t_sji_idx_txt, /ALIGN_RIGHT)
		  IF dt_set THEN $
        t_sji_real_val = WIDGET_LABEL(params_sji_base, VALUE=t_sji_real_txt, $
          /ALIGN_RIGHT)
      IF (raster_time_fb OR refraster_time_fb) THEN $
        t_sji_raster_real_val = WIDGET_LABEL(params_sji_base, VALUE='N/A', $
          /ALIGN_RIGHT)

  ; Data value parameters
      no_label    = WIDGET_LABEL(verlabel_base, VALUE='Data values', /ALIGN_LEFT)
      real_label  = WIDGET_LABEL(verlabel_base, VALUE='Value', /ALIGN_RIGHT)
      main_label  = WIDGET_LABEL(params_main_base, VALUE=' ', /ALIGN_RIGHT)
      dataval_real_txt = STRING(((*hdr.imdata)[lp_start])[x_start,y_start], $
        FORMAT='(E10.4)')
		  dataval_real_val = WIDGET_LABEL(params_main_base, VALUE=dataval_real_txt, $
        /ALIGN_RIGHT)
      ; dataval_unit_txt = WIDGET_LABEL(params_main_base, VALUE='['+hdr.bunit+']', /ALIGN_RIGHT)
      ref_label   = WIDGET_LABEL(params_ref_base, VALUE=' ', /ALIGN_RIGHT)
      IF hdr.showref THEN BEGIN
        dataval_ref_real_txt = STRING(((*hdr.refdata)[lp_ref_start])[x_start,y_start], $
          FORMAT='(E10.4)')
        ; dataval_ref_unit_txt = '['+hdr.refbunit+']'
      ENDIF ELSE BEGIN
        dataval_ref_real_format = ''
        dataval_ref_real_txt = 'N/A'
        ; dataval_ref_unit_txt = ' '
      ENDELSE
		  dataval_ref_real_val = WIDGET_LABEL(params_ref_base, VALUE=dataval_ref_real_txt, $
        /ALIGN_RIGHT)
      ;dataval_ref_unit_txt = WIDGET_LABEL(params_ref_base, VALUE=dataval_ref_unit_txt, $
      ;  /ALIGN_CENTER)
      sji_label   = WIDGET_LABEL(params_sji_base, VALUE=' ', /ALIGN_RIGHT)
      ;dataval_sji_unit_txt = ' '
      IF hdr.sjifile THEN BEGIN
        dataval_sji_real_txt = STRING(((*hdr.sjidata)[0])[x_start,y_start], $
          FORMAT='(E10.4)')
        ;IF (STRTRIM(hdr.sjibunit,2) NE '0') THEN dataval_sji_unit_txt = '['+hdr.sjibunit+']'
      ENDIF ELSE BEGIN
        dataval_sji_real_format = ''
        dataval_sji_real_txt = 'N/A'
      ENDELSE
      dataval_sji_real_val = WIDGET_LABEL(params_sji_base, VALUE=dataval_sji_real_txt, $
        /ALIGN_RIGHT)
      ;dataval_sji_unit_txt = WIDGET_LABEL(params_sji_base, VALUE=dataval_sji_unit_txt, $
      ;  /ALIGN_CENTER)

   param_base = WIDGET_BASE(control_panel, /COLUMN)
    ; Column 1 of parameters overview containing cursor x,y and zoomfactor
;		disp = WIDGET_BASE(param_base, /COLUMN, /FRAME)
    ; Zommfactor info
		zoom_base = WIDGET_BASE(param_base, /ROW, /FRAME)
		zoom_txt = WIDGET_LABEL(zoom_base, VALUE = 'Zoom:')
		zoom_val = WIDGET_LABEL(zoom_base, $
      VALUE = STRING(zoomfactors[0]*100.,FORMAT='(I4)')+'%', /DYNAMIC_RESIZE)
;      VALUE = STRTRIM(LONG(zoomfactors[0]),2), /DYNAMIC_RESIZE)

	bg = WIDGET_BASE(cpanel, EVENT_PRO = 'CRISPEX_PB_BG')
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(initialising control panel)', /WIDGET, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Initializing control panel... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- SETTING UP DATA POINTERS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(setting up data pointers)', /WIDGET, /OVER
	feedback_text = [feedback_text,'> Setting up data pointers... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	index = INTARR(hdr.nx,hdr.ny,2)
	xarr = INDGEN(hdr.nx)
	yarr = INDGEN(hdr.ny)
	FOR i=0,hdr.nx-1 DO BEGIN
		index[i,*,1] = yarr
	ENDFOR
	FOR j=0,hdr.ny-1 DO BEGIN
		index[*,j,0] = xarr
	ENDFOR
	indexmap= PTR_NEW(index, /NO_COPY)
	indices = PTR_NEW(INTARR(hdr.nx,hdr.ny,2))

	xyslice	= PTR_NEW(BYTARR(hdr.nx,hdr.ny))
	sel_xyslice	= PTR_NEW(BYTARR(hdr.nx,hdr.ny))
	dopslice= PTR_NEW(BYTARR(hdr.nx,hdr.ny))
	sel_dopslice= PTR_NEW(BYTARR(hdr.nx,hdr.ny))
	emptydopslice= PTR_NEW(BYTARR(hdr.nx,hdr.ny))
	maskslice= PTR_NEW(BYTARR(hdr.nx,hdr.ny))

  phiscan = PTR_NEW(MAKE_ARRAY(hdr.nx,hdr.ny,hdr.nlp, TYPE=hdr.imtype))
  sspscan = PTR_NEW(MAKE_ARRAY(hdr.nx,hdr.ny,hdr.nlp*hdr.ns, TYPE=hdr.imtype))
	phislice= PTR_NEW(BYTARR(hdr.nlp,nphi))
	IF ((hdr.spfile EQ 1) OR (hdr.single_cube[0] GE 1)) THEN BEGIN
		loopslab= PTR_NEW(0)
		loopslice = PTR_NEW(0)
		refloopslab= PTR_NEW(0)
		refloopslice = PTR_NEW(0)
		crossloc= PTR_NEW(INTARR(nphi))
		exact_loopslab= PTR_NEW(0)
		exact_loopslice = PTR_NEW(BYTARR(nphi,hdr.mainnt))
		exact_crossloc= PTR_NEW(INTARR(nphi))
		rest_loopslab = PTRARR(nphi,/ALLOCATE_HEAP)
		rest_loopslice = PTRARR(nphi,/ALLOCATE_HEAP)
		rest_crossloc = PTRARR(nphi,/ALLOCATE_HEAP)
		FOR i=0,nphi-1 DO BEGIN
			*rest_loopslab[i]= PTR_NEW(0)
			*rest_loopslice[i] = PTR_NEW(0)
			*rest_crossloc[i]= PTR_NEW(0)
		ENDFOR
		det_loopslab= PTR_NEW(0)
		det_loopslice = PTR_NEW(BYTARR(nphi,hdr.mainnt))
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
	
	IF hdr.showref THEN $
		sel_refslice= PTR_NEW(BYTARR(hdr.refnx,hdr.refny)) $
  ELSE $
    sel_refslice = 0
  IF hdr.sjifile THEN $
		sel_sjislice= PTR_NEW(BYTARR(hdr.sjinx,hdr.sjiny)) $
  ELSE $
    sel_sjislice = 0
	
  imwintitle = 'CRISPEX'+instance_label+': Main image'
;	CRISPEX_WINDOW, imwinx, imwiny, control_panel, imwintitle, imwin, xywid, DRAWID = xydrawid, $
;		DRAWBASE = drawbase, 0, 0, /SCROLL, XSCROLL=xpos_slider, YSCROLL=ypos_slider;XSCROLL=imwinx, YSCROLL=imwiny
    ;, RESIZING = 1, RES_HANDLER = 'CRISPEX_DISPLAYS_XYREF_RESIZE'
  ; Create location for main image including scroll bars
  main = WIDGET_BASE(cpanel,/COLUMN)
  draw_verslid_base = WIDGET_BASE(main,/ROW)
  draw_horslid_base = WIDGET_BASE(main,/ROW)
	xydraw = WIDGET_DRAW(draw_verslid_base, XSIZE=imwinx, YSIZE=imwiny, RETAIN = 2)
  ypos_slider = WIDGET_SLIDER(draw_verslid_base,VALUE=0,MIN=0,MAX=1,/SUPPRESS,/DRAG,$
                          EVENT_PRO='CRISPEX_SLIDER_YPOS',/VERTICAL, YSIZE=imwiny)
  xpos_slider = WIDGET_SLIDER(draw_horslid_base,VALUE=0,MIN=0,MAX=1,/SUPPRESS,/DRAG,$
                          EVENT_PRO='CRISPEX_SLIDER_XPOS', XSIZE=imwinx)
	WIDGET_CONTROL, cpanel, /REALIZE, TLB_GET_SIZE=cpanel_size
  ; Determine window offsets based on realised control panel size and position
  ; If reference cube present, check if it would fit next to main image
  refxoffset = cpanel_size[0]
  refyoffset = ydelta
  IF hdr.showref THEN BEGIN
    windows_xextent = cpanel_size[0]+ spwinx + imwinx + lswinx + xdelta*3
    IF (windows_xextent LE x_scr_size) THEN refxoffset += xdelta 
  ENDIF 
 
  spxoffset = refxoffset+(hdr.showref*imwinx)+xdelta
  spyoffset = lswiny + ydelta
  lsxoffset = spxoffset
  
  WIDGET_CONTROL, xydraw, GET_VALUE=xydrawid

	WIDGET_CONTROL, xydraw, EVENT_PRO = 'CRISPEX_CURSOR', /SENSITIVE, /DRAW_MOTION_EVENTS, $
    /TRACKING_EVENTS, /DRAW_BUTTON_EVENTS
	
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(setting up data pointers)', /WIDGET, /OVER, /DONE
	IF startupwin THEN BEGIN
		WSET, startupwid
		WSHOW, startupwid
		feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Setting up data pointers... done!']
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	ENDIF

;--------------------------------------------------------------------------------- DEFINE INFO POINTER
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(defining main pointer)', /WIDGET, /OVER
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
    master_time_ids:master_time_ids, time_offset_slider:time_offset_slid, $
		lower_lp_text:lower_lp_text, upper_lp_text:upper_lp_text, $	
		reset_lprange_but:reset_lprange_but, lp_slider:lp_slid, $				
		lp_speed_slider:lp_speed_slid, lp_blink_slider:lp_blink_slid, $	
		lp_blink_but:lp_blink_but, lp_ref_but:lp_ref_but, lp_ref_slider:lp_ref_slid, $
		x_slider:x_slid, y_slider:y_slid, lock_button:lockbut, unlock_button:unlockbut, $
		zoom_button_ids:zoom_button_ids, xpos_slider:xpos_slider, ypos_slider:ypos_slider, $			
    stokes_button_ids:stokes_button_ids, stokes_spbutton_ids:stokes_spbutton_ids, $
    specwin_buts:specwin_buts, refspecwin_buts:refspecwin_buts, $
    specwin_button_ids:specwin_button_ids, refspecwin_button_ids:refspecwin_button_ids, $
		detspect_label:detspect_label, scale_detspect_but:scale_detspect_but, $
		detspect_im_but:detspect_im_but, detspect_ref_but:detspect_ref_but, $
		ls_toggle_but:ls_toggle_but, subtract_but:subtract_but, $		
		lower_y_text:lower_y_text, upper_y_text:upper_y_text, $	
		sp_toggle_but:sp_toggle_but, refsp_toggle_but:refsp_toggle_but, int_toggle_but:int_toggle_but, $
		phis_toggle_but:phis_toggle_but, reference_but:reference_but, doppler_but:doppler_but, $						
    scaling_cbox:scaling_cbox, imagescale_cbox:imagescale_cbox, $
    histo_opt_txt:histo_opt_txt, $
    gamma_label:gamma_label, gamma_slider:gamma_slider, scalemin_slider:scalemin_slider, $
    scalemax_slider:scalemax_slider, $
    scaling_reset_button:scaling_reset_but, scaling_reset_all_but:scaling_reset_all_but, $
    diagscale_label:diagscale_label, ls_mult_cbox:ls_mult_cbox, ls_mult_txt:ls_mult_txt, $
		phi_slider:phi_slid, nphi_slider:nphi_slid, loop_feedb_but:loop_feedb_but, $	
		bwd_move_slit:bwd_move_slit, fwd_move_slit:fwd_move_slit, $
		loop_slit_but:loop_slit_but, rem_loop_pt_but:rem_loop_pt_but, loop_slice_but:loop_slice_but, $	
		overlay_but:overlay_but, loop_overlay_all:loop_overlay_al, loop_overlay_sav:loop_overlay_sav, $		
		linestyle_0:linestyle_0, linestyle_1:linestyle_1, linestyle_2:linestyle_2, $			
		measure_but:measure_but, apix_label:apix_label, apix_unit:apix_unit, apix_text:apix_text, $
		measure_asec_lab:measure_asec_lab, measure_asec_text:measure_asec_text, $
		measure_km_lab:measure_km_lab, measure_km_text:measure_km_text, $					
		mask_button_ids:mask_button_ids, masks_overlay_ct_cbox:masks_overlay_ct_cbox, $
		masks_overlay_col_slid:masks_overlay_col_slid, raster_button:raster_but, $
    dispwid:dispwid, clear_current_inst:clear_current_inst, $
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
;--------------------------------------------------------------------------------- HEADER CONTROLS 
	ctrlshdr = { $
		header_select:0, header_txt:0 $
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
    xycoord_val:xycoord_val, xycoord_real_val:xycoord_real_val,  $
    refxycoord_val:refxycoord_val, refxycoord_real_val:refxycoord_real_val,  $
    sjixycoord_val:sjixycoord_val, sjixycoord_real_val:sjixycoord_real_val,  $
    lp_idx_val:lp_idx_val, lp_real_val:lp_real_val, lp_vdop_val:lp_vdop_val, $
    lp_ref_idx_val:lp_ref_idx_val, lp_ref_real_val:lp_ref_real_val, $
    lp_ref_vdop_val:lp_ref_vdop_val, $
    t_idx_val:t_idx_val, t_real_val:t_real_val, $ 
    t_raster_real_val:t_raster_real_val, $
    t_ref_idx_val:t_ref_idx_val, t_ref_real_val:t_ref_real_val, $
    t_raster_ref_real_val:t_raster_ref_real_val, $
    t_sji_idx_val:t_sji_idx_val, t_sji_real_val:t_sji_real_val, $
    dataval_real_val:dataval_real_val, dataval_ref_real_val:dataval_ref_real_val, $
    dataval_sji_real_val:dataval_sji_real_val, $
    zoom_val:zoom_val $
	}
;--------------------------------------------------------------------------------- PREFERENCE BUTTONS
	ctrlspref = { $
		startup_autopl:0, startup_win:0, displays_bgcols:0, displays_plcols:0, displays_interp:0, $
		displays_phislice:0, displays_slices:0, histo_opt_txt:0, gamma_slid:0, gamma_label:0,  $	
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
		imrefdetspect:0, lp_ref_lock:lp_ref_lock, bwd_insensitive:1, fwd_insensitive:1 $
	}
;--------------------------------------------------------------------------------- CURSOR 
	curs = { $
		sx:sx_start, sy:sy_start, sxlock:sx_start, sylock:sx_start, $		
    sxsji:sxsji_start, sysji:sysji_start, $
		xlock:x_start, ylock:y_start, lockset:0 $				
	}
;--------------------------------------------------------------------------------- DATA 
	data = { $
    imagedata:hdr.imdata, xyslice:xyslice, refdata:hdr.refdata, refslice:hdr.refslice, $
    maskdata:hdr.maskdata, maskslice:maskslice, dopslice:dopslice, dopplerscan:dopplerscan, $
    spdata:hdr.spdata, sspscan:sspscan, refspdata:hdr.refspdata, $
    refscan:hdr.refscan, refsspscan:hdr.refsspscan, spslice:PTR_NEW(0), refspslice:PTR_NEW(0), $
		emptydopslice:emptydopslice, scan:hdr.scan, phiscan:phiscan, phislice:phislice, $				
    sjidata:hdr.sjidata, sjislice:hdr.sjislice, $
		indexmap:indexmap, indices:indices, ratio:ratio, $
		lunsp:hdr.lunsp, lunim:hdr.lunim, lunrefim:hdr.lunrefim, lunrefsp:hdr.lunrefsp, lunmask:hdr.lunmask $
	}
;--------------------------------------------------------------------------------- DATA PARAMETERS
	dataparams = { $
    ; Filenames
		imfilename:hdr.imfilename, spfilename:hdr.spfilename, refimfilename:hdr.refimfilename, $
    refspfilename:hdr.refspfilename, maskfilename:hdr.maskfilename, $	
    ; Headers
    hdrs:[PTR_NEW(hdr.hdrs_main),PTR_NEW(hdr.hdrs_ref),PTR_NEW(hdr.hdrs_sji)], $
    next:[N_ELEMENTS(hdr.hdrs_main),N_ELEMENTS(hdr.hdrs_ref),N_ELEMENTS(hdr.hdrs_sji)], $
    ; Spatial dimensions
		x:x_start, y:y_start, d_nx:hdr.nx, d_ny:hdr.ny, nx:hdr.nx, ny:hdr.ny, $							
    sjinx:hdr.sjinx, sjiny:hdr.sjiny, sjidx:hdr.sjidx, sjidy:hdr.sjidy, $
    xsji:xsji_start, ysji:ysji_start, tarr_sji:hdr.tarr_sji, $;tarr_main:hdr.tarr_main, tarr_ref:hdr.tarr_ref, $
    tarr_raster_main:hdr.tarr_raster_main, tarr_raster_ref:hdr.tarr_raster_ref, $
		lc:hdr.lc, lp:lp_start, lp_ref:lp_ref_start, lp_dop:lp_start, nlp:hdr.nlp, refnlp:hdr.refnlp,$	
    ns:hdr.ns, s:0L, $					
		lps:hdr.lps, ms:hdr.ms, spec:hdr.mainspec, $
		reflps:hdr.reflps, refms:hdr.refms, refspec:hdr.refspec, $
		nt:hdr.mainnt, mainnt:hdr.mainnt, refnt:hdr.refnt, masknt:hdr.masknt, sjint:hdr.sjint, $		
    default_toffset_main:hdr.toffset_main, default_toffset_ref:hdr.toffset_ref, $
    dx:hdr.dx, dy:hdr.dy, $
    bunit:[hdr.bunit,hdr.refbunit], lpunit:[hdr.lpunit,hdr.reflpunit], xunit:hdr.xunit,$
    yunit:hdr.yunit, tunit:hdr.tunit $
	}
;--------------------------------------------------------------------------------- DATA SWITCH
	dataswitch = { $
		onecube:hdr.onecube, reffile:hdr.showref, refspfile:hdr.refspfile, spfile:hdr.spfile, $
    maskfile:hdr.maskfile, sjifile:hdr.sjifile  $							
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
		t_first:t_first, t_last:t_last, t_range:hdr.mainnt, t_low:t_first, t_upp:t_last, $					
		x_first:x_first, x_last:x_last, y_first:y_first, y_last:y_last, $				
		lp_first:lp_first, lp_last:lp_last, lp_range:hdr.nlp, lp_low:lp_first, lp_upp:lp_last, $
		lp_ref_first:lp_ref_first, lp_ref_last:lp_ref_last, lp_ref_low:lp_ref_first, $
    lp_ref_upp:(hdr.refnlp-1), lp_ref_range:hdr.refnlp, $
		nlpreb:nlpreb, ntreb:ntreb, refnlpreb:nlpreb, refntreb:refntreb, refloopnlxreb:nlpreb, refloopntreb:ntreb, $
		loopnlxreb:nlpreb, loopntreb:loopntreb, restloopnlxreb:nlpreb, restloopntreb:loopntreb, $
		retrdetnlxreb:nlpreb, retrdetntreb:loopntreb, phisnlpreb:nlpreb, nphireb:nphireb, $					
		xi:hdr.xi, yi:hdr.yi, xo:hdr.xo, yo:hdr.yo, xi_ref:hdr.xi_ref, yi_ref:hdr.yi_ref, $
    xo_ref:hdr.xo_ref, yo_ref:hdr.yo_ref, phisxtri:PTR_NEW(0), phisytri:PTR_NEW(0), $
    phisxi:PTR_NEW(0), phisyi:PTR_NEW(0), phisxo:PTR_NEW(0), phisyo:PTR_NEW(0), $
    xyrastersji:hdr.xyrastersji, sjix0:hdr.sjix0, sjiy0:hdr.sjiy0, $
    xsji_first:0L, xsji_last:(hdr.sjinx-1), ysji_first:0L, ysji_last:(hdr.sjiny-1),  $
    rastercont:hdr.rastercont, $
		interpspslice:interpspslice, phislice_update:phislice_update, slices_imscale:slices_imscale, $
    tsel_main:PTR_NEW(hdr.tsel_main), tsel_ref:PTR_NEW(hdr.tsel_ref), $
    tsel_sji:PTR_NEW(hdr.tsel_sji), master_time:0, $
    tarr_main:PTR_NEW(hdr.tarr_main[hdr.tsel_main]), $
    tarr_ref:PTR_NEW(hdr.tarr_ref[hdr.tsel_ref]), $
    tarr_sji:PTR_NEW(hdr.tarr_sji[hdr.tsel_sji]), $
    t:t_start, t_main:hdr.tsel_main[0], t_ref:hdr.tsel_ref[0], t_sji:hdr.tsel_sji[0], $
    t_low_main:hdr.tarr_main[0], t_upp_main:hdr.tarr_main[hdr.mainnt-1], $
    t_low_ref:hdr.tarr_ref[0], t_upp_ref:hdr.tarr_ref[hdr.refnt-1], $
    toffset_main:hdr.toffset_main, toffset_ref:hdr.toffset_ref $
	}
;--------------------------------------------------------------------------------- DATA DISPLAY SWITCHES
	dispswitch = { $
		restricted_t_range:PTR_NEW(0), restricted_lp_range:PTR_NEW(0), restricted_lp_ref_range:PTR_NEW(0), $
		exts:exts_set, refexts:refexts_set, warpspslice:hdr.warpspslice, warprefspslice:hdr.warprefspslice, $
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
		sel_diagnostics:PTR_NEW(hdr.sel_diagnostics), lines_diagnostics:PTR_NEW(hdr.lines_diagnostics), $
	  linlab_diagnostics:['Solid', 'Dotted', 'Dashed', 'Dash Dot', 'Dash Dot Dot', 'Long Dashes'],$
	  colors_diagnostics:[0,200,135,120,100,90,230,40],$
	  collab_diagnostics:['Black', 'Red', 'Pink', 'Purple', 'Blue', 'Turquoise', 'Grey', 'Green'],$
    selcol_diagnostics:PTR_NEW(hdr.selcol_diagnostics), $
    disp_diagnostics:REPLICATE(1,hdr.ndiagnostics), ndisp_diagnostics:hdr.ndiagnostics, $
    disp_refdiagnostics:REPLICATE(1,hdr.nrefdiagnostics), ndisp_refdiagnostics:hdr.nrefdiagnostics, $
		diagnostics:hdr.diagnostics, ndiagnostics:hdr.ndiagnostics, lock_t:1, $ 
    diag_start:hdr.diag_start, diag_width:hdr.diag_width, $
    diag_starts:PTR_NEW(hdr.diag_start), diag_widths:PTR_NEW(hdr.diag_width), $
    wheredispdiag:PTR_NEW(LONARR(hdr.ndiagnostics)), $
    refdiagnostics:hdr.refdiagnostics, nrefdiagnostics:hdr.nrefdiagnostics, $
    refdiag_start:hdr.refdiag_start, refdiag_width:hdr.refdiag_width, $
    refdiag_starts:PTR_NEW(hdr.refdiag_start), refdiag_widths:PTR_NEW(hdr.refdiag_width), $
    wheredisprefdiag:PTR_NEW(LONARR(hdr.nrefdiagnostics)), $
    lp_diag_all:0, lp_ref_diag_all:0 $
	}
;--------------------------------------------------------------------------------- I/O PARAMS
	ioparams = { $
		hdr:hdr $
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
		arcsecpix:hdr.dx, spatial_measurement:0, np:0, $					
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
		det_overlay_all:1, loopslit:0, overlalways:1, looppath_feedback:1, mask:hdr.maskfile, $
    maskim:[hdr.maskfile,hdr.showref,0], sjiraster:1 $		
	}
;--------------------------------------------------------------------------------- PARAMETER WINDOW CONTROLS 
	paramparams = { $
		wav_h:wav_h, sp_h:sp_h, scale_cubes:scale_cubes_vals, $
    xcoord_format:xcoord_format, ycoord_format:ycoord_format, $
    xcoord_real_format:xcoord_real_format, ycoord_real_format:ycoord_real_format, $
    refxcoord_format:refxcoord_format, refycoord_format:refycoord_format, $
    refxcoord_real_format:refxcoord_real_format, $
    refycoord_real_format:refycoord_real_format, $
    sjixcoord_format:sjixcoord_format, sjiycoord_format:sjiycoord_format, $
    sjixcoord_real_format:sjixcoord_real_format, $
    sjiycoord_real_format:sjiycoord_real_format, $
    lp_idx_format:lp_idx_format, lp_real_format:lp_real_format, $
    lp_vdop_format:lp_vdop_format, $
    lp_ref_idx_format:lp_ref_idx_format, lp_ref_real_format:lp_ref_real_format, $
    lp_ref_vdop_format:lp_ref_vdop_format, $
    t_idx_format:t_idx_format, t_real_format:t_real_format, $
    t_raster_real_format:t_raster_real_format, $
    t_ref_idx_format:t_ref_idx_format, t_ref_real_format:t_ref_real_format, $
    t_raster_ref_real_format:t_raster_ref_real_format, $
    t_sji_idx_format:t_sji_idx_format, t_sji_real_format:t_sji_real_format $
	}
;--------------------------------------------------------------------------------- PARAM SWITCHES
	paramswitch = { $
    dt_set:dt_set, t_raster:raster_time_fb $ 
	}
;--------------------------------------------------------------------------------- PATHS AND DIRECTORIES
	paths = { $
		ipath:ipath, opath:opath, opath_write:opath_write, hostname:hostname, $
    dir_settings:dir_settings, dir_settings_write:dir_settings_write, $
		dir_aux:dir_aux, dir_tanat:dir_aux, tanat_repointed:0 $
	}
;--------------------------------------------------------------------------------- PLAYBACK
	pbparams = { $
		t_step:t_step, t_speed:t_speed, direction:direction, bg:bg, mode:'PAUSE', lmode:'LOOP', $				
		imrefmode:0, lp_blink:lp_start, lp_blink_init:lp_start, spmode:0, spdirection:1 $						
	}
;--------------------------------------------------------------------------------- PHI SLIT VARIABLES 
	phiparams = { $
		d_nphi_set:nphi, nphi:nphi, angle:angle, nphi_set:LONG(hdr.ny/3.)-1, sphi:0, $			
		x_pts:PTR_NEW(0.), y_pts:PTR_NEW(0.), nw_cur:LONG(hdr.ny/3.), curindex:0, maxindex:0 $				
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
    xtickinterval:0, xdoptickinterval:0, xreftickinterval:0, xrefdoptickinterval:0, $
		ls_low_y:ls_low_y, ls_upp_y:ls_upp_y, ls_yrange:ls_yrange, $		
		ls_low_y_ref:ls_low_y_ref, ls_upp_y_ref:ls_upp_y_ref, ls_yrange_ref:ls_yrange_ref, $
		int_low_y:int_low_y, int_upp_y:int_upp_y, int_low_t:t_first, int_upp_t:t_last, $
		dt:hdr.dt, v_dop:hdr.v_dop, v_dop_ref:hdr.v_dop_ref, $		
    diag_ratio:PTR_NEW(hdr.diag_width/FLOAT(hdr.diag_width)), $
    diag_range_sp:PTR_NEW((hdr.diag_width/FLOAT(TOTAL(hdr.diag_width)))*xplspw), $
    refdiag_ratio:PTR_NEW(hdr.refdiag_width/FLOAT(hdr.refdiag_width)), $
    refdiag_range_sp:PTR_NEW((hdr.refdiag_width/FLOAT(TOTAL(hdr.refdiag_width)))*xplspw), $
    diag_range_phis:PTR_NEW((hdr.diag_width/FLOAT(TOTAL(hdr.diag_width)))*phisxplspw), $
    phis_yrange:[-1.,1.] $
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
		heightset:heightset, refheightset:refheightset, multichannel:hdr.multichannel, scalestokes:scalestokes, $
		v_dop_set:hdr.v_dop_set[0], v_dop_set_ref:hdr.v_dop_set[1], subtract:0, ref_subtract:0 $						
	}
;--------------------------------------------------------------------------------- PLOT TITLES
	plottitles = { $
		spxtitle:hdr.spxtitle, spytitle:spytitle, spwintitle:spwintitle, $
		refspxtitle:hdr.refspxtitle, refspwintitle:refspwintitle, $
		lsytitle:lsytitle, reflsytitle:reflsytitle, $
		lswintitle:lswintitle, reflswintitle:reflswintitle, phiswintitle:phiswintitle $
	}
;--------------------------------------------------------------------------------- PREFERENCE SETTINGS
	prefs = { $
		autoplay:autoplay, startupwin:startupwin, defsaveid:defsaveid,  $
		defipath:defipath, prefipath:prefipath, defopath:defopath, prefopath:prefopath, $
		bgplotcol_old:bgplotcol, plotcol_old:plotcol, interpspslice_old:interpspslice, $
		slices_imscale_old:slices_imscale, histo_opt_val:histo_opt_val, gamma_val:gamma_val, $
		tmp_autoplay:autoplay, tmp_startupwin:startupwin, tmp_defsaveid:defsaveid, $
		tmp_bgplotcol:bgplotcol, tmp_plotcol:plotcol, tmp_defipath:defipath, tmp_prefipath:prefipath, $
		tmp_defopath:defopath, tmp_prefopath:prefopath, tmp_interpspslice:interpspslice, $
		tmp_phislice_update:phislice_update, tmp_slices_imscale:slices_imscale, $								
    tmp_histo_opt_val:histo_opt_val, tmp_gamma_val:gamma_val, $
		default_autoplay:default_autoplay, default_startupwin:default_startupwin, $
		default_bgplotcol:default_bgplotcol, default_plotcol:default_plotcol, $
		default_defipath:default_defipath, default_prefipath:default_prefipath, $
		default_defopath:default_defopath, default_prefopath:default_prefopath, $
		default_defsaveid:default_defsaveid, default_interpspslice:default_interpspslice, $
		default_phislice_update:default_phislice_update, default_slices_imscale:default_slices_imscale, $
    default_histo_opt_val:default_histo_opt_val, default_gamma_val:default_gamma_val, $
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
		imagescale:imagescale, imagemin:immin[hdr.lc], imagemax:immax[hdr.lc], $
;    scale_range:PTR_NEW(REPLICATE(0.,4)), scale_minimum:PTR_NEW(REPLICATE(0.,4)), $	
;		scale_max_val:PTR_NEW(REPLICATE(255.,4)), scale_min_val:PTR_NEW(REPLICATE(0.,4)), $				
;		rel_scale_max_val:PTR_NEW(REPLICATE(100.,4)), rel_scale_min_val:PTR_NEW(REPLICATE(0.,4)), $			
		refmin:refmin[hdr.reflc], refmax:refmax[hdr.reflc], imrefscaling:0, relative:relative_scaling, $	
		scalestokes_max:hdr.scalestokes_max, dopmin:dopplermin[hdr.lc], dopmax:dopplermax[hdr.lc], $
    sjimin:sjimin, sjimax:sjimax, t_current:t_start, $
    imagemin_curr:immin[hdr.lc], imagemax_curr:immax[hdr.lc], $
		refmin_curr:refmin[hdr.reflc], refmax_curr:refmax[hdr.reflc], $
		dopmin_curr:dopplermin[hdr.lc], dopmax_curr:dopplermax[hdr.lc], $
    sjimin_curr:sjimin, sjimax_curr:sjimax, $
    spslice_min:REPLICATE(immin[hdr.lc],2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    spslice_max:REPLICATE(immax[hdr.lc],2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    phislice_min:REPLICATE(immin[hdr.lc],2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    phislice_max:REPLICATE(immax[hdr.lc],2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    gamma:REPLICATE(gamma_val,2*hdr.ndiagnostics+hdr.nrefdiagnostics+1), $
;    contrast:REPLICATE(0,2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    minimum:REPLICATE(0,2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    maximum:REPLICATE(100,2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    idx:0, diagscale_label_vals:diagscale_label_vals, mult_diag:0, $
    histo_opt_val:REPLICATE(histo_opt_val,2*hdr.ndiagnostics+hdr.nrefdiagnostics+1),  $
    mult_val:[main_mult_val,ref_mult_val], $
    sel_xyslice:sel_xyslice, sel_refslice:sel_refslice, sel_dopslice:sel_dopslice, $
    sel_sjislice:sel_sjislice $
	}
;--------------------------------------------------------------------------------- STOKES PARAMS
	stokesparams = { $
		labels:hdr.stokes_labels, button_labels:stokes_button_labels, select_sp:hdr.stokes_select_sp, $	
		prev_select_sp:hdr.stokes_select_sp $
	}
;--------------------------------------------------------------------------------- VERSION INFO
	versioninfo = { $
		version_number:version_number, revision_number:revision_number $				
	}
;--------------------------------------------------------------------------------- WINDOW IDs
	winids = { $
;		root:control_panel, imtlb:imwin, imwid:xywid, $		
		root:cpanel, imwid:xydrawid, $		
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
    sjiwid:0, sjitlb:0, sjidrawid:0, sjidrawbase:0, $
		savetlb:0, detsavetlb:0, restoretlb:0, preftlb:0, $							
		estimatetlb:0, savewintlb:0, saveoptwintlb:0, restsestlb:0, paramtlb:0, $					
		feedbacktlb:0, abouttlb:0, errtlb:0, warntlb:0, restsesfeedbtlb:0, $
    shorttlb:0, headertlb:0, $
		imwintitle:imwintitle, spwintitle:'',lswintitle:'',refwintitle:'',refspwintitle:'',reflswintitle:'', $
		imrefwintitle:'',dopwintitle:'',phiswintitle:'',restloopwintitle:PTR_NEW(''),retrdetwintitle:'',$
		loopwintitle:'',refloopwintitle:'',intwintitle:'', sjiwintitle:'' $
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
    sjiwinx:sjiwinx, sjiwiny:sjiwiny, $
		xdelta:xdelta, ydelta:ydelta,  $
    spxoffset:spxoffset, spyoffset:spyoffset, $
    lsxoffset:lsxoffset, lsyoffset:0, $
    refxoffset:refxoffset, refyoffset:refyoffset, aboutxoffset:startup_xpos, $
    aboutyoffset:startup_ypos $
		}
;--------------------------------------------------------------------------------- WINDOW SWITCHES
	winswitch = { $
		showls:0, showsp:0, showrefsp:hdr.refspfile, showrefls:showrefls, showimref:0, $
		estimate_win:0, showphis:0, showref:hdr.showref, showparam:0, showint:0, $
		showloop:0, showrefloop:0, showrestloop:0, showretrdet:0, showdop:0, dispwids:0, $
    showsji:hdr.sjifile $
	}
;--------------------------------------------------------------------------------- ZOOMING 
	zooming = { $
		factor:zoomfactors[0], factorswitch:factorswitch, factors:zoomfactors, xpos:0., ypos:0. $;, $
;    handle_extreme:handle_extreme $								
	}
;--------------------------------------------------------------------------------- DEFINE INFO POINTER
	info = { $
		ctrlscp:PTR_NEW(ctrlscp, /NO_COPY), $
		ctrlsdet:PTR_NEW(ctrlsdet, /NO_COPY), $
		ctrlsfeedb:PTR_NEW(ctrlsfeedb, /NO_COPY), $
		ctrlshdr:PTR_NEW(ctrlshdr, /NO_COPY), $
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
		ioparams:PTR_NEW(ioparams, /NO_COPY), $
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
	WIDGET_CONTROL, cpanel, SET_UVALUE = info

;	WIDGET_CONTROL, imwin, SET_UVALUE = info

	pseudoevent = { WIDGET_BUTTON, id:cpanel, top:cpanel, handler:0L, select:1 }
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(defining main pointer)', /WIDGET, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Defining main pointer... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text

;--------------------------------------------------------------------------------- DETERMINE DISPLAY OF PLOTS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(determining display of plots)', /WIDGET, /OVER
	feedback_text = [feedback_text,'> Determining display of plots... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	IF ((*(*info).dataswitch).spfile EQ 1) THEN BEGIN
		WIDGET_CONTROL, sp_toggle_but, SET_BUTTON = 1
		(*(*info).winswitch).showsp = 1
	ENDIF ELSE BEGIN
		IF (hdr.single_cube[0] EQ 1) THEN BEGIN
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
			IF (hdr.single_cube[0] GE 1) THEN WIDGET_CONTROL, approxmenu, SENSITIVE = 0
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
;		IF ((onecube EQ 0) OR (hdr.nt EQ 1)) THEN BEGIN
		IF (hdr.mainnt EQ 1) THEN BEGIN
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
			WIDGET_CONTROL, save_as_png_all, SENSITIVE = 0
			WIDGET_CONTROL, save_as_mpeg, SENSITIVE = 0
			WIDGET_CONTROL, sh_fbwd_button, SENSITIVE = 0
			WIDGET_CONTROL, sh_backward_button, SENSITIVE = 0
			WIDGET_CONTROL, sh_pause_button, SENSITIVE = 0
			WIDGET_CONTROL, sh_forward_button, SENSITIVE = 0
			WIDGET_CONTROL, sh_ffwd_button, SENSITIVE = 0
			WIDGET_CONTROL, timeslicemenu, SENSITIVE = 0
      WIDGET_CONTROL, det_file_loop, SENSITIVE=0
			*(*(*info).data).scan = (*(*(*info).data).scan)[0]
			*(*(*info).data).phiscan = *(*(*info).data).scan
		ENDIF 
		WIDGET_CONTROL, loop_slit_but, SENSITIVE = exts_set
		WIDGET_CONTROL, loop_feedb_but, SENSITIVE = exts_set
;		WIDGET_CONTROL, save_as_menu, SENSITIVE = exts_set
		WIDGET_CONTROL, timeslicemenu, SENSITIVE = 0 ;exts_set
	ENDELSE
	IF (N_ELEMENTS(hdr.single_cube[0]) EQ 1) THEN BEGIN
		WIDGET_CONTROL, ls_toggle_but, SET_BUTTON = (hdr.single_cube[0] NE 1) 
		(*(*info).winswitch).showls = (hdr.single_cube[0] NE 1) 
	ENDIF ELSE BEGIN
		WIDGET_CONTROL, ls_toggle_but, SET_BUTTON = 1
		(*(*info).winswitch).showls = 1
	ENDELSE
  WIDGET_CONTROL, (*(*info).ctrlscp).lower_lp_text, SENSITIVE=((*(*info).intparams).ndiagnostics LE 1) 
  WIDGET_CONTROL, (*(*info).ctrlscp).upper_lp_text, SENSITIVE=((*(*info).intparams).ndiagnostics LE 1) 

	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(determining display of plots)', /WIDGET, /OVER, /DONE
	feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Determining display of plots... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
;--------------------------------------------------------------------------------- OPENING WINDOWS
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(realising widget)', /WIDGET, /OVER
	feedback_text = [feedback_text,'> Realising widget... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	CRISPEX_MASK_BUTTONS_SET, pseudoevent
  set_zoomfac = CRISPEX_BGROUP_ZOOMFAC_SET(pseudoevent, /NO_DRAW, SET_FACTOR=0)
  IF (hdr.mainnt GT 1) THEN CRISPEX_DISPRANGE_T_RANGE, pseudoevent, /NO_DRAW
;  set_stokes = CRISPEX_DISPLAYS_STOKES_SELECT_XY(pseudoevent, /NO_DRAW, SET_STOKES=0)
	IF (*(*info).winswitch).showsp THEN BEGIN
		(*(*info).winswitch).showsp = 0
    CRISPEX_DRAW_GET_SPECTRAL_AXES, pseudoevent, /MAIN
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
    CRISPEX_DRAW_GET_SPECTRAL_AXES, pseudoevent, /REFERENCE
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
		CRISPEX_DISPLAYS_PHIS_TOGGLE, pseudoevent, /NO_DRAW
		IF startupwin THEN WSHOW, startupwid
	ENDIF
  IF (*(*info).winswitch).showsji THEN BEGIN
    (*(*info).winswitch).showsji = 0
		CRISPEX_DISPLAYS_SJI_TOGGLE, pseudoevent, /NO_DRAW
		IF startupwin THEN WSHOW, startupwid
  ENDIF

	CRISPEX_FIND_CLSAV, pseudoevent
	IF ((*(*info).retrparams).clfilecount NE 0) THEN BEGIN
		WIDGET_CONTROL, (*(*info).ctrlscp).sel_saved_loop, SENSITIVE = 1
		WIDGET_CONTROL, (*(*info).ctrlscp).all_saved_loop, SENSITIVE = 1
	ENDIF

	CRISPEX_UPDATE_T, pseudoevent
;	CRISPEX_DISPLAYS_PARAM_OVERVIEW_TOGGLE, pseudoevent
	IF ((*(*info).dataswitch).spfile EQ 1) OR (*(*info).dataswitch).onecube THEN $
    CRISPEX_UPDATE_SLICES, pseudoevent, /NO_DRAW
	IF showrefls THEN BEGIN
		IF (ref_detspect_scale EQ 0) THEN BEGIN
			CRISPEX_DISPRANGE_LS_SCALE_REF, pseudoevent
			CRISPEX_DISPRANGE_LS_RANGE, pseudoevent, /NO_DRAW
		ENDIF
		CRISPEX_DISPLAYS_DETSPECT_IM_SELECT, pseudoevent
	ENDIF
	IF (detspect_scale EQ 0) THEN BEGIN
		CRISPEX_DISPRANGE_LS_SCALE_MAIN, pseudoevent
		CRISPEX_DISPRANGE_LS_RANGE, pseudoevent, /NO_DRAW
	ENDIF
  CRISPEX_SCALING_APPLY_SELECTED, pseudoevent
	CRISPEX_DRAW, pseudoevent
	
;	IF ((*(*info).winswitch).showsp OR (*(*info).winswitch).showphis OR (*(*info).winswitch).showrefsp) THEN spwset = 1 ELSE spwset = 0
	spwset = ((*(*info).winswitch).showsp OR (*(*info).winswitch).showphis OR $
            (*(*info).winswitch).showrefsp) 
	IF (*(*info).winswitch).showls THEN $
    lsoffset = (*(*info).winswitch).showls * ((*(*info).winsizes).lswiny + $
    (*(*info).winsizes).ydelta) + (*(*info).dataswitch).refspfile * (*(*info).winsizes).ydelta $
   ELSE $
		lsoffset = (*(*info).winswitch).showrefls * (reflswiny + (*(*info).winsizes).ydelta)
;	WIDGET_CONTROL, cpanel, /REALIZE;, $
;    TLB_SET_XOFFSET = (*(*info).winsizes).xywinx + (*(*info).winsizes).refxoffset + $
;    2*(*(*info).winsizes).xdelta + spwset * ((*(*info).winsizes).xdelta + (*(*info).winsizes).spwinx), $
;		TLB_SET_YOFFSET = lsoffset
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(realising widget)', /WIDGET, /OVER, /DONE
	IF startupwin THEN BEGIN
		WSET, startupwid
		feedback_text = [feedback_text[0:N_ELEMENTS(feedback_text)-2],'> Realising widget... done!']
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	ENDIF

;--------------------------------------------------------------------------------- START MANAGING
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(start mangaing)', /WIDGET, /OVER
	feedback_text = [feedback_text,'> Start managing... ']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	XMANAGER, 'CRISPEX', cpanel, /NO_BLOCK
;	XMANAGER, 'CRISPEX', imwin, /NO_BLOCK
	WSHOW, (*(*info).winids).imwid
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, $
                                        '(start mangaing)', /WIDGET, /OVER, /DONE
	feedback_text = [feedback_text[0]+'done!',feedback_text[1:N_ELEMENTS(feedback_text)-2],$
    '> Start managing... done!']
	IF startupwin THEN CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, feedback_text
	WAIT,0.1
	IF (*(*info).prefs).autoplay THEN CRISPEX_PB_FORWARD, pseudoevent
	IF resave_preferences THEN CRISPEX_PREFERENCES_SAVE_SETTINGS, pseudoevent, /RESAVE
	IF (TOTAL(verbosity[0:1]) GE 1) THEN CRISPEX_UPDATE_STARTUP_SETUP_FEEDBACK, 'Set-up done!'
	IF startupwin THEN BEGIN
		WSET, startupwid
		CRISPEX_UPDATE_STARTUP_FEEDBACK, startup_im, xout, yout, 'Set-up done!'
		WAIT,0.25
		WIDGET_CONTROL, startuptlb, /DESTROY
	ENDIF
	CRISPEX_VERBOSE_SET_BUTTONS, pseudoevent

  ; Issue last warning/error messages before finishing setup
  IF (extreme_aspect AND (hdr.nx EQ 1)) THEN $
    CRISPEX_WINDOW_OK, pseudoevent, 'Warning', $
      'Extreme aspect ratio detected with NX = 1. Blowing ',$
      'up x-dimension for better visualisation: note that ',$
      'the image pixel aspect ratio is now inaccurate.', $
			OK_EVENT='CRISPEX_CLOSE_EVENT_WINDOW', /BLOCK
END
