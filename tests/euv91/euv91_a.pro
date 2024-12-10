;pro euv91
;************************************************************************
;	This program creates the solar EUV irradiance for 18-1050 A based on
;	the multiple linear regression model for solar EUV flux. Reference:
;	Tobiska, W.K., Revised solar extreme ultraviolet flux model, JATP,
;	(accepted) 1991. This paper was presented as XXVIII COSPAR, 
;	The Hague, STP I.2-P1, 1990.
;	The model values produced in this program include:
;		E_mod_ph	= photon flux (39 wlg, # days modeled)
;		E_mod		= energy flux (39 wlg, # days modeled)
;		specfluxm	= photon flux (21 5nm bins, # days modeled)
;	Files required for model to run:
;		euv91coef.txt	= ASCII file of model coefficients
;		euv91index.dat	= ASCII file of proxy values
;	Files created by running this model:
;		[Note: if these files are represented as arrays, then the 
;		 rows are the values of a given wavelength interval for all
;		 dates (starting with the label of the wavelength) and the 
;		 columns are the values of a given date for all wavelengths.]
;
;		<filename>.dat	= binary file (assoc(unit,fltarr(moddays+1))
;				  22 records: record            item
;						0		dates
;						1-21		5nm ph flux
;		<filename>_F.dat= binary file (assoc(unit,fltarr(moddays+1))
;				  40 records: record            item
;						0		dates
;						1-39		photon flux
;						40		daily flux total
;		<filename>_E.dat= binary file (assoc(unit,fltarr(moddays+1))
;				  40 records: record            item
;						0		dates
;						1-39		energy flux
;						40		daily flux total
;
;	K. Tobiska              Baseline version: 7-31-90      Rev: 4-30-91
;************************************************************************

print,'Initializing variables and reading files'
moddays=1 & yyddd=0.				; number of modeled days
c=fltarr(12,52)					; array of model coefficients
indices=fltarr(10,3750)				; array of proxy values
dates=fltarr(7500)				; array of dates available
lya=dates & he1=dates & f10=dates & f81=dates	; more arrays for indices
ij=fltarr(52)+1.				; identity vector
specwave=findgen(21)*50+25&specwave(0)=35.	; spectrum wavelengths
filename=' '					; filename for saving data
wave1=[18.62,30.02,50.52,100.54,150.10,$
        200.02,256.32,284.15,251.10,303.31,303.78,303.31,368.07,$
	356.01,401.14,465.22,453.00,500.00,$
	554.37,584.33,554.37,609.76,629.73,609.76,650.30,$
	703.36,701.00,765.15,770.41,787.71,750.01,801.00,851.00,$
	901.00,977.02,951.00,1025.72,1031.91,1001.00]
wave2=[29.52,49.22,99.99,148.40,198.58,249.18,$
	256.32,284.15,299.50,303.31,303.78,349.85,368.07,$
	399.82,436.70,465.22,499.37,550.00,$
	554.37,584.33,599.60,609.76,629.73,644.10,700.00,$
	703.36,750.00,765.15,770.41,787.71,800.00,850.00,900.00,$
	950.00,977.02,1000.00,1025.72,1031.91,1050.00]
wave=(wave1+wave2)/2.				; average wavelength interval
temp1=fltarr(39)&temp2=temp1&temp=fltarr(6)	; temporary arrays

openr,1,'euv91coef.txt'				; get the coefficients
  line=' ' & for i=0,2 do readf,1,line & readf,1,c
openr,2,'euv91index2.dat'			; get the indices
  for i=0,1 do readf,2,line & readf,2,indices
close,1 & close,2
  dates(0)=transpose(indices(0,*)) & dates(3750)=transpose(indices(5,*))
  lya(0)=transpose(indices(1,*)) & lya(3750)=transpose(indices(6,*))
  he1(0)=transpose(indices(2,*)) & he1(3750)=transpose(indices(7,*))
  f10(0)=transpose(indices(3,*)) & f10(3750)=transpose(indices(8,*))
  f81(0)=transpose(indices(4,*)) & f81(3750)=transpose(indices(9,*))
read,' How many days will be modeled? ',moddays
  ik=fltarr(moddays)+1.				; identity vector
  E_mod_ph=fltarr(39,moddays)			; 39 wlgs (photon flux)
  specfluxm=fltarr(21,moddays)			; consolidate to 5 nm bins
getdate:
read,' What is the starting date (yyddd: 68172-88366)? ',yyddd
  ind1=where(dates eq yyddd)			; find the starting date
  if (ind1(0) eq -1) then goto,getdate		; date out of range, try again
date=dates(ind1(0):ind1(0)+moddays-1)		; dates to be modeled
i1=lya(ind1(0):ind1(0)+moddays-1)		; lyman alpha
i2=he1(ind1(0):ind1(0)+moddays-1)		; helium I 10830
i3=f10(ind1(0):ind1(0)+moddays-1)		; daily f10.7
i4=f81(ind1(0):ind1(0)+moddays-1)		; 81-day smoothed f10.7
  for i=0,moddays-1 do begin			; substitute for missing data
    if (lya(ind1(0)+i)eq 0)and(he1(ind1(0)+i)gt 0) then i1(i) = i2(i) ; Lya=HeI
    if (lya(ind1(0)+i)eq 0) then i1(i) = (8.7e8)*i3(i) + 1.9e11       ; Lya=F10
    if (he1(ind1(0)+i)eq 0) then i2(i) = i1(i)                        ; HeI=Lya
  end

w1=1+(transpose(c(08,*))#ik)*(ij#exp(-i1*1e-10)); He I scaling factor
w2=1+(transpose(c(09,*))#ik)*(ij#exp(-i2*1e-10)); Lyman alpha scaling factor
w3=1+(transpose(c(10,*))#ik)*(ij#exp(-i3*.7))	; F81 scaling factor
w4=1+(transpose(c(11,*))#ik)*(ij#exp(-i4*.7))	; F10.7 scaling factor
model=((transpose(c(3,*))#ik)+$			; create the model 
                    (transpose(c(4,*))#ik)*(ij#i1)+$
                    (transpose(c(5,*))#ik)*(ij#i2)+$
                    (transpose(c(6,*))#ik)*(ij#i3)+$
                    (transpose(c(7,*))#ik)*(ij#i4))     *w1*w2*w3*w4

print,'Making 5 nm bins and 39 wavelength groups'
E_mod_ph(0,0)=model(0,*)
E_mod_ph(1,0)=model(1,*)+model(2,*)
E_mod_ph(2,0)=model(3,*)+model(4,*)
E_mod_ph(3,0)=model(5,*)+model(6,*)
E_mod_ph(4,0)=model(7,*)+model(8,*)
E_mod_ph(5,0)=model(9,*)+model(10,*)
E_mod_ph(6,0)=model(11,*)
E_mod_ph(7,0)=model(12,*)
E_mod_ph(8,0)=model(13,*)+model(14,*)
E_mod_ph(9,0)=model(15,*)
E_mod_ph(10,0)=model(16,*)
E_mod_ph(11,0)=model(17,*)
E_mod_ph(12,0)=model(18,*)
E_mod_ph(13,0)=model(19,*)+model(20,*)
E_mod_ph(14,0)=model(21,*)+model(22,*)
E_mod_ph(15,0)=model(23,*)
E_mod_ph(16,0)=model(24,*)+model(25,*)
E_mod_ph(17,0)=model(26,*)+model(27,*)
E_mod_ph(18,0)=model(28,*)
E_mod_ph(19,0)=model(29,*)
E_mod_ph(20,0)=model(30,*)
E_mod_ph(21,0)=model(31,*)
E_mod_ph(22,0)=model(32,*)
E_mod_ph(23,0)=model(33,*)+model(34,*)
E_mod_ph(24,0)=model(35,*)+model(36,*)
E_mod_ph(25,0)=model(37,*)
E_mod_ph(26,0)=model(38,*)
E_mod_ph(27,0)=model(39,*)
E_mod_ph(28,0)=model(40,*)
E_mod_ph(29,0)=model(41,*)
E_mod_ph(30,0)=model(42,*)+model(43,*)
E_mod_ph(31,0)=model(44,*)
E_mod_ph(32,0)=model(45,*)
E_mod_ph(33,0)=model(46,*)
E_mod_ph(34,0)=model(47,*)
E_mod_ph(35,0)=model(48,*)
E_mod_ph(36,0)=model(49,*)
E_mod_ph(37,0)=model(50,*)
E_mod_ph(38,0)=model(51,*)
E_mod=(12400.*1.6022e-12)*(E_mod_ph*((1./wave)#ik))	; energy flux

specfluxm(0,0)=E_mod_ph(0,*)+E_mod_ph(1,*)
specfluxm(1,0)=E_mod_ph(2,*)
specfluxm(2,0)=E_mod_ph(3,*)
specfluxm(3,0)=E_mod_ph(4,*)
specfluxm(4,0)=E_mod_ph(5,*)
specfluxm(5,0)=E_mod_ph(6,*)+E_mod_ph(7,*)+E_mod_ph(8,*)
specfluxm(6,0)=E_mod_ph(9,*)+E_mod_ph(10,*)+E_mod_ph(11,*)
specfluxm(7,0)=E_mod_ph(12,*)+E_mod_ph(13,*)
specfluxm(8,0)=E_mod_ph(14,*)
specfluxm(9,0)=E_mod_ph(15,*)+E_mod_ph(16,*)
specfluxm(10,0)=E_mod_ph(17,*)
specfluxm(11,0)=E_mod_ph(18,*)+E_mod_ph(19,*)+E_mod_ph(20,*)
specfluxm(12,0)=E_mod_ph(21,*)+E_mod_ph(22,*)+E_mod_ph(23,*)
specfluxm(13,0)=E_mod_ph(24,*)
specfluxm(14,0)=E_mod_ph(25,*)+E_mod_ph(26,*)
specfluxm(15,0)=E_mod_ph(27,*)+E_mod_ph(28,*)+E_mod_ph(29,*)+E_mod_ph(30,*)
specfluxm(16,0)=E_mod_ph(31,*)
specfluxm(17,0)=E_mod_ph(32,*)
specfluxm(18,0)=E_mod_ph(33,*)
specfluxm(19,0)=E_mod_ph(34,*)+E_mod_ph(35,*)
specfluxm(20,0)=E_mod_ph(36,*)+E_mod_ph(37,*)+E_mod_ph(38,*)

read,' What is the name of the file for saving the values: ',filename
  print,'opening '+filename+'.dat'
  openw,18,filename+'.dat',4*moddays+1		; 5 nm binned photon flux
  e=assoc(18,fltarr(moddays+1))
  print,'opening '+filename+'_F.dat'
  openw,19,filename+'_F.dat',4*moddays+1	; photon flux
  ef=assoc(19,fltarr(moddays+1))
  print,'opening '+filename+'_E.dat'
  openw,20,filename+'_E.dat',4*moddays+1	; energy flux
  ee=assoc(20,fltarr(moddays+1))
  for ii=0,21 do e(ii)=fltarr(moddays+1)	; initialize data files
  for ii=0,40 do begin
    ef(ii)=fltarr(moddays+1)
    ee(ii)=fltarr(moddays+1)
  end
  for ii=1,21 do begin				; insert the wavelengths
    temp(0)=specwave(ii-1)
    e(ii)=temp
  end
  for ii=1,39 do begin				; insert the wavelengths
    temp1(0)=wave1(ii-1)&temp2(0)=wave2(ii-1)
    ef(ii)=temp1
    ee(ii)=temp2
  end

print,'Inserting modeled values into data files'
  for i=1,moddays do begin
    print,date(i-1)
    temp=e(0)&temp(i)=date(i-1)&e(0)=temp
    for ii=1,21 do begin			; 5 nm bins - photon flux
      temp=e(ii)&temp(i)=specfluxm(ii-1,i-1)&e(ii)=temp
    end
    temp=ef(0)&temp(i)=date(i-1)&ef(0)=temp
    for ii=1,39 do begin
      temp=ef(ii)&temp(i)=E_mod_ph(ii-1,i-1)&ef(ii)=temp
    end
    temp=ef(40)&temp(i)=total(E_mod_ph(*,i-1))&ef(40)=temp
    temp=ee(0)&temp(i)=date(i-1)&ee(0)=temp
    for ii=1,39 do begin 
      temp=ee(ii)&temp(i)=E_mod(ii-1,i-1)&ee(ii)=temp
    end
    temp=ee(40)&temp(i)=total(E_mod(*,i-1))&ee(40)=temp
  end

close,18,19,20
end
