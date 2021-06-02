#/usr/bin/python3
import os
import numpy as np
import rpy2.robjects as robjects
rpf=os.path.split(os.path.realpath(__file__))[0]

def _input_sza(sza):
    sza=np.array(sza,dtype=np.float)
    sza=np.deg2rad(sza)
    sza=robjects.FloatVector(sza)
    return sza

def _input_date(date):
    #day of year
    doy=(date-date.astype("datetime64[Y]")).astype(int)+1

    if doy>365:
        doy=365
    Dayth=robjects.FloatVector([doy])
    #month of year
    moy=(date.astype("datetime64[M]")-date.astype("datetime64[Y]")).astype(int)+1
    Month=robjects.FloatVector([moy])
    #year
    year=date.astype("datetime64[Y]").astype(int)+1970
    Year=robjects.FloatVector([year])
    return Dayth,Month,Year

def _input_floats(X,N):
    X=np.array(X)
    if X.size==1:
        X=np.ones(N)*X
    elif X.size!=N:
        raise ValueError("size of array must be 1 or same as sza.")
    else:
        X=X
    X=robjects.FloatVector(X)
    return X

def _ang_turb(AOD550,ang):
    """Calculate Angstroem turbidity. From optical depth and Angstroem exponent.
    """
    return AOD550 * 0.55**ang

def _AM(sza):
    """Calculate absolute airmass
    """
    sza=np.array(sza) # [rad]
    AM=1/(np.cos(sza)+0.50572*((6.07995+(np.pi/2-sza)*180/np.pi)**(-1.6364)))
    return AM

def _Eext(date,Esc=1361.1):
    """
    Calculate extraterrestrial irradiance Eext with the Spencer Equation
        + Spencer,J.: Fourier series representation of the position of the sun. Search 1971;2(5). 172-172.
    
    Parameters
    ----------
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    Esc: float [Wm-2], optional
        Solar constant, default: 1361.1 Wm-2
    
    Returns
    -------
    Eext: float [Wm-2]
        extraterrestrial irradiance 
    """
    Dayth,Month,Year = _input_date(date)
    doy=float(Dayth[0])
    D = 2.*np.pi*(doy-1)/365.
    Eext = 1.00011 
    Eext+= 0.034221*np.cos(D)
    Eext+= 0.001280*np.sin(D)
    Eext+= 0.000719*np.cos(2*D)
    Eext+= 0.000077*np.sin(2*D)
    Eext=Eext*Esc
    return Eext

def _TL2_R(sza,wv,AOD550,ang):
    """
    Linke Turbidity at 2 am  - Remund 2003
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
        
    Returns
    -------
    TL2_R: array_like len(sza)
        Linke Turbidity at 2 air mass -  Remund 2003 version
    """
    ## work on inputs
    N=len(sza)
    wv =    np.array(_input_floats(wv,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    ang =     np.array(_input_floats(ang,N))
    
    #calculate angstroem turbidity
    b=_ang_turb(AOD550,ang)
    #calculate Linke Turbidity
    TL2R =   (1.8494 + 0.2425*wv - 0.0203*wv**2)
    TL2R+= b*(15.427 + 0.3153*wv - 0.0254*wv**2)
    return TL2R 

def _TL2_D(sza,wv,AOD550,ang):
    """
    Linke Turbidity at 2 am - Dogniaux 1972
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
        
    Returns
    -------
    TL2_D: array_like len(sza)
        Linke Turbidity at 2 air mass -  Dogniaux 1972 version    
    """
    ## work on inputs
    N=len(sza)
    wv =    np.array(_input_floats(wv,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    ang =     np.array(_input_floats(ang,N))
    
    #calculate angstroem turbidity
    b=_ang_turb(AOD550,ang)
    #
    TL2D = b*(16. + 0.22*wv)
    TL2D+= 0.1 
    TL2D+= (175. - sza) / (39.5*np.exp(-1.*wv) + 47.1)
    TL2D = TL2D/0.8662
    return TL2D

def _TL2_I(sza,pressure,wv,AOD550):
    """
    Linke Turbidity at 2 am - Ineichen 2008
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
        
    Returns
    -------
    TL2_I: array_like len(sza)
        Linke Turbidity at 2 air mass - Ineichen 2008 version
    """
    ## work on inputs
    N=len(sza)
    press = np.array(_input_floats(pressure,N))
    wv =    np.array(_input_floats(wv,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    
    
    p0=1013.25
    TL2I = 3.91*AOD550*np.exp(0.689*p0/press)
    TL2I+= 0.376*np.log(wv)
    TL2I+= 2.+ 0.54*p0/press - 0.5*(p0/press)**2 + 0.16*(p0/press)**3
    
    return TL2I


def _TL2_Gu(sza,pressure,wv,u_O,u_Ns,u_Nt,AOD550,ang):
    """
    Linke Turbidity at 2 am - Gueymard 1997
    + DOI: 10.1175/1520-0450(1998)037<0414:TDFBIM>2.0.CO;2
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    u_O: float or array_like  len(sza), [atm cm]
        total amount of O3 - typical values [0.1,0.5]
    u_Ns: float or array_like  len(sza), [atm cm]
        total amount of stratospheric NO2  - typical values [0,0.0005]
    u_Nt: float or array_like  len(sza), [atm cm]
        total amount of tropospheric NO2  - typical values [0,0.02]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of aerosol optical depth
        
    Returns
    -------
    TL2_Gu: array_like len(sza)
        Linke Turbidity at 2 air mass -  Gueymard1997 version
    """
    ## work on inputs
    N=len(sza)
    press = np.array(_input_floats(pressure,N))
    wv =    np.array(_input_floats(wv,N))
    u_O =    np.array(_input_floats(u_O,N))
    u_Ns =    np.array(_input_floats(u_Ns,N))
    u_Nt =    np.array(_input_floats(u_Nt,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    ang =     np.array(_input_floats(ang,N))
    
    ## calculate  all the stuff
    p0=1013.25 # hPa
    q=1.-(press/p0)
    
    #calculate optical mass Eq.(8)
    #rayleigh
    mR=1./(np.cos(np.deg2rad(sza))+0.45665*(sza**0.07)*((96.4836-sza)**(-1.697)))
    #watervapour
    mW=1./(np.cos(np.deg2rad(sza))+0.031141*(sza**0.1)*((92.471-sza)**(-1.3814)))
    #aerosol
    mA=mW
    #N02
    mN=mW
    
    #calculate broadband dry clean atmosphere optical depth (Eq(9))
    a0 = 1.00000 - 0.98173*q
    a1 = 0.18164 - 0.24259*q + 0.0507390*q**2
    a2 = 0.18164 - 0.17005*q - 0.0084949*q**2
    
    b0 = -0.0080617 + 0.028303*u_O - 0.014055*u_O**2
    b1 = +0.0113180 - 0.041018*u_O + 0.023471*u_O**2
    b2 = -0.0044577 + 0.016728*u_O - 0.010910*u_O**2
    
    c0 = 0.0036916 + 0.047361*u_O + 0.0058324*u_O**2
    c1 = 0.0154710 + 0.061662*u_O - 0.0440220*u_O**2
    c2 = 0.0399040 - 0.038633*u_O + 0.0548990*u_O**2
    
    f1 = (a0+a1*mR) / (1.+a2*mR)
    f2 = b0 + b1*mR**0.25 + b2*np.log(mR)
    f3 = (0.19758 + 0.00088585*mR - 0.097557*mR**0.2) / (1. + 0.0044767*mR)
    f4 = (c0 + c1*mR**(-0.72))/(np.exp(1.+c2*mR))
    f5 = u_Ns*(2.8669 - 0.078633*np.log(mR)**2.36)
    
    tau_c=f1*(f2+f3)+f4+f5
    
    # calculate broadband water vapor optical depth Eq.(A1)
    M=(1.7135 + 0.10004*mW + 0.00053986*mW**2)/(1.7149 + 0.097294*mW + 0.002567*mW**2)
    
    a0 = (1.72800 - 2.14510*q) / (1. - 0.96212*q)
    a1 = (0.37042 + 0.64537*q) / (1. + 0.94528*q)
    a2 = (3.51450 - 0.12483*q) / (1. - 0.34018*q)
    
    b0 = (0.63889 - 0.81121*q) / (1. - 0.799880*q)
    b1 = (0.06836 + 0.49008*q) / (1. + 4.723400*q)
    b2 = (2.15670 + 1.45460*q) / (1. + 0.038808*q)
    
    c0 = (-0.185700 + 0.23871*q) / (1. - 0.841110*q)
    c1 = (-0.022344 - 0.19312*q) / (1. + 6.216900*q)
    c2 = ( 2.170900 + 1.64230*q) / (1. + 0.062545*q)
    
    d0 = 3.3704 + 6.8096*q
    d1 = (12.487 - 18.517*q - 0.4089*q**2) / (1. - 1.4104*q)
    d2 = (2.5024 - 0.56834*q - 1.4623*q**2)/ (1. - 1.0252*q)
    d3 = (-0.030833 - 1.172*q - 0.98878*q**2) / (1. + 31.546*q)
    
    g1 = (a0*wv + a1*wv**1.6)/(1.+a2*wv)
    g2 = (b0*wv + b1*wv**1.6)/(1.+b2*wv)
    g3 = (c0*wv + c1*wv**1.6)/(1.+c2*wv)
    g4 = (d0*wv + d1*wv**0.62)/(1. + d2*wv + d3*wv**2)
    
    tau_wv = M*(g1 + g2*M*mW + g3*(M*mW)**1.28)/(1. + g4*M*mW)
    
    # calculate broadband aerosol optical depth Eq.14
    
    a0=(1.6685 + 4.1257*wv + 0.018748*wv**2)/(1. + 2.336*wv)
    a1=(0.075379 + 0.066532*wv - 0.0042634*wv**2)/(1. + 1.9477*wv)
    a2=(0.12867 + 0.24264*wv - 0.0087874*wv**2)/(1. + 3.3566*wv)
    
    b0=(-0.032335 - 0.0060424*wv) / (1.+0.023563*wv)
    b1=(-0.38229 - 0.0009926*wv) / (1.+0.044137*wv**0.594)
    b2=(-0.0059467 + 0.0054054*wv) / (1.+0.91487*wv)
    b3=(0.21989 + 0.041897*wv) / (1.+0.35717*wv)
    
    n = (1.3211 + 2.2036*wv) / (1. + 1.9367*wv)
    
    s1 = (a0 +a1*mA) / (1. + a2*mA)
    s2 = (b0 + b1*mA + b2*mA**2) / (1. + b3*mA**n)
     
    b=_ang_turb(AOD550,ang)
    tau_A=b*(s1 + s2*b)
    
    #calculate broadband N02 optical depth Eq(11)
    tau_N= u_Nt*(2.8669 - 0.078633*(np.log(mN)**2.36))
    
    TL2Gu = 1. + (mA/mR)*(tau_wv + tau_A + tau_N)/tau_c 
    TL2Gu = TL2Gu/0.8662
    return TL2Gu
    
def _TL2_M(sza,wv,altitude,AOD550,ang):
    """
    Linke Turbidity at 2 am - Molineaux 1995
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    altitude: float, [m]
        location altitude above sea level
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
        
    Returns
    -------
    TL2_M: array_like len(sza)
        Linke Turbidity at 2 air mass -  Molineaux 1995 version
    """
    ## work on inputs
    N=len(sza)
    wv =    np.array(_input_floats(wv,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    ang =     np.array(_input_floats(ang,N))
    h=altitude
    b=_ang_turb(AOD550,ang)
    
    f=np.exp(-0.12*h/1000.)/(np.cos(np.deg2rad(sza)) + 0.50572*(96.07995 - sza)**(-1.6364))
    TL2M = 1.5 + 12.4*b + 0.5*wv**(1./3.)
    TL2M+= 4.*(b-0.1)*np.log(f)
    return TL2M
    
def _TL2_Gr(sza,wv,AOD550,ang):
    """
    Linke Turbidity at 2 am - Grenier 1994
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle
    wv: float or array_like len(sza), [atm.cm]
        total columnar (pricipitable) water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
        
    Returns
    -------
    TL2_Gr: array_like len(sza)
        Linke Turbidity at 2 air mass -  Grenier 1994 version
    """
    ## work on inputs
    N=len(sza)
    wv =    np.array(_input_floats(wv,N))
    AOD550 = np.array(_input_floats(AOD550,N))
    ang =     np.array(_input_floats(ang,N))
    b=_ang_turb(AOD550,ang)
    
    TL2Gr = b - (-0.10545 - 0.02005*wv + 0.0050689*wv**2 - 0.0005202*wv**3)
    TL2Gr = TL2Gr / (0.073554 - 0.0029011*wv + 0.00075553*wv**2 - 0.000078281*wv**3) 
    
    return TL2Gr
    
def model_01_TJ(sza,date):
    """
    ==============
    ### 01-TJ 1957
    ==============
    Clear-sky irradiance model by Threkeld and Jordan 1957.
    
    References
    ----------
        +Threlkeld, J. L., & Jordan, R. C. (1957). 
            Direct solar radiation available on clear days. Heat.,
            Piping Air Cond., 29(12).
        +Masters, G. M. (2013).
            Renewable and efficient electric power systems. 
            John Wiley & Sons.
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)

    #
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/1-TJ.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)

    out=np.array(list(r.IrradianceTJ()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_02_schulze(sza):
    """
    ==============
    ### 02-Schulze 1975
    ==============
    Clear-sky irradiance model by Schulze 1975.
    
    References
    ----------
        +Schulze, R. E. (1976).
            A physically based method of estimating solar radiation from suncards.
            Agricultural Meteorology, 16(1), 85-101.
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 

    Returns
    -------
    DNI: array.shape(sza.shape,date.shape), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.shape,date.shape), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.shape,date.shape), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)

    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/2-Schulze.R"))
    r.assign("sza",sza)
    
    out=np.array(list(r.IrradianceSchulze()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI
    

def model_03_dpp(sza):
    """
    ==============
    ### 03-DPP 1978
    ==============
    Clear-sky irradiance model by Schulze 1975.
    
    References
    ----------
        +Daneshyar, M. (1978).
            Solar radiation statistics for Iran. 
            Sol. Energy;(United States), 21(4).
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). 
            Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. 
            Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 

    Returns
    -------
    DNI: array.shape(sza.shape,date.shape), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.shape,date.shape), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.shape,date.shape), [Wm-2]
        Global horizontal irradiance
    """
        #Inputs
    sza              = _input_sza(sza)
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/3-DPP.R"))
    r.assign("sza",sza)
    
    out=np.array(list(r.IrradianceDPP()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_04_adnot(sza):
    """
    ==============
    ### 04-Adnot 1979
    ==============
    Clear-sky irradiance model by Adnot et al 1979.
    
    References
    ----------
        +Adnot, J., Bourges, B., Campana, D., & Gicquel, R. (1979). 
            Utilisation de courbes de fréquences cumulées d'irradiation solaire globale pour le calcul des installations solaires.
        +Barbieri, F., Rifflart, C., Vo, B. T., Rajakaruna, S., & Ghosh, A. (2016, July). 
            A comparative study of clear-sky irradiance models for Western Australia. 
            In 2016 IEEE Power and Energy Society General Meeting (PESGM) (pp. 1-5). IEEE.
            
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 

    Returns
    -------
    DNI: np.nan
    DHI: np.nan
    GHI: array.shape(sza.shape,date.shape), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/4-Adnot.R"))
    r.assign("sza",sza)
    
    out=np.array(list(r.IrradianceAdnot()))
    GHI=out
    
    GHI=np.array(GHI)
    DNI=np.zeros(GHI.shape)*np.nan
    DHI=np.zeros(GHI.shape)*np.nan
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_05_biga(sza):
    """
    ==============
    ### 05-Biga 1976
    ==============
    Clear-sky irradiance model by Biga 1979.
    
    References
    ----------
        +Schulze, R. E. (1976). 
            A physically based method of estimating solar radiation from suncards. 
            Agricultural Meteorology, 16(1), 85-101.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 

    Returns
    -------
    DNI: array.shape(sza.shape,date.shape), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.shape,date.shape), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.shape,date.shape), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/5-Biga.R"))
    r.assign("sza",sza)
    
    out=np.array(list(r.IrradianceBiga()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_06_ashrae(sza,date):
    """
    ==============
    ### 06-ashrae 1985
    ==============
    Clear-sky irradiance model by Handbook 1985.
    
    References
    ----------
        +Handbook, A. F. (1985). American society of heating, refrigerating and air-conditioning engineers. Inc.:
            Atlanta, GA, USA.
        +El Mghouchi, Y., Ajzoul, T., Taoukil, D., & El Bouardi, A. (2016). 
            The most suitable prediction model of the solar intensity, on horizontal plane, at various weather conditions in a specified location in Morocco. 
            Renewable and Sustainable Energy Reviews, 54, 84-98.
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the month of the year from.

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)

    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/6-ASHRAE.R"))
    r.assign("sza",sza)
    r.assign("Month",Month)

    out=np.array(list(r.IrradianceASHRAE()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_07_sharma(sza,date,Esc=1366.1):
    """
    ==============
    ### 07-sharma 1965
    ==============
    Clear-sky irradiance model by Sharma 1965.
    
    References
    ----------
        +Sharma, M. R., & Pal, R. S. (1965). 
            Interrelationships between total, direct, and diffuse solar radiation in the tropics. 
            Solar Energy, 9(4), 183-192.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). 
            Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. 
            Renewable Energy, 55, 85-103.
    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    Esc              = _input_floats(Esc,len(sza))

    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/7-Sharma.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceSharma()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_08_elmghouchi(sza,date,Esc=1367):
    """
    ==============
    ### 08-El Mghouchi 2014
    ==============
    Clear-sky irradiance model by El Mghouchi 2014.
    
    References
    ----------
        +El Mghouchi, Y., El Bouardi, A., Choulli, Z., & Ajzoul, T. (2014). 
            New model to estimate and evaluate the solar radiation. 
            International Journal of Sustainable Built Environment, 3(2), 225-234.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    Esc              = _input_floats(Esc,len(sza))
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/8-ElMghouchi.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceElMghouchi()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_09_yangwalsh(sza,date,Esc=1362):
    """
    ==============
    ### 09-YangWalsh
    ==============
    Clear-sky irradiance model by Yang and Walsh 2014.
    
    References
    ----------
        +Yang, D., Walsh, W. M., & Jirutitijaroen, P. (2014). 
            Estimation and applications of clear sky global horizontal irradiance at the equator. 
            Journal of Solar Energy Engineering, 136(3), 034505.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    Esc: float [Wm-2], optional
        Solar constant, default: 1362 Wm-2

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    Esc              = _input_floats(Esc,len(sza))
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/9-YangWalsh.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceYangWalsh()))
    GHI=out
    
    GHI=np.array(GHI)
    DNI=np.zeros(GHI.shape)*np.nan
    DHI=np.zeros(GHI.shape)*np.nan
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_10_HLJ(sza,date,altitude,Esc=1353):
    """
    ==============
    ### 10-HLJ
    ==============
    Clear-sky irradiance model by Hottel 1976.
    
    References
    ----------
        +Hottel, H. C. (1976). 
            A simple model for estimating the transmittance of direct solar radiation through clear atmospheres. 
            Solar energy, 18(2), 129-134.
        +Gueymard, C. A. (2012). 
            Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. 
            Solar Energy, 86(8), 2145-2169. 

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        locations altitude
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    Esc              = _input_floats(Esc,len(sza))
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/10-HLJ.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceHLJ(altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_11_kumar(sza,date,pressure,Esc=1353):
    """
    ==============
    ### 11-Kumar 1997
    ==============
    Clear-sky irradiance model by Kumar 1997.
    
    References
    ----------
        +Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    Esc              = _input_floats(Esc,len(sza))
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/11-Kumar.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("press",press)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceKumar()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_12_campbell(sza,date,pressure,Esc=1366.1):
    """
    ==============
    ### 12-Campbell
    ==============
    Clear-sky irradiance model by Campbell and Normann 1998.
    
    References
    ----------
        +Campbell, G. S., & Norman, J. M. (1998). An introduction to environmental biophysics. An introduction to environmental biophysics., (Ed. 2).

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    Esc              = _input_floats(Esc,len(sza))
    
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/12-Campbell.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("press",press)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceCampbell()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_13_furich(sza,date,altitude,Esc=1367):
    """
    ==============
    ### 13- Fu and Rich 1999 (Release v9.2 2008)
    ==============
    Clear-sky irradiance model by Fu and Rich 1999.
    
    References
    ----------
        +Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        locations altitude
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    Esc              = _input_floats(Esc,len(sza))

    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/13-FuRich.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)

    out=np.array(list(r.IrradianceFurich(altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_14_atwaterball_1(sza,date,pressure,wv,Esc=1353):
    """
    ==============
    ### 14-Atwater and Ball_1 1981
    ==============
    Clear-sky irradiance model by Atwater and Ball 1981
    
    References
    ----------
        +Atwater, M. A., & Ball, J. T. (1981). A surface solar radiation model for cloudy atmospheres. Monthly weather review, 109(4), 878-888.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: nan
    DHI: nan
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/14-AtwaterBall1.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("Esc",Esc)

    # run R code
    out=np.array(list(r.IrradianceAtwater1()))
    GHI=np.array(out)
    DNI=np.zeros(GHI.shape)*np.nan
    DHI=np.zeros(GHI.shape)*np.nan
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_15_kasm(sza,date,pressure,wv,Esc=1367):
    """
    ==============
    ### 15-Kasm 1983
    ==============
    Clear-sky irradiance model by Kasten 1984
    
    References
    ----------
        +Kasten, F. (1984). Parametrisierung der Globalstahlung durch Bedeckungsgrad und Trubungsfaktor. Annalen der Meteorologie Neue, 20, 49-50.
        +Badescu, V. (1997). Verification of some very simple clear and cloudy sky models to evaluate global solar irradiance. Solar Energy, 61(4), 251-264.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source((os.path.join(pf,"R/15-KASM.R")))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("Esc",Esc)

    # run R code
    res=r.IrradianceKASM()
    # print(np.array(list(res)))
    out=np.array(list(r.IrradianceKASM()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_16_capderou(sza,date,pressure,altitude,latitude,Esc=1367):
    """
    ==============
    ### 16-Capderou 1987
    ==============
    Clear-sky irradiance model by Capderou 1987
    
    References
    ----------
        +Capderou, M. (1987). Theoretical and experimental models solar atlas of Algeria (in French) Tome 1 and 2. Algeria: University Publications Office.
        +Marif, Y., Chiba, Y., Belhadj, M. M., Zerrouki, M., & Benhammou, M. (2018). A clear sky irradiation assessment using a modified Algerian solar atlas model in Adrar city. Energy Reports, 4, 84-90.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure 
    altitude: float, [m]
        location altitude
    latitude: float, [deg]
        location latitude [-90,90]
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/16-Capderou.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("press",press)
    r.assign("Esc",Esc)
    
    # run R code
    out=np.array(list(r.IrradianceCapderou(latitude,altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_17_22_kasten(sza,date,altitude,TL2,Esc = 1367.):
    """
    ==============
    ### 17-22 Kasten 1984
    ==============
    Clear-sky irradiance model by Kasten
    
    References
    ----------
        +Ineichen, P., & Perez, R. (2002). A new airmass independent formulation for the Linke turbidity coefficient. Solar Energy, 73(3), 151-157.

    
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        location altitude
    TL2: array_like len(sza)
        Linke Turbidity at 2 air mass 
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: nan
    DHI: nan
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    
    sza=np.array(_input_sza(sza))

    Eext=_Eext(date,Esc)
    AM=_AM(sza)
    fh1 = np.exp(-altitude/8000)
    fh2 = np.exp(-altitude/1250)

    GHI=0.84*Eext*np.cos(sza)*np.exp(-0.027*AM*(fh1+fh2*(TL2-1.)))
    GHI[GHI<0]=0
    DNI=np.zeros(GHI.shape)*np.nan
    DHI=np.zeros(GHI.shape)*np.nan
    
    return DNI,DHI,GHI

def model_23_28_Heliosat1(sza,date,pressure,TL2,Esc=1367.):
    """
    ==============
    ### 23-28 Heliosat-1 1996
    ==============
    Clear-sky irradiance model by Hammer
    
    References
    ----------
        +Hammer, A., Heinemann, D., Hoyer, C., Kuhlemann, R., Lorenz, E., Müller, R., & Beyer, H. G. (2003). Solar energy assessment using remote sensing technologies. Remote Sensing of Environment, 86(3), 423-432.
        +Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    TL2: array_like len(sza)
        Linke Turbidity at 2 air mass 
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    sza=np.array(_input_sza(sza))
    press=np.array(_input_floats(pressure,len(sza)))
    
    Eext=_Eext(date,Esc)
    AM=_AM(sza)
    ama=AM*press/1013.25
    
    #rayleigh optical depth
    theta=1./(6.62960+1.7513*ama-0.1202*ama**2+0.0065*ama**3-0.00013*ama**4)
    theta[ama>20]=1./(10.4 + 0.718*ama[ama>20])
    
    DNI=Eext*np.exp(-ama*theta*TL2*0.8662)
    DNI[DNI<0]=0
    DHI=Eext*(0.0065+(0.0646*TL2-0.045)*np.cos(sza)-(0.0327*TL2-0.014)*(np.cos(sza))**2)
    DHI[DHI<0]=0
    GHI=DNI*np.cos(sza)+DHI
    GHI[GHI<0]=0
    return DNI,DHI,GHI

def model_29_34_esra(sza,date,altitude,TL2,Esc=1367.):
    """
    ==============
    ### 29-34 ESRA 2000
    ==============
    Clear-sky irradiance model for Esra
    
    References
    ----------
        +Rigollier, C., Bauer, O., & Wald, L. (2000). On the clear sky model of the ESRA—European Solar Radiation Atlas—with respect to the Heliosat method. Solar energy, 68(1), 33-48.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        location altitude
    TL2: array_like len(sza)
        Linke Turbidity at 2 air mass 
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    sza=np.array(_input_sza(sza))

    Eext=_Eext(date,Esc)
    
    alpha=np.pi/2.-sza
    alphatrue=alpha+0.061359*(0.1594+1.123*alpha+0.065656*alpha**2)/(1+28.9344*alpha+277.3971*alpha**2)   
    ame=np.exp(-altitude/8434.5)/(np.sin(alphatrue)+0.50572*(alphatrue/np.pi*180+6.07995)**-1.6364)
    
    #rayleigh optical depth
    rot =1./ (6.6296+1.7513*ame-0.1202*(ame**2)+0.0065*ame**3-0.00013*ame**4)
    rot[ame>20]=1./(10.4 + 0.718*ame[ame>20])


    #direct beam irradiance
    DNI = Eext*np.exp(-0.8662*TL2*ame*rot)
   
    #the diffuse transmission function at zenith
    TRD = (-1.5843e-2)+(3.0543e-2)*TL2+(3.797e-4)*TL2**2
    #the diffuse angular function
    a0 = (2.6463e-1)-(6.1581e-2)*TL2+(3.1408e-3)*TL2**2
    a1 = 2.0402+(1.8945e-2)*TL2-(1.1161e-2)*TL2**2
    a2 = -1.3025+(3.9231e-2)*TL2+(8.5079e-3)*TL2**2
    a0[a0*TRD<2e-3] = (2e-3)/TRD[a0*TRD<2e-3]
    FD = a0+a1*np.sin(alpha)+a2*(np.sin(alpha))**2
    #diffuse horizontal irradiance
    DHI = Eext*TRD*FD    
    
    #global horizontal irradiance
    GHI= DNI*np.cos(sza)+DHI

    DNI[DNI<0]=0
    DHI[DHI<0]=0
    GHI[GHI<0]=0
    return DNI,DHI,GHI


def model_35_40_Heliosat2(sza,date,altitude,TL2,Esc=1367.13):
    """
    ==============
    ### 35-40 Heliosat-2 2002
    ==============
    Clear-sky irradiance model by Lefevre
    
    References
    ----------
        +Lefevre, M., Albuisson, M., & Wald, L. (2002). Joint report on interpolation scheme ‘meteosat’and database ‘climatology’i (meteosat). SoDa Deliverable D3-8 and D5-1-4. Internal document.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        location altitude
    TL2: array_like len(sza)
        Linke Turbidity at 2 air mass 
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    sza=np.array(_input_sza(sza))
    N=len(sza)
    
    Eext=_Eext(date,Esc)
    
    alpha=np.pi/2.-sza
    alphatrue=alpha+0.061359*(0.1594+1.123*alpha+0.065656*alpha**2)/(1+28.9344*alpha+277.3971*alpha**2)   
    am=1./(np.sin(alphatrue)+0.50572*(alphatrue/np.pi*180.+6.07995)**-1.6364)
    
    #corr_rot is the correction of the integral Rayleigh optical thickness due to the elevation of the site     
    corr_rot075=1.248174-0.011997*am+0.00037*am**2
    corr_rot05=1.68219-0.03059*am+0.00089*am**2
    Y=np.vstack((corr_rot05,corr_rot075))
    Y=np.vstack((Y,np.ones(N)))
    
    #piecewise linear interpolation
    X=np.array([[0.5 ,1.],
                [0.75,1.],
                [1.  ,1.]])

    lin=np.linalg.lstsq(X,Y,rcond=None)[0]
    m=lin[0,:]
    y0=lin[1,:]
    pchup0=np.exp(-altitude/8434.5)     #all sites satisfy p/p0>=0.5
    corr_rot=m*pchup0+y0
    
    
    #Rayleigh Optical Thickness (rot)    
    rot = (corr_rot**-1)*(6.625928+1.92969*am-0.170073*am**2+0.011517*am**3-0.000285*am**4)**-1
    rot[am>20]=(10.4 + 0.718*am[am>20]*pchup0)**-1
    
    pi=np.pi
    L00=np.zeros(sza.size)    #define the length of L00 to be the same as the length of alpha
    L00[alpha>pi/6.]=-1.7349e-2
    L00[(alpha>pi/12)*(alpha<=pi/6)]=-8.2193e-3
    L00[alpha<=pi/12]=-1.1656e-3
    L01=np.zeros(sza.size)    #define the length
    L01[alpha>pi/6]=-5.8985e-3
    L01[(alpha>pi/12)*(alpha<=pi/6)]=4.5643e-4
    L01[alpha<=pi/12]=1.8408e-4
    L02=np.zeros(sza.size)   #define the length
    L02[alpha>pi/6.]=6.8868e-4
    L02[(alpha>pi/12)*(alpha<=pi/6)]=6.7916e-5
    L02[alpha<=pi/12]=-4.8754e-7
       
    L10=np.zeros(sza.size)    #define the length
    L10[alpha>pi/6.]=1.0258
    L10[(alpha>pi/12)*(alpha<=pi/6)]=8.9233e-1
    L10[alpha<=pi/12]=7.4095e-1
    L11=np.zeros(sza.size)    #define the length
    L11[alpha>pi/6.]=-1.2196e-1
    L11[(alpha>pi/12)*(alpha<=pi/6)]=-1.9991e-1
    L11[alpha<=pi/12]=-2.2427e-1
    L12=np.zeros(sza.size)    #define the length
    L12[alpha>pi/6.]=1.9299e-3
    L12[(alpha>pi/12)*(alpha<=pi/6)]=9.9741e-3
    L12[alpha<=pi/12]=1.5314e-2
    
    L20=np.zeros(sza.size)    #define the length
    L20[alpha>pi/6.]=-7.2178e-3
    L20[(alpha>pi/12)*(alpha<=pi/6)]=2.5428e-1
    L20[alpha<=pi/12]=3.4959e-1
    L21=np.zeros(sza.size)    #define the length
    L21[alpha>pi/6.]=1.3086e-1
    L21[(alpha>pi/12)*(alpha<=pi/6)]=2.6140e-1
    L21[alpha<=pi/12]=7.2313e-1
    L22=np.zeros(sza.size)    #define the length
    L22[alpha>pi/6.]=-2.8405e-3
    L22[(alpha>pi/12)*(alpha<=pi/6)]=-1.7020e-2
    L22[alpha<=pi/12]=-1.2305e-1
    L23=np.zeros(sza.size)    #define the length
    L23[alpha>pi/6.]=0
    L23[(alpha>pi/12)*(alpha<=pi/6)]=0
    L23[alpha<=pi/12]=5.9194e-3

    
    #direct beam irradiance
    Trb=np.exp(-0.8662*TL2*np.exp(altitude/8434.5)*rot)
    c0=L00+L01*TL2*pchup0+L02*(TL2*pchup0)**2
    c1=L10+L11*TL2*pchup0+L12*(TL2*pchup0)**2
    c2=L20+L21*TL2*pchup0+L22*(TL2*pchup0)**2+L23*(TL2*pchup0)**3
    Fb=c0+c1*np.sin(alpha)+c2*(np.sin(alpha))**2
    DNI=Eext*Trb*Fb/np.sin(alpha)
    DNI[DNI<0]=0
    
    #diffuse horizontal irradiance
    TRD = (-1.5843e-2)+(3.0543e-2)*TL2+(3.797e-4)*TL2**2
    a01 = (2.6463e-1)-(6.1581e-2)*TL2+(3.1408e-3)*TL2**2
    a11 = 2.0402+(1.8945e-2)*TL2-(1.1161e-2)*TL2**2
    a21 = -1.3025+(3.9231e-2)*TL2+(8.5079e-3)*TL2**2
    xiabiao=a01*TRD<(2e-3)
    a01[xiabiao] = (2e-3)/TRD[xiabiao]
    FD = a01+a11*np.sin(alpha)+a21*(np.sin(alpha))**2
    DHI = Eext*TRD*FD
    DHI[DHI<0]=0
    
    #global horizontal irradiance
    GHI=DNI*np.cos(sza)+DHI
    GHI[GHI<0]=0 
    
    return DNI,DHI,GHI

def model_41_46_ineichen_perez(sza,date,altitude,TL2,Esc = 1367.):
    """
    ==============
    ### 41-46 Ineichen & Perez 2002
    ==============
    Clear-sky irradiance model by Ineichen and Perez 2002
    
    References
    ----------
        +Ineichen, P., & Perez, R. (2002). A new airmass independent formulation for the Linke turbidity coefficient. Solar Energy, 73(3), 151-157.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float, [m]
        location altitude
    TL2: array_like len(sza)
        Linke Turbidity at 2 air mass 
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: nan
    DHI: nan
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    
    sza=np.array(_input_sza(sza))
    
    Eext=_Eext(date,Esc)
    AM=_AM(sza)
    
    fh1 = np.exp(-altitude/8000)
    fh2 = np.exp(-altitude/1250)
    cg1 = (5.09e-5)*(altitude)+0.868
    cg2 = (3.92e-5)*(altitude)+0.0387
    
    GHI= cg1*Eext*np.cos(sza)*np.exp(-cg2*AM*(fh1+fh2*(TL2-1)))
    GHI[GHI<0]=0
    DNI=np.zeros(GHI.size)*np.nan
    DHI=np.zeros(GHI.size)*np.nan
    
    return DNI,DHI,GHI


def model_47_CLS(sza,date,pressure,wv,albedo,Esc=1353):
    """
    ==============
    ### 47-CLS 1976
    ==============
    Clear-sky irradiance model by Suckling and Hay 1976
    
    References
    ----------
        +Suckling, P. W., & Hay, J. E. (1976). Modelling direct, diffuse, and total solar radiation for cloudless days. Atmosphere, 14(4), 298-308.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/47-CLS.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)

    out=np.array(list(r.IrradianceCLS()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI
    
def model_48_king(sza,date,pressure,wv,albedo,AOD550,ang,Esc=1366.1):
    """
    ==============
    ### 48-King 1979
    ==============
    Clear-sky irradiance model by King 1979
    
    References
    ----------
        +King, R., & Buckius, R. O. (1979). Direct solar transmittance for a clear sky. Solar Energy, 22(3), 297-301.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))
    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/48-King.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)

    out=np.array(list(r.IrradianceKing()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_49_josefsson(sza,date,pressure,wv,albedo,Esc=1366.1):
    """
    ==============
    ### 49-Josefsson 1985
    ==============
    Clear-sky irradiance model by Josefsson 1985
    
    References
    ----------
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/49-Josefsson.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)


    out=np.array(list(r.IrradianceJosefsson()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_50_badescu(sza,date,pressure,wv,u_O,Esc=1366.1):
    """
    ==============
    ### 50-Badescu 2008
    ==============
    Clear-sky irradiance model by Badescu 2008
    
    References
    ----------
        +Badescu, V. (2008). Use of sunshine number for solar irradiance time series generation. In Modeling Solar Radiation at the Earth’s Surface (pp. 327-355). Springer, Berlin, Heidelberg.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    u_O              = _input_floats(u_O,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/50-Badescu.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",u_O)


    out=np.array(list(r.IrradianceBadescu()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_51_simple_solis(sza,date,pressure,wv,AOD700,Esc=1367.):
    """
    ==============
    ### 51-Simplified Soils 2008
    ==============
    Clear-sky irradiance model for Solis 2008
    
    References
    ----------
        +Ineichen, P. (2008). A broadband simplified version of the Solis clear sky model. Solar Energy, 82(8), 758-762.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    AOD700: float or array_like len(sza)
        aerosol optical depth at 700nm
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    AOD700           = _input_floats(AOD700,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/51-Simplified_Solis.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("aod700",AOD700)


    out=np.array(list(r.IrradianceSolis()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_52_advanced_solis(sza,date,pressure,wv,AOD550,Esc=1367):
    """
    ==============
    ### 52-Advanced Soils 2018
    ==============
    Clear-sky irradiance model for Solis 2018
    
    References
    ----------
        +Ineichen, P. (2018). High turbidity solis clear sky model: development and validation. Remote Sensing, 10(3), 435.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/52-Advanced_Solis.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("aod550",AOD550)


    out=np.array(list(r.IrradianceSolis18()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_53_perrin(sza,date,pressure,wv,ozone,AOD550,ang,Esc=1367):
    """
    ==============
    ### 53-Perrin 1975
    ==============
    Clear-sky irradiance model by e Brichambaut 1975
    
    References
    ----------
        +e Brichambaut, C. P. (1975). Estimation des Ressources Energétiques en France. Cahiers de l’AFEDES, (1).
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like  len(sza), [atm cm]
        total amount of O3 - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        atmospherical optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent of optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))


    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/53-Perrin.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradiancePerrin()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_54_CEM(sza,date,pressure,wv,albedo,AOD0,Esc=1353.):
    """
    ==============
    ### 54-CEM 1978
    ==============
    Clear-sky irradiance model by Atwater and Ball 1978
    
    References
    ----------
        +Atwater, M. A., & Ball, J. T. (1978). A numerical solar radiation model based on standard meteorological observations. Solar Energy, 21(3), 163-170.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD0: float or array_like len(sza)
        broadband aerosol optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD0             = _input_floats(AOD0,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/54-CEM.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("broadbandaod",AOD0)


    out=np.array(list(r.IrradianceCEM()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_55_atwaterball2(sza,date,pressure,wv,albedo,AOD0,Esc=1353.):
    """
    ==============
    ### 55-Atwater and Ball_2 1981 
    ==============
    Clear-sky irradiance model by Atwater and Ball 1981
    
    References
    ----------
        +Bird, R., & Hulstrom, R. (1981). A simplified clear sky model for direct and diffuse insolation on horizontal surfaces. SERI. Solar Energy Research Institute, Golden, Colorado, 642-661.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD0: float or array_like len(sza)
        broadband aerosol optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD0             = _input_floats(AOD0,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/55-AtwaterBall2.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("broadbandaod",AOD0)


    out=np.array(list(r.IrradianceAtwater2()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_56_RSC(sza,date,pressure,wv,albedo,AOD550,ang,Esc=1371):
    """
    ==============
    ### 56-RSC 1985 
    ==============
    Clear-sky irradiance model by Carroll 1985
    
    References
    ----------
        +Carroll, J. J. (1985). Global transmissivity and diffuse fraction of solar radiation for clear and cloudy skies as measured and as predicted by bulk transmissivity models. Solar Energy, 35(2), 105-118.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1371 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))


    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/56-RSC.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceRSC()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_57_PSIM(sza,date,pressure,wv,albedo,AOD550,ang,Esc=1373):
    """
    ==============
    ### 57-PSIM 1993
    ==============
    Clear-sky irradiance model by Gueymard 1993
    
    References
    ----------
        +Gueymard, C. (1993). Mathermatically integrable parameterization of clear-sky beam and global irradiances and its use in daily irradiation applications. Solar Energy, 50(5), 385-397.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1373 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))


    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/57-PSIM.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradiancePSIM()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_58_bashahu(sza,date,pressure,wv,AOD550,ang,Esc=1367):
    """
    ==============
    ### 57-PSIM 1993
    ==============
    Clear-sky irradiance model by Gueymard 1993
    
    References
    ----------
        +Gueymard, C. (1993). Mathermatically integrable parameterization of clear-sky beam and global irradiances and its use in daily irradiation applications. Solar Energy, 50(5), 385-397.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/58-Bashahu.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceBashahu()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_59_MMAC(sza,date,pressure,wv,albedo,AOD0,Esc=1353.):
    """
    ==============
    ### 59-MMAC 1993,2003
    ==============
    Clear-sky irradiance model by Atwater and Ball 1981
    
    References
    ----------
        +Gueymard, C. (1993). Critical analysis and performance assessment of clear sky solar irradiance models using theoretical and measured data. Solar Energy, 51(2), 121-138.
        +Gueymard, C. A. (2003). Direct solar transmittance and irradiance predictions with broadband models. Part I: detailed theoretical performance assessment. Solar Energy, 74(5), 355-379.
        +Davies, J. A., & McKay, D. C. (1982). Estimating solar irradiance and components. Solar Energy, 29(1), 55-64.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD0: float or array_like len(sza)
        broadband aerosol optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD0             = _input_floats(AOD0,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/59-MMAC.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("broadbandaod",AOD0)


    out=np.array(list(r.IrradianceMMac()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_60_yang(sza,date,pressure,wv,ozone,AOD550,ang,Esc=1361.1):
    """
    ==============
    ### 60-Yang 2005
    ==============
    Clear-sky irradiance model by Yang 2005
    
    References
    ----------
        +Yang, K., & Koike, T. (2005). A general model to estimate hourly and daily solar radiation for hydrological studies. Water Resources Research, 41(10).

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1361.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/60-Yang.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceYang()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_61_calinoiu(sza,date,wv,ozone,no2,AOD550,ang,Esc=1366.1):
    """
    ==============
    ### 61-Calinoiu 2018
    ==============
    Clear-sky irradiance model by Calinoiu 2018
    
    References
    ----------
        +Calinoiu, D., Stefu, N., Boata, R., Blaga, R., Pop, N., Paulescu, E., ... & Paulescu, M. (2018). Parametric modeling: A simple and versatile route to solar irradiance. Energy Conversion and Management, 164, 175-187.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    no2: float or array_like len(sza), [atm.cm]
        total columnar NO2 - typical values [0.,0.02]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    no2              = _input_floats(no2,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/61-Calinoiu.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("NO2",no2)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceCalinoiu()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_62_Hoyt(sza,date,pressure,wv,ozone,albedo,AOD550,ang,Esc=1367):
    """
    ==============
    ### 62-Hoyt 1978
    ==============
    Clear-sky irradiance model by Hoyt 1987
    
    References
    ----------
        +Iqbal, M. (2012). An introduction to solar radiation. Elsevier.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))


    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/62-Hoyt.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceHoyt()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_63_mac2(sza,date,pressure,wv,albedo,AOD550,ang,Esc=1353):
    """
    ==============
    ### 63-MAC2 1982
    ==============
    Clear-sky irradiance model by Davies and McKay 1982
    
    References
    ----------
        +Davies, J. A., & McKay, D. C. (1982). Estimating solar irradiance and components. Solar Energy, 29(1), 55-64.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/58-Bashahu.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("albedo",albedo)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceBashahu()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_64_metstat(sza,date,pressure,wv,ozone,albedo,AOD0,Esc = 1367):
    """
    ==============
    ### 64-METSTAT 1998
    ==============
    Clear-sky irradiance model by Maxwell 1998
    
    References
    ----------
        +Maxwell, E. L. (1998). METSTAT—The solar radiation model used in the production of the National Solar Radiation Data Base (NSRDB). Solar Energy, 62(4), 263-279.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    albedo: float or array_like len(sza)
        surface albedo
    AOD0: float or array_like len(sza)
        broadband aerosol optical depth
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD0             = _input_floats(AOD0,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/64-METSTAT.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("albedo",albedo)
    r.assign("broadbandaod",AOD0)


    out=np.array(list(r.IrradianceMetstat()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_65_pr(sza,date,pressure,wv,ozone,albedo,AOD550,ang,Esc=1373):
    """
    ==============
    ### 65-PR 2000
    ==============
    Clear-sky irradiance model by Psiloglou 2000
    
    References
    ----------
        +Gueymard, C. A. (2003). Direct solar transmittance and irradiance predictions with broadband models. Part I: detailed theoretical performance assessment. Solar Energy, 74(5), 355-379.
        +Psiloglou, B. E., Santamouris, M., & Asimakopoulos, D. N. (2000). Atmospheric broadband model for computation of solar radiation at the Earth’s surface. Application to Mediterranean climate. pure and applied geophysics, 157(5), 829-860.
        +Badescu, V., Gueymard, C. A., Cheval, S., Oprea, C., Baciu, M., Dumitrescu, A., ... & Rada, C. (2013). Accuracy analysis for fifty-four clear-sky solar radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy, 55, 85-103.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1373 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/65-PR.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradiancePR()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_66_paulescu_schlett(sza,date,pressure,altitude,wv,ozone,AOD550,ang,Esc=1369):
    """
    ==============
    ### 66-Paulescu & Schlett 2003
    ==============
    Clear-sky irradiance model by Paulescu & Schlett 2003
    
    References
    ----------
        +Paulescu, M., & Schlett, Z. (2003). A simplified but accurate spectral solar irradiance model. Theoretical and applied climatology, 75(3-4), 203-212.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    altitude: float or array_like len(sza), [m]
        location altitude
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1369 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    altitude         = _input_floats(altitude,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/66-PaulescuSchlett.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradiancePaulescu(altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_67_mrm_v5(sza,date,pressure,wv,ozone,albedo,AOD550,ang,Esc=1366.1):
    """
    ==============
    ### 67-MRMv5 2008
    ==============
    Clear-sky irradiance model by Kambezidis 2008
    
    References
    ----------
        +Kambezidis, H. D., & Psiloglou, B. E. (2008). The meteorological radiation model (MRM): advancements and applications. In Modeling solar radiation at the earth’s surface (pp. 357-392). Springer, Berlin, Heidelberg.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/67-MRMv5.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("albedo",albedo)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceMRM5()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_68_mrm_v61(sza,date,pressure,wv,ozone,albedo,AOD550,aerosol_file=None,Esc=1360.8):
    """
    ==============
    ### 68-MRMv6.1 2017
    ==============
    Clear-sky irradiance model by Kambezidis 2017
    
    References
    ----------
        +Kambezidis, H. D., & Psiloglou, B. E. (2008). The meteorological radiation model (MRM): advancements and applications. In Modeling solar radiation at the earth’s surface (pp. 357-392). Springer, Berlin, Heidelberg.
        +Kambezidis, H. D., Psiloglou, B. E., Karagiannis, D., Dumka, U. C., & Kaskaoutis, D. G. (2016). Recent improvements of the Meteorological Radiation Model for solar irradiance estimates under all-sky conditions. Renewable Energy, 93, 142-158.
        +Kambezidis, H. D., Psiloglou, B. E., Karagiannis, D., Dumka, U. C., & Kaskaoutis, D. G. (2017). Meteorological Radiation Model (MRM v6. 1): Improvements in diffuse radiation estimates and a new approach for implementation of cloud products. Renewable and Sustainable Energy Reviews, 74, 616-637.
        +Jamie Bright's fortran code
        
    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    albedo: float or array_like len(sza)
        surface albedo
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    Esc: float [Wm-2], optional
        Solar constant, default: 1360.8 Wm-2    
    aerosol_file: string, optional
        path to aerosol-table provided by Jamie Brigth including AOD and SSA at different wavelength

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    press            = _input_floats(pressure,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    Esc              = _input_floats(Esc,len(sza))
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    
    pf = os.path.dirname(os.path.realpath(__file__))
    if aerosol_file==None:
        aerosol_file=os.path.join(pf,"R/AEROSOL-TABLE.txt")
    aerosol_table    = np.loadtxt(aerosol_file,delimiter=',')

    # R code import and assignment
    r=robjects.r
    
    r.source(os.path.join(pf,"R/68-MRMv61.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",press)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("albedo",albedo)
    r.assign("aod550",AOD550)
    r.assign("aerosol_table",aerosol_table)

    out=np.array(list(r.IrradianceMRM61()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    rpy2.robjects.numpy2ri.deactivate()
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_69_janjai(sza,date,altitude,wv,ozone,AOD550,ang,Esc=1366.1):
    """
    ==============
    ### 69-Janjai 2011
    ==============
    Clear-sky irradiance model by Janjai 2011
    
    References
    ----------
        +Janjai, S., Sricharoen, K., & Pattarapanitchai, S. (2011). Semi-empirical models for the estimation of clear sky solar global and direct normal irradiances in the tropics. Applied Energy, 88(12), 4749-4755.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float or array_like len(sza), [m]
        location altitude
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1366.1 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    altitude         = _input_floats(altitude,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/69-Janjai.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceJanjai(altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_70_bird(sza,date,pressure,albedo,wv,ozone,AOD550,ang,Esc=1353):
    """
    ==============
    ### 70-Bird 1981
    ==============
    Clear-sky irradiance model by Bird 1981
    
    References
    ----------
        +Bird, R. E., & Hulstrom, R. L. (1981). Simplified clear sky model for direct and diffuse insolation on horizontal surfaces (No. SERI/TR-642-761). Solar Energy Research Inst., Golden, CO (USA).

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    albedo: float or array_like len(sza)
        surface albedo
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1353 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    pressure         = _input_floats(pressure,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/70-Bird.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",pressure)
    r.assign("albedo",albedo)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceBird()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def model_71_iqbalC(sza,date,pressure,albedo,wv,ozone,AOD550,ang,Esc=1367):
    """
    ==============
    ### 71-Iqbal-C 1983
    ==============
    Clear-sky irradiance model by Iqbal 1983
    
    References
    ----------
        +Iqbal, M. (2012). An introduction to solar radiation. Elsevier.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    albedo: float or array_like len(sza)
        surface albedo
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    pressure         = _input_floats(pressure,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/71-IqbalC.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",pressure)
    r.assign("albedo",albedo)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceIqbalc()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_72_mod_iqbalC(sza,date,altitude,albedo,wv,ozone,AOD550,ang,Esc=1367):
    """
    ==============
    ### 72-Modified Iqbal-C 2012
    ==============
    Clear-sky irradiance model by Iqbal 1983 modified 2012
    
    References
    ----------
        +Bird, R. E., & Hulstrom, R. L. (1981). Simplified clear sky model for direct and diffuse insolation on horizontal surfaces (No. SERI/TR-642-761). Solar Energy Research Inst., Golden, CO (USA).
        +Gueymard, C. A. (2012). Clear-sky irradiance predictions for solar resource mapping and large-scale applications: Improved validation methodology and detailed performance analysis of 18 broadband radiative models. Solar Energy, 86(8), 2145-2169.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    altitude: float or array_like len(sza), [m]
        location altitude
    albedo: float or array_like len(sza)
        surface albedo
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    altitude         = _input_floats(altitude,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/72-Modified-IqbalC.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("albedo",albedo)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceMIqbalc(altitude)))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI

def model_73_rest2v5(sza,date,pressure,albedo,wv,ozone,no2,AOD550,ang,Esc=1367):
    """
    ==============
    ### 71-Iqbal-C 1983
    ==============
    Clear-sky irradiance model by Iqbal 1983
    
    References
    ----------
        +Iqbal, M. (2012). An introduction to solar radiation. Elsevier.

    Parameters
    ----------
    sza: array_like, [deg]
        solar zenith angle 
    date: np.datetime64
        date retrieve the day of the year from.
        day of the year is restricted to [1,365]
    pressure: float or array_like len(sza), [mb] [hPa]
        local barometric pressure
    albedo: float or array_like len(sza)
        surface albedo
    wv: float or array_like len(sza), [atm.cm]
        total columnar water vapour, typical values [0,6]
    ozone: float or array_like len(sza), [atm.cm]
        total columnar ozone - typical values [0.1,0.5]
    no2: float or array_like  len(sza), [atm cm]
        total columnar amount of NO2  - typical values [0,0.02]
    AOD550: float or array_like len(sza)
        aerosol optical depth at 550nm
    ang: float or array_like len(sza)
        angstroem exponent
    Esc: float [Wm-2], optional
        Solar constant, default: 1367 Wm-2    

    Returns
    -------
    DNI: array.shape(sza.size), [Wm-2]
        Direct normal irradiance
    DHI: array.shape(sza.size), [Wm-2]
        Diffuse horizontal irradiance
    GHI: array.shape(sza.size), [Wm-2]
        Global horizontal irradiance
    """
    #Inputs
    sza              = _input_sza(sza)
    Dayth,Month,Year = _input_date(date)
    pressure         = _input_floats(pressure,len(sza))
    albedo           = _input_floats(albedo,len(sza))
    wv               = _input_floats(wv,len(sza))
    ozone            = _input_floats(ozone,len(sza))
    no2              = _input_floats(no2,len(sza))
    AOD550           = _input_floats(AOD550,len(sza))
    ang              = _input_floats(ang,len(sza))
    Esc              = _input_floats(Esc,len(sza))

    ang_beta=_input_floats(_ang_turb(np.array(AOD550),np.array(ang)),len(sza))

    # R code import and assignment
    r=robjects.r
    pf = os.path.dirname(os.path.realpath(__file__))
    r.source(os.path.join(pf,"R/73-REST2v5.R"))
    r.assign("sza",sza)
    r.assign("Dayth",Dayth)
    r.assign("Year",Year)
    r.assign("Esc",Esc)
    r.assign("press",pressure)
    r.assign("albedo",albedo)
    r.assign("wv",wv)
    r.assign("ozone",ozone)
    r.assign("NO2",no2)
    r.assign("ang_alpha",ang)
    r.assign("ang_beta",ang_beta)


    out=np.array(list(r.IrradianceRest2v5()))
    DNI,DHI,GHI=out
    DNI=np.array(DNI)
    DHI=np.array(DHI)
    GHI=np.array(GHI)
    
    #remove r variables from memory
    r("rm(list=ls())")
    return DNI,DHI,GHI


def clear_sky_model_with_aod(aod,name,sza,date,pressure,altitude,albedo,wv,ozone,Esc):
    if name == 'MMAC':
        kwargs=dict(AOD0 = aod,
                    sza = sza,
                    date = date,
                    pressure = pressure,
                    wv = wv,
                    albedo = albedo)
        return model_59_MMAC(**kwargs)
    elif name == 'MRM61':
        kwargs=dict(sza = sza,
                    date = date,
                    pressure = pressure,
                    wv = wv,
                    ozone = ozone,
                    albedo = albedo,
                    AOD550 = aod)
        return model_68_mrm_v61(**kwargs)
    elif name == 'CEM':
        ### 54 CEM
        kwargs = dict(sza = sza,
                    date = date,
                    pressure = pressure,
                    wv = wv,
                    albedo = albedo,
                    AOD0 = aod)
        return model_54_CEM(**kwargs)
    elif name == 'KastenI':
        ### Kasten - IneichenTurbidity
        TL2I=_TL2_I(sza=sza,pressure=pressure,wv=wv,AOD550=aod)
        kwargs = dict(sza = sza,
                      date = date,
                      altitude = altitude,
                      TL2 = TL2I)
        return model_17_22_kasten(**kwargs)
    elif name == 'Heliosat1I':
        TL2I=_TL2_I(sza=sza,pressure=pressure,wv=wv,AOD550=aod)
        kwargs = dict(sza = sza,
                      date = date,
                      pressure = pressure,
                      TL2 = TL2I)
        return model_23_28_Heliosat1(**kwargs)
    elif name == 'ESRA':
        TL2I=_TL2_I(sza=sza,pressure=pressure,wv=wv,AOD550=aod)
        kwargs = dict(sza = sza,
                      date = date,
                      altitude = altitude,
                      TL2 = TL2I)
        return model_29_34_esra(**kwargs)
    elif name == 'Heliosat2I':
        TL2I=_TL2_I(sza=sza,pressure=pressure,wv=wv,AOD550=aod)
        kwargs = dict(sza = sza,
                      date = date,
                      altitude = altitude,
                      TL2 = TL2I,
                      Esc = Esc)
        return model_35_40_Heliosat2(**kwargs)
    elif name == 'IneichenPerez':
        TL2I=_TL2_I(sza=sza,pressure=pressure,wv=wv,AOD550=aod)
        kwargs = dict(sza = sza,
                      date = date,
                      altitude = altitude,
                      TL2 = TL2I)
        return model_41_46_ineichen_perez(**kwargs)    
    elif name == 'SOLISsimple':
        kwargs = dict(sza = sza,
                      date = date,
                      pressure = pressure,
                      wv = wv,
                      AOD700 = aod)
        return model_51_simple_solis(**kwargs)
    elif name == 'SOLISadvanced':
        kwargs = dict(sza = sza,
                      date = date,
                      pressure = pressure,
                      wv = wv,
                      AOD550 = aod,
                      Esc = Esc)
        return model_52_advanced_solis(**kwargs)
    elif name == 'AtwaterBall2':
        kwargs = dict(sza = sza,
                      date = date,
                      pressure = pressure,
                      wv = wv,
                      albedo = albedo,
                      AOD0 = aod)
        return model_55_atwaterball2(**kwargs)
    elif name == 'METSTAT':
        kwargs = dict(sza = sza,
                      date = date,
                      pressure = pressure,
                      wv = wv,
                      ozone = ozone,
                      albedo = albedo,
                      AOD0 = aod)
        return model_64_metstat(**kwargs)
    else:
        raise ValueError(str("Model name must be one of:"+
                             " MRM61, MMAC, METSTAT, AtwaterBall2, "+
                             "SOLISsimple, SOLISadvanced, IneichenPerez, "+
                             "Heliosat1I, Heliosat2I, ESRA, KastenI, CEM"))
    
    


    
    
    
