"""
This example runs all clear sky models with some standard atmospheric assumptions.
The parameters are all float values, which means they are assumend to be constant
over the day. But the input could also be an array of length of the sza input array.
Therefore, time is controlled via the solar zenith angle. 
"""
import numpy as np
import matplotlib.pyplot as plt

import clear_sky_models.models as csm

if __name__=='__main__':
    ### Inputs
    dates=np.datetime64("2014-04-01")
    S0=1367 # solar constant
    szas=np.arange(0,90,5) # solar zenith angle[deg]
    altitude=200 # altitude [m]
    pressure=1013.25 # local barometric pressure [hPa] [mb]
    latitude=50 # location latitude [-90,90]
    albedo=0.2 # surface albedo
    
    # gases
    wv=5 # total columnar amount of water vapour [atm.cm], typical values [0,6]
    u_O=0.35 # total amount of O3 [atm cm] - typical values [0.1,0.5]
    u_Nt=0.002 # total amount of tropospheric NO2 [atm cm]  - typical values [0,0.02]
    u_Ns=0.0001 # total amount of stratospheric NO2 [atm cm] - typical values [0,0.0005]
    u_N=u_Nt+u_Ns # total amount of columnar NO2 [atmc m] - typical values [0,0.0205]

    #aerosol
    ang=1.3 # angstroem exponent of atmospheric optical depth
    AOD550=0.25 # aerosol optical depth at 550nm
    AOD700=AOD550*(700./550.)**(-ang) # aerosol optical depth at 700nm
    AOD0=0.1 # broadband aerosol optical depth
    
    
    ### converted inputs
    Eext=csm._Eext(date=dates,Esc=S0)
    ang_turb=csm._ang_turb(AOD550,ang)
    
    TL2R=csm._TL2_R(sza=szas,wv=wv,AOD550=AOD550,ang=ang)
    TL2D=csm._TL2_D(sza=szas,wv=wv,AOD550=AOD550,ang=ang)
    TL2I=csm._TL2_I(sza=szas,pressure=pressure,wv=wv,AOD550=AOD550)
    TL2Gu=csm._TL2_Gu(sza=szas,pressure=pressure,wv=wv,u_O=u_O,u_Ns=u_Ns,u_Nt=u_Nt,AOD550=AOD550,ang=ang)
    TL2M=csm._TL2_M(sza=szas,wv=wv,altitude=altitude,AOD550=AOD550,ang=ang)
    TL2Gr=csm._TL2_Gr(sza=szas,wv=wv,AOD550=AOD550,ang=ang)

    ### solar irradiance models
    ### 01 TJ
    DNI01,DHI01,GHI01=csm.model_01_TJ(sza=szas,
                            date=dates)
    
    ### 02 Schulze
    DNI02,DHI02,GHI02=csm.model_02_schulze(sza=szas)
    
    ### 03 DPP
    DNI03,DHI03,GHI03=csm.model_03_dpp(sza=szas)
    
    ### 04 Adnot
    DNI04,DHI04,GHI04=csm.model_04_adnot(sza=szas)
    
    ### 05 Biga
    DNI05,DHI05,GHI05=csm.model_05_biga(sza=szas)
    
    ### 06 Ashrae
    DNI06,DHI06,GHI06=csm.model_06_ashrae(sza=szas,
                                date=dates)
    
    ### 07 Sharma
    DNI07,DHI07,GHI07=csm.model_07_sharma(sza=szas,
                                date=dates,
                                Esc=S0)
    
    ### 08 El Mghouchi
    DNI08,DHI08,GHI08=csm.model_08_elmghouchi(sza=szas,
                                    date=dates,
                                    Esc=S0)
    
    ### 09 Yangwalsh
    DNI09,DHI09,GHI09=csm.model_09_yangwalsh(sza=szas,
                                  date=dates,
                                  Esc=S0)
    
    ### 10 HLJ
    DNI10,DHI10,GHI10=csm.model_10_HLJ(sza=szas,
                            date=dates,
                            altitude=altitude,
                            Esc=S0)
    
    ### 11 Kumar
    DNI11,DHI11,GHI11=csm.model_11_kumar(sza=szas,
                              date=dates,
                              pressure=pressure,
                              Esc=S0)
    
    ### 12 Campbell
    DNI12,DHI12,GHI12=csm.model_12_campbell(sza=szas,
                                  date=dates,
                                  pressure=pressure,
                                  Esc=S0)
    
    ### 13 Furich
    DNI13,DHI13,GHI13=csm.model_13_furich(sza=szas,
                                date=dates,
                                altitude=altitude,
                                Esc=S0)
    
    ### 14 Atwaterball
    DNI14,DHI14,GHI14=csm.model_14_atwaterball_1(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    Esc=S0)
    
    ### 15 KASM
    DNI15,DHI15,GHI15=csm.model_15_kasm(sza=szas,
                             date=dates,
                             pressure=pressure,
                             wv=wv,
                             Esc=S0)
    
    ### 16 Capderou
    DNI16,DHI16,GHI16=csm.model_16_capderou(sza=szas,
                                 date=dates,
                                 pressure=pressure,
                                 altitude=altitude,
                                 latitude=latitude,
                                 Esc=S0)
    
    ### 17 - 22 Kasten
    DNI17,DHI17,GHI17=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2R,
                                         Esc=S0)

    DNI18,DHI18,GHI18=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2D,
                                         Esc=S0)
    DNI19,DHI19,GHI19=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2I,
                                         Esc=S0)
    DNI20,DHI20,GHI20=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2Gu,
                                         Esc=S0)
    DNI21,DHI21,GHI21=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2M,
                                         Esc=S0)
    DNI22,DHI22,GHI22=csm.model_17_22_kasten(sza=szas,
                                         date=dates,
                                         altitude=altitude,
                                         TL2=TL2Gr,
                                         Esc=S0)
    ### 23 - 28 Heliosat1
    DNI23,DHI23,GHI23=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2R,
                                             Esc=S0)
    DNI24,DHI24,GHI24=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2D,
                                             Esc=S0)
    DNI25,DHI25,GHI25=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2I,
                                             Esc=S0)
    DNI26,DHI26,GHI26=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2Gu,
                                             Esc=S0)
    DNI27,DHI27,GHI27=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2M,
                                             Esc=S0)
    DNI28,DHI28,GHI28=csm.model_23_28_Heliosat1(sza=szas,
                                             date=dates,
                                             pressure=pressure,
                                             TL2=TL2Gr,
                                             Esc=S0)
    
    ### 29 - 34 ESRA
    DNI29,DHI29,GHI29=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2R,
                                             Esc=S0)
    DNI30,DHI30,GHI30=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2D,
                                             Esc=S0)
    DNI31,DHI31,GHI31=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2I,
                                             Esc=S0)
    DNI32,DHI32,GHI32=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gu,
                                             Esc=S0)
    DNI33,DHI33,GHI33=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2M,
                                             Esc=S0)
    DNI34,DHI34,GHI34=csm.model_29_34_esra(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gr,
                                             Esc=S0)
    
    ### 35 - 40 Heliosat2
    DNI35,DHI35,GHI35=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2R,
                                             Esc=S0)
    DNI36,DHI36,GHI36=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2D,
                                             Esc=S0)
    DNI37,DHI37,GHI37=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2I,
                                             Esc=S0)
    DNI38,DHI38,GHI38=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gu,
                                             Esc=S0)
    DNI39,DHI39,GHI39=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2M,
                                             Esc=S0)
    DNI40,DHI40,GHI40=csm.model_35_40_Heliosat2(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gr,
                                             Esc=S0)
    
    ### 41 - 46 Ineichen Perez
    DNI41,DHI41,GHI41=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2R,
                                             Esc=S0)
    DNI42,DHI42,GHI42=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2D,
                                             Esc=S0)
    DNI43,DHI43,GHI43=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2I,
                                             Esc=S0)
    DNI44,DHI44,GHI44=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gu,
                                             Esc=S0)
    DNI45,DHI45,GHI45=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2M,
                                             Esc=S0)
    DNI46,DHI46,GHI46=csm.model_41_46_ineichen_perez(sza=szas,
                                             date=dates,
                                             altitude=altitude,
                                             TL2=TL2Gr,
                                             Esc=S0)
    
    ### 47 CLS
    DNI47,DHI47,GHI47=csm.model_47_CLS(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    Esc=S0)
    ### 48 King
    DNI48,DHI48,GHI48=csm.model_48_king(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)

    ### 49 Josefsson
    DNI49,DHI49,GHI49=csm.model_49_josefsson(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    Esc=S0)
    
    ### 50 Badescu
    DNI50,DHI50,GHI50=csm.model_50_badescu(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    u_O=u_O,
                                    Esc=S0)
    
    ### 51 simplified solis
    DNI51,DHI51,GHI51=csm.model_51_simple_solis(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    AOD700=AOD700,
                                    Esc=S0)
    
    ### 52 advanced solis
    DNI52,DHI52,GHI52=csm.model_52_advanced_solis(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    AOD550=AOD550,
                                    Esc=S0)
    
    ### 53 Perrin
    DNI53,DHI53,GHI53=csm.model_53_perrin(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 54 CEM
    DNI54,DHI54,GHI54=csm.model_54_CEM(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD0=AOD0,
                                    Esc=S0)
    
    ### 55 Atwater and Ball 2
    DNI55,DHI55,GHI55=csm.model_55_atwaterball2(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD0=AOD0,
                                    Esc=S0)
    
    ### 56 RSC
    DNI56,DHI56,GHI56=csm.model_56_RSC(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 57 PSIM
    DNI57,DHI57,GHI57=csm.model_57_PSIM(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 58 Bashahu
    DNI58,DHI58,GHI58=csm.model_58_bashahu(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 59 MMAC
    DNI59,DHI59,GHI59=csm.model_59_MMAC(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD0=AOD0,
                                    Esc=S0)
    
    ### 60 Yang
    DNI60,DHI60,GHI60=csm.model_60_yang(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 61 Calinoiu
    DNI61,DHI61,GHI61=csm.model_61_calinoiu(sza=szas,
                                    date=dates,                                    
                                    wv=wv,
                                    ozone=u_O,
                                    no2=u_N,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 62 Hoyt
    DNI62,DHI62,GHI62=csm.model_62_Hoyt(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 63 MAC2
    DNI63,DHI63,GHI63=csm.model_63_mac2(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 64 METSTAT
    DNI64,DHI64,GHI64=csm.model_64_metstat(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    albedo=albedo,
                                    AOD0=AOD0,
                                    Esc=S0)
    
    ### 65 PR
    DNI65,DHI65,GHI65=csm.model_65_pr(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 66 Paulesc and Schlett
    DNI66,DHI66,GHI66=csm.model_66_paulescu_schlett(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    altitude=altitude,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 67 MRM v5
    DNI67,DHI67,GHI67=csm.model_67_mrm_v5(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 68 MRM v6.1
    DNI68,DHI68,GHI68=csm.model_68_mrm_v61(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    wv=wv,
                                    ozone=u_O,
                                    albedo=albedo,
                                    AOD550=AOD550,
                                    Esc=S0)
    
    ### 69 Janjai
    DNI69,DHI69,GHI69=csm.model_69_janjai(sza=szas,
                                    date=dates,
                                    altitude=altitude,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 70 Bird
    DNI70,DHI70,GHI70=csm.model_70_bird(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    albedo=albedo,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 71 Iqbal C
    DNI71,DHI71,GHI71=csm.model_71_iqbalC(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    albedo=albedo,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 72 modified Iqbal C
    DNI72,DHI72,GHI72=csm.model_72_mod_iqbalC(sza=szas,
                                    date=dates,
                                    altitude=altitude,
                                    albedo=albedo,
                                    wv=wv,
                                    ozone=u_O,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)
    
    ### 73 REST2 v5
    DNI73,DHI73,GHI73=csm.model_73_rest2v5(sza=szas,
                                    date=dates,
                                    pressure=pressure,
                                    albedo=albedo,
                                    wv=wv,
                                    ozone=u_O,
                                    no2=u_N,
                                    AOD550=AOD550,
                                    ang=ang,
                                    Esc=S0)

    ### make plots
    for i in range(73):
        N = i+1
        if i%10==0:
            if i!=0:
                axs[0].legend()
                plt.tight_layout()
                plt.savefig(f"example_plots/example{int(i/10):d}.png")
            # plot only 10 models in one plot
            fig,axs = plt.subplots(3,1,figsize=(7,15))
            for ax in axs:
                ax.set_xlabel('sza [deg]')
                ax.grid(True)
            axs[0].set_ylabel('GHI [Wm-2]')
            axs[1].set_ylabel('DHI [Wm-2]')
            axs[2].set_ylabel('DNI [Wm-2]')

        axs[0].plot(szas,globals()[f'GHI{N:02d}'],label=f'CSM{N:02d}')
        axs[1].plot(szas,globals()[f'DHI{N:02d}'])
        axs[2].plot(szas,globals()[f'DNI{N:02d}'])
    axs[0].legend()
    plt.tight_layout()
    plt.savefig(f"example_plots/example{int(i/10):d}.png")

