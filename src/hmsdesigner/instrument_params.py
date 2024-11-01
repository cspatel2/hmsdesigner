# %%
from __future__ import annotations
import os
from typing import Dict, List, Optional, SupportsFloat as Numeric
from dataclasses import dataclass
import tosholi


@dataclass
class HmsSysParam:
    """## Instrument Optics Parameters
    """
    FlCollimator: float  # collimator focal length (mm)
    FlPrimeCamera: float  # grating camera focal length (mm)
    SlitLengthmm: float  # slit length (mm)
    sigma: float  # groove density (grooves/mm)
    relSlitPositionmm: float  # distance between the two slits (mm)
    # distance between slit closest to mosaic and the farest edge of the mosaic (mm)
    SlitA2FarEdgemm: float
    # distance between slit closest to mosaic and the closest edge of the mosaic (mm)
    SlitA2CloseEdgemm: float
    MosaicWidthmm: float  # mosaic width (mm)
    MosaicHeightmm: float  # mosaic height (mm)
    MosaicWindowWidthmm: float  # mosaic window width (mm)
    MosaicWindowHeightmm: float  # mosaic window height (mm)
    MosaicFilters: List[List[str]]  # mosaic filters


@dataclass
class HmsWlParam:
    """## Instrument interest wavelength parameters
    """
    wl: float  # wavelength (nm)
    color: str  # color
    SlitNum: int  # number of slits
    DiffractionOrder: int  # diffraction order
    PanelLetter: str  # panel letter
    PanelWindowWidthmm: float  # panel window width (mm)
    PanelWidthmm: float  # panel width (mm) < PanelWindowWidthmm
    PanelWindowHeightmm: float  # panel window height (mm)
    PanelHeightmm: float  # panel height (mm) < PanelWindowHeightmm


@dataclass
class HmsInstr:
    """## Instrument Adjustment Parameters
    """
    alpha: float  # angle of incidence (deg)
    max_ord: int  # maximum diffraction order
    mgamma_deg: float  # grating incidence angle
    imgsz: int  # detector size
    imgrot: float  # detector rotation
    slitwidth: Optional[float] = None  # slit width (mm)


@dataclass
class HmsParams:
    """## Instrument Parameters
    """
    system: str  # Instrument name
    SysParam: HmsSysParam  # Instrument Optics Parameters
    WlParam: Dict[str, HmsWlParam]  # Instrument interest wavelength parameters
    InstParam: Optional[HmsInstr] = None  # Instrument Adjustment Parameters


# %%
if __name__ == '__main__':
    hmsAOrigin_ParamDict = {
        'FlCollimator': 400,
        'FlPrimeCamera': 443.401,
        'SlitLengthmm': 55.29,
        'sigma': 98.76,
        'relSlitPositionmm': 25,
        'SlitA2FarEdgemm': 58.35,
        'SlitA2CloseEdgemm': 58.35 - 52.50,
        'MosaicWidthmm': 52.50,
        'MosaicHeightmm': 55.19,
        'MosaicWindowWidthmm': 49.52,
        'MosaicWindowHeightmm': 50.05,
        'MosaicFilters': [['4278', '6300', '5577', '7774'],
                          ['6563', '4861']]
    }

    hmsAOrigin_wlParamDict = {
        '5577': {'wl': 557.7,
                 'color': 'green',
                 'SlitNum': 2,
                 'DiffractionOrder': 33,
                 'PanelLetter': 'b',
                 'PanelWindowWidthmm': 10.10,
                 'PanelWidthmm': 10.10,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6300': {'wl': 630.0,
                 'color': 'Red',
                 'SlitNum': 2,
                 'DiffractionOrder': 29,
                 'PanelLetter': 'c',
                 'PanelWindowWidthmm': 9.70,
                 'PanelWidthmm': 9.70,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '7774': {'wl': 777.4,
                 'color': 'plum',
                 'SlitNum': 2,
                 'DiffractionOrder': 24,
                 'PanelLetter': 'a',
                 'PanelWindowWidthmm': 12.22,
                 'PanelWidthmm': 12.22,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9
                 },
        '4278': {'wl': 427.8,
                 'color': 'blue',
                 'SlitNum': 1,
                 'DiffractionOrder': 43,
                 'PanelLetter': 'd',
                 'PanelWindowWidthmm': 17.50,
                 'PanelWidthmm': 17.50 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6563': {'wl': 656.3,
                 'color': 'Darkred',
                 'SlitNum': 3,
                 'DiffractionOrder': 28,
                 'PanelLetter': 'f',
                 'PanelWindowWidthmm': 18.83,
                 'PanelWidthmm': 18.83 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
        '4861': {'wl': 486.1,
                 'color': 'cyan',
                 'SlitNum': 4,
                 'DiffractionOrder': 38,
                 'PanelLetter': 'e',
                 'PanelWindowWidthmm': 30.69,
                 'PanelWidthmm': 30.69,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
    }

    hmsA_ParamDict = {
        'hmsVersion': 'A',
        'FlCollimator': 400,
        'FlPrimeCamera': 443.401,
        'SlitLengthmm': 55.29,
        'sigma': 98.76,
        'relSlitPositionmm': 21.43,
        'SlitA2FarEdgemm': 63.35,
        'SlitA2CloseEdgemm': 13.83,
        'MosaicWidthmm': 52.50,
        'MosaicHeightmm': 55.19,
        'MosaicWindowWidthmm': 49.52,
        'MosaicWindowHeightmm': 50.05,
        'MosaicFilters': [['4278', '6300', '5577', '7774'],
                          ['6563', '4861']]
    }

    hmsA_wlParamDict = {
        '5577': {'wl': 557.7,
                 'color': 'green',
                 'SlitNum': 2,
                 'DiffractionOrder': 33,
                 'PanelLetter': 'b',
                 'PanelWindowWidthmm': 10.10,
                 'PanelWidthmm': 10.10,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6300': {'wl': 630.0,
                 'color': 'Red',
                 'SlitNum': 2,
                 'DiffractionOrder': 29,
                 'PanelLetter': 'c',
                 'PanelWindowWidthmm': 9.70,
                 'PanelWidthmm': 9.70,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '7774': {'wl': 777.4,
                 'color': 'plum',
                 'SlitNum': 1,
                 'DiffractionOrder': 24,
                 'PanelLetter': 'a',
                 'PanelWindowWidthmm': 12.22,
                 'PanelWidthmm': 12.22,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9
                 },
        '4278': {'wl': 427.8,
                 'color': 'blue',
                 'SlitNum': 1,
                 'DiffractionOrder': 43,
                 'PanelLetter': 'd',
                 'PanelWindowWidthmm': 17.50,
                 'PanelWidthmm': 17.50 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6563': {'wl': 656.3,
                 'color': 'Darkred',
                 'SlitNum': 3,
                 'DiffractionOrder': 28,
                 'PanelLetter': 'f',
                 'PanelWindowWidthmm': 18.83,
                 'PanelWidthmm': 18.83 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
        '4861': {'wl': 486.1,
                 'color': 'cyan',
                 'SlitNum': 4,
                 'DiffractionOrder': 38,
                 'PanelLetter': 'e',
                 'PanelWindowWidthmm': 30.69,
                 'PanelWidthmm': 30.69,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
    }

    hmsAEclipse_wlParamDict = {
        '5577': {'wl': 557.7,
                 'color': 'green',
                 'SlitNum': 2,
                 'DiffractionOrder': 33,
                 'PanelLetter': 'b',
                 'PanelWindowWidthmm': 10.10,
                 'PanelWidthmm': 10.10,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6300': {'wl': 630.0,
                 'color': 'Red',
                 'SlitNum': 2,
                 'DiffractionOrder': 29,
                 'PanelLetter': 'c',
                 'PanelWindowWidthmm': 9.70,
                 'PanelWidthmm': 9.70,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '7774': {'wl': 777.4,
                 'color': 'plum',
                 'SlitNum': 1,
                 'DiffractionOrder': 24,
                 'PanelLetter': 'a',
                 'PanelWindowWidthmm': 12.22,
                 'PanelWidthmm': 12.22,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9
                 },
        '4278': {'wl': 427.8,
                 'color': 'blue',
                 'SlitNum': 1,
                 'DiffractionOrder': 43,
                 'PanelLetter': 'd',
                 'PanelWindowWidthmm': 17.50,
                 'PanelWidthmm': 17.50 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+2.9,
                 },
        '6563': {'wl': 656.3,
                 'color': 'Darkred',
                 'SlitNum': 3,
                 'DiffractionOrder': 28,
                 'PanelLetter': 'f',
                 'PanelWindowWidthmm': 18.83,
                 'PanelWidthmm': 18.83 + 2.68,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
        '4861': {'wl': 486.1,
                 'color': 'cyan',
                 'SlitNum': 4,
                 'DiffractionOrder': 38,
                 'PanelLetter': 'e',
                 'PanelWindowWidthmm': 30.69,
                 'PanelWidthmm': 30.69,
                 'PanelWindowHeightmm': 25.02,
                 'PanelHeightmm': 25.02+1.94,
                 },
    }

    hmsAEclipse_ParamDict = {
        'hmsVersion': 'A - ECLIPSE',
        'FlCollimator': 400,
        'FlPrimeCamera': 429.420,
        'SlitLengthmm': 55.29,
        'sigma': 98.76,
        'relSlitPositionmm': 21.43,
        'SlitA2FarEdgemm': 63.35,
        'SlitA2CloseEdgemm': 13.83,
        'MosaicWidthmm': 52.50,
        'MosaicHeightmm': 55.19,
        'MosaicWindowWidthmm': 49.52,
        'MosaicWindowHeightmm': 50.05,
        'MosaicFilters': [['4278', '6300', '5577', '7774'],
                          ['6563', '4861']]
    }

    hmsB_ParamDict = {
        'hmsVersion': 'B',  # hitmis stayed here for solar eclipse
        'FlCollimator': 400,  # focal length of collimator, mm.
        # focal length of collimator (grating -> mosaic), mm.
        'FlPrimeCamera': 365.891,
        'SlitLengthmm': 55.30,  # slight length, mm.
        # slit length = slitlen/focallength,Deg.
        'sigma': 1/(10245.208387*100/80.5),  # measured grating density
        'relSlitPositionmm': 20.04,  # Distance between the two slits, mm.
        # distance between slit closest to mosaic and the farest edge of the mosaic, mm.
        'SlitA2FarEdgemm': 71.35,
        # distance between slit closest to mosaic and the closest edge of the mosaic, mm.
        'SlitA2CloseEdgemm': 15.83,
        'MosaicWidthmm': 52.52,
        'MosaicHeightmm': 55.10,
        # Blue tap facing away from grating, if we are looking from grating to mosaic
        'MosaicWindowWidthmm': 49.96,  # The window that the camera images, mm.
        # the window height that the camera images, mm.
        'MosaicWindowHeightmm': 50.15,
        'MosaicFilters': [['4278', '7774', '6300', '5577'],  # top panel from L -> R
                          ['4861', '6563']]  # bottom panel from L -> R
    }

    hmsB_wlParamDict = {
        '5577': {'wl': 557.7,
                 'color': 'green',
                 'SlitNum': 2,
                 'DiffractionOrder': 33,
                 'PanelLetter': 'a',
                 'PanelWindowWidthmm': 18.18 + 1.44,
                 'PanelWidthmm': 18.18 + 1.44,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 },
        '6300': {'wl': 630.0,
                 'color': 'Red',
                 'SlitNum': 2,
                 'DiffractionOrder': 29,
                 'PanelLetter': 'b',
                 'PanelWindowWidthmm': 7.09+1.44,
                 'PanelWidthmm': 7.09+1.44,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 },
        '7774': {'wl': 777.4,
                 'color': 'plum',
                 'SlitNum': 2,
                 'DiffractionOrder': 24,
                 'PanelLetter': 'c',
                 'PanelWindowWidthmm': 8.39+1.92,
                 'PanelWidthmm': 8.39+1.92,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 },
        '4278': {'wl': 427.8,
                 'color': 'blue',
                 'SlitNum': 1,
                 'DiffractionOrder': 43,
                 'PanelLetter': 'd',
                 'PanelWindowWidthmm': 12.89,  # measure the window to + or -
                 'PanelWidthmm': 12.89,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 },
        '6563': {'wl': 656.3,
                 'color': 'Darkred',
                 'SlitNum': 3,
                 'DiffractionOrder': 28,
                 'PanelLetter': 'e',
                 'PanelWindowWidthmm': 28.89 + 1.96,
                 'PanelWidthmm': 28.89 + 1.96,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 },
        '4861': {'wl': 486.1,
                 'color': 'cyan',
                 'SlitNum': 4,
                 'DiffractionOrder': 38,
                 'PanelLetter': 'f',
                 'PanelWindowWidthmm': 20.61,  # measure the window to + or -
                 'PanelWidthmm': 20.61,
                 'PanelWindowHeightmm': 26.69,  # measure the window to + or -
                 'PanelHeightmm': 26.69,
                 }
    }

    hmsBOrigin_ParamDict = {
        'hmsVersion': 'B-ORIGIN',
        'FlCollimator': 400,  # focal length of collimator, mm.
        # focal length of collimator (grating -> mosaic), mm.
        'FlPrimeCamera': 365.891,
        'SlitLengthmm': 55.30,  # slight length, mm.
        # slit length = slitlen/focallength,Deg.
        'sigma': 1 / 10245.208387,  # measured grating density
        #     'sigma': 10131.7122594, #new grating
        'relSlitPositionmm': 20.04,  # Distance between the two slits, mm.
        # distance between slit closest to mosaic and the farest edge of the mosaic, mm.
        'SlitA2FarEdgemm': 71.35,
        # distance between slit closest to mosaic and the closest edge of the mosaic, mm.
        'SlitA2CloseEdgemm': 15.83,
        'MosaicWidthmm': 52.52,
        'MosaicHeightmm': 55.10,
        # Blue tap facing away from grating, if we are looking from grating to mosaic
        'MosaicWindowWidthmm': 49.96,  # The window that the camera images, mm.
        # the window height that the camera images, mm.
        'MosaicWindowHeightmm': 50.15,
        'MosaicFilters': [['7841', '4278', '6300', '5577'],  # top panel from L -> R
                          ['6563', '4861']]  # ,'6544']]  #bottom panel from L -> R
    }

    hmsBOrigin_wlParamDict = {
        '5577': {'wl': 557.7,
                 'color': 'green',
                 'SlitNum': 2,
                 'DiffractionOrder': 33,
                 'PanelLetter': 'a',
                 'PanelWindowWidthmm': 14.75,
                 'PanelWidthmm': 14.75,
                 'PanelWindowHeightmm': 24.775,  # measure the window to + or -
                 'PanelHeightmm': 24.775+2.64,
                 },
        '6300': {'wl': 630.0,
                 'color': 'Red',
                 'SlitNum': 2,
                 'DiffractionOrder': 29,
                 'PanelLetter': 'b',
                 'PanelWindowWidthmm': 6.5,  # 1.44 ,
                 'PanelWidthmm': 6.5,
                 'PanelWindowHeightmm': 24.775,  # measure the window to + or -
                 'PanelHeightmm': 24.775+2.64,
                 },
        '4278': {'wl': 427.8,
                 'color': 'blue',
                 'SlitNum': 1,
                 'DiffractionOrder': 43,
                 'PanelLetter': 'd',
                 'PanelWindowWidthmm': 6,  # measure the window to + or -
                 'PanelWidthmm': 6,
                 'PanelWindowHeightmm': 24.775,  # measure the window to + or -
                 'PanelHeightmm': 24.775+2.64,
                 },
        '7841': {'wl': 784.1,
                 'color': 'hotpink',
                 'SlitNum': 2,
                 'DiffractionOrder': 23,
                 'PanelLetter': 'd',
                 # measure the window to + or -
                 'PanelWindowWidthmm': 14.5+8.21,
                 'PanelWidthmm': 14.5 + 8.21 + 2.56,
                 'PanelWindowHeightmm': 24.775,  # measure the window to + or -
                 'PanelHeightmm': 24.775+2.64,
                 },

        #     '7774':{'wl':777.4,
        #             'color':'plum',
        #             'SlitNum':2,
        #             'DiffractionOrder': 23,
        #             'PanelLetter':'c',
        #             'PanelWindowWidthmm':8.21,
        #             'PanelWidthmm': 8.21 + 2.56,
        #             'PanelWindowHeightmm':24.775, #measure the window to + or -
        #             'PanelHeightmm':24.775+2.64,
        #             },
        '4861': {'wl': 486.1,
                 'color': 'cyan',
                 'SlitNum': 4,
                 'DiffractionOrder': 38,
                 'PanelLetter': 'f',
                 'PanelWindowWidthmm': 21,  # measure the window to + or -
                 'PanelWidthmm': 21,
                 'PanelWindowHeightmm': 25.375,  # measure the window to + or -
                 'PanelHeightmm': 25.375 + 2.31,
                 },
        '6563': {'wl': 656.3,
                 'color': 'Darkred',
                 'SlitNum': 3,
                 'DiffractionOrder': 28,
                 'PanelLetter': 'e',
                 'PanelWindowWidthmm': 28.96,
                 'PanelWidthmm': 28.96 + 2.56,
                 'PanelWindowHeightmm': 25.375,  # measure the window to + or -
                 'PanelHeightmm': 25.375 + 2.31,
                 },

    }

    sys = HmsSysParam(**hmsAOrigin_ParamDict)
    wls = {k: HmsWlParam(**v) for (k, v) in hmsAOrigin_wlParamDict.items()}
    hmsparam = HmsParams(sys.hmsVersion, sys, wls)
    # print(hash(sys))
    dir = os.path.dirname(__file__)
    with open(os.path.join(dir, 'hmsa_origin.toml'), 'wb') as ofile:
        tosholi.dump(hmsparam, ofile)

    sys = HmsSysParam(**hmsA_ParamDict)
    wls = {k: HmsWlParam(**v) for (k, v) in hmsA_wlParamDict.items()}
    hmsparam = HmsParams(sys.hmsVersion, sys, wls)
    # print(hash(sys))
    # print(hmsparam.to_dict())
    with open(os.path.join(dir, 'hmsa_aurora.toml'), 'wb') as ofile:
        tosholi.dump(hmsparam, ofile)

    sys = HmsSysParam(**hmsAEclipse_ParamDict)
    wls = {k: HmsWlParam(**v) for (k, v) in hmsAEclipse_wlParamDict.items()}
    hmsparam = HmsParams(sys.hmsVersion, sys, wls)
    with open(os.path.join(dir, 'hmsa_eclipse.toml'), 'wb') as ofile:
        tosholi.dump(hmsparam, ofile)

    sys = HmsSysParam(**hmsB_ParamDict)
    wls = {k: HmsWlParam(**v) for (k, v) in hmsB_wlParamDict.items()}
    hmsparam = HmsParams(sys.hmsVersion, sys, wls)
    with open(os.path.join(dir, 'hmsb_old.toml'), 'wb') as ofile:
        tosholi.dump(hmsparam, ofile)

    sys = HmsSysParam(**hmsBOrigin_ParamDict)
    wls = {k: HmsWlParam(**v) for (k, v) in hmsBOrigin_wlParamDict.items()}
    hmsparam = HmsParams(sys.hmsVersion, sys, wls)
    with open(os.path.join(dir, 'hmsb_origin.toml'), 'wb') as ofile:
        tosholi.dump(hmsparam, ofile)

# %%
