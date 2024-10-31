# %%
from __future__ import annotations
from typing import Dict, Optional, SupportsFloat as Numeric
import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from collections.abc import Iterable
import tosholi
import astropy.io.fits as pf
from scipy.optimize import curve_fit
from skimage.exposure import equalize_hist

from .instrument_params import HmsParams, HmsSysParam, HmsWlParam, HmsInstr
# %%
# %%


class HMS_ImagePredictor:
    """## HMS Image Predictor class.
    """

    def __init__(self, configfile: str, alpha: Optional[Numeric] = None, *, num_orders: Optional[int] = None, mgammadeg: Optional[Numeric] = None, pix: Optional[int] = None, img_rot: Optional[Numeric] = None, slitwidth: Optional[Numeric] = None):
        """## Create an instance of the HMS_ImagePredictor class.
        Applies the grating equation (sinβ = mλ/(d*sinγ) - sinα) to predict the location of spectral lines on the detector.

        ### Args:
            - `configfile (str)`: Configuration TOML file for the instrument.
            - `alpha (Optional[Numeric], optional)`: Override α. Defaults to None.
            - `num_orders (Optional[int], optional)`: Maximum order to calculate. Defaults to None.
            - `mgammadeg (Optional[Numeric], optional)`: Override γ. Defaults to None.
            - `pix (Optional[int], optional)`: Detector pixels. Defaults to None.
            - `img_rot (Optional[Numeric], optional)`: Image plane rotation. Defaults to None.
            - `slitwidth (Optional[Numeric], optional)`: Slit width override, in mm. Defaults to None.

        ### Raises:
            - `FileNotFoundError`: Configuration file not found.
            - `TypeError`: Invalid file extension.
            - `ValueError`: Alpha, num_orders, mgammadeg, and pix must be set.
        """
        self.hmsVersion = ''
        self.slitwidth: Optional[Numeric] = None
        if not os.path.exists(configfile):
            raise FileNotFoundError(f"File {configfile} not found.")
        
        self.hmsParamDict: HmsSysParam = None
        self.wlParamDict: Dict[str, HmsWlParam] = None

        ext = os.path.splitext(configfile)[-1].lower()
        if ext == '.toml':
            with open(configfile, 'rb') as ifile:
                params = tosholi.load(HmsParams, ifile)
                self.hmsParamDict = params.SysParam
                self.wlParamDict = params.WlParam
                self.hmsVersion = params.system
                if params.InstParam is not None:
                    instparam: HmsInstr = params.InstParam
                    self.alpha = instparam.alpha
                    self.num_orders = instparam.max_ord
                    self.mgammadeg = instparam.mgamma_deg
                    self.pix = instparam.imgsz
                    self.img_rot = instparam.imgrot
                    self.slitwidth = instparam.slitwidth
        else:
            raise TypeError(
                f"Invalid file extension {ext}. Please provide a .toml file.")

        # Load params if they are some
        if alpha is not None:
            self.alpha = alpha
        if num_orders is not None:
            self.num_orders = num_orders
        if mgammadeg is not None:
            self.mgammadeg = mgammadeg
        if pix is not None:
            self.pix = pix
        if img_rot is not None:
            self.img_rot = img_rot
        if slitwidth is not None:
            self.slitwidth = slitwidth

        if self.alpha is None or self.num_orders is None or self.mgammadeg is None or self.pix is None:
            raise ValueError(
                "Alpha, num_orders, mgammadeg, and pix must be set.")
        # Initialize other parameters
        # focal length of collimator slit -> grating
        self.f = self.hmsParamDict.FlCollimator
        # focal length of collimator grating -> mosaic
        self.fprime = self.hmsParamDict.FlPrimeCamera
        self.sigma = 1e6 / self.hmsParamDict.sigma # mm
        self.alphas = self.relative_alphas(self.alpha)  # deg
        self.alpha_slitA = self.alphas[0]  # deg
        self.alpha_slitB = self.alphas[1]  # deg
        slitlengthdeg = np.arctan(
            self.hmsParamDict.SlitLengthmm*0.5 / self.f)*360/np.pi  # deg
        self.gamma = self.mgammadeg + \
            ((slitlengthdeg / 2) * np.linspace(-1, 1, 2 * 50 + 1))  # deg
        # all int orders
        self.orders = np.arange(-self.num_orders,
                                self.num_orders + 1, 1, dtype=int)
        self.MosaicWindowHeightmm = self.hmsParamDict.MosaicWindowHeightmm
        self.MosaicWindowWidthmm = self.hmsParamDict.MosaicWindowWidthmm
        self.MosaicHeightmm = self.hmsParamDict.MosaicHeightmm
        self.MosaicWidthmm = self.hmsParamDict.MosaicWidthmm

        self.lw = 0.4  # line width
        self.ms = 0.2  # maker size
        self.ls = '--o'  # line style
        self.lc = 'black'  # line color

    def get_beta(self, wl: Numeric, m: int, *, gamma: Numeric = 90) -> Numeric:
        """## Calculate the diffracted angle β using the grating equation: \n sinβ = mλ/(d*sinγ) - sinα

        ### Args:
            - `wl (Numeric)`: Wavelength in nm.
            - `m (int)`: Diffraction order.
            - `gamma (Numeric, optional)`: Incident angle parallel to the grating. Defaults to 90.

        ### Raises:
            - `ValueError`: α is not set.

        ### Returns:
            - `Numeric`: Diffracted angle β.
        """
        if self.alpha is None:
            raise ValueError("Alpha is not set. Please set the alpha value.")
        const = m*wl/(self.sigma*np.sin(np.deg2rad(gamma)))  # unitless
        sinb = const - np.sin(np.deg2rad(self.alpha))  # unitless
        beta = np.arcsin(sinb)  # rad
        return np.rad2deg(beta)

    def resolving_power(self, wl: Numeric, m: int, *, gamma: Numeric = 90) -> Numeric:
        """## Calculate the resolving power of the HMS using the formula: \n R = mλ/(σ*sin γ*cos β*B) \n where B is the slit width in radians.

        ### Args:
            - `wl (Numeric)`: Wavelength in nm.
            - `m (int)`: Diffraction order.
            - `gamma (Numeric, optional)`: Incident angle parallel to the grating. Defaults to 90.

        ### Raises:
            - `ValueError`: α is not set.
            - `ValueError`: Slit width is not set.

        ### Returns:
            - `Numeric`: Resolving power (λ/Δλ).
        """
        if self.alpha is None:
            raise ValueError("Alpha is not set. Please set the alpha value.")
        if self.slitwidth is None:
            raise ValueError(
                "Slit width is not set. Please set the slit width.")
        blur = self.slitwidth/self.f  # rad
        const = m*wl/(self.sigma*np.sin(np.deg2rad(gamma)))  # unitless
        sinb = const - np.sin(np.deg2rad(self.alpha))  # unitless
        cosb = np.sqrt(1 - sinb**2)  # unitless
        res = (const/cosb) / blur
        return res

    # Spatial to angular coordinates
    def mm2deg(self, len_mm: float, fl: float) -> float:
        """converts linear postion (mm) to angluar postion (degrees)

        Args:
            len_mm (float): linear postion in mm. 
            fl (float): focal length of the lens in mm.

        Returns:
            float: angluar postion in degrees.
        """
        return np.rad2deg(len_mm / fl)

    # Angular to spatial coordinates
    def deg2mm(self, beta_deg: float, fl: float) -> float:
        """ converts angluar postion (degrees) to linear postion (mm)

        Args:
            beta_deg (float): angluar postion in degrees.
            fl (float): focal length of the lens in mm.

        Returns:
            float: linear postion in mm.
        """
        return fl * np.deg2rad(beta_deg)

    # Angular to detector coordinates
    def deg2pix(self, deg: float, axis: int, fl: float, totalpix: int = None) -> float:
        """ Converts angluar postion (degrees) to pixel postion (pix)
        Args:
            deg (float): angluar postion in degrees.
            axis (int): 0 for x-axis, 1 for y-axis.
            fl (float): focal length of the lens in mm.
            totalpix (int, optional): toal number of pixels on the detector along the axis. if None-it defualts to the detector size used to iniitalize the class. Defaults to None.

        Raises:
            ValueError: Axis has to be 0 [x-axis] or 1 [y-axis]

        Returns:
            float: pixel postion (pix)
        """
        if totalpix is None:
            totalpix = self.pix
        if axis == 0:
            mlen = self.MosaicWindowWidthmm  # axis 0 = x, mlen = width of mosaic window of HMS
        elif axis == 1:
            # axis 1 = y, mlen = height of mosaic window of HMS
            mlen = self.MosaicWindowHeightmm
        else:
            raise ValueError("Axis can be 0 [x-axis] or 1 [y-axis]")
        return self.deg2mm(deg, fl) * totalpix / mlen

    # Spatial to pixel coordinates
    def mm2pix(self, len_mm: float, axis: int, totalpix: int = None) -> float:
        if totalpix is None:
            totalpix = self.pix

        if axis == 0:
            mlen = self.MosaicWindowWidthmm  # axis 0 = x, mlen = width of mosaic window of HMS  # noqa: E701
        elif axis == 1:
            # axis 1 = y, mlen = height of mosaic window of HMS
            mlen = self.MosaicWindowHeightmm
        else:
            raise ValueError("Axis can be 0 [x-axis] or 1 [y-axis]")
        return len_mm * totalpix / mlen

    def relative_alphas(self, alpha):
        alpha_slitA = alpha  # deg
        alpha_slitB = alpha_slitA + \
            self.mm2deg(self.hmsParamDict.relSlitPositionmm, self.f)  # deg
        return [alpha_slitA, alpha_slitB]  # deg

    def annotate_width(self, ax, lenmm, x_start, y0, wl, linewidth, linestyle, color, measurement: bool, toppanel: bool = True):
        if toppanel:
            ymin = 0.5
            ymax = 1
            ha = 'left'
            va = 'bottom'
            y = y0-12

        else:
            ymin = 0
            ymax = 0.5
            ha = 'left'
            va = 'bottom'
            y = y0 - 12
        ax.axvline(x=x_start, ymin=ymin, ymax=ymax,
                   linewidth=linewidth, linestyle=linestyle, color=color)
        plt.text(x_start+lenmm-.5, y, f'Panel:{wl} A', color='black', fontsize=7, ha=ha, va=va, rotation=90,
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.3'))
        xy = (x_start, y0)
        xytext = (x_start + lenmm, y0)
        if measurement:
            plt.annotate('', xy=xy, xytext=xytext,
                         arrowprops=dict(arrowstyle='<->', color='black', lw=1))
            plt.text(x_start + lenmm/2, y0, f'{lenmm:.2f} mm',
                     color='black', fontsize=7, ha='center', va='center',
                     bbox=dict(facecolor='white', alpha=1, edgecolor='none', boxstyle='round,pad=0.3'))

    def annotate_height(self, ax, lenmm, y_start, x0, linewidth, linestyle, color, measurement: bool):
        ax.axhline(y=y_start+lenmm, linewidth=linewidth,
                   linestyle=linestyle, color=color)
        xy = (x0, y_start)
        xytext = (x0, y_start+lenmm)
        if measurement:
            plt.annotate('', xy=xy, xytext=xytext,
                         arrowprops=dict(arrowstyle='<->', color='black', lw=1))
            plt.text(x0, y_start + lenmm/5, f'{lenmm:.2f} mm',
                     color='black', fontsize=7, ha='center', va='center', rotation=270,
                     bbox=dict(facecolor='white', alpha=1, edgecolor='none', boxstyle='round,pad=0.3'))

    def Linear(self, x: float, m: float, b: float): return m*x+b

    def single_defraction(self, alpha: float, gamma: float, m_order: int, sigma: float, wl: float) -> float:
        """ Calculates difraction angle beta from a diffraction grating using the grating equation. \n sinβ = mλ/(d*sinγ) - sinα 

        Args:
            alpha (float): angle of incidence perpendicluar to groves, α [Degrees (0-360) or Radians(0-2π)].
            gamma (float): angle of incidence parallel to groves,γ [Degrees (0-360) or Radians(0-2π)]. 
            m_order (int):  difraction order, m.
            sigma (float): grating grove spacing (nm). or 1e6/Grovedensity(grooves/mm)
            wl (float): Wavelength, λ.

        Returns:
            float: angle of difraction, β [Degrees (0-360)].
        """
        if isinstance(m_order, Iterable):
            m_iter = list(m_order)
        else:
            m_iter = [m_order]

        if isinstance(gamma, Iterable):
            beta_arr = []
            for m in m_iter:
                betas = [self.single_defraction(
                    alpha, g_, m, sigma, wl) for g_ in gamma]
                beta_arr.append(betas)
            return np.asanyarray(beta_arr)

        alpha_rad = alpha*np.pi/180
        gamma_rad = gamma*np.pi/180  # noqa: F405
        beta_rad = np.arcsin(-np.sin(alpha_rad) +
                             ((m_order*wl)/(self.sigma*np.sin(gamma_rad))))

        return np.rad2deg(beta_rad)

    def Grating(self, alpha: float, gamma: float, m_order: int, sigma: float, wl: float) -> np.asanyarray:
        """ Calulates all possible diffraction angles for orders up to -m to +m using the grating equation: \n sinβ = mλ/(d*sinγ) - sinα 

        Args:
            alpha (float): angle of incidence perpendicluar to groves, α [Degrees (0-360) or Radians(0-2π)].
            gamma (float): angle of incidence parallel to groves,γ [Degrees (0-360) or Radians(0-2π)]. 
            m_order (int): difraction order, m.
            sigma (float): grating grove spacing (nm). or 1e6/Grovedensity(grooves/mm)
            wl (float): Wavelength, λ.

        Returns:
            np.asanyarray: diffraction angles, array of shape (len(wl), len(alpha), 2*m_order + 1, len(gamma))
        """
        if isinstance(alpha, Iterable):
            alpha_iter = list(alpha)
        else:
            alpha_iter = [alpha]

        if isinstance(gamma, Iterable):
            gamma_iter = list(gamma)
        else:
            gamma_iter = [gamma]

        if isinstance(wl, Iterable):
            wl_iter = list(wl)
        else:
            wl_iter = [wl]

        results = []
        for wl in wl_iter:
            res_per_wl = []
            for alpha in alpha_iter:
                beta = self.single_defraction(
                    alpha, gamma_iter, m_order, sigma, wl)
                res_per_wl.append(beta)
            results.append(res_per_wl)
        # shape(result) = (len(wl), len(alpha), 2*m_order + 1, len(gamma))
        return np.asanyarray(results)

    def plot_spectral_lines(self, ImageAt: str, Tape2Grating: bool, mosaic: bool = True, wls: Iterable = [486.1, 427.8, 557.7, 630.0, 656.3, 777.4], measurement: bool = True, fprime: float = None, *, fig_args = (), fig_kwargs = dict(figsize=(7, 6), dpi=300, tight_layout=True)) -> plt.Figure:

        self.wls = wls
        ImageAt = ImageAt.lower()

        if fprime is not None:
            self.fprime = fprime

        # Calculate beta values using grating equation
        betas = self.Grating(self.alphas, self.gamma,
                             self.orders, self.sigma, wls)  # deg
        betas = self.deg2mm(betas, fl=self.fprime)  # mm
        gamma = self.deg2mm(self.gamma, fl=self.fprime)  # mm

        fig, ax = plt.subplots(*fig_args, **fig_kwargs)
        plt.rcParams['font.family'] = 'monospace'
        plt.rcParams['font.monospace'] = 'Andale Mono'
        matplotlib.rc('text', usetex=False)
        matplotlib.rc('xtick', labelsize=10)
        matplotlib.rc('ytick', labelsize=10)
        version = self.hmsVersion.upper()

        ax.set_title(
            f"{version} \nDiffracted Spectral Lines at {ImageAt.upper()} \nAlpha = {self.alpha} Deg")

        if ImageAt == 'mosaic':
            ax.set_xlabel('Beta [mm]')
            ax.set_ylabel('Gamma [mm]')
            self._plot_on_mosaic(ax, betas, gamma, wls, Tape2Grating)

        elif ImageAt == 'mosaicwindow':
            ax.set_xlabel('Beta [mm]')
            ax.set_ylabel('Gamma [mm]')
            self._plot_on_mosaicwindow(
                ax, betas, gamma, wls, Tape2Grating, mosaic, measurement)

        elif ImageAt == 'detector':
            ax.set_xlabel('Detector X [pix]')
            ax.set_ylabel('Detector Y [pix]')
            self._plot_on_detector(ax, betas, gamma, wls, Tape2Grating, fprime)

        else:
            ax.set_xlabel('Beta [mm]')
            ax.set_ylabel('Gamma [mm]')
            self._plot_all_lines_on_mosaic(ax, betas, gamma, wls, Tape2Grating)

        return fig

    def _plot_all_lines_on_mosaic(self, ax, betas, gamma, wls, Tape2Grating):
        for widx, wl in enumerate(wls):
            w = str(int(wl * 10))  # nm -> A
            if w not in self.wlParamDict.keys():
                w = self.find_nearest(wl)
                wdict = self.wlParamDict.get(w)
                color = 'Black'
                gamma_plot = gamma
                for aidx, a in enumerate(self.alphas):
                    for midx, m in enumerate(self.orders):
                        beta_plot = betas[widx][aidx][midx]
                        ax.plot(beta_plot, gamma_plot, self.ls,
                                markersize=self.ms, color=color, linewidth=self.lw)
                        idx = int(0.1*len(beta_plot))
                        xy = (beta_plot[idx], gamma_plot[30])
                        xytext = (beta_plot[idx]-1, gamma_plot[idx]+1)
                        ax.annotate(f'[{aidx},{m}]', xy, xytext, rotation=270)
                        idx = int(0.8*len(beta_plot))
                        xy = (beta_plot[idx], gamma_plot[-30])
                        xytext = (beta_plot[idx]-1, gamma_plot[idx]+1)
                        ax.annotate(f'[{aidx},{m}]', xy, xytext, rotation=270)
            else:
                wdict = self.wlParamDict.get(w)
                color = wdict.color
                if wdict.SlitNum in [1, 3]:
                    aidx = 1
                elif wdict.SlitNum in [2, 4]:
                    aidx = 0
                else:
                    raise ValueError("Slit Number must be 1, 2, 3, or 4.")

                # set gamma filter accoding to SlitNum
                self.g0 = self.deg2mm(self.mgammadeg, fl=self.fprime)
                if wdict.SlitNum in [1, 2]:
                    gmask = (gamma >= self.g0)
                else:
                    gmask = (gamma <= self.g0)
                gamma_plot = gamma[gmask]

                for midx, m in enumerate(self.orders):
                    beta_plot = betas[widx][aidx][midx][gmask]
                    ax.plot(beta_plot, gamma_plot, self.ls,
                            markersize=self.ms, color=color, linewidth=self.lw)

        self.g0 = self.deg2mm(self.mgammadeg, fl=self.fprime)
        # To plot just the Mosaic
        X1 = self.deg2mm(self.alpha, fl=self.fprime) - \
            self.hmsParamDict.SlitA2FarEdgemm  # mm
        X2 = X1 + self.MosaicWindowWidthmm  # mm

        Y1 = self.g0 + self.MosaicWindowHeightmm/2  # mm
        Y2 = self.g0 - self.MosaicWindowHeightmm/2  # mm

        # ax.axhline(y = 92, color = 'Black',linewidth = 1)

        # limit from L-> R if looking from grating to mosaic.
        ax.set_xlim(X2, X1)
        print(X2, X1)
        ax.set_ylim(Y2, Y1)
        print(Y2, Y1)

    def find_nearest(self, wl):
        wlnum = np.array([int(w) / 10 for w in self.wlParamDict.keys()])
        nearest_idx = np.argmin(np.abs(wlnum - wl))
        return list(self.wlParamDict.keys())[nearest_idx]

    def _plot_on_mosaic(self, ax, betas, gamma, wls, Tape2Grating):
        for widx, wl in enumerate(wls):
            w = str(int(wl * 10))  # nm -> A

            if w not in self.wlParamDict.keys():
                w = self.find_nearest(wl)
                wdict = self.wlParamDict.get(w)
                color = 'Black'
            else:
                wdict = self.wlParamDict.get(w)
                color = wdict.color

            # set alpha idex
            if wdict.SlitNum in [1, 3]:
                aidx = 1
            elif wdict.SlitNum in [2, 4]:
                aidx = 0
            else:
                raise ValueError("Slit Number must be 1, 2, 3, or 4.")

            # set gamma filter accoding to SlitNum
            self.g0 = self.deg2mm(self.mgammadeg, fl=self.fprime)
            if wdict.SlitNum in [1, 2]:
                gmask = (gamma >= self.g0)
            else:
                gmask = (gamma <= self.g0)
            gamma_plot = gamma[gmask]

            for midx, m in enumerate(self.orders):
                beta_plot = betas[widx][aidx][midx][gmask]
                ax.plot(beta_plot, gamma_plot, self.ls,
                        markersize=self.ms, color=color, linewidth=self.lw)

                idx = int(0.8*len(beta_plot))
                xy = (beta_plot[idx], gamma_plot[-30])
                xytext = (beta_plot[idx]-.75, gamma_plot[idx]+4)
                slitnum = int(wdict.SlitNum)
                ax.annotate(f'[{slitnum}, {int(wl*10)} A]', xy, xytext,
                            rotation=270, va='top', ha='center', fontsize=8)

         # To plot just the Mosaic
        X1 = self.deg2mm(self.alpha, fl=self.fprime) - \
            self.hmsParamDict.SlitA2FarEdgemm  # mm
        X2 = X1 + self.MosaicWindowWidthmm  # mm

        Y1 = self.g0 + self.MosaicWindowHeightmm/2  # mm
        Y2 = self.g0 - self.MosaicWindowHeightmm/2  # mm

        # ax.axhline(y = 92, color = 'Black',linewidth = 1)

        # limit from L-> R if looking from grating to mosaic.
        ax.set_xlim(X2, X1)
        ax.set_ylim(Y2, Y1)

    def _plot_on_mosaicwindow(self, ax, betas, gamma, wls, Tape2Grating, Mosaic=True, Measurements: bool = True, savefig: bool = True):
        for widx, wl in enumerate(wls):
            w = str(int(wl * 10))  # nm -> A

            if w not in self.wlParamDict.keys():
                w = self.find_nearest(wl)
                wdict = self.wlParamDict.get(w)
                color = 'Black'
            else:
                wdict = self.wlParamDict.get(w)
                color = wdict.color

            # set alpha idex
            if wdict.SlitNum in [1, 3]:
                aidx = 1
            elif wdict.SlitNum in [2, 4]:
                aidx = 0
            else:
                raise ValueError("Slit Number must be 1, 2, 3, or 4.")

            # set gamma filter accoding to SlitNum
            self.g0 = self.deg2mm(self.mgammadeg, fl=self.fprime)
            if wdict.SlitNum in [1, 2]:
                gmask = (gamma >= self.g0)
            else:
                gmask = (gamma <= self.g0)
            gamma_plot = gamma[gmask]

            midx = list(self.orders).index(wdict.DiffractionOrder)
            beta_plot = betas[widx][aidx][midx][gmask]
            ax.plot(beta_plot, gamma_plot, self.ls,
                    markersize=self.ms, color=color, linewidth=self.lw)
            idx = int(0.8*len(beta_plot))
            xy = (beta_plot[idx], gamma_plot[-30])
            xytext = (beta_plot[idx]-.75, gamma_plot[idx]+4)
            slitnum = int(wdict.SlitNum)
            ax.annotate(f'[{slitnum}, {int(wl*10)} A]', xy, xytext,
                        rotation=270, va='top', ha='center', fontsize=8)

            # ax.annotate(f'[{int(wl*10)} A,{self.orders[midx]}]',xy,xytext,rotation = 270,va='top', ha = 'center',fontsize = 8)

        # To plot just the Mosaic
        X1 = self.deg2mm(self.alpha, fl=self.fprime) - \
            self.hmsParamDict.SlitA2FarEdgemm  # mm
        X2 = X1 + self.MosaicWindowWidthmm  # mm

        Y1 = self.g0 + self.MosaicWindowHeightmm/2  # mm
        Y2 = self.g0 - self.MosaicWindowHeightmm/2  # mm
        # plot whole mosaic
        X1_w = X1 - (self.MosaicWidthmm-self.MosaicWindowWidthmm)
        X2_w = X2

        # Y1_w = Y1 + 3.2
        # Y2_w = Y2 -1.94
        # y postion adjusted by ledge height
        toppanel_wls = list(self.hmsParamDict.MosaicFilters[0])
        wdict = self.wlParamDict[toppanel_wls[0]]
        Y1_w = Y1 + (wdict.PanelHeightmm - wdict.PanelWindowHeightmm)

        # y postion adjusted by ledge height
        bottompanel_wls = list(self.hmsParamDict.MosaicFilters[1])
        wdict = self.wlParamDict[bottompanel_wls[0]]
        Y2_w = Y2 - (wdict.PanelHeightmm - wdict.PanelWindowHeightmm)

        # limit from L-> R if looking from grating to mosaic.
        ax.set_xlim(X2_w, X1_w)
        ax.set_ylim(Y2_w, Y1_w)  # limist from Bottom to top

        ls = '--'
        lw = 1
        c = 'silver'
        ax.axvline(X1, linewidth=lw, linestyle=ls, color=c)
        ax.axhline(Y1, linewidth=lw, linestyle=ls, color=c)
        ax.axhline(Y2, linewidth=lw, linestyle=ls, color=c)

        c = 'Black'
        ls = '-'
        lw = 1.5
        ax.axhline(
            Y2_w + self.wlParamDict['6563'].PanelHeightmm, linewidth=lw, linestyle=ls, color=c)
        # test for the difference between horizontal line at gamma = 90 and the middle seam of mosiac. they should be right on top of eachother.
        # ax.axhline(self.g0,linewidth = lw, linestyle = ls,color ='orange')
        # print(Y2_w + self.wlParamDict['6563'].PanelHeightmm - self.g0)
        if Mosaic:
            # Bottom Panel of Mosaic------------------------------------
            idx = 1
            # mm, height at which the panel width label should be at.
            widthlabel = self.g0 - self.MosaicWindowHeightmm/4
            # mm, width at which the panel height label should be at.
            heightlabel = self.deg2mm(
                self.alpha, fl=self.fprime) - self.hmsParamDict.SlitA2FarEdgemm*0.95
            x_bottom = X2_w  # left edge
            y_bottom = Y2_w  # bottom edge
            wls = list(self.hmsParamDict.MosaicFilters[idx])
            wls.reverse()
            height = self.wlParamDict[wls[idx]].PanelHeightmm
            self.annotate_height(ax, height, y_bottom,
                                 heightlabel, lw, ls, c, Measurements)
            for wl in wls:
                wdict = self.wlParamDict[wl]
                lenmm = float(wdict.PanelWidthmm)
                x_bottom -= lenmm
                self.annotate_width(
                    ax, lenmm, x_bottom, widthlabel, wl, lw, ls, c, Measurements, False)

            # top Panel of Mosaic---------------------------------------
            idx = 0
            # mm, height at which the panel width label should be at.
            widthlabel = self.g0 + self.MosaicWindowHeightmm/4
            # mm, width at which the panel height label should be at.
            heightlabel = self.deg2mm(
                self.alpha, fl=self.fprime) - self.hmsParamDict.SlitA2FarEdgemm*0.95
            X2 = X1 + self.MosaicWindowWidthmm  # mm
            Y1 = self.g0 + self.MosaicWindowHeightmm/2  # mm
            x_top = X2_w  # left edge
            y_bottom = Y2_w + \
                self.wlParamDict[wls[idx]].PanelHeightmm  # top edge
            wls = list(self.hmsParamDict.MosaicFilters[idx])
            wls.reverse()
            height = self.wlParamDict[wls[idx]].PanelHeightmm
            self.annotate_height(ax, height, y_bottom,
                                 heightlabel, lw, ls, c, Measurements)
            for wl in wls:
                wdict = self.wlParamDict[wl]
                lenmm = float(wdict.PanelWidthmm)
                x_top -= lenmm
                self.annotate_width(
                    ax, lenmm, x_top, widthlabel, wl, lw, ls, c, Measurements)

        if savefig:
            plt.savefig(
                f"hms{self.hmsVersion.upper()}_Mosaic_Measurements.png")

    def _plot_on_detector(self, ax, betas, gamma, wls, Tape2Grating, fprime=None):
        if fprime is not None:
            self.fprime = fprime
        # Transform to pixel coordinates
        self.wls = wls
        X1 = self.deg2mm(self.alpha, fl=self.fprime) - \
            self.hmsParamDict.SlitA2FarEdgemm  # mm
        X2 = X1 + self.MosaicWindowWidthmm  # mm
        X1, X2 = self.mm2pix(X1, 0, self.pix), self.mm2pix(
            X2, 0, self.pix)  # pix
        self.mx, self.bx = 1, -np.min([X1, X2])  # slope,intercept
        X1, X2 = self.Linear(X1, self.mx, self.bx), self.Linear(
            X2, self.mx, self.bx)  # pix -> detectorpix

        self.g0 = self.deg2mm(self.mgammadeg, fl=self.fprime)
        Y1 = self.g0 + self.MosaicWindowHeightmm/2  # mm
        Y2 = self.g0 - self.MosaicWindowHeightmm/2  # mm
        Y1, Y2 = self.mm2pix(Y1, 1, self.pix), self.mm2pix(
            Y2, 1, self.pix)  # pix
        self.my, self.by = 1, -np.min([Y1, Y2])  # slope,intercept
        Y1, Y2 = self.Linear(Y1, self.my, self.by), self.Linear(
            Y2, self.my, self.by)  # pix -> detectorpix

        # both axis flip reflecting off of the fold mirror
        ax.set_xlim(X1, X2)
        ax.set_ylim(Y1, Y2)

        betas = self.Linear(self.mm2pix(betas, 0, self.pix), self.mx, self.bx)
        gamma = self.Linear(self.mm2pix(gamma, 1, self.pix), self.my, self.by)

        for widx, wl in enumerate(wls):
            w = str(int(wl * 10))  # nm -> A
            if w not in self.wlParamDict.keys():
                w = self.find_nearest(wl)
                wdict = self.wlParamDict.get(w)
                color = 'Black'
            else:
                wdict = self.wlParamDict.get(w)
                color = wdict.color

            # set alpha idex
            if wdict.SlitNum in [1, 3]:
                aidx = 1
            elif wdict.SlitNum in [2, 4]:
                aidx = 0
            else:
                raise ValueError("Slit Number must be 1, 2, 3, or 4.")

            # set gamma filter accoding to SlitNum
            # g0 = self.Linear(self.deg2pix(self.mgammadeg,1,self.f,self.pix),mx,bx)
            g0 = self.deg2mm(self.mgammadeg, self.fprime)  # deg -> mm
            g0 = self.mm2pix(g0, 1)  # mm -> pix
            self.g0 = self.Linear(g0, self.my, self.by)  # pix -> dpix
            if wdict.SlitNum in [1, 2]:
                gmask = (gamma >= self.g0)
            else:
                gmask = (gamma <= self.g0)

            gamma_plot = gamma[gmask]
            for midx, m in enumerate(self.orders):
                beta_plot = betas[widx][aidx][midx][gmask]
                ax.plot(beta_plot, gamma_plot, self.ls,
                        markersize=self.ms, color=color, linewidth=self.lw)
                idx = int(len(beta_plot)/4)
                xy = (beta_plot[idx], gamma_plot[idx+10])
                xytext = (beta_plot[idx]+1, gamma_plot[idx+5])
                slitnum = wdict.SlitNum
                if self.slitwidth is not None:
                    res = self.resolving_power(wl, m)*1e-3
                    beta = self.get_beta(wl, m) - self.alpha
                    ax.annotate(f'[{wl},{slitnum},{m},{res:.0f}k,{round(beta, 0):.0f}°]',
                                xy, xytext, rotation=270)
                else:
                    ax.annotate(f'[{wl},{slitnum},{m}]',
                                xy, xytext, rotation=270)
        # Additional plotting and detector setup based on Tape2Grating

    def clip_image_percentiles(self, image: np.ndarray, lower_percentile: float = 0, upper_percentile: float = 99) -> np.ndarray:
        lower_value = np.percentile(image, lower_percentile)
        upper_value = np.percentile(image, upper_percentile)

        clipped_image = np.clip(image, lower_value, upper_value)
        return clipped_image

    def overlay_on_Image(self, impath: str, wls: Iterable = [486.1, 427.8, 557.7, 630.0, 656.3, 777.4], vmin: float = 0, vmax: float = np.nan, fprime=None):
        self.wls = wls
        if '.png' in impath:
            data = plt.imread(impath)
        elif '.fit' in impath:
            with pf.open(impath) as hdul:
                data = hdul[1].data.astype(np.float64)
        else:
            raise ValueError("Unsupported file format")
        # data = self.clip_image_percentiles(data)
        self.pix = data.shape[1]
        fig = self.plot_spectral_lines(
            ImageAt='detector', Tape2Grating=True, wls=wls)
        data = equalize_hist(data)
        cmap = 'gray'
        if np.isnan(vmax):
            plt.imshow(data, cmap=cmap)
        else:
            plt.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap)
        plt.colorbar()

    def line_spacing(self, wls: Iterable, alphas: Iterable = None, cidx: Iterable = [0], observed: bool = False):
        self.wls = wls
        betas = []
        if observed:
            if alphas is None:
                raise ValueError(
                    "Alphas must be provided for observed line spacing.")
            elif len(alphas) != len(wls):
                raise ValueError(
                    "An Alpha must be provided for each wl.This is the alpha where the model at that wl matches observation.")
            paired = sorted(zip(wls, alphas), reverse=False,
                            key=lambda x: x[0])
            wls, alphas = zip(*paired)
        else:
            wls = sorted(wls, reverse=False)

        for widx, wl in enumerate(wls):
            w = str(int(wl * 10))  # nm -> A
            if w not in self.wlParamDict.keys():
                w = self.find_nearest(wl)
                wdict = self.wlParamDict.get(w)
                color = 'Black'
            else:
                wdict = self.wlParamDict.get(w)
                color = wdict.color

            # Determine alpha based on the observed parameter
            if observed:
                if alphas is None:
                    raise ValueError(
                        "Alphas must be provided for observed line spacing.")
                rel_alphas = self.relative_alphas(alphas[widx])
            else:
                rel_alphas = self.alphas

            if wdict.SlitNum in [1, 3]:
                alpha = rel_alphas[1]
            else:
                alpha = rel_alphas[1]

            morder = int(wdict.DiffractionOrder)
            beta_wl = []
            for i in cidx:
                print(alpha, self.gamma[i], morder, self.sigma, wl)
                beta_wl.append(self.single_defraction(
                    alpha, self.gamma[i], morder, self.sigma, wl))
            betas.append(np.array(beta_wl))

        return np.array(self.deg2mm((betas[0] - betas[1]), fl=self.f))

    def fit_fprime(self, wls: str, alp: Iterable):
        """_summary_

        Args:
            element (str): Oxygen lines or Hydrogen lines
            alp (Iterable): Alphas at which the model matched observed img order of 'oxygen': wl= [630.0,557.7]
            and hydrogen':wl = [656.3, 486.1]
        """
        plt.figure()
        plt.title(f"HiT&MIS {self.hmsVersion}")
        cidx = np.arange(50, 90, 1)  # make it 101 after
        # True line spacing
        ds = self.line_spacing(wls, cidx=cidx, observed=False)
        print(ds)
        # Observed Line Spacing
        ds_prime = self.line_spacing(wls, alp, cidx, True)

        def fitfunc(ds: float, fprime: float): return (fprime*ds/400)

        popt, pcov = curve_fit(fitfunc, ds, ds_prime)

        plt.scatter(ds, ds_prime, marker='x', color='red', s=5)
        plt.plot(ds, fitfunc(ds, popt[0]), label=f"Slope (f') = {popt[0]:.3f}")
        plt.xlabel("True line Spacing, ds [mm]")
        plt.ylabel("observed line Spacing, ds' [mm]")
        plt.legend()


# %%
# predictor = HMS_ImagePredictor(hmsVersion='bo',mgammadeg=90,alpha=66.45)

# %%
# img = predictor.plot_spectral_lines('MosaicWindow',True,wls = np.arange(774.1,788.4), mosaic=True,measurement=False)

# img = predictor.plot_spectral_lines('Mosaicwindow',True,wls = [557.7, 630.0,427.8, 784.1,777.4,486.1,656.3,656.8,481,644,786.0,782.1,780.8,652.2,654.4,653.3, 774.4], mosaic=True,measurement=True)

# %%
# img = predictor.plot_spectral_lines('Mosaicwindow',True,wls = [557.7, 630.0,427.8, 784.1,777.4,486.1,656.3])
# %%
# if __name__ == "__main__":
    # predictor = HMS_ImagePredictor(hmsVersion='b',mgammadeg=90,alpha=85)
    # imgdir = 'Images/hmsA_testimg/ccdi_20240811_003321.fits'
    # predictor.overlay_on_Image(imgdir)
    # predictor.plot_spectral_lines('', True,False,wls=[486.1, 656.3,844.6,557.7])
    # predictor.plot_spectral_lines('detector',True)
# %%


# %%
