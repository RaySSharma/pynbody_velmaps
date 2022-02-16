#!/usr/bin/env python
# coding: utf-8

import pynbody
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.stats import gaussian_fwhm_to_sigma


class VelocityMap:
    def __init__(
        self,
        halo,
        z,
        cosmo,
        image_width_kpc,
        aperture_kpc,
        pixel_scale_arcsec=0.05,
        fwhm_arcsec=None,
    ):
        """Generates a velocity map for the particles in `halo`, following the prescription in Koudmani et al (2021). Defaults chosen to match the MANGA survey.

        Args:
            halo (pynbody.snapshot.SimSnap): Pynbody halo object, already centered + aligned.
            z (float): Redshift of halo.
            cosmo (astropy.cosmology.core.FlatLambdaCDM): Astropy cosmology object for simulation cosmology.
            image_width_kpc (float): Width of final image in kpc.
            aperture_kpc (float): Aperture radius in kpc within which to generate map. Typically 1.5 * effective radius.
            pixel_scale_arcsec (float): Desired pixel scale of final image, in arcseconds. Defaults to 0.05".
            fwhm_arcsec (float): Desired FWHM of final image, in arcseconds. Skip FWHM by setting to `None`. Defaults to None.
        """

        self.pixel_scale_arcsec = pixel_scale_arcsec
        self.image_width_kpc = image_width_kpc
        self.fwhm_arcsec = fwhm_arcsec
        self.aperture_radius = aperture_kpc

        self.kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).to("kpc arcsec**-1").value
        self.npixels = self.calc_npixels(self.image_width_kpc, self.pixel_scale_arcsec)
        self.kpc_per_pixel = self.image_width_kpc / self.npixels

        self.particles = self.restrict_to_aperture(halo)

        self.data = self.generate_los_map()
        self.data.mask = self.mask_aperture()
        if self.fwhm_arcsec is not None:
            self.data = self.convolve_fwhm()

    def calc_npixels(self, image_width_kpc, pixel_scale_arcsec):
        """Converts image width and pixel scale to an image width in pixels.

        Args:
            image_width_kpc (float): Desired image width in kpc units.
            pixel_scale_arcsec (float): Desired pixel scale in arcsec.

        Returns:
            int: Image width in pixels.
        """
        image_width_arcsec = image_width_kpc / self.kpc_per_arcsec
        npixels = int(image_width_arcsec / pixel_scale_arcsec)
        return npixels

    def restrict_to_aperture(self, halo):
        """Restricts halo to a 3d radius from the center.

        Args:
            halo (pynbody.snapshot.SimSnap): Pynbody halo object.

        Returns:
            pynbody.snapshot.SimSnap: Cropped halo object.
        """
        if self.aperture_radius is not None:
            inside_aperture = pynbody.filt.Sphere(self.aperture_radius, cen=(0,0,0))
            return halo[inside_aperture]
        else:
            return halo

    def mask_aperture(self):
        """Mask data outside the aperture.

        Returns:
            array-like: Boolean mask marking pixels within the aperture.
        """
        w, h = self.data.shape
        center = (int(w / 2), int(h / 2))
        radius = self.aperture_radius / self.kpc_per_arcsec / self.pixel_scale_arcsec

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

        mask = dist_from_center <= radius
        return mask

    def generate_los_map(self):
        """Creates line-of-sight velocity map from self.particles.

        Returns:
            array-like: pixel-by-pixel velocity map
        """
        im = pynbody.plot.sph.image(
            self.particles,
            qty="vz",
            av_z="mass",
            units="km s**-1",
            log=False,
            width=self.image_width_kpc,
            resolution=self.npixels,
            noplot=True,
            show_cbar=False,
        )[::-1, :]
        return np.ma.masked_array(im)

    def convolve_fwhm(self):
        """Convolve a observational PSF with full width half max of self.fwhm_arcsec.

        Returns:
            array-like: Masked numpy array containing the image, as well as the mask outlining the aperture.
        """
        sigma_arcsec = self.fwhm_arcsec * gaussian_fwhm_to_sigma
        sigma_pixels = sigma_arcsec / self.pixel_scale_arcsec
        im = gaussian_filter(self.data.data * self.data.mask, sigma=sigma_pixels)
        return np.ma.masked_array(im, mask=~self.data.mask)
