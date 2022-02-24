#!/usr/bin/env python
# coding: utf-8

import pynbody
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.stats import gaussian_fwhm_to_sigma


class Aperture(pynbody.filt.Filter):

    """
    Return particles that are within an x-y aperture with `radius', centered on the point `cen`.

    Inputs:
    -------
    *radius* : extent of the aperture. Can be a number or a string specifying the units.
    *cen* : center of the aperture, in (x, y). default = (0,0)
    """

    def __init__(self, radius, cen=(0, 0)):
        self._descriptor = "aperture"
        self.cen = np.asarray(cen)
        if self.cen.shape != (2,):
            raise ValueError("Centre must be length 2 array")

        if isinstance(radius, str):
            radius = pynbody.units.Unit(radius)

        self.radius = radius

    def __call__(self, sim):
        radius = self.radius
        wrap = -1.0

        with sim.immediate_mode:
            pos = sim["pos"]

        if pynbody.units.is_unit_like(radius):
            radius = float(radius.in_units(pos.units, **pos.conversion_context()))

        if "boxsize" in sim.properties:
            wrap = sim.properties["boxsize"]

        if pynbody.units.is_unit_like(wrap):
            wrap = float(wrap.in_units(pos.units, **pos.conversion_context()))

        cen = self.cen
        if pynbody.units.has_units(cen):
            cen = cen.in_units(pos.units)

        distance = ((sim["pos"][:, :2] - self.cen) ** 2).sum(axis=1)
        return distance < radius ** 2

    def __repr__(self):
        if pynbody.units.is_unit(self.radius):

            return "Aperture('%s', %s)" % (str(self.radius), repr(self.cen))
        else:
            return "Aperture(%.2e, %s)" % (self.radius, repr(self.cen))


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
        halo_max_size=None,
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

        if halo_max_size is None:
            self.halo_max_size = aperture_kpc
        else:
            self.halo_max_size = halo_max_size

        self.kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).to("kpc arcsec**-1").value
        self.npixels = self.calc_npixels(self.image_width_kpc, self.pixel_scale_arcsec)
        self.kpc_per_pixel = self.image_width_kpc / self.npixels

        self.particles = self.restrict_halo(halo)
        self.mask = self.mask_aperture()
        self.raw = self.generate_los_map()

        self.data = np.asarray(self.raw)

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

    def restrict_halo(self, halo):
        """Restricts halo to a 3d radius from the center.

        Args:
            halo (pynbody.snapshot.SimSnap): Pynbody halo object.

        Returns:
            pynbody.snapshot.SimSnap: Cropped halo object.
        """
        inside_aperture = Aperture(self.halo_max_size, cen=(0, 0))
        return halo[inside_aperture]

    def mask_aperture(self):
        """Mask data outside the aperture.

        Returns:
            array-like: Boolean mask marking pixels within the aperture.
        """
        if self.aperture_radius is not None:
            w, h = self.npixels, self.npixels
            center = (int(w / 2), int(h / 2))
            radius = (
                self.aperture_radius / self.kpc_per_arcsec / self.pixel_scale_arcsec
            )

            Y, X = np.ogrid[:h, :w]
            dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

            mask = dist_from_center <= radius
        else:
            mask = np.ones((self.npixels, self.npixels))
        return mask

    def generate_los_map(self):
        """Creates line-of-sight velocity map from self.particles.

        Returns:
            array-like: pixel-by-pixel velocity map
        """
        im = pynbody.plot.sph.image(
            self.particles,
            qty="map_qty",
            log=False,
            units=self.particles["map_qty"].units,
            width=self.image_width_kpc,
            resolution=self.npixels,
            noplot=True,
            show_cbar=False,
        )
        im2 = pynbody.plot.sph.image(
            self.particles,
            qty="map_weights",
            log=False,
            units=self.particles["map_weights"].units,
            width=self.image_width_kpc,
            resolution=self.npixels,
            noplot=True,
            show_cbar=False,
        )
        im = im / im2
        return im.in_units("km s**-1")

    def create_masked_image(self, im, mask):
        return np.ma.masked_array(im[::-1, :], mask=mask)

    def convolve_fwhm(self):
        """Convolve a observational PSF with full width half max of self.fwhm_arcsec.

        Returns:
            array-like: Masked numpy array containing the image, as well as the mask outlining the aperture.
        """
        sigma_arcsec = self.fwhm_arcsec * gaussian_fwhm_to_sigma
        sigma_pixels = sigma_arcsec / self.pixel_scale_arcsec
        im = gaussian_filter(self.data.data, sigma=sigma_pixels)
        return im