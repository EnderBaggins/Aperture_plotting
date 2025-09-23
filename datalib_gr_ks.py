#!/usr/bin/env python

import h5py
import numpy as np
import toml
from pathlib import Path
import os
import re
from datalib import Data, flag_to_species
from datalib_logsph import DataSph

def Sigma(r,th,a):
  return r**2+a**2*np.cos(th)**2
def Delta(r,a):
  return r**2-2.0*r+a**2
def AA(r,th,a):
  return (r**2+a**2)**2-Delta(r,a)*a**2*np.sin(th)**2
def alpha(r,th,a):
  return 1.0/np.sqrt(1.0+2.0*r/Sigma(r,th,a))
def beta1u(r,th,a):
  return 2.0*r/(Sigma(r,th,a)+2.0*r)
def gmsqrt(r,th,a):
  return Sigma(r,th,a)*np.sin(th)*np.sqrt(1.0+2.0*r/Sigma(r,th,a))
def gsqrt(r,th,a):
  return Sigma(r,th,a)*np.sin(th)
def beta1d(r,th,a):
  return 2.0*r/Sigma(r,th,a)
def beta3d(r,th,a):
  return -2.0*a*r*np.sin(th)**2/Sigma(r,th,a)

def gd00(r,th,a):
  return -(1-2.0*r/Sigma(r,th,a))
def gd01(r,th,a):
  return 2.0*r/Sigma(r,th,a)
def gd03(r,th,a):
  return -2.0*a*r*np.sin(th)**2/Sigma(r,th,a)
def gd11(r,th,a):
  return 1.0+2.0*r/Sigma(r,th,a)
def gd13(r,th,a):
  return -a*(1.0+2.0*r/Sigma(r,th,a))*np.sin(th)**2
def gd22(r,th,a):
  return Sigma(r,th,a)
def gd33(r,th,a):
  return AA(r,th,a)*np.sin(th)**2/Sigma(r,th,a)

# These are gamma metric functions
def gmd11(r,th,a):
  return gd11(r,th,a)
def gmd13(r,th,a):
  return gd13(r,th,a)
def gmd22(r,th,a):
  return gd22(r,th,a)
def gmd33(r,th,a):
  return gd33(r,th,a)

def gmu11(r,th,a):
  return (a**2+r**2)/Sigma(r,th,a)-2.0*r/(Sigma(r,th,a)+2.0*r)
def gmu13(r,th,a):
  return a/Sigma(r,th,a)
def gmu22(r,th,a):
  return 1.0/Sigma(r,th,a)
def gmu33(r,th,a):
  return 1.0/Sigma(r,th,a)/np.sin(th)**2

def gu00(r,th,a):
    return -(1.0+2.0*r/Sigma(r,th,a))
def gu01(r,th,a):
    return 2.0*r/Sigma(r,th,a)
def gu11(r,th,a):
    return Delta(r,a)/Sigma(r,th,a)
def gu13(r,th,a):
    return a/Sigma(r,th,a)
def gu22(r,th,a):
    return 1.0/Sigma(r,th,a)
def gu33(r,th,a):
    return 1.0/Sigma(r,th,a)/np.sin(th)**2

# Outer event horizon
def rs_o(a):
    return 1.0+np.sqrt(1.0-a**2)
# Inner event horizon
def rs_i(a):
    return 1.0-np.sqrt(1.0-a**2)
def t_BL(t_KS,r,a):
    return t_KS-1.0/np.sqrt(1.0-a**2)*(rs_o(a)*np.log(np.fabs(r/rs_o(a)-1.0))
                                       -rs_i(a)*np.log(np.fabs(r/rs_i(a)-1.0)))

# Initialize the 3D and 4D Levi-Civita symbols
levi_civita3 = np.zeros((3, 3, 3))
levi_civita4 = np.zeros((4, 4, 4, 4))

for i in range(3):
    for j in range(3):
        for k in range(3):
            if i == j or i == k or j == k:
                levi_civita3[i, j, k] = 0
            elif (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
                levi_civita3[i, j, k] = 1
            else:
                levi_civita3[i, j, k] = -1

for i in range(4):
    for j in range(4):
        for k in range(4):
            for l in range(4):
                if i == j or i == k or i == l or j == k or j == l or k == l:
                    levi_civita4[i, j, k, l] = 0
                elif (i, j, k, l) in [(0, 1, 2, 3), (0, 3, 1, 2), (0, 2, 3, 1), (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0), (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1), (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0)]:
                    levi_civita4[i, j, k, l] = 1
                else:
                    levi_civita4[i, j, k, l] = -1

def reduce_by_tiles_reduceat(arr, tile_size):
    N = arr.shape[0]
    # Sum along rows first
    row_sums = np.add.reduceat(arr, np.arange(0, N, tile_size), axis=0)
    # Then sum along columns
    return np.add.reduceat(row_sums, np.arange(0, N, tile_size), axis=1)


class DataKerrSchild(DataSph):
  _mesh_loaded = False
  _metric_ready = False

  def __init__(self, path, tile_size=4):
    super().__init__(path)
    self.a = self._conf["bh_spin"]
    self.rH = rs_o(self.a)
    self.extra_fld_keys = ["fluxB", "Dd1", "Dd2", "Dd3", "D", "Bd1", "Bd2", "Bd3", "B", 
                           "Ed1", "Ed2", "Ed3", "Hd1", "Hd2", "Hd3",
                           "sigma", "flux_upper", "flux_lower", "n_proper", 
                           "fluid_u_upper", "fluid_u_lower", "fluid_b_upper", 
                           "stress_e", "stress_p", "frf_transform", "frf_transform_inv", "frf_T_munu",
                           "plasma_temp", "pressure_para", "pressure_perp", "plasma_beta", "frf_B",
                           "stress_reduced", "flux_upper_reduced", "flux_lower_reduced", "n_proper_reduced",
                           "fluid_u_upper_reduced", "fluid_u_lower_reduced", "fluid_b_upper_reduced"]
    self.tile_size = tile_size
    self.reduced_shape = (self.x1.shape[0] // tile_size, self.x1.shape[1] // tile_size)
    self._rs_reduced = np.exp(
      np.linspace(
        0,
        self._conf["size"][0],
        self._conf["N"][0]
        // self._conf["downsample"] // self.tile_size,
      ) + self._conf["lower"][0]
    )
    self._ths_reduced = np.linspace(
          0,
          self._conf["size"][1],
          self._conf["N"][1] // self._conf["downsample"] // self.tile_size,
        ) + self._conf["lower"][1]

    self._rv_reduced, self._thv_reduced = np.meshgrid(self._rs_reduced, self._ths_reduced)
    self.x1_reduced = self._rv_reduced * np.sin(self._thv_reduced)
    self.x2_reduced = self._rv_reduced * np.cos(self._thv_reduced)
    self.compute_metrics()
    self.reload()


  def save_diagnostics(self, name, array):
    path = os.path.join(self._path, f"diagnostics.{self._current_fld_step:05d}.h5")
    with h5py.File(path, "a") as f:
      if not name in f:
        f[name] = array

  def compute_metrics(self):
    if self._metric_ready:
      return
    self.g_upper = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))

    self.g_upper[:, :, 0, 0] = gu00(self._rv, self._thetav, self.a)
    self.g_upper[:, :, 0, 1] = self.g_upper[:, :, 1, 0] = gu01(self._rv, self._thetav, self.a)
    self.g_upper[:, :, 1, 1] = gu11(self._rv, self._thetav, self.a)
    self.g_upper[:, :, 2, 2] = gu22(self._rv, self._thetav, self.a)
    self.g_upper[:, :, 3, 3] = gu33(self._rv, self._thetav, self.a)
    self.g_upper[:, :, 1, 3] = self.g_upper[:, :, 3, 1] = gu13(self._rv, self._thetav, self.a)

    self.g_lower = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))

    self.g_lower[:, :, 0, 0] = gd00(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 0, 1] = self.g_lower[:, :, 1, 0] = gd01(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 0, 3] = self.g_lower[:, :, 3, 0] = gd03(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 1, 1] = gd11(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 2, 2] = gd22(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 3, 3] = gd33(self._rv, self._thetav, self.a)
    self.g_lower[:, :, 1, 3] = self.g_lower[:, :, 3, 1] = gd13(self._rv, self._thetav, self.a)

    self.g_upper_reduced = np.zeros((self.reduced_shape[0], self.reduced_shape[1], 4, 4))
    self.g_lower_reduced = np.zeros((self.reduced_shape[0], self.reduced_shape[1], 4, 4))

    self.g_upper_reduced[:, :, 0, 0] = gu00(self._rv_reduced, self._thv_reduced, self.a)
    self.g_upper_reduced[:, :, 0, 1] = self.g_upper_reduced[:, :, 1, 0] = gu01(self._rv_reduced, self._thv_reduced, self.a)
    self.g_upper_reduced[:, :, 1, 1] = gu11(self._rv_reduced, self._thv_reduced, self.a)
    self.g_upper_reduced[:, :, 2, 2] = gu22(self._rv_reduced, self._thv_reduced, self.a)
    self.g_upper_reduced[:, :, 3, 3] = gu33(self._rv_reduced, self._thv_reduced, self.a)
    self.g_upper_reduced[:, :, 1, 3] = self.g_upper_reduced[:, :, 3, 1] = gu13(self._rv_reduced, self._thv_reduced, self.a)

    self.g_lower_reduced[:, :, 0, 0] = gd00(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 0, 1] = self.g_lower_reduced[:, :, 1, 0] = gd01(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 0, 3] = self.g_lower_reduced[:, :, 3, 0] = gd03(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 1, 1] = gd11(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 2, 2] = gd22(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 3, 3] = gd33(self._rv_reduced, self._thv_reduced, self.a)
    self.g_lower_reduced[:, :, 1, 3] = self.g_lower_reduced[:, :, 3, 1] = gd13(self._rv_reduced, self._thv_reduced, self.a)

    self.sqrt_g = gsqrt(self._rv, self._thetav, self.a)
    self.sqrt_g_reduced = gsqrt(self._rv_reduced, self._thv_reduced, self.a)
    self._metric_ready = True

  def raise_4d_vec(self, vec_lower):
    if not self._metric_ready:
      self.compute_metrics()
    return np.einsum("ijab,ijb->ija", self.g_upper, vec_lower)

  def lower_4d_vec(self, vec_upper):
    if not self._metric_ready:
      self.compute_metrics()
    return np.einsum("ijab,ijb->ija", self.g_lower, vec_upper)
  
  def dot_4d(self, vec_upper, vec_lower):
    return np.einsum("ija, ija->ij", vec_upper, vec_lower)

  def inner_product_4d_contravariant(self, vec1, vec2):
    if not self._metric_ready:
      self.compute_metrics()
    return np.einsum("ijab,ija,ijb->ij", self.g_lower, vec1, vec2)

  def inner_product_4d_covariant(self, vec1, vec2):
    if not self._metric_ready:
      self.compute_metrics()
    return np.einsum("ijab,ija,ijb->ij", self.g_upper, vec1, vec2)

  def _load_fld_quantity(self, key):
    path = os.path.join(self._path, f"fld.{self._current_fld_step:05d}.h5")
    if key == "fluxB":
      self._load_sph_mesh()
      self.__dict__[key] = np.cumsum(self.B1 * gmsqrt(self._rv, self._thetav, self.a) * self._dtheta, axis=0)
    # Lower components of D
    elif key == "Dd1":
      self.__dict__[key] = gmd11(self._rv, self._thetav, self.a) * self.E1 + gmd13(self._rv, self._thetav, self.a) * self.E3
    elif key == "Dd2":
      self.__dict__[key] = gmd22(self._rv, self._thetav, self.a) * self.E2
    elif key == "Dd3":
      self.__dict__[key] = gmd13(self._rv, self._thetav, self.a) * self.E1 + gmd33(self._rv, self._thetav, self.a) * self.E3
    elif key == "D":
      self.__dict__[key] = np.sqrt(self.E1 * self.Dd1 + self.E2 * self.Dd2 + self.E3 * self.Dd3)
    # Lower components of B
    elif key == "Bd1":
      self.__dict__[key] = gmd11(self._rv, self._thetav, self.a) * self.B1 + gmd13(self._rv, self._thetav, self.a) * self.B3
    elif key == "Bd2":
      self.__dict__[key] = gmd22(self._rv, self._thetav, self.a) * self.B2
    elif key == "Bd3":
      self.__dict__[key] = gmd13(self._rv, self._thetav, self.a) * self.B1 + gmd33(self._rv, self._thetav, self.a) * self.B3
    elif key == "B":
      self.__dict__[key] = np.sqrt(self.B1 * self.Bd1 + self.B2 * self.Bd2 + self.B3 * self.Bd3)
    # Auxiliary field E lower
    elif key == "Ed1":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Dd1
    elif key == "Ed2":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Dd2 - gmsqrt(self._rv, self._thetav, self.a) * beta1u(self._rv, self._thetav, self.a) * self.B3
    elif key == "Ed3":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Dd3 + gmsqrt(self._rv, self._thetav, self.a) * beta1u(self._rv, self._thetav, self.a) * self.B2
    # Auxiliary field H lower
    elif key == "Hd1":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Bd1
    elif key == "Hd2":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Bd2 + gmsqrt(self._rv, self._thetav, self.a) * beta1u(self._rv, self._thetav, self.a) * self.E3
    elif key == "Hd3":
      self.__dict__[key] = alpha(self._rv, self._thetav, self.a) * self.Bd3 - gmsqrt(self._rv, self._thetav, self.a) * beta1u(self._rv, self._thetav, self.a) * self.E2
    elif key == "sigma": # this is the cold sigma, computed in the FRF
      # self.__dict__[key] = self.B**2 / (self.Rho_p - self.Rho_e + 1e-6)
      b2 = self.inner_product_4d_contravariant(self.frf_B, self.frf_B)
      self.__dict__[key] = b2 / self.n_proper
    elif key == "flux_upper":
      self.__dict__[key] = self.raise_4d_vec(self.flux_lower)
    elif key == "flux_lower":
      self.__dict__[key] = np.stack([self.num_e + self.num_p, self.flux_e1 + self.flux_p1,
                                     self.flux_e2 + self.flux_p2, self.flux_e3 + self.flux_p3], axis=-1)
    elif key == "n_proper":
      # self.__dict__[key] = np.sqrt(np.abs(inner_product_4d_covariant(self.flux_lower, self.flux_lower, self._rv, self._thetav, self.a)))
      # Technically we know that flux is timelike and we should just insert a negative sign. However abs is more robust
      self.__dict__[key] = np.sqrt(np.abs(self.dot_4d(self.flux_upper, self.flux_lower)))
    elif key == "fluid_u_upper":
      indices = np.where(self.n_proper > 0)
      u_upper = np.copy(self.flux_upper)
      u_upper[indices] = self.flux_upper[indices] / self.n_proper[indices][..., np.newaxis]
      self.__dict__[key] = u_upper
    elif key == "fluid_u_lower":
      indices = np.where(self.n_proper > 0)
      u_lower = np.copy(self.flux_lower)
      u_lower[indices] = self.flux_lower[indices] / self.n_proper[indices][..., np.newaxis]
      self.__dict__[key] = u_lower
    elif key == "frf_B":
      # compute b vector in the fluid rest frame
      B = np.stack([np.zeros_like(self.B1), self.B1, self.B2, self.B3], axis=-1)
      E = np.stack([np.zeros_like(self.Ed1), self.Ed1, self.Ed2, self.Ed3], axis=-1)
      u_lower = self.fluid_u_lower
      alpha_val = alpha(self._rv, self._thetav, self.a)
      sqrt_gm = gmsqrt(self._rv, self._thetav, self.a)
      b0 = (u_lower[...,1] * B[...,1] + u_lower[...,2] * B[...,2] + u_lower[...,3] * B[...,3]) / alpha_val
      b1 = -u_lower[...,0] * B[...,1] / alpha_val - (u_lower[...,2] * E[...,3] - u_lower[...,3] * E[...,2]) / alpha_val / sqrt_gm
      b2 = -u_lower[...,0] * B[...,2] / alpha_val - (u_lower[...,3] * E[...,1] - u_lower[...,1] * E[...,3]) / alpha_val / sqrt_gm
      b3 = -u_lower[...,0] * B[...,3] / alpha_val - (u_lower[...,1] * E[...,2] - u_lower[...,2] * E[...,1]) / alpha_val / sqrt_gm
      # b_upper = np.array([b0, b1, b2, b3])
      self.__dict__[key] = np.stack([b0, b1, b2, b3], axis=-1)
    elif key == "fluid_b_upper":
      bnorm = np.sqrt(self.inner_product_4d_contravariant(self.frf_B, self.frf_B))
      # bnorm = np.sqrt(inner_product_4d_contravariant(self.frf_B, self.frf_B, self._rv, self._thetav, self.a))
      self.__dict__[key] = self.frf_B / bnorm[..., np.newaxis]
    elif key == "stress_e":
      stress_e = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))
      stress_e[:, :, 0, 0] = self.stress_e00
      stress_e[:, :, 0, 1] = stress_e[:, :, 1, 0] = self.stress_e01
      stress_e[:, :, 0, 2] = stress_e[:, :, 2, 0] = self.stress_e02
      stress_e[:, :, 0, 3] = stress_e[:, :, 3, 0] = self.stress_e03
      stress_e[:, :, 1, 1] = self.stress_e11
      stress_e[:, :, 1, 2] = stress_e[:, :, 2, 1] = self.stress_e12
      stress_e[:, :, 1, 3] = stress_e[:, :, 3, 1] = self.stress_e13
      stress_e[:, :, 2, 2] = self.stress_e22
      stress_e[:, :, 2, 3] = stress_e[:, :, 3, 2] = self.stress_e23
      stress_e[:, :, 3, 3] = self.stress_e33
      # stress_e = np.stack([[self.stress_e00, self.stress_e01, self.stress_e02, self.stress_e03],
      #                      [self.stress_e01, self.stress_e11, self.stress_e12, self.stress_e13],
      #                      [self.stress_e02, self.stress_e12, self.stress_e22, self.stress_e23],
      #                      [self.stress_e03, self.stress_e13, self.stress_e23, self.stress_e33]
      #                     ], axis=-1)  # shape: (4, 4, grid_y, grid_x)
      # stress_e = np.moveaxis(stress_e, 0, -2)  # shape: (grid_y, grid_x, 4, 4)
      self.__dict__[key] = stress_e
    elif key == "stress_p":
      stress_p = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))
      stress_p[:, :, 0, 0] = self.stress_p00
      stress_p[:, :, 0, 1] = stress_p[:, :, 1, 0] = self.stress_p01
      stress_p[:, :, 0, 2] = stress_p[:, :, 2, 0] = self.stress_p02
      stress_p[:, :, 0, 3] = stress_p[:, :, 3, 0] = self.stress_p03
      stress_p[:, :, 1, 1] = self.stress_p11
      stress_p[:, :, 1, 2] = stress_p[:, :, 2, 1] = self.stress_p12
      stress_p[:, :, 1, 3] = stress_p[:, :, 3, 1] = self.stress_p13
      stress_p[:, :, 2, 2] = self.stress_p22
      stress_p[:, :, 2, 3] = stress_p[:, :, 3, 2] = self.stress_p23
      stress_p[:, :, 3, 3] = self.stress_p33
      # stress_p = np.stack([[self.stress_p00, self.stress_p01, self.stress_p02, self.stress_p03],
      #                      [self.stress_p01, self.stress_p11, self.stress_p12, self.stress_p13],
      #                      [self.stress_p02, self.stress_p12, self.stress_p22, self.stress_p23],
      #                      [self.stress_p03, self.stress_p13, self.stress_p23, self.stress_p33]
      #                     ], axis=-1)  # shape: (4, 4, grid_y, grid_x)
      # stress_p = np.moveaxis(stress_p, 0, -2)  # shape: (grid_y, grid_x, 4, 4)
      self.__dict__[key] = stress_p
    elif key == "stress_reduced":
      stress = (self.stress_e + self.stress_p) * self.sqrt_g[..., np.newaxis, np.newaxis]
      stress_reduced = reduce_by_tiles_reduceat(stress, self.tile_size) / self.sqrt_g_reduced[..., np.newaxis, np.newaxis] / self.tile_size**2
      self.__dict__[key] = stress_reduced
    elif key == "flux_lower_reduced":
      flux = self.flux_lower * self.sqrt_g[..., np.newaxis]
      flux_reduced = reduce_by_tiles_reduceat(flux, self.tile_size) / self.sqrt_g_reduced[..., np.newaxis] / self.tile_size**2
      self.__dict__[key] = flux_reduced
    elif key == "flux_upper_reduced":
      self.__dict__[key] = np.einsum('ijab,ijb->ija', self.g_upper_reduced, self.flux_lower_reduced)
    elif key == "n_proper_reduced":
      self.__dict__[key] = np.sqrt(np.abs(self.dot_4d(self.flux_upper_reduced, self.flux_lower_reduced)))
    elif key == "frf_transform":
      t_vec = np.array([0, 1, 0, 0])
      e2_vec = np.einsum('abcd,ijb,ijc,d->ija', levi_civita4, self.fluid_u_upper, self.fluid_b_upper, t_vec)
      e2_vec = raise_4d_vec(e2_vec, self._rv, self._thetav, self.a)
      e2_vec /= np.sqrt(np.abs(self.inner_product_4d_contravariant(e2_vec, e2_vec)))[..., np.newaxis]
      # e2_vec /= np.sqrt(np.abs(inner_product_4d_contravariant(e2_vec, e2_vec, self._rv, self._thetav, self.a)))[..., np.newaxis]
      e3_vec = np.einsum('abcd,ijb,ijc,ijd->ija', levi_civita4, self.fluid_u_upper, self.fluid_b_upper, e2_vec)
      e3_vec = raise_4d_vec(e3_vec, self._rv, self._thetav, self.a)
      e3_vec /= np.sqrt(np.abs(self.inner_product_4d_contravariant(e3_vec, e3_vec)))[..., np.newaxis]
      # e3_vec /= np.sqrt(np.abs(inner_product_4d_contravariant(e3_vec, e3_vec, self._rv, self._thetav, self.a)))[..., np.newaxis]
      Rs = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))
      # This is dx/dx_hat
      Rs[:, :, 0, 0] = self.fluid_u_upper[..., 0]
      Rs[:, :, 0, 1] = e2_vec[..., 0]
      Rs[:, :, 0, 2] = e3_vec[..., 0]
      Rs[:, :, 0, 3] = self.fluid_b_upper[..., 0]
      Rs[:, :, 1, 0] = self.fluid_u_upper[..., 1]
      Rs[:, :, 1, 1] = e2_vec[..., 1]
      Rs[:, :, 1, 2] = e3_vec[..., 1]
      Rs[:, :, 1, 3] = self.fluid_b_upper[..., 1]
      Rs[:, :, 2, 0] = self.fluid_u_upper[..., 2]
      Rs[:, :, 2, 1] = e2_vec[..., 2]
      Rs[:, :, 2, 2] = e3_vec[..., 2]
      Rs[:, :, 2, 3] = self.fluid_b_upper[..., 2]
      Rs[:, :, 3, 0] = self.fluid_u_upper[..., 3]
      Rs[:, :, 3, 1] = e2_vec[..., 3]
      Rs[:, :, 3, 2] = e3_vec[..., 3]
      Rs[:, :, 3, 3] = self.fluid_b_upper[..., 3]
      self.__dict__[key] = Rs
    elif key == "frf_transform_inv":
      fluid_b_lower = self.lower_4d_vec(self.fluid_b_upper)
      e2_vec_lower = self.lower_4d_vec(self.frf_transform[:, :, :, 1])
      e3_vec_lower = self.lower_4d_vec(self.frf_transform[:, :, :, 2])
      Rs_inv = np.zeros((self.x1.shape[0], self.x1.shape[1], 4, 4))
      # This is dx_hat/dx
      Rs_inv[:, :, 0, 0] = self.fluid_u_lower[..., 0]
      Rs_inv[:, :, 0, 1] = self.fluid_u_lower[..., 1]
      Rs_inv[:, :, 0, 2] = self.fluid_u_lower[..., 2]
      Rs_inv[:, :, 0, 3] = self.fluid_u_lower[..., 3]
      Rs_inv[:, :, 1, 0] = e2_vec_lower[..., 0]
      Rs_inv[:, :, 1, 1] = e2_vec_lower[..., 1]
      Rs_inv[:, :, 1, 2] = e2_vec_lower[..., 2]
      Rs_inv[:, :, 1, 3] = e2_vec_lower[..., 3]
      Rs_inv[:, :, 2, 0] = e3_vec_lower[..., 0]
      Rs_inv[:, :, 2, 1] = e3_vec_lower[..., 1]
      Rs_inv[:, :, 2, 2] = e3_vec_lower[..., 2]
      Rs_inv[:, :, 2, 3] = e3_vec_lower[..., 3]
      Rs_inv[:, :, 3, 0] = fluid_b_lower[..., 0]
      Rs_inv[:, :, 3, 1] = fluid_b_lower[..., 1]
      Rs_inv[:, :, 3, 2] = fluid_b_lower[..., 2]
      Rs_inv[:, :, 3, 3] = fluid_b_lower[..., 3]
      self.__dict__[key] = Rs_inv
    elif key == "frf_T_munu":
      T_munu = self.stress_e + self.stress_p
      T_munu_frf = np.einsum('...ki,...lj,...kl->...ij', self.frf_transform, self.frf_transform, T_munu)
      self.__dict__[key] = T_munu_frf
    elif key == "plasma_temp":
      indices = np.where(self.n_proper > 0)
      pressure = np.zeros_like(self.n_proper)
      pressure[indices] = (self.frf_T_munu[:, :, 1, 1] + self.frf_T_munu[:, :, 2, 2] + self.frf_T_munu[:, :, 3, 3])[indices] / 3.0
      self.__dict__[key] = pressure / self.n_proper
    elif key == "pressure_para":
      self.__dict__[key] = self.frf_T_munu[:, :, 3, 3]
    elif key == "pressure_perp":
      self.__dict__[key] = (self.frf_T_munu[:, :, 1, 1] + self.frf_T_munu[:, :, 2, 2]) / 2.0
    elif key == "plasma_beta":
      self.__dict__[key] = self.plasma_temp * self.n_proper / (0.5 * self.inner_product_4d_contravariant(self.frf_B, self.frf_B) + 1e-6)
      # self.__dict__[key] = self.plasma_temp * self.n_proper / (0.5 * inner_product_4d_contravariant(self.frf_B, self.frf_B, self._rv, self._thetav, self.a) + 1e-6)
    # elif key
    # elif key == "J":
    #   self._J = np.sqrt(self.J1 * self.J1 + self.J2 * self.J2 + self.J3 * self.J3)
    # elif key == "EdotB":
    #   self._EdotB = self.E1 * self.B1 + self.E2 * self.B2 + self.E3 * self.B3
    # elif key == "JdotB":
    #   self._JdotB = self.J1 * self.B1 + self.J2 * self.B2 + self.J3 * self.B3
      # elif key == "EdotB":
      #     setattr(self, "_" + key, data["EdotBavg"][()])
    else:
      data = h5py.File(path, "r")
      self.__dict__[key] = data[key][()]
      data.close()

# Compute the local flux 4 vector of electrons
def compute_fluid_4flux_e_upper(data):
  a = data.conf["bh_spin"]
  u0 = gu00(data._rv, data._thetav, a) * data.num_e + gu01(data._rv, data._thetav, a) * data.flux_e1
  u1 = gu11(data._rv, data._thetav, a) * data.flux_e1 + gu01(data._rv, data._thetav, a) * data.num_e + gu13(data._rv, data._thetav, a) * data.flux_e3
  u2 = gu22(data._rv, data._thetav, a) * data.flux_e2
  u3 = gu33(data._rv, data._thetav, a) * data.flux_e3 + gu13(data._rv, data._thetav, a) * data.flux_e1
  # Stack the components along the last axis to get shape (..., 4)
  return np.stack([u0, u1, u2, u3], axis=-1)

# Compute the local flux 4 vector of positrons
def compute_fluid_4flux_p_upper(data):
  a = data.conf["bh_spin"]
  u0 = gu00(data._rv, data._thetav, a) * data.num_p + gu01(data._rv, data._thetav, a) * data.flux_p1
  u1 = gu11(data._rv, data._thetav, a) * data.flux_p1 + gu01(data._rv, data._thetav, a) * data.num_p + gu13(data._rv, data._thetav, a) * data.flux_p3
  u2 = gu22(data._rv, data._thetav, a) * data.flux_p2
  u3 = gu33(data._rv, data._thetav, a) * data.flux_p3 + gu13(data._rv, data._thetav, a) * data.flux_p1
  # Stack the components along the last axis to get shape (..., 4)
  return np.stack([u0, u1, u2, u3], axis=-1)

# Compute the local flux 4 vector of electrons and positrons as a single fluid
def compute_fluid_4flux_upper(data):
  flux_upper_e = compute_fluid_4flux_e_upper(data)
  flux_upper_p = compute_fluid_4flux_p_upper(data)
  return flux_upper_e + flux_upper_p

# Compute the proper density of electrons
def compute_fluid_proper_density_e(data):
  flux_upper_e = compute_fluid_4flux_e_upper(data)
  # flux_upper_p = compute_fluid_4flux_p_upper(data)
  n_e = np.sqrt(np.abs(flux_upper_e[:,:,0]*data.num_e + flux_upper_e[:,:,1]*data.flux_e1 + flux_upper_e[:,:,2]*data.flux_e2 + flux_upper_e[:,:,3]*data.flux_e3))
  # n_p = np.sqrt(np.maximum(0.0, flux_upper_p[:,:,0]*data.num_p + flux_upper_p[:,:,1]*data.flux_p1 + flux_upper_p[:,:,2]*data.flux_p2 + flux_upper_p[:,:,3]*data.flux_p3))
  return n_e

# Compute the proper density of positrons
def compute_fluid_proper_density_p(data):
  flux_upper_p = compute_fluid_4flux_p_upper(data)
  n_p = np.sqrt(np.abs(flux_upper_p[:,:,0]*data.num_p + flux_upper_p[:,:,1]*data.flux_p1 + flux_upper_p[:,:,2]*data.flux_p2 + flux_upper_p[:,:,3]*data.flux_p3))
  return n_p

# Compute the proper density of electrons and positrons as a single fluid
def compute_fluid_proper_density(data):
  flux_upper = compute_fluid_4flux_upper(data)
  n_e = np.sqrt(np.abs(flux_upper[:,:,0]*data.num_e + flux_upper[:,:,1]*data.flux_e1 + flux_upper[:,:,2]*data.flux_e2 + flux_upper[:,:,3]*data.flux_e3))
  n_p = np.sqrt(np.abs(flux_upper[:,:,0]*data.num_p + flux_upper[:,:,1]*data.flux_p1 + flux_upper[:,:,2]*data.flux_p2 + flux_upper[:,:,3]*data.flux_p3))
  return n_e + n_p

# Inner product of two 4d contravariant vectors
def inner_product_4d_contravariant(v1, v2, r, th, a):
  g00 = gd00(r, th, a)
  g01 = gd01(r, th, a)
  g03 = gd03(r, th, a)
  g11 = gd11(r, th, a)
  g13 = gd13(r, th, a)
  g22 = gd22(r, th, a)
  g33 = gd33(r, th, a)
  return (g00 * v1[...,0] * v2[...,0] +
          2 * g01 * v1[...,0] * v2[...,1] +
          2 * g03 * v1[...,0] * v2[...,3] +
          g11 * v1[...,1] * v2[...,1] +
          2 * g13 * v1[...,1] * v2[...,3] +
          g22 * v1[...,2] * v2[...,2] +
          g33 * v1[...,3] * v2[...,3])

# Inner product of two 4d covariant vectors
def inner_product_4d_covariant(v1, v2, r, th, a):
  g00 = gu00(r, th, a)
  g01 = gu01(r, th, a)
  g11 = gu11(r, th, a)
  g13 = gu13(r, th, a)
  g22 = gu22(r, th, a)
  g33 = gu33(r, th, a)
  return (g00 * v1[...,0] * v2[...,0] +
          2 * g01 * v1[...,0] * v2[...,1] +
          g11 * v1[...,1] * v2[...,1] +
          2 * g13 * v1[...,1] * v2[...,3] +
          g22 * v1[...,2] * v2[...,2] +
          g33 * v1[...,3] * v2[...,3])

# Inner product of two 3d contravariant vectors
def inner_product_3d_contravariant(v1, v2, r, th, a):
  g11 = gmd11(r, th, a)
  g13 = gmd13(r, th, a)
  g22 = gmd22(r, th, a)
  g33 = gmd33(r, th, a)
  return (g11 * v1[...,0] * v2[...,0] +
          g22 * v1[...,1] * v2[...,1] +
          g33 * v1[...,2] * v2[...,2]
          + 2.0 * g13 * v1[...,0] * v2[...,2]
         )

# Inner product of two 3d covariant vectors
def inner_product_3d_covariant(v1, v2, r, th, a):
  g11 = gmu11(r, th, a)
  g13 = gmu13(r, th, a)
  g22 = gmu22(r, th, a)
  g33 = gmu33(r, th, a)
  return (g11 * v1[...,0] * v2[...,0] +
          g22 * v1[...,1] * v2[...,1] +
          g33 * v1[...,2] * v2[...,2]
          + 2.0 * g13 * v1[...,0] * v2[...,2]
         )

# Raise a 4d covariant vector to a contravariant vector 
def raise_4d_vec(v, r, th, a):
  g00 = gu00(r, th, a)
  g01 = gu01(r, th, a)
  g11 = gu11(r, th, a)
  g13 = gu13(r, th, a)
  g22 = gu22(r, th, a)
  g33 = gu33(r, th, a)
  u0 = g00 * v[...,0] + g01 * v[...,1]
  u1 = g01 * v[...,0] + g11 * v[...,1] + g13 * v[...,3]
  u2 = g22 * v[...,2]
  u3 = g13 * v[...,1] + g33 * v[...,3]
  return np.stack([u0, u1, u2, u3], axis=-1)

# Lower a 4d contravariant vector to a covariant vector
def lower_4d_vec(v, r, th, a):
  g00 = gd00(r, th, a)
  g01 = gd01(r, th, a)
  g03 = gd03(r, th, a)
  g11 = gd11(r, th, a)
  g13 = gd13(r, th, a)
  g22 = gd22(r, th, a)
  g33 = gd33(r, th, a)
  l0 = g00 * v[...,0] + g01 * v[...,1] + g03 * v[...,3]
  l1 = g01 * v[...,0] + g11 * v[...,1] + g13 * v[...,3]
  l2 = g22 * v[...,2]
  l3 = g13 * v[...,1] + g33 * v[...,3] + g03 * v[...,0]
  return np.stack([l0, l1, l2, l3], axis=-1)
