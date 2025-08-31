#!/usr/bin/env python

import h5py
import numpy as np
import toml
from pathlib import Path
import os
import re
from datalib import Data, flag_to_species
from datalib_logsph import DataSph

def rho2(a, r, th):
  return r * r + (a * np.cos(th))**2

def sqrt_gamma(a, r, th):
  r2 = rho2(a, r, th)
  return r2 * np.sin(th) * np.sqrt(1.0 + 2.0 * r / r2)

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


class DataKerrSchild(DataSph):
  _mesh_loaded = False

  def __init__(self, path):
    super().__init__(path)
    self.a = self._conf["bh_spin"]
    self.rH = rs_o(self.a)
    self.extra_fld_keys = ["fluxB", "Dd1", "Dd2", "Dd3", "Bd1", "Bd2", "Bd3", "Ed1", "Ed2", "Ed3", "Hd1", "Hd2", "Hd3", "sigma", "flux_upper", "flux_lower",
                           "n_proper", "fluid_u_upper", "fluid_u_lower", "fluid_b_upper"]
    self.reload()

#   def load_sph_mesh(self):
#      return super().load_sph_mesh()

#   def __getattr__(self, key):
#     if key not in self.__dict__:
#       if key in self._fld_keys:
#         self._load_fld_quantity(key)
#       elif key in self._ptc_keys:
#         self._load_ptc_quantity(key)
#       elif key in self._mesh_keys:
#         self._load_mesh_quantity(key)
#       elif key == "keys":
#         self.__dict__[key] = self._fld_keys + self._ptc_keys + self._mesh_keys
#       elif key == "conf":
#         self.__dict__[key] = self._conf
#       else:
#         return None
#     return self.__dict__[key]

  def _load_fld_quantity(self, key):
    path = os.path.join(self._path, f"fld.{self._current_fld_step:05d}.h5")
    if key == "fluxB":
      self._load_sph_mesh()
      self.__dict__[key] = np.cumsum(self.B1 * sqrt_gamma(self.a, self._rv, self._thetav) * self._dtheta, axis=0)
    # Lower components of D
    elif key == "Dd1":
      self.__dict__[key] = gmd11(self._rv, self._thetav, self.a) * self.E1 + gmd13(self._rv, self._thetav, self.a) * self.E3
    elif key == "Dd2":
      self.__dict__[key] = gmd22(self._rv, self._thetav, self.a) * self.E2
    elif key == "Dd3":
      self.__dict__[key] = gmd13(self._rv, self._thetav, self.a) * self.E1 + gmd33(self._rv, self._thetav, self.a) * self.E3
    # Lower components of B
    elif key == "Bd1":
      self.__dict__[key] = gmd11(self._rv, self._thetav, self.a) * self.B1 + gmd13(self._rv, self._thetav, self.a) * self.B3
    elif key == "Bd2":
      self.__dict__[key] = gmd22(self._rv, self._thetav, self.a) * self.B2
    elif key == "Bd3":
      self.__dict__[key] = gmd13(self._rv, self._thetav, self.a) * self.B1 + gmd33(self._rv, self._thetav, self.a) * self.B3
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
    elif key == "sigma": # this is the cold sigma
      self.__dict__[key] = (self.B1 * self.Bd1 + self.B2 * self.Bd2 + self.B3 * self.Bd3) / (self.Rho_p - self.Rho_e + 1e-6)
    elif key == "flux_upper":
      self.__dict__[key] = compute_fluid_4flux_upper(self)
    elif key == "flux_lower":
      self.__dict__[key] = np.stack([self.num_e + self.num_p, self.flux_e1 + self.flux_p1,
                                     self.flux_e2 + self.flux_p2, self.flux_e3 + self.flux_p3], axis=-1)
    elif key == "n_proper":
      self.__dict__[key] = np.sqrt(np.abs(inner_product_4d_covariant(self.flux_lower, self.flux_lower, self._rv, self._thetav, self.a)))
    elif key == "fluid_u_upper":
      indices = np.where(self.n_proper > 0)
      u_upper = self.flux_upper
      u_upper[indices] = self.flux_upper[indices] / self.n_proper[indices][..., np.newaxis]
      self.__dict__[key] = u_upper
    elif key == "fluid_u_lower":
      indices = np.where(self.n_proper > 0)
      u_lower = self.flux_lower
      u_lower[indices] = self.flux_lower[indices] / self.n_proper[indices][..., np.newaxis]
      self.__dict__[key] = u_lower
    elif key == "fluid_b_upper":
      # compute b vector in the fluid rest frame
      B = np.stack([np.zeros_like(data.B1), data.B1, data.B2, data.B3], axis=-1)
      E = np.stack([np.zeros_like(data.Ed1), data.Ed1, data.Ed2, data.Ed3], axis=-1)
      u_lower = self.fluid_u_lower
      # B = np.array([0.0, data.B1[Nth, Nr], data.B2[Nth, Nr], data.B3[Nth, Nr]])
      # E = np.array([0.0, data.Ed1[Nth, Nr], data.Ed2[Nth, Nr], data.Ed3[Nth, Nr]])
      alpha_val = alpha(self._rv, self._thetav, self.a)
      sqrt_gm = gmsqrt(self._rv, self._thetav, self.a)
      b0 = (u_lower[...,1] * B[...,1] + u_lower[...,2] * B[...,2] + u_lower[...,3] * B[...,3]) / alpha_val
      b1 = -u_lower[...,0] * B[...,1] / alpha_val - (u_lower[...,2] * E[...,3] - u_lower[...,3] * E[...,2]) / alpha_val / sqrt_gm
      b2 = -u_lower[...,0] * B[...,2] / alpha_val - (u_lower[...,3] * E[...,1] - u_lower[...,1] * E[...,3]) / alpha_val / sqrt_gm
      b3 = -u_lower[...,0] * B[...,3] / alpha_val - (u_lower[...,1] * E[...,2] - u_lower[...,2] * E[...,1]) / alpha_val / sqrt_gm
      b_upper = np.array([b0, b1, b2, b3])
      bnorm = np.sqrt(inner_product_4d_covariant(b_upper, b_upper, self._rv, self._thetav, self.a))
      self.__dict__[key] = np.stack([b0, b1, b2, b3], axis=-1) / bnorm

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