import numpy as np
import xarray as xr

"""
z系での鉛直微分を計算する

引数
var -> 3次元 (bottom_top, south_north, west_east) のxarray配列
z -> 圧力の3次元配列
var_name -> 鉛直微分の変数名
*出力ファイルは時間および格子点ごとに鉛直座標の値が変わるので水平各格子点ごとに計算が必要
"""

def zdiff(var, z, var_name):
  dvar_dz=np.empty_like(var.values)
  nz, ny, nx=var.shape

#各格子点での鉛直微分を計算
  for j in range(ny):
    for i in range(nx):
      zcoord=z[:, j, i].values
      dvar_dz[:, j, i]=np.gradient(var[:, j, i].values, zcoord)

#xarrayに戻す
  dvar_dz = xr.DataArray(
    dvar_dz,
    coords=var.coords,
    dims=var.dims,
    name=var_name
  )

  return dvar_dz


"""
p系での鉛直微分を計算する

引数
var -> 3次元 (bottom_top, south_north, west_east) のxarray配列
p -> 圧力の3次元配列
var_name -> 鉛直微分の変数名
"""

def pdiff(var, p, var_name):
  dvar_dp=np.empty_like(var.values)
  nz, ny, nx=var.shape

#鉛直微分を計算
  for j in range(ny):
    for i in range(nx):
      pcoord=p[:, j, i].values
      dvar_dp[:, j, i]=np.gradient(var[:, j, i].values, pcoord)

#xarrayに戻す
  dvar_dp = xr.DataArray(
    dvar_dp,
    coords=var.coords,
    dims=var.dims,
    name=var_name
  )

  return dvar_dp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              