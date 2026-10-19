import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vertcross, CoordPair,get_cartopy,latlon_coords

"""
wrfinput/outファイルの(緯度,経度)にマッチにする(y, x)の組を取得する
引数
filepass -> 位置を決定したいwrfinput/outのファイルパス 基本絶対パス
lat, lon -> 変換したい緯度, 経度
"""
def latlon_to_yx(filepass, lat, lon):
  input=Dataset(filepass)
  lat_all=getvar(input, "lat")
  lon_all=getvar(input, "lon")
  dist=np.sqrt((lat_all - lat)**2 + (lon_all - lon)**2)
  y, x=np.unravel_index(np.argmin(dist.values), dist.shape)
  print(f"(lat, lon) = ({lat}, {lon}) to (y, x)= ({y}, {x})")

  return [y, x]

"""
任意のwrfoutファイルから鉛直断面を取り出す
引数
wrfout -> wrfoutのデータセット Dataset(ファイルパス)で取り出せる
var, unit -> 変数名, 単位
start, end -> 断面の開始と終了位置 [y, x]のリストで上記のlatlon_to_yxによって求める
vcord -> 鉛直座標 z(高度) or p(気圧)

*wrfoutの鉛直座標は時々刻々変化するので固定の鉛直座標を作ってそこに線形補完する
"""
def get_vertical_cross_section(wrfout, var, unit, start, end, vcord):
  start_point=CoordPair(x=start[1], y=start[0])
  end_point=CoordPair(x=end[1], y=end[0])

  ds=wrfout
  var_all=getvar(ds, var, units=unit)

  if vcord == "z":
    vert=getvar(ds,"z")
    var_vert=vertcross(var_all, vert, start_point=start_point, end_point=end_point, latlon=True)
#いったん断面を取得

    z_fixed = np.arange(0, 20001, 200)
    var_vert = var_vert.interp(vertical=z_fixed)
#固定座標系を参照して線形補間

  if vcord == "p":
    vert=getvar(ds, "pressure")
    var_vert=vertcross(var_all, vert, start_point=start_point, end_point=end_point, latlon=True)

    p_fixed=np.arange(1010, 101, -10)
    var_vert = var_vert.interp(vertical=p_fixed)

  return var_vert

"""
getvarでの単位指定がない変数
"""
def get_vertical_cross_section_nounit(wrfout, var, start, end, vcord):
  start_point=CoordPair(x=start[1], y=start[0])
  end_point=CoordPair(x=end[1], y=end[0])

  ds=wrfout
  var_all=getvar(ds, var)

  if vcord == "z":
    vert=getvar(ds,"z")
    var_vert=vertcross(var_all, vert, start_point=start_point, end_point=end_point, latlon=True)
#いったん断面を取得

    z_fixed = np.arange(0, 20001, 200)
    var_vert = var_vert.interp(vertical=z_fixed)
#固定座標系を参照して線形補間

  if vcord == "p":
    vert=getvar(ds, "pressure")
    var_vert=vertcross(var_all, vert, start_point=start_point, end_point=end_point, latlon=True)

    p_fixed=np.arange(1010, 101, -10)
    var_vert = var_vert.interp(vertical=p_fixed)

  return var_vert
"""
断面に平行な風ベクトルを取り出す？
"""