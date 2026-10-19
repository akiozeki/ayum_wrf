import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

"""
水平面図におけるグリッドを表示
"""
def draw_grid(ax):
  gl = ax.gridlines(draw_labels=True, 
                  linewidth=1, 
                  color='gray', 
                  alpha=0.5, 
                  linestyle='--')
# ラベルの表示位置を制御
  gl.xlines = True         # 経度線を描く
  gl.ylines = True         # 緯度線を描く
  gl.top_labels = False
  gl.right_labels = False
  gl.x_inline = False      # ラベルを図の内側（インライン）に書かない設定
  gl.y_inline = False      # 緯度も念のため設定
  gl.ylabel_style = {'rotation': 0}
  gl.xlabel_style = {'rotation': 0}
  gl.xpadding = 10


"""
2次元マップを描画する
引数
lon, lat -> 経度, 経度の一次元配列
field -> 描画したい変数の2次元配列
projection -> WRFの地図投影法, get_cartopy()で取得する
domain -> 描画範囲[], 0にすると自動
cmap -> カラーマップ
vmin, vmax , vint -> カラーバーの最小，最大，間隔
extend -> カラーバーの表示形式, "min" or "max", or "both"
fname, unit -> 変数名, 単位
fig_title -> 図のタイトル
"""
def draw_latlon_contourf(lon, lat, field, projection, domain, cmap, vmin, vmax, vint, extend, fname, unit, fig_title):
  fig=plt.figure()
  ax=plt.axes(projection=projection)
  if domain != 0:
    ax.set_extent(domain)

  shade=ax.contourf(
    lon, lat, field,
    cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
    levels=np.arange(vmin, vmax+vint, vint), extend=extend,
    transform=ccrs.PlateCarree()
  )

  cbar=plt.colorbar(
    shade,orientation="vertical",
    shrink=0.8
    )
  
  cbar.set_label(f"{fname} [ {unit} ]", rotation=0, labelpad=20)

# グリッド線の設定
  draw_grid(ax)

  ax.coastlines()
#  ax.add_feature(cfeature.LAND,color="gray")
  ax.set_title(fig_title,fontsize=18)


"""
有意な所のみ描画する
sig_map -> 有意な所に1が入っている2次元配列 サイズはフィールドと一致している必要がある
"""
def draw_latlon_contourf_onlysig(lon, lat, field, projection, domain, cmap, vmin, vmax, vint, extend, fname, unit, fig_title, sig_map):
  sig_field=np.where(sig_map == 1, field, np.nan)
#有意でないところはNanマスク
  fig=plt.figure()
  ax=plt.axes(projection=projection)
  if domain != 0:
    ax.set_extent(domain)

  shade=ax.contourf(
    lon, lat, sig_field,
    cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
    levels=np.arange(vmin, vmax+vint, vint), extend=extend,
    transform=ccrs.PlateCarree()
  )

  cbar=plt.colorbar(
    shade,orientation="vertical",
    shrink=0.8
    )
  
  cbar.set_label(f"{fname} [ {unit} ]", rotation=0, labelpad=20)

# グリッド線の設定
  draw_grid(ax)
  ax.coastlines()
#  ax.add_feature(cfeature.LAND,color="gray")
  ax.set_title(fig_title,fontsize=18)


"""
2次元マップを描画する(等値線)
引数
lon, lat -> 経度, 経度の一次元配列
field -> 描画したい変数の2次元配列
projection -> WRFの地図投影法, get_cartopy()で取得する
domain -> 描画範囲[], 0にすると自動
vmin, vmax , vint -> カラーバーの最小，最大，間隔
fig_title -> 図のタイトル
"""
def draw_latlon_contour(lon, lat, field, projection, domain, vmin, vmax, vint, fig_title):
  fig=plt.figure(figsize=(8,8))
  ax=plt.axes(projection=projection)
  if domain != 0:
    ax.set_extent(domain)

  contour=ax.contour(
    lon, lat, field,
    levels=np.arange(vmin, vmax+vint, vint), 
    transform=ccrs.PlateCarree()
  )
  ax.clabel(contour)

# グリッド線の設定
  draw_grid(ax)
  ax.coastlines()
#  ax.add_feature(cfeature.LAND,color="gray")
  ax.set_title(fig_title,fontsize=18)


  """
鉛直断面図を描画する
引数
idx, vert -> 水平, 鉛直座標 wrf-pythonから取得する(ert=var_vert.vertical.values, idx=var_vert.cross_line_idx.values)
field -> 描画したい変数の鉛直断面
vcord -> 鉛直座標 z(高度) or p(気圧)
vbottom, vtop
cmap -> カラーマップ
vmin, vmax , vint -> カラーバーの最小，最大，間隔
extend -> カラーバーの表示形式, "min" or "max", or "both"
fname, unit -> 変数名, 単位
fig_title -> 図のタイトル
nan -> 欠損値配列 fieldと同じ次元
"""
def draw_vertical_cross(idx, vert, field, vcord, vbottom, vtop, cmap, vmin, vmax, vint, extend, fname, unit, fig_title, nan):
  plt.rcParams["font.size"]=12
  fig=plt.figure()
  ax=plt.axes()
  ax.set_xticks([])      
  ax.set_xlabel("Horizontal")
  ax.set_ylabel(f"Vertical ({vcord})")
  ax.set_ylim(vbottom, vtop)

  ax.contourf(idx,vert,nan,levels=[0.5,1.5],colors=["black"])

  shade=ax.contourf(
    idx, vert, field,
    cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
    levels=np.arange(vmin, vmax+vint, vint), extend=extend
  )

  cbar=plt.colorbar(
      shade,orientation="vertical",
      shrink=0.8
      )
    
  cbar.set_label(f"{fname} [ {unit} ]", rotation=0, labelpad=20)
   
  contour=ax.contour(idx, vert, field, levels=np.arange(vmin, vmax+vint, vint), colors="black",linewidths=1.0)
  ax.clabel(contour)

  ax.set_title(fig_title)
  plt.tight_layout()

#風ベクトルも描きたい  
# xstep=4
# vstep=4
# qx,qy,qk=0.9,1.03,10
# vector=ax.quiver(
#  idx[::xstep],
#  vert[::vstep],
#  u1_mean[::vstep,::xstep],
#  vdummy[::vstep,::xstep],
#  scale=200,
#  width=0.008,
#  color="green"
#  )
# ax.quiverkey(vector,qx,qy,qk,f"{qk}m/s")


"""
有意な所のみ描画する
sig_map -> 有意な所に1が入っている2次元配列 サイズはフィールドと一致している必要がある
"""
def draw_vertical_cross_onlysig(idx, vert, field, vcord, vbottom, vtop, cmap, vmin, vmax, vint, extend, fname, unit, fig_title, nan, sig_map):
  sig_field=np.where(sig_map == 1, field, np.nan)

  plt.rcParams["font.size"]=12
  fig=plt.figure()
  ax=plt.axes()
  ax.set_xticks([])      
  ax.set_xlabel("Horizontal")
  ax.set_ylabel(f"Vertical ({vcord})")
  ax.set_ylim(vbottom, vtop)

  ax.contourf(idx,vert,nan,levels=[0.5,1.5],colors=["black"])

  shade=ax.contourf(
    idx, vert, sig_field,
    cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
    levels=np.arange(vmin, vmax+vint, vint), extend=extend
  )

  cbar=plt.colorbar(
      shade,orientation="vertical",
      shrink=0.8
      )
    
  cbar.set_label(f"{fname} [ {unit} ]", rotation=0, labelpad=20)
   
  contour=ax.contour(idx, vert, sig_field, levels=np.arange(vmin, vmax+vint, vint), colors="black",linewidths=1.0)
  ax.clabel(contour)

  ax.set_title(fig_title)
  plt.tight_layout()