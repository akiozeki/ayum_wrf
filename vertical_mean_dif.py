###任意の鉛直断面を取った変数の平均と実験差を計算，描画する
from __future__ import print_function, division
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from wrf import getvar, vertcross, CoordPair,get_cartopy,latlon_coords
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.vertical_cross_section import latlon_to_yx, get_vertical_cross_section
from myFunction.draw_field import draw_vertical_cross, draw_vertical_cross_onlysig
#from myFunction.ttest_Nstar import ttest_vertical_cross
from myFunction.ttest import ttest_vertical_cross


n=1
#アンサンブルメンバ, アンサンブルを取っていない場合は0

case="2025winter"
ex1="CTL"
ex2="ME00"

var="z"
unit="m"

#vcord="z"
#vunit="m"
vcord="p"
vunit="hPa"

vtop=300
vbottom=1000
#鉛直座標

have_test="yes"
#検定の有無 yes or no

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir1=f"{wrf_dir}/wrfout/{case}/{ex1}"
wrfout_dir2=f"{wrf_dir}/wrfout/{case}/{ex2}"
fig_dir=f"{wrf_dir}/fig/{case}"

start_date=datetime(2024, 12, 1, 0)
end_date=datetime(2025, 2, 28, 18)
dh=3
print(f"Set Time : {start_date}----{end_date}")

var_list1=[]
var_list2=[]

#断面図の位置設定
# start_lat,start_lon=40.5,122
# end_lat,end_lon=40.5,137
start_lat,start_lon=50,120
end_lat,end_lon=30,135
if n == 0:
  wrfinput=f"{wrfout_dir1}/wrfinput_d01"
if n != 0:
  wrfinput=f"{wrfout_dir1}/n{n}/wrfinput_d01"

start=latlon_to_yx(wrfinput, start_lat, start_lon)
end=latlon_to_yx(wrfinput, end_lat, end_lon)

date=start_date
while date <= end_date:
  y=date.year
  m=date.month
  d=date.day
  h=date.hour

  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}"
  if n == 0:
    wrfout1=f"{wrfout_dir1}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/wrfout_d01_{time_str}:00:00"

  if n != 0:
    wrfout1=f"{wrfout_dir1}/n{n}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/n{n}/wrfout_d01_{time_str}:00:00"
#ここで「n〇」とすれば任意のアンサンブルメンバーを1つ取り出せる

  print(f"Read wrfout File {time_str}")

  ds1=Dataset(wrfout1)
  ds2=Dataset(wrfout2)

#断面データを取り出す
  var_vert1=get_vertical_cross_section(ds1, var, unit, start, end, vcord)
  var_vert2=get_vertical_cross_section(ds2, var, unit, start, end, vcord)
 # print(var_vert1)
  #配列は(鉛直, 水平)

#座標情報を取りだす
  if date == end_date:
#ここでどの断面を取ったかを表す図も描いておく
    ter=getvar(ds1,"ter",units="m")
    cart_proj=get_cartopy(ter)
    lats,lons=latlon_coords(ter)
    fig=plt.figure()
    ax=plt.axes(projection=cart_proj)
    ax.plot(
      [start_lon,end_lon],
      [start_lat,end_lat],
      color="blue",
      linewidth=3,
      transform=ccrs.PlateCarree()
    )

    ax.contourf(
      lons,lats,ter,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
      transform=ccrs.PlateCarree()
    )
    # グリッド線の設定
    gl = ax.gridlines(draw_labels=True, 
                  linewidth=1, 
                  color='gray', 
                  alpha=0.5, 
                  linestyle='--')
# ラベルの表示位置を制御
    gl.xlines = True         # 経度線を描く
    gl.ylines = True         # 緯度線を描く
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator=FixedLocator([127,130,133,136,139])
    #gl.ylocator=FixedLocator()
    gl.x_inline = False      # ラベルを図の内側（インライン）に書かない設定
    gl.y_inline = False      # 緯度も念のため設定
    gl.ylabel_style = {'rotation': 0}
    gl.xlabel_style = {'rotation': 0}
    gl.xpadding = 10
    ax.coastlines()

    plt.savefig(f"{fig_dir}/vert_cross_line{start}to({end}.png")
    plt.close("all")

    vert=var_vert1.vertical.values
    idx=var_vert1.cross_line_idx.values



  var_list1.append(var_vert1)
  var_list2.append(var_vert2)


  ds1.close()
  ds2.close()

  date+=timedelta(hours=dh)

var_array1=np.stack(var_list1, axis=0)
var_array2=np.stack(var_list2, axis=0)

##平均値の差の検定
if have_test == "yes":
  sign_map=ttest_vertical_cross(var_array1, var_array2)

#平均と差を計算
var_mean1=np.nanmean(var_array1, axis=0)
var_mean2=np.nanmean(var_array2, axis=0)
var_dif=var_mean1 - var_mean2

nan1=np.isnan(var_mean1)
nan2=np.isnan(var_mean2)
dummy=np.zeros_like(var_mean1)
#欠損値の配列も作る


# #Plot 
draw_vertical_cross(idx, vert, var_mean1, vcord, vbottom, vtop, 
                     "viridis", 0, 10000, 500, "max", "Z", "m", "", nan1)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_{ex1}_mean_n={n}.png")
plt.close("all")

draw_vertical_cross(idx, vert, var_mean2, vcord, vbottom, vtop, 
                     "viridis", 0, 10000, 500, "max", "Z", "m", "", nan2)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_{ex2}_mean_n={n}.png")
plt.close("all")

draw_vertical_cross(idx, vert, var_dif, vcord, vbottom, vtop,
                     "seismic", -18, 18, 4, "both", "Z", "m", "", nan1)
plt.contourf(idx ,vert, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"])

#draw_vertical_cross_onlysig(idx, vert, var_dif, vcord, vbottom, vtop,
#                     "seismic", -18, 18, 4, "both", "Z", "m", "", nan1, sign_map)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_Dif_mean_n={n}.png")
plt.close("all")

# levs=np.arange(260,324,4)
# levs_dif=np.arange(-7,9,2.0)
# #levc=np.arange(0,10000,500)
# #levc_dif=np.arange(-18,22,4)
# levc=np.arange(260,324,4)
# levc_dif=np.arange(-7,9,2.0)
# #levs=np.arange(0,10000,500)
# #levs_dif=np.arange(-18,22,4)