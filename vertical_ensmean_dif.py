###任意の鉛直断面を取った変数の平均と実験差を計算，描画する
# -*- coding: utf-8 -*-
from datetime import datetime,timedelta
from netCDF4 import Dataset
from wrf import getvar,get_cartopy,latlon_coords,geo_bounds,interplevel
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
from scipy.stats import ttest_ind
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.vertical_cross_section import latlon_to_yx, get_vertical_cross_section
from myFunction.draw_field import draw_vertical_cross, draw_vertical_cross_onlysig
from myFunction.ttest_Nstar import ttest_vertical_cross_ensmean
from myFunction.ttest import ttest_vertical_cross

N=5
print("N (Ensemble Member) = ",N)
#アンサンブル数

case="2025winter"
ex1="CTL"
ex2="ME00"

var="z"
unit="m"

# vcord="z"
# vunit="m"
vcord="p"
vunit="hPa"

vtop=300
vbottom=1000
#鉛直座標


#断面図の位置設定
# start_lat,start_lon=40.5,122
# end_lat,end_lon=40.5,137
start_lat,start_lon=50,120
end_lat,end_lon=30,135

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir1=f"{wrf_dir}/wrfout/{case}/{ex1}"
wrfout_dir2=f"{wrf_dir}/wrfout/{case}/{ex2}"
fig_dir=f"{wrf_dir}/fig/{case}"

wrfinput=f"{wrfout_dir1}/n1/wrfinput_d01"

start=latlon_to_yx(wrfinput, start_lat, start_lon)
end=latlon_to_yx(wrfinput, end_lat, end_lon)

##検定
have_test="yes"
#平均値の差の検定の有無 yes or no
use_Nstar="yes"
#有効サンプル数を使った検定をするか
#yes ->各メンバでの有効サンプル数を計算しその和を使って自由度を計算する
#時間がかかるのでsign_mapは保存し，書き直す時をそれを使う
#no ->アンサンブル数を使って自由度を計算する
if have_test == "yes" and use_Nstar =="yes":
  test="TTest_Nstar"
  print("test -> ", test)
elif have_test == "yes" and use_Nstar == "no":
  test="TTest"
  print("test -> ", test)
else:
  test=""   
  print("No Test")
  sign_file=f"{wrfout_dir2}/sign_map_vertcross_{var}_mean_ens(N={N}).nc"
  sign_map=np.load(sign_file)

start_date=datetime(2024, 12, 1, 0)
end_date=datetime(2025, 2, 28, 18)
dh=3
print(f"Set Time : {start_date}----{end_date}")


mean_list1=[]
mean_list2=[]
#各メンバの期間平均値を保持するリスト

if use_Nstar == "yes":
  ens_list1=[]
  ens_list2=[]

n=1
while n <= N:
  date=start_date

  var_list1=[]
  var_list2=[]
#各時刻の値を保持するリスト

  while date <= end_date:
    y=date.year
    m=date.month
    d=date.day
    h=date.hour

    time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}"
    print(f"Read wrfout File {time_str} (N={n})")

#ファイルオープンとデータ取得
    wrfout1=f"{wrfout_dir1}/n{n}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/n{n}/wrfout_d01_{time_str}:00:00"
    ds1=Dataset(wrfout1)
    ds2=Dataset(wrfout2)
  #変数取り出す
    
#断面データを取り出す
    var_vert1=get_vertical_cross_section(ds1, var, unit, start, end, vcord)
    var_vert2=get_vertical_cross_section(ds2, var, unit, start, end, vcord)
 #   print(var_vert1)
#配列は(鉛直, 水平)

#座標情報を取りだす
    if n == N and date == end_date:
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
    #グリッド線の設定
      gl = ax.gridlines(draw_labels=True, 
                    linewidth=1, 
                    color='gray', 
                    alpha=0.5, 
                    linestyle='--')
#ラベルの表示位置を制御
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
    
    ds1.close()
    ds2.close()

    var_list1.append(var_vert1)
    var_list2.append(var_vert2)

    date+=timedelta(hours=dh)

#鉛直座標が動的であることが確認できる
  # vref = var_list1[0].vertical.values

  # for i, v in enumerate(var_list1):
  #   diff = np.max(np.abs(vref - v.vertical.values))
  #   print(i, diff)

  # for i in [0, 6, 7, 8]:
  #   print(f"i={i}")
  #   print(var_list1[i].vertical.values)

  darray1=np.stack(var_list1, axis=0)
  darray2=np.stack(var_list2, axis=0)

  if test == "TTest_Nstar":
    ens_list1.append(darray1)
    ens_list2.append(darray2)
  
  else: 
    mean1=np.nanmean(darray1, axis=0)
    mean2=np.nanmean(darray2, axis=0)    
    mean_list1.append(mean1)
    mean_list2.append(mean2)

  n+=1

if test == "TTest_Nstar":
  enstime_darray1=np.stack(ens_list1,axis=0)
  enstime_darray2=np.stack(ens_list2,axis=0)
#4次元配列(member, time, south_north, west_east)

else:
  ens_darray1=np.stack(mean_list1, axis=0)
  ens_darray2=np.stack(mean_list2, axis=0)
 
#3次元配列(member, south_north, west_east)
#print(ens_darray1)
#darray[member, south_north, west_east]

##ココで母平均の差の検定
if have_test == "yes" :
  if use_Nstar == "yes":
    sign_map=ttest_vertical_cross_ensmean(enstime_darray1, enstime_darray2, 0.95)
    sign_path=f"{wrfout_dir2}/sign_map_vertcross_{var}_mean_ens.npy"
    np.save(sign_path, sign_map)
    print("Save -> ",sign_path)
    ens_darray1=enstime_darray1.mean(dim="time", skipna=True)
    ens_darray2=enstime_darray2.mean(dim="time", skipna=True)
#時間方向に圧縮   
 
  if use_Nstar == "no":
    sign_map=ttest_vertical_cross(ens_darray1, ens_darray2, 0.95)

#平均値と差を計算
ens_mean1=np.nanmean(ens_darray1, axis=0)
ens_mean2=np.nanmean(ens_darray2, axis=0)
ens_dif=ens_mean1 - ens_mean2

nan1=np.isnan(ens_mean1)
nan2=np.isnan(ens_mean2)
dummy=np.zeros_like(ens_mean1)
#欠損値の配列も作る

##Plot 
draw_vertical_cross(idx, vert, ens_mean1, vcord, vbottom, vtop, 
                     "viridis", 0, 10000, 500, "max", "Z", "m", "", nan1)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_{ex1}_mean_ens(N={N}).png")
plt.close("all")

draw_vertical_cross(idx, vert, ens_mean2, vcord, vbottom, vtop, 
                     "viridis", 0, 10000, 500, "max", "Z", "m", "", nan2)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_{ex2}_mean_ens(N={N}).png")
plt.close("all")

draw_vertical_cross(idx, vert, ens_dif, vcord, vbottom, vtop,
                     "seismic", -18, 18, 4, "both", "Z", "m", "", nan1)
plt.contourf(idx ,vert, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"])

# draw_vertical_cross_onlysig(idx, vert, ens_dif, vcord, vbottom, vtop,
#                     "seismic", -18, 18, 4, "both", "Z", "m", "", nan1, sign_map)
plt.savefig(f"{fig_dir}/Vertical_Cross_{var}_Dif_mean_ens(N={N}).png")
plt.close("all")

# levs=np.arange(260,324,4)
# levs_dif=np.arange(-7,9,2.0)
# #levc=np.arange(0,10000,500)
# #levc_dif=np.arange(-18,22,4)
# levc=np.arange(260,324,4)
# levc_dif=np.arange(-7,9,2.0)
# #levs=np.arange(0,10000,500)
# #levs_dif=np.arange(-18,22,4)