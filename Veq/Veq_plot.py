###別ファイルで作ったnetCDFファイルを使って指定した期間の平均を求める
from datetime import datetime
from netCDF4 import Dataset
import xarray as xr
import numpy as np
from wrf import getvar,interplevel,get_cartopy,latlon_coords
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys
from pathlib import Path

script_dir=Path(__file__).resolve().parent.parent
# myFunction と Veq の親ディレクトリを取得

sys.path.insert(0, str(script_dir))
# Python のモジュール検索パスに追加

from myFunction.ttest import ttest_2Dmap
from myFunction.ttest_Nstar import ttest_2Dmap_ensmean
from myFunction.draw_field import draw_grid

"""各種設定"""
type="ens"
#ens or singl

N=5
print("N (Ensemble Member) = ",N)
#アンサンブル数 (type == "ens") or アンサンブルメンバ (type == "singl")

#高度指定
lev=750

##変数指定
var="zeta_smooth"
#var="dtzeta_dt"
#var="Adv.zeta_x"
#var="Adv.zeta_y"
#var="Adv.zeta_h"
#var="Adv.zeta_p"
#var="Adv.f"
#var="Stretching"
#var="Tilting"
#var="Sum"
#var="Res"
#var=0

case="2025winter"
ex1="CTL"
ex2="ME00"

n_member=1
dx=25
sigma=6
#gtime="start"
gtime="end"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig/{case}/Veq"

domain=[120,140,32,50]
#図の描画範囲
#factor=0
#factor=5
factor=6
#factor=9
# factor=10
# factor=11
#描画のための倍率(10のfactor乗)


levc_dif=np.arange(-8.8,10.2,1.6)
levc_plus=np.arange(0.6,16.6,1.6)
levc_minus=np.arange(-15,0,1.6)

levc_plus=np.arange(1.6,17.7,3.2)
levc_minus=np.arange(-17.6,0,3.2)

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
  sign_file=f"{wrfout_dir}/{ex2}/sign_map_{var}_mean_ens(N={N}).npy"
  sign_map=np.load(sign_file)



"""座標系情報を得るために任意のwrfoutファイルを開く"""
sdate=datetime(2024,12,1,0)
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour

coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

ter1=getvar(coord_ds1, "ter")
ter2=getvar(coord_ds2, "ter")
land=getvar(coord_ds1,"LANDMASK")
cart_proj=get_cartopy(land)
print("projection:",cart_proj)
lats,lons=latlon_coords(land)
#print(lats)
coord_ds1.close()
coord_ds2.close()

"""データの読み込み"""
if type == "ens":
  dfile1=f"{wrfout_dir}/{ex1}/{ex1}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}_ensemble.nc"
  dfile2=f"{wrfout_dir}/{ex2}/{ex2}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}_ensemble.nc"


elif type == "singl":
  dfile1=f"{wrfout_dir}/{ex1}/n1/{ex1}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"
  dfile2=f"{wrfout_dir}/{ex2}/n1/{ex2}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"

ds1=xr.open_dataset(dfile1)
ds2=xr.open_dataset(dfile2)

#print(ds1)
#print(ds2)


#変数指定 -> var[member, time, south_north, west_east]
if var == "Adv.zeta_h" :
  var1=ds1["Adv.zeta_x"] + ds1["Adv.zeta_y"]
  var2=ds2["Adv.zeta_x"] + ds2["Adv.zeta_y"]


elif var == "Sum"  or var == "Res":
  sum1=ds1["Adv.zeta_x"] + ds1["Adv.zeta_y"] + ds1["Adv.zeta_p"] + ds1["Adv.f"] + ds1["Stretching"] + ds1["Tilting"]
  sum2=ds2["Adv.zeta_x"] + ds2["Adv.zeta_y"] + ds2["Adv.zeta_p"] + ds2["Adv.f"] + ds2["Stretching"] + ds2["Tilting"]
  
  if var == "Sum":
    var1=sum1
    var2=sum2

  elif var == "Res":
    var1=ds1["dzeta_dt"] - sum1
    var2=ds2["dzeta_dt"] - sum2

else:
  var1=ds1[var]
  var2=ds2[var]
 
#print(var1)

ds1.close()
ds2.close()


"""母平均の差の検定と平均"""
if type == "ens":
  if have_test == "yes" :
    if use_Nstar == "yes":
      sign_map=ttest_2Dmap_ensmean(var1, var2, 0.95)
      sign_path=f"{wrfout_dir}/{ex2}/sign_map_{var}_mean_ens(N={N}).npy"
      np.save(sign_path, sign_map)
      print("Save -> ",sign_path)

  var1=var1.mean(dim="member", skipna=True)
  var2=var2.mean(dim="member", skipna=True)
#時間方向に圧縮   
  if have_test == "yes" and use_Nstar == "no":
      sign_map=ttest_2Dmap(var1, var2, 0.95)
 
  
#平均値と差を計算
mean1=var1.mean(dim="time", skipna=True)
mean2=var2.mean(dim="time", skipna=True)
dif= mean1 - mean2

"""Plot""" 
#EX1
fig=plt.figure()
ax=plt.axes(projection=cart_proj)
ax.set_extent(domain)
ax.contourf(
      lons,lats,ter1,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
      transform=ccrs.PlateCarree()
    )

contour_plus=ax.contour(
  lons,lats,mean1*10**factor,
  levels=levc_plus,
  colors="black",linewidths=2.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus, fontsize=12)

contour_minus=ax.contour(
  lons,lats,mean1*10**factor,
  levels=levc_minus,
  colors="black",linewidths=2.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_minus, fontsize=12)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex1}-{var}{lev}[1e-{factor}]-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex1}{var}_{lev}hpa_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")

#EX2
fig=plt.figure()
ax=plt.axes(projection=cart_proj)
ax.set_extent(domain)
ax.contourf(
      lons,lats,ter2,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
      transform=ccrs.PlateCarree()
    )
 
contour_plus=ax.contour(
  lons,lats,mean2*10**factor,
  levels=levc_plus,
  colors="black",linewidths=2.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus, fontsize=12)

contour_minus=ax.contour(
  lons,lats,mean2*10**factor,
  levels=levc_minus,
  colors="black",linewidths=2.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_minus, fontsize=12)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex2}-{var}{lev}[1e-{factor}]-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex2}{var}_{lev}_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")


##Differemce
fig=plt.figure()
ax=plt.axes(projection=cart_proj)
ax.set_extent(domain)

contour_plus=ax.contour(
  lons,lats,dif*10**factor,
  levels=levc_plus,
  colors="red",linewidths=2.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus, fontsize=12)

contour_minus=ax.contour(
  lons,lats,dif*10**factor,
  levels=levc_minus,
  colors="blue",linewidths=2.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus, fontsize=12)

ax.contourf(
      lons,lats,ter1,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
      transform=ccrs.PlateCarree()
    )
plt.contourf(lons ,lats, sign_map,
             levels=[0.5,1.5], colors="none", hatches=["oo"],
             transform=ccrs.PlateCarree()
            )

# shade=ax.contourf(
#   lons,lats,dif*10**factor,
#   levels=levc_dif,
#   cmap="bwr",
#   extend="both",
#   transform=ccrs.PlateCarree()
#   )

# plt.colorbar(shade)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"Dif-{var}{lev}[1e-{factor}]-dx=dy={dx}")

plt.savefig(f"{fig_dir}/Dif_{var}_{lev}hpa_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")

print("End Program")