###別ファイルで作ったnetCDFファイルを使って指定した期間の平均を求める
import os
from datetime import datetime,timedelta
from netCDF4 import Dataset
from wrf import getvar,get_cartopy,latlon_coords,geo_bounds,interplevel,ALL_TIMES
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.ticker import FixedLocator

import sys
from pathlib import Path
#scriptディレクトリを検索パスに追加
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from myFunction.draw_field import draw_latlon_contourf, draw_latlon_contourf_onlysig


####################################
#各パラメータを設定
var_name="vT_bpass"
var_unit="m/s * K"
vmin=-12
vmax=24
vint=4
dif_vmin=-7.5
dif_vmax=7.5
dif_vint=1
cmap="coolwarm"
extend="both"

# var_name="0.5VV_bpass"
# var_unit="(m/s)^2"
# vmin=0
# vmax=48
# vint=6
# dif_vmin=-15
# dif_vmax=15
# dif_vint=2
# cmap="coolwarm"
# extend="max"

p=850
mday_high=2
mday_low=7

case="2025winter"
ex1="CTL"
ex2="ME00"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout"
fig_dir=f"{wrf_dir}/fig/{case}/disturb"

#アンサンブルメンバー数 / アンサンブル平均したい場合はそのサイズ
ensmean="yes"
#"yes" or "no" noの場合以下の要素数は必ず1個
target_members=[1, 2, 3, 4, 5]
#取り出すメンバ
N=len(target_members)

if ensmean == "yes":
  input_dir1=f"{wrfout_dir}/{case}/{ex1}"
  input_dir2=f"{wrfout_dir}/{case}/{ex2}"

elif ensmean == "no":
  input_dir1=f"{wrfout_dir}/{case}/{ex1}/n{target_members[0]}"
  input_dir2=f"{wrfout_dir}/{case}/{ex2}/n{target_members[0]}"
  
have_test="yes"
#平均値の差の検定の有無 yes or no
#noは既に計算している場合も含む
if have_test == "yes":
  test="TTest_Nstar"
  print("test -> ", test)
else:
  test=""   
  print("No Test")
  if ensmean == "yes":
    sign_file=f"{wrfout_dir}/{case}/{ex2}/sign_map_{var_name}_mean_ens(N={N}).nc"
  elif ensmean == "no":
    sign_file=f"{wrfout_dir}/{case}/{ex2}/n{target_members[0]}/sign_map_{var_name}_mean.nc"

sdate=datetime(2024,12,1,0)
edate=datetime(2025,2,28,18)
print(f"Used Data : {sdate}----{edate}")
domain=0
#domain=[120,150,30,45]

############################

##座標系情報を得るために任意のwrfoutファイルを開く
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour

#座標情報を取りだす
coord_ds1=Dataset(f"{wrfout_dir}/{case}/{ex1}/n{target_members[0]}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
coord_ds2=Dataset(f"{wrfout_dir}/{case}/{ex2}/n{target_members[0]}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

land1=getvar(coord_ds1,"LANDMASK")
ter1=getvar(coord_ds1,"ter",units="m")
ter1=np.ma.masked_where(land1 != 1,ter1)
land2=getvar(coord_ds2,"LANDMASK")
ter2=getvar(coord_ds2,"ter",units="m")
ter2=np.ma.masked_where(land2 != 1,ter2)

cart_proj=get_cartopy(land1)
print("projection:",cart_proj)
lats,lons=latlon_coords(land1)
coord_ds1.close()
coord_ds2.close()

##データの取り出しと計算
if var_name == "vT_bpass":
  vfile1_high=f"{input_dir1}/{ex1}_va{p}_{mday_high}dMean.nc"
  vfile2_high=f"{input_dir2}/{ex2}_va{p}_{mday_high}dMean.nc"
  Tfile1_high=f"{input_dir1}/{ex1}_tk{p}_{mday_high}dMean.nc"
  Tfile2_high=f"{input_dir2}/{ex2}_tk{p}_{mday_high}dMean.nc"

  v_bar1_high=xr.open_dataset(vfile1_high)["va_bar"]
  v_bar2_high=xr.open_dataset(vfile2_high)["va_bar"]
  T_bar1_high=xr.open_dataset(Tfile1_high)["tk_bar"]
  T_bar2_high=xr.open_dataset(Tfile2_high)["tk_bar"]

  vfile1_low=f"{input_dir1}/{ex1}_va{p}_{mday_low}dMean.nc"
  vfile2_low=f"{input_dir2}/{ex2}_va{p}_{mday_low}dMean.nc"
  Tfile1_low=f"{input_dir1}/{ex1}_tk{p}_{mday_low}dMean.nc"
  Tfile2_low=f"{input_dir2}/{ex2}_tk{p}_{mday_low}dMean.nc"

  v_bar1_low=xr.open_dataset(vfile1_low)["va_bar"]
  v_bar2_low=xr.open_dataset(vfile2_low)["va_bar"]
  T_bar1_low=xr.open_dataset(Tfile1_low)["tk_bar"]
  T_bar2_low=xr.open_dataset(Tfile2_low)["tk_bar"]

  v_bar1_high, v_bar1_low=xr.align(v_bar1_high, v_bar1_low, join="inner")
  v_bar2_high, v_bar2_low=xr.align(v_bar2_high, v_bar2_low, join="inner")
  T_bar1_high, T_bar1_low=xr.align(T_bar1_high, T_bar1_low, join="inner")
  T_bar2_high, T_bar2_low=xr.align(T_bar2_high, T_bar2_low, join="inner")
#移動平均の広さによって時間次元数が異なるので小さいほうにそろえる

  var1=(v_bar1_high - v_bar1_low) * (T_bar1_high - T_bar1_low)
  var2=(v_bar2_high - v_bar2_low) * (T_bar2_high - T_bar2_low)

  
elif var_name == "0.5VV_bpass":
  ufile1_high=f"{input_dir1}/{ex1}_ua{p}_{mday_high}dMean.nc"
  ufile2_high=f"{input_dir2}/{ex2}_ua{p}_{mday_high}dMean.nc"
  vfile1_high=f"{input_dir1}/{ex1}_va{p}_{mday_high}dMean.nc"
  vfile2_high=f"{input_dir2}/{ex2}_va{p}_{mday_high}dMean.nc"

  u_bar1_high=xr.open_dataset(ufile1_high)["ua_bar"]
  u_bar2_high=xr.open_dataset(ufile2_high)["ua_bar"]
  v_bar1_high=xr.open_dataset(vfile1_high)["va_bar"]
  v_bar2_high=xr.open_dataset(vfile2_high)["va_bar"]

  ufile1_low=f"{input_dir1}/{ex1}_ua{p}_{mday_low}dMean.nc"
  ufile2_low=f"{input_dir2}/{ex2}_ua{p}_{mday_low}dMean.nc"
  vfile1_low=f"{input_dir1}/{ex1}_va{p}_{mday_low}dMean.nc"
  vfile2_low=f"{input_dir2}/{ex2}_va{p}_{mday_low}dMean.nc"

  u_bar1_low=xr.open_dataset(ufile1_low)["ua_bar"]
  u_bar2_low=xr.open_dataset(ufile2_low)["ua_bar"]
  v_bar1_low=xr.open_dataset(vfile1_low)["va_bar"]
  v_bar2_low=xr.open_dataset(vfile2_low)["va_bar"]

  u_bar1_high, u_bar1_low=xr.align(u_bar1_high, u_bar1_low, join="inner")
  u_bar2_high, u_bar2_low=xr.align(u_bar2_high, u_bar2_low, join="inner")
  v_bar1_high, v_bar1_low=xr.align(v_bar1_high, v_bar1_low, join="inner")
  v_bar2_high, v_bar2_low=xr.align(v_bar2_high, v_bar2_low, join="inner")
#移動平均の広さによって時間次元数が異なるので小さいほうにそろえる

  var1=0.5 * ((u_bar1_high - u_bar1_low)**2 + (v_bar1_high - v_bar1_low)**2)
  var1=0.5 * ((u_bar2_high - u_bar2_low)**2 + (v_bar2_high - v_bar2_low)**2)

else:
  raise ValueError("Set Valid Var")    


##ココで母平均の差の検定
if ensmean == "yes":
  from myFunction.ttest_Nstar import ttest_2Dmap_ensmean

  if have_test == "yes" :
    sign_map=ttest_2Dmap_ensmean(var1, var2, 0.95)
    sign_path=f"{wrfout_dir}/{case}/{ex2}/sign_map_{var_name}_mean_ens.npy"
    np.save(sign_path, sign_map)
    print("Save -> ",sign_path)

  
  ens_mean1=var1.mean(dim="member", skipna=True)
  ens_mean2=var2.mean(dim="member", skipna=True)
  mean1=ens_mean1.mean(dim="time", skipna=True)
  mean2=ens_mean2.mean(dim="time", skipna=True)
#メンバーおよび時間方向に平均   

elif ensmean == "no":
  from myFunction.ttest_Nstar import ttest_2Dmap

  if have_test == "yes" :
    sign_map=ttest_2Dmap(var1, var2, 0.95)
    sign_path=f"{wrfout_dir}/{case}/{ex2}/sign_map_{var_name}_mean_ens.npy"
    np.save(sign_path, sign_map)
    print("Save -> ",sign_path)

  mean1=var1.mean(dim="time", skipna=True)
  mean2=var2.mean(dim="time", skipna=True)
#時間方向に平均

##陸面マスク
# mean1=np.ma.masked_where(land1 == 1,mean1)
# mean2=np.ma.masked_where(land2 == 1,mean2)
#sign_map=np.ma.masked_where(land1 == 1,sign_map)

dif=mean1 - mean2

##Plot
##Plot
#引数：lon, lat, field, projection, domain, cmap, vmin, vmax, vint, extend, fname, unit, fig_title
#EX1
draw_latlon_contourf(lons, lats, mean1, cart_proj, domain,
                     cmap, vmin, vmax, vint, extend, "", var_unit, f"{ex1} {var_name}")

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var_name}_{ex1}_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var_name}_{ex1}_mean(n={target_members[0]}).png")

plt.close("all")

#EX2
draw_latlon_contourf(lons, lats, mean2, cart_proj, domain,
                     cmap, vmin, vmax, vint, extend, "", var_unit, f"{ex2} {var_name}")

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var_name}_{ex2}_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var_name}_{ex2}_mean(n={target_members[0]}).png")
plt.close("all")

#Dif
draw_latlon_contourf(lons, lats, dif, cart_proj, domain,
                     "seismic", dif_vmin, dif_vmax, dif_vint, "both", "", var_unit, f"Dif {var_name}")
# draw_latlon_contourf_onlysig(lons, lats, ens_dif, cart_proj, 
#                      "seismic", -19, 21, 2, "both", "Z", "m", f"Dif {var_label}", sign_map)
plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
#plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["///"], transform=ccrs.PlateCarree())

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var_name}_Dif_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var_name}_Dif_mean(n={target_members[0]}).png")
plt.close("all")

print("End Program")