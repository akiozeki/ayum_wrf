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
# var="dtheta_dx"
# varlabel="K/100km"
# factor=100000
# #単位変化のための係数
# vmin=-2.0
# vmax=2.0
# vint=0.2
# dif_vmin=-0.55
# dif_vmax=0.55
# dif_vint=0.1
# cmap="coolwarm"
# extend="both"

# var="dtheta_dy"
# varlabel="K/100km"
# factor=100000
# vmin=-2.0
# vmax=2.0
# vint=0.2
# dif_vmin=-0.55
# dif_vmax=0.55
# dif_vint=0.1
# cmap="coolwarm"
# extend="min"

# var="grad_theta"
# varlabel="K/100km"
# factor=100000
# vmin=0
# vmax=2.0
# vint=0.2
# dif_vmin=-0.7
# dif_vmax=0.7
# dif_vint=0.2
# cmap="coolwarm"
# extend="max"

# var="B.V._Frequency"
# #描画はその逆数
# varlabel="1e+2 s"
# factor=0.01
# vmin=0.4
# vmax=1.2
# vint=0.1
# dif_vmin=-0.35
# dif_vmax=0.35
# dif_vint=0.1
# cmap="coolwarm"
# extend="max"

var="Eady_Growth_Rate"
varlabel="/Day"
factor=86400
vmin=0
vmax=1.8
vint=0.2
dif_vmin=-0.7
dif_vmax=0.7
dif_vint=0.2
cmap="coolwarm"
extend="max"

p=850
mday=7

case="2025winter"
ex1="CTL"
ex2="ME00"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig/{case}/egr"

#アンサンブルメンバー数 / アンサンブル平均したい場合はそのサイズ
ensmean="yes"
#"yes" or "no" noの場合以下の要素数は必ず1個
target_members=[1, 2, 3, 4, 5]
#取り出すメンバ
N=len(target_members)

  
have_test="no"
#平均値の差の検定の有無 yes or no
#noは既に計算している場合も含む
if have_test == "yes":
  test="TTest_Nstar"
  print("test -> ", test)
else:
  test=""   
  print("No Test")
  if ensmean == "yes":
    sign_file=f"{wrfout_dir}/{ex2}/sign_map_{var}_mean_ens.npy"
  elif ensmean == "no":
    sign_file=f"{wrfout_dir}/{ex2}/n{target_members[0]}/sign_map_{var}_mean.npy"

sdate=datetime(2024,12,1,0)
edate=datetime(2025,2,28,18)
print(f"Set Time : {sdate}----{edate}")
#domain=0
domain=[120,150,30,45]

############################

##座標系情報を得るために任意のwrfoutファイルを開く
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour

#座標情報を取りだす
coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/n{target_members[0]}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/n{target_members[0]}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

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

if ensmean == "yes":
  dfile1=f"{wrfout_dir}/{ex1}/{ex1}_EGR_Seo{p}_{mday}dMean.nc"
  dfile2=f"{wrfout_dir}/{ex2}/{ex2}_EGR_Seo{p}_{mday}dMean.nc"
elif ensmean == "no":
  dfile1=f"{wrfout_dir}/{ex1}/n{target_members[0]}/{ex1}_EGR_Seo{p}_{mday}dMean.nc"
  dfile2=f"{wrfout_dir}/{ex2}/n{target_members[0]}/{ex2}_EGR_Seo{p}_{mday}dMean.nc"

ds1=xr.open_dataset(dfile1)

if var == "grad_theta":
  var1=np.sqrt(ds1["dtheta_dx"]**2 + ds1["dtheta_dy"]**2) * factor

elif var == "B.V._Frequency":
  var1=(1/ds1["B.V._Frequency"]) * factor

else:
  var1=ds1[f"{var}"] * factor

ds1.close()


ds2=xr.open_dataset(dfile2)

if var == "grad_theta":
  var2=np.sqrt(ds2["dtheta_dx"]**2 + ds2["dtheta_dy"]**2) * factor

elif var == "B.V._Frequency":
  var2=(1/ds2["B.V._Frequency"]) * factor
  var="revN"

else:
  var2=ds2[f"{var}"] * factor
 
#print(var2) 
ds2.close()

##ココで母平均の差の検定
if ensmean == "yes":
  from myFunction.ttest_Nstar import ttest_2Dmap_ensmean

  if have_test == "yes" :
    sign_map=ttest_2Dmap_ensmean(var1, var2, 0.95)
    sign_path=f"{wrfout_dir}/{ex2}/sign_map_{var}_mean_ens.npy"
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
    sign_path=f"{wrfout_dir}/{ex2}/sign_map_{var}_mean_ens.npy"
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
log2dif=np.log2(mean1) - np.log2(mean2)


##Plot
##Plot
#引数：lon, lat, field, projection, domain, cmap, vmin, vmax, vint, extend, fname, unit, fig_title
#EX1
draw_latlon_contourf(lons, lats, mean1, cart_proj, domain,
                     cmap, vmin, vmax, vint, extend, "", varlabel, f"{ex1} {var}")

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var}_{ex1}_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var}_{ex1}_mean(n={target_members[0]}).png")

plt.close("all")

#EX2
draw_latlon_contourf(lons, lats, mean2, cart_proj, domain,
                     cmap, vmin, vmax, vint, extend, "", varlabel, f"{ex2} {var}")
if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var}_{ex2}_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var}_{ex2}_mean(n={target_members[0]}).png")
plt.close("all")

#Dif
draw_latlon_contourf(lons, lats, dif, cart_proj, domain,
                     "seismic", dif_vmin, dif_vmax, dif_vint, "both", "", varlabel, f"Dif {var}")
# draw_latlon_contourf_onlysig(lons, lats, ens_dif, cart_proj, 
#                      "seismic", -19, 21, 2, "both", "Z", "m", f"Dif {var_label}", sign_map)
sign_map=np.load(sign_file)
plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
#plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["///"], transform=ccrs.PlateCarree())

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var}_Dif_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var}_Dif_mean(n={target_members[0]}).png")
plt.close("all")


#log2dif
draw_latlon_contourf(lons, lats, log2dif, cart_proj, domain,
                     "seismic", -1.5, 1.5, 0.2, "both", "", "", f"log2Dif {var}")
# draw_latlon_contourf_onlysig(lons, lats, ens_dif, cart_proj, 
#                      "seismic", -19, 21, 2, "both", "Z", "m", f"Dif {var_label}", sign_map)
# plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
# plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["///"], transform=ccrs.PlateCarree())

if ensmean == "yes":
  plt.savefig(f"{fig_dir}/{var}_log2Dif_mean_ens(N={N}).png")
elif ensmean == "no":  
  plt.savefig(f"{fig_dir}/{var}_log2Dif_mean(n={target_members[0]}).png")
plt.close("all")

print("End Program")