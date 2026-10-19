###特定気圧面の変数の平均と実験差を計算，描画する
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

from myFunction.draw_field import draw_latlon_contourf, draw_latlon_contourf_onlysig
from myFunction.ttest import ttest_2Dmap
from myFunction.ttest_Nstar import ttest_2Dmap_ensmean


N=5
print("N (Ensemble Member) = ",N)
#アンサンブル数

var="z"
lev=850
var_label="Z850"
var_unit="m"
domain=0
#domain=[120,140,32,50]
vmin=1000
vmax=1600
vint=50
dif_vmin=-18
dif_vmax=18
dif_vint=4

case="2025winter"
ex1="CTL"
ex2="ME00"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir1=f"{wrf_dir}/wrfout/{case}/{ex1}"
wrfout_dir2=f"{wrf_dir}/wrfout/{case}/{ex2}"
fig_dir=f"{wrf_dir}/fig/{case}"

start_date=datetime(2024, 12, 1, 0)
end_date=datetime(2025, 2, 28, 18)
dh=3
print(f"Set Time : {start_date}----{end_date}")


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
  sign_file=f"{wrfout_dir2}/sign_map_{var_label}_mean_ens(N={N}).nc"
  sign_map=np.load(sign_file)

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

    if n == N and date == end_date :
      ds_ref=Dataset(f"{wrfout_dir1}/n1/wrfout_d01_{time_str}:00:00")
 #     ds_ref=Dataset(f"{wrfout_dir1}/wrfout_d01_{time_str}:00:00")
      ter=getvar(ds_ref,"ter")
      land=getvar(ds_ref,"LANDMASK")
      cart_proj=get_cartopy(land)
      print("projection:",cart_proj)
      lats,lons=latlon_coords(land)

#ファイルオープンとデータ取得
    wrfout1=f"{wrfout_dir1}/n{n}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/n{n}/wrfout_d01_{time_str}:00:00"
    ds1=Dataset(wrfout1)
    ds2=Dataset(wrfout2)
  #変数取り出す
    p1=getvar(ds1,"p",units="hpa")
    var1=getvar(ds1,var,units=var_unit)
    pvar1=interplevel(var1,p1,lev)
    
    p2=getvar(ds2,"p",units="hpa")
    var2=getvar(ds2,var,units=var_unit)
    pvar2=interplevel(var2,p2,lev)
    
    
    ds1.close()
    ds2.close()

    var_list1.append(pvar1)
    var_list2.append(pvar2)

    date+=timedelta(hours=dh)

  darray1=xr.concat(var_list1,dim="time")
  darray2=xr.concat(var_list2,dim="time")

  if test == "TTest_Nstar":
    ens_list1.append(darray1)
    ens_list2.append(darray2)
  
  else: 
    mean1=darray1.mean(dim="time", skipna=True)
    mean2=darray2.mean(dim="time", skipna=True)

    mean_list1.append(mean1)
    mean_list2.append(mean2)

  print("************************************")
  n+=1


if test == "TTest_Nstar":
  enstime_darray1=xr.concat(ens_list1,dim="member")
  enstime_darray2=xr.concat(ens_list2,dim="member")
#4次元配列(member, time, south_north, west_east)

else:
  ens_darray1=xr.concat(mean_list1,dim="member")
  ens_darray2=xr.concat(mean_list2,dim="member")
#3次元配列(member, south_north, west_east)
#print(ens_darray1)


##ココで母平均の差の検定
if have_test == "yes" :
  if use_Nstar == "yes":
    sign_map=ttest_2Dmap_ensmean(enstime_darray1, enstime_darray2, 0.95)
    sign_path=f"{wrfout_dir2}/sign_map_{var_label}_mean_ens.npy"
    np.save(sign_path, sign_map)
    print("Save -> ",sign_path)
    ens_darray1=enstime_darray1.mean(dim="time", skipna=True)
    ens_darray2=enstime_darray2.mean(dim="time", skipna=True)
#時間方向に圧縮   
 
  if use_Nstar == "no":
    sign_map=ttest_2Dmap(ens_darray1, ens_darray2, 0.95)

#平均値と差を計算
ens_mean1=ens_darray1.mean(dim="member", skipna=True)
ens_mean2=ens_darray2.mean(dim="member", skipna=True)
ens_dif=ens_mean1 - ens_mean2

##Plot
#引数：lon, lat, field, projection, domain, cmap, vmin, vmax, vint, extend, fname, unit, fig_title
draw_latlon_contourf(lons, lats, ens_mean1, cart_proj, domain,
                     "viridis", vmin, vmax, vint, "max", "Z", "m", f"{ex1} {var_label}")
plt.savefig(f"{fig_dir}/{var_label}_{ex1}_mean_ens(N={N}).png")
plt.close("all")

draw_latlon_contourf(lons, lats, ens_mean2, cart_proj, domain,
                     "viridis", vmin, vmax, vint, "max", "Z", "m", f"{ex2} {var_label}")
plt.savefig(f"{fig_dir}/{var_label}_{ex2}_mean_ens(N={N}).png")
plt.close("all")

draw_latlon_contourf(lons, lats, ens_dif, cart_proj, domain,
                     "seismic", dif_vmin, dif_vmax, dif_vint, "both", "Z", "m", f"Dif {var_label}")
# draw_latlon_contourf_onlysig(lons, lats, ens_dif, cart_proj, 
#                      "seismic", -19, 21, 2, "both", "Z", "m", f"Dif {var_label}", sign_map)
plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
#plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["///"], transform=ccrs.PlateCarree())

plt.savefig(f"{fig_dir}/{var_label}_Dif_mean_ens(N={N})_{test}.png")
plt.close("all")