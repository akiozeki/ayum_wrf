###特定気圧面の変数の平均と実験差を計算，描画する **アンサンブル平均は別ファイル->
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

from myFunction.draw_field import draw_latlon_contourf
from myFunction.ttest_Nstar import ttest_2Dmap
#from myFunction.ttest import ttest_2Dmap

n=1
#アンサンブルメンバ, アンサンブルを取っていない場合は0

varn="z"
lev=750
var_unit="m"
domain=[120,140,32,50]
vmin=2000
vmax=2600
vint=50
dif_vmin=-18
dif_vmax=18
dif_vint=4


var_label="Z750"
case="2025winter"
ex1="CTL"
ex2="ME00"

have_test="no"
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

date=start_date

while date <= end_date:
  y=date.year
  m=date.month
  d=date.day
  h=date.hour

  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}"
  print(f"Read wrfout File {time_str}")

  if date == end_date :
    ds_ref=Dataset(f"{wrfout_dir1}/n1/wrfout_d01_{time_str}:00:00")
 #   ds_ref=Dataset(f"{wrfout_dir1}/wrfout_d01_{time_str}:00:00")
    ter=getvar(ds_ref,"ter")
    land=getvar(ds_ref,"LANDMASK")
    cart_proj=get_cartopy(land)
    print("projection:",cart_proj)
    lats,lons=latlon_coords(land)


  if n == 0:
    wrfout1=f"{wrfout_dir1}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/wrfout_d01_{time_str}:00:00"

  if n != 0:
    wrfout1=f"{wrfout_dir1}/n{n}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir2}/n{n}/wrfout_d01_{time_str}:00:00"
#ここで「n〇」とすれば任意のアンサンブルメンバーを1つ取り出せる


#ファイルオープンとデータ取得
  ds1=Dataset(wrfout1)
  ds2=Dataset(wrfout2)
  #変数取り出す
  p1=getvar(ds1,"p",units="hpa")
  if var_unit == "" :
    var1=getvar(ds1,varn)
  else :  
    var1=getvar(ds1,varn,units=var_unit)

  pvar1=interplevel(var1,p1,lev)
    
  p2=getvar(ds2,"p",units="hpa")
  if var_unit == "" :
    var2=getvar(ds2,varn)
  else :
    var2=getvar(ds2,varn,units=var_unit)

  pvar2=interplevel(var2,p2,lev)
    
    
  ds1.close()
  ds2.close()

  var_list1.append(pvar1)
  var_list2.append(pvar2)

  date+=timedelta(hours=dh)

darray1=xr.concat(var_list1,dim="time")
darray2=xr.concat(var_list2,dim="time")
#print(darray1)
#darray[time, south_north, west_east]

##ココで母平均の差の検定
if have_test == "yes":
  sign_map=ttest_2Dmap(darray1, darray2)

#平均値と差を計算
mean1=darray1.mean(dim="time", skipna=True)
mean2=darray2.mean(dim="time", skipna=True)
dif=mean1 - mean2


##Plot
draw_latlon_contourf(lons, lats, mean1, cart_proj, domain,
                     "viridis", vmin, vmax, vint, "both", varn, "", f"{ex1} {var_label}")
plt.savefig(f"{fig_dir}/{var_label}_{ex1}_mean_n={n}.png")
plt.close("all")

draw_latlon_contourf(lons, lats, mean2, cart_proj, domain,
                     "viridis", vmin, vmax, vint, "both", varn, "", f"{ex2} {var_label}")
plt.savefig(f"{fig_dir}/{var_label}_{ex2}_mean_n={n}.png")
plt.close("all")

draw_latlon_contourf(lons, lats, dif, cart_proj, domain,
                     "seismic", dif_vmin, dif_vmax, dif_vint, "both", varn, "", f"Dif {var_label}")
if have_test == "yes":
  plt.contourf(lons ,lats, sign_map, levels=[0.5,1.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
plt.savefig(f"{fig_dir}/{var_label}_Dif_mean_n={n}.png")
plt.close("all")