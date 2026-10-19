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
from scipy.stats import ttest_ind

from myFunction.draw_field import draw_latlon_contourf, draw_latlon_contourf_onlysig
from myFunction.ttest import ttest_2Dmap
#from myFunction.ttest_Nstar import ttest_2Dmap

#各パラメータを設定
#変数EGR or revN or Ushear
var="EGR"
var_unit="[ /Day ]"
factor=86400
vmin, vmax, vint=0, 2.8, 0.4
vmin_dif, vmax_dif, vint_dif=-0.9, 1.1, 0.2
#factorは単位変化のための係数

#var="N"
##描画はその逆数
#var_unit="[ 10^2 s ]"
#factor=0.01
#vmin, vmax, vint=0.5,1.2,0.1
#vmin_dif, vmax_dif, vint_dif=-0.25,0.35,0.1

#var="Ushear"
#var_unit="[ m/s /km ]"
#factor=1000
#vmin, vmax, vint=0,14,2
#vmin_dif, vmax_dif, vint_dif=-4.5,5.5,1

vmin_log2, vmax_log2, vint_log2=-1.4, 1.8, 0.4

p=850
mday=7
key=f"{var}{p}_{mday}dMean"

case="2025DJF"
ex1="CTL"
ex2="ME00"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig/{case}"

target_members=[1]
#計算するアンサンブルメンバn〇，calcのコードと対応
ensmean="no"
#"yes" or "no"

if ensmean == "yes":
  input_file1=f"{wrfout_dir}/{ex1}/{ex1}_{key}.nc"
  input_file2=f"{wrfout_dir}/{ex2}/{ex2}_{key}.nc"

elif ensmean == "no":
  n = target_members[0]
  input_file1=f"{wrfout_dir}/{ex1}/n{n}/{ex1}_{key}.nc"
  input_file2=f"{wrfout_dir}/{ex2}/n{n}/{ex2}_{key}.nc"

# os.makedirs(f"{fig_dir}/{ex1}/egr",exist_ok=True)
# os.makedirs(f"{fig_dir}/{ex2}/egr",exist_ok=True)

have_test="yes"
#平均値の差の検定の有無 yes or no

#作図の描画設定
domain=[125,140,35,44]
landmask="no"
"yes" or "no"

sdate=datetime(2024,12,1,0)
edate=datetime(2025,2,28,18)
time_step=timedelta(hours=6)
n_step=(edate - sdate) // time_step
print(f"Set Time : {sdate}----{edate},n_step={n_step+1}")
############################

##座標系情報を得るために任意のwrfoutファイルを開く

sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour

if ensmean == "no":
    coord_member = target_members[0]
    coord_file1 = (
        f"{wrfout_dir}/{ex1}/n{coord_member}/"
        f"wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00"
    )
    coord_file2 = (
        f"{wrfout_dir}/{ex2}/n{coord_member}/"
        f"wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00"
    )
else:
    coord_file1 = (
        f"{wrfout_dir}/{ex1}/"
        f"wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00"
    )
    coord_file2 = (
        f"{wrfout_dir}/{ex2}/"
        f"wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00"
    )

coord_ds = Dataset(coord_file1)
land = getvar(coord_ds, "LANDMASK")
cart_proj = get_cartopy(land)
ter1 = getvar(coord_ds, "ter", units="m")
ter1 = np.ma.masked_where(land != 1, ter1)
coord_ds.close()

coord_ds2 = Dataset(coord_file2)
land2 = getvar(coord_ds2, "LANDMASK")
ter2 = getvar(coord_ds2, "ter", units="m")
ter2 = np.ma.masked_where(land2 != 1, ter2)
coord_ds2.close()


ds1=xr.open_dataset(input_file1)
var1=ds1[f"{var}"]

if var == "N":
 var1=1/var1

lats,lons=latlon_coords(var1)
ds1.close()

ds2=xr.open_dataset(input_file2)
var2=ds2[f"{var}"]
if var == "N":
 var2=1/var2
 var="revN"
#print(var2) 
ds2.close()

if ensmean == "yes":
  var1=var1.mean(dim="member",skipna=True)
  var2=var2.mean(dim="member",skipna=True)

##ココで母平均の差の検定
if have_test == "yes":
  if ensmean == "no":
    sign_map=ttest_2Dmap(var1, var2)
  elif ensmean == "yes":
    sign_map=ttest_2Dmap(var1, var2)

##平均計算
#時間方向へ平均(member,ny,nx)
mean1=var1.mean(dim="time",skipna=True)
mean2=var2.mean(dim="time",skipna=True)


##陸面マスク
mean1=np.ma.masked_where(land == 1,mean1)
mean2=np.ma.masked_where(land == 1,mean2)
#sign=np.ma.masked_where(land == 1,sign)

#print(mean1)
#print(mean2)
dif=mean1 - mean2
log2dif=np.log2(mean1) - np.log2(mean2)

key=f"{var}{p}-{mday}dMean"

##Plot
#EX1
draw_latlon_contourf(lons, lats, mean1, cart_proj, 
                     "Reds", vmin, vmax, vint, "max", var, var_unit, f"{ex1} {var} {var_unit}")
#plt.scatter(lon_point,lat_point,marker="^",s=500,color="green",transform=ccrs.PlateCarree())
plt.savefig(f"{fig_dir}/{ex1}Mean_{key}.png")
plt.close("all")

#EX2
draw_latlon_contourf(lons, lats, mean2, cart_proj, 
                     "Reds", vmin, vmax, vint, "max", var, var_unit, f"{ex2} {var} {var_unit}")
#plt.scatter(lon_point,lat_point,marker="^",s=500,color="green",transform=ccrs.PlateCarree())
plt.savefig(f"{fig_dir}/{ex2}Mean_{key}.png")
plt.close("all")

#Dif
draw_latlon_contourf_onlysig(lons, lats, dif, cart_proj, 
                     "seismic", vmin_dif, vmax_dif, vint_dif, "both", var, var_unit, f"Dif {var} {var_unit}", sign_map)
#plt.plot([122,137],[40.5,40.5],color="blue",linewidth=3.0,transform=ccrs.PlateCarree())
#plt.scatter(lon_point,lat_point,marker="^",s=500,color="green",transform=ccrs.PlateCarree())
#plt.contourf(lons ,lats, sign_map, levels=[-0.5,0.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())
#plt.contourf(lons ,lats, sign_map, levels=[-0.5,0.5], colors="none", hatches=["oo"], transform=ccrs.PlateCarree())

plt.savefig(f"{fig_dir}/Dif_Mean_{key}.png")
plt.close("all")

#Log2Dif
draw_latlon_contourf_onlysig(lons, lats, log2dif, cart_proj, 
                     "seismic", vmin_log2, vmax_log2, vint_log2, "both", var, var_unit, f"log2Dif {var} {var_unit}", sign_map)
#plt.plot([122,137],[40.5,40.5],color="blue",linewidth=3.0,transform=ccrs.PlateCarree())
#plt.scatter(lon_point,lat_point,marker="^",s=500,color="green",transform=ccrs.PlateCarree())
plt.savefig(f"{fig_dir}/log2Dif_Mean_{key}.png")
plt.close("all")

print("End Program")