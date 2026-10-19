###別ファイルで作ったnetCDFファイルを使って指定した期間の平均を求める
from datetime import datetime
import calendar
from netCDF4 import Dataset
import xarray as xr
import numpy as np
from wrf import getvar,interplevel,get_cartopy,latlon_coords
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.draw_field import draw_grid 

#月を設定
month=2
month_name=calendar.month_abbr[month]

#高度指定
lev=750

##変数指定
#var="theta"
#var="dtheta_dt"
#var="Adv.theta_x"
#var="Adv.theta_y"
#var="Adv.theta_h"
#var="Adv.theta_p"
#var="Adv.f"
#var="Stretching"
#var="Tilting"
#var="Sum"
var="Res"
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
#factor=6
#factor=9
factor=10
#factor=11
#描画のための倍率(10のfactor乗)


levc_dif=np.arange(-8.8,10.2,1.6)
levc_plus=np.arange(0.6,16.6,1.6)
levc_minus=np.arange(-15,0,1.6)

#levc_plus=np.arange(0.2,16.3,3.2)
#levc_minus=np.arange(-16.2,0,3.2)


##座標系情報を得るために任意のwrfoutファイルを開く
sdate=datetime(2024,12,1,0)
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour

if n_member >= 1 :
 coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
 coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
else :
 coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
 coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

ter1=getvar(coord_ds1, "ter")
ter2=getvar(coord_ds2, "ter")
land=getvar(coord_ds1,"LANDMASK")
cart_proj=get_cartopy(land)
print("projection:",cart_proj)
lats,lons=latlon_coords(land)
#print(lats)
coord_ds1.close()
coord_ds2.close()


dfile1=f"{wrfout_dir}/{ex1}/n1/{ex1}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"
dfile2=f"{wrfout_dir}/{ex2}/n1/{ex2}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"

ds_all1=xr.open_dataset(dfile1)
ds_all2=xr.open_dataset(dfile2)

ds_all1=xr.decode_cf(ds_all1)
ds_all2=xr.decode_cf(ds_all2)

#指定の月のデータのみ取り出す
def select_month(ds_all, year, month):

    start = np.datetime64(
        f"{year:04d}-{month:02d}-01T00:00:00"
    )
    # 翌月
    if month == 12:
        next_year = year + 1
        next_month = 1
    else:
        next_year = year
        next_month = month + 1

    end = np.datetime64(
        f"{next_year:04d}-{next_month:02d}-01T00:00:00"
    )

    ds = ds_all.where(
        (ds_all["Time"] >= start) &
        (ds_all["Time"] < end),
        drop=True
    )

    print("----------------------------------------")
    print(f"Selected month : {year}-{month:02d}")
    print(f"Start time     : {ds['Time'].values[0]}")
    print(f"End time       : {ds['Time'].values[-1]}")
    print(f"Number of data : {ds.sizes['time']}")
    print("----------------------------------------")

    return ds


#print(ds1)
#print(ds2)

ds1 = select_month(ds_all1, 2025 if month != 12 else 2024, month)
ds2 = select_month(ds_all2, 2025 if month != 12 else 2024, month)


if var == "Adv.theta_h" :
  var1=ds1["Adv.theta_x"].mean(dim="time",skipna=True) + ds1["Adv.theta_y"].mean(dim="time",skipna=True)
  var2=ds2["Adv.theta_x"].mean(dim="time",skipna=True) + ds2["Adv.theta_y"].mean(dim="time",skipna=True)


elif var == "Sum"  or var == "Res":
  sum1=(
  ds1["Adv.theta_x"].mean(dim="time",skipna=True)
  +ds1["Adv.theta_y"].mean(dim="time",skipna=True)
  +ds1["Adv.theta_p"].mean(dim="time",skipna=True)
  +ds1["Adv.f"].mean(dim="time",skipna=True)
  +ds1["Stretching"].mean(dim="time",skipna=True)
  +ds1["Tilting"].mean(dim="time",skipna=True)
  )

  sum2=(
  ds2["Adv.theta_x"].mean(dim="time",skipna=True)
  +ds2["Adv.theta_y"].mean(dim="time",skipna=True)
  +ds2["Adv.theta_p"].mean(dim="time",skipna=True)
  +ds2["Adv.f"].mean(dim="time",skipna=True)
  +ds2["Stretching"].mean(dim="time",skipna=True)
  +ds2["Tilting"].mean(dim="time",skipna=True)
  )

  if var == "Sum":
    var1=sum1
    var2=sum2

  elif var == "Res":
    var1=ds1["dtheta_dt"].mean(dim="time",skipna=True) - sum1
    var2=ds2["dtheta_dt"].mean(dim="time",skipna=True) - sum2

else:
  var1=ds1[var].mean(dim="time",skipna=True)
  var2=ds2[var].mean(dim="time",skipna=True)
 
#print(var1)

ds1.close()
ds2.close()

dif=var1 - var2

##Plot 
##EX1
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
  lons,lats,var1*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,var1*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_minus)

draw_grid(ax)

ax.coastlines()
ax.set_title(f"{ex1}-{var}{lev}[1e-{factor}]-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{month_name}_{ex1}{var}_{lev}hpa_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")

##EX2
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
  lons,lats,var2*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,var2*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_minus)

draw_grid(ax)

ax.coastlines()
ax.set_title(f"{ex2}-{var}{lev}[1e-{factor}]-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{month_name}_{ex2}{var}_{lev}_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")


##Differemce
fig=plt.figure()
ax=plt.axes(projection=cart_proj)
ax.set_extent(domain)

contour_plus=ax.contour(
  lons,lats,dif*10**factor,
  levels=levc_plus,
  colors="red",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,dif*10**factor,
  levels=levc_minus,
  colors="blue",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

ax.contourf(
      lons,lats,ter1,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
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

plt.savefig(f"{fig_dir}/{month_name}_Dif_{var}_{lev}hpa_dxdy{dx}_{gtime}sigma{sigma}.png")

plt.close("all")

print("End Program")