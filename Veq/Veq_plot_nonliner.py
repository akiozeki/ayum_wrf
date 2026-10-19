###別ファイルで作ったnetCDFファイルを使って指定した期間の平均を求める
from datetime import datetime
from netCDF4 import Dataset
import xarray as xr
import numpy as np
from wrf import getvar,interplevel,get_cartopy,latlon_coords
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.draw_field import draw_grid

def gfilter(var,gs):
  var_smooth=xr.DataArray(
    gaussian_filter(var,sigma=(gs,gs)),
    coords=var.coords,
    dims=var.dims
  )
  return(var_smooth)


#高度
lev=750

##変数指定
#var="Adv.zeta_h"
var="Stretching"
#var="Adv.p"未実装
#var="Tilting"未実装


case="2025winter"
ex1="CTL"
ex2="ME00"

n_member=1
dx=25
dy=25
#km単位
sigma=6
gtime="end"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig/{case}/Veq"

domain=[120,140,32,50]
#図の描画範囲
#factor=11
factor=10

levc_dif=np.arange(-8.8,10.2,1.6)
levc_plus=np.arange(0.6,16.6,1.6)
levc_minus=np.arange(-15,0,1.6)

sdate=datetime(2024,12,1,0)


##座標系情報を得るために任意のwrfoutファイルを開く
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

omega=7.292*10**-5
lat=getvar(coord_ds1, "lat")
f=2*omega*np.sin(np.radians(lat))
f=f.assign_coords({
    "west_east": np.arange(lat.west_east.size) * dx*1000,
    "south_north": np.arange(lat.south_north.size) * dy*1000
  })  
    
ter1=getvar(coord_ds1, "ter")
ter2=getvar(coord_ds2, "ter")
land=getvar(coord_ds1,"LANDMASK")
cart_proj=get_cartopy(land)
#print("projection:",cart_proj)
lats,lons=latlon_coords(land)
#print(lats)
coord_ds1.close()
coord_ds2.close()


dfile1=f"{wrfout_dir}/{ex1}/n1/{ex1}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"
dfile2=f"{wrfout_dir}/{ex2}/n1/{ex2}Voleq{lev}hPa_dxdy{dx}_{gtime}sigma{sigma}.nc"

ds1=xr.open_dataset(dfile1)
ds2=xr.open_dataset(dfile2)

if gtime == "start":
  u1=ds1["u_smooth"]
  u2=ds2["u_smooth"]
  v1=ds1["v_smooth"]
  v2=ds2["v_smooth"]
  # omg1=ds1["omg_smooth"]
  # omg2=ds2["omg_smooth"]
  zeta1=ds1["zeta_smooth"]
  zeta2=ds2["zeta_smooth"]

if gtime == "end":
  u1=ds1["u"]
  u2=ds2["u"]
  v1=ds1["v"]
  v2=ds2["v"]
  zeta1=ds1["zeta"]
  zeta2=ds2["zeta"]

u_mean1=u1.mean(dim="time", skipna=True)
u_mean2=u2.mean(dim="time", skipna=True)
v_mean1=v1.mean(dim="time", skipna=True)
v_mean2=v2.mean(dim="time", skipna=True)
zeta_mean1=zeta1.mean(dim="time", skipna=True)
zeta_mean2=zeta2.mean(dim="time", skipna=True)

u_ano1=u1 - u_mean1
u_ano2=u2 - u_mean2
v_ano1=v1 - v_mean1
v_ano2=v2 - v_mean2
zeta_ano1=zeta1 - zeta_mean1
zeta_ano2=zeta2 - zeta_mean2

# print(u_mean1)
# print(u_ano1)
# print(theta_mean1)
# print(theta_ano1)

if var == "Adv.zeta_h":
   stationary1=-(u_mean1 * zeta_mean1.differentiate("west_east") + v_mean1 * zeta_mean1.differentiate("south_north"))
   stationary2=-(u_mean2 * zeta_mean2.differentiate("west_east") + v_mean2 * zeta_mean2.differentiate("south_north"))
   non_stationary1=-(u_ano1 * zeta_ano1.differentiate("west_east") + (v_ano1 * zeta_ano1.differentiate("south_north"))).mean(dim="time", skipna=True)
   non_stationary2=-(u_ano2 * zeta_ano2.differentiate("west_east") + (v_ano2 * zeta_ano2.differentiate("south_north"))).mean(dim="time", skipna=True)

elif var == "Stretching":
  stationary1=-(zeta_mean1 + f)*(u_mean1.differentiate("west_east") + v_mean1.differentiate("south_north"))
  stationary2=-(zeta_mean2 + f)*(u_mean2.differentiate("west_east") + v_mean2.differentiate("south_north"))
  non_stationary1=-(zeta_ano1*(u_ano1.differentiate("west_east")+v_ano1.differentiate("south_north"))).mean(dim="time", skipna=True)
  non_stationary2=-(zeta_ano2*(u_ano2.differentiate("west_east")+v_ano2.differentiate("south_north"))).mean(dim="time", skipna=True)

if gtime == "end":
  stationary1=gfilter(stationary1, sigma)
  stationary2=gfilter(stationary2, sigma)
  non_stationary1=gfilter(non_stationary1, sigma)
  non_stationary2=gfilter(non_stationary2, sigma)

ds1.close()
ds2.close()

stationary_dif=stationary1 - stationary2
non_stationary_dif=non_stationary1 - non_stationary2

# print(stationary_dif)
# print(non_stationary_dif)


# if var == "Adv.theta_h":
#     adv_theta_mean1 = ds1["Adv.theta_h"].mean(dim="time")
#     residual = adv_theta_mean1 - (stationary1 + non_stationary1)
#     print("Residual max =", np.nanmax(np.abs(residual)))

if var == "Stretching":
    stretch_mean1 = ds1["Stretching"].mean(dim="time")
    residual = stretch_mean1 - (stationary1 + non_stationary1)
    print("Residual max =", np.nanmax(np.abs(residual)))


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
  lons,lats,stationary1*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,stationary1*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex1}-Stationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex1}_Stationary{var}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

plt.close("all")


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
  lons,lats,non_stationary1*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,non_stationary1*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex1}-nonStationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex1}_nonStationary_{var}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

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
  lons,lats,stationary2*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,stationary2*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex2}-Stationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex2}_Stationary{var}_{lev}_dxdy{dx}_sigma{sigma}.png")

plt.close("all")


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
  lons,lats,non_stationary2*10**factor,
  levels=levc_plus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,non_stationary2*10**factor,
  levels=levc_minus,
  colors="black",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex2}-nonStationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/{ex2}_nonStationary{var}_{lev}_dxdy{dx}_sigma{sigma}.png")

plt.close("all")


##Differemce
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
  lons,lats,stationary_dif*10**factor,
  levels=levc_plus,
  colors="red",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,stationary_dif*10**factor,
  levels=levc_minus,
  colors="blue",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"Dif-Stationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}//Dif_Stationary_{var}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

plt.close("all")


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
  lons,lats,non_stationary_dif*10**factor,
  levels=levc_plus,
  colors="red",linewidths=1.0,
  transform=ccrs.PlateCarree()
  )
ax.clabel(contour_plus)

contour_minus=ax.contour(
  lons,lats,non_stationary_dif*10**factor,
  levels=levc_minus,
  colors="blue",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
ax.clabel(contour_minus)

draw_grid(ax)
ax.coastlines()
ax.set_title(f"Dif-nonStationary-{var}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}//Dif_nonStationary_{var}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

plt.close("all")

print("End Program")