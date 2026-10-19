###別ファイルで作ったnetCDFファイルを使って指定した期間の平均を求める
from datetime import datetime
from netCDF4 import Dataset
import xarray as xr
import numpy as np
from wrf import getvar,get_cartopy,latlon_coords
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.draw_field import draw_grid


##変数指定
#shading
#svar="thg"
#svar="dthg_dt"
#svar="Adv.thgx"
#svar="Adv.thgy"
# svar="Adv.thgh"
# svar="Adv.f"
# svar="Stretching"
# svar="Str_res"
#svar="Sum"
svar=0

#contour
#cvar="Z"
#cvar="thg"
cvar="dzetag_dt"
# cvar="Adv.zetagx"
# cvar="Adv.zetagy"
#cvar="Adv.zetagh"
#cvar="Adv.f"
#cvar="Stretching"
#cvar="Res"
#cvar="Sum"
#cvar=0

#vector
#uvar="Ug"
#vvar="Vg"
#uvar="Ua"
#vvar="Va"
uvar=0
vvar=0
#0だとベクトル地衡風の描画なし

##実験，計算設定
case="2025winter"
ex1="CTL"
ex2="ME00"

n_member=1
dx=25
lev=750
sigma=6
gtime="end"

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig"
input_dir=wrfout_dir

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

sdate=datetime(2024,12,1,0)


##座標系情報を得るために任意のwrfoutファイルを開く
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour



coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

land=getvar(coord_ds1,"LANDMASK")
cart_proj=get_cartopy(land)
ter1=getvar(coord_ds1,"ter")
ter2=getvar(coord_ds2,"ter")
print("projection:",cart_proj)
lats,lons=latlon_coords(land)
#print(lats)
coord_ds1.close()
coord_ds2.close()


dfile1=f"{input_dir}/{ex1}/n1/{ex1}{lev}hPa_QGVeq_dxdy{dx}_{gtime}sigma{sigma}.nc"
dfile2=f"{input_dir}/{ex2}/n1/{ex2}{lev}hPa_QGVeq_dxdy{dx}_{gtime}sigma{sigma}.nc"

ds1=xr.open_dataset(dfile1)
ds2=xr.open_dataset(dfile2)

#print(ds1)
#print(ds2)
if cvar == 0 :
  cvar1=cvar2=0
  print("No Shading Var")

elif cvar == "Adv.zetagh" :
  cvar1=ds1["Adv.zetagx"].mean(dim="time",skipna=True) + ds1["Adv.zetagy"].mean(dim="time",skipna=True)
  cvar2=ds2["Adv.zetagx"].mean(dim="time",skipna=True) + ds2["Adv.zetagy"].mean(dim="time",skipna=True)


elif cvar == "Sum" :
  cvar1=ds1["Adv.zetagx"].mean(dim="time",skipna=True) + ds1["Adv.zetagy"].mean(dim="time",skipna=True) + ds1["Adv.f"].mean(dim="time",skipna=True) + ds1["Stretching"].mean(dim="time",skipna=True)
  cvar2=ds2["Adv.zetagx"].mean(dim="time",skipna=True) + ds2["Adv.zetagy"].mean(dim="time",skipna=True) + ds2["Adv.f"].mean(dim="time",skipna=True) + ds2["Stretching"].mean(dim="time",skipna=True)

elif cvar == "Res" :
  cvar1=ds1["dzetag_dt"].mean(dim="time",skipna=True) - ds1["Adv.zetagx"].mean(dim="time",skipna=True) - ds1["Adv.zetagy"].mean(dim="time",skipna=True) - ds1["Adv.f"].mean(dim="time",skipna=True) - ds1["Stretching"].mean(dim="time",skipna=True)
  cvar2=ds2["dzetag_dt"].mean(dim="time",skipna=True) - ds2["Adv.zetagx"].mean(dim="time",skipna=True) - ds2["Adv.zetagy"].mean(dim="time",skipna=True) - ds2["Adv.f"].mean(dim="time",skipna=True) - ds2["Stretching"].mean(dim="time",skipna=True)

else:
  cvar1=ds1[cvar].mean(dim="time",skipna=True)
  cvar2=ds2[cvar].mean(dim="time",skipna=True)

cdif=cvar1 - cvar2

# if cvar == 0 :
#   print("No Contour Var")

# elif cvar == "Adv.thgh" :
#   cvar1=ds1["Adv.thgx"].mean(dim="time",skipna=True) + ds1["Adv.thgy"].mean(dim="time",skipna=True)
#   cvar2=ds2["Adv.thgx"].mean(dim="time",skipna=True) + ds2["Adv.thgy"].mean(dim="time",skipna=True)
#   cdif=cvar1 - cvar2

# elif cvar == "Str_res" :
#   cvar1=0 - ds1["Adv.thgx"].mean(dim="time",skipna=True) - ds1["Adv.thgy"].mean(dim="time",skipna=True) - ds1["Adv.f"].mean(dim="time",skipna=True)
#   cvar2=0 - ds2["Adv.thgx"].mean(dim="time",skipna=True) - ds2["Adv.thgy"].mean(dim="time",skipna=True) - ds2["Adv.f"].mean(dim="time",skipna=True)
#   cdif=cvar1 - cvar2

# else :
#    cvar1=ds1[cvar].mean(dim="time",skipna=True)
#    cvar2=ds2[cvar].mean(dim="time",skipna=True)
#    cdif=cvar1 - cvar2
    
if vvar == 0 :
  print("No Vector Var")

else :
  u1=ds1[uvar].mean(dim="time",skipna=True)
  v1=ds1[vvar].mean(dim="time",skipna=True)
  u2=ds2[uvar].mean(dim="time",skipna=True)
  v2=ds2[vvar].mean(dim="time",skipna=True)
  udif=u1 -u2
  vdif=v1 -v2

ds1.close()
ds2.close()

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

# if svar != 0 :
#   shade=ax.contourf(
#     lons,lats,svar1,
#     levels=levs,
#     cmap="bwr",
#     extend="both",
#     transform=ccrs.PlateCarree()
#     )

#   cbar=plt.colorbar(
#     shade,orientation="horizontal",
#     shrink=0.6,
#     aspect=25
#     )
#   cbar.ax.tick_params(labelsize=10)

if cvar != 0:
  contour_plus=ax.contour(
    lons,lats,cvar1*10**factor,
    levels=levc_plus,
    colors="black",linewidths=1.0,
    transform=ccrs.PlateCarree()
  )
  ax.clabel(contour_plus)

  contour_minus=ax.contour(
    lons,lats,cvar1*10**factor,
    levels=levc_minus,
    colors="black",linewidths=1.0,
    transform=ccrs.PlateCarree()
  )
  ax.clabel(contour_minus)

if uvar == "Ug" or uvar == "Ua" :
  step=12
  qx,qy,qk=1.1,-0.1,20
  vector=ax.quiver(
  lons.values[::step,::step],
  lats.values[::step,::step],
  u1.values[::step,::step],
  v1.values[::step,::step],
  scale=400,
  width=0.005,
  color="black",
  transform=ccrs.PlateCarree()
    )

  ax.quiverkey(vector,qx,qy,qk,f"{qk}m/s")

# グリッド線の設定
draw_grid(ax)

ax.coastlines()

ax.set_title(f"{ex1}-{svar}-{cvar}-{vvar}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/QGVeq{ex1}{svar}_{cvar}_{uvar}{vvar}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

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


# if svar != 0 :
#   shade=ax.contourf(
#     lons,lats,svar2,
#     levels=levs,
#     cmap="bwr",
#     extend="both",
#     transform=ccrs.PlateCarree()
#     )

#   cbar=plt.colorbar(
#     shade,orientation="horizontal",
#     shrink=0.6,
#     aspect=25
#     )
  
#   cbar.ax.tick_params(labelsize=10)

if cvar != 0:
  contour_plus=ax.contour(
    lons,lats,cvar2*10**factor,
    levels=levc_plus,
    colors="black",linewidths=1.0,
    transform=ccrs.PlateCarree()
  )
  ax.clabel(contour_plus)

  contour_minus=ax.contour(
    lons,lats,cvar2*10**factor,
    levels=levc_minus,
    colors="black",linewidths=1.0,
    transform=ccrs.PlateCarree()
  )
  ax.clabel(contour_minus)

if uvar == "Ug" or uvar == "Ua" :
  step=12
  qx,qy,qk=1.1,-0.1,20
  vector=ax.quiver(
  lons.values[::step,::step],
  lats.values[::step,::step],
  u2.values[::step,::step],
  v2.values[::step,::step],
  scale=400,
  width=0.005,
  color="black",
  transform=ccrs.PlateCarree()
   )

  ax.quiverkey(vector,qx,qy,qk,f"{qk}m/s")

draw_grid(ax)
ax.coastlines()
ax.set_title(f"{ex1}-{svar}-{cvar}-{vvar}{lev}-dx=dy={dx}")

plt.savefig(f"{fig_dir}/QGVeq{ex2}{svar}_{cvar}_{uvar}{vvar}_{lev}_dxdy{dx}_sigma{sigma}.png")

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

# if svar != 0 :
#   shade=ax.contourf(
#     lons,lats,sdif,
#     levels=levs_dif,
#     cmap="bwr",
#     extend="both",
#     transform=ccrs.PlateCarree()
#     )

#   cbar=plt.colorbar(
#     shade,orientation="horizontal",
#     shrink=0.6,
#     aspect=25
#     )
  
#   cbar.ax.tick_params(labelsize=10)


if cvar != 0:
  contour_plus=ax.contour(
  lons,lats,cdif*10**factor,
  levels=levc_plus,
  colors="red",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
  ax.clabel(contour_plus)

  contour_minus=ax.contour(
  lons,lats,cdif*10**factor,
  levels=levc_minus,
  colors="blue",linewidths=1.0,
  transform=ccrs.PlateCarree()
    )
  ax.clabel(contour_minus)

  
if uvar == "Ug" or uvar == "Ua" :
  step=12
  qx,qy,qk=1.1,-0.1,5
  vector=ax.quiver(
  lons.values[::step,::step],
  lats.values[::step,::step],
  udif.values[::step,::step],
  vdif.values[::step,::step],
  scale=50,
  width=0.005,
  color="black",
  transform=ccrs.PlateCarree()
   )

  ax.quiverkey(vector,qx,qy,qk,f"{qk}m/s")

draw_grid(ax)
ax.coastlines()
#ax.set_title(f"Dif-{svar}-{cvar}-{vvar}{lev}-dx=dy={dx}")


plt.savefig(f"{fig_dir}/QGVeqDif_{svar}_{cvar}_{uvar}{vvar}_{lev}hpa_dxdy{dx}_sigma{sigma}.png")

plt.close("all")

print("End Program")
