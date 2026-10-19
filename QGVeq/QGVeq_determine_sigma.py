###ここでは実験ごとに渦度方程式の右辺各項(Adv.theta_g, Adv.f, Stretching)を求めてnetCDFファイルに保存する
from scipy.ndimage import gaussian_filter
from datetime import datetime,timedelta
import math
from netCDF4 import Dataset
from wrf import getvar,interplevel,get_cartopy,latlon_coords
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from PIL import Image
import glob
import pprint

from myFunction.draw_field import draw_latlon_contour
from myFunction.animation import create_animation_from_png

s_sigma=0
e_sigma=30
delta_sigma=1
#σの初期値, 終了値, 刻み幅
sig_list=np.arange(s_sigma, e_sigma+0.01, delta_sigma)


def gfilter(var,gs):
  var_smooth=xr.DataArray(
    gaussian_filter(var,sigma=(gs,gs)),
    coords=var.coords,
    dims=var.dims
  )
  return(var_smooth)

##実験設定
case="2025winter"
ex="CTL"

map_proj="lambert"
#"lambert" or "mercator" 
ref_lat=42.5
#メルカトルであれば使用する

lev=800
delta_p=10
top=lev - delta_p
bottom=lev + delta_p

start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=6
dsec=dh*3600
print(f"Set Time : {start_date}----{end_date} dh={dh}")

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout"
fig_dir=f"{wrf_dir}/fig/{case}/QGVeq"

##時間方向に一定な定数
a=6370*10**3
g=9.81
omega=7.292*10**-5
#Holton,p12参照
ref_lat=42.5
f0=2*omega*math.sin(math.radians(ref_lat))

dx=25*(10**3)
dy=25*(10**3)

print("a =",a)
print("g =",g)
print("Ω =",omega)
print("ref_lat =",ref_lat)
print("f0 =",f0)
print("dx=dy=",dx)

#for sigma in sig_list:
code=0
sigma=s_sigma
list_error=[]
list_sigma=[]

while code <= 5:
  list_z=[]
  list_thg=[]
  list_left=[]
  list_right=[]  

  date=start_date
  while date <= end_date:
    y=date.year
    m=date.month
    d=date.day
    h=date.hour
    time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}:00:00"
#    print(time_str)

##変数取り出し
    ds=Dataset(f"{wrfout_dir}/{case}/{ex}/n1/wrfout_d01_{time_str}")
  
    lat=getvar(ds,"lat")
    p=getvar(ds,"p",units="hpa")

    ua=getvar(ds,"ua",units="ms-1")
    va=getvar(ds,"va",units="ms-1")
    za=getvar(ds,"z",units="m")
    wa=getvar(ds, "omega")

    u=interplevel(ua,p,lev)
    v=interplevel(va,p,lev)
    z=interplevel(za,p,lev)
    wtop=interplevel(wa,p,top)
    wbottom=interplevel(wa,p,bottom)

    # u=u.isel(
    #   south_north=slice(5, -5),
    #   west_east=slice(5, -5)
    # )

    # v=v.isel(
    #       south_north=slice(5, -5),
    #       west_east=slice(5, -5)
    # )

    # z=z.isel(
    #       south_north=slice(5, -5),
    #       west_east=slice(5, -5)
    # )


#計算前に各変数へガウシアンフィルタをかける
    z=gfilter(z,sigma)
    if date == start_date:
      cart_proj=get_cartopy(lat)
      lats,lons=latlon_coords(lat)
      draw_latlon_contour(lons, lats, z, cart_proj, 
                     1600, 2000, 50, f"Sigma = {sigma}")
      plt.savefig(f"{fig_dir}/InitialZ_sigma={sigma:02d}.png")
      plt.close("all")

    u=gfilter(u,sigma)
    v=gfilter(v,sigma)
    wtop=gfilter(wtop,sigma)
    wbottom=gfilter(wbottom,sigma)

    ds.close()

##計算
    if map_proj == "lambert":
##Zの変数へユークリッド距離に応じた座標を付与
      z = z.assign_coords({
        "west_east": np.arange(z.west_east.size) * dx,
        "south_north": np.arange(z.south_north.size) * dy
      })

      u = u.assign_coords({
              "west_east": np.arange(u.west_east.size) * dx,
              "south_north": np.arange(u.south_north.size) * dy
            })
      
      v = v.assign_coords({
              "west_east": np.arange(v.west_east.size) * dx,
              "south_north": np.arange(v.south_north.size) * dy
            })

#微分計算    
      dzdx=z.differentiate("west_east")
      dzdy=z.differentiate("south_north")
      dz2dx2=dzdx.differentiate("west_east")
      dz2dy2=dzdy.differentiate("south_north")
      dz3dx3=dz2dx2.differentiate("west_east")
      dz3dy3=dz2dy2.differentiate("south_north")
      dz3dxdy2=dz2dy2.differentiate("west_east")
      dz3dx2dy=dz2dx2.differentiate("south_north")

      dwdp=(wtop - wbottom)/(2*delta_p*100)

#地衡風と渦度
      f0=omega*np.cos(ref_lat*np.pi/180)
      ug=-(g/f0)*dzdy
      vg=(g/f0)*dzdx
      thg=(g/f0)*(dz2dx2 + dz2dy2)

#非地衡風成分
      ua=u - ug
      va=v - vg
      duadx=ua.differentiate("west_east")
      dvady=va.differentiate("south_north")
#      duadx=u.differentiate("west_east")
#      dvady=v.differentiate("south_north")

#相対渦度移流
      adv_thgx= -(g**2/f0**2)*(-dzdy*(dz3dx3 + dz3dxdy2))
      adv_thgy= - (g**2/f0**2)*(dzdx*(dz3dx2dy + dz3dy3))
    

#惑星渦度移流
      beta=2*omega*np.cos(ref_lat*np.pi/180)/a
      adv_f=-beta*vg

#ストレッチング項
      strt=-f0*(duadx + dvady)
#      strt=f0*dwdp


    if map_proj == "mercator":
      rate_cos=np.cos(ref_lat*np.pi/180)/np.cos(lat*np.pi/180)
 ##メルカトルは未完成   
    right_hand=adv_thgx + adv_thgy + adv_f + strt

    list_z.append(z)
    list_thg.append(thg)
    list_right.append(right_hand)

    date+=timedelta(hours=dh)

  da=xr.concat(list_z, dim="time")
  thg=xr.concat(list_thg, dim="time")
#ここでオイラー微分項を計算
  left=thg.differentiate("time") / dsec

  right=xr.concat(list_right, dim="time")
  

  error=abs(left - right)


  # 右辺を時間平均->この値が各格子点で0になってほしい
  error_index = error.mean(
    dim="time",
    skipna=True
  )

# 空間平均
  error_index = error_index.mean(
    dim=["south_north", "west_east"],
    skipna=True
  ).item()

#  error_index=error.mean(dim=["time", "south_north", "west_east"], skipna=True).item()
  print(f"Calc G.F. (sigma={sigma}) {error_index} ")
  if error_index < 10**-10 :
    code+=1
    if code == 1:
      target_sigma=sigma
      for t in range(1, len(da["time"])+1):
        draw_latlon_contour(lons, lats, da[t-1], cart_proj, 
                             1600, 2000, 50, f"t={t}")
        plt.savefig(f"{fig_dir}/Z_sigma={target_sigma}_t={t:03d}.png")
        plt.close("all")       

  list_error.append(error_index)

  list_sigma.append(sigma)


  sigma+=delta_sigma

print(list_error)
plt.rcParams["font.size"]=14
plt.plot(list_sigma, list_error, color="green", lw=5)
plt.yscale("log")
plt.xlabel("Sigma")
plt.ylabel("Error")
plt.axhline(y=10**-10, color="red", lw=2, linestyle="--")
plt.tight_layout()
plt.savefig(f"{fig_dir}/Sigma_Error.png")


##補足のアニメーション
#各σにおける初期時刻のz分布
create_animation_from_png(fig_dir, "InitialZ_sigma=", 1)
create_animation_from_png(fig_dir, f"Z_sigma={target_sigma}", 1)
