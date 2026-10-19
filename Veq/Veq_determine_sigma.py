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
e_sigma=20
delta_sigma=1
#σの初期値, 刻み幅

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

start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
dsec=dh*3600
print(f"Set Time : {start_date}----{end_date} dh={dh}")

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout"
fig_dir=f"{wrf_dir}/fig/{case}/Veq"

##時間方向に一定な定数
a=6370*10**3
g=9.81
omega=7.292*10**-5
#Holton,p12参照

dx=25*(10**3)
dy=25*(10**3)

print("a =",a)
print("g =",g)
print("Ω =",omega)
print("ref_lat =",ref_lat)
print("dx=dy=",dx)

#for sigma in sig_list:
list_absolute_error=[]
list_relative_error=[]
list_sigma=[]


print("Calculate G.F.")
print("sigma  absolute_error_index  relative_error_index")
min_relative_error = np.inf
best_sigma = None
for sigma in range(s_sigma, e_sigma+1, delta_sigma):
  list_theta=[]
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
    lat=getvar(ds,"lat")
    f=2*omega*np.sin(np.radians(lat))
    f = f.assign_coords({
      "west_east": np.arange(lat.west_east.size) * dx,
      "south_north": np.arange(lat.south_north.size) * dy
    }) 
  
    ua=getvar(ds,"ua",units="ms-1")
    va=getvar(ds,"va",units="ms-1")

    u=interplevel(ua,p,lev)
    v=interplevel(va,p,lev)


#計算前に各変数へガウシアンフィルタをかける
    u=gfilter(u,sigma)
    v=gfilter(v,sigma)

    ds.close()

##計算
    if map_proj == "lambert":
##変数へユークリッド距離に応じた座標を付与
      u = u.assign_coords({
        "west_east": np.arange(u.west_east.size) * dx,
        "south_north": np.arange(u.south_north.size) * dy
      })
      
      v = v.assign_coords({
         "west_east": np.arange(v.west_east.size) * dx,
         "south_north": np.arange(v.south_north.size) * dy
            })

#微分計算  
    dudx=u.differentiate("west_east")
    dudy=u.differentiate("south_north")
    dvdx=v.differentiate("west_east")
    dvdy=v.differentiate("south_north")

#渦度
    theta=dvdx - dudy
    if date == start_date:
          cart_proj=get_cartopy(lat)
          lats,lons=latlon_coords(lat)
          draw_latlon_contour(lons, lats, theta*10**5, cart_proj, 
                         -9.6, 9.7, 1.6, f"Sigma = {sigma} [1e-5 /s]")
          plt.savefig(f"{fig_dir}/Initial_Theta_sigma={sigma:02d}.png")
          plt.close("all")
    
#相対渦度移流
    adv_thx=-u*theta.differentiate("west_east")
    adv_thy=-v*theta.differentiate("south_north")
    
#惑星渦度移流
    adv_f=-v*f.differentiate("south_north")

#ストレッチング項
    strt=-(theta + f)*(dudx + dvdy)


    if map_proj == "mercator":
      rate_cos=np.cos(ref_lat*np.pi/180)/np.cos(lat*np.pi/180)
 ##メルカトルは未完成   


    right_hand=adv_thx + adv_thy + adv_f + strt

    list_theta.append(theta)
    list_right.append(right_hand)

    date+=timedelta(hours=dh)

#ここでオイラー微分項を計算
  array_theta=xr.concat(list_theta, dim="time")
  array_theta=array_theta.assign_coords({
      "time": np.arange(array_theta.Time.size)*dsec
    })
  left=array_theta.differentiate("time") 

  right=xr.concat(list_right, dim="time")
  

  absolute_error=abs(left - right)
  scale = abs(left) + abs(right)
#収支があっているほどこれらの値は0に近づく

  absolute_error_index =absolute_error.mean(dim="time", skipna=True).mean(
        dim=["south_north", "west_east"],
        skipna=True
    ) .item()

  scale_index =scale.mean(dim="time", skipna=True).mean(
        dim=["south_north", "west_east"],
        skipna=True
    ) .item()

  relative_error_index=absolute_error_index / scale_index
  
  print(f"{sigma}  {absolute_error_index}  {relative_error_index}")

  if relative_error_index < min_relative_error:
     min_relative_error = relative_error_index
     best_sigma = sigma      

     best_theta = array_theta.copy(deep=True)


  list_absolute_error.append(absolute_error_index)
  list_relative_error.append(relative_error_index)

  list_sigma.append(sigma)

print("--------------------------------")
print(f"Minimum Relative Error = {min_relative_error} Sigma={best_sigma}")
print("--------------------------------")



plt.rcParams["font.size"]=14
fig=plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(list_sigma, list_absolute_error, color="blue", lw=3)
ax2=ax1.twinx()
ax2.plot(list_sigma, list_relative_error, color="red", lw=5)
ax2.scatter(best_sigma, min_relative_error, color="black", s=80, zorder=5)

ax1.set_yscale("log")
ax1.set_xlabel("Sigma")
ax1.set_ylabel("Absolute Error", color="blue")
ax2.set_ylabel("Relative Error", color="red")
#plt.axhline(y=10**-10, color="red", lw=2, linestyle="--")
plt.tight_layout()
plt.savefig(f"{fig_dir}/Sigma_Error.png")


##補足のアニメーション
#各σにおける初期時刻の相対渦度分布
create_animation_from_png(fig_dir, "Initial_Theta_sigma=", 1)
#相対誤差が最も小さいシグマでの相対渦度時系列を書く
for t in range(best_theta.sizes["time"]):
  draw_latlon_contour(
    lons, lats, best_theta[t]*10**5, cart_proj,
    -9.6, 9.7, 1.6,
     f"t={t+1} (sigma={best_sigma}) [1e-5 /s]"
    )
  plt.savefig(f"{fig_dir}/Theta_sigma={best_sigma}_t={t+1:03d}.png")
  plt.close("all")
create_animation_from_png(fig_dir, f"Theta_sigma={best_sigma}", 1)
