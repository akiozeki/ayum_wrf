##ココでは温位と西風をプロット
from __future__ import print_function, division
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from wrf import getvar, vertcross, CoordPair,get_cartopy,latlon_coords
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from myFunction.vertical_cross_section import latlon_to_yx, get_vertical_cross_section, get_vertical_cross_section_nounit
from myFunction.draw_field import draw_vertical_cross, draw_vertical_cross_onlysig
#from myFunction.ttest_Nstar import ttest_vertical_cross
from myFunction.ttest import ttest_vertical_cross

n=1
case="2025winter"
ex1="CTL"
ex2="ME00"

vcord="p"
vunit="hPa"
left=31
right=60
vtop=500
vbottom=1000
xstep=1
vstep=8
#ベクトルの間引き間隔

levs=np.arange(-0.14,0,0.02)
levs_dif=np.arange(-7.5*10**-2,8.5*10**-2,1*10**-2)
levc_dif=np.arange(-4.5,5.5,1)
levc=np.arange(260,324,4)
# levc_dif=np.arange(-5,6,1.0)


wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
fig_dir=f"{wrf_dir}/fig/{case}"


start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
print(f"Set Time : {start_date}----{end_date}")

theta_list1=[]
theta_list2=[]
dtheta_dp_list1=[]
dtheta_dp_list2=[]
omg_list1=[]
omg_list2=[]

#断面図の位置設定
##例のごとくlatlon->yx変換
start_lat,start_lon=50,120
end_lat,end_lon=30,135
if n == 0:
  wrfinput=f"{wrfout_dir}/{ex1}/wrfinput_d01"
elif n != 0:
  wrfinput=f"{wrfout_dir}/{ex1}/n{n}/wrfinput_d01"

start=latlon_to_yx(wrfinput, start_lat, start_lon)
end=latlon_to_yx(wrfinput, end_lat, end_lon)

date=start_date
while date <= end_date:
  y=date.year
  m=date.month
  d=date.day
  h=date.hour

  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}"
  if n == 0:
    wrfout1=f"{wrfout_dir}/{ex1}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir}/{ex2}/wrfout_d01_{time_str}:00:00"

  if n != 0:
    wrfout1=f"{wrfout_dir}/{ex1}/n{n}/wrfout_d01_{time_str}:00:00"
    wrfout2=f"{wrfout_dir}/{ex2}/n{n}/wrfout_d01_{time_str}:00:00"
#ここで「n〇」とすれば任意のアンサンブルメンバーを1つ取り出せる

  print(f"Read wrfout File {time_str}")

  ds1=Dataset(wrfout1)
  ds2=Dataset(wrfout2)

#断面データを取り出す
  theta_vert1=get_vertical_cross_section(ds1, "theta", "K", start, end, vcord)
  theta_vert2=get_vertical_cross_section(ds2, "theta", "K", start, end, vcord)
  omg_vert1=get_vertical_cross_section_nounit(ds1, "omega", start, end, vcord)
  omg_vert2=get_vertical_cross_section_nounit(ds2, "omega", start, end, vcord)

  vert=theta_vert1.vertical.values
  dtheta_dp1=xr.DataArray(
      np.gradient(theta_vert1.values,vert,axis=0),
      coords=theta_vert1.coords,
      dims=theta_vert1.dims
  )
  dtheta_dp2=xr.DataArray(
      np.gradient(theta_vert2.values,vert,axis=0),
      coords=theta_vert2.coords,
      dims=theta_vert2.dims
  )
 # print(var_vert1)
  #配列は(鉛直, 水平)

#座標情報を取りだす
  if date == end_date:
#ここでどの断面を取ったかを表す図も描いておく
    ter=getvar(ds1,"ter",units="m")
    cart_proj=get_cartopy(ter)
    lats,lons=latlon_coords(ter)
    fig=plt.figure()
    ax=plt.axes(projection=cart_proj)
    ax.plot(
      [start_lon,end_lon],
      [start_lat,end_lat],
      color="blue",
      linewidth=3,
      transform=ccrs.PlateCarree()
    )

    ax.contourf(
      lons,lats,ter,
      cmap="Greys",
      extend="max",
      levels=np.arange(0,2000,100),
      transform=ccrs.PlateCarree()
    )
    # グリッド線の設定
    gl = ax.gridlines(draw_labels=True, 
                  linewidth=1, 
                  color='gray', 
                  alpha=0.5, 
                  linestyle='--')
# ラベルの表示位置を制御
    gl.xlines = True         # 経度線を描く
    gl.ylines = True         # 緯度線を描く
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator=FixedLocator([127,130,133,136,139])
    #gl.ylocator=FixedLocator()
    gl.x_inline = False      # ラベルを図の内側（インライン）に書かない設定
    gl.y_inline = False      # 緯度も念のため設定
    gl.ylabel_style = {'rotation': 0}
    gl.xlabel_style = {'rotation': 0}
    gl.xpadding = 10
    ax.coastlines()

    plt.savefig(f"{fig_dir}/vert_cross_line{start}to({end}.png")
    plt.close("all")

    
    idx=theta_vert1.cross_line_idx.values



  theta_list1.append(theta_vert1)
  theta_list2.append(theta_vert2)

  dtheta_dp_list1.append(dtheta_dp1)
  dtheta_dp_list2.append(dtheta_dp2)

  omg_list1.append(omg_vert1)
  omg_list2.append(omg_vert2)

  ds1.close()
  ds2.close()

  date+=timedelta(hours=dh)

theta_array1=np.stack(theta_list1,axis=0)
theta_mean1=np.nanmean(theta_array1,axis=0)
#print(theta1_array)
#print(theta1_mean)

theta_array2=np.stack(theta_list2,axis=0)
theta_mean2=np.nanmean(theta_array2,axis=0)

omg_array1=np.stack(omg_list1,axis=0)
omg_mean1=np.nanmean(omg_array1,axis=0)
#print(omg1_mean)

omg_array2=np.stack(omg_list2,axis=0)
omg_mean2=np.nanmean(omg_array2,axis=0)


dtheta_dp_array1=np.stack(dtheta_dp_list1,axis=0)
dtheta_dp_mean1=np.nanmean(dtheta_dp_array1,axis=0)

dtheta_dp_array2=np.stack(dtheta_dp_list2,axis=0)
dtheta_dp_mean2=np.nanmean(dtheta_dp_array2,axis=0)

theta_dif=theta_mean1 - theta_mean2
omg_dif=omg_mean1 - omg_mean2
dtheta_dp_dif=dtheta_dp_mean1 - dtheta_dp_mean2

print(theta_dif)

xdummy=np.zeros_like(omg_mean1)
nan1=np.isnan(theta_mean1)
nan2=np.isnan(theta_mean2)
#x方向は描かないので，0埋め配列を置く


#Plot 
#ex1
plt.rcParams["font.size"]=16
fig=plt.figure()
ax=plt.axes()
ax.set_xticks([])      
ax.set_xlabel("Horizontal")
ax.set_ylabel(f"Vertical [{vunit}]")
ax.set_xlim(left,right)
ax.set_ylim(vbottom,vtop)
ax.contourf(idx,vert,nan1,levels=[0.5,1.5],colors=["black"])

shade=ax.contourf(idx,vert,dtheta_dp_mean1,levels=levs,extend="min",cmap="Reds_r")
cbar=plt.colorbar(shade)
#cbar.ax.set_xlabel("K",rotation=90,labelpad=20)
#cbar.ax.tick_params(labelsize=10)

contour=ax.contour(idx,vert,theta_mean1,levels=levc,colors="black",linewidths=1.0)
ax.clabel(contour)

qx,qy,qk=0.9,-0.1,1
vector=ax.quiver(
 idx[::xstep],
 vert[::vstep],
 xdummy[::vstep,::xstep],
 -omg_mean1[::vstep,::xstep],
 scale=10,
 width=0.008,
 color="orange"
 )

ax.quiverkey(vector,qx,qy,qk,"1 Pa/s")



#ax.set_title(f"{ex1}-Vertical")
plt.tight_layout()
#if vcord == "pressure" : ax.invert_yaxis()#鉛直軸反転
plt.savefig(f"{fig_dir}/{ex1}-{case}_atheta_ap_vertical.png")
plt.close("all")

#ex2
plt.rcParams["font.size"]=16
fig=plt.figure()
ax=plt.axes()
ax.set_xticks([])      
ax.set_xlabel("Horizontal")
ax.set_ylabel(f"Vertical [{vunit}]")
ax.set_xlim(left,right)
ax.set_ylim(vbottom,vtop)
ax.contourf(idx,vert,nan2,levels=[0.5,1.5],colors=["black"])

shade=ax.contourf(idx,vert,dtheta_dp_mean2,levels=levs,extend="min",cmap="Reds_r")
plt.colorbar(shade)
    
contour=ax.contour(idx,vert,theta_mean2,levels=levc,colors="black",linewidths=1.0)
ax.clabel(contour)


qx,qy,qk=0.9,-0.1,1
vector=ax.quiver(
  idx[::xstep],
  vert[::vstep],
  xdummy[::vstep,::xstep],
  -omg_mean2[::vstep,::xstep],
  scale=10,
  width=0.008,
  color="orange"
  )

ax.quiverkey(vector,qx,qy,qk,"1 Pa/s")

#ax.set_title(f"{ex2}-Vertical")
plt.tight_layout()
#if vcord == "pressure" : ax.invert_yaxis()#鉛直軸反転不要
plt.savefig(f"{fig_dir}/{ex2}-{case}vertical_atheta_ap.png")
plt.close("all")

#difference
plt.rcParams["font.size"]=16
fig=plt.figure()
ax=plt.axes()
ax.set_xticks([])      
ax.set_xlabel("Horizontal")
ax.set_ylabel(f"Vertical [{vunit}]")
ax.set_xlim(left,right)
ax.set_ylim(vbottom,vtop)
ax.contourf(idx,vert,nan1,levels=[0.5,1.5],colors=["black"])

shade=ax.contourf(idx,vert,dtheta_dp_dif,levels=levs_dif,cmap="bwr",extend="both")
plt.colorbar(shade)
    
contour=ax.contour(idx,vert,theta_dif,levels=levc_dif,colors="black",linewidths=1.0)


ax.clabel(contour)

qx,qy,qk=0.9,-0.1,1
vector=ax.quiver(
  idx[::xstep],
  vert[::vstep],
  xdummy[::vstep,::xstep],
  -omg_dif[::vstep,::xstep],
  scale=5,
  width=0.008,
  color="green"
  )

ax.quiverkey(vector,qx,qy,qk,"1 Pa/s")

#ax.set_title(f"{ex1} - {ex2}-Vertical")
plt.tight_layout()
#if vcord == "pressure" : ax.invert_yaxis()#鉛直軸反転->不要
plt.savefig(f"{fig_dir}//Dif-{case}vertical_atheta_ap.png")
plt.close("all")

print("End Program")    
