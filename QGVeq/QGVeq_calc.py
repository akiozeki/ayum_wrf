###ここでは実験ごとに渦度方程式の右辺各項(Adv.theta_g, Adv.f, Stretching)を求めてnetCDFファイルに保存する
from scipy.ndimage import gaussian_filter
from datetime import datetime,timedelta
import math
from netCDF4 import Dataset
from wrf import getvar,interplevel
import numpy as np
import xarray as xr


def gfilter(var,gs):
  var_smooth=xr.DataArray(
    gaussian_filter(var,sigma=(gs,gs)),
    coords=var.coords,
    dims=var.dims
  )
  return(var_smooth)

##実験設定
case="2025winter"
#ex="CTL"
ex="ME00"

map_proj="lambert"
#"lambert" or "mercator" 
ref_lat=42.5
#メルカトルであれば使用する

lev=750
gsigma=6
#ガウシアンフィルタで用いる標準偏差(0でナシ)
#gtime="start"
gtime="end"
#ガウシアンフィルタで用いる標準偏差(0でナシ)

start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
dsec=dh*3600
print(f"Set Time : {start_date}----{end_date} dh={dh}")

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}/{ex}/n1"

##定数
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

listz=[]
listug=[]
listvg=[]
listua=[]
listva=[]
list=[]
listx=[]
listy=[]
listf=[]
lists=[]

date=start_date
while date <= end_date:
  y=date.year
  m=date.month
  d=date.day
  h=date.hour
  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}:00:00"
  print(time_str)

##変数取り出し
  ds=Dataset(f"{wrfout_dir}/wrfout_d01_{time_str}")
  
  lat=getvar(ds,"lat")
  p=getvar(ds,"p",units="hpa")

  ua=getvar(ds,"ua",units="ms-1")
  va=getvar(ds,"va",units="ms-1")
  za=getvar(ds,"z",units="m")

  u=interplevel(ua,p,lev)
  v=interplevel(va,p,lev)
  z=interplevel(za,p,lev)


#計算前に各変数へガウシアンフィルタをかける
  if gtime == "start" :
    z=gfilter(z,gsigma)
    u=gfilter(u,gsigma)
    v=gfilter(v,gsigma)

#  print(z)

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

#地衡風と渦度
    ug=-(g/f0)*dzdy
    vg=(g/f0)*dzdx
    zetag=(g/f0)*(dz2dx2 + dz2dy2)

#非地衡風成分
    ua=u - ug
    va=v - vg
    duadx=ua.differentiate("west_east")
    dvady=va.differentiate("south_north")

#相対渦度移流
    adv_zetagx= -(g**2/f0**2)*(-dzdy*(dz3dx3 + dz3dxdy2))
    adv_zetagy= - (g**2/f0**2)*(dzdx*(dz3dx2dy + dz3dy3))
    

#惑星渦度移流
    beta=2*omega*np.cos(np.radians(ref_lat))/a
    adv_f=-beta*vg

#ストレッチング項
    str=-f0*(duadx + dvady)


  if map_proj == "mercator":
    rate_cos=np.cos(ref_lat*np.pi/180)/np.cos(lat*np.pi/180)
 ##メルカトルは未完成   

  if gtime == "end" :
    zetag=gfilter(zetag, gsigma)
    adv_zetagx=gfilter(adv_zetagx, gsigma)
    adv_zetagy=gfilter(adv_zetagy, gsigma)
    adv_f=gfilter(adv_f, gsigma)
    str=gfilter(str, gsigma)


  listz.append(z)  
  listug.append(ug)
  listvg.append(vg)
  listua.append(ua)
  listva.append(va)
  list.append(zetag)
  listx.append(adv_zetagx)
  listy.append(adv_zetagy)
  listf.append(adv_f)
  lists.append(str)

  date+=timedelta(hours=dh)

daz=xr.concat(listz,dim="time")
daug=xr.concat(listug,dim="time")
davg=xr.concat(listvg,dim="time")
daua=xr.concat(listua,dim="time")
dava=xr.concat(listva,dim="time")
da=xr.concat(list,dim="time")
dax=xr.concat(listx,dim="time")
day=xr.concat(listy,dim="time")
daf=xr.concat(listf,dim="time")
das=xr.concat(lists,dim="time")

daz=daz.assign_coords({
      "time": np.arange(daz.Time.size)*dh*3600 
        })
daug=daug.assign_coords({
      "time": np.arange(daug.Time.size)*dh*3600 
        })
davg=davg.assign_coords({
      "time": np.arange(davg.Time.size)*dh*3600 
        })
daua=daua.assign_coords({
      "time": np.arange(daua.Time.size)*dh*3600 
        })
dava=dava.assign_coords({
      "time": np.arange(dava.Time.size)*dh*3600 
        })
da=da.assign_coords({
      "time": np.arange(da.Time.size)*dh*3600 
        })
dax=dax.assign_coords({
      "time": np.arange(dax.Time.size)*dh*3600 
        })
day=day.assign_coords({
      "time": np.arange(day.Time.size)*dh*3600 
        })
daf=daf.assign_coords({
      "time": np.arange(daf.Time.size)*dh*3600 
        })
das=das.assign_coords({
      "time": np.arange(das.Time.size)*dh*3600 
        })

#ここでオイラー微分項を計算
da_dt=da.differentiate("time")

daz=daz.rename("Z")
daug=daug.rename("Ug")
davg=davg.rename("Vg")
daua=daua.rename("Ua")
dava=dava.rename("Va")
da=da.rename("zetag")
da_dt=da_dt.rename("dzetag_dt")
dax=dax.rename("Adv.zetagx")
day=day.rename("Adv.zetagy")
daf=daf.rename("Adv.f")
das=das.rename("Stretching")

ds_out=xr.merge([daz,daug,davg,daua,dava,da,da_dt,dax,day,daf,das])
print(ds_out)

#ランベルト図法はこれを表す変数がnetCDF形式に対応しないので削除
ds_out.attrs.pop("projection", None)
for var in ds_out.data_vars:
    ds_out[var].attrs.pop("projection", None)


output=f"{wrfout_dir}/{ex}{lev}hPa_QGVeq_dxdy{int(dy/1000)}_{gtime}sigma{gsigma}.nc"

ds_out.to_netcdf(output)
print("Save :",output)