### アンサンブルメンバーごとに渦度方程式の各項を計算し、member 軸を追加して保存する
from pathlib import Path
from scipy.ndimage import gaussian_filter
from datetime import datetime, timedelta
from netCDF4 import Dataset
from wrf import getvar, interplevel
import numpy as np
import xarray as xr

"""水平ガウシアンフィルタ。時間方向にはフィルタをかけない。"""
def gfilter(var, gs):
    if var.ndim == 2:
        sigma = (gs, gs)
    elif var.ndim == 3:
        sigma = (0, gs, gs)
    else:
        raise ValueError(f"gfilter は 2次元または3次元を想定しています: ndim={var.ndim}")

    return xr.DataArray(
        gaussian_filter(var, sigma=sigma),
        coords=var.coords,
        dims=var.dims,
        attrs=var.attrs,
        name=var.name,
    )


"""各種設定"""
case="2025winter"
#ex="CTL"
ex="ME00"

map_proj = "lambert"  
#"lambert" or "mercator"（mercator には未対応）
ref_lat=42.5

lev=750
top=700
bottom=800

gsigma=6
# ガウシアンフィルタで用いる標準偏差（0 でフィルタなし）
# gtime = "start"
gtime="end"  
#"start" または "end"->基本endが良い

start_date=datetime(2024, 12, 1, 0)
end_date=datetime(2025, 2, 28, 18)
dh=3
dsec=dh * 3600
print(f"Set Time: {start_date} ---- {end_date}, dh={dh} h")

wrf_dir="/DATA/USER/ozeki/MyWRF"
# 各メンバーのディレクトリ: .../wrfout/{case}/{ex}/n1, n2, ...
wrfout_base=f"{wrf_dir}/wrfout/{case}/{ex}"

# 使用するアンサンブルメンバー
total_member=5
members = [f"n{i}" for i in range(1, total_member+1)]

## 定数
a = 6370 * 10**3
g = 9.81
OMEGA = 7.292 * 10**-5
ref_lat = 42.5

dx = 25 * 10**3
dy = 25 * 10**3

print("a =", a)
print("g =", g)
print("Ω =", OMEGA)
print("ref_lat =", ref_lat)
print("dx=dy =", dx)


"""指定されたアンサンブルメンバー n1, n2, ... のデータを計算する関数を定義する"""
def calc_member(member):

  wrfout_dir = f"{wrfout_base}/{member}"

  print(f"\n===== Processing member: {member} =====")
  listu=[]
  listusm=[]
  listv=[]
  listvsm=[]
  listomg=[]
  listomgsm=[]
  list=[]
  listsm=[]
  listx=[]
  listy=[]
  listp=[]
  listf=[]
  lists=[]
  listt=[]

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
    f=2*OMEGA*np.sin(np.radians(lat))
    f = f.assign_coords({
      "west_east": np.arange(lat.west_east.size) * dx,
      "south_north": np.arange(lat.south_north.size) * dy
    })      

    p=getvar(ds,"p",units="hpa")

    ua=getvar(ds,"ua",units="ms-1")
    va=getvar(ds,"va",units="ms-1")
    omga=getvar(ds, "omega") / 100.0
    #Pa/s -> hPa/s 

    u=interplevel(ua, p, lev)
    utop=interplevel(ua, p, top)
    ubottom=interplevel(ua, p, bottom)

    v=interplevel(va, p, lev)
    vtop=interplevel(va, p, top)
    vbottom=interplevel(va, p, bottom)

    omg=interplevel(omga, p, lev)


    #先にガウシアンフィルタする場合はここで
    if gtime == "start" :
      u=gfilter(u, gsigma)
      utop=gfilter(utop, gsigma)
      ubottom=gfilter(ubottom, gsigma)

      v=gfilter(v, gsigma)
      vtop=gfilter(vtop, gsigma)
      vbottom=gfilter(vbottom, gsigma)

      omg=gfilter(omg, gsigma)

    ds.close()

    ##計算
    if map_proj == "lambert":
    ##変数へユークリッド距離 (m) に応じた座標を付与
      coord = {
        "west_east": np.arange(u.west_east.size)*dx,
        "south_north": np.arange(u.south_north.size)*dy
        }

      u=u.assign_coords(coord)
      utop=utop.assign_coords(coord)
      ubottom=ubottom.assign_coords(coord)

      v=v.assign_coords(coord)
      vtop=vtop.assign_coords(coord)
      vbottom=vbottom.assign_coords(coord)

      omg=omg.assign_coords(coord)

#水平微分計算    
      dudx=u.differentiate("west_east")
      dudy=u.differentiate("south_north")
      dudytop=utop.differentiate("south_north")
      dudybottom=ubottom.differentiate("south_north")

      dvdx=v.differentiate("west_east")
      dvdxtop=vtop.differentiate("west_east")
      dvdxbottom=vbottom.differentiate("west_east")
      dvdy=v.differentiate("south_north")
      domgdx=omg.differentiate("west_east")
      domgdy=omg.differentiate("south_north")

#鉛直微分計算
      dudp=(utop - ubottom) / (top - bottom)
      dvdp=(vtop - vbottom) / (top - bottom)

#渦度
      zeta=dvdx - dudy
      zetatop=dvdxtop - dudytop
      zetabottom=dvdxbottom - dudybottom
        
#水平相対渦度移流
      adv_zetax=-u*zeta.differentiate("west_east")
      adv_zetay=-v*zeta.differentiate("south_north")

#惑星渦度移流
      adv_f=-v*f.differentiate("south_north")

#鉛直相対渦度移流
      dzetaadp=(zetatop - zetabottom) / (top - bottom)
      adv_zetap=-omg*dzetaadp

#ストレッチング項
      str=-(zeta + f)*(dudx + dvdy)

#ティルティング
      tilt=dudp*domgdy - dvdp*domgdx

    # if map_proj == "mercator":
    #   rate_cos=np.cos(ref_lat*np.pi/180)/np.cos(lat*np.pi/180)
    #  ##メルカトルは未完成

    if gtime == "end" :
      u_sm=gfilter(u, gsigma)
      v_sm=gfilter(v, gsigma)
      omg_sm=gfilter(omg,gsigma)
      zeta_sm=gfilter(zeta, gsigma)
      adv_zetax=gfilter(adv_zetax, gsigma)
      adv_zetay=gfilter(adv_zetay, gsigma)
      adv_zetap=gfilter(adv_zetap, gsigma)
      adv_f=gfilter(adv_f, gsigma)
      str=gfilter(str, gsigma)
      tilt=gfilter(tilt, gsigma)

    listu.append(u)
    listv.append(v)
    listomg.append(omg)
    list.append(zeta)
    listx.append(adv_zetax)
    listy.append(adv_zetay)
    listp.append(adv_zetap)
    listf.append(adv_f)
    lists.append(str)
    listt.append(tilt)

    if gtime == "end":
      listusm.append(u_sm)
      listvsm.append(v_sm)
      listomgsm.append(omg_sm)
      listsm.append(zeta_sm)

    date+=timedelta(hours=dh)


  dau = xr.concat(listu, dim="time")  
  dav = xr.concat(listv, dim="time")
  daomg=xr.concat(listomg, dim="time")
  da = xr.concat(list, dim="time")
  dax=xr.concat(listx, dim="time")
  day=xr.concat(listy, dim="time")
  dap=xr.concat(listp, dim="time")
  daf=xr.concat(listf, dim="time")
  das=xr.concat(lists, dim="time")
  dat=xr.concat(listt, dim="time")

  dau=dau.assign_coords({
    "time": np.arange(dau.Time.size)*dh*3600 
    })
  dav=dav.assign_coords({
    "time": np.arange(dav.Time.size)*dh*3600 
    })
  daomg=daomg.assign_coords({
    "time": np.arange(daomg.Time.size)*dh*3600 
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
  dap=dap.assign_coords({
    "time": np.arange(dap.Time.size)*dh*3600 
    })
  daf=daf.assign_coords({
    "time": np.arange(daf.Time.size)*dh*3600 
    })
  das=das.assign_coords({
    "time": np.arange(das.Time.size)*dh*3600 
    })
  dat=dat.assign_coords({
    "time": np.arange(dat.Time.size)*dh*3600 
    })

  if gtime == "end":
    dau_sm = xr.concat(listusm, dim="time")  
    dav_sm = xr.concat(listvsm, dim="time")
    daomg_sm=xr.concat(listomgsm, dim="time")
    da_sm = xr.concat(listsm, dim="time")
    dau_sm=dau_sm.assign_coords({
      "time": np.arange(dau_sm.Time.size)*dh*3600 
      })
    dav_sm=dav_sm.assign_coords({
      "time": np.arange(dav_sm.Time.size)*dh*3600 
      })
    daomg_sm=daomg_sm.assign_coords({
      "time": np.arange(daomg_sm.Time.size)*dh*3600 
      })
    da_sm=da_sm.assign_coords({
      "time": np.arange(da_sm.Time.size)*dh*3600 
      })

#ここでオイラー微分項を計算
  da_dt=da.differentiate("time")
  if gtime == "end":
    da_dt=gfilter(da_dt, gsigma)
    #print(da_dt)

  dau=dau.rename("u")
  dav=dav.rename("v")
  daomg=daomg.rename("omg")
  da=da.rename("zeta")
  da_dt=da_dt.rename("dzeta_dt")
  dax=dax.rename("Adv.zeta_x")
  day=day.rename("Adv.zeta_y")
  dap=dap.rename("Adv.zeta_p")
  daf=daf.rename("Adv.f")
  das=das.rename("Stretching")
  dat=dat.rename("Tilting")

  if gtime == "start":
    ds_out=xr.merge([dau,dav,daomg,da,da_dt,dax,day,dap,daf,das,dat])

  elif gtime == "end":
    dau_sm=dau_sm.rename("u_smooth")
    dav_sm=dav_sm.rename("v_smooth")
    daomg_sm=daomg_sm.rename("omg_smooth")
    da_sm=da_sm.rename("zeta_smooth")
    ds_out=xr.merge([dau,dau_sm,dav,dav_sm,daomg,daomg_sm,da,da_sm,da_dt,dax,day,dap,daf,das,dat])
     
  #print(ds_out)
  return ds_out


"""各メンバーを計算し、最後に member 軸で結合する"""
member_datasets = []
for member in members:
    member_ds = calc_member(member)
    member_datasets.append(member_ds)

# 出力変数の次元は (member, time, south_north, west_east) となる
# member 座標には n1, n2, ... を格納する
ds_out = xr.concat(
    member_datasets,
    dim=xr.IndexVariable("member", range(len(members))),
    join="exact",
    compat="equals",
)

print("\n===== Ensemble dataset =====")
print(ds_out)

# save = 0
save = 1

if save == 1:
# Lambert 図法の projection 属性は netCDF 保存時の問題を避けるため削除
  ds_out.attrs.pop("projection", None)
  for var in ds_out.data_vars:
    ds_out[var].attrs.pop("projection", None)

  output = (f"{wrfout_base}/{ex}Voleq{lev}hPa_dxdy{int(dy / 1000)}_{gtime}sigma{gsigma}_ensemble.nc")
  ds_out.to_netcdf(output)
  print("Save:", output)

print("End Program")
