from netCDF4 import Dataset
import xarray as xr
from wrf import getvar,interplevel
from datetime import datetime

"""
ある時刻の特定気圧面の場についてアンサンブル平均を計算する

引数
dir -> n1, n2....というディレクトリが並んでいるディレクトリ
ensN -> アンサンブルサイズ
time -> 時刻 datetime(y, m, d ,h)
var, unit -> 変数名とその単位 getvarで選択できるもの
level -> 気圧面(hpa)
"""
def ensmean_horizontal_P(dir, ensN, time, var, unit, level):
  y=time.year
  m=time.month
  d=time.day
  h=time.hour
  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}"
  print(f"Calculate Ensemble Mean(N={ensN}) var={var} time={time_str} level={level}")
  
  var_list=[]
  for N in range(1, ensN+1):
    wrfout=f"{dir}/n{N}/wrfout_d01_{time_str}:00:00"
    ds=Dataset(wrfout)
    p=getvar(ds, "p", units="hpa")
    var=getvar(ds, var, units=unit)
    pvar=interplevel(var, p, level)

    var_list.append(pvar)

  darray=xr.concat(var_list, dim="member")
  ens_mean=darray.mean(dim="member", skipna=True)

  return ens_mean