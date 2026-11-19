###変数の移動平均場と偏差場を計算する　→修正版："_imp.py
from datetime import datetime,timedelta
from netCDF4 import Dataset
from wrf import getvar, interplevel, ALL_TIMES
import numpy as np
import xarray as xr

######################################
#変数
var_name="va"
#wrf-pythonで取り出すときの変数名
var_unit="ms-1"
#var_unit=0
#wrf-pythonで取り出す変数の単位(単位指定がない場合は0を指定)

var_type="pres"
p=850
#計算する気圧面

#var_type="surf"
#気圧面の値か表面の値か


#ファイルパスと時間
case="2025winter"
ex="CTL"
wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
output_dir=wrfout_dir
start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
#dh=12
print(f"Set Time : {start_date}----{end_date} dh={dh}")


##移動平均の設定##
moving_day=7
step_per_day=24 // dh
step_all=moving_day*step_per_day + 1
half_step=step_all // 2

#アンサンブルメンバー数 / アンサンブル軸がある場合はそのサイズ
ensmean="yes"
#"yes" or "no" noの場合以下の要素数は必ず1個
target_members=[1, 2, 3, 4, 5]
#取り出すメンバ
if ensmean == "no" and len(target_members) != 1:
  raise ValueError("If ensmean=no, target_member's size =1")

if ensmean =="yes":
  ens_size=len(target_members)
  saving_name=f"{output_dir}/{ex}/{ex}_{var_name}{p}_{moving_day}dMean.nc"

elif ensmean == "no":
  ens_size=1
  saving_name=f"{output_dir}/{ex}/n{target_members[0]}/{ex}_{var_name}{p}_{moving_day}dMean.nc"

#解像度
dx=25*1000
dy=25*1000
#######################

slist_bar=[]
slist_prime=[]
#保存前の配列を作るリスト

for member in target_members:
  print(f"n={member}")

  mlist_bar=[]
  mlist_prime=[]
#あるメンバについての(nt, ny, nx)の配列を作るためのリスト
  date=start_date
  while date <= end_date:

##指定した時間の範囲の中で必要な時刻を取り出す
    net_start=date - timedelta(hours = dh * half_step)
    net_end=date + timedelta(hours = dh *half_step) 

    if net_start < start_date or net_end > end_date:
      print(f"Skip:{date}")
      date+=timedelta(hours=dh)
      continue

    y=date.year
    m=date.month
    d=date.day
    h=date.hour
    time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}:00:00"
    print(f"time:{time_str} (n={member})")
     

##移動平均を取る範囲のデータを配列としてまとめる
    move_date=date - timedelta(hours=dh * half_step)
    print(f"{date}_Calculate_{moving_day}days_Moving_Mean")
    print("from",move_date)
    moving_list=[]
    for t in range(step_all):
      my=move_date.year
      mm=move_date.month
      md=move_date.day
      mh=move_date.hour
      moving_time_str=f"{my}-{mm:02d}-{md:02d}_{mh:02d}:00:00"
  
      
      moving_file=Dataset(f"{wrfout_dir}/{ex}/n{member}/wrfout_d01_{moving_time_str}")
      moving_list.append(moving_file)
      
      move_date+=timedelta(hours=dh)
#ココで移動平均内での時間更新      
    print("to", move_date - timedelta(hours=dh))
    
##変数取り出し
    if var_type == "surf":
      if var_unit != 0:
        var=getvar(moving_list, var_name ,timeidx=ALL_TIMES,method="join", units=var_unit)
      elif var_unit == 0:
        var=getvar(moving_list, var_name ,timeidx=ALL_TIMES,method="join")
      else:
        raise ValueError("Not Set Unit")

    if var_type == "pres":
      if var_unit != 0:
        var_all=getvar(moving_list, var_name ,timeidx=ALL_TIMES,method="join", units=var_unit)
      elif var_unit == 0:
        var_all=getvar(moving_list, var_name ,timeidx=ALL_TIMES,method="join")
      else:
        raise ValueError("Not Set Unit")

      p_3d=getvar(moving_list,"p",timeidx=ALL_TIMES,method="join",units="hpa")
      var=interplevel(var_all, p_3d, p)

##ファイルclose
    for ds in moving_list:
      ds.close()
  
##変数の移動平均(bar:ここで次元圧縮:file*ny*nx -> 1*ny*nx)
    var_mmean=var.mean(dim="file",skipna=True)
    mlist_bar.append(var_mmean)

##改めて対称時刻のwrfoutを開いて偏差を計算する(prime)
    current_file=Dataset(f"{wrfout_dir}/{ex}/n{member}/wrfout_d01_{time_str}")
    if var_type == "surf":
      if var_unit != 0:
        current_var=getvar(current_file, var_name , units=var_unit)
      elif var_unit == 0:
        current_var=getvar(current_file, var_name)
      else:
        raise ValueError("Not Set Unit")

    if var_type == "pres":
      if var_unit != 0:
        current_var_all=getvar(current_file, var_name ,units=var_unit)
      elif var_unit == 0:
        current_var_all=getvar(current_file, var_name)
      else:
        raise ValueError("Not Set Unit")

      p_3d=getvar(current_file, "p", units="hpa")
      var=interplevel(current_var_all, p_3d, p)

    var_anom=var - var_mmean
    mlist_prime.append(var_anom)

    date+=timedelta(hours=dh)

  
  mdarray_bar=xr.concat(mlist_bar, dim="time")
  mdarray_prime=xr.concat(mlist_prime, dim="time")

  slist_bar.append(mdarray_bar)
  slist_prime.append(mdarray_prime)


sdarray_bar=xr.concat(slist_bar, dim="member")
sdarray_prime=xr.concat(slist_prime, dim="member")

sdarray_bar=sdarray_bar.rename(f"{var_name}_bar")
sdarray_prime=sdarray_prime.rename(f"{var_name}_prime")

output=xr.merge([sdarray_bar, sdarray_prime])
output = output.assign_coords(
    member=("member", target_members)
)
#output_save=output.transpose("member", "time", "south_north", "west_east")
#->datasetの次元の順序入れ替え，反応せず

if output.sizes["member"] == 1:   
  output.squeeze(dim="member")
#アンサンブルナシの場合member軸を削除
output.to_netcdf(saving_name)
print(output)

print(f"Save {saving_name}")
print("End Program")