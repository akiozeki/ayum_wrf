###変数の移動平均場および偏差場を計算する
###Claudeによる高速化版: 全時刻を一度だけ読み込み、rollingで移動平均を計算
from datetime import datetime,timedelta
from numbers import Number
from netCDF4 import Dataset
from wrf import getvar, interplevel
import numpy as np
import xarray as xr

######################################
#変数
var_name="ua"
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
ex="ME00"
wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
output_dir=wrfout_dir
start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
#dh=12
print(f"Set Time : {start_date}----{end_date} dh={dh}")


##移動平均の設定##
moving_day=2
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


def load_field(path):
  """1つのwrfoutから対象変数の2次元場(ny, nx)を取り出す"""
  nc=Dataset(path)
  try:
    kw={} if var_unit == 0 else {"units": var_unit}

    if var_type == "surf":
      v=getvar(nc, var_name, **kw)

    elif var_type == "pres":
      v_all=getvar(nc, var_name, **kw)
      p_3d=getvar(nc, "p", units="hpa")
      v=interplevel(v_all, p_3d, p)

    else:
      raise ValueError("var_type must be 'surf' or 'pres'")

    #時刻ごとに値が違うスカラー座標は連結時に邪魔になるので落とす
    v=v.drop_vars(["Time", "XTIME"], errors="ignore").load()
  finally:
    nc.close()
  return v


#全時刻のリスト(start_date〜end_date)
all_dates=[]
date=start_date
while date <= end_date:
  all_dates.append(date)
  date+=timedelta(hours=dh)

n_time=len(all_dates)
if n_time <= 2*half_step:
  raise ValueError("期間が移動平均の窓に対して短すぎます")

#移動平均の窓が期間内に収まる時刻だけを使う(元コードのSkip条件と同じ)
valid=slice(half_step, n_time - half_step)
valid_dates=all_dates[valid]
print(f"Valid Time : {valid_dates[0]}----{valid_dates[-1]} ({len(valid_dates)} steps)")

slist_bar=[]
slist_prime=[]

for member in target_members:
  print(f"n={member}")

##全時刻のデータを一度だけ読み込む(各ファイル1回)
  fields=[]
  for date in all_dates:
    time_str=f"{date.year}-{date.month:02d}-{date.day:02d}_{date.hour:02d}:00:00"
    fields.append(load_field(f"{wrfout_dir}/{ex}/n{member}/wrfout_d01_{time_str}"))

  var_all=xr.concat(fields, dim="time", coords="minimal", compat="override")
  var_all=var_all.assign_coords(time=np.array(all_dates, dtype="datetime64[ns]"))
  del fields

##移動平均(bar): 中心化した窓、NaNは無視(元コードのskipna=Trueに対応)
  var_bar=var_all.rolling(time=step_all, center=True, min_periods=1).mean()

##窓が不完全な両端を除外
  var_bar=var_bar.isel(time=valid)
  var_cur=var_all.isel(time=valid)

##偏差(prime)
  var_prime=var_cur - var_bar

  slist_bar.append(var_bar)
  slist_prime.append(var_prime)
  del var_all, var_bar, var_cur, var_prime


sdarray_bar=xr.concat(slist_bar, dim="member")
sdarray_prime=xr.concat(slist_prime, dim="member")

sdarray_bar=sdarray_bar.rename(f"{var_name}_bar")
sdarray_prime=sdarray_prime.rename(f"{var_name}_prime")

output=xr.merge([sdarray_bar, sdarray_prime])
output = output.assign_coords(
    member=("member", target_members)
)

if output.sizes["member"] == 1:
  output=output.squeeze(dim="member", drop=True)
#アンサンブルナシの場合member軸を削除

##netCDFに保存できない型の属性(wrf-pythonのprojectionなど)を削除
valid_attr_types=(str, Number, np.ndarray, np.number, list, tuple, bytes)
def drop_invalid_attrs(attrs, label):
  for k in [k for k, v in attrs.items() if not isinstance(v, valid_attr_types)]:
    print(f"Drop attr '{k}' from {label}")
    del attrs[k]

drop_invalid_attrs(output.attrs, "dataset")
for name, v in output.variables.items():
  drop_invalid_attrs(v.attrs, name)

output.to_netcdf(saving_name)
print(output)

print(f"Save {saving_name}")
print("End Program")