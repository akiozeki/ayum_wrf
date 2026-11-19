###Seo_et_al.(2021)で定義されたeady_growth_rateを計算する
###(高速化版: 全時刻を一度だけ読み込み、rollingで移動平均を計算)
from datetime import datetime,timedelta
from numbers import Number
from netCDF4 import Dataset
from wrf import getvar,interplevel
import numpy as np
import xarray as xr

#計算する気圧面
p=850
p_t=800
p_w=900

#ファイルパスと時間###
case="2025winter"
ex="CTL"
wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout/{case}"
output_dir=wrfout_dir
start_date=datetime(2024,12,1,0)
end_date=datetime(2025,2,28,18)
dh=3
print(f"Set Time : {start_date}----{end_date} dh={dh}")


##移動平均の設定##
moving_day=7
step_per_day=24 // dh
step_all=moving_day*step_per_day + 1
half_step=step_all // 2

#アンサンブルメンバー数 / アンサンブル平均したい場合はそのサイズ
ensmean="no"
#"yes" or "no" noの場合以下の要素数は必ず1個
target_members=[1]
#取り出すメンバ
if ensmean == "no" and len(target_members) != 1:
  raise ValueError("If ensmean=no, target_member's size =1")

if ensmean =="yes":
  ens_size=len(target_members)
  saving_name=f"{output_dir}/{ex}/{ex}_EGR_Seo{p}_{moving_day}dMean.nc"

elif ensmean == "no":
  ens_size=1
  saving_name=f"{output_dir}/{ex}/n{target_members[0]}/{ex}_EGR_Seo{p}_{moving_day}dMean.nc"

#######################
g=9.81
Rd=287
dx=25*1000
dy=25*1000


def clean_coords(v):
  """XLAT, XLONG と次元座標以外の座標(Time, XTIMEなど)を落とす"""
  drop=[c for c in v.coords if c not in v.dims and c not in ("XLAT", "XLONG")]
  return v.drop_vars(drop)


def load_fields(path):
  """1つのwrfoutから必要な2次元場(ny, nx)をまとめて取り出す"""
  nc=Dataset(path)
  try:
    p_3d=getvar(nc, "p", units="hpa")
    theta=getvar(nc, "theta", units="K")
    tk=getvar(nc, "tk")

    ds=xr.Dataset({
      "theta_t": clean_coords(interplevel(theta, p_3d, p_t)).load(),
      "theta_m": clean_coords(interplevel(theta, p_3d, p)).load(),
      "theta_w": clean_coords(interplevel(theta, p_3d, p_w)).load(),
      "temp":    clean_coords(interplevel(tk, p_3d, p)).load(),
    })
  finally:
    nc.close()
  return ds


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

slist=[]

for member in target_members:
  print(f"n={member}")

##全時刻のデータを一度だけ読み込む(各ファイル1回)
  fields=[]
  for date in all_dates:
    time_str=f"{date.year}-{date.month:02d}-{date.day:02d}_{date.hour:02d}:00:00"
    fields.append(load_fields(f"{wrfout_dir}/{ex}/n{member}/wrfout_d01_{time_str}"))

  ds_all=xr.concat(fields, dim="time", coords="minimal", compat="override")
  ds_all=ds_all.assign_coords(time=np.array(all_dates, dtype="datetime64[ns]"))
  del fields

##x,y座標(m)を付与(differentiate用)
  xy_coord = {
        "west_east": np.arange(ds_all.sizes["west_east"])*dx,
        "south_north": np.arange(ds_all.sizes["south_north"])*dy
        }
  ds_all=ds_all.assign_coords(xy_coord)

##移動平均: 中心化した窓、NaNは無視(元コードのskipna=Trueに対応)
  ds_mean=ds_all.rolling(time=step_all, center=True, min_periods=1).mean()
##窓が不完全な両端を除外
  ds_mean=ds_mean.isel(time=valid)
  del ds_all

  theta_t_mmean=ds_mean["theta_t"]
  theta_m_mmean=ds_mean["theta_m"]
  theta_w_mmean=ds_mean["theta_w"]
  temp_mmean=ds_mean["temp"]

###計算(全時刻まとめて)
#傾圧性
  dtheta_dx=theta_m_mmean.differentiate("west_east")
  dtheta_dy=theta_m_mmean.differentiate("south_north")
  grad_theta=np.sqrt(dtheta_dx**2 + dtheta_dy**2)

#静的安定度
  dtheta_dp=(theta_t_mmean - theta_w_mmean) / ((p_t - p_w)*100)
  const=(g**2 * p * 100) / Rd
  N=np.sqrt(const * (1 / (temp_mmean * theta_m_mmean)) * (-dtheta_dp))

#eady_growth_rate
  egr=0.31 * g / (N * theta_m_mmean) * grad_theta

  print("EGR Space&Time Mean = ",np.nanmean(egr.values) * 86400, "[Day-1]")

  ds_out=xr.Dataset({
    "dtheta_dx": dtheta_dx,
    "dtheta_dy": dtheta_dy,
    "B.V._Frequency": N,
    "Eady_Growth_Rate": egr,
  })
  slist.append(ds_out)
  del ds_mean, ds_out


output=xr.concat(slist, dim="member")
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