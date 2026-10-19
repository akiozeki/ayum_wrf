##指定領域のフルード数(Fr=U/NH)の領域平均値の高度分布を計算する
from datetime import datetime,timedelta
from wrf import getvar,interplevel,get_cartopy,latlon_coords
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from myFunction.vertical_cross_section import latlon_to_yx

##定数
H=1000
#山の高さスケール
g=9.81

##実験設定
case="2025winter"
ex="CTL"

map_proj="lambert"
#"lambert" or "mercator" 
ref_lat=42.5
#メルカトルであれば使用する

start_date=datetime(2024,12,1,0)
end_date=datetime(2024,12,31,18)
dh=3
dsec=dh*3600
print(f"Set Time : {start_date}----{end_date} dh={dh}")

wrf_dir="/DATA/USER/ozeki/MyWRF"
wrfout_dir=f"{wrf_dir}/wrfout"
fig_dir=f"{wrf_dir}/fig/{case}"

##領域設定
wrfinput=f"{wrfout_dir}/{ex}/n1/wrfinput_d01"
sw_lat, sw_lon=
ne_lat, ne_lon=
sw_point=latlon_to_yx
ne_point=latlon_to_yx

date=start_date
while date <= end_date:
  y=date.year
  m=date.month
  d=date.day
  h=date.hour
  time_str=f"{y}-{m:02d}-{d:02d}_{h:02d}:00:00"
#   print(time_str)

##変数取り出し
  ds=Dataset(f"{wrfout_dir}/{case}/{ex}/n1/wrfout_d01_{time_str}")
  print(ds)