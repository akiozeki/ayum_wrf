from datetime import datetime, timedelta
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
from wrf import getvar,interplevel,get_cartopy,latlon_coords
from netCDF4 import Dataset

from myFunction.draw_field import draw_grid
from myFunction.animation import create_animation_from_png

# =========================
# 設定
# =========================
case = "2025winter"
ex1 = "CTL"
ex2 = "ME00"
lev=750

wrfout_dir = f"/DATA/USER/ozeki/MyWRF/wrfout/{case}"
fig_dir=f"/DATA/USER/ozeki/MyWRF/fig/{case}/Veq"

fname1 = f"{wrfout_dir}/{ex1}/n1/CTLVoleq{lev}hPa_dxdy25_endsigma6.nc"
fname2 = f"{wrfout_dir}/{ex2}/n1/ME00Voleq{lev}hPa_dxdy25_endsigma6.nc"

varname = "zeta_smooth"
ex="ex1"
#ex1 or ex2 or dif
key = f"{varname}_{ex}"
persec = 0.5
#1秒に何フレームかを表す

# コンターの設定
factor=5
levc_plus = np.arange(0.2,17,0.8)
levc_minus = np.arange(-16.2,0,0.8)

##座標系情報を得るために任意のwrfoutファイルを開く
sdate=datetime(2024,12,1,0)
sy=sdate.year
sm=sdate.month
sd=sdate.day
sh=sdate.hour
dh=3

coord_ds1=Dataset(f"{wrfout_dir}/{ex1}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")
coord_ds2=Dataset(f"{wrfout_dir}/{ex2}/n1/wrfout_d01_{sy}-{sm:02d}-{sd:02d}_{sh:02d}:00:00")

ter1=getvar(coord_ds1, "ter")
ter2=getvar(coord_ds2, "ter")
land=getvar(coord_ds1,"LANDMASK")
cart_proj=get_cartopy(land)
print("projection:",cart_proj)
lats,lons=latlon_coords(land)
#print(lats)
coord_ds1.close()
coord_ds2.close()


# =========================
# データ読み込み
# =========================
ds1 = xr.open_dataset(fname1)
ds2 = xr.open_dataset(fname2)

da1 = ds1[varname].transpose(
    "time", "south_north", "west_east"
)
da2 = ds2[varname].transpose(
    "time", "south_north", "west_east"
)

# 時刻・配列サイズの確認
if not np.array_equal(ds1["time"].values, ds2["time"].values):
  raise ValueError("CTLとME00の時刻が一致していません")

if da1.shape != da2.shape:
  raise ValueError("CTLとME00の配列サイズが一致していません")

# =========================
# 差を計算
# =========================
diff = da1 - da2


# コンター描画用データ
if ex == "ex1":
  data=da1 * 10**factor

elif ex == "ex2":
   data=da2 * 10**factor

elif ex == "dif":
  data = diff.values * 10**factor

else:
  raise ValueError("Data Not Detected")

times = ds1["time"].values
nt = len(times)

# =========================
# 描画領域
# =========================
fig = plt.figure(figsize=(10, 8))
date = sdate
for i in range(len(times)):

  fig = plt.figure(figsize=(10, 8))
  ax = fig.add_subplot(1, 1, 1, projection=cart_proj)

  ax.set_extent([120, 150, 30, 55], crs=ccrs.PlateCarree())
  ax.contour(
    lons, lats, data[i],
    levels=levc_plus,
    colors="red",
    linewidths=1.0,
    transform=ccrs.PlateCarree()
    )

    # 負の値：青
  ax.contour(
    lons, lats, data[i],
    levels=levc_minus,
    colors="blue",
    linewidths=1.0,
    transform=ccrs.PlateCarree()
    )

  draw_grid(ax)
  ax.coastlines()

  ax.set_title(f"{date}")

# PNG保存
  png_path = f"{fig_dir}/{key}{i:04d}.png"
  fig.savefig(png_path, dpi=150, bbox_inches="tight")

  plt.close(fig)

  print(f"Saved: {png_path}")

  date+=timedelta(hours=dh)

# PNG画像をGIFに変換
create_animation_from_png(
  fig_dir=fig_dir,
  key=key,
  persec=persec
)