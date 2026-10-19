import numpy as np
import xarray as xr

"""
指定した日数の移動平均を計算する。
da -> xarray.DataArray, time次元を持つデータ
days -> 移動平均の期間 [日]
center -> 
    True  : 移動平均の中心を元の時刻に合わせる
    False : 後方移動平均
"""

def moving_average_day(da, days, center=True):
 
    # 時間間隔を確認
    dt = da.time.diff("time").median().values
    dt_hours = dt / np.timedelta64(1, "h")

    # 指定日数をデータ点数に変換
    window = int(round(days * 24 / dt_hours))

    if window < 1:
        raise ValueError("移動平均の期間が短すぎます")

    # 移動平均
    da_ma = da.rolling(
        time=window,
        center=center,
        min_periods=window
    ).mean()

    return da_ma