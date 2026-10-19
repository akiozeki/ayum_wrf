import numpy as np
from scipy import signal
from scipy.stats import t
import math

"""
時系列データ(x)の平均, 標準偏差, サンプル数を計算する
"""
def calc_statistics(x):
  x_bar=np.mean(x)
  sigma=np.std(x)
  N=len(x)
  #r=np.corrcoef(x[0:-1], x[1:])[0,1]
  #x[0:-1]は配列の初期要素～最終1つ手前，x[1:]は2要素目～最終
  #N_star=N*(1-r)/(1+r)
#  print("x_bar, sigma, r, N, N_star =")
#  print([x_bar, sigma, N, N_star])

  return x_bar, sigma, N


"""
2つの時系列データ(x1, x2)について母平均の差のt検定を行う
出力は有意(True) or 有意でない(False)と有効自由度
3つ目の引数は信頼係数 両側95% -> 0.975, 90% -> 0.95
"""
def ttest(x1, x2, cc):
  x_bar1, sigma1, N1=calc_statistics(x1)
  x_bar2, sigma2, N2=calc_statistics(x2)

  t_sample=(x_bar1 - x_bar2) / (np.sqrt(sigma1**2/N1 + sigma2**2/N2))
#検定量tが求まる

  dof=N1 + N2 - 2
  t_border=t.ppf(cc, dof)
#棄却域が求まる

  test_result=(np.abs(t_sample) > t_border)

  return test_result

#実行テスト
# a1=[0.183, 0.927, 0.451, 0.068, 0.792, 0.514, 0.339, 0.875, 0.226, 0.641]
# a2=[0.58, 0.104, 0.886, 0.473, 0.295, 0.967, 0.532, 0.147, 0.681, 0.359]
# print("a1=", a1)
# print("a2=", a2)
# print(ttest(a1, a2, 0.975))


"""
2つの(時間, 南北, 東西)の3次元配列(z1, z2)について各格子点(南北, 東西)で母平均の差のt検定を行う
出力は, 有意->1, 有意でない->0の2次元マップ
"""
def ttest_2Dmap(z1, z2, cc):
  nt, nj, ni=z1.shape
  print(f"sample size : (nt, ni, nj)=",(nt, nj, ni))

  confidence_coefficient=cc

  sign_map=np.zeros((nj,ni))
  #dof_map=np.zeros((nj,ni))

  for i in range(0, ni):
    for j in range(0, nj):
      x1=z1[:, j, i]
      x2=z2[:, j, i]

      result=ttest(x1, x2, confidence_coefficient)
      sign_map[j, i]=result.astype(int)
 #     dof_map[j, i]=dof
 # print("有効自由度の領域平均：", np.nanmean(dof_map))
  print(f"有意 (片側 {confidence_coefficient}%) な格子点/全格子点数:{np.sum(sign_map)}/{sign_map.size}")
  return sign_map


"""
2つの(時間, 鉛直, 水平)の3次元配列(z1, z2)について各格子点(鉛直, 水平)で母平均の差のt検定を行う
変数名以外はttest_2Dmapと同じ
出力は, 有意->1, 有意でない->0の2次元マップ
"""
def ttest_vertical_cross(z1, z2, cc):
  nt, nk, nh=z1.shape
  print(f"sample size : (nt, nk, nh)=",(nt, nk, nh))

  confidence_coefficient=cc

  sign_map=np.zeros((nk,nh))
  dof_map=np.zeros((nk,nh))

  for h in range(0, nh):
    for k in range(0, nk):
      x1=z1[:, k, h]
      x2=z2[:, k, h]

      result=ttest(x1, x2, confidence_coefficient)
      sign_map[k, h]=result.astype(int)
  #    dof_map[k, h]=dof
  #print("有効自由度の領域平均：", np.nanmean(dof_map))
  print(f"有意 (片側 {confidence_coefficient}%) な格子点/全格子点数:{np.sum(sign_map)}/{sign_map.size}")
  
  return sign_map