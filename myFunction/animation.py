#surf_bin.pyで出力されたpng画像データをつなげてgif形式に保存する
#参考https://emotionexplorer.blog.fc2.com/blog-entry-356.html
from PIL import Image
import glob
import pprint
import os

"""
pngファイルをつなげてgifアニメーションを作成する
fig_dir -> pngファイルがあるディレクトリ(アニメーションの出力先)
key -> アニメーションにするファイルの共通接頭辞
persec -> アニメーションの送り間隔(秒)
"""
def create_animation_from_png(fig_dir, key, persec):

  files=sorted(glob.glob(f"{fig_dir}/{key}*.png"))
  pprint.pprint(files)

  images=list(map(lambda file : Image.open(file),files))
  images[0].save(
    f"{fig_dir}/animation_{key}.gif",
    save_all=True,append_images=images[1:],optimize=True,
    duration=persec*1000,loop=0)

  for file in files:
    os.remove(file)
  
  print("Create Animation")
