import os
from glob import glob
import imageio

def save_as_gif():
	images = [imageio.imread(x) for x in glob(os.path.join(os.getcwd(), '*.png'))]
	framerate = float(1.0/30)
	imageio.mimsave(os.path.join(os.getcwd(), "spanish_castle_illusion.gif"), images, fps=framerate)

if __name__ == '__main__':
	save_as_gif()