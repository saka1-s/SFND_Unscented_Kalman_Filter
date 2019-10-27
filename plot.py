from matplotlib import pyplot
import numpy as np
import pandas as pd
rmse = pd.read_csv('rmse.csv')
pyplot.plot(rmse["t"],rmse["px"],label="px")
pyplot.plot(rmse["t"],rmse["py"],label="py")
pyplot.plot(rmse["t"],rmse["vx"],label="vx")
pyplot.plot(rmse["t"],rmse["vy"],label="vy")
pyplot.legend()
pyplot.show()