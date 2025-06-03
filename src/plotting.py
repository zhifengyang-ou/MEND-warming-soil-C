import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fp = './photosyn/'
df = pd.read_csv(fp+'cflux.out')

fig, axs = plt.subplots(2, 1)
axs[0].plot(df.Ac)
#axs[0].set_xlim(0, 2)
#axs[0].set_xlabel('Time')
#axs[0].set_ylabel('s1 and s2')
axs[0].grid(True)

axs[1].plot(df['An'])
axs[1].grid(True)
fig.tight_layout()
plt.show()
