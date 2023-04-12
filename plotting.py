import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("output/FullSD3_VectorRotationData.csv")
# data = pd.read_csv("C:/Users/Arman/Downloads/FullSD3_VectorRotationData.csv")

# data = pd.read_csv("Z:/Undergrads/Arman 2022/SD3 AND NGLN MESH EXPORTS/SD3 3/SD3 3_VectorRotationData.csv")

data.reset_index(drop=True)

# data["Z Diff"] = data["Z Diff"].abs()
# data["Y Diff"] = data["Y Diff"].abs()

# print(data)

# sns.scatterplot(data=data, x="IOU", y="VectorDegrees")
# sns.set_style("ticks")

mpfiColor = ["#007167"]

sns.set_palette(sns.color_palette(mpfiColor))

h = sns.jointplot(x=data["IOU"], y=data["VectorAngle"], kind='scatter')

plt.xlabel('Intersection Over Union', fontsize=15)
plt.ylabel('Skew (Degrees)', fontsize=15)
h.fig.suptitle('Skew Degrees vs Intersection Over Union', fontsize=20)

h.fig.tight_layout()
h.fig.subplots_adjust(top=0.90)

h.fig.set_figwidth(8)
h.fig.set_figheight(8)

plt.show()