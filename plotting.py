import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("output/FullSD3_VectorRotationData.csv")
data.reset_index(drop=True)

# data["Z Diff"] = data["Z Diff"].abs()
# data["Y Diff"] = data["Y Diff"].abs()

# print(data)

# sns.scatterplot(data=data, x="IOU", y="VectorDegrees")
# sns.set_style("ticks")
h = sns.jointplot(x=data["IOU"], y=data["VectorAngle"], kind='scatter')

h.set_axis_labels('Intersection Over Union', 'Skew (Degrees)')
h.fig.suptitle('Skew Degrees vs Intersection Over Union')

h.fig.tight_layout()
h.fig.subplots_adjust(top=0.90)

plt.show()