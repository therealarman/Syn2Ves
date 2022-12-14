import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("output/RotationData.csv")
data.reset_index(drop=True)

data["Z Diff"] = data["Z Diff"].abs()
data["Y Diff"] = data["Y Diff"].abs()

print(data)

sns.scatterplot(data=data, x="Z Diff", y="Y Diff", size="IOU")
plt.show()