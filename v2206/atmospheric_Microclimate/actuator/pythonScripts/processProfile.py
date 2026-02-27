import numpy as np
import pandas as pd

df = pd.read_csv("data_uProfile.csv", usecols=["Points:1", "UMean:0", "UMean:1", "UMean:2"])
df.to_csv("uInletProfile.csv", index=False)

# counters = []
# y = []
# u = []
# yPrev = 0
# uSum = 0
# counter = 0
# for row in df.itertuples():
#     if row[1] == yPrev:
#         counter += 1
#         uSum += row[2]
#     else:
#         # log prev position
#         counters.append(counter)
#         uMean = uSum / counter
#         u.append(uMean)
#         y.append(yPrev)
#         # reset for current position
#         counter = 1
#         yPrev = row[1]
#         uSum = row[2]
    