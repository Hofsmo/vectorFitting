import pandas as pd
import numpy as np
from datetime import datetime

df = pd.read_csv('data.csv')

# Timeformatstring
timef = "%Y/%m/%d %H:%M:%S.%f"
t0 = datetime.strptime(df.Timestamp[0], timef)

t = [(datetime.strptime(Timestamp, timef)-t0).total_seconds()
     for Timestamp in df.Timestamp]
