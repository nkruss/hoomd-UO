import numpy as np
import pandas as pd

def read_pos(fname: str):
    df = pd.read_csv(fname)
    print(df)
