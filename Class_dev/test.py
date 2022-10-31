import pandas as pd

galaxy_df = pd.read_csv("test1.csv")
valid_galaxies = galaxy_df[galaxy_df['mass']<9.5]
valid_galaxies.to_csv("verify.csv")