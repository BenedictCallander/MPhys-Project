import panadas as pd
MS1df = pd.read_csv('csv/tng33MAIN.csv')

MS1df.met = 10*((MS1df.met-MS1df.met.min())/(MS1df.met.max()-MS1df.met.min()))
