import pandas as pd
import numpy as np

# create test dataframe
raw_data = {'name': ['Willard Morris', 'Al Jennings', 'Omar   Mullins', 'Spencer McDaniel'],
        'age': [20, 19, 22, 21],
        'favorite_color': ['blue', 'red', 'yellow', "green"],
        'grade': [88, 92, 95, 70]}
raw_df = pd.DataFrame(raw_data)

print(raw_df)

# make a list of all column names
colnames = raw_df.columns
wanted = ['age', 'grade']

# list comprehensions
final_colnames = [col for col in colnames if col in wanted] + [col for col in colnames if col not in wanted]

# feed it my list of names
raw_df = raw_df[final_colnames]
print(raw_df)
