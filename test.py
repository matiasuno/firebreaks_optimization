import pandas as pd

df = pd.read_csv("ignitions_log.csv", header=None, names=['points'])

# Now, you can access the single column DataFrame
single_column_df = df['points']

# Convert the DataFrame column to a list
column_values_list = single_column_df.tolist()

# Print the list
print(column_values_list)