from LHCO_reader import LHCO_reader
import os
import matplotlib.pyplot as plt
background = "background/"
events1 = LHCO_reader.Events(f_name=background+"ttbar.lhco") # 100000
events2 = LHCO_reader.Events(f_name=background+"ttbar_largejet.lhco")

# print(len(events1)) # 100000
# print(len(events2)) # 100000

dirpath = os.path.join(os.path.dirname(__file__), "sphaleron")
path_list = []
for _, _, filename in os.walk(dirpath):
    path_list = ["{}\{}".format(dirpath,file) for file in filename]    

for path in path_list:
    print(path + ": ", len(LHCO_reader.Events(f_name=path)))

matrix = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
]

# Transpose the matrix using list comprehension
transposed_matrix = [[row[i] for row in matrix] for i in range(len(matrix[0]))]

# Print original matrix
print("Original Matrix:")
for row in matrix:
    print(row)

# Print transposed matrix
print("\nTransposed Matrix:")
for row in transposed_matrix:
    print(row)

fig, ax = plt.subplot(2, 3)

plt.scatter()