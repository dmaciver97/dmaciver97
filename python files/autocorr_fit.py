import argparse
import csv


parser = argparse.ArgumentParser()
parser.add_argument('--label', type=str, default='2.5A.TDdat', help='TDdat text file')
args = parser.parse_args()

label = args.label

def row_reader(row):
    print(row)
    for num in row:
        lst = (num.split(' '))
        if len(num)==37:
            item1= lst[0]
            item2=lst[1]
        else:
            item3=lst[1]
            item4=lst[2]

lead_path = 'C:\\Users\\Administrator\\programs\\git repos\\dmaciver97\\python files\\TDdat files\\'
#material, power, size, radius= get_label_path(label)

path = lead_path+label

with open(path, 'r') as in_file:
    reader = csv.reader(in_file.readlines()[5:])
    
    with open('log.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerows(row_reader(row) for row in reader)