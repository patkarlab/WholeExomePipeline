#!/usr/bin/env python3

import os
import sys
import pandas as pd
import csv
import openpyxl


input_xlsx_file = sys.argv[1]
csv_file1 = sys.argv[2]

wb = openpyxl.load_workbook(input_xlsx_file)
if os.path.getsize(csv_file1) != 0:
	worksheet = wb.create_sheet('pharma_markers')
	read_csv1 = csv.reader(open(csv_file1, 'r', encoding='utf-8'), delimiter='\t')
	for rows in read_csv1:
			worksheet.append(rows)

wb.save(input_xlsx_file)
