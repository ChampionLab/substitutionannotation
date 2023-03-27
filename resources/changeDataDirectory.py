#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
package_dir = os.path.dirname(os.path.split(os.path.abspath(__file__))[0])
file_path = os.path.join(package_dir, 'temp', 'data_directory.txt')

with open(file_path, "w") as f:
    inputdir = input('Please enter the directory with the data:')
    f.write(inputdir)

print('Done')