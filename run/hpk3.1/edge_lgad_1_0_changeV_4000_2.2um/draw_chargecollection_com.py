#!/usr/bin/env python

import os,sys,re,traceback
from datetime import datetime
from string import Template

A = [80,100,120,140,160,180,200]
for i in A:
	template_file = open(r'./plot_charge_collection.c','r')
	tmpl = Template(template_file.read())
	lines = []
	lines.append(tmpl.substitute(VBIAS=i))
	python_file_name = 'temp_plot_charge_collection.c'
	python_file = open(python_file_name,'w')
	python_file.writelines(lines)
	python_file.close()
	os.system('root -q -b -l temp_plot_charge_collection.c')
