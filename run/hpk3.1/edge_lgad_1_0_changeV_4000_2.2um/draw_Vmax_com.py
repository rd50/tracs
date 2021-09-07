#!/usr/bin/env python

import os,sys,re,traceback
from datetime import datetime
from string import Template

A = [80,90,100,110,120,130,140,150,160,170,180,190,200]
for i in A:
	template_file = open(r'./plt_Vmax.c','r')
	tmpl = Template(template_file.read())
	lines = []
	lines.append(tmpl.substitute(VBIAS=i))
	python_file_name = 'temp_plt_Vmax.c'
	python_file = open(python_file_name,'w')
	python_file.writelines(lines)
	python_file.close()
	os.system('root -q -b -l temp_plt_Vmax.c')
