#!/usr/bin/env python

import os,sys,re,traceback
from datetime import datetime
from string import Template
A = [0.005,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045]
for i in A:
	for j in range(13):
		template_file = open(r'./plt_waf_lgad.c','r')
		tmpl = Template(template_file.read())
		lines = []
		q =  i + 11.953
		a = 80 + j * 10
		lines.append(tmpl.substitute(VBIAS=a,Z=i,Z_DATA=q))
		python_file_name = 'temp_plt_waf_lgad.c'
		python_file = open(python_file_name,'w')
		python_file.writelines(lines)
		python_file.close()
		os.system('root -q -b -l temp_plt_waf_lgad.c')
