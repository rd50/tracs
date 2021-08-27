#!/usr/bin/env python

import os,sys,re,traceback
from datetime import datetime
from string import Template
import time

i=sys.argv[1]
template_file = open(r'./W2_HPK3_1','r')
tmpl = Template(template_file.read())
lines = []
lines.append(tmpl.substitute(Vbias=i))
p = './Vbias'+str(i)+'V/'
if not os.path.exists(p):
    os.makedirs(p)
config_file_name = p + 'W2_HPK3_1'
config_file = open(config_file_name,'w')
config_file.writelines(lines)
config_file.close()
print('generate config file : ',config_file_name)
os.system('cp ./etct.carriers ' + p + 'etct.carriers')
os.system('./runall.sh '+p+' &')


