#!/usr/bin/env python

import os,sys,re,traceback
from datetime import datetime
from string import Template
import time

os.system('rm *um -rf')

Vbias = [80,90,100,110,120,130,140,150,160,170,180,190,200]
for V in Vbias:
    A = ' '
    for i in range(51):
        template_file = open(r'./W2_HPK3_1','r')
        tmpl = Template(template_file.read())
        lines = []
        lines.append(tmpl.substitute(Z_POSITION=i,Vbias='${Vbias}'))
        p = './'+str(i)+'um/'
        if not os.path.exists(p):
            os.makedirs(p)
        python_file_name = p + 'W2_HPK3_1'
        python_file = open(python_file_name,'w')
        python_file.writelines(lines)
        python_file.close()
        print('generate config file : ',python_file_name)
        if V == 80:
            os.system('cp ./etct.carriers ' + p + 'etct.carriers')
            os.system('cp ./config_Vbias.py ' + p + 'config_Vbias.py')
            os.system('cp ./runall.sh ' + p + 'runall.sh')
        os.system('cd ' + p +' && python '+'config_Vbias.py '+str(V) +' -print >'+'run'+str(i)+'um'+ '.dt 2>&1 &')

    for j in range(51):
        B = str(j) + 'um/'+'Vbias'+str(V)+'V/NOirrad_dt1ps_4pF_tNOtrappingns_dz1um_dy5dV20V_0nns_edge_0_rc.hetct.root'
        A = A + B + ' '

    index = 0
    while index == 0:
        time.sleep(5)
        index = 1
        for j in range(51):
            D = str(j) + 'um/'+'Vbias'+str(V)+'V/NOirrad_dt1ps_4pF_tNOtrappingns_dz1um_dy5dV20V_0nns_edge_0_rc.hetct.root'
            if not os.path.exists(D):
                index = 0
        if index == 1:
            os.system('./mergeall.sh '+ 'edge_lgad_'+str(V)+'V.root' +A)


