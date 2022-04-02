# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:26:03 2021

@author: Sarah

Run in 
"""
from pyesgf.logon import LogonManager # if getting keyerror need to change 'HOME' to 'USERPROFILE' if using windows in \lib\site-packages\pyesgf\logon.py; line 58, 62

lm = LogonManager()
lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
lm.is_logged_on()