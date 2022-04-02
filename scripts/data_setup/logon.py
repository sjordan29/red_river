from pyesgf.logon import LogonManager # if getting keyerror need to change 'HOME' to 'USERPROFILE' if using windows in \lib\site-packages\pyesgf\logon.py; line 58, 62
lm = LogonManager()
lm.logoff()
lm.logon_with_openid('https://esgf-data.dkrz.de/esgf-idp/openid/smj5vup', password='76Crooked!', bootstrap=True)
# lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
print(lm.is_logged_on())
