#!/usr/bin/python

import mechanize
import sys

for f in sys.argv:
    if f != sys.argv[1]:
        br = mechanize.Browser()
        br.set_handle_robots(False)
        br.set_handle_refresh(False)
        br.open('http://dude.docking.org/generate')
        br.select_form(nr=0)
        
        br.set_all_readonly(False)
        br.form["email"] = "jss97@pitt.edu"
        br.form["turing"] = "never been better!"
        br.form.add_file(open(f), 'text/plain', f)
        req = br.submit()
        
        print req.read()
