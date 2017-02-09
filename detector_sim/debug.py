# detects any Non-ASCII characters in a file, or alternatively
# just add "# -*- coding: utf-8 -*-" to top of code

with open("pmt_res.py") as fp:
    for i, line in enumerate(fp):
        if "\xe2" in line:
            print i, repr(line)