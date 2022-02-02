## Short Shell script that converts code to pretty-looking (ish) pdf


enscript -E -q -Z -p - -f Courier8 problem_1_hw4.py --color | ps2pdf - out.pdf
