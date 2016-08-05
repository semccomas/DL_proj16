#just made this to make sure each sliding table window has its own OHE, it does
import sys
f= open(sys.argv[1]).read().splitlines()
o= open(sys.argv[2]).read().splitlines()
count= 0
n= 0
for line in f:
	print line
	count= count +1
	if count == 20:
                print o[n]
		n= n+1
		print "--------------------------------------------------", n
                count= 0
