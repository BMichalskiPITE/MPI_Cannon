import random
import sys

size = int(sys.argv[1])
for i in range(1,3):
	with open("matrix{}_{}.txt".format(size,i),"w+") as m:
		m.write("{} {} ".format(size, size))
		for _ in range(size*size):
			m.write("{} ".format(random.randint(0,10)))
