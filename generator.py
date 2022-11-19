import sys, random

numtogen = int(sys.argv[1])

print(f"count={numtogen}")
for i in range(numtogen):
    print(f"{i+1} {random.randrange(100, 990)} {random.randrange(100, 990)}")