cluster: cluster.c
	gcc -std=c99 -Wall -Wextra -Werror -DNDEBUG cluster.c -o cluster -lm