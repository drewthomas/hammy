nwgks: nwgks.[ch] data.[ch] Makefile
	gcc -Os -s -Wall -Wextra -ansi -pedantic -std=c99 \
		-DMAKE_STANDALONE_NWGKS nwgks.c data.c -lm -o nwgks
