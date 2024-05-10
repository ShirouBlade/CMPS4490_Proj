all: aerial

aerial: aerial.cpp
	g++ aerial.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o aerial

clean:
	rm -f aerial

