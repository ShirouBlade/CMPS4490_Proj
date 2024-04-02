all: lab2 fps prototype

lab2: lab2.cpp timers.cpp
	g++ lab2.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o lab2
prototype: proj_prototype.cpp
	g++ proj_prototype.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o prototype
fps: fps.cpp
	g++ fps.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o fps
clean:
	rm -f lab2

