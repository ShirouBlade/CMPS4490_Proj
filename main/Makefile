all: fps prototype

prototype: proj_prototype.cpp
	g++ proj_prototype.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o prototype
	
fps: fps.cpp
	g++ fps.cpp timers.cpp libggfonts.a -Wall -lX11 -lGL -lGLU -lm -o fps
clean:
	rm -f prototype

