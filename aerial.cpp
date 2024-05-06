//program: proj_prototype.cpp
//author:  Danny Simpson
//date:    April 02, 2024
//
//
//
//Prototype for aerial obstacle course.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
//X11 utilities not currently needed.
//#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "defs.h"
#include "log.h"
#include "fonts.h"
#include <unistd.h>
#include <time.h>


typedef float Flt;
typedef Flt Vec[3];
typedef Flt Vec4[4];
typedef Flt Matrix[4][4];

#define M_PI 3.14159265358979323846
//some macros
const Vec upv = {0.0, 1.0, 0.0};
const int MAX_SMOKES = 400;
//Class for a vector object.
class Myvec {
public:
	Flt x, y, z;
	Myvec(Flt a, Flt b, Flt c) {
		x = a;
		y = b;
		z = c;
	}
	void make(Flt a, Flt b, Flt c) {
		x = a;
		y = b;
		z = c;
	}
	void negate() {
		x = -x;
		y = -y;
		z = -z;
	}
	void zero() {
		x = y = z = 0.0f;
	}
	Flt dot(Myvec v) {
		return (x*v.x + y*v.y + z*v.z);
	}
	Flt lenNoSqrt() {
		return (x*x + y*y + z*z);
	}
	Flt len() {
		return sqrtf(lenNoSqrt());
	}
	void copy(Myvec b) {
		b.x = x;
		b.y = y;
		b.z = z;
	}
	void add(Myvec b) {
		x = x + b.x;
		y = y + b.y;
		z = z + b.z;
	}
	void sub(Myvec b) {
		x = x - b.x;
		y = y - b.y;
		z = z - b.z;
	}
	void scale(Flt s) {
		x *= s;
		y *= s;
		z *= s;
	}
	void addS(Myvec b, Flt s) {
		x = x + b.x * s;
		y = y + b.y * s;
		z = z + b.z * s;
	}
    // Function to find
    // cross product of two vector array.
    //void crossProduct(int vect_A[], int vect_B[], int cross_P[]){
    //    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    //    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    //    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    //}
	Flt normalize() {
		return 0.0f;
	}
};
class Smoke {
    public:
        Vec pos;
        Vec vert[16];
        Flt radius;
        int n;
        float Distance;
        struct timespec tstart;
        Flt maxtime;
        Flt alpha;
        Smoke() { }
};
class Rings{
    public:
        Vec pos;
        Rings() { }
};
class Cloud {
    public:
        int x;
        int y;
        int random;
        float cloudmap32[32*32];
        float cloudmap256[256*256];
        Flt alpha;
};
//-----------------------------------------------------------------------------
//Setup timers
const double physicsRate = 1.0 / 60.0;
const double oobillion = 1.0 / 1e9;
extern struct timespec timeStart, timeCurrent;
extern struct timespec timePause;
extern double physicsCountdown;
extern double timeSpan;
extern double timeDiff(struct timespec *start, struct timespec *end);
extern void timeCopy(struct timespec *dest, struct timespec *source);
//-----------------------------------------------------------------------------





class Global {
public:
	int xres, yres;
    int menu, pause, optionSelected;
	Flt aspectRatio;
	Vec cameraPosition;
    float initialz;
    float cameraDistancePlayer;
    float angleAroundPlayer;
    float pitch;
    float yaw;
    float roll;
    Matrix cameraMat;
    Vec cameraDir;
	GLfloat lightPosition[4];
    struct timespec smokeStart, smokeTime;
    struct timespec cloudStart, cloudTime;
    struct timespec introStart, introPause, introTime;
    Smoke *smoke;
    Smoke *cloud;
    Rings *ring;
    Cloud *background;
    int nsmokes;
    int nrings;
    int nclouds;
    int vsync;
    int gamestart;
    float offset[3];
    int intro;
    int fps;
    Vec planePos;
    Vec plane2Pos;
    Vec planeAngle;
    Vec plane2Dir;
    Vec4 plane2Joystick;
    ~Global() {
        if (smoke)
        delete [] smoke;
        if (cloud)
        delete [] cloud;
    }
	Global() {
		//constructor
        pitch = 20;
        yaw = 0;
		xres=640;
		yres=480;
        menu = 0;
        pause = 0; // 0 -> Game not paused, 1 -> Game paused
        optionSelected = 0; // 0 - Resume, 1 - Options, 2 - Main Menu
		aspectRatio = (GLfloat)xres / (GLfloat)yres;
		MakeVector(0.0, 2.5, 15.5, cameraPosition);
        initialz = cameraPosition[2];
        MakeVector(0.0, 0.0, -1.0, cameraDir);
        MakeVector(0.0, 0.0, 0.0, plane2Dir);
        MakeVector(0.0, 0.0, 0.0, planeAngle);
		//light is up high, right a little, toward a little
		MakeVector(100.0f, 240.0f, 40.0f, lightPosition);
		lightPosition[3] = 1.0f;
        MakeVector(0.0, 1.0, 8.0, plane2Pos);
		//init_opengl();
        nsmokes = 0;
        nclouds = 0;
        nrings = 0;
        smoke = new Smoke[MAX_SMOKES];
        cloud = new Smoke[(MAX_SMOKES * 6)];
        ring = new Rings[MAX_SMOKES];
        background = new Cloud;
        gamestart = 0;
        offset[2] = 0.0f;
        intro = 1;
        vsync = 1;
        cameraDistancePlayer = 20;
        angleAroundPlayer = 0;
	}
	void init_opengl();
	void init();
	void check_mouse(XEvent *e);
	int check_keys(XEvent *e);
	void physics();
	void render();
} g;
// Function to find
// cross product of two vector array.
void crossProduct(const Vec v1, const Vec v2, Vec v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v3[0];
}

class X11_wrapper {
private:
	Window win;
	GLXContext glc;
public:
	Display *dpy;
	X11_wrapper() {
		//Look here for information on XVisualInfo parameters.
		//http://www.talisman.org/opengl-1.1/Reference/glXChooseVisual.html
		Window root;
		GLint att[] = { GLX_RGBA,
						GLX_STENCIL_SIZE, 2,
						GLX_DEPTH_SIZE, 24,
						GLX_DOUBLEBUFFER, None };
		Colormap cmap;
		XSetWindowAttributes swa;
		setup_screen_res(640, 480);
		dpy = XOpenDisplay(NULL);
		if (dpy == NULL) {
			printf("\n\tcannot connect to X server\n\n");
			exit(EXIT_FAILURE);
		}
		root = DefaultRootWindow(dpy);
		XVisualInfo *vi = glXChooseVisual(dpy, 0, att);
		if (vi == NULL) {
			printf("\n\tno appropriate visual found\n\n");
			exit(EXIT_FAILURE);
		} 
		cmap = XCreateColormap(dpy, root, vi->visual, AllocNone);
		swa.colormap = cmap;
		swa.event_mask = ExposureMask | KeyPressMask | KeyReleaseMask |
							StructureNotifyMask | SubstructureNotifyMask;
		win = XCreateWindow(dpy, root, 0, 0, g.xres, g.yres, 0,
								vi->depth, InputOutput, vi->visual,
								CWColormap | CWEventMask, &swa);
		set_title();
		glc = glXCreateContext(dpy, vi, NULL, GL_TRUE);
		glXMakeCurrent(dpy, win, glc);
	}
	~X11_wrapper() {
		XDestroyWindow(dpy, win);
		XCloseDisplay(dpy);
	}
	void setup_screen_res(const int w, const int h) {
		g.xres = w;
		g.yres = h;
		g.aspectRatio = (GLfloat)g.xres / (GLfloat)g.yres;
	}
	void check_resize(XEvent *e) {
		//The ConfigureNotify is sent by the
		//server if the window is resized.
		if (e->type != ConfigureNotify)
			return;
		XConfigureEvent xce = e->xconfigure;
		if (xce.width != g.xres || xce.height != g.yres) {
			//Window size did change.
			reshape_window(xce.width, xce.height);
		}
	}
	void reshape_window(int width, int height) {
		//window has been resized.
		setup_screen_res(width, height);
		//
		glViewport(0, 0, (GLint)width, (GLint)height);
		glMatrixMode(GL_PROJECTION); glLoadIdentity();
		glMatrixMode(GL_MODELVIEW); glLoadIdentity();
		glOrtho(0, g.xres, 0, g.yres, -1, 1);
		set_title();
	}
	void set_title() {
		//Set the window title bar.
		XMapWindow(dpy, win);
		XStoreName(dpy, win, "fps framework");
	}
	bool getXPending() {
		return XPending(dpy);
	}
	XEvent getXNextEvent() {
		XEvent e;
		XNextEvent(dpy, &e);
		return e;
	}
	void swapBuffers() {
		glXSwapBuffers(dpy, win);
	}
} x11;

void OverlapOctaves(float  *map32, float  *map256);
void ExpFilter(float  *map);
int main()
{
	g.init_opengl();
	srand(time(NULL));
	int done = 0;
    int fps = 0;
    g.menu = 0;
    clock_gettime(CLOCK_REALTIME, &timePause);
    clock_gettime(CLOCK_REALTIME, &timeStart);
    struct timespec fpsStart;
    struct timespec fpsCurr;
    physicsCountdown = 0.0;
    clock_gettime(CLOCK_REALTIME, &fpsStart);
	while (!done) {
		while (x11.getXPending()) {
			XEvent e = x11.getXNextEvent();
			x11.check_resize(&e);
			g.check_mouse(&e);
			done = g.check_keys(&e);
		}
        clock_gettime(CLOCK_REALTIME, &timeCurrent);
        timeSpan = timeDiff(&timeStart, &timeCurrent);
        timeCopy(&timeStart, &timeCurrent);
	physicsCountdown += timeSpan;
	    while (physicsCountdown >= physicsRate){
            g.physics();
            physicsCountdown -= physicsRate;
        }
        ++fps;
        clock_gettime(CLOCK_REALTIME, &fpsCurr);
        double diff = timeDiff(&fpsStart, &fpsCurr);
        if (diff >= 1.0) {
            g.fps = fps;
            fps = 0;
            timeCopy(&fpsStart, &fpsCurr);
        }
		g.render();
		x11.swapBuffers();
        OverlapOctaves(g.background->cloudmap32, g.background -> cloudmap256);
        ExpFilter(g.background -> cloudmap256);
	}
	cleanup_fonts();
	return 0;
}

void SetNoise(float  *map);
void Global::init() {
	//Place general program initializations here.
    SetNoise(background -> cloudmap32);
}

void Global::init_opengl()
{
	//OpenGL initialization
	glClearColor(0.0f, 0.1f, 0.3f, 0.0f);
	//Enable surface rendering priority using a Z-buffer.
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	//Enable Phong shading of surfaces.
	glShadeModel(GL_SMOOTH);
	//Enable this so material colors are the same as vertex colors.
	glEnable(GL_COLOR_MATERIAL);
	//Enable diffuse lighting of surfaces with normals defined.
	glEnable(GL_LIGHTING);
	//Turn on a light
	glLightfv(GL_LIGHT0, GL_POSITION, g.lightPosition);
	glEnable(GL_LIGHT0);
	//Do this to allow fonts
	glEnable(GL_TEXTURE_2D);
	initialize_fonts();
}


void Global::check_mouse(XEvent *e)
{
	//Did the mouse move?
	//Was a mouse button clicked?
	static int savex = 0;
	static int savey = 0;
	//
	if (e->type == ButtonRelease) {
		return;
	}
	if (e->type == ButtonPress) {
		if (e->xbutton.button==1) {
			//Left button is down
		}
		if (e->xbutton.button==3) {
			//Right button is down
		}
	}
	if (savex != e->xbutton.x || savey != e->xbutton.y) {
		//Mouse moved
		savex = e->xbutton.x;
		savey = e->xbutton.y;
	}
}
void identity33(Matrix m);
void yy_transform(const Vec rotate, Matrix a);
void trans_vector(Matrix mat, const Vec in, Vec out);
void matrixFromAxisAngle(const Vec v, Flt ang, Matrix m);
void screenShot();

int Global::check_keys(XEvent *e)
{
	//Was there input from the keyboard?
	if (e->type == KeyPress) {
		int key = (XLookupKeysym(&e->xkey, 0) & 0x0000ffff);
        if(key == XK_p){
            pause ^= 1;
            return 0;
        }
        float speed = 0.5;
        float dist = 0.0;
        Vec up = {0.0, 1.0, 0.0};
        const float sensitivity = 0.1f;
        //int qdown = 0;
        if (key == XK_p) {
            pause ^= 1; // Toggle pause state
            return 0;
        }
        if (pause) {
            if (key == XK_Down) {
                optionSelected = (optionSelected + 1) % 3; // Cycle through options
            }
            if (key == XK_Up) {
                optionSelected = (optionSelected + 2) % 3; // Cycle through options 
            }
            if (key == XK_Return || key == XK_space) {
                switch (optionSelected) {
                    case 0:
                        pause = 0; // Resume game
                        break;
                    case 1:
                        // Option menu
                        break;
                    case 2: // Main menu
                        {
                            gamestart = 0;
                            pause = 0;
                            break;
                        }
                }
            }
            return 0; // ignore other keys when paused
        } else{
		switch(key) {
			case XK_1:
				break;
            case XK_m:
                g.gamestart = 1;
                break;
			case XK_Right:
				//g.cameraPosition[0] += 1.0;
                {    
                    Vec v = {0.0, -0.1, 0.0};
                    Matrix m;
                    identity33(m);
                    yy_transform(v, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                    //trans_vector(m, g.plane2Pos, g.plane2Pos);
                    g.planeAngle[0]--;
                    if(g.planeAngle[1] < 30)
                        g.planeAngle[1]++;
                    //Vec v2 = {0.0, -0.01, 0.0};
                    //Matrix n;
                    //identity33(n);
                    //yy_transform(v2, n);
                    //trans_vector(n, g.cameraPosition, g.cameraPosition);
                }
				break;
			case XK_Left:
				//g.cameraPosition[0] -= 1.0;
                {
                    Vec v = {0.0, 0.1, 0.0};
                    Matrix m;
                    identity33(m);
                    yy_transform(v, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                    //trans_vector(m, g.plane2Dir, g.plane2Dir);
                    //trans_vector(m, g.plane2Pos, g.plane2Pos);
                    g.planeAngle[0]++;
                    if(planeAngle[1] > -30)
                        g.planeAngle[1]--;
                    //Vec v2 = {0.0, 0.01, 0.0};
                    //Matrix n;
                    //identity33(n);
                    //yy_transform(v2, n);
                    //trans_vector(n, g.cameraPosition, g.cameraPosition);
                }
				break;
			case XK_Up:
				//g.cameraPosition[1] += 0.2;
                //tilt camera down
                //need to rotate camera vector
                {
                    //Vec v = {0.1, 0.0, 0.0};
                    Matrix m;
                    //identity33(m);
                    //yy_transform(v, m);
                    Vec x;
                    crossProduct(g.cameraDir, up, x);
                    g.plane2Dir[0] = g.cameraDir[0];
                    g.plane2Dir[1] = g.cameraDir[1];
                    g.plane2Dir[2] = g.cameraDir[2];
                    matrixFromAxisAngle(x, 0.1, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                    //trans_vector(m, g.plane2Dir, g.plane2Dir);
                    g.planeAngle[2]++;
                    //Matrix n;
                    //Vec y;
                    //crossProduct(g.cameraPosition, up, y);
                    //matrixFromAxisAngle(y, 0.001, n);
                    //trans_vector(n, g.cameraPosition, g.cameraPosition);
                }
				break;
			case XK_Down:
				//g.cameraPosition[1] -= 0.2;
                {
                    //Vec v = {-0.1, 0.0, 0.0};
                    Matrix m;
                    //identity33(m);
                    //yy_transform(v, m);
                    Vec x;
                    crossProduct(g.cameraDir, up, x);
                    g.plane2Dir[0] = g.cameraDir[0];
                    g.plane2Dir[1] = g.cameraDir[1];
                    g.plane2Dir[2] = g.cameraDir[2];
                    x[0] = -x[0];
                    x[1] = -x[1];
                    x[2] = -x[2];
                    matrixFromAxisAngle(x, 0.1, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                    //trans_vector(m, g.plane2Dir, g.plane2Dir);
                    //Matrix n;
                    //Vec y;
                    //crossProduct(g.cameraPosition, up, y);
                    //matrixFromAxisAngle(y, 0.001, n);
                    //trans_vector(n, g.cameraPosition, g.cameraPosition);
                    g.planeAngle[2]--;
                }
				break;
			case XK_f:
				//g.cameraPosition[2] -= 1.0;
                // add the directio vector to the camera position vector
                for(int i=0; i < 3; i++){
                    g.plane2Pos[i] += g.cameraDir[i] * speed;
                    //g.cameraPosition[i] += g.cameraDir[i] * speed;
                }
				break;
			case XK_b:
				//g.cameraPosition[2] += 1.0;
                for(int i=0; i < 3; i++){
                    g.plane2Pos[i] -= g.cameraDir[i] * speed;
                    //g.cameraPosition[i] -= g.cameraDir[i] * speed;
                }
				break;
            case XK_w:
                //g.cameraPosition[1] += 0.2;
                Vec heading;
                heading[0] = 0.5 * sin((g.planeAngle[0] * M_PI)/180);
                heading[2] = 0.5 * cos((g.planeAngle[0] * M_PI)/180);
                heading[1] = -0.1 * g.planeAngle[2];
                for(int i=0; i < 3; i++){
                    g.plane2Pos[i] += heading[i];
                    //g.cameraPosition[i] += g.cameraDir[i] * speed;
                }
                if(g.planeAngle[1] > 0.0f)
                    g.planeAngle[1]--;
                if(g.planeAngle[1] < 0.0f)
                    g.planeAngle[1]++;
                break;
            case XK_s:
                //g.cameraPosition[1] -= 0.2;
                //for(int i=0; i < 3; i++){
                //    g.plane2Pos[i] -= g.cameraDir[i] * speed;
                //    g.cameraPosition[i] -= g.cameraDir[i] * speed;
                //}
                break;
            case XK_a:
                /*
                dist = -1.0;
                Vec left;
                crossProduct(g.cameraDir, up, left);
                for(int i = 0; i < 3; i++){
                    g.plane2Pos[i] += left[i] * dist;
                    //g.cameraPosition[i] += left[i] * dist;
                }
                g.planeAngle[0]--;
                */
                break;
            case XK_d:
                /*
                dist = 1.0;
                Vec right;
                crossProduct(g.cameraDir, up, right);
                for(int i = 0; i < 3; i++){
                    g.plane2Pos[i] += right[i] * dist;
                    //g.cameraPosition[i] += right[i] * dist;
                }
                g.planeAngle[0]++;
                */
                break;
            case XK_q:
                //qdown ^= qdown;
                //if(qdown != 0){
                //    speed = 10.0;
                //} else {
                //    speed = 0.1;
                //}
                system("convert -loop 0 -coalesce -layers OptimizeFrame -delay 20 ./images/img*.jpg abc.gif");
                break;
            case XK_space:
                g.plane2Pos[1] += 0.2;
                //g.cameraPosition[1] += 0.2;
                break;
            case XK_Tab:
                g.plane2Pos[1] -= 0.2;
                //g.cameraPosition[1] -= 0.2;
                break;
            case XK_p:
                screenShot();
                break;
	    case XK_v: {
		       g.vsync ^= 1;
		       static PFNGLXSWAPINTERVALEXTPROC glXSwapIntervalEXT = NULL;
                        glXSwapIntervalEXT =
                                (PFNGLXSWAPINTERVALEXTPROC)glXGetProcAddressARB(
                                (const GLubyte *)"glXSwapIntervalEXT");
                        GLXDrawable drawable = glXGetCurrentDrawable();
                        if (g.vsync) {
                                glXSwapIntervalEXT(x11.dpy, drawable, 1);
                        } else {
                                glXSwapIntervalEXT(x11.dpy, drawable, 0);
                        }
			break;
	    }
            case XK_l:
                g.menu ^= 1;
                break;
			case XK_Escape:
				return 1;
		}
	}
    }
	return 0;
}

float Noise(int x, int y, int random)
{
    int n = x + y * 57 + random * 131;
    n = (n<<13) ^ n;
    return (1.0f - ( (n * (n * n * 15731 + 789221) +
            1376312589)&0x7fffffff)* 0.000000000931322574615478515625f);
}

void SetNoise(float  *map)
{
  float temp[34][34];

  int random=rand() % 5000;

  for (int y=1; y<33; y++)
  for (int x=1; x<33; x++)
  {
    temp[x][y] = 128.0f + Noise(x,  y,  random)*128.0f;
  }

  for (int x=1; x<33; x++)
  {
    temp[0][x] = temp[32][x];
    temp[33][x] = temp[1][x];
    temp[x][0] = temp[x][32];
    temp[x][33] = temp[x][1];
  }
  temp[0][0] = temp[32][32];
  temp[33][33] = temp[1][1];
  temp[0][33] = temp[32][1];
  temp[33][0] = temp[1][32];

  for (int y=1; y<33; y++)
    for (int x=1; x<33; x++)
    {
      float center = temp[x][y]/4.0f;
      float sides = (temp[x+1][y] + temp[x-1][y] + temp[x][y+1] + temp[x][y-1])/8.0f;
      float corners = (temp[x+1][y+1] + temp[x+1][y-1] + temp[x-1][y+1] + temp[x-1][y-1])/16.0f;

      map[((x-1)*32) + (y-1)] = center + sides + corners;
    }
}
float Interpolate(float x, float y, float  *map)
{
  int Xint = (int)x;
  int Yint = (int)y;

  float Xfrac = x - Xint;
  float Yfrac = y - Yint;

  int X0 = Xint % 32;
  int Y0 = Yint % 32;
  int X1 = (Xint + 1) % 32;
  int Y1 = (Yint + 1) % 32;

  float bot = map[X0*32 + Y0] + Xfrac * (map[X1*32 + Y0] - map[X0*32 + Y0]);
  float top = map[X0*32 + Y1] + Xfrac * (map[X1*32 +  Y1] - map[X0*32 + Y1]);

  return (bot + Yfrac * (top - bot));
}

void OverlapOctaves(float  *map32, float  *map256)
{
  for (int x=0; x<256*256; x++)
  {
    map256[x] = 0;
  }

  for (int octave=0; octave<4; octave++)
    for (int x=0; x<256; x++)
      for (int y=0; y<256; y++)
      {
        float scale = 1 / pow(2, 3-octave);
        float noise = Interpolate(x*scale, y*scale , map32);

        map256[(y*256) + x] += noise / pow(2, octave);
      }
}

void ExpFilter(float  *map)
{
  float cover = 20.0f;
  float sharpness = 0.95f;

  for (int x=0; x<256*256; x++)
  {
    float c = map[x] - (255.0f-cover);
    if (c<0)     c = 0;
    map[x] = 255.0f - ((float)(pow(sharpness, c))*255.0f);
  }
}

void identity33(Matrix m)
{
	m[0][0] = m[1][1] = m[2][2] = 1.0f;
	m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] = m[2][1] = 0.0f;
}

void yy_transform(const Vec rotate, Matrix a)
{
	//This function applies a rotation to a matrix.
	//It actually concatenates a transformation to the matrix.
	//Call this function first, then call trans_vector() to apply the
	//rotations to an object or vertex.
	//
	if (rotate[0] != 0.0f) {
		Flt ct = cos(rotate[0]), st = sin(rotate[0]);
		Flt t10 = ct*a[1][0] - st*a[2][0];
		Flt t11 = ct*a[1][1] - st*a[2][1];
		Flt t12 = ct*a[1][2] - st*a[2][2];
		Flt t20 = st*a[1][0] + ct*a[2][0];
		Flt t21 = st*a[1][1] + ct*a[2][1];
		Flt t22 = st*a[1][2] + ct*a[2][2];
		a[1][0] = t10;
		a[1][1] = t11;
		a[1][2] = t12;
		a[2][0] = t20;
		a[2][1] = t21;
		a[2][2] = t22;
		return;
	}
	if (rotate[1] != 0.0f) {
		Flt ct = cos(rotate[1]), st = sin(rotate[1]);
		Flt t00 = ct*a[0][0] - st*a[2][0];
		Flt t01 = ct*a[0][1] - st*a[2][1];
		Flt t02 = ct*a[0][2] - st*a[2][2];
		Flt t20 = st*a[0][0] + ct*a[2][0];
		Flt t21 = st*a[0][1] + ct*a[2][1];
		Flt t22 = st*a[0][2] + ct*a[2][2];
		a[0][0] = t00;
		a[0][1] = t01;
		a[0][2] = t02;
		a[2][0] = t20;
		a[2][1] = t21;
		a[2][2] = t22;
		return;
	}
	if (rotate[2] != 0.0f) {
		Flt ct = cos(rotate[2]), st = sin(rotate[2]);
		Flt t00 = ct*a[0][0] - st*a[1][0];
		Flt t01 = ct*a[0][1] - st*a[1][1];
		Flt t02 = ct*a[0][2] - st*a[1][2];
		Flt t10 = st*a[0][0] + ct*a[1][0];
		Flt t11 = st*a[0][1] + ct*a[1][1];
		Flt t12 = st*a[0][2] + ct*a[1][2];
		a[0][0] = t00;
		a[0][1] = t01;
		a[0][2] = t02;
		a[1][0] = t10;
		a[1][1] = t11;
		a[1][2] = t12;
		return;
	}
}
void matrixFromAxisAngle(const Vec v, Flt ang, Matrix m)
{
    // arguments
    // v   = vector indicating the axis
    // ang = amount of rotation
    // m   = matrix to be updated
    // This source was used during research...
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/
    //
    struct Axisangle {
        Flt angle;
        Flt x,y,z;
    } a1;
    a1.x = v[0];
    a1.y = v[1];
    a1.z = v[2];
    a1.angle = ang;
    //
    Flt c = cos(a1.angle);
    Flt s = sin(a1.angle);
    Flt t = 1.0 - c;
    m[0][0] = c + a1.x * a1.x * t;
    m[1][1] = c + a1.y * a1.y * t;
    m[2][2] = c + a1.z * a1.z * t;
    //
    Flt tmp1 = a1.x * a1.y * t;
    Flt tmp2 = a1.z * s;
    m[1][0] = tmp1 + tmp2;
    m[0][1] = tmp1 - tmp2;
    //
    tmp1 = a1.x * a1.z * t;
    tmp2 = a1.y * s;
    m[2][0] = tmp1 - tmp2;
    m[0][2] = tmp1 + tmp2;
    tmp1 = a1.y * a1.z * t;
    tmp2 = a1.x * s;
    m[2][1] = tmp1 + tmp2;
    m[1][2] = tmp1 - tmp2;
}
void screenShot()
{
    //A capture of the OpenGL window client area.
    static int inc = 0;
    //Get pixels...
    unsigned char *data = new unsigned char [g.xres * g.yres * 3];
    glReadPixels(0, 0, g.xres, g.yres, GL_RGB, GL_UNSIGNED_BYTE, data);
    //Write a PPM file...
    char ts[256], tj[256];
    sprintf(ts, "./images/img%03i.ppm", inc);
    sprintf(tj, "./images/img%03i.jpg", inc);
    ++inc;
    FILE *fpo = fopen(ts, "w");
    fprintf(fpo, "P6\n");
    fprintf(fpo, "%i %i\n", g.xres, g.yres);
    fprintf(fpo, "255\n");
    //Image is upside-down.
    //Go backwards a row at a time...
    unsigned char *p = data;
    p = p + ((g.yres-1) * g.xres * 3);
    unsigned char *start = p;
    for (int i=0; i<g.yres; i++) {
        for (int j=0; j<g.xres*3; j++) {
            fprintf(fpo, "%c", *p);
            ++p;
        }
        start = start - (g.xres*3);
        p = start;
    }
    fclose(fpo);
    delete [] data;
    char t2[2560];
    sprintf(t2, "convert %s %s", ts, tj);
    system(t2);
    unlink(ts);
}
void trans_vector(Matrix mat, const Vec in, Vec out)
{
	Flt f0 = mat[0][0] * in[0] + mat[1][0] * in[1] + mat[2][0] * in[2];
	Flt f1 = mat[0][1] * in[0] + mat[1][1] * in[1] + mat[2][1] * in[2];
	Flt f2 = mat[0][2] * in[0] + mat[1][2] * in[1] + mat[2][2] * in[2];
	out[0] = f0;
	out[1] = f1;
	out[2] = f2;
}

void cube(float w1, float h1, float d1)
{
	float w = w1 * 0.5;
	float d = d1 * 0.5;
	float h = h1 * 0.5;
	//notice the normals being set
	glBegin(GL_QUADS);
		// top
		glNormal3f( 0.0f, 1.0f, 0.0f);
		glVertex3f( w, h,-d);
		glVertex3f(-w, h,-d);
		glVertex3f(-w, h, d);
		glVertex3f( w, h, d);
		// bottom
		glNormal3f( 0.0f, -1.0f, 0.0f);
		glVertex3f( w,-h, d);
		glVertex3f(-w,-h, d);
		glVertex3f(-w,-h,-d);
		glVertex3f( w,-h,-d);
		// front
		glNormal3f( 0.0f, 0.0f, 1.0f);
		glVertex3f( w, h, d);
		glVertex3f(-w, h, d);
		glVertex3f(-w,-h, d);
		glVertex3f( w,-h, d);
		// back
		glNormal3f( 0.0f, 0.0f, -1.0f);
		glVertex3f( w,-h,-d);
		glVertex3f(-w,-h,-d);
		glVertex3f(-w, h,-d);
		glVertex3f( w, h,-d);
		// left side
		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-w, h, d);
		glVertex3f(-w, h,-d);
		glVertex3f(-w,-h,-d);
		glVertex3f(-w,-h, d);
		// right side
		glNormal3f( 1.0f, 0.0f, 0.0f);
		glVertex3f( w, h,-d);
		glVertex3f( w, h, d);
		glVertex3f( w,-h, d);
		glVertex3f( w,-h,-d);
		glEnd();
	glEnd();
}

void tube(int n, float rad, float len)
{
	//Tube is centered at the origin.
	//Base of tube is at { 0, 0, 0 }.
	//Top of tube is at { 0, len, 0 }.
	//
	const int MAX_POINTS = 100;
	static float pts[MAX_POINTS][3];
	static int firsttime = 1;
	if (firsttime) {
		firsttime = 0;
		double angle = 0.0;
		double inc = (3.1415926535 * 2.0) / (double)n;
		for (int i=0; i<n; i++) {
			pts[i][0] = cos(angle) * rad;
			pts[i][2] = sin(angle) * rad;
			pts[i][1] = 0.0f;
			angle += inc;
		}
	}
	glBegin(GL_QUADS);
		for (int i=0; i<n; i++) {
			int j = (i+1) % n;
			glNormal3f(pts[i][0], 0.0, pts[i][2]);
			glVertex3f(pts[i][0], 0.0, pts[i][2]);
			glVertex3f(pts[i][0], len, pts[i][2]);
			glNormal3f(pts[j][0], 0.0, pts[j][2]);
			glVertex3f(pts[j][0], len, pts[j][2]);
			glVertex3f(pts[j][0], 0.0, pts[j][2]);
		}
	glEnd();
}


void drawGround()
{
	int n = 11;
	float w = 10.0;
	float d = 10.0;
	float w2 = w*.49;
	float d2 = d*.49;
	float x = -(n/2)*w;
	float y =    0.0;
	float z = -(n/2)*d;
	float xstart = x;
	static int firsttime = 1;
	if (firsttime) {
		firsttime = 0;
		printf("x: %f\n", x);
	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			glPushMatrix();
			glTranslatef(x, y, z);
			glColor3ub(j*20, 100, 120-i*5);
			glBegin(GL_QUADS);
				glNormal3f( 0.0f, 1.0f, 0.0f);
				glVertex3f( w2, 0.0, -d2);
				glVertex3f(-w2, 0.0, -d2);
				glVertex3f(-w2, 0.0,  d2);
				glVertex3f( w2, 0.0,  d2);
			glEnd();
			glPopMatrix();
			x += w;
		}
		x = xstart;
		z += w;
	}
}

void make_view_matrix(const Vec p1, const Vec p2, Matrix m)
{
    //Line between p1 and p2 form a LOS Line-of-sight.
    //A rotation matrix is built to transform objects to this LOS.
    //Diana Gruber  http://www.makegames.com/3Drotation/
    m[0][0]=m[1][1]=m[2][2]=1.0f;
    m[0][1]=m[0][2]=m[1][0]=m[1][2]=m[2][0]=m[2][1]=0.0f;
    Vec out = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
    //
    Flt l1, len = out[0]*out[0] + out[1]*out[1] + out[2]*out[2];
    if (len == 0.0f) {
        MakeVector(0.0f,0.0f,1.0f,out);
    } else {
        l1 = 1.0f / sqrtf(len);
        out[0] *= l1;
        out[1] *= l1;
        out[2] *= l1;
    }
    m[2][0] = out[0];
    m[2][1] = out[1];
    m[2][2] = out[2];
    Vec up = { -out[1] * out[0], upv[1] - out[1] * out[1], -out[1] * out[2] };
    //
    len = up[0]*up[0] + up[1]*up[1] + up[2]*up[2];
    if (len == 0.0f) {
        MakeVector(0.0f,0.0f,1.0f,up);
    }
    else {
        l1 = 1.0f / sqrtf(len);
        up[0] *= l1;
        up[1] *= l1;
        up[2] *= l1;
    }
    m[1][0] = up[0];
    m[1][1] = up[1];
    m[1][2] = up[2];
    //make left vector.
    VecCross(up, out, m[0]);
}
void vecNormalize(Vec v)
{
    Flt len = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if (len == 0.0)
        return;
    len = 1.0 / sqrt(len);
    v[0] *= len;
    v[1] *= len;
    v[2] *= len;
}

void vecScale(Vec v, Flt s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
}

void swap(struct Smoke* xp, struct Smoke* yp)
{
    struct Smoke temp;
    temp = *xp;
    *xp = *yp;
    *yp = temp;
}
void vecSub(Vec in1, Vec in2, Vec out){
    out[0] = in2[0] - in1[0];
    out[1] = in2[1] - in1[1];
    out[2] = in2[2] - in1[2];
}
void drawSmoke()
{
    
    for(int i = 0; i<g.nsmokes;i++){
        Flt dx = g.smoke[i].pos[0] - g.cameraPosition[0];
        Flt dy = g.smoke[i].pos[1] - g.cameraPosition[1];
        Flt dz = g.smoke[i].pos[2] - g.cameraPosition[2];
        g.smoke[i].Distance = dx*dx + dy*dy + dz*dz;
        for(int j =0; j<g.nsmokes - i;j++){
        if(g.smoke[j].Distance < g.smoke[j+1].Distance){
            swap(&g.smoke[j], &g.smoke[j+1]);
        }
        }
    }
    
    //
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for (int i=0; i<g.nsmokes; i++) {
    glPushMatrix();
    glTranslatef(g.smoke[i].pos[0],g.smoke[i].pos[1],g.smoke[i].pos[2]);
    
        Vec v;
        vecSub(g.smoke[i].pos, g.cameraPosition, v);
        Vec z = {0.0f, 0.0f, 0.0f};
        make_view_matrix(z, v, g.cameraMat);
        float mat[16];
        mat[0] = g.cameraMat[0][0];
        mat[1] = g.cameraMat[0][1];
        mat[2] = g.cameraMat[0][2];
        mat[4] = g.cameraMat[1][0];
        mat[5] = g.cameraMat[1][1];
        mat[6] = g.cameraMat[1][2];
        mat[8] = g.cameraMat[2][0];
        mat[9] = g.cameraMat[2][1];
        mat[10] = g.cameraMat[2][2];
        mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.0f;
        mat[15] = 1.0f;
        glMultMatrixf(mat);
    
    glColor4ub(25, 25, 25, (unsigned char)g.smoke[i].alpha);
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0.0, 0.0, 1.0);
    for (int j=0; j<g.smoke[i].n; j++) {
        //each vertex of the smoke
        //glVertex3fv(g.smoke[i].vert[j]);
        vecNormalize(g.smoke[i].vert[j]);
        vecScale(g.smoke[i].vert[j], g.smoke[i].radius);
        glVertex3fv(g.smoke[i].vert[j]);
    }
    glEnd();
    glPopMatrix();
    }
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}
void drawCloud()
{
    
    for(int i = 0; i<g.nclouds;i++){
        Flt dx = g.cloud[i].pos[0] - g.cameraPosition[0];
        Flt dy = g.cloud[i].pos[1] - g.cameraPosition[1];
        Flt dz = g.cloud[i].pos[2] - g.cameraPosition[2];
        g.cloud[i].Distance = dx*dx + dy*dy + dz*dz;
        for(int j =0; j<g.nclouds - i;j++){
        if(g.cloud[j].Distance < g.cloud[j+1].Distance){
            swap(&g.cloud[j], &g.cloud[j+1]);
        }
        }
    }
    
    //
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for (int i=0; i<g.nclouds; i++) {
    g.cloud[i].pos[0] += g.offset[0];
    g.cloud[i].pos[1] += g.offset[1];
    glPushMatrix();
    glTranslatef(g.cloud[i].pos[0],g.cloud[i].pos[1],g.cloud[i].pos[2]);
        
        Vec v;
        vecSub(g.cloud[i].pos, g.cameraPosition, v);
        Vec z = {0.0f, 0.0f, 0.0f};
        make_view_matrix(z, v, g.cameraMat);
        float mat[16];
        mat[0] = g.cameraMat[0][0];
        mat[1] = g.cameraMat[0][1];
        mat[2] = g.cameraMat[0][2];
        mat[4] = g.cameraMat[1][0];
        mat[5] = g.cameraMat[1][1];
        mat[6] = g.cameraMat[1][2];
        mat[8] = g.cameraMat[2][0];
        mat[9] = g.cameraMat[2][1];
        mat[10] = g.cameraMat[2][2];
        mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.0f;
        mat[15] = 1.0f;
        glMultMatrixf(mat);
    
    glColor4ub(255, 255, 255, (unsigned char)g.cloud[i].alpha);
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0.0, 0.0, 1.0);
    for (int j=0; j<g.cloud[i].n; j++) {
        //each vertex of the smoke
        //glVertex3fv(g.smoke[i].vert[j]);
        vecNormalize(g.cloud[i].vert[j]);
        vecScale(g.cloud[i].vert[j], g.cloud[i].radius);
        glVertex3fv(g.cloud[i].vert[j]);
    }
    glEnd();
    glPopMatrix();
    }
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void make_a_smoke()
{
    if (g.nsmokes < MAX_SMOKES) {
    Smoke *s = &g.smoke[g.nsmokes];
    s->pos[0] = rnd() * 5.0 - 2.5;
    s->pos[2] = rnd() * 5.0 - 2.5;
    s->pos[1] = rnd() * 0.1 + 0.1;
    s->radius = rnd() * 1.0 + 0.5;
    s->n = rand() % 5 + 5;
    Flt angle = 0.0;
    Flt inc = (PI*2.0) / (Flt)s->n;
    for (int i=0; i<s->n; i++) {
        s->vert[i][0] = cos(angle) * s->radius;
        s->vert[i][1] = sin(angle) * s->radius;
        s->vert[i][2] = 0.0;
        angle += inc;
    }
    s->maxtime = 8.0;
    s->alpha = 100.0;
    clock_gettime(CLOCK_REALTIME, &s->tstart);
    ++g.nsmokes;
    }
}
void make_a_cloud()
{
    if (g.nclouds < MAX_SMOKES * 6) {
    Smoke *c = &g.cloud[g.nclouds];
    c->pos[0] = rnd() * 5.0 - 2.5;
    c->pos[2] = rnd() * 5.0 - 2.5;
    c->pos[1] = rnd() * 0.1 + 0.1;
    c->radius = (rnd() * 1.0 + 0.5)/4;
    c->n = rand() % 5 + 5;
    Flt angle = 0.0;
    Flt inc = (PI*2.0) / (Flt)c->n;
    for (int i=0; i<c->n; i++) {
        c->vert[i][0] = cos(angle) * c->radius;
        c->vert[i][1] = sin(angle) * c->radius;
        c->vert[i][2] = 0.0;
        angle += inc;
    }
    c->maxtime = 8.0;
    c->alpha = 100.0;
    clock_gettime(CLOCK_REALTIME, &c->tstart);
    ++g.nclouds;
    }

}

void make_a_plane()
{   
    float w = 0.5 * 0.5;
    float d = 0.5 * 0.5;
    float h = 0.5 * 0.5;
	//notice the normals being set
	glBegin(GL_QUADS);
		// top
		glNormal3f( 0.0f, 1.0f, 0.0f);
		glVertex3f( w, h,-d);
		glVertex3f(-w, h,-d);
		glVertex3f(-w, h, d);
		glVertex3f( w, h, d);
		// bottom
		glNormal3f( 0.0f, -1.0f, 0.0f);
		glVertex3f( w*2,-h, d*2);
		glVertex3f(-w*2,-h, d*2);
		glVertex3f(-w*2,-h,-d*2);
		glVertex3f( w*2,-h,-d*2);
		// front
		glNormal3f( 0.0f, 0.0f, 1.0f);
		glVertex3f( w, h, d);
		glVertex3f(-w, h, d);
		glVertex3f(-w*2,-(((h*30)/2)-1.0), d*2);
		glVertex3f(w*2,-(((h*30)/2)-1.0), d*2);
		// back
		glNormal3f( 0.0f, 0.0f, -1.0f);
		glVertex3f( w*2,-h,-d*2);
		glVertex3f(-w*2,-h,-d*2);
		glVertex3f(-w, h,-d);
		glVertex3f( w, h,-d);
		// left side
		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-w, h, d);
		glVertex3f(-w, h,-d);
		glVertex3f(-w*2,-h,-d*2);
		glVertex3f(-w*2,-h, d*2);
		// right side
		glNormal3f( 1.0f, 0.0f, 0.0f);
		glVertex3f( w, h,-d);
		glVertex3f( w, h, d);
		glVertex3f( w*2,-h, d*2);
		glVertex3f( w*2,-h,-d*2);
        // front plane body top
       /* glNormal3f( 0.0f, 0.0f, 1.0f);
        glVertex3f( -w*2,-h, d*2);
        glVertex3f( w*2,-h, d*2);
        glVertex3f( w*2,-(((h*30)/2)-1.0), d*2 );
        glVertex3f( -w*2,-(((h*30)/2)-1.0), d*2 );*/
        // plane cockpit bottom
        glColor3ub(255,255,255);
        glNormal3f( 0.0f, 0.0f, 1.0f);
        glVertex3f( -w*2,-(((h*30)/2)-1.0), d*2 );
        glVertex3f( w*2,-(((h*30)/2)-1.0), d*2 );
        glVertex3f( w*2,-(((h*30)/2)+1.0), d*2 );
        glVertex3f( -w*2,-(((h*30)/2)+1.0), d*2 );
		glEnd();
	glEnd();
}
void make_a_plane2()
{  
    int n = 10;
    float rad = 0.5;
    float len = 0.5 * 0.5 * 30;
   const int MAX_POINTS = 100;
        static float pts[MAX_POINTS][3];
        static int firsttime = 1;
        if (firsttime) {
                firsttime = 0;
                double angle = 0.0;
                double inc = (3.1415926535 * 2.0) / (double)n;
                for (int i=0; i<n; i++) {
                        pts[i][0] = cos(angle) * rad;
                        pts[i][2] = sin(angle) * rad;
                        pts[i][1] = 0.0f;
                        angle += inc;
                }
        }
        glBegin(GL_QUADS);
        // Plane body 
                for (int i=0; i<n; i++) {
                        int j = (i+1) % n;
                        glNormal3f(pts[i][0], pts[i][2], 0.0);
                        glVertex3f(pts[i][0], pts[i][2], 0.0);
                        glVertex3f(pts[i][0], pts[i][2], len);
                        glNormal3f(pts[j][0], pts[j][2], 0.0);
                        glVertex3f(pts[j][0], pts[j][2], len);
                        glVertex3f(pts[j][0], pts[j][2], 0.0);
                }
		// Plane body top
		for (int i=0; i<n+1; i++) {
                        int j = (i+1) % n;
                        glNormal3f(0.0, 0.0, 1.0);
                        glVertex3f(0.0, 0.0, len);
			glVertex3f(pts[j][0], pts[j][2], len);
			glVertex3f(pts[i][0], pts[i][2], len);
			glVertex3f(0.0, 0.0, len);
                }
        // Tail Body
		for (int i=0; i<n+1; i++) {
                        int j = (i+1) % n;
                        glNormal3f(pts[i][0], pts[i][2], pts[i][2]);
                        glVertex3f(0.0, 0.5, -3.0);
                        glVertex3f(pts[j][0], pts[j][2], 0.0);
                        glVertex3f(pts[i][0], pts[i][2], 0.0);
                        glVertex3f(0.0, 0.5, -3.0);
                }
            // Wings
    
        glNormal3f(1.0, 0.0, 0.0);
        // Right wing top
        glColor3f(0.5, 0.5, 0.5); // Gray color for the wings
        glVertex3f(0.0, 0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, 0.25, len/2-1.0); // Wing root
        glVertex3f(5.0, 0.0, len/2-2.5); // Wing edge
        glVertex3f(5.0, 0.0, len/2-1.0);
        glNormal3f(-1.0, 0.0, 0.0);
        // Right wing bottom
        glColor3f(0.5, 0.5, 0.5); // Gray color for the wings
        glVertex3f(0.0, -0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, -0.25, len/2-1.0); // Wing root
        glVertex3f(5.0, 0.0, len/2-2.5); // Wing edge
        glVertex3f(5.0, 0.0, len/2-1.0);

        // Right wing right side bottom-top
        glNormal3f(0.5, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, 0.0, len/2-1.5); // Wing root
        glVertex3f(0.0, 0.25, len/2-1.0); // Wing edge
        glVertex3f(5.0, 0.0, len/2 - 2.5);
        glVertex3f(5.0, 0.0, len/2-2.5);
        // Right wing right side bottom-bottom
        glNormal3f(-0.5, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, 0.0, len/2-1.5); // Wing root
        glVertex3f(0.0, -0.25, len/2-1.0); // Wing edge
        glVertex3f(5.0, 0.0, len/2 - 2.5);
        glVertex3f(5.0, 0.0, len/2-2.5);
        // Right wing right side top-top
        glNormal3f(0.0, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, -0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, 0.25, len/2+2.5); // Wing edge
        glVertex3f(5.0, 0.0, len/2-1.0);
        glVertex3f(5.0, 0.0, len/2-1.0);
        glNormal3f(1.0, 0.0, 0.0);
        // Left wing top
        glColor3f(0.5, 0.5, 0.5); // Gray color for the wings
        glVertex3f(0.0, 0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, 0.25, len/2-1.0); // Wing root
        glVertex3f(-5.0, 0.0, len/2-2.5); // Wing edge
        glVertex3f(-5.0, 0.0, len/2-1.0);
        glNormal3f(-1.0, 0.0, 0.0);
        // Left wing bottom
        glColor3f(0.5, 0.5, 0.5); // Gray color for the wings
        glVertex3f(0.0, -0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, -0.25, len/2-1.0); // Wing root
        glVertex3f(-5.0, 0.0, len/2-2.5); // Wing edge
        glVertex3f(-5.0, 0.0, len/2-1.0);
        // Left wing right side bottom-top
        glNormal3f(0.5, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, 0.0, len/2-1.5); // Wing root
        glVertex3f(0.0, 0.25, len/2-1.0); // Wing edge
        glVertex3f(-5.0, 0.0, len/2 - 2.5);
        glVertex3f(-5.0, 0.0, len/2-2.5);
        // Left wing right side bottom-bottom
        glNormal3f(-0.5, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, 0.0, len/2-1.5); // Wing root
        glVertex3f(0.0, -0.25, len/2-1.0); // Wing edge
        glVertex3f(-5.0, 0.0, len/2 - 2.5);
        glVertex3f(-5.0, 0.0, len/2-2.5);
        // Left wing right side top-top
        glNormal3f(0.0, 1.0, 0.0);
        glColor3f(1.0, 0.0, 0.0); // red color for the wings
        glVertex3f(0.0, -0.25, len/2+2.5); // Wing root
        glVertex3f(0.0, 0.25, len/2+2.5); // Wing edge
        glVertex3f(-5.0, 0.0, len/2-1.0);
        glVertex3f(-5.0, 0.0, len/2-1.0);

        glColor3f(0.25, 0.25, 0.25);
        glNormal3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 0.25, 0.0);
        glVertex3f(1.5, 0.25, -1.5);
        glVertex3f(1.5, 0.25, -1.5);
        glVertex3f(0.0, 0.25, -1.5);
        glNormal3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 0.25, 0.0);
        glVertex3f(-1.5, 0.25, -1.5);
        glVertex3f(-1.5, 0.25, -1.5);
        glVertex3f(0.0, 0.25, -1.5);
        glNormal3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.25, 0.0);
        glVertex3f(0.0, 1.5, -1.5);
        glVertex3f(0.0, 1.5, -1.5);
        glVertex3f(0.0, 0.25, -1.5);
        
        glNormal3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.25, 0.0);
        glVertex3f(1.5, 0.25, -1.5);
        glVertex3f(0.0, 0.30, 0.0);
        glVertex3f(1.5, 0.25, -1.5);




    //glVertex3f(0.0, len / 2, 0.0); // Wing root
    //glVertex3f(1.0, len / 2, 0.0); // Wing tip
    //glVertex3f(0.0, len / 2, 5.0); // Wing edge
        
    //glColor3f(0.5, 0.5, 0.5); // Gray color for the tail
    //glVertex3f(0.0, len, 0.0); // Tail root
    //glVertex3f(0.0, len + 3.0, 0.0); // Tail tip
    //glVertex3f(0.0, len, 5.0); // Tail edge
    glEnd();
}

void Global::physics()
{
//	g.cameraPosition[2] -= 0.1;
//	g.cameraPosition[0] = 1.0 + sin(g.cameraPosition[2]*0.3);
    if(g.gamestart > 0){
        clock_gettime(CLOCK_REALTIME, &g.smokeTime);
        double d = timeDiff(&g.smokeStart, &g.smokeTime);
        if (d > 0.01) {
        //time to make another smoke particle
        make_a_smoke();
        timeCopy(&g.smokeStart, &g.smokeTime);
        }
        //move smoke particles
        for (int i=0; i<g.nsmokes; i++) {
        //smoke rising
        g.smoke[i].pos[1] += 0.015;
        g.smoke[i].pos[1] += ((g.smoke[i].pos[1]*0.24) * (rnd() * 0.075));
        //expand particle as it rises
        g.smoke[i].radius += g.smoke[i].pos[1]*0.002;
        //wind might blow particle
        if (g.smoke[i].pos[1] > 10.0) {
            //experiment here with different values
            g.smoke[i].pos[0] -= rnd() * 0.1;
        }
        }
        //this is where a smoke particle will fade away as it lingers
        int i=0;
        while (i < g.nsmokes) {
        struct timespec bt;
        clock_gettime(CLOCK_REALTIME, &bt);
        double d = timeDiff(&g.smoke[i].tstart, &bt);
        if (d > g.smoke[i].maxtime - 3.0) {
            g.smoke[i].alpha *= 0.95;
            if (g.smoke[i].alpha < 1.0)
            g.smoke[i].alpha = 1.0;
        }
        if (d * 3 > (g.smoke[i].maxtime)) {
            //delete this smoke
            --g.nsmokes;
            g.smoke[i] = g.smoke[g.nsmokes];
            continue;
        }
        ++i;
        }
    } else {
        clock_gettime(CLOCK_REALTIME, &g.cloudTime);
        double d = timeDiff(&g.cloudStart, &g.cloudTime);
        if (d > 0.01) {
        //time to make another smoke particle
        make_a_cloud();
        timeCopy(&g.cloudStart, &g.cloudTime);
        }
        //move smoke particles
        for (int i=0; i<g.nclouds; i++) {
        //smoke rising
        g.cloud[i].pos[0] += 0.015;
        g.cloud[i].pos[0] += ((g.cloud[i].pos[0]*0.24) * (rnd() * 0.075));
        //expand particle as it rises
        //g.cloud[i].radius += g.cloud[i].pos[0]*0.002;
        //wind might blow particle
        g.cloud[i].pos[1] += 0.015;
        if (g.cloud[i].pos[0] > 10.0) {
            //experiment here with different values
            //g.smoke[i].pos[1] += rnd() * 0.01;
        }
        }
        //this is where a smoke particle will fade away as it lingers
     /*   int i=0;
        while (i < g.nclouds) {
        struct timespec bt;
        clock_gettime(CLOCK_REALTIME, &bt);
        double d = timeDiff(&g.cloud[i].tstart, &bt);
        if (d > g.cloud[i].maxtime - 3.0) {
            g.cloud[i].alpha *= 0.95;
            if (g.cloud[i].alpha < 1.0)
            g.cloud[i].alpha = 1.0;
        }
        if (d > g.cloud[i].maxtime) {
            //delete this smoke
            --g.nclouds;
            g.cloud[i] = g.cloud[g.nclouds];
            continue;
        }
        ++i;
        
    }*/

}
 
}

void Global::render()
{
    if(g.gamestart > 0){
        static int justonce = 1;
        if(justonce == 1){
            g.cameraPosition[2] = g.initialz;
            justonce = 0;
        }
        Rect r;
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        if (pause) {
            // Setup for 2D rendering
            glViewport(0, 0, xres, yres);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluOrtho2D(0, xres, 0, yres);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();

            // Disable lighting for clear text rendering
            glDisable(GL_LIGHTING);
            glDisable(GL_DEPTH_TEST);

            Rect r;
            r.bot = yres / 2 + 50;
            r.left = xres / 2;
            r.center = 1;
            ggprint12(&r, 16, 0x00ff0000, "Game Paused");

            // Options for the pause menu
            const char *options[] = { "Resume", "Options", "Main Menu" };
            for (int i = 0; i < 3; i++) {
                r.bot -= 20;
                unsigned int color;
                if (optionSelected == i) {
                    color = 0x00ff0000; // selected: red
                } else {
                    color = 0xffffffff; // nonselected: white
                }
                ggprint8b(&r, 16, color, options[i]);
            }

            // Re-enable anything disabled before returning
            glEnable(GL_LIGHTING);
            glEnable(GL_DEPTH_TEST);
            return;
        }
        //
        //3D mode
        //
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0f, g.aspectRatio, 0.5f, 1000.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        //for documentation...
        Vec up = { 0.0, 1.0, 0.0 };
        float horDist = g.cameraDistancePlayer * cos((pitch * M_PI)/180);
        float virtDist = g.cameraDistancePlayer * sin((pitch * M_PI)/180);
        float theta = planeAngle[0] + g.angleAroundPlayer;;
        float offsetX = horDist * sin((theta*M_PI)/180);
        float offsetZ = horDist * cos((theta*M_PI)/180);
        g.cameraPosition[1] = g.plane2Pos[1] + virtDist;
        g.cameraPosition[0] = g.plane2Pos[0] - offsetX;
        g.cameraPosition[2] = g.plane2Pos[2] - offsetZ;
        gluLookAt(
            g.cameraPosition[0], g.cameraPosition[1], g.cameraPosition[2],
            g.plane2Pos[0], 
            g.plane2Pos[1], 
            g.plane2Pos[2]-1.0,
            up[0], up[1], up[2]);
        glLightfv(GL_LIGHT0, GL_POSITION, g.lightPosition);
        //
        drawGround();
        drawSmoke();
        //
        //Draw tube
        glPushMatrix();
        glRotatef(45.0, 0.0, 1.0, 0.0);
        glTranslatef(-8.0, 2.0, -8.0);
        glRotatef(-90.0, 1.0, 0.0, 0.0);
        tube(10, 2.0, 20.0);
        glPopMatrix();
        //
        static float angle = 0.0;
        glPushMatrix();
        glTranslatef(-8.0, 2.0, -8.0);
        glRotatef(angle, 0.0, 1.0, 0.0);
        angle += 1.1;
        cube(0.5, 0.5, 4.0);
        glPopMatrix();

        //glPushMatrix();
       // glTranslatef(g.cameraDir[0], g.cameraDir[1], g.cameraDir[2]-5);
       // make_a_plane();
       // glPopMatrix();
        glPushMatrix();
        glTranslatef(g.plane2Pos[0], g.plane2Pos[1], g.plane2Pos[2]);
	//glRotatef(90.0, 0.0, 0.0, 1.0);
    //glRotatef(-90.0, 1.0, 0.0, 0.0);
    //glRotatef(90.0, 0.0, 0.0, 1.0);
    //glRotatef(90.0, 1.0, 0.0, 0.0);
    glRotatef(g.planeAngle[1], 0.0, 0.0, 1.0);
    glRotatef(g.planeAngle[0], 0.0, 1.0, 0.0);
    glRotatef(g.planeAngle[2], 1.0, 0.0, 0.0);
    make_a_plane2();
        glPopMatrix();
        
        //
        //
        //switch to 2D mode
        //
        glViewport(0, 0, g.xres, g.yres);
        glMatrixMode(GL_MODELVIEW);   glLoadIdentity();
        glMatrixMode (GL_PROJECTION); glLoadIdentity();
        gluOrtho2D(0, g.xres, 0, g.yres);
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        r.bot = g.yres - 20;
        r.left = 10;
        r.center = 0;

        ggprint8b(&r, 16, 0x00887766, "fps framework");
        ggprint8b(&r, 16, 0x00990000, "For help press 'l'");
        if(g.menu){
            ggprint8b(&r, 16, 0x00990000, "w: Move Foward");
            ggprint8b(&r, 16, 0x00990000, "s: Move Back");
            ggprint8b(&r, 16, 0x00990000, "a: Move Left");
            ggprint8b(&r, 16, 0x00990000, "s: Move Right");
            ggprint8b(&r, 16, 0x00990000, "Up: look down");
            ggprint8b(&r, 16, 0x00990000, "Down: look up");
            ggprint8b(&r, 16, 0x00990000, "left: look Left");
            ggprint8b(&r, 16, 0x00990000, "right: look Right");
            ggprint8b(&r, 16, 0x00990000, "p: take screenshot");

        }
        glPopAttrib();
    } else {
        /*
        double d;
        clock_gettime(CLOCK_REALTIME, &g.introTime);
        d = timeDiff(&introStart, &introTime);
        timeCopy(&introStart, &introTime);
        */
        
        Rect r;
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        //
        //3D mode
        //
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        GLuint texture[256][256][3];
        for(int i = 0; i < 256; i++) {
            for(int j = 0; j < 256; j++) {
                float color = background -> cloudmap256[i*256+j];
                texture[i][j][0] = color;
                texture[i][j][1] = color;
                texture[i][j][2] = color;
            }
        }

        gluPerspective(45.0f, g.aspectRatio, 0.5f, 1000.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        //for documentation...
        Vec up = { 0.0, 1.0, 0.0 };
        gluLookAt(
            g.cameraPosition[0], g.cameraPosition[1], g.cameraPosition[2],
            g.cameraPosition[0] + g.cameraDir[0],
            g.cameraPosition[1] + g.cameraDir[1],
            g.cameraPosition[2] + g.cameraDir[2],
            up[0], up[1], up[2]);
        glLightfv(GL_LIGHT0, GL_POSITION, g.lightPosition);
        drawCloud();
        glViewport(0, 0, g.xres, g.yres);
        glMatrixMode(GL_MODELVIEW);   glLoadIdentity();
        glMatrixMode (GL_PROJECTION); glLoadIdentity();
        gluOrtho2D(0, g.xres, 0, g.yres);
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        r.bot = g.yres - 20;
        r.left = 10;
        r.center = 0;
        //ggprint8b(&r, 16, 0x00990000, "to start press 'm'");
        glPopAttrib();
	char buff[50];
        sprintf(buff, "fps: %d", g.fps);
        ggprint8b(&r, 16, 0xffffffff, buff);
	sprintf(buff, "vsync: %d", g.vsync);
        ggprint8b(&r, 16, 0xffffffff, buff);

        if(g.cameraPosition[2] > (5.0-2.5))
            g.cameraPosition[2] -= 0.05;
        if(g.cameraPosition[2] <= (5.0-2.5)){
            r.bot = g.yres/2;
            r.left = g.xres/2 - g.xres/10;
            r.center = 0;
            ggprint16(&r, 16, 0x00990000, "Aerial Obstacle");
            r.left = g.xres/2 - g.xres/10 + 10;
            ggprint12(&r, 16, 0x00990000, "to start game");
            r.left = g.xres/2 - g.xres/10 + 20;
            ggprint12(&r, 16, 0x00990000, "press 'm'");
        }

    }
}



