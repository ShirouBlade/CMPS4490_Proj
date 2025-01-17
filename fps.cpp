//
//program: fps.cpp
//author:  Gordon Griesel
//date:    Winter 2020
//
//Modified by: Danny Simpson
//Date: January 30, 2024
//
//Framework for a 3D game.
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
#include "defs.h"
#include "log.h"
#include "fonts.h"

typedef float Flt;
typedef Flt Vec[3];
typedef Flt Matrix[4][4];

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
	Flt normalize() {
		return 0.0f;
	}
};

class Global {
public:
	int xres, yres;
	Flt aspectRatio;
	Vec cameraPosition;
    Vec cameraDir;
	GLfloat lightPosition[4];
	Global() {
		//constructor
		xres=640;
		yres=480;
		aspectRatio = (GLfloat)xres / (GLfloat)yres;
		MakeVector(0.0, 1.0, 8.0, cameraPosition);
        MakeVector(0.0, 0.0, -1.0, cameraDir);
		//light is up high, right a little, toward a little
		MakeVector(100.0f, 240.0f, 40.0f, lightPosition);
		lightPosition[3] = 1.0f;
		//init_opengl();
	}
	void init_opengl();
	void init();
	void check_mouse(XEvent *e);
	int check_keys(XEvent *e);
	void physics();
	void render();
} g;


class X11_wrapper {
private:
	Display *dpy;
	Window win;
	GLXContext glc;
public:
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


int main()
{
	g.init_opengl();
	int done = 0;
	while (!done) {
		while (x11.getXPending()) {
			XEvent e = x11.getXNextEvent();
			x11.check_resize(&e);
			g.check_mouse(&e);
			done = g.check_keys(&e);
		}
		g.physics();
		g.render();
		x11.swapBuffers();
	}
	cleanup_fonts();
	return 0;
}

void Global::init() {
	//Place general program initializations here.
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

int Global::check_keys(XEvent *e)
{
	//Was there input from the keyboard?
	if (e->type == KeyPress) {
		int key = (XLookupKeysym(&e->xkey, 0) & 0x0000ffff);
        float speed = 0.1;
        //int qdown = 0;
		switch(key) {
			case XK_1:
				break;
			case XK_Right:
				//g.cameraPosition[0] += 1.0;
                {    
                    Vec v = {0.0, -0.1, 0.0};
                    Matrix m;
                    identity33(m);
                    yy_transform(v, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
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
                }
				break;
			case XK_Up:
				//g.cameraPosition[1] += 0.2;
                //tilt camera down
                //need to rotate camera vector
                {
                    Vec v = {0.1, 0.0, 0.0};
                    Matrix m;
                    identity33(m);
                    yy_transform(v, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                }
				break;
			case XK_Down:
				//g.cameraPosition[1] -= 0.2;
                {
                    Vec v = {-0.1, 0.0, 0.0};
                    Matrix m;
                    identity33(m);
                    yy_transform(v, m);
                    trans_vector(m, g.cameraDir, g.cameraDir);
                }
				break;
			case XK_f:
				//g.cameraPosition[2] -= 1.0;
                // add the directio vector to the camera position vector
                for(int i=0; i < 3; i++){
                    g.cameraPosition[i] += g.cameraDir[i] * speed;
                }
				break;
			case XK_b:
				//g.cameraPosition[2] += 1.0;
                for(int i=0; i < 3; i++){
                    g.cameraPosition[i] -= g.cameraDir[i] * speed;
                }
				break;
            case XK_w:
                g.cameraPosition[1] += 0.2;
                break;
            case XK_s:
                g.cameraPosition[1] += 0.2;
                break;
            case XK_a:
                g.cameraPosition[0] -= 1.0;
                break;
            case XK_d:
                g.cameraPosition[0] += 1.0;
                break;
            case XK_q:
                //qdown ^= qdown;
                //if(qdown != 0){
                //    speed = 10.0;
                //} else {
                //    speed = 0.1;
                //}
                break;
			case XK_Escape:
				return 1;
		}
	}
	return 0;
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

void Global::physics()
{
//	g.cameraPosition[2] -= 0.1;
//	g.cameraPosition[0] = 1.0 + sin(g.cameraPosition[2]*0.3);
}

void Global::render()
{
	Rect r;
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
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
	gluLookAt(
		g.cameraPosition[0], g.cameraPosition[1], g.cameraPosition[2],
		g.cameraPosition[0] + g.cameraDir[0], 
        g.cameraPosition[1] + g.cameraDir[1], 
        g.cameraPosition[2] + g.cameraDir[2],
		up[0], up[1], up[2]);
	glLightfv(GL_LIGHT0, GL_POSITION, g.lightPosition);
	//
	drawGround();
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
	glPopAttrib();
}



