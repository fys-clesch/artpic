#include "artpic.h"
#include "msg.h"
#include "alloc.h"
#include "color.h"
#include "control.h"
#include "shapes.h"
#include "draw.h"
#include "viewer.h"

/**
 * Good OpenGL tutorial: http://www.3dsource.de/faq/viewing.htm.
 */
uint MAIN_WIDTH,
     MAIN_HEIGHT,
     PARENT_WIN,
     SUB_WIN_1;
const uint INIT_MAIN_WIDTH = 800,
           INIT_MAIN_HEIGHT = 600;
const double NEARCLIPPING = -1000.,
             FARCLIPPING = 1000.;
double ROT_Y = 0.,
       NOR_Y = 0.,
       ROT_X = 0.,
       MOV_Y = 0.,
       MOV_X = 0.,
       ZOOM = DEFAULT_ZOOM;
uchar USE_LIGHT = 0;

extern const uchar MOUSE_L_DOWN,
                   MOUSE_R_DOWN,
                   STOP_ROT,
                   TRANSLATE,
                   ESC_CLICKED;
extern const uint DRAW_RAY_N;
extern const double BIN_SPHERE3_ALPHA;

/** \brief Initiates the OpenGL viewer.
 *
 * \param argc int The number of arguments taken from the main function.
 * \param argv char** The list of arguments taken from the main function.
 * \return void
 *
 */
void go_freeglut(int argc, char **argv)
{
    PARENT_WIN = init_main(argc, argv);
    glutDisplayFunc(&maindisplay);
    glutReshapeFunc(&reshape);
    glutCloseFunc(&close_all);
    glutMotionFunc(&motion);
    glutPassiveMotionFunc(&passivemotion);
    glutMouseFunc(&mouse);
    glutKeyboardFunc(&keyboard);
    glutSpecialFunc(&arrow_keys);
    glutVisibilityFunc(&visibility);
    glutWindowStatusFunc(&win_state);
    SUB_WIN_1 = glutCreateSubWindow(PARENT_WIN, 5, 5, MAIN_WIDTH / 4, MAIN_HEIGHT / 5);
    glutDisplayFunc(&subdisplay);
    glutKeyboardFunc(&keyboard);
    glutSpecialFunc(&arrow_keys);
    glutMainLoop();
}

/** \brief Initializes the OpenGL mode
 *
 * \param argc int The first argument from main.
 * \param argv char** The second argument from main.
 * \return uint The identifier of the main window.
 *
 */
uint init_main(int argc, char **argv)
{
    fncyprint("\ninitializing graphics...", 50);
    MAIN_WIDTH = INIT_MAIN_WIDTH;
    MAIN_HEIGHT = INIT_MAIN_HEIGHT;
    glutInit(&argc, argv);
    glutInitWindowSize(MAIN_WIDTH, MAIN_HEIGHT);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutCreateWindow("artpic");
    glutSetCursor(GLUT_CURSOR_CROSSHAIR);
    /* glEnable(GL_DEPTH_TEST); */ /**< If enabled, one can not see through a transparent sphere. */
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHT0);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    /* glEnable(GL_NORMALIZE); */
    glEnable(GL_COLOR_MATERIAL);
    glLightfv(GL_LIGHT0, GL_AMBIENT, LIGHT_AMBIENT);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, LIGHT_DIFFUSE);
    glLightfv(GL_LIGHT0, GL_SPECULAR, LIGHT_SPECULAR);
    glLightfv(GL_LIGHT0, GL_POSITION, LIGHT_POSITION);
    glMaterialfv(GL_FRONT, GL_AMBIENT, MAT_AMBIENT);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MAT_DIFFUSE);
    glMaterialfv(GL_FRONT, GL_SPECULAR, MAT_SPECULAR);
    glMaterialfv(GL_FRONT, GL_SHININESS, HIGH_SHININESS);
    glClearColor(0., 0., 0., 0.);
    glClearDepth(1.);
    glViewport(0, 0, MAIN_WIDTH, MAIN_HEIGHT);
    if(glutGet(GLUT_DISPLAY_MODE_POSSIBLE))
        fprintf(stdout, " done\n");
    else
        fprintf(stdout, " ! can not establish the correct display mode. flawless display is unlikely.\n");
    fprintf(stdout,
            "rendering with\n"
            "OGL version    %s\nrenderer       %s\nvendor         %s\n",
            glGetString(GL_VERSION), glGetString(GL_RENDERER), glGetString(GL_VENDOR));
    return glutGetWindow();
}

/** \brief Handles the updating functionality.
 *
 * \param void
 * \return void
 *
 * This is called by glutIdleFunc.
 */
void animate(void)
{
    static uint avg = 50,
                avgtime = 0,
                frams = 0;
    if(frams++ < avg)
        avgtime = glutGet(GLUT_ELAPSED_TIME);
    else
    {
        static char title[64];
        static uint offset = 0;
        frams = 0;
        snprintf(title, 64, "artpic - %.1f fps", 1e3 * avg / (avgtime - offset));
        glutSetWindowTitle(title);
        avgtime = offset = glutGet(GLUT_ELAPSED_TIME);
    }
    glutSetWindow(PARENT_WIN);
    glutPostRedisplay();
    glutSetWindow(SUB_WIN_1);
    glutPostRedisplay();
    usleep(10000);
    if(!MOUSE_L_DOWN && !STOP_ROT)
    {
        static const double vrot = .02;
        ROT_Y += vrot * (glutGet(GLUT_ELAPSED_TIME) - avgtime); /**< Constant rotation. */
    }
    ROT_Y = fmod(ROT_Y, 360.);
    ROT_X = fmod(ROT_X, 360.);
}

/** \brief Changes the drawing mode in case the window is not visible.
 *
 * \param vis const int The visibility status.
 * \return void
 *
 * This is called by glutVisibilityFunc.
 */
void visibility(const int vis)
{
    if(vis == GLUT_VISIBLE) glutIdleFunc(animate);
    else
    {
        glutSetWindowTitle("artpic - idle");
        glutIdleFunc(NULL);
    }
}

/** \brief Changes the drawing mode in case the window is covered or inactive.
 *
 * \param state const int The visibility status.
 * \return void
 *
 * This is called by glutWindowStatusFunc.
 */
void win_state(const int state)
{
    if(state == GLUT_FULLY_RETAINED) glutIdleFunc(animate);
    else
    {
        glutSetWindowTitle("artpic - idle");
        glutIdleFunc(NULL);
    }
}

/** \brief Handles the main window of the application.
 *
 * \param void
 * \return void
 *
 * This is called by the first instance of glutDisplayFunc.
 */
void maindisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScaled(ZOOM, ZOOM, ZOOM);
    glTranslated(MOV_X, MOV_Y, 0.);
    glRotated(ROT_X, 1., 0., 0.);
    glRotated(ROT_Y, NOR_Y, 1., 0.);
    //handle_glray((const glray *)NULL,0,DRAW_EM_ALL,1);
    handle_glray((const glray *)NULL, 0, DRAW_SOME, 1, DRAW_RAY_N);
    if(USE_LIGHT)
        glEnable(GL_LIGHTING);
    //handle_prtcls_boxed((const sphrcl_prtcl const *)NULL,(const boundingbox *)NULL,0,DRAW_EM_ALL);
    handle_prtcls((const sphrcl_prtcl *const)NULL, DRAW_EM_ALL);
    handle_bin_sphere3((const bin_hit_screen *const)NULL, DRAW_EM_ALL, SINCOS_MAP, INTENSITY);
    if(USE_LIGHT)
        glDisable(GL_LIGHTING);
    float width;
    glGetFloatv(GL_LINE_WIDTH, &width);
    glLineWidth(2.);
    const double xs[3] = { -10., 0., 0.},
                 xf[3] = {10., 0., 0.},
                 ys[3] = {0., -10., 0.},
                 yf[3] = {0., 10., 0.},
                 zs[3] = {0., 0., -10.},
                 zf[3] = {0., 0., 10.};
    glColor4dv(CWHITE);
    draw_arrowv(xs, xf, .15, .1, 6);
    draw_arrowv(ys, yf, .15, .1, 6);
    draw_arrowv(zs, zf, .15, .1, 6);
    draw_coord_ov();
    glLineWidth(width);
    draw_cb(0);
    glutSwapBuffers();
}

/** \brief A subdisplay handler.
 *
 * \param void
 * \return void
 *
 * This can be called after an instance of glutCreateSubWindow.
 */
void subdisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0., 4., 0., 4.);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glColor3dv(CWHITE);
    /*drawblockstring2d("these keys may have the one or the other effect: "
                      "F1, F11, F12, Esc, R, m, Mouse Buttons. "
                      "this list might be incomplete.",
                      23, .2, 3.6, 0., -.4);*/
    char info[256];
    snprintf(info, 256, "printing ray %u", DRAW_RAY_N);
    drawblockstring2d(info, 23, .2, 3.6, 0., -.4);
    snprintf(info, 256, "screen alpha %.2f", BIN_SPHERE3_ALPHA);
    drawblockstring2d(info, 23, .2, 3.2, 0., -.4);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glutSwapBuffers();
}

/** \brief Takes care of proper resizing the graphics after the window size is changed.
 *
 * \param width const int The new width of the window.
 * \param height const int The new height of the window.
 * \return void
 *
 */
void reshape(const int width, const int height)
{
    if((uint)width != MAIN_WIDTH || (uint)height != MAIN_HEIGHT)
    {
        MAIN_WIDTH = width;
        MAIN_HEIGHT = height;
        glViewport(0, 0, MAIN_WIDTH, MAIN_HEIGHT);
    }
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-(double)MAIN_WIDTH / 2., (double)MAIN_WIDTH / 2.,
            -(double)MAIN_HEIGHT / 2., (double)MAIN_HEIGHT / 2.,
            NEARCLIPPING, FARCLIPPING);
}

/** \brief Handles deallocation of the memory allocated during the drawing processes.
 *
 * \param void
 * \return void
 *
 */
void free_drawing_mem(void)
{
    fprintf(stdout, "freeing drawing memory... ");
    fflush(stdout);
    handle_glray((const glray *)NULL, 0, FREE_EM_ALL, 0);
    handle_sphere3(0, 0, 0., (const double *)NULL, FREE_EM_ALL);
    handle_bin_sphere3((const bin_hit_screen *const)NULL, FREE_EM_ALL, 0, 0);
    handle_prtcls_boxed((const sphrcl_prtcl *const)NULL, (const boundingbox *)NULL, 0, FREE_EM_ALL);
    draw_cb(1);
    fprintf(stdout, "done\n");
}

/** \brief Handles proper closing of the windows and deallocation of the memory.
 *
 * \param void
 * \return void
 *
 * This is called by glutCloseFunc.
 */
void close_all(void)
{
    static uchar done = 0;
    if(!done) /**< No chance to free the memory twice. */
    {
        if(glutGetWindow() && ESC_CLICKED)
        {
            glutIdleFunc(NULL);
            free_drawing_mem();
            glutDestroyWindow(PARENT_WIN);
            glutLeaveMainLoop();
        }
        else
        {
            free_drawing_mem();
            glutLeaveMainLoop();
        }
        done = 1;
    }
}
