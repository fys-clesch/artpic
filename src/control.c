#include "artpic.h"
#include "viewer.h"
#include "control.h"

uchar MOUSE_L_DOWN = 0,
      MOUSE_R_DOWN = 0,
      STOP_ROT = 0,
      TRANSLATE = 0,
      ESC_CLICKED = 0;
uint DRAW_RAY_N = 0;
double BIN_SPHERE3_ALPHA = DEFAULT_BIN_SPHERE3_ALPHA;

extern const uint MAIN_WIDTH,
                  MAIN_HEIGHT,
                  INIT_MAIN_WIDTH,
                  INIT_MAIN_HEIGHT;
extern double ROT_Y,
              NOR_Y,
              ROT_X,
              MOV_X,
              MOV_Y,
              ZOOM;

/** \brief Handles keyboard events.
 *
 * \param key const uchar An identifier of the key.
 * \param x int The current x position of the cursor.
 * \param y int The current y position of the cursor.
 * \return void
 *
 */
void keyboard(const uchar key, int x, int y)
{
    switch(key)
    {
    case 27: /**< Esc. */
        ESC_CLICKED = 1;
        close_all();
        break;
    case 'm':
        TRANSLATE = !TRANSLATE;
        if(TRANSLATE)
            glutSetCursor(GLUT_CURSOR_INFO);
        else
            glutSetCursor(GLUT_CURSOR_CROSSHAIR);
        break;
    case 'R':
        ROT_Y = NOR_Y = ROT_X = MOV_X = MOV_Y = 0.;
        ZOOM = DEFAULT_ZOOM;
        break;
    default:
        ESC_CLICKED = x = y = 0;
    }
}

/** \brief Handles special keyboard events.
 *
 * \param a_keys const int An identifier of the special key.
 * \param x int The current x position of the cursor.
 * \param y int The current y position of the cursor.
 * \return void
 *
 */
void arrow_keys(const int a_keys, int x, int y)
{
    const double alpha_inc = .1;
    switch(a_keys)
    {
    case GLUT_KEY_F1:
        STOP_ROT = !STOP_ROT;
        break;
    case GLUT_KEY_F12:
        glutFullScreen();
        break;
    case GLUT_KEY_F11:
        glutReshapeWindow(INIT_MAIN_WIDTH, INIT_MAIN_HEIGHT);
        glutPositionWindow(25, 25);
        break;
    case GLUT_KEY_UP:
        if(BIN_SPHERE3_ALPHA <= 1. - alpha_inc)
            BIN_SPHERE3_ALPHA += .1;
        break;
    case GLUT_KEY_DOWN:
        if(BIN_SPHERE3_ALPHA >= 0. + alpha_inc)
            BIN_SPHERE3_ALPHA -= .1;
        break;
    case GLUT_KEY_LEFT:
        DRAW_RAY_N++;
        break;
    case GLUT_KEY_RIGHT:
        if(DRAW_RAY_N > 0)
            DRAW_RAY_N--;
        break;
    default:
        x = y = 0;
    }
}

/** \brief Handles mouse events.
 *
 * \param mode const int The mode of the mouse, i.e. a motion with a pressed button.
 * \param button const int An identifier of the mousebutton.
 * \param state const int An identifier of the state of a mousebutton.
 * \param x const int The current x position of the cursor.
 * \param y const int The current y position of the cursor.
 * \return void
 *
 */
void trackball(const int mode, const int button, const int state, const int x, const int y)
{
    static double startMX = 0., startMY = 0.;
    double deltaMX = 0., deltaMY = 0.;
    switch(mode)
    {
    case KNOWMOUSEBUTTON: /**< Handle the buttons. */
        if(!button && !state && !TRANSLATE && !MOUSE_R_DOWN) /**< Rotate. */
        {
            MOUSE_L_DOWN = 1;
            startMY = x;
            startMY -= ROT_Y;
            startMX = y;
            startMX -= ROT_X;
        }
        else if(!button && !state && TRANSLATE && !MOUSE_R_DOWN) /**< Translate. */
        {
            MOUSE_L_DOWN = 1;
            startMY = y;
            startMY += MOV_Y * MAIN_HEIGHT / 4;
            startMX = x;
            startMX -= MOV_X * MAIN_WIDTH / 4.;
        }
        else if(!button && state)
            MOUSE_L_DOWN = 0;
        else if(button == 2 && !state && !MOUSE_L_DOWN) /**< Zoom. */
        {
            MOUSE_R_DOWN = 1;
            startMY = y;
        }
        else if(button == 2 && state)
            MOUSE_R_DOWN = 0;
        break;
    case MOUSEMOTION: /**< Handle the effect of a mouse motion. */
        if(MOUSE_L_DOWN && !TRANSLATE && !MOUSE_R_DOWN) /**< Rotate. */
        {
            deltaMX = y - startMX;
            deltaMY = x - startMY;
            ROT_X = deltaMX;
            ROT_Y = deltaMY;
        }
        else if(MOUSE_L_DOWN && TRANSLATE && !MOUSE_R_DOWN) /**< Translate. */
        {
            deltaMX = x - startMX;
            deltaMY = y - startMY;
            MOV_X = deltaMX / MAIN_WIDTH * 4.;
            MOV_Y = -deltaMY / MAIN_HEIGHT * 4.;
        }
        else if(MOUSE_R_DOWN && !MOUSE_L_DOWN) /**< Zoom. */
        {
            deltaMY = y - startMY;
            const double minzoom = .02;
            if(ZOOM >= minzoom)
                ZOOM *= (1. + deltaMY / 200.);
            else
                ZOOM = minzoom;
#if DEBUG_MOTION
            fprintf(stdout, "zoom: %f\n", ZOOM);
#endif
            startMY = y;
        }
        break;
    case RESET:
        startMX = startMY = deltaMX = deltaMY = 0.;
        break;
    }
}

/** \brief Handles the pure motion of the mouse for glutMotionFunc.
 *
 * \param x const int The current x position of the cursor.
 * \param y const int The current y position of the cursor.
 * \return void
 *
 */
void motion(const int x, const int y)
{
#if DEBUG_MOTION
    fprintf(stdout, "mouse motion: %d %d\n", x, y);
#endif
    trackball(MOUSEMOTION, 0, 0, x, y);
}

/** \brief Handles the pure motion of the mouse for glutPassiveMotionFunc.
 *
 * \param x const int The current x position of the cursor.
 * \param y const int The current y position of the cursor.
 * \return void
 *
 * Not necessary, but it creates a valid handle for glutPassiveMotionFunc.
 */
void passivemotion(const int x, const int y)
{
#if DEBUG_MOTION
    fprintf(stdout, "passive mouse motion: %d, %d\n", x, y);
#endif
}

/** \brief Handles and initiates the interaction with the mouse.
 *
 * \param button const int An identifier of the button pressed.
 * \param state const int An identifier of the state of the button.
 * \param x const int The current x position of the cursor.
 * \param y const int The current y position of the cursor.
 * \return void
 *
 */
void mouse(const int button, const int state, const int x, const int y)
{
#if DEBUG_MOTION
    fprintf(stdout, "button = %d, state = %d %d %d\n", button, state, x, y);
#endif
    trackball(KNOWMOUSEBUTTON, button, state, x, y);
}
