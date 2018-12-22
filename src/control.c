#include "artpic.h"
#include "viewer.h"
#include "control.h"

uchar lmouse_down = 0, rmouse_down = 0, stop_rot = 0, translate = 0, esc_pressed = 0;
uint draw_ray_n = 0;
double bin_sphere3_alpha = DEFAULT_BIN_SPHERE3_ALPHA;

extern const uint main_width, main_height, init_main_width, init_main_height;
extern double rotY, norY, rotX, movX, movY, zoom;

/** \brief Handles keyboard events
 *
 * \param key const uchar An identifier of the key
 * \param x int The current x position of the cursor
 * \param y int The current y position of the cursor
 * \return void
 *
 */
void keyboard(const uchar key, int x, int y)
{
    switch(key)
    {
    case 27: /**< Esc */
        esc_pressed = 1;
        close_all();
        break;
    case 'm':
        translate = !translate;
        if(translate) glutSetCursor(GLUT_CURSOR_INFO);
        else glutSetCursor(GLUT_CURSOR_CROSSHAIR);
        break;
    case 'R':
        rotY = norY = rotX = movX = movY = 0.;
        zoom = DEFAULT_ZOOM;
        break;
    default:
        esc_pressed = x = y = 0;
    }
}

/** \brief Handles special keyboard events
 *
 * \param a_keys const int An identifier of the special key
 * \param x int The current x position of the cursor
 * \param y int The current y position of the cursor
 * \return void
 *
 */
void arrow_keys(const int a_keys, int x, int y)
{
    const double alpha_inc = .1;
    switch(a_keys)
    {
    case GLUT_KEY_F1:
        stop_rot = !stop_rot;
        break;
    case GLUT_KEY_F12:
        glutFullScreen();
        break;
    case GLUT_KEY_F11:
        glutReshapeWindow(init_main_width, init_main_height);
        glutPositionWindow(25, 25);
        break;
    case GLUT_KEY_UP:
        if(bin_sphere3_alpha <= 1. - alpha_inc) bin_sphere3_alpha += .1;
        break;
    case GLUT_KEY_DOWN:
        if(bin_sphere3_alpha >= 0. + alpha_inc) bin_sphere3_alpha -= .1;
        break;
    case GLUT_KEY_LEFT:
        draw_ray_n++;
        break;
    case GLUT_KEY_RIGHT:
        if(draw_ray_n > 0) draw_ray_n--;
        break;
    default:
        x = y = 0;
    }
}

/** \brief Handles mouse events
 *
 * \param mode const int The mode of the mouse, i.e. a motion with a pressed button
 * \param button const int An identifier of the mousebutton
 * \param state const int An identifier of the state of a mousebutton
 * \param x const int The current x position of the cursor
 * \param y const int The current y position of the cursor
 * \return void
 *
 */
void trackball(const int mode, const int button, const int state, const int x, const int y)
{
    static double startMX = 0., startMY = 0.;
    double deltaMX = 0., deltaMY = 0.;
    switch(mode)
    {
    case KNOWMOUSEBUTTON: /**< Handle the buttons */
        if(!button && !state && !translate && !rmouse_down) /**< Rotate */
        {
            lmouse_down = 1;
            startMY = x;
            startMY -= rotY;
            startMX = y;
            startMX -= rotX;
        }
        else if(!button && !state && translate && !rmouse_down) /**< Translate */
        {
            lmouse_down = 1;
            startMY = y;
            startMY += movY * main_height / 4;
            startMX = x;
            startMX -= movX * main_width / 4.;
        }
        else if(!button && state) lmouse_down = 0;
        else if(button == 2 && !state && !lmouse_down) /**< Zoom */
        {
            rmouse_down = 1;
            startMY = y;
        }
        else if(button == 2 && state) rmouse_down = 0;
        break;
    case MOUSEMOTION: /**< Handle the effect of a mouse motion */
        if(lmouse_down && !translate && !rmouse_down) /**< Rotate */
        {
            deltaMX = y - startMX;
            deltaMY = x - startMY;
            rotX = deltaMX;
            rotY = deltaMY;
        }
        else if(lmouse_down && translate && !rmouse_down) /**< Translate */
        {
            deltaMX = x - startMX;
            deltaMY = y - startMY;
            movX = deltaMX / main_width * 4.;
            movY = -deltaMY / main_height * 4.;
        }
        else if(rmouse_down && !lmouse_down) /**< Zoom */
        {
            deltaMY = y - startMY;
            const double minzoom = .02;
            if(zoom >= minzoom)
                zoom *= (1. + deltaMY / 200.);
            else
                zoom = minzoom;
#if DEBUG_MOTION
            fprintf(stdout, "zoom: %f\n", zoom);
#endif
            startMY = y;
        }
        break;
    case RESET:
        startMX = startMY = deltaMX = deltaMY = 0.;
        break;
    }
}

/** \brief Handles the pure motion of the mouse for glutMotionFunc
 *
 * \param x const int The current x position of the cursor
 * \param y const int The current y position of the cursor
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

/** \brief Handles the pure motion of the mouse for glutPassiveMotionFunc
 *
 * \param x const int The current x position of the cursor
 * \param y const int The current y position of the cursor
 * \return void
 *
 * Not necessary, but it creates a valid handle for glutPassiveMotionFunc
 */
void passivemotion(const int x, const int y)
{
#if DEBUG_MOTION
    fprintf(stdout, "passive mouse motion: %d, %d\n", x, y);
#endif
}

/** \brief Handles and initiates the interaction with the mouse
 *
 * \param button const int An identifier of the button pressed
 * \param state const int An identifier of the state of the button
 * \param x const int The current x position of the cursor
 * \param y const int The current y position of the cursor
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
