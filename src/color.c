#include "artpic.h"
#include "msg.h"
#include "auxf.h"

/**
 * A table to convert wavelengths from 380 to 780 nm to RGB values corresponding to a CIE 1964.
 * Standard observer and a CIE system with IlluminantE.
 */
const double LIGHTTOCOLOR[401][3] =
{
    {0.0729182, 0, 1}, {0.0468786, 0.0022235, 1}, {0.0436796, 0.00703644, 1},
    {0.0441658, 0.00696042, 1}, {0.0444142, 0.00691748, 1}, {0.0443842, 0.00699808, 1},
    {0.0442733, 0.00699036, 1}, {0.0441636, 0.00684462, 1}, {0.0439286, 0.0066874, 1},
    {0.0435627, 0.00667337, 1}, {0.0430219, 0.00675771, 1}, {0.0424597, 0.00683335, 1},
    {0.0418853, 0.00688686, 1}, {0.0412837, 0.00691953, 1}, {0.0406319, 0.00694846, 1},
    {0.0399816, 0.00698257, 1}, {0.0393239, 0.00701798, 1}, {0.0386318, 0.00705692, 1},
    {0.0378973, 0.00710174, 1}, {0.0371471, 0.00714259, 1}, {0.036375, 0.00719247, 1},
    {0.0355807, 0.00725098, 1}, {0.0347599, 0.00731813, 1}, {0.0339001, 0.00740081, 1},
    {0.0329977, 0.00750304, 1}, {0.0320592, 0.00762257, 1}, {0.031087, 0.00775328, 1},
    {0.0300575, 0.00789467, 1}, {0.0289466, 0.00804553, 1}, {0.0277435, 0.00820507, 1},
    {0.0264415, 0.00837577, 1}, {0.0250492, 0.00855768, 1}, {0.0235772, 0.00875545, 1},
    {0.0220305, 0.00897498, 1}, {0.0204182, 0.00921805, 1}, {0.0187449, 0.00948817, 1},
    {0.017018, 0.00978925, 1}, {0.0152488, 0.0101287, 1}, {0.0134455, 0.0105149, 1},
    {0.011612, 0.0109553, 1}, {0.00974829, 0.0114544, 1}, {0.00785198, 0.0120148, 1},
    {0.00592294, 0.0126345, 1}, {0.00395913, 0.0133117, 1}, {0.00195854, 0.0140444, 1},
    {0, 0.0149148, 1}, {0, 0.0178041, 1}, {0, 0.0207648, 1},
    {0, 0.0237763, 1}, {0, 0.0268202, 1}, {0, 0.0298781, 1},
    {0, 0.0329447, 1}, {0, 0.0360394, 1}, {0, 0.0391928, 1},
    {0, 0.0424392, 1}, {0, 0.0458097, 1}, {0, 0.0493238, 1},
    {0, 0.0529749, 1}, {0, 0.0567468, 1}, {0, 0.0606263, 1},
    {0, 0.0646056, 1}, {0, 0.0686942, 1}, {0, 0.0729281, 1},
    {0, 0.0773545, 1}, {0, 0.0820217, 1}, {0, 0.0869737, 1},
    {0, 0.0922327, 1}, {0, 0.0977796, 1}, {0, 0.103572, 1},
    {0, 0.109565, 1}, {0, 0.115715, 1}, {0, 0.122023, 1},
    {0, 0.128564, 1}, {0, 0.13546, 1}, {0, 0.142847, 1},
    {0, 0.150859, 1}, {0, 0.159583, 1}, {0, 0.169006, 1},
    {0, 0.179056, 1}, {0, 0.189648, 1}, {0, 0.2007, 1},
    {0, 0.212194, 1}, {0, 0.22425, 1}, {0, 0.237068, 1},
    {0, 0.250877, 1}, {0, 0.265914, 1}, {0, 0.282363, 1},
    {0, 0.30029, 1}, {0, 0.319715, 1}, {0, 0.340677, 1},
    {0, 0.363245, 1}, {0, 0.387488, 1}, {0, 0.413392, 1},
    {0, 0.440851, 1}, {0, 0.469687, 1}, {0, 0.499657, 1},
    {0, 0.530612, 1}, {0, 0.562708, 1}, {0, 0.596283, 1},
    {0, 0.631701, 1}, {0, 0.669272, 1}, {0, 0.709085, 1},
    {0, 0.75084, 1}, {0, 0.794032, 1}, {0, 0.838126, 1},
    {0, 0.882635, 1}, {0, 0.927301, 1}, {0, 0.97226, 1},
    {0, 1, 0.982467}, {0, 1, 0.939492}, {0, 1, 0.899091},
    {0, 1, 0.860968}, {0, 1, 0.825036}, {0, 1, 0.791264},
    {0, 1, 0.759596}, {0, 1, 0.729942}, {0, 1, 0.702206},
    {0, 1, 0.676301}, {0, 1, 0.652141}, {0, 1, 0.629635},
    {0, 1, 0.608673}, {0, 1, 0.589109}, {0, 1, 0.570731},
    {0, 1, 0.553308}, {0, 1, 0.536624}, {0, 1, 0.520489},
    {0, 1, 0.504791}, {0, 1, 0.489541}, {0, 1, 0.474812},
    {0, 1, 0.460681}, {0, 1, 0.447217}, {0, 1, 0.434423},
    {0, 1, 0.422197}, {0, 1, 0.410376}, {0, 1, 0.398794},
    {0, 1, 0.387287}, {0, 1, 0.375732}, {0, 1, 0.36407},
    {0, 1, 0.352278}, {0, 1, 0.340329}, {0, 1, 0.328191},
    {0, 1, 0.315819}, {0, 1, 0.303166}, {0, 1, 0.290179},
    {0, 1, 0.276808}, {0, 1, 0.263003}, {0, 1, 0.248724},
    {0, 1, 0.233947}, {0, 1, 0.218669}, {0, 1, 0.202898},
    {0, 1, 0.186652}, {0, 1, 0.169931}, {0, 1, 0.152682},
    {0, 1, 0.134809}, {0, 1, 0.116185}, {0, 1, 0.0966631},
    {0, 1, 0.0760761}, {0, 1, 0.0542541}, {0, 1, 0.0310102},
    {0, 1, 0.00612291}, {0.0216832, 1, 0.00146317}, {0.0476887, 1, 0.000404092},
    {0.0750332, 1, 0}, {0.103439, 1, 0}, {0.132518, 1, 0},
    {0.162263, 1, 0}, {0.192697, 1, 0}, {0.223904, 1, 0},
    {0.256, 1, 0}, {0.289093, 1, 0}, {0.323276, 1, 0},
    {0.35861, 1, 0}, {0.395105, 1, 0}, {0.43276, 1, 0},
    {0.471588, 1, 0}, {0.511624, 1, 0}, {0.552935, 1, 0},
    {0.595622, 1, 0}, {0.639813, 1, 0}, {0.68565, 1, 0},
    {0.733277, 1, 0}, {0.782842, 1, 0}, {0.834526, 1, 0},
    {0.888453, 1, 0}, {0.944723, 1, 0}, {1, 0.99659, 0},
    {1, 0.939282, 0}, {1, 0.886177, 0}, {1, 0.836922, 0},
    {1, 0.791185, 0}, {1, 0.748656, 0}, {1, 0.709035, 0},
    {1, 0.672038, 0}, {1, 0.637415, 0}, {1, 0.604953, 0},
    {1, 0.574469, 0}, {1, 0.545798, 0}, {1, 0.518788, 0},
    {1, 0.493304, 0}, {1, 0.469237, 0}, {1, 0.446498, 0},
    {1, 0.425014, 0}, {1, 0.40471, 0}, {1, 0.385509, 0},
    {1, 0.367332, 0}, {1, 0.350097, 0}, {1, 0.33374, 0},
    {1, 0.31822, 0}, {1, 0.30351, 0}, {1, 0.289582, 0},
    {1, 0.276404, 0}, {1, 0.263938, 0}, {1, 0.25213, 0},
    {1, 0.240926, 0}, {1, 0.230278, 0}, {1, 0.220148, 0},
    {1, 0.210502, 0}, {1, 0.201313, 0}, {1, 0.192559, 0},
    {1, 0.184218, 0}, {1, 0.176272, 0}, {1, 0.168704, 0},
    {1, 0.161494, 0}, {1, 0.154623, 0}, {1, 0.148073, 0},
    {1, 0.141829, 0}, {1, 0.135874, 0}, {1, 0.130197, 0},
    {1, 0.124784, 0}, {1, 0.119623, 0}, {1, 0.114706, 0},
    {1, 0.110022, 0}, {1, 0.10556, 0}, {1, 0.10131, 0},
    {1, 0.0972619, 0}, {1, 0.0934088, 0}, {1, 0.0897455, 0},
    {1, 0.0862715, 0}, {1, 0.0829903, 0}, {1, 0.0799075, 0},
    {1, 0.0770295, 0}, {1, 0.0743573, 0}, {1, 0.0718778, 0},
    {1, 0.0695702, 0}, {1, 0.0674101, 0}, {1, 0.0653713, 0},
    {1, 0.0634294, 0}, {1, 0.0615654, 0}, {1, 0.059762, 0},
    {1, 0.0579999, 0}, {1, 0.0562596, 0}, {1, 0.0545262, 0},
    {1, 0.0527985, 0}, {1, 0.0510852, 0}, {1, 0.0494008, 0},
    {1, 0.0477618, 0}, {1, 0.0461815, 0}, {1, 0.0446603, 0},
    {1, 0.0431919, 0}, {1, 0.041768, 0}, {1, 0.040378, 0},
    {1, 0.0390157, 0}, {1, 0.037682, 0}, {1, 0.0363827, 0},
    {1, 0.0351252, 0}, {1, 0.0339159, 0}, {1, 0.0327623, 0},
    {1, 0.0316726, 0}, {1, 0.0306559, 0}, {1, 0.0297249, 0},
    {1, 0.0288921, 0}, {1, 0.0281638, 0}, {1, 0.0275363, 0},
    {1, 0.0269994, 0}, {1, 0.0265378, 0}, {1, 0.02613, 0},
    {1, 0.0257569, 0}, {1, 0.0254104, 0}, {1, 0.0250822, 0},
    {1, 0.0247649, 0}, {1, 0.0244526, 0}, {1, 0.0241389, 0},
    {1, 0.0238239, 0}, {1, 0.0235076, 0}, {1, 0.0231983, 0},
    {1, 0.0229027, 0}, {1, 0.022624, 0}, {1, 0.0223667, 0},
    {1, 0.022126, 0}, {1, 0.0219045, 0}, {1, 0.0217008, 0},
    {1, 0.0215083, 0}, {1, 0.021326, 0}, {1, 0.0211569, 0},
    {1, 0.020999, 0}, {1, 0.0208489, 0}, {1, 0.0207085, 0},
    {1, 0.0205786, 0}, {1, 0.0204621, 0}, {1, 0.020352, 0},
    {1, 0.02025, 0}, {1, 0.0201563, 0}, {1, 0.0200789, 0},
    {1, 0.0200136, 0}, {1, 0.0199458, 0}, {1, 0.0198934, 0},
    {1, 0.0198548, 0}, {1, 0.0198104, 0}, {1, 0.0197635, 0},
    {1, 0.0197343, 0}, {1, 0.0197043, 0}, {1, 0.0196648, 0},
    {1, 0.0196201, 0}, {1, 0.019583, 0}, {1, 0.0195336, 0},
    {1, 0.0194991, 0}, {1, 0.0194792, 0}, {1, 0.0194674, 0},
    {1, 0.0194376, 0}, {1, 0.0194198, 0}, {1, 0.0194365, 0},
    {1, 0.0194656, 0}, {1, 0.0194712, 0}, {1, 0.0194739, 0},
    {1, 0.0195032, 0}, {1, 0.01952, 0}, {1, 0.0195231, 0},
    {1, 0.0195016, 0}, {1, 0.0194592, 0}, {1, 0.0194627, 0},
    {1, 0.0195096, 0}, {1, 0.019553, 0}, {1, 0.0195173, 0},
    {1, 0.0195274, 0}, {1, 0.0195643, 0}, {1, 0.019645, 0},
    {1, 0.0197231, 0}, {1, 0.0198551, 0}, {1, 0.0199479, 0},
    {1, 0.0200346, 0}, {1, 0.0199948, 0}, {1, 0.0199501, 0},
    {1, 0.0199991, 0}, {1, 0.0201587, 0}, {1, 0.0201773, 0},
    {1, 0.0201283, 0}, {1, 0.020129, 0}, {1, 0.0201023, 0},
    {1, 0.0203239, 0}, {1, 0.0204411, 0}, {1, 0.0204403, 0},
    {1, 0.0205791, 0}, {1, 0.0206992, 0}, {1, 0.0209655, 0},
    {1, 0.0211085, 0}, {1, 0.0211556, 0}, {1, 0.0211232, 0},
    {1, 0.0211336, 0}, {1, 0.0214109, 0}, {1, 0.0219011, 0},
    {1, 0.0219328, 0}, {1, 0.0219441, 0}, {1, 0.0226185, 0},
    {1, 0.0227448, 0}, {1, 0.0227483, 0}, {1, 0.0227454, 0},
    {1, 0.0221541, 0}, {1, 0.0219891, 0}, {1, 0.0220109, 0},
    {1, 0.0213101, 0}, {1, 0.0208097, 0}, {1, 0.020094, 0},
    {1, 0.0201024, 0}, {1, 0.0216923, 0}, {1, 0.0226198, 0},
    {1, 0.0227534, 0}, {1, 0.0224353, 0}, {1, 0.0229141, 0},
    {1, 0.0229933, 0}, {1, 0.0233349, 0}, {1, 0.0252097, 0},
    {1, 0.0267292, 0}, {1, 0.027483, 0}, {1, 0.02809, 0},
    {1, 0.0303177, 0}, {1, 0.0296452, 0}, {1, 0.0289016, 0},
    {1, 0.0244465, 0}, {1, 0.0221457, 0}, {1, 0.0182286, 0},
    {1, 0, 0.0356377}, {1, 0, 0.146942}, {1, 0, 0.179703},
    {1, 0, 0.181322}, {1, 0, 0.181333}, {1, 0, 0.181333},
    {1, 0, 0.181333}, {1, 0, 0.181333}, {1, 0, 0.181333},
    {1, 0, 0.181333}, {1, 0, 0.181333}
};

const double CWHITE[4] = {1., 1., 1., 1.},
             CBLACK[4] = {0., 0., 0., 1.},
             CTRANS[4] = {1., 1., 1., .01},
             CGREY[4] = {.8, .8, .8, .5},
             CLIGHTGREY[4] = {.9, .9, .9, .5},
             CRED[4] = {1., .01, .05, 1.},
             CGREEN[4] = {.01, 1., .01, 1.},
             CBLUE[4] = {.01, .1, 1., 1.};

const float LIGHT_AMBIENT[4] = {0., 0., 0., .8},
            LIGHT_DIFFUSE[4] = {.8, .8, .8, .8},
            LIGHT_SPECULAR[4] = {.8, .8, .8, .8},
            LIGHT_POSITION[4] = {2., 5., 5., 0.},
            MAT_AMBIENT[4] = {.7, .7, .7, .8},
            MAT_DIFFUSE[4] = {.8, .8, .8, .8},
            MAT_SPECULAR[4] = {1., 1., 1., .8},
            HIGH_SHININESS[1] = {100.};

/** \brief Norms a double array of RGB values to the maximum.
 *
 * \param rgb double* The RGB array.
 * \return void
 *
 */
void norm_rgb(double *rgb)
{
    double max = MAXOF(rgb[0], MAXOF(rgb[1], rgb[2]));
    if(max > 0.)
    {
        max = 1. / max;
        rgb[0] *= max;
        rgb[1] *= max;
        rgb[2] *= max;
    }
}

/** \brief Makes the minimum value of a double array non-negative.
 *
 * \param rgb double* The RGB array.
 * \return void
 *
 */
void constrain_rgb(double *rgb)
{
    double w;
    w = (rgb[0] < 0.) ? rgb[0] : 0.;
    w = (rgb[1] < w) ? rgb[1] : w;
    w = (rgb[2] < w) ? rgb[2] : w;
    if(w != 0.)
    {
        rgb[0] -= w;
        rgb[1] -= w;
        rgb[2] -= w;
    }
}

/** \brief Converts a RGB value to a HSV value.
 *
 * \param rgb double* An array containing the RGB values.
 * \param hsv double* An array to store the HSV values.
 * \return void
 *
 * This function works for different scopes, that is RGB and HSV can overlap.
 * RGB might be changed due to constrain_rgb and norm_rgb.
 */
void rgbtohsv(double *rgb, double *hsv)
{
    double max, min;
    if(rgb[0] < 0. || rgb[1] < 0. || rgb[2] < 0.)
        constrain_rgb(rgb);
    if(rgb[0] > 1. || rgb[1] > 1. || rgb[2] > 1.)
        norm_rgb(rgb);
    max = MAXOF(rgb[0], MAXOF(rgb[1], rgb[2]));
    min = MINOF(rgb[0], MINOF(rgb[1], rgb[2]));
    if(max == min)
        hsv[0] = 0.;
    else if(max == rgb[0])
        hsv[0] = 60.*((rgb[1] - rgb[2]) / (max - min));
    else if(max == rgb[1])
        hsv[0] = 60.*(2. + (rgb[2] - rgb[0]) / (max - min));
    else if(max == rgb[2])
        hsv[0] = 60.*(4. + (rgb[0] - rgb[1]) / (max - min));
    else
        error_msg("this should not happen", ERR_ARG);
    if(hsv[0] < 0.)
        hsv[0] += 360.;
    if(max == 0.)
        hsv[1] = 0.;
    else
        hsv[1] = (max - min) / max;
    if(hsv[0] != 0. && hsv[1] != 0.)
        hsv[2] = max;
    else
        hsv[2] = 0.;
}

/** \brief Converts a HSV value to a RGB value.
 *
 * \param hsv double* An array containing the HSV values.
 * \param rgb double* An array to store the RGB values.
 * \return void
 *
 * This function works for different scopes, that is RGB and HSV can overlap.
 */
void hsvtorgb(double *hsv, double *rgb)
{
    double f, p, q, t;
    uint h;
    h = (uint)floor(hsv[0] / 60.);
    f = hsv[0] / 60. - h;
    p = hsv[2] * (1. - hsv[1]);
    q = hsv[2] * (1. - hsv[1] * f);
    t = hsv[2] * (1. - hsv[1] * (1. - f));
    if(h == 0 || h == 6)
    {
        rgb[0] = hsv[0];
        rgb[1] = t;
        rgb[2] = p;
    }
    else if(h == 1)
    {
        rgb[0] = q;
        rgb[1] = hsv[0];
        rgb[2] = p;
    }
    else if(h == 2)
    {
        rgb[0] = p;
        rgb[1] = hsv[0];
        rgb[2] = t;
    }
    else if(h == 3)
    {
        rgb[0] = p;
        rgb[1] = q;
        rgb[2] = hsv[0];
    }
    else if(h == 4)
    {
        rgb[0] = t;
        rgb[1] = p;
        rgb[2] = hsv[0];
    }
    else if(h == 5)
    {
        rgb[0] = hsv[0];
        rgb[1] = p;
        rgb[2] = q;
    }
    else
    {
        error_msg("can not convert your hsv value. normed to degrees?", ERR_ARG);
        rgb[0] = rgb[1] = rgb[2] = 0.;
    }
}

/** \brief Colours the given bin_hit_screen variable with a colour function.
 *
 * \param ptc patch3* The patches that have to be coloured.
 * \param nptc const uint The number of patches.
 * \param *bhs const bin_hit_screen const The variable from which the results are taken which are going to used for colouring.
 * \param cfun const colourfun The function to be used for mapping.
 * \param alpha const double The alpha channel of the RGB colour.
 * \param type const bin_hit_print_type specifies The quantity that should be mapped to RGB space.
 * \return void
 *
 * The patch and the bin size have to be the same (this is tested)
 * Using a double const for bhs to make sure that whether the pointer nor the 'data of it' will be changed.
 */
void colour_bin_patch3(patch3 *ptc, const uint nptc, const bin_hit_screen *const bhs,
                       const colourfun cfun, const double alpha, const bin_hit_print_type type)
{
    uint i;
    double max, min = DBL_MAX, diff = 0xDEAD,
                t_logminval = 1., /**< A temporary minimum value for the log-scaling.  */
                **fptr;
    uchar logmod;
    switch(type)
    {
    case E_FIELD_AMPLITUDE:
        fptr = (*bhs).amp;
        max = (*bhs).amp_max;
        logmod = 1;
        break;
    case POLARISATION_DENSITY:
        fptr = (*bhs).pol_dens;
        max = (*bhs).pol_max;
        logmod = 0;
        break;
    case MOD_E_FIELD_AMPLITUDE:
        fptr = (*bhs).mod_amp;
        max = (*bhs).mod_amp_max;
        logmod = 1;
        break;
    case INTENSITY:
        fptr = (*bhs).ray_int;
        max = (*bhs).ray_int_max;
        logmod = 1;
        break;
    default:
        error_msg("unrecognised option. choosing e-field amplitude.", ERR_ARG);
        fptr = (*bhs).amp;
        max = (*bhs).amp_max;
        logmod = 1;
    }
    for(i = 0; i < (*bhs).nbin; i++)
        if(min > fptr[(*bhs).idx[i].ia][(*bhs).idx[i].ib])
            min = fptr[(*bhs).idx[i].ia][(*bhs).idx[i].ib];
    diff = fabs(max - min);
    if(diff == 0.)
        error_msg("could not detect any event on the screen", ERR_ARG);
    if(logmod)
    {
        if(min < 0.)
        {
            error_msg("log-mode is not possible", ERR_ARG);
            logmod = 0;
        }
        else
        {
            const double logminval = 1e-99; /**< The minimum value used for the log-scaling. */
            if(min > logminval) t_logminval = min;
            else t_logminval = logminval;
            min = (min > t_logminval) ? log10(min) : log10(t_logminval);
            diff = log10(max) - min;
        }
    }
    if(cfun == SINCOS_MAP)
        for(i = 0; i < nptc; i++)
        {
            double nval = fptr[(*bhs).idx[i].ia][(*bhs).idx[i].ib];
            if(logmod)
                nval = (nval > t_logminval) ? log10(nval) : min;
            nval = M_PI * (nval - min) / diff;
            if(nval >= 0. && nval < M_PIh)
            {
                ptc[i].rgba[0] = 0.;
                sincosd(nval, &ptc[i].rgba[2], &ptc[i].rgba[1]);
                ptc[i].rgba[3] = alpha;
            }
            else if(nval >= M_PIh && nval <= M_PI)
            {
                sincosd(nval, &ptc[i].rgba[0], &ptc[i].rgba[1]);
                ptc[i].rgba[0] = -ptc[i].rgba[0];
                ptc[i].rgba[2] = 0.;
                ptc[i].rgba[3] = alpha;
            }
            else
                memcpy(ptc[i].rgba, CGREY, 4 * sizeof(double));
        }
    else if(cfun == VISIBLE_LIGHT)
    {
        for(i = 0; i < nptc; i++)
        {
            double nval = fptr[(*bhs).idx[i].ia][(*bhs).idx[i].ib];
            if(logmod)
                nval = (nval > t_logminval) ? log10(nval) : min;
            nval = (nval - min) / diff;
            if(nval >= 0. && nval <= 1.)
            {
                uint t1 = (uint)floor(nval * 400.);
                ptc[i].rgba[0] = LIGHTTOCOLOR[t1][0];
                ptc[i].rgba[1] = LIGHTTOCOLOR[t1][1];
                ptc[i].rgba[2] = LIGHTTOCOLOR[t1][2];
                ptc[i].rgba[3] = alpha;
            }
            else
                memcpy(ptc[i].rgba, CGREY, 4 * sizeof(double));
        }
    }
    else
        error_msg("unknown colourfun option", ERR_ARG);
    if(type == POLARISATION_DENSITY)
        for(i = 0; i < nptc; i++)
            if(fptr[(*bhs).idx[i].ia][(*bhs).idx[i].ib] < DBL_EPSILON2)
            {
                ptc[i].rgba[0] = ptc[i].rgba[1] = ptc[i].rgba[2] = 0.;
                ptc[i].rgba[3] = alpha;
            }
}

/** \brief Maps a single double value to a RGBA value.
 *
 * \param rgba double* The RGBA value to be produced.
 * \param nval double The value to be mapped to the RGBA space.
 * \param min const double The minimum value to be mapped.
 * \param max const double The maximum value to be mapped.
 * \param cfun const colourfun The function to be used for mapping.
 * \return void
 *
 */
void map_cfun(double *rgba, double nval, const double min, const double max, const colourfun cfun)
{
    double diff = max - min;
    if(cfun == SINCOS_MAP)
    {
        nval = M_PI * (nval - min) / diff;
        if(nval >= 0. && nval < M_PIh)
        {
            rgba[0] = 0.;
            sincosd(nval, &rgba[2], &rgba[1]);
            rgba[3] = 1.;
        }
        else if(nval >= M_PIh && nval <= M_PI)
        {
            sincosd(nval, &rgba[0], &rgba[1]);
            rgba[0] = -rgba[0];
            rgba[2] = 0.;
            rgba[3] = 1.;
        }
        else
            memcpy(rgba, CGREY, 4 * sizeof(double));
    }
    else if(cfun == VISIBLE_LIGHT)
    {
        nval = (nval - min) / diff;
        if(nval >= 0. && nval <= 1.)
        {
            uint t1 = (uint)floor(nval * 400.);
            rgba[0] = LIGHTTOCOLOR[t1][0];
            rgba[1] = LIGHTTOCOLOR[t1][1];
            rgba[2] = LIGHTTOCOLOR[t1][2];
            rgba[3] = 1.;
        }
        else
            memcpy(rgba, CGREY, 4 * sizeof(double));
    }
    else
        error_msg("wrong identifier", ERR_ARG);
}
