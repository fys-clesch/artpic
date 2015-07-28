#include "artpic.h"
#include "msg.h"

const cdoub n_vac = 1. + 0.i;
const double mu_vac = 1.;

/** \brief computes the refractive index of h2o
 *
 * \param lam double the wavelength
 * \param s const char* the unit of lam
 * \return cdoub the output
 *
 * http://refractiveindex.info/?group=LIQUIDS&material=Water
 * reference for the data: M. Daimon and A. Masumura
 * Measurement of the refractive index of distilled water from the near-infrared region to the ultraviolet region
 * Appl. Opt. 46, 3811-3820 (2007)
 * doi:10.1364/AO.46.003811
 * conditions: 20.0 deg
 * valid for 182 to 1129 nm
 */
cdoub rfrct_indx_h2o(double lam, const char *s)
{
	if(!strcmp(s, "pm")) lam *= 1e-6;
	else if(!strcmp(s, "nm")) lam *= 1e-3;
	else if(!strcmp(s, "um")) {}
	else if(!strcmp(s, "mm")) lam *= 1e3;
	else
	{
		error_msg("wrong unit. setting to 500 nm.", ERR_ARG);
		lam = 500e-3;
	}
	if(lam < 182e-3 || lam > 1129e-3) fprintf(stdout, "rfrct_indx_h2o only valid for wavelengths (182 to 1129) nm");
	double t = lam * lam, n;
	n = 5.684027565e-1 / (t - 5.101829712e-3) +
	    1.726177391e-1 / (t - 1.821153936e-2) +
	    2.086189578e-2 / (t - 2.620722293e-2) +
	    1.130748688e-1 / (t - 1.069792721e1);
	n *= t;
	n += 1.;
	if(n < 0.)
	{
		error_msg("imaginary value, returning 1.", ERR_ARG);
		return 1. + 0.i;
	}
	else return (sqrt(n) + 1.1e-7 * 1.i);
}

/** \brief computes the refractive index of air
 *
 * \param lam double the wavelength
 * \param s const char* the unit of lam
 * \return cdoub the output
 *
 * reference for the data: Philip E. Ciddor
 * Refractive index of air: new equations for the visible and near infrared
 * Appl. Optics 35, 1566-1573 (1996)
 * doi:10.1364/AO.35.001566
 * conditions: 15.0 deg, 101 325 Pa
 */
cdoub rfrct_indx_air(double lam, const char *s)
{
	if(!strcmp(s, "pm")) lam *= 1e-6;
	else if(!strcmp(s, "nm")) lam *= 1e-6;
	else if(!strcmp(s, "um")) {}
	else if(!strcmp(s, "mm")) lam *= 1e3;
	else
	{
		error_msg("wrong unit. setting to 500 nm.", ERR_ARG);
		lam = 500e-3;
	}
	if(lam < 200e-3) fprintf(stdout, "rfrct_indx_air only valid for wavelengths >= 200 nm");
	double t = 1. / (lam * lam), n;
	n = 1. +
	    5792105e-8 / (238.0185 - t) +
	    167917e-8 / (57.362 - t);
	if(n < 0.)
	{
		error_msg("imaginary value, returning 1.", ERR_ARG);
		return 1. + 0.i;
	}
	else return (n + 0.i);
}
