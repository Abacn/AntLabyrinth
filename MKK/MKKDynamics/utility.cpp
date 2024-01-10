#include "utility.h"
#define _USE_MATH_DEFINES
#include <math.h>

// global random generator
std::mt19937_64 rg;
std::uniform_real_distribution <> rrand;

// get diameter range (1-x) - (1+x)
double diameterRange(double polydispersity)
{
  if (polydispersity == 0.0) return 0.0;
  else
  {
    // approximate as K^2 = x^2 / 3 + x^4/5 + x^6/7
    double K = polydispersity * polydispersity;
    double dtmp = pow(-427.0 - 3375.0 * K + 5.0 * sqrt(5.0 * (2765.0 + 23058.0 * K + 91125.0 * K * K)), 1.0 / 3.0);
    double xsqr = -(7.0 / 15.0) + (6.0 * pow(2.0 * 7.0 * 7.0, 1.0 / 3.0)) / (5.0 * dtmp)
                    - 1.0 / 15.0 * pow(7.0 / 2.0, 1.0 / 3.0) * dtmp;
    return sqrt(xsqr);
  }
}

// get average volume of sphere in radius (1-x) to (1+x)
// and distribution P(r) ~ r^(-3)
double averagesprvolume(double rx)
{
	double factor = 1.0;
	if (0.0 != rx)
	{
#if (DIM == 2)
	  factor = ((rx * rx - 1) * (rx * rx - 1) * atan(rx)) / (rx);
#else
	  factor = (pow(1 + rx, DIM) * (1 - rx) * (1 - rx) - pow(1 - rx, DIM) * (1 + rx) * (1 + rx)) / (2 * rx * (DIM - 2));
#endif
	}
	return VOLUMESPHERE*factor;
}

// get r from packing fraction
double gerRfrompfrac(double vspheres, double rx, int N)
{
  // This formula has finite error in finite N
  // return pow(vspheres / N / averagesprvolume(rx), 1.0 / DIM);
  // Note. Phi = 2^d * pfrac

  double rtom2min = 1.0 / ((1.0 + rx) * (1.0 + rx));
  double rtom2step = (1.0 / ((1.0 - rx) * (1.0 - rx)) - 1.0 / ((1.0 + rx) * (1.0 + rx))) / (N - 1);
  double totV = 0;
  for (int i = 0; i < N; i++)
  {
    totV += pow(rtom2min + i * rtom2step, -0.5*DIM );
  }
  return pow(vspheres / (VOLUMESPHERE * totV), 1.0 / DIM);
}

// lattice point on Dn lattice
void DnLatticePoint(double p[DIM], int32_t mirror[DIM])
{
  int rp, gind = 0, summr = 0;
  double maxdelta = 0, dtmp;
  for(rp=0; rp<DIM; ++rp)
  {
    // fast round eqv mirror[rp] = round(p[rp]);
    dtmp = p[rp] + 6755399441055744.0;
    mirror[rp] = reinterpret_cast<int32_t&>(dtmp);
    summr += mirror[rp];
    p[rp] -= mirror[rp];
  }
  if(summr%2 != 0)
  {
    for(rp=0; rp<DIM; ++rp)
    {
      dtmp = fabs(p[rp]);
      if(maxdelta < dtmp)
      {
        maxdelta = dtmp;
        gind = rp;
      }
    }
    // use next nearest neighbor image
    if(p[gind] < 0)
    {
      --mirror[gind];
      ++p[gind];
    }
    else
    {
      ++mirror[gind];
      --p[gind];
    }
  }
}

void setminimg(double x[DIM])
{
#if (1==BOX_TYPE)
  double dtmp;
  for (int k = 0; k < DIM; ++k)
  {
    // fast round
    dtmp = x[k] + 6755399441055744.0;
    x[k] -= reinterpret_cast<int32_t&>(dtmp);
  }
#else
  int32_t itmp[DIM];
  DnLatticePoint(x, itmp);
#endif
}
