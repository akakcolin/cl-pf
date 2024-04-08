#ifndef SCALAR_H
// Scalar Arithmetic Operations and Macros

#define SCALAR_H

#include<math.h>


/** Mathematical Constants **/

#define HALF      (0.5)
#define THIRD     (1./3.)
#define FOURTH    (0.25)
#define FIFTH     (0.2)
#define SIXTH     (1./6.)
#define SEVENTH   (1./7.)
#define EIGHTH    (0.125)
#define NINTH     (1./9.)
#define TENTH     (1./10.)
#define ELEVENTH  (1./11.)
#define TWELFTH   (1./12.)

/* extended precision to 60 digits */
#define SQRT2 (1.41421356237309504880168872420969807856967187537694807317668)
#define SQRT3 (1.73205080756887729352744634150587236694280525381038062805581)
#define SQRT5 (2.23606797749978969640917366873127623544061835961152572427090)
#define SQRT6 (2.44948974278317809819728407470589139196594748065667012843269)
#define SQRT7 (2.64575131106459059050161575363926042571025918308245018036833)
#define SQRT8 (2.82842712474619009760337744841939615713934375075389614635336)
#define CBRT2 (1.25992104989487316476721060727822835057025146470150798008198)
#define CBRT3 (1.44224957030740838232163831078010958839186925349935057754642)
#define CBRT4 (1.58740105196819947475170563927230826039149332789985300980829)
#define CBRT5 (1.70997594667669698935310887254386010986805511054305492438286)
#define CBRT6 (1.81712059283213965889121175632726050242821046314121967148133)
#define CBRT7 (1.91293118277238910119911683954876028286243905034587576621065)
#define CBRT9 (2.08008382305190411453005682435788538633780534037326210969759)
/* (-1+sqrt(5.))/2 */
#define GOLDEN_RATIO \
  (.618033988749894848204586834365638117720309179805762862135450)
#define GOLDEN_RATIO_RECIPROCAL (GOLDEN_RATIO + 1)

#ifndef PI
#define PI    (3.14159265358979323846264338327950288419716939937510582097494)
#endif
#define RADIAN_TO_DEGREE(r) ((r)/PI*180.)
#define DEGREE_TO_RADIAN(d) ((d)/180.*PI)

#define SIN(x) sin(DEGREE_TO_RADIAN(x))
#define COS(x) cos(DEGREE_TO_RADIAN(x))
#define TAN(x) tan(DEGREE_TO_RADIAN(x))
#define ASIN(x) RADIAN_TO_DEGREE(asin(x))
#define ACOS(x) RADIAN_TO_DEGREE(acos(x))
#define ATAN(x) RADIAN_TO_DEGREE(atan(x))

#ifndef EPS
/* IEEE 754 floating-point arithmetic relative error bound */
#define EPS (2.220446e-16)
#endif


#define EQ(x,y)  ( (x) == (y) )
#define GE(x,y)  ( (x) >= (y) )
#define GT(x,y)  ( (x) >  (y) )
#define LE(x,y)  ( (x) <= (y) )
#define LT(x,y)  ( (x) <  (y) )
#define ABS(x)  ((x)>0?(x):-(x))

#ifndef MIN 
#define MIN(x,y)      ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y)      ((x)>(y)?(x):(y))
#endif

#define SQUARE(i)     ((i)*(i))
#define CUBE(i)       ((i)*(i)*(i))
#define QUAD(i)       ((i)*(i)*(i)*(i))
#define CEIL(i,j)     (((i)<=0)?0:((i)-1)/(j)+1)
#define SEIL(i,j)     ((size_t)CEIL(i,j))
#define ROUNDUP_TO(i,j)         (CEIL(i,j)*(j))
#define SEPARATION(x,y)         (fabs((x)-(y)))
#define SMALLSEPARATION(x,y,z)  (fabs((x)-(y))<z)


#define SPHERE_VOLUME(radius) (4.*PI/3.*CUBE(radius))
#define safe_avg(sum,count) ((double)(sum)/(((count)!=0)?(count):1))


#define DISTANCE2(x,y,z) (SQUARE(x)+SQUARE(y)+SQUARE(z))
#define DISTANCE(x,y,z)   sqrt(DISTANCE2(x,y,z))
#define DISTANCED2(x,y)  (SQUARE(x)+SQUARE(y))
#define DISTANCE2D(x,y)   sqrt(DISTANCED2(x,y))

        
#endif
