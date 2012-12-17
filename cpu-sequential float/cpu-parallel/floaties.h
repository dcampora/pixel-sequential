
#include <math.h>

#ifndef FLOATIES
#define FLOATIES 1

float 	f_m_maxXSlope			= 0.400;
float 	f_m_maxYSlope			= 0.300;
float 	f_m_maxZForRBeamCut		= 200.0;
float 	f_m_maxR2Beam			= 1.0;
int 	f_m_maxMissed			= 4;
float 	f_m_extraTol			= 0.150;
float 	f_m_maxChi2ToAdd		= 100.0;
float 	f_m_maxChi2SameSensor	= 16.0;
float   f_m_maxChi2Short		= 6.0 ;
float   f_m_maxChi2PerHit		= 16.0;
int 	f_m_sensNum				= 48;
float   f_w						= 0.050 / sqrt( 12. );

#endif 
