
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <map>

#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"

using namespace std;
using namespace tbb;

#ifndef TYPES
#define TYPES 1

typedef struct {
	int startPosition;
	int hitsNum;
	float z;
} f_sensorInfo;

struct frack {
	float m_x0;
	float m_tx;
	float m_y0;
	float m_ty;

	float m_s0;
	float m_sx;
	float m_sz;
	float m_sxz;
	float m_sz2;

	float m_u0;
	float m_uy;
	float m_uz;
	float m_uyz;
	float m_uz2;
	
	int trackHitsNum;

	vector<int> hits;
};

// TODO: Create a file with floats instead of doubles

extern int* f_numHitsEvents;
extern int* f_sensor_hitStarts;
extern int* f_sensor_hitsNums;
extern float* f_sensor_Zs;
extern int* f_hit_ids;
extern int* f_track_ids;
extern int* f_hit_sensorNums;
extern int* f_hit_isUseds;
extern float* f_hit_Xs;
extern float* f_hit_Ys;
extern float* f_hit_Zs;
extern float* f_hit_Ws;

extern int num_events;
extern int max_hits;
extern int sens_num;

extern float 	f_m_maxXSlope;
extern float 	f_m_maxYSlope;
extern float 	f_m_maxZForRBeamCut;
extern float 	f_m_maxR2Beam;
extern int 		f_m_maxMissed;
extern float 	f_m_extraTol;
extern float 	f_m_maxChi2ToAdd;
extern float 	f_m_maxChi2SameSensor;
extern float	f_m_maxChi2Short;
extern float	f_m_maxChi2PerHit;
extern int 		f_m_sensNum;
extern float	f_w;

extern map<string, float> calcResults(vector<float> times);
extern void printStatistics(map<string, float> statistics, string pre);

extern vector< vector<frack> > parallel_fracks_vector;

class f_TBBSearchByPair {
public:
	void operator() (const blocked_range<int>& r) const;
};

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

float f_zBeam(frack *tr);
float f_r2AtZ( float z , frack *tr);
void f_solve (frack *tr);
float f_chi2Hit( float x, float y, float hitX, float hitY, float hitW);
float f_xAtHit(frack *tr, float z );
float f_yAtHit( frack *tr, float z  );
float f_chi2Track(frack *tr, int offset);
float f_chi2(frack *t);
bool f_addHitsOnSensor( f_sensorInfo *sensor, float xTol, float maxChi2, frack *tr, int threadId );

// TODO: add the remove functions
void f_removeHit(frack *t, int worstHitOffset);
void f_removeWorstHit(frack* t);

bool f_all3SensorsAreDifferent(frack* t);
int f_nbUnused(frack* t, char* usedHits, int event_hit_displ);

void f_addHit ( frack *tr, int offset);
void f_setTrack(frack *tr, int hit0offset, int hit1offset);

void printValuableStatistics();

#endif