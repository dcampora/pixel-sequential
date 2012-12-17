
//#include "stdafx.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include <stdlib.h>
#include <math.h>

#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"

extern int* numHitsEvents;
extern int* sensor_hitStarts;
extern int* sensor_hitsNums;
extern double* sensor_Zs;
extern int* hit_ids;
extern int* track_ids;
extern int* hit_sensorNums;
extern int* hit_isUseds;
extern double* hit_Xs;
extern double* hit_Ys;
extern double* hit_Zs;
extern double* hit_Ws;

extern int num_events;
extern int max_hits;
extern int sens_num;

using namespace std;
using namespace tbb;

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif


typedef struct {
	int startPosition;
	int hitsNum;
	double z;
} sensorInfo;

struct track {
	double m_x0;
	double m_tx;
	double m_y0;
	double m_ty;

	double m_s0;
	double m_sx;
	double m_sz;
	double m_sxz;
	double m_sz2;

	double m_u0;
	double m_uy;
	double m_uz;
	double m_uyz;
	double m_uz2;
	
	int trackHitsNum;

	// int internalId;
	// int firstHit;
	// bool m_backward;

	// for testing purposes
	vector<int> hits;
};

extern vector< vector<track> > parallel_tracks_vector;

class TBBSearchByPair {
public:
	void operator() (const blocked_range<int>& r) const;
};

// void searchByPair_tbb(const blocked_range<int>& r);

double zBeam(track *tr);
double r2AtZ( double z , track *tr);
void solve (track *tr);
double chi2Hit( double x, double y, double hitX, double hitY, double hitW);
double xAtHit(track *tr, double z );
double yAtHit( track *tr, double z  );
double chi2Track(track *tr, int offset);
double chi2(track *t);
bool addHitsOnSensor( sensorInfo sensor, double xTol, double maxChi2,
							 track *tr, int threadId );
void removeHit(track *t, int worstHitOffset);
void removeWorstHit(track* t);
bool all3SensorsAreDifferent(track* t);
int nbUnused(track* t, char* usedHits, int event_hit_displ);

void addHit ( track *tr, int offset);
void setTrack(track *tr, int hit0offset, int hit1offset);

void searchByPair(int event_no, vector<track>& vector_tracks);
