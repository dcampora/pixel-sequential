
// #include "stdafx.h"
#include "SeqAlg.h"

double 	m_maxXSlope			= 0.400;
double 	m_maxYSlope			= 0.300;
double 	m_maxZForRBeamCut	= 200.0;
double 	m_maxR2Beam			= 1.0;
int 	m_maxMissed			= 4;
double 	m_extraTol			= 0.150;
double 	m_maxChi2ToAdd		= 100.0;
double 	m_maxChi2SameSensor	= 16.0;
double  m_maxChi2Short		= 6.0 ;
double  m_maxChi2PerHit		= 16.0;
int 	m_sensNum			= 48;

double  w                   = 0.050 / sqrt( 12. );

double zBeam(track *tr) {
	return -( tr->m_x0 * tr->m_tx + tr->m_y0 * tr->m_ty ) / ( tr->m_tx * tr->m_tx + tr->m_ty * tr->m_ty );
}

double r2AtZ( double z , track *tr) {
    double xx = tr->m_x0 + z * tr->m_tx;
    double yy = tr->m_y0 + z * tr->m_ty;
    return xx*xx + yy * yy;
 }

void solve (track *tr) {
	double den = ( tr->m_sz2 * tr->m_s0 - tr->m_sz * tr->m_sz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_tx     = ( tr->m_sxz * tr->m_s0  - tr->m_sx  * tr->m_sz ) / den;
	tr->m_x0     = ( tr->m_sx  * tr->m_sz2 - tr->m_sxz * tr->m_sz ) / den;

	den = ( tr->m_uz2 * tr->m_u0 - tr->m_uz * tr->m_uz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_ty     = ( tr->m_uyz * tr->m_u0  - tr->m_uy  * tr->m_uz ) / den;
	tr->m_y0     = ( tr->m_uy  * tr->m_uz2 - tr->m_uyz * tr->m_uz ) / den;
}

void addHit ( track *tr, int offset) {

	// track_ids[offset] = tr->internalId;
	tr->trackHitsNum++;

	double z = hit_Zs[offset];
	double x = hit_Xs[offset];
	// double w = hit_Ws[offset];

	double wz = w * z;

	tr->m_s0  += w;
	tr->m_sx  += w * x;
	tr->m_sz  += wz;
	tr->m_sxz += wz * x;
	tr->m_sz2 += wz * z;

	double y = hit_Ys[offset];

	tr->m_u0  += w;
	tr->m_uy  += w * y;
	tr->m_uz  += wz;
	tr->m_uyz += wz * y;
	tr->m_uz2 += wz * z;

	if( tr->trackHitsNum > 1 ) solve(tr);

	tr->hits.push_back(offset);
}

void setTrack(track *tr, int hit0_offset, int hit1_offset){

	// track_ids[hit0_offset] = tr->internalId;
	tr->trackHitsNum = 1;

	double z = hit_Zs[hit0_offset];
	double x = hit_Xs[hit0_offset];
	// double w = hit_Ws[hit0_offset];

	double wz = w * z;

	tr->m_s0  = w;
	tr->m_sx  = w * x;
	tr->m_sz  = wz;
	tr->m_sxz = wz * x;
	tr->m_sz2 = wz * z;

	double y = hit_Ys[hit0_offset];

	tr->m_u0  = w;
	tr->m_uy  = w * y;
	tr->m_uz  = wz;
	tr->m_uyz = wz * y;
	tr->m_uz2 = wz * z;

	// TODO: Remove when not needed
	tr->hits.push_back(hit0_offset);

	addHit (tr, hit1_offset);
}

double chi2Hit( double x, double y, double hitX, double hitY, double hitW){
	double dx = x - hitX;
	double dy = y - hitY;
	return dx * dx * (hitW) + dy * dy * (hitW);
}
double xAtHit(track *tr, double z )
{
	return tr->m_x0 + tr->m_tx * z;
}
double yAtHit( track *tr, double z  )
{
	return tr->m_y0 + tr->m_ty * z;
}
double chi2Track(track *tr, int offset)
{
	double z = hit_Zs[offset];
	return chi2Hit( xAtHit( tr, z ), yAtHit(tr, z ), hit_Xs[offset], hit_Ys[offset], w);
}
double chi2(track *t)
{
	double ch = 0.0;
	int nDoF  = -4;
	int hitNumber;
	for (int i=0; i<t->hits.size(); i++){
		hitNumber = t->hits[i];
		ch += chi2Track(t, hitNumber);
		nDoF += 2;
	}
	return ch/nDoF;
}

bool addHitsOnSensor( sensorInfo *sensor, double xTol, double maxChi2,
							 track *tr, int eventId ) {
	
	if (sensor->hitsNum == 0) return false;
	int offset = eventId * max_hits;

	double xGuess = xAtHit(tr, sensor->z) - xTol - 1;
	int lastHit = sensor->startPosition + sensor->hitsNum - 1;
	if(hit_Xs[offset + lastHit] < xGuess) return false;

	int hitStart = sensor->startPosition;
	unsigned int step = sensor->hitsNum;
	while ( step > 2 ) {
		step = step/2;
		if (hit_Xs[offset + hitStart + step] < xGuess) hitStart += step;
	}

	bool added = false;
	int tmpOffset = 0;
	double xPred;
	for(int iH=hitStart; iH<=lastHit; ++iH){
		tmpOffset = offset + iH;
		xPred = xAtHit(tr, hit_Zs[tmpOffset]);
		if ( hit_Xs[tmpOffset] + xTol < xPred ) continue;
		if ( hit_Xs[tmpOffset] - xTol > xPred ) break;
		if ( chi2Track(tr, tmpOffset) < maxChi2 ) {
			addHit(tr, tmpOffset);
			// *usedHit = tmpOffset; - Used hits are tagged by the end of the algorithm, not before.
			added = true;
		}
	}
	return added;
}

void removeHit(track *tr, int worstHitOffset){
	tr->trackHitsNum--;

	double z = hit_Zs[worstHitOffset];
	// double w = hit_Ws[worstHitOffset];
	double x = hit_Xs[worstHitOffset];
	
	double wz = w * z;

	tr->m_s0  -= w;
	tr->m_sx  -= w * x;
	tr->m_sz  -= wz;
	tr->m_sxz -= wz * x;
	tr->m_sz2 -= wz * z;

	double y = hit_Ys[worstHitOffset];

	tr->m_u0  -= w;
	tr->m_uy  -= w * y;
	tr->m_uz  -= wz;
	tr->m_uyz -= wz * y;
	tr->m_uz2 -= wz * z;

	vector<int>::iterator it = find(tr->hits.begin(), tr->hits.end(), worstHitOffset);

	tr->hits.erase(it);

	if( tr->trackHitsNum > 1 ) solve(tr);
}

//== Remove the worst hit until all chi2 are good
void removeWorstHit(track* tr)
{
	double topChi2 = 1.e9;
	int worstHitOffset;
	while( topChi2 > m_maxChi2PerHit ) {
	    topChi2 = 0.0;

	    // This for loop gets the worst hit
		for (int i=0; i<tr->hits.size(); i++){
			double myChi2 = chi2Track(tr, tr->hits[i]);
			if (myChi2 > topChi2){
				topChi2 = myChi2;
				worstHitOffset = tr->hits[i];
			}
		}

	    // If it's bad, we remove it
	    if ( topChi2 > m_maxChi2PerHit ) {
	      // usedHits[worstHitOffset] = 0;
		  // It has still not been added to isUseds, no need to do this :)

	      removeHit(tr, worstHitOffset);
		  // This changes the chi2 of the track, which is why 
	    }

	    // And the algorithm goes on... ?
	    // Every hit with chi2 > maxChi2 will be removed... is this the desired behaviour?
		// -> yes, read description above
	}
}

bool all3SensorsAreDifferent(track *t) {
    double s0 = hit_sensorNums[t->hits[0]];
    double s1 = hit_sensorNums[t->hits[1]];
    double s2 = hit_sensorNums[t->hits[2]];
	
    if ( s0 == s1 ) return false;
    if ( s0 == s2 ) return false;
    if ( s1 == s2 ) return false;
    return true;
}

int nbUnused(track *t, char* usedHits, int event_hit_displ) {
	int nn = 0;
	for (vector<int>::iterator it = t->hits.begin(); it != t->hits.end(); ++it){
		if (!usedHits[(*it) - event_hit_displ])
			++nn;
	}
	return nn;
}

// Note: parallel_tracks_vector should have n items allocated already.
void TBBSearchByPair::operator() (const blocked_range<int>& r) const{
	for (int i = r.begin(); i != r.end(); ++i){
		searchByPair(i, parallel_tracks_vector[i]);
	}
}

void searchByPair(int eventRealId, vector<track>& tracks_vector) {
	int lastSensor  = sens_num - 1;
	int firstSensor = 2;

	int eventId = eventRealId % num_events;
	char* usedHits = (char*) calloc(max_hits, sizeof(char));

	track m_track;

	// Helping variables
	int event_sensor_displ = eventId * sens_num;
	int event_hit_displ = eventId * max_hits;

	int sens0, sens1, first1, hit0_no, hit1_no;
	sensorInfo sensor0, sensor1, extra_sensor;
	double dxMax, dyMax;

	// Iterate from the last until the first+1 (including it)
	for ( sens0 = lastSensor; firstSensor <= sens0; sens0 -= 1 ) {
 	 // sens1 is next sensor in same side
 	 sens1 = sens0 - 2;

 	 sensor0.startPosition = sensor_hitStarts[event_sensor_displ + sens0];
 	 sensor0.hitsNum	= sensor_hitsNums[event_sensor_displ + sens0];
 	 sensor0.z = sensor_Zs[event_sensor_displ + sens0];
 	 sensor1.startPosition = sensor_hitStarts[event_sensor_displ + sens1];
 	 sensor1.hitsNum	= sensor_hitsNums[event_sensor_displ + sens1];
 	 sensor1.z = sensor_Zs[event_sensor_displ + sens1];

     // Indicator of max dx, dy permissible over next hit
 	 dxMax = m_maxXSlope * fabs( sensor1.z - sensor0.z );
 	 dyMax = m_maxYSlope * fabs( sensor1.z - sensor0.z );

	 first1 = 0;
 	 for (hit0_no = 0; hit0_no < sensor0.hitsNum; hit0_no++){
 		int hit0_offset = event_hit_displ + sensor0.startPosition + hit0_no;
 		if ( usedHits[hit0_offset - event_hit_displ] ){
			continue;
		}
 		double x0  = hit_Xs[hit0_offset];
 		double y0  = hit_Ys[hit0_offset];

		// cout << " hit " << hit_ids[hit0_offset] << " against:" << endl;

		// Min and max x permissible over next hit
		double xMin = x0 - dxMax;
		double xMax = x0 + dxMax;

		// Iterate hits in s1
		for (hit1_no = first1; hit1_no < sensor1.hitsNum; hit1_no++) {

			int hit1_offset = event_hit_displ + sensor1.startPosition + hit1_no;
			double x1 = hit_Xs[hit1_offset];

			// Require second hit not to be used, and to be within x limits
			if ( x1 < xMin ) {
				first1 = hit1_no + 1; // Start on the item (hit1_no+1) for next iteration
				continue;
			}
			if ( x1 > xMax ){
				break; // hits are ordered by x (clever)
			}
			if ( usedHits[hit1_offset - event_hit_displ] ){
				continue;
			}

			// Check hit is within y limits
			double y1  = hit_Ys[hit1_offset];
			if ( fabs( y1 - y0 ) > dyMax ){
				continue;
			}

			// Creates a track starting at itH0 (first hit) and containing itH1 (second hit - addHit)
			m_track.hits.clear();
			setTrack(&m_track, hit0_offset, hit1_offset);

			// cout << endl << "hit" << hit_ids[hit0_offset] << " and hit" << hit_ids[hit1_offset] << " are compatible, creating track" << endl;

			//== Cut on R2Beam if needed : backward tracks, i.e zBeam > first hit
			if ( sensor0.z < m_maxZForRBeamCut ) {
				double z_beam  = zBeam(&m_track);
				if ( z_beam > sensor0.z ) {
					double r2Beam = r2AtZ( z_beam, &m_track );
					if ( r2Beam > m_maxR2Beam ){
						continue;
					}
				}
			}

			//== Extend downstream, on both sides of the detector as soon as one hit is missed
			int extraStep = 2;
			int extraSens = sens1-extraStep;

			int nbMissed = 0;

			double lastZ = sensor1.z;

			while ( extraSens >= 0 ) {
				extra_sensor.startPosition = sensor_hitStarts[event_sensor_displ + extraSens];
				extra_sensor.hitsNum	= sensor_hitsNums[event_sensor_displ + extraSens];
				extra_sensor.z = sensor_Zs[event_sensor_displ + extraSens];

				double tol     = m_extraTol;
				double maxChi2 = m_maxChi2ToAdd;
				if ( extra_sensor.z < lastZ - 100.0 ) {
					tol     = 2 * tol;
					maxChi2 = 2 * maxChi2;
				}

				bool added = addHitsOnSensor(&extra_sensor, tol, maxChi2, &m_track, eventId);
				if ( added ) {
					nbMissed = 0;
					lastZ = extra_sensor.z;
				} else {
					nbMissed += extraStep;
					extraStep = 1;
				}
				if ( m_maxMissed < nbMissed ){
					break;
				}

				extraSens -= extraStep;
			}

			//== Try upstream if almost forward tracks
			if ( sensor0.z > m_maxZForRBeamCut ) {
				extraStep = 1;
				extraSens = sens0 + 3;  // + 2 already tried...
				nbMissed = 2;

				while ( extraSens <= lastSensor ) {
					extra_sensor.startPosition = sensor_hitStarts[event_sensor_displ + extraSens];
					extra_sensor.hitsNum	= sensor_hitsNums[event_sensor_displ + extraSens];
					extra_sensor.z = sensor_Zs[event_sensor_displ + extraSens];

                    bool added = addHitsOnSensor(&extra_sensor, m_extraTol, m_maxChi2ToAdd, &m_track, eventId);

					if ( added ) {
						nbMissed = 0;
					} else {
						nbMissed += extraStep;
					}
					if ( m_maxMissed < nbMissed ){
						break;
					}
					extraSens += extraStep;
				}

			}

			removeWorstHit(&m_track);

			if ( m_track.trackHitsNum < 3 ){
				// cout << "Track only has " << m_track.trackHitsNum << " hits!" << endl;
				continue;
			}

			//== Add compatible hits in sens0 and sens1.
			int next_hit = hit0_no + 1 + sensor0.startPosition;
			if ( next_hit != sensor0.startPosition + sensor0.hitsNum ) {
				if ( chi2Track(&m_track, event_hit_displ + next_hit) < m_maxChi2SameSensor ) {
					hit0_no++;
					addHit(&m_track, event_hit_displ + next_hit);
				}
			}

			next_hit = hit1_no + 1 + sensor1.startPosition;
			if ( next_hit != sensor1.startPosition + sensor1.hitsNum ) {
				if ( chi2Track(&m_track, event_hit_displ + next_hit) < m_maxChi2SameSensor ) {
					hit1_no++;
					addHit(&m_track, event_hit_displ + next_hit);
				}
			}

			//== Final check: if only 3 hits, all should be unused and chi2 good.
			if ( m_track.trackHitsNum == 3 ) {
				if ( !all3SensorsAreDifferent(&m_track) ) {
					// cout << "Not all three sensors are different" << endl;
					continue;
				}
				if( nbUnused(&m_track, usedHits, event_hit_displ) != 3){
					// cout << "There is less than 3 non-used hits" << endl;
					continue;
				}
				if( chi2(&m_track) > m_maxChi2Short){
					// cout << "Chi2 test not passed" << endl;
					continue;
				}
			} else {
				if ( nbUnused(&m_track, usedHits, event_hit_displ) < .6 * m_track.trackHitsNum ) {
					// cout << "More than 60% of the hits are used already" << endl;
					continue;
				}
			}

			// cout << endl << "++ Writing track" << endl;
			tracks_vector.push_back(m_track);

			// Tag used hits IF the track is bigger than 3!!
			if (m_track.trackHitsNum > 3){
				for (int i=0; i<m_track.hits.size(); i++){
					usedHits[m_track.hits[i] - event_hit_displ] = 1;
				}
				break;
			}
		} // itH1


 	 } // itH0
   } // sens0

   free(usedHits);
 }

