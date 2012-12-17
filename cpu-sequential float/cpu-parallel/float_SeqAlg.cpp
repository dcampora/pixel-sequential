
#include "stdafx.h"
#include "float_SeqAlg.h"

// Note: parallel_tracks_vector should have n items allocated already.
void f_TBBSearchByPair::operator() (const blocked_range<int>& r) const{
	for (int i = r.begin(); i != r.end(); ++i){
		f_searchByPair(i, parallel_fracks_vector[i]);
	}
}

void f_searchByPair(int eventRealId, vector<frack>& tracks_vector) {
	int lastSensor  = sens_num - 1;
	int firstSensor = 2;

	frack m_track;

	int eventId = eventRealId % num_events;
	char* usedHits = (char*) calloc(max_hits, sizeof(char));

	// Helping variables
	int event_sensor_displ = eventId * sens_num;
	int event_hit_displ = eventId * max_hits;

	int sens0, sens1, first1, hit0_no, hit1_no;
	f_sensorInfo sensor0, sensor1, extra_sensor;
	float dxMax, dyMax;

	// Iterate from the last until the first+1 (including it)
	for ( sens0 = lastSensor; firstSensor <= sens0; sens0 -= 1 ) {
 	 // sens1 is next sensor in same side
 	 sens1 = sens0 - 2;

 	 sensor0.startPosition = f_sensor_hitStarts[event_sensor_displ + sens0];
 	 sensor0.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens0];
 	 sensor0.z = f_sensor_Zs[event_sensor_displ + sens0];
 	 sensor1.startPosition = f_sensor_hitStarts[event_sensor_displ + sens1];
 	 sensor1.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens1];
 	 sensor1.z = f_sensor_Zs[event_sensor_displ + sens1];

     // Indicator of max dx, dy permissible over next hit
 	 dxMax = f_m_maxXSlope * fabs( sensor1.z - sensor0.z );
 	 dyMax = f_m_maxYSlope * fabs( sensor1.z - sensor0.z );

	 first1 = 0;
 	 for (hit0_no = 0; hit0_no < sensor0.hitsNum; hit0_no++){
 		int hit0_offset = event_hit_displ + sensor0.startPosition + hit0_no;
 		if ( usedHits[hit0_offset - event_hit_displ] ){
			continue;
		}
 		float x0  = f_hit_Xs[hit0_offset];
 		float y0  = f_hit_Ys[hit0_offset];

		// cout << " hit " << f_hit_ids[hit0_offset] << " against:" << endl;

		// Min and max x permissible over next hit
		float xMin = x0 - dxMax;
		float xMax = x0 + dxMax;

		// Iterate hits in s1
		for (hit1_no = first1; hit1_no < sensor1.hitsNum; hit1_no++) {

			int hit1_offset = event_hit_displ + sensor1.startPosition + hit1_no;
			float x1 = f_hit_Xs[hit1_offset];

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
			float y1  = f_hit_Ys[hit1_offset];
			if ( fabs( y1 - y0 ) > dyMax ){
				continue;
			}

			// Creates a frack starting at itH0 (first hit) and containing itH1 (second hit - addHit)
			m_track.hits.clear();
			f_setTrack(&m_track, hit0_offset, hit1_offset);

			// cout << endl << "hit" << f_hit_ids[hit0_offset] << " and hit" << f_hit_ids[hit1_offset] << " are compatible, creating frack" << endl;

			//== Cut on R2Beam if needed : backward tracks, i.e zBeam > first hit
			if ( sensor0.z < f_m_maxZForRBeamCut ) {
				float z_beam  = f_zBeam(&m_track);
				if ( z_beam > sensor0.z ) {
					float r2Beam = f_r2AtZ( z_beam, &m_track );
					if ( r2Beam > f_m_maxR2Beam ){
						continue;
					}
				}
			}

			//== Extend downstream, on both sides of the detector as soon as one hit is missed
			int extraStep = 2;
			int extraSens = sens1-extraStep;

			int nbMissed = 0;

			float lastZ = sensor1.z;

			// Test current frack against tracks placed behind
			while ( extraSens >= 0 ) {
				extra_sensor.startPosition = f_sensor_hitStarts[event_sensor_displ + extraSens];
				extra_sensor.hitsNum	= f_sensor_hitsNums[event_sensor_displ + extraSens];
				extra_sensor.z = f_sensor_Zs[event_sensor_displ + extraSens];

				float tol     = f_m_extraTol;
				float maxChi2 = f_m_maxChi2ToAdd;
				if ( extra_sensor.z < lastZ - 100.0 ) {
					tol     = 2 * tol;
					maxChi2 = 2 * maxChi2;
				}

				bool added = f_addHitsOnSensor(&extra_sensor, tol, maxChi2, &m_track, eventId);
				if ( added ) {
					nbMissed = 0;
					lastZ = extra_sensor.z;
				} else {
					nbMissed += extraStep;
					extraStep = 1;
				}
				if ( f_m_maxMissed < nbMissed ){
					break;
				}

				extraSens -= extraStep;
			}

			//== Try upstream if almost forward tracks
			if ( sensor0.z > f_m_maxZForRBeamCut ) {
				extraStep = 1;
				extraSens = sens0 + 3;  // + 2 already tried...
				nbMissed = 2;

				while ( extraSens <= lastSensor ) {
					extra_sensor.startPosition = f_sensor_hitStarts[event_sensor_displ + extraSens];
					extra_sensor.hitsNum	= f_sensor_hitsNums[event_sensor_displ + extraSens];
					extra_sensor.z = f_sensor_Zs[event_sensor_displ + extraSens];

                    bool added = f_addHitsOnSensor(&extra_sensor, f_m_extraTol, f_m_maxChi2ToAdd, &m_track, eventId);

					if ( added ) {
						nbMissed = 0;
					} else {
						nbMissed += extraStep;
					}
					if ( f_m_maxMissed < nbMissed ){
						break;
					}
					extraSens += extraStep;
				}

			}

			f_removeWorstHit(&m_track);

			if ( m_track.trackHitsNum < 3 ){
				// cout << "Track only has " << m_track.trackHitsNum << " hits!" << endl;
				continue;
			}

			//== Add compatible hits in sens0 and sens1.
			int next_hit = hit0_no + 1 + sensor0.startPosition;
			if ( next_hit != sensor0.startPosition + sensor0.hitsNum ) {
				if ( f_chi2Track(&m_track, event_hit_displ + next_hit) < f_m_maxChi2SameSensor ) {
					hit0_no++;
					f_addHit(&m_track, event_hit_displ + next_hit);
				}
			}

			next_hit = hit1_no + 1 + sensor1.startPosition;
			if ( next_hit != sensor1.startPosition + sensor1.hitsNum ) {
				if ( f_chi2Track(&m_track, event_hit_displ + next_hit) < f_m_maxChi2SameSensor ) {
					hit1_no++;
					f_addHit(&m_track, event_hit_displ + next_hit);
				}
			}

			//== Final check: if only 3 hits, all should be unused and chi2 good.
			if ( m_track.trackHitsNum == 3 ) {
				if ( !f_all3SensorsAreDifferent(&m_track) ) {
					// cout << "Not all three sensors are different" << endl;
					continue;
				}
				if( f_nbUnused(&m_track, usedHits, event_hit_displ) != 3){
					// cout << "There is less than 3 non-used hits" << endl;
					continue;
				}
				if( f_chi2(&m_track) > f_m_maxChi2Short){
					// cout << "Chi2 test not passed" << endl;
					continue;
				}
			} else {
				if ( f_nbUnused(&m_track, usedHits, event_hit_displ) < .6 * m_track.trackHitsNum ) {
					// cout << "More than 60% of the hits are used already" << endl;
					continue;
				}
			}

			// cout << endl << "++ Writing frack" << endl;
			tracks_vector.push_back(m_track);



			// Tag used hits IF the frack is bigger than 3!!
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