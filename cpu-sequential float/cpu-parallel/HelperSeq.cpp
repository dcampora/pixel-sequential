
#include "stdafx.h"
#include "HelperSeq.h"

float f_zBeam(frack *tr) {
	return -( tr->m_x0 * tr->m_tx + tr->m_y0 * tr->m_ty ) / ( tr->m_tx * tr->m_tx + tr->m_ty * tr->m_ty );
}

float f_r2AtZ( float z , frack *tr) {
    float xx = tr->m_x0 + z * tr->m_tx;
    float yy = tr->m_y0 + z * tr->m_ty;
    return xx*xx + yy * yy;
 }

void f_solve (frack *tr) {
	float den = ( tr->m_sz2 * tr->m_s0 - tr->m_sz * tr->m_sz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_tx     = ( tr->m_sxz * tr->m_s0  - tr->m_sx  * tr->m_sz ) / den;
	tr->m_x0     = ( tr->m_sx  * tr->m_sz2 - tr->m_sxz * tr->m_sz ) / den;

	den = ( tr->m_uz2 * tr->m_u0 - tr->m_uz * tr->m_uz );
	if ( fabs(den) < 10e-10 ) den = 1.;
	tr->m_ty     = ( tr->m_uyz * tr->m_u0  - tr->m_uy  * tr->m_uz ) / den;
	tr->m_y0     = ( tr->m_uy  * tr->m_uz2 - tr->m_uyz * tr->m_uz ) / den;
}

void f_addHit ( frack *tr, int offset) {

	// track_ids[offset] = tr->internalId;
	tr->trackHitsNum++;

	float z = f_hit_Zs[offset];
	float x = f_hit_Xs[offset];
	// float f_w = f_hit_Ws[offset];

	float wz = f_w * z;

	tr->m_s0  += f_w;
	tr->m_sx  += f_w * x;
	tr->m_sz  += wz;
	tr->m_sxz += wz * x;
	tr->m_sz2 += wz * z;

	float y = f_hit_Ys[offset];

	tr->m_u0  += f_w;
	tr->m_uy  += f_w * y;
	tr->m_uz  += wz;
	tr->m_uyz += wz * y;
	tr->m_uz2 += wz * z;

	if(tr->trackHitsNum > 1) f_solve(tr);

	tr->hits.push_back(offset);
}

void f_setTrack(frack *tr, int hit0_offset, int hit1_offset){

	// track_ids[hit0_offset] = tr->internalId;
	tr->trackHitsNum = 1;

	float z = f_hit_Zs[hit0_offset];
	float x = f_hit_Xs[hit0_offset];
	// float f_w = f_hit_Ws[hit0_offset];

	float wz = f_w * z;

	tr->m_s0  = f_w;
	tr->m_sx  = f_w * x;
	tr->m_sz  = wz;
	tr->m_sxz = wz * x;
	tr->m_sz2 = wz * z;

	float y = f_hit_Ys[hit0_offset];

	tr->m_u0  = f_w;
	tr->m_uy  = f_w * y;
	tr->m_uz  = wz;
	tr->m_uyz = wz * y;
	tr->m_uz2 = wz * z;

	// TODO: Remove when not needed
	tr->hits.push_back(hit0_offset);

	f_addHit (tr, hit1_offset);
}

float f_chi2Hit( float x, float y, float hitX, float hitY, float hitW){
	float dx = x - hitX;
	float dy = y - hitY;
	return dx * dx * (hitW) + dy * dy * (hitW);
}
float f_xAtHit(frack *tr, float z )
{
	return tr->m_x0 + tr->m_tx * z;
}
float f_yAtHit( frack *tr, float z  )
{
	return tr->m_y0 + tr->m_ty * z;
}
float f_chi2Track(frack *tr, int offset)
{
	float z = f_hit_Zs[offset];
	return f_chi2Hit( f_xAtHit( tr, z ), f_yAtHit(tr, z ), f_hit_Xs[offset], f_hit_Ys[offset], f_w);
}
float f_chi2(frack *t)
{
	float ch = 0.0;
	int nDoF  = -4;
	int hitNumber;
	for (int i=0; i<t->hits.size(); i++){
		hitNumber = t->hits[i];
		ch += f_chi2Track(t, hitNumber);
		nDoF += 2;
	}
	return ch/nDoF;
}

bool f_addHitsOnSensor( f_sensorInfo *sensor, float xTol, float maxChi2,
							 frack *tr, int eventId ) {
	
	if (sensor->hitsNum == 0) return false;
	int offset = eventId * max_hits;

	float xGuess = f_xAtHit(tr, sensor->z) - xTol - 1;
	int lastHit = sensor->startPosition + sensor->hitsNum - 1;
	if(f_hit_Xs[offset + lastHit] < xGuess) return false;

	int hitStart = sensor->startPosition;
	unsigned int step = sensor->hitsNum;
	while ( step > 2 ) {
		step = step/2;
		if (f_hit_Xs[offset + hitStart + step] < xGuess) hitStart += step;
	}

	bool added = false;
	int tmpOffset = 0;
	float xPred;
	for(int iH=hitStart; iH<=lastHit; ++iH){
		tmpOffset = offset + iH;
		xPred = f_xAtHit(tr, f_hit_Zs[tmpOffset]);
		if ( f_hit_Xs[tmpOffset] + xTol < xPred ) continue;
		if ( f_hit_Xs[tmpOffset] - xTol > xPred ) break;
		if ( f_chi2Track(tr, tmpOffset) < maxChi2 ) {
			f_addHit(tr, tmpOffset);
			// *usedHit = tmpOffset; - Used hits are tagged by the end of the algorithm, not before.
			added = true;
		}
	}
	return added;
}

void f_removeHit(frack *tr, int worstHitOffset){
	tr->trackHitsNum--;

	float z = f_hit_Zs[worstHitOffset];
	// float f_w = f_hit_Ws[worstHitOffset];
	float x = f_hit_Xs[worstHitOffset];
	
	float wz = f_w * z;

	tr->m_s0  -= f_w;
	tr->m_sx  -= f_w * x;
	tr->m_sz  -= wz;
	tr->m_sxz -= wz * x;
	tr->m_sz2 -= wz * z;

	float y = f_hit_Ys[worstHitOffset];

	tr->m_u0  -= f_w;
	tr->m_uy  -= f_w * y;
	tr->m_uz  -= wz;
	tr->m_uyz -= wz * y;
	tr->m_uz2 -= wz * z;

	vector<int>::iterator it = find(tr->hits.begin(), tr->hits.end(), worstHitOffset);

	tr->hits.erase(it);

	if( tr->trackHitsNum > 1 ) f_solve(tr);
}

//== Remove the worst hit until all chi2 are good
void f_removeWorstHit(frack* tr)
{
	float topChi2 = 1.e9;
	int worstHitOffset;

	while( topChi2 > f_m_maxChi2PerHit ) {
	    topChi2 = 0.0;


	    // This for loop gets the worst hit
		for (int i=0; i<tr->hits.size(); i++){

			float myChi2 = f_chi2Track(tr, tr->hits[i]);
			if (myChi2 > topChi2){
				topChi2 = myChi2;
				worstHitOffset = tr->hits[i];
			}
		}

	    // If it's bad, we remove it
	    if ( topChi2 > f_m_maxChi2PerHit ) {
	      // f_hit_isUseds[worstHitOffset] = 0;
		  // It has still not been added to isUseds, no need to do this :)

	      f_removeHit(tr, worstHitOffset);
		  // This changes the chi2 of the frack, which is why 
	    }

	    // And the algorithm goes on... ?
	    // Every hit with chi2 > maxChi2 will be removed... is this the desired behaviour?
		// -> yes, read description above
	}
}

bool f_all3SensorsAreDifferent(frack *t) {
    float s0 = f_hit_sensorNums[t->hits[0]];
    float s1 = f_hit_sensorNums[t->hits[1]];
    float s2 = f_hit_sensorNums[t->hits[2]];
	
    if ( s0 == s1 ) return false;
    if ( s0 == s2 ) return false;
    if ( s1 == s2 ) return false;
    return true;
}

int f_nbUnused(frack *t, char* usedHits, int event_hit_displ) {
	int nn = 0;
	for (vector<int>::iterator it = t->hits.begin(); it != t->hits.end(); ++it){
		if (!usedHits[(*it) - event_hit_displ])
			++nn;
	}
	return nn;
}


void printValuableStatistics(){
	vector<float> means_hits, means_two_hits, means_three_hits;

	for(int event_no = 0; event_no < num_events; event_no++){

		int sens0, sens1, sens2;
		int i, j, k;
		f_sensorInfo sensor0, sensor1, sensor2;

		int event_sensor_displ = event_no * sens_num;
		int event_hit_displ = event_no * max_hits;

		vector<float> hits;
		for(i=0; i<sens_num; i++)
			hits.push_back(f_sensor_hitsNums[event_sensor_displ + i]);
	
		vector<float> hits_two_cardinal;
		vector<float> hits_three_cardinal;

		// Iterate in pairs
		for(sens0 = sens_num - 1; sens0 >= 2; sens0 = sens0 - 1){
			sens1 = sens0 - 2;
		
			sensor0.startPosition = f_sensor_hitStarts[event_sensor_displ + sens0];
			sensor0.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens0];
			sensor0.z = f_sensor_Zs[event_sensor_displ + sens0];
			sensor1.startPosition = f_sensor_hitStarts[event_sensor_displ + sens1];
			sensor1.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens1];
			sensor1.z = f_sensor_Zs[event_sensor_displ + sens1];
		
			hits_two_cardinal.push_back((sensor0.hitsNum * sensor1.hitsNum));
		}
	
		// Iterate in trios
		for(sens0 = sens_num - 2; sens0 >= 4; sens0 = sens0 - 1){
			sens1 = sens0 - 2;
			sens2 = sens1 - 2;
		
			sensor0.startPosition = f_sensor_hitStarts[event_sensor_displ + sens0];
			sensor0.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens0];
			sensor0.z = f_sensor_Zs[event_sensor_displ + sens0];
			sensor1.startPosition = f_sensor_hitStarts[event_sensor_displ + sens1];
			sensor1.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens1];
			sensor1.z = f_sensor_Zs[event_sensor_displ + sens1];
			sensor2.startPosition = f_sensor_hitStarts[event_sensor_displ + sens2];
			sensor2.hitsNum	= f_sensor_hitsNums[event_sensor_displ + sens2];
			sensor2.z = f_sensor_Zs[event_sensor_displ + sens2];
		
			hits_three_cardinal.push_back((sensor2.hitsNum * sensor1.hitsNum + sensor2.hitsNum * sensor0.hitsNum));
		}
	
		map<string, float> results1 = calcResults(hits);
		map<string, float> results2 = calcResults(hits_two_cardinal);
		map<string, float> results3 = calcResults(hits_three_cardinal);

		means_hits.push_back(results1["mean"]);
		means_two_hits.push_back(results2["mean"]);
		means_three_hits.push_back(results3["mean"]);

		cout << "event " << event_no << ":" << endl
			 << "hits:" << endl;
		printStatistics(results1, string(" "));

		cout << "pairs:" << endl;
		printStatistics(results2, string(" "));

		cout << "trios:" << endl;
		printStatistics(results3, string(" "));

		cout << endl;
	}

	map<string, float> results1 = calcResults(means_hits);
	map<string, float> results2 = calcResults(means_two_hits);
	map<string, float> results3 = calcResults(means_three_hits);

	cout << endl << "Statistics for all events:" << endl
		 << "Means of hits:" << endl;
	printStatistics(results1, string(" "));
	
	cout << endl;
	printStatistics(results2, string(" "));
	
	cout << endl;
	printStatistics(results3, string(" "));

	cout << endl << "Have a nice day! :)" << endl;
}