// tests.cpp : Defines the entry point for the console application.
//

// #include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>

#include <cmath>

#include "tbb/tick_count.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include "floaties.h"

#include "SeqAlg.h"
#include "float_SeqAlg.h"

using namespace std;
using namespace tbb;

template<class T>
void printResultTracks(vector<T> tracks, int event_no, string track_folder);

bool areTracksEqual(vector<track> tracks_1, vector<track> tracks_2);

template<class T, class Z>
void compareTrackVectors(vector<vector<T> > correct_tracks_vector,
	vector<vector<Z> > new_tracks_vector);

template<class T, class Z>
int compareTracks(vector<T> tracks_1, vector<Z> tracks_2);


void readFile(string filename, char*& input, bool do_malloc);
map<string, float> calcResults(vector<float> times);
float infinite();

template<class T>
void clean_vector(vector<vector<T> >& _tracks_vector);

// Prepare the data types (globally accesible)
int* numHitsEvents;
int* sensor_hitStarts;
int* sensor_hitsNums;
double* sensor_Zs;
int* hit_ids;
int* track_ids;
int* hit_sensorNums;
int* hit_isUseds;
double* hit_Xs;
double* hit_Ys;
double* hit_Zs;
double* hit_Ws;

int* f_numHitsEvents;
int* f_sensor_hitStarts;
int* f_sensor_hitsNums;
float* f_sensor_Zs;
int* f_hit_ids;
int* f_track_ids;
int* f_hit_sensorNums;
int* f_hit_isUseds;
float* f_hit_Xs;
float* f_hit_Ys;
float* f_hit_Zs;
float* f_hit_Ws;

vector<vector<track> > parallel_tracks_vector;
vector<vector<frack> > parallel_fracks_vector;
vector<vector<track> > sequential_tracks_vector;

int num_events, sens_num, max_hits;
int num_events_to_process;

template <class T>
static std::string toString(T t){
	std::stringstream ss;
	std::string s;
	ss << t;
	ss >> s;
	return s;
}

int main(int argc, char* argv[])
{
	char c1, c2, *optarg;
	num_events = 200;
	sens_num = 48;
	max_hits = 3792;

	num_events_to_process = 200;
	int repeat_n_times = 1;
	int max_threads = 1;
	int compute = 7;

	int last_opt = 0;
	while ((argc - last_opt) > 1){
		c1 = argv[last_opt+1][0];
		c2 = argv[last_opt+1][1];
		optarg = argv[last_opt+2];
		last_opt += 2;

		/* if(c2)
			printf(" %c%c %s", c1, c2, optarg);
		else
			printf(" %c %s", c1, optarg); */

		switch(c1){
		case '-':
			switch(c2){
			case 'n':
				num_events_to_process = atoi(optarg);
				break;
			case 'r':
				repeat_n_times = atoi(optarg);
				break;
			case 't':
				max_threads = atoi(optarg);
				break;
			case 'c':
				compute = atoi(optarg);
				break;
			case 'h':
				cout << "Use:" << endl
					 << " cpu-parallel -n <num_events> -r <repeat_times> -t <max_threads> -c <configuration>" << endl;
				return 0;
			}
            break;
		}
	}

	cout << "Setup:" << endl
		 << " Number of events to process: " << num_events_to_process << endl
		 << " Max number of threads: " << max_threads << endl
		 << " Number of repetitions per experiment: " << repeat_n_times << endl << endl;
	
	// Give me them datas!!11!
	string c = "input.dump";
	char* input;

	// Read file
	readFile(c, input, true);

	// Pointers, pointers, pointers!1!!eleven!
	int dataOffset = max_hits * num_events;

	numHitsEvents = (int*) &input[0];
	sensor_hitStarts = (int*) (numHitsEvents + num_events);
	sensor_hitsNums = (int*) (sensor_hitStarts + num_events * sens_num);
	sensor_Zs = (double*) (sensor_hitsNums + num_events * sens_num);
	
	hit_ids = (int*) (sensor_Zs + num_events * sens_num);
	track_ids = hit_ids +  dataOffset;
	hit_sensorNums = track_ids + dataOffset;
	hit_isUseds = hit_sensorNums + dataOffset;
	hit_Xs = (double*) (hit_isUseds + dataOffset);
	hit_Ys = hit_Xs + dataOffset;
	hit_Zs = hit_Ys + dataOffset;
	hit_Ws = hit_Zs + dataOffset;

	
	// Give me them datas!!11!
	string c_float = "input_float.dump";
	char* input_float;

	// Read file
	readFile(c_float, input_float, true);

	// Pointers, pointers, pointers!1!!eleven!
	f_numHitsEvents = (int*) &input_float[0];
	f_sensor_hitStarts = (int*) (f_numHitsEvents + num_events);
	f_sensor_hitsNums = (int*) (f_sensor_hitStarts + num_events * sens_num);
	f_sensor_Zs = (float*) (f_sensor_hitsNums + num_events * sens_num);
	
	f_hit_ids = (int*) (f_sensor_Zs + num_events * sens_num);
	f_track_ids = f_hit_ids +  dataOffset;
	f_hit_sensorNums = f_track_ids + dataOffset;
	f_hit_isUseds = f_hit_sensorNums + dataOffset;
	f_hit_Xs = (float*) (f_hit_isUseds + dataOffset);
	f_hit_Ys = f_hit_Xs + dataOffset;
	f_hit_Zs = f_hit_Ys + dataOffset;
	f_hit_Ws = f_hit_Zs + dataOffset;


	vector<float> times;

	// Allocate space for tracks vectors
	parallel_tracks_vector.clear();
	parallel_fracks_vector.clear();
	sequential_tracks_vector.clear();
	for (int i=0; i<num_events_to_process; ++i){
		vector<track> t1, t2;
		vector<frack> t3;

		parallel_tracks_vector.push_back(t1);
		sequential_tracks_vector.push_back(t2);
		parallel_fracks_vector.push_back(t3);
	}

	// Prints awesome statistics of the data
	// printValuableStatistics();


	// Start sequential algorithm
	cout << "Starting searchByPair..." << endl;
	for (int t=0; t<repeat_n_times; t++){
		if(t!=0) readFile(c, input, false);

		clean_vector(sequential_tracks_vector);

		tick_count tc_start = tick_count::now();
		for (int i=1; i<=num_events_to_process; i++){
			searchByPair(i-1, sequential_tracks_vector[i-1]);
		}
		tick_count tc_end = tick_count::now();

		times.push_back((tc_end - tc_start).seconds());
	}
	map<string, float> sequential_results = calcResults(times);
	
	map<string, float> f_sequential_results;
	vector<map<string, float> > parallel_experiments;
	vector<map<string, float> > f_parallel_experiments;

	if((compute >> 1) & 0x01){
		// Start sequential float-algorithm
		cout << "Starting f_searchByPair..." << endl;
		for (int t=0; t<repeat_n_times; t++){
		
			readFile(c_float, input_float, false);
			clean_vector(parallel_fracks_vector);

			tick_count tc_start = tick_count::now();
			for (int i=1; i<=num_events_to_process; i++){
				f_searchByPair(i-1, parallel_fracks_vector[i-1]);
			}
			tick_count tc_end = tick_count::now();

			times.push_back((tc_end - tc_start).seconds());
		}
		f_sequential_results = calcResults(times);
	}


	if((compute >> 2) & 0x01){
		cout << "Starting parallel (by events) f_searchByPair..." << endl;

		// Perform parallel experiments!
		for(int i=1; i<=max_threads; i=i*2){
			cout << " Performing " << i << "-thread experiment" << endl;
		
			task_scheduler_init init(i);

			// Perform experiment!
			times.clear();
			for (int t=0; t<repeat_n_times; t++){
				// Re-read file
				readFile(c_float, input_float, false);
				clean_vector(parallel_fracks_vector);

				tick_count parallel_start = tick_count::now();
				parallel_for(blocked_range<int>(0, num_events_to_process),
							 f_TBBSearchByPair());
				tick_count parallel_end = tick_count::now();

				times.push_back((parallel_end - parallel_start).seconds());
			}

			// Store
			f_parallel_experiments.push_back(calcResults(times));
		}
		// */
	}


	if(compute & 0x01){
		cout << "Starting parallel (by events) searchByPair..." << endl;

		// Perform parallel experiments!
		for(int i=1; i<=max_threads; i=i*2){
			cout << " Performing " << i << "-thread experiment" << endl;
		
			task_scheduler_init init(i);

			// Perform experiment!
			times.clear();
			for (int t=0; t<repeat_n_times; t++){
				// Re-read file
				readFile(c, input, false);

				clean_vector(parallel_tracks_vector);

				tick_count parallel_start = tick_count::now();
				parallel_for(blocked_range<int>(0, num_events_to_process),
							 TBBSearchByPair());
				tick_count parallel_end = tick_count::now();

				times.push_back((parallel_end - parallel_start).seconds());
			}

			// Store
			parallel_experiments.push_back(calcResults(times));
		}
	}

	
	/*
	compareTrackVectors(sequential_tracks_vector, parallel_fracks_vector);
	
	cout << "Writing results..." << endl;
	for (int i=0; i<num_events_to_process; i++)
		printResultTracks(parallel_fracks_vector[i], i, "tracks");

	*/

	
	// compareTrackVectors(sequential_tracks_vector, parallel_fracks_vector);
	
	/*
	cout << "Writing results..." << endl;
	for (int i=0; i<num_events_to_process; i++)
		printResultTracks(parallel_fracks_vector[i], i, "tracks"); */

	cout << std::setprecision(3);

	cout << "Algorithm execution finished!" << endl
		 << "Sequential: " << sequential_results["mean"] << " seconds (mean), "
		 << sequential_results["deviation"] << " deviation, "
		 << sequential_results["variance"] << " variance, "
		 << sequential_results["min"] << " min, "
		 << sequential_results["max"] << " max" << endl;
	
	if((compute >> 1) & 0x01){
		cout << "f_Sequential: " << f_sequential_results["mean"] << " seconds (mean), "
			 << f_sequential_results["deviation"] << " deviation, "
			 << f_sequential_results["variance"] << " variance, "
			 << f_sequential_results["min"] << " min, "
			 << f_sequential_results["max"] << " max, "
			 << sequential_results["mean"] / f_sequential_results["mean"] << "x" << endl;
	}

	if((compute >> 2) & 0x01){
		cout << "f_Parallel by event: " << endl;
		for (int i=1, j=0; i<=max_threads; i=i*2, j++){
			cout << i << " threads: "
				 << f_parallel_experiments[j]["mean"] << " seconds (mean), "
				 << f_parallel_experiments[j]["deviation"] << " deviation, "
				 << f_parallel_experiments[j]["variance"] << " variance, "
				 << f_parallel_experiments[j]["min"] << " min, "
				 << f_parallel_experiments[j]["max"] << " max, "
				 << sequential_results["mean"] / f_parallel_experiments[j]["mean"] << "x" << endl;
		}
	}

	if(compute & 0x01){
		cout << "Parallel by event: " << endl;
		for (int i=1, j=0; i<=max_threads; i=i*2, j++){
			cout << i << " threads: "
				 << parallel_experiments[j]["mean"] << " seconds (mean), "
				 << parallel_experiments[j]["deviation"] << " deviation, "
				 << parallel_experiments[j]["variance"] << " variance, "
				 << parallel_experiments[j]["min"] << " min, "
				 << parallel_experiments[j]["max"] << " max, "
				 << sequential_results["mean"] / parallel_experiments[j]["mean"] << "x" << endl;
		}
	}

	/*
	cout << "Comparing results (with last parallel experiment)..." << endl;
	bool equal = 1;
	for (int i=0; i<num_events_to_process; i++)
		equal = equal && areTracksEqual(sequential_tracks_vector[i], parallel_tracks_vector[i]);
	if (equal)
		cout << endl << "Results are identical! Good job!" << endl;
	else
		cout << endl << "Results differ. Hmm." << endl;

	cout << "Writing results..." << endl;
	for (int i=0; i<num_events; i++)
		printResultTracks(sequential_tracks_vector[i], i, "tracks");
	
	if(!equal)
		for (int i=0; i<num_events; i++)
			printResultTracks(parallel_tracks_vector[i], i, "parallel_tracks");
	*/

	// compareTrackVectors(sequential_tracks_vector, sequential_tracks_vector);

	
	cout << "Done!" << endl;

	// getchar();

}

template<class T>
void clean_vector(vector<vector<T> >& _tracks_vector){
	for(typename vector<vector<T> >::iterator it = _tracks_vector.begin(); it != _tracks_vector.end(); it++){
		(*it).clear();
	}
}

map<string, float> calcResults(vector<float> times){
	// sqrt ( E( (X - m)2) )
	map<string, float> results;
	float deviation = 0.0, variance = 0.0, mean = 0.0, min = infinite(), max = 0.0;

	int n = 0;
	float seconds;
	for(vector<float>::iterator it = times.begin(); it != times.end(); it++){
		n++;
		seconds = (*it);
		mean = (mean * (n - 1) + seconds) / n;
		variance += seconds * seconds;

		if (seconds < min) min = seconds;
		if (seconds > max) max = seconds;
	}

	variance = (variance - times.size() * mean * mean) / times.size();
	deviation = sqrt(variance);

	results["variance"] = variance;
	results["deviation"] = deviation;
	results["mean"] = mean;
	results["min"] = min;
	results["max"] = max;

	return results;
}

void printStatistics(map<string, float> statistics, string pre){
	// Assuming format created before
	cout << pre << "Mean: " << statistics["mean"] << endl
		 << pre << "Deviation, variance: " << statistics["deviation"] << ", " << statistics["variance"] << endl
		 << pre << "Min, max: " << statistics["min"] << ", " << statistics["max"] << endl;
}

float infinite() { 
	static const int value = 0x7f800000; 
	return *(float*)& value;
}

void readFile(string filename, char*& input, bool do_malloc ){
	int size;
	
	// Give me them datas!!11!
	ifstream infile (filename.c_str(), ifstream::binary);

	// get size of file
	infile.seekg(0, ifstream::end);
	size = infile.tellg();
	infile.seekg(0);

	// read content of infile with pointers
	if (do_malloc)
		input = (char*) malloc(size);
	infile.read (input, size);
	infile.close();
}

template<class T, class Z>
void compareTrackVectors(vector<vector<T> > correct_tracks_vector, vector<vector<Z> > new_tracks_vector){
	float mean1 = 0.0, mean2 = 0.0;
	for(int i=0; i<num_events_to_process; i++){
		cout << "event " << i << endl;
		int c = compareTracks<T, Z>(correct_tracks_vector[i], new_tracks_vector[i]);
		cout << endl;

		mean1 = (mean1 * i + (new_tracks_vector[i].size() != 0 ? (float) c / (float) new_tracks_vector[i].size() : 1.0)) / (float) (i+1);
		mean2 = (mean2 * i + (correct_tracks_vector[i].size() != 0 ? (float) c / (float) correct_tracks_vector[i].size() : 1.0)) / (float) (i+1);
	}

	cout << "Mean 1 (correct / no_new_tracks): " << mean1 * 100.0 << "%" << endl
		 << "Mean 2 (correct / no_total_tracks): " << mean2 * 100.0 << "%" << endl;
}

template<class T, class Z>
int compareTracks(vector<T> correct_tracks, vector<Z> new_tracks){
	cout << "Track sizes: " << correct_tracks.size() << ", " << new_tracks.size();
	if( correct_tracks.size() != new_tracks.size() )
		cout << " (differ)";
	cout << endl;

	int num_tracks_equal = 0;
	bool equal = 1;
	for(typename vector<T>::iterator it = correct_tracks.begin(); it != correct_tracks.end(); it++){
		vector<int> hits_correct = (*it).hits;

		for(typename vector<Z>::iterator it2 = new_tracks.begin(); it2 != new_tracks.end(); it2++){
			vector<int> hits_new = (*it2).hits;

			if (hits_correct.size() == hits_new.size()){
				equal = 1;
				for(vector<int>::iterator it3 = hits_correct.begin(); it3 != hits_correct.end(); it3++){
					if (find(hits_new.begin(), hits_new.end(), (*it3)) == hits_new.end()){
						equal = 0;
						break;
					}
				}

				if (equal) num_tracks_equal++;
			}
			if (equal){
				equal = 0;
				break;
			}
		}
	}

	cout << "Equal tracks: " << num_tracks_equal << " ("
		 << (new_tracks.size() != 0 ? (((float) num_tracks_equal) / ((float) new_tracks.size())) * 100.0 : 100) << "%, "
		 << (correct_tracks.size() != 0 ? (((float) num_tracks_equal) / ((float) correct_tracks.size())) * 100.0 : 100)
		 << "%)" << endl;
	return num_tracks_equal;
}

bool areTracksEqual(vector<track> tracks_1, vector<track> tracks_2){
	if(tracks_1.size() != tracks_2.size()){
		cout << "x (tracks size)";
		return false;
	}

	for(int i=0; i<tracks_1.size(); ++i){
		track track_1 = tracks_1[i];
		track track_2 = tracks_2[i];
		if(track_1.hits.size() != track_2.hits.size()){
			cout << "x (hits size)";
			return false;
		}

		for(int j=0; j<track_1.hits.size(); ++j){
			if(track_1.hits[j] != track_2.hits[j]){
				cout << "x (hits)";
				return false;
			}
		}
	}
	cout << ".";
	return true;
}

template<class T>
void printResultTracks(vector<T> tracks, int event_no, string track_folder_container){
	string track_filename = track_folder_container + "//tracks_" + toString<int>(event_no) + ".txt";
	ofstream track_file(track_filename.c_str());
	track_file << std::setprecision(3);

	int z;
	
	int t_no = 0;
	for (typename vector<T>::iterator it = tracks.begin(); it != tracks.end(); it++){
		// info() << format( "Dist%8.3f chi%7.3f ", track.distance( *itH ), track.chi2( *itH ) );
		track_file << "track " << t_no++ /*<< " chi2 " << chi2(&(*it))*/ << std::endl;
		for (vector<int>::iterator ith = (*it).hits.begin(); ith != (*it).hits.end(); ith++){
			z = (*ith);
			track_file << "hit " << hit_ids[(*ith)] << " s " << hit_sensorNums[(*ith)] << " ("
					   << hit_Xs[(*ith)] << ", " << hit_Ys[(*ith)] << ", " << hit_Zs[(*ith)] << ")" << endl;
		}
		track_file << endl;
	}

	track_file.close();
}
