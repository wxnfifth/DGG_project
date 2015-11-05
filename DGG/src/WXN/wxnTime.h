#ifndef _WXNTIME_H_
#define _WXNTIME_H_

#include <Windows.h>
#include <ctime>
#include <string>
using std::string;

class ElapasedTime{
	double PCFreq;
	__int64 CounterStart;
	double previous_time;

public:
	ElapasedTime(bool startNow = true){
		if( startNow == false){
			PCFreq = 0.0;
			CounterStart = 0.0;
			return;
		}
		LARGE_INTEGER li;
		if( !QueryPerformanceFrequency(&li) ){
			fprintf(stderr,  "QueryPerformanceFrequency failed!\n");
			return;
		}
		PCFreq = double(li.QuadPart) / 1000.0;
		//printf( "inipcfreq %lf\n" , (double)PCFreq );
		QueryPerformanceCounter(&li);
		CounterStart = li.QuadPart;
		previous_time = 0;
	}
	void start(){
		LARGE_INTEGER li;
		if( !QueryPerformanceFrequency(&li) ){
			fprintf(stderr,  "QueryPerformanceFrequency failed!\n");
			return;
		}
		PCFreq = double(li.QuadPart) / 1000.0;
		QueryPerformanceCounter(&li);
		CounterStart = li.QuadPart;
		previous_time = 0;
	}
	double getTime(){
	    LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		return double(li.QuadPart-CounterStart)/PCFreq/1000.0;
	}

	void printEstimateTime(double time_interval, double percent) {//percent in 0-1
		double current_time = getTime();
		double past_time = current_time - previous_time;
		if (past_time > time_interval) {
			previous_time = current_time;
			double remain_time = current_time / percent * (1 - percent);
			if (remain_time > 3600) {
				int hour = remain_time / 3600.;
				double seconds = remain_time - hour * 3600.;
				printf("Computed %.0lf percent, time %lf, estimate_remain_time %d h %lf s\n",
					percent * 100.0, current_time, hour, seconds );
			} else{
				printf("Computed %.0lf percent, time %lf, estimate_remain_time %lf\n",
					percent * 100.0, current_time, remain_time);
			}
		}
	}


	void printTime(string tmpMessage = ""){
	    LARGE_INTEGER li;
		QueryPerformanceCounter(&li);
		//printf( "li %lf\n" , (double)li.QuadPart);
		//printf( "pcfreq %lf\n" , PCFreq );
		if( tmpMessage.size() != 0 ) 
			tmpMessage = tmpMessage +",";
		fprintf(stderr, "%s time : %lf seconds\n\r" ,  tmpMessage.c_str() , double(li.QuadPart-CounterStart)/PCFreq/1000.0 );
		return;
	}
};


#endif

