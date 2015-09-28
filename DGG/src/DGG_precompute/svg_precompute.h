#ifndef _SVG_PRECOMPUTE_H_
#define _SVG_PRECOMPUTE_H_
#include "stdafx.h"




void svg_precompute(const string& input_obj_name, const int fixed_k, string& svg_file_name); 

void svg_precompute_fix_neighbor(const string& input_obj_name, const int fixed_k, string& svg_file_name);

void svg_precompute_mmp(const string& input_obj_name, const int fixed_k, string& svg_file_name);

void svg_precompute_hy(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);


void svg_precompute_jiajun(const string& input_obj_name, double eps_vg, string& svg_file_name);
    
void svg_precompute_jiajun_output(const string& input_obj_name, double eps_vg, string& svg_file_name);

void svg_precompute_hy_pruning(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);

void  svg_precompute_hy_multithread(const string& input_obj_name, double eps_vg, string& output_filename, double const_for_theta, int thread_num);

void svg_precompute_hy_fast(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);

void svg_precompute_ich(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);


#endif
