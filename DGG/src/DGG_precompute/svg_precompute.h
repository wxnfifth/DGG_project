#ifndef _SVG_PRECOMPUTE_H_
#define _SVG_PRECOMPUTE_H_
#include "stdafx.h"



string get_DGG_filename(const string& input_obj_name, const string& method_name, double eps_vg, double const_for_theta);

void svg_precompute(const string& input_obj_name, const int fixed_k, string& svg_file_name); 

void svg_precompute_debug(const string& input_obj_name, const int fixed_k);


void svg_precompute_fix_neighbor(const string& input_obj_name, const int fixed_k, string& svg_file_name);

void svg_precompute_mmp(const string& input_obj_name, const int fixed_k, string& svg_file_name);

void svg_precompute_hy(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);

void svg_precompute_jiajun(const string& input_obj_name, double eps_vg, string& svg_file_name);
    
void svg_precompute_jiajun_output(const string& input_obj_name, double eps_vg, string& svg_file_name);

void svg_precompute_hy_pruning(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);

void  svg_precompute_hy_multithread(const string& input_obj_name, double eps_vg, string& output_filename, double const_for_theta, int thread_num);

void svg_precompute_hy_fast(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta);

void svg_precompute_ich(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, bool is_debug_mode);

void svg_precompute_ich_multithread_before_pruning(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta,
	 int thread_num, double& ich_multi_time);

void svg_precompute_ich_multithread(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, int thread_num);


void svg_precompute_LiuYongjin_fixing(const string& input_file_name, double eps_vg, double const_for_theta, const string& svg_file_name);

void svg_precompute_ich_debug(const string& input_obj_filename, const string& debug_svg_filename, const string& neighbor_filename);

#endif
