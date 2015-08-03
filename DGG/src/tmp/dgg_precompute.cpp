
#include "stdafx.h"
#include "svg_precompute.h"


int main(int argc, char** argv)
{
  string input_file_name;
  int fixed_K = 50;// a parameter for our algorithm, default is 50
  double eps_vg;
  string method = "n";
  double const_for_theta = 5;

  if (argc >= 4) {
    method = argv[3];
    if (method == "n" || method == "f") {
      input_file_name = argv[1];
      fixed_K = atoi(argv[2]);
      method = argv[3];
    } else if (method == "j" || method == "y") {
      input_file_name = argv[1];
      eps_vg = atof(argv[2]);
      method = argv[3];
    } else if (method == "h" || method == "p" || method == "d") {
      input_file_name = argv[1];
      eps_vg = atof(argv[2]);
      method = argv[3];
      if (argc == 5) {
        const_for_theta = atoi(argv[4]);
      } else {
        printf("wrong argument, usage example 'SVG_precompute.exe bunny.obj 50  h(hy's method) 10(const num)");
        exit(1);
      }
    } else if (method == "n") {
      input_file_name = argv[1];
      fixed_K = atoi(argv[2]);
      method = argv[3];
    } 
  } else {
    printf("wrong argument, usage example 'SVG_precompute.exe bunny.obj 50 [n(normal) f(fix_neighbor) j(jiajun's method)' y(jiajun's method and jiajun's output) h(hy's method) p(hy's method with prunning) d(dgg's method with prunning and fast)");
    exit(1);
  }

  string svg_file_name;
  if (method == "n") {
    svg_precompute(input_file_name, fixed_K, svg_file_name);
  } else if (method == "f") {
    svg_precompute_fix_neighbor(input_file_name, fixed_K, svg_file_name);
  } else if (method == "m") {
    svg_precompute_mmp(input_file_name, fixed_K, svg_file_name);
  } else if (method == "j") {
    svg_precompute_jiajun(input_file_name, eps_vg, svg_file_name);
  } else if (method == "y") {
    svg_precompute_jiajun_output(input_file_name, eps_vg, svg_file_name);
  } else if (method == "h") {
    svg_precompute_hy(input_file_name, eps_vg, svg_file_name, const_for_theta);
  } else if (method == "p") {
    svg_precompute_hy_pruning(input_file_name, eps_vg, svg_file_name, const_for_theta);
  } else if (method == "d") {
    svg_precompute_hy_fast(input_file_name, eps_vg, svg_file_name, const_for_theta);


  }

         


  return 0;
}

