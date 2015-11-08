import os
import sys
import math
import numpy

obj_name = sys.argv[1]
model_name = obj_name[:-4]
noise_percent_list = numpy.arange(0,0.4,0.05)

k = 500

print noise_percent_list

error_file = '%s_svg_noise_data.txt' % (model_name)

for noise_percent in noise_percent_list:
    cmd_line = '..\\..\\bin\\figure_sphere.exe %s %lf' % (obj_name, noise_percent)
    print cmd_line
    os.system(cmd_line)
    output_obj_name = '%s_noise%.2lf.obj' % (model_name, noise_percent)
    output_model_name = output_obj_name[:-4]
    precompute_log_filename = 'dgg_precompute_%s_k%d.log' % (output_model_name,k)
    svg_precompute_cmd_line = '..\\..\\bin\\dgg_precompute.exe %s %d n 2> %s' % (output_obj_name,k,precompute_log_filename)
    os.system(svg_precompute_cmd_line)
    svg_binary_filename = '%s_SVG_k%d.binary' % (output_model_name,k)
    svg_log_filename=svg_binary_filename[:-7] + '_dij.log'
    svg_dij_cmd_line = '..\\..\\bin\\dgg_lc.exe %s %s dij 2> %s' %(output_obj_name,svg_binary_filename,svg_log_filename)
    print svg_dij_cmd_line
    os.system(svg_dij_cmd_line)
    with open(svg_log_filename,'r') as f:
        for line in f:
            lst = line.split()
            if lst[0] == 'total_average_error':
                average_error = float(lst[1])
            elif lst[0] == 'average_neigh':
                average_neigh = float(lst[1])
            elif lst[0] == 'total_average_running_time':
                average_running_time = float(lst[1])                
    with open(error_file,'a') as f:
        f.write( '%lf\t%lf\t%lf\t%lf\n' % (noise_percent,average_neigh,average_error,average_running_time))
    
    

exit(1)

k_start = int(sys.argv[2])
k_end = int(sys.argv[3])
sample_num = int(sys.argv[4])
#eps_start > eps_end
k_list = []
k_list.append(k_start)
log_interval = (math.log10(k_end) - math.log10(k_start)) / (sample_num -1)
for i in range(1,sample_num-1):
    k_list.append(int(math.pow(10,math.log10(k_start) + log_interval * i)))
    
k_list.append(k_end)

print k_list
    
#obj_path = model_name + '.obj'
model_name = obj_name[:-4]
obj_path = obj_name

error_file = '%s_k_error_file_time.txt' % (model_name)

for k in k_list:
    precompute_log_filename = 'dgg_precompute_%s_k%d.log' % (obj_path,k)
    svg_precompute_cmd_line = 'c:/util/dgg_precompute.exe %s %d n > %s' % (obj_path,k,precompute_log_filename)
    svg_binary_filename = '%s_SVG_k%d.binary' % (model_name,k)
    if not os.path.isfile(svg_binary_filename):
        print svg_precompute_cmd_line
        os.system(svg_precompute_cmd_line)
    svg_log_filename=svg_binary_filename[:-7] + '_dij.log'
    svg_dij_cmd_line = 'c:/util/dgg_lc.exe %s %s dij > %s' %(obj_path,svg_binary_filename,svg_log_filename)
    print svg_dij_cmd_line
    os.system(svg_dij_cmd_line)
    
    with open(svg_log_filename,'r') as f:
        for line in f:
            lst = line.split()
            if lst[0] == 'total_average_error':
                average_error = float(lst[1])
            elif lst[0] == 'average_neigh':
                average_neigh = float(lst[1])
            elif lst[0] == 'total_average_running_time':
                average_running_time = float(lst[1])                
    with open(error_file,'a') as f:
        f.write( '%lf\t%lf\t%lf\n' % (average_neigh,average_error,average_running_time))
    

#os.system('dir >tmp.txt')
#rt = subprocess.check_output(*popenargs, **kwargs)  






