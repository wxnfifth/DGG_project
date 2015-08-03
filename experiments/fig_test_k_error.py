import os
import sys
import math

#print sys.argv

obj_path = sys.argv[1]

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
model_name = os.path.basename(obj_path)[:-4]
dir_name = os.path.dirname(obj_path)

error_file = os.path.join(dir_name,'%s_k_error_file_time.txt' % (model_name))

for k in k_list:
    precompute_log_filename = os.path.join(dir_name,'dgg_precompute_%s_k%d.log' % (model_name,k))
    svg_precompute_cmd_line = r'..\bin\dgg_precompute.exe %s %d n 2> %s' % (obj_path,k,precompute_log_filename)
    svg_binary_filename = os.path.join(dir_name,'%s_SVG_k%d.binary' % (model_name,k))
    if not os.path.isfile(svg_binary_filename):
        print svg_precompute_cmd_line
        os.system(svg_precompute_cmd_line)
    svg_log_filename=svg_binary_filename[:-7] + '_dij.log'
    svg_dij_cmd_line = r'..\bin\dgg_lc.exe %s %s dij 2> %s' %(obj_path,svg_binary_filename,svg_log_filename)
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
        f.write( '%.10f\t%.10f\t%.10f\t%.10f\n' % (average_neigh,average_error,average_running_time,k))
    

#os.system('dir >tmp.txt')
#rt = subprocess.check_output(*popenargs, **kwargs)  






