import os
import sys
import math

#print sys.argv

obj_name = sys.argv[1]

eps_start = float(sys.argv[2])
eps_end = float(sys.argv[3])
sample_num = int(sys.argv[4])
method = sys.argv[5]
input_const = int(sys.argv[6])
const_method = sys.argv[7]
if const_method != 'angle' and const_method != 'choose':
    print 'sys.argv[7] must be angle or choose'
    quit()

#eps_start > eps_end
eps_list = []
eps_list.append(eps_end)
if sample_num > 1:
    log_interval = (math.log10(eps_start) - math.log10(eps_end)) / (sample_num -1)
    for i in range(1,sample_num-1):
        eps_list.append(math.pow(10,log_interval * i + math.log10(eps_end)))
    eps_list.append(eps_start)
    eps_list.reverse()
    
print eps_list
##for e in eps_list:
##    print e
##quit()

#obj_path = model_name + '.obj'
model_name = obj_name[:-4]
obj_path = obj_name

if method == 'fan':
    error_file = '%s_error_file_fan_cc%.0f.txt' % (model_name,input_const)
elif method == 'lc':
    error_file = '%s_error_file_lc_cc%.0f.txt' % (model_name,input_const)
else:
    print 'method is fan or lc'
    quit()

    
for eps in eps_list:
    if const_method == 'angle':
        constant = int(input_const / 180.0 * math.pi / math.asin(math.sqrt(eps)) + 1)
    elif const_method == 'choose':
        constant = input_const
    #constant = 5
    precompute_log_filename = 'dgg_precompute_%s_%f_%d.log' % (obj_path,eps,constant)
    svg_precompute_cmd_line = r'..\bin\dgg_precompute.exe %s %f p %d 2> %s 1>&2' % (obj_path,eps,constant,precompute_log_filename)
    svg_binary_filename = '%s_DGG%f_c%d_pruning.binary' % (model_name,eps,constant)
    if not os.path.isfile(svg_binary_filename):
        print svg_precompute_cmd_line
        os.system(svg_precompute_cmd_line)
        flag_preprocess = True
        #svg_tmp_binary_filename = '%s_DGG%f_c%d.binary' % (model_name,eps,constant)
        #os.system('del ' + svg_tmp_binary_filename)
    else:
        flag_preprocess = False
    svg_log_filename=svg_binary_filename[:-7] + '_hy.log'
    #if not os.path.isfile(svg_log_filename):
    if method == 'fan':
        svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s hy 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
    elif method == 'lc':
        svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s lll 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
        
    print svg_lll_cmd_line
    os.system(svg_lll_cmd_line)
    print svg_log_filename
    with open(svg_log_filename,'r') as f:
        for line in f:
            lst = line.split()
            if len(lst) == 0:
                continue
            if lst[0] == 'total_average_error':
                average_error = float(lst[1])
            elif lst[0] == 'average_neigh':
                average_neigh = float(lst[1])
            elif lst[0] == 'total_average_running_time':
                average_running_time = float(lst[1])
    if flag_preprocess:
        with open(precompute_log_filename,'r') as f:
            for line in f:
                lst = line.split()
                #print lst
                if lst[0] == 'Average_degree_before':
                    average_neigh_before = float(lst[1])       
        with open(error_file,'a') as f:
            f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps,average_neigh_before))
    else:
        with open(error_file,'a') as f:
            f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps))
        

#os.system('dir >tmp.txt')
#rt = subprocess.check_output(*popenargs, **kwargs)  






