import os
import math
import argparse
#print sys.argv

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_name", type=str, help="input file name",default='')
parser.add_argument('input_eps',type=float,nargs='+',help = 'list of eps')
parser.add_argument('-m','--method',type=str,choices=['fan','lc','dggdij'], help='chois is fan or lc or dggdij')
parser.add_argument('-o','--input_const',type=int,help = 'input chost')
parser.add_argument('-c','--const_method',type=str,choices=['angle','choose'], help='const_method is choose or angle')
#parser.add_argument("-f",'--face_num',type=int,help='number of face ',default=0)
args = parser.parse_args()


#obj_name = sys.argv[1]
obj_name = args.input_name

#method = sys.argv[5]
#input_const = int(sys.argv[6])
#const_method = sys.argv[7]
eps_list = args.input_eps
print eps_list
method = args.method
input_const = args.input_const
const_method = args.const_method
#eps_start > eps_end
#eps_list = []
#eps_list.append(eps_end)
#if sample_num > 1:
#    log_interval = (math.log10(eps_start) - math.log10(eps_end)) / (sample_num -1)
#    for i in range(1,sample_num-1):
#        eps_list.append(math.pow(10,log_interval * i + math.log10(eps_end)))
#    eps_list.append(eps_start)
#    eps_list.reverse()
##for e in eps_list:
##    print e
##quit()

#obj_path = model_name + '.obj'
model_name = os.path.basename(obj_name)[:-4]
dir_name = os.path.dirname(obj_name)
obj_path = obj_name

if method == 'fan':
    error_file = os.path.join(dir_name,'%s_error_file_ich_fan_cc%.0f.txt' % (model_name,input_const))
elif method == 'lc':
    error_file = os.path.join(dir_name,'%s_error_file_ich_lc_cc%.0f.txt' % (model_name,input_const))
elif method == 'dggdij':
    error_file = os.path.join(dir_name,'%s_error_file_ich_dggdij_cc%.0f.txt' % (model_name,input_const))
else:
    print 'method is fan or lc'
    quit()

print error_file

    
for eps in eps_list:
    if const_method == 'angle':
        constant = int(input_const / 180.0 * math.pi / math.asin(math.sqrt(eps)) + 1)
    elif const_method == 'choose':
        constant = input_const
    #constant = 5
    precompute_log_filename = os.path.join(dir_name,'dggich_precompute_%s_%.10f_%d.log' % (model_name,eps,constant))
    svg_precompute_cmd_line = r'..\bin\dgg_precompute.exe %s %.10f im %d 4 2> %s 1>&2' % (obj_path,eps,constant,precompute_log_filename)
    svg_binary_filename = os.path.join(dir_name,'%s_DGGICH%.10f_c%d_pruning.binary' % (model_name,eps,constant))
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
        svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s hy dbl 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
    elif method == 'lc':
        svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s lll dbl 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
    elif method == 'dggdij':
        svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s dggdij dbl 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
	        
    print svg_lll_cmd_line
    os.system(svg_lll_cmd_line)
    print svg_log_filename
    average_error = 0
    average_neigh = 0
    average_running_time = 0
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
        average_neigh_before = 0
        preprocessing_time = 0
        with open(precompute_log_filename,'r') as f:
            for line in f:
                lst = line.split()
                #print lst
                if len(lst) == 0:
                    continue
                elif lst[0] == 'Average_degree_before':
                    average_neigh_before = float(lst[1])
                elif lst[0] == 'total_time_and_pruning':
                    preprocessing_time = float(lst[1])
        with open(error_file,'a') as f:
            f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\t%lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps,average_neigh_before,preprocessing_time))
    else:
        with open(error_file,'a') as f:
            f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps))
        

#os.system('dir >tmp.txt')
#rt = subprocess.check_output(*popenargs, **kwargs)  






