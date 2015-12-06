import os
import math
import argparse
#print sys.argv

def main():
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
        svg_precompute_cmd_line = r'..\bin\dgg_precompute.exe %s %f im %d 4 2> %s 1>&2' % (obj_path,eps,constant,precompute_log_filename)
        svg_binary_filename = os.path.join(dir_name,'%s_DGGICH%f_c%d_pruning.binary' % (model_name,eps,constant))
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
        elif method == 'dggdij':
            svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s dggdij 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
    	        
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
            ich_multi_time = 0
            pruning_time = 0
            with open(precompute_log_filename,'r') as f:
                for line in f:
                    lst = line.split()
                    #print lst
                    if len(lst) < 2:
                        continue
                    elif lst[0] == 'Average_degree_before':
                        average_neigh_before = float(lst[1])
                    elif lst[0] == 'ich_multi_time,':
                        ich_multi_time = float(lst[3])
                    elif lst[0] == 'prunning' and lst[1] == 'time':
                        pruning_time = float(lst[2])
            preprocessing_time = ich_multi_time * 3.3 + pruning_time
            with open(error_file,'a') as f:
                f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\t%lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps,average_neigh_before,preprocessing_time))
        else:
            with open(error_file,'a') as f:
                f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps))
        if average_error == 0:
            continue
        eps_error = average_error
        x = math.log(eps_error)
        estimate_k = int(-2.676 * x * x * x - 45.47 * x * x - 279.9 * x - 584.6)
        print estimate_k
        low_k = estimate_k - 50
        if low_k < 20:
            low_k = 20
        high_k = estimate_k + 50
        while True:
            low_k_error = find_k_error(dir_name,model_name,obj_path,low_k,error_file)
            if low_k_error < eps_error:
                low_k -= 50
                if low_k < 30:
                    low_k = 30
                    break
            else:
                break

        while True:
            high_k_error = find_k_error(dir_name,model_name,obj_path,high_k,error_file)
            if high_k_error > eps_error:
                high_k += 50
            else:
                break
        print low_k, high_k
        while True:
            mid_k = int((low_k+high_k) / 2)
            mid_k_error = find_k_error(dir_name,model_name,obj_path,mid_k,error_file)
            if eps_error > mid_k_error:
                high_k = mid_k
            else:
                low_k = mid_k
            if high_k - low_k <= 2:
                break
            
            
    
    #os.system('dir >tmp.txt')
    #rt = subprocess.check_output(*popenargs, **kwargs)  

def find_k_error(dir_name, model_name, obj_path, k,error_file):
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
            if len(lst) == 0:
                continue
            if lst[0] == 'total_average_error':
                average_error = float(lst[1])
            elif lst[0] == 'average_neigh':
                average_neigh = float(lst[1])
            elif lst[0] == 'total_average_running_time':
                average_running_time = float(lst[1])                
    with open(error_file,'a') as f:
        f.write( '%.10f\t%.10f\t%.10f\t%.10f\n' % (average_neigh,average_error,average_running_time,k))
    return average_error


main()



