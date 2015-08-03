import os
import sys
import math

def main():
    obj_name = sys.argv[1]

    num_eps = int(sys.argv[2])

    eps_list = []
    for i in range(3,3+num_eps):
        eps_list.append(float(sys.argv[i]))

    model_name = os.path.basename(obj_name)[:-4]
    dir_name = os.path.dirname(obj_name)
    obj_path = obj_name

    if len(sys.argv) > 3 + num_eps:
        constant = float(sys.argv[3+num_eps])
    else:
        constant = 20

    error_file = os.path.join(dir_name, '%s_eps_with_k_cc%.0f.txt' % (model_name,constant))

    for eps in eps_list:
        precompute_log_filename = os.path.join(dir_name, 'dgg_precompute_%s_%f_%d.log' % (model_name,eps,constant))
        svg_precompute_cmd_line = r'..\bin\dgg_precompute.exe %s %f p %d 2> %s 1>&2' % (obj_path,eps,constant,precompute_log_filename)
        svg_binary_filename = os.path.join(dir_name,'%s_DGG%f_c%d_pruning.binary' % (model_name,eps,constant))
        if not os.path.isfile(svg_binary_filename):
            print svg_precompute_cmd_line
            os.system(svg_precompute_cmd_line)
            flag_preprocess = True
        else:
            flag_preprocess = False
        svg_log_filename=svg_binary_filename[:-7] + '_hy.log'
        method = 'fan'
        if method == 'fan':
            svg_lll_cmd_line = r'..\bin\dgg_lc.exe %s %s hy 2> %s 1>&2' %(obj_path,svg_binary_filename,svg_log_filename)
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
    
        with open(error_file,'a') as f:
            f.write( '%lf\t%.10lf\t%lf\t%.10lf\t%lf\n' % (average_neigh,average_error,average_running_time,eps,average_error/eps))
        eps_error = average_error
        x = math.log(eps_error)
        estimate_k = int(-2.676 * x * x * x - 45.47 * x * x - 279.9 * x - 584.6)
        low_k = estimate_k - 50
        high_k = estimate_k + 50
        while True:
            low_k_error = find_k_error(dir_name,model_name,obj_path,low_k,error_file)
            if low_k_error > eps_error:
                low_k -= 100
            else:
                break

        while True:
            high_k_error = find_k_error(dir_name,model_name,obj_path,high_k,error_file)
            if high_k_error < eps_error:
                high_k += 100
            else:
                break
        print low_k, high_k
        while True:
            mid_k = int((low_k+high_k) / 2)
            mid_k_error = find_k_error(dir_name,model_name,obj_path,mid_k,error_file)
            if mid_k_error < eps_error:
                low_k = mid_k
            else:
                high_k = mid_k
            if high_k - low_k <= 2:
                break
            



    

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




