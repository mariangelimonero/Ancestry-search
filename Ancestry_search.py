import os
import subprocess
import argparse
import multiprocessing as mp
import pandas as pd
from operator import add, truediv

parser = argparse.ArgumentParser(description="Ancestry analysis using Admixture, Structure and RFMix.\n")
parser.add_argument("--locfolder", help="Name for location file.", required=True)
parser.add_argument("--fbfolder", help="Name for forward-backward file.", required=True)
parser.add_argument("--o", help="Name of the output string used for the pipeline run.", required=True)
args = parser.parse_args()

def run_ancestry_finder(fb,locfile,sampleindex,samplename):
    _fb_map = []
    _item = []
    _loc_map = []

    print("Entering " + samplename + " data matrix...")
    print("Gathering ForwardBackward file data...")
    with open(fb,'r') as fb_file:
        for item in fb_file:
            _item = item.split(' ')
            _fb_map.append(_item)
        fb_file.close()

    print("Gathering Location file data...")
    with open(locfile,'r') as loc_file:
        for item in loc_file:
            _loc_map.append(item.strip('\n'))
        loc_file.close()

    print("Gathering ancestry data...")
    ploidy = 2
    num_pops = 3
    col_num = int(sampleindex) * num_pops * ploidy - num_pops * ploidy
    sample_filename = samplename + "_all_chr.anc"
    with open(sample_filename,'a') as out:
        print('Writing ' + sample_filename + ' file...')
        for i,rs in enumerate(_loc_map):
            _tmp_out = [_fb_map[i][col_num],_fb_map[i][col_num+1],_fb_map[i][col_num+2]]
            out.write(rs + '\t' + str(_tmp_out) + '\n')
        out.close()
        print('Closing ' + sample_filename + ' file.')
    return

def pool_runner(_fb,_loc,ind,name):
    for i in range(0,22):
        run_ancestry_finder(_fb[i],_loc[i],ind,name)
    return

if __name__ == '__main__':
    print("Opening data files...")
    _sample_names = []
    with open('Results/' + args.o + '/Sample21_30_03_29_22.txt','r') as rsn:
        for name in rsn:
            _sample_names.append(name.strip('\n'))

    _fb_file_list = []
    for i in range(1,23):
        for file in os.listdir(args.fbfolder):
            if file.endswith('.1.ForwardBackward.txt') and ('chr'+str(i)+'.') in file:
                _fb_file_list.append(os.path.join(args.fbfolder,file))

    _loc_file_list = []
    for i in range(1,23):
        for file in os.listdir(args.locfolder):
            if ('chr'+str(i)+'.') in file:
                _loc_file_list.append(os.path.join(args.locfolder,file))

    if os.cpu_count() > 2:
        process_pool = mp.Pool(os.cpu_count()-1)
    else:
        process_pool = mp.Pool(os.cpu_count()-1)

    print("Matching ancestries with samples...")

    _args = []
    for n in range(0,len(_sample_names)):
        _args.append((_fb_file_list,_loc_file_list,n+1,_sample_names[n]))

    process_pool.starmap(pool_runner,_args)
    process_pool.close()
    process_pool.join()

    print("Getting rsnumber list...")
    _rs_names = ['rsnumber']
    _chr = ['chr']
    with open("Results/" + args.o + '/' + args.o + '.bim') as o_bim:
        for line in o_bim:
            line = line.strip('\n')
            _rs_names.append(line.split('\t')[1])
            _chr.append(line.split('\t')[0])
        o_bim.close()
    _rs_names = pd.DataFrame(_rs_names)
    _chr = pd.DataFrame(_chr)
    print("Concatenating rsnumber data with chromosome number...")
    _rs_names = pd.concat([_rs_names,_chr],axis=1)

    print("Getting positions list...")
    _all_data = [['POS','cM']]
    try:
        with open("Results/" + args.o + "/All_chr.positions") as fpos:
            for pos in fpos:
                pos = pos.strip('\n')
                _all_data.append(pos.split('\t'))
            fpos.close()
        _all_data = pd.DataFrame(_all_data)
        print("Concatenating position data with rsnumber data...")
        _all_data = pd.concat([_rs_names,_all_data],axis=1)
        _anc_data = pd.DataFrame()
        _si = 1
        _avg_anc = []
    except Exception as error :
        print('Error: '+error)

    print("Setting base average to 0...")
    for i in range(len(_rs_names)-1):
        _avg_anc.append([0,0,0])

    print("Starting data concatenation...")
    for file in os.listdir(os.getcwd()):
        if file.endswith('.anc'):
            with open(file,'r') as f:
                print('Concatenating ' + str(_si) + '/' + str(len(_sample_names)) + ' sample with lists')
                print ("\033[A\033[A")
                _tmp = []
                _tmp.append(str(file.strip('_all_chr.anc')))
                _i = 0
                for row in f:
                    row = row.strip("\n").replace("'", "").replace(" ","")
                    _tmp_row = row.split('\t')[-1]
                    _tmp.append(_tmp_row)
                    _tmp_row = _tmp_row.strip('[]').split(',')
                    _tmp_row = [float(x) for x in _tmp_row]
                    for i in range(0,2):
                        _avg_anc[_i][i] += _tmp_row[i]
                    _i += 1
                _anc_data = pd.concat([_anc_data,pd.DataFrame(_tmp)],axis=1)
                _si+=1
                f.close()

    _all_avgs = ['average_ancestry']
    for item in _avg_anc:
        _all_avgs.append([round(item[x] / len(_sample_names),5) for x in range(len(item))])
    _all_data = pd.concat([_all_data,pd.DataFrame(_all_avgs)],axis=1)
    _all_data = pd.concat([_all_data,_anc_data],axis=1)
    print('\rSaving all data to csv...')
    _all_data.to_csv('Results/'+ args.o + '/rsnumber_ancestry_all_3.csv',index=False,header=False)
    print("Cleaning up...")
    subprocess.call('rm -f *.anc',shell=True)
    print("-----Finished-----")
