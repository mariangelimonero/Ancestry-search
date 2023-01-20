#Bracket_search: Determinar posiciones cercanas (+-0.5cM) entre la lista de CVD y la lista de GWAS-Pablo.
#Funciona
#fecha: 06-12-2022
#Output: Genera 3 documentos que son el args.o (intersección entre CVD y GWAS), bracket_search.csv (match de marcadores entre 0.5cM) y bracket_search_final.csv(match de marcadores entre 0.5cM y sus respectivas ancestrías).
#Linea de comando: python3 Parte3.py --rs_list 60marcadoresgeneticos_copy.txt --ancestry_file rsnumber_ancestry_all.csv --o Parte3_06_12_22.csv

import argparse
import pandas as pd
import csv

parser = argparse.ArgumentParser(description="Ancestry analysis using Admixture, Structure and RFMix.\n")
parser.add_argument("--rs_list", help="Name for literature and GWAS genetic markers positions file (.txt file).", required=True)
parser.add_argument("--ancestry_file", help="Name for ancestry file (.csv file).", required=True)
parser.add_argument("--o", help="Name for output file.", required=True)
args = parser.parse_args()

rsnumber_list = set()
random_rs_list = []


with open(args.rs_list,'r') as file:
    for rs in file:
        rsnumber_list.add(rs.strip('\n'))
        print("rs_list read")

with open(args.ancestry_file,'r') as anc, open(args.o,'w') as out:
    print("Searching ancestries for queried rsnumbers...")
    for row in anc:
        if "rsnumber" in row:
            out.write(row)
        else:
            rsnumber = [row.strip('\n').split(',')[2]]
            random_rs_list.append(row)
            if set(rsnumber).intersection(rsnumber_list):
                out.write(row)
                out.close
                print("done")

with open ("bracket_search_results.csv", 'w+') as info:
    df = pd.read_csv(args.o)
    print("csv file read")
    df.reset_index(drop=True, inplace=True)
    writer = csv.writer(info,lineterminator = "\n")
    head = ["POS","POS2","CM1","CM2"]
    writer.writerow(head)
    for i in range(len(df)) :
        x = df.loc[i, "POS"]
        y = df.loc[i, "cM"]
        y2 = df.loc[i, "cM"] + 0.5
        y3 = df.loc[i, "cM"] - 0.5
        print("bracket search")
        for i in range(len(df)) :
            z = df.loc[i, "cM"]
            n2 = df.loc[i, "POS"]
            if (z <= y2 and z >= y3) and ( n2 != x) :
                row = [x,df.loc[i, "POS"],y, z]
                writer.writerow(row)
                print(row)

df_pep =  pd.read_csv("bracket_search_results.csv")

df_pep2 = pd.merge(df_pep,df,on = "POS")

df_pep2.to_csv("bracket_search_final.csv",index = False)
