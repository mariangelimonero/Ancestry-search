import argparse
import random

parser = argparse.ArgumentParser(description="Ancestry analysis using Admixture, Structure and RFMix.\n")
parser.add_argument("--rs_list", help="Name for rsnumbers file.", required=True)
parser.add_argument("--ancestry_file", help="Name for ancestry file.", required=True)
parser.add_argument("--random", help="Take N random rsnumbers from whole genome.\n")
parser.add_argument("--o", help="Name for output file.", required=True)
args = parser.parse_args()

rsnumber_list = set()
intersection_rs_list = []
random_rs_list = []

with open(args.rs_list,'r') as file:
    for rs in file:
        rsnumber_list.add(rs.strip('\n'))

with open(args.ancestry_file,'r') as anc, open(args.o,'w') as out:
    print("Searching ancestries for queried rsnumbers...")
    for row in anc:
        if "rsnumber" in row:
            out.write(row)
            print("checked if")
        else:
            rsnumber = [row.strip('\n').split(',')[0]]
            random_rs_list.append(row)
            print("rsnumber in random_rs_list")
            if set(rsnumber).intersection(rsnumber_list):
                intersection_rs_list.append(row)
                print("rsnumber intersection in intersection_rs_list")
            else:
                random_rs_list.append(row)
                print("Random rsnumbers in random_rs_list")
    if args.random:
        print("Randomly searching " + args.random + " rsnumbers ancestries from whole genome...")
        random.shuffle(random_rs_list)
        for i in range(int(args.random)):
            out.write(random_rs_list.pop())
