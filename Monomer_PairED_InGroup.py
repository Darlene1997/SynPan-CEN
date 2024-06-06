import pandas as pd
import argparse
import Levenshtein
import os


def calculate_ED_in_groups(ref, que, each_chr, Dir):
    os.system('echo "############## 1 calculate_ED_in_groups #####################" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    if os.path.exists(str(Dir) + "01_AllMatchList/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_all.matchList"):
        os.system('echo "Do %s_%s_%s_group_pairED_all.matchList exists!" >> %s_%s_%s_test.log' % (each_chr, ref, que, each_chr, ref, que))
    else:
        df_ref_monomer = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt" , sep="\t", header=None, names=["Chr", "Start", "End", "Strand", "New_id", "Group"])
        df_ref_seq = pd.read_table("../../0_monomer_anno/CR_CENtype_repeats_" + str(ref) + ".csv", sep=",") ##start,end,length,seq,strand,class,filename,Chr,edit.distance,repetitiveness
        df_que_monomer = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["Chr", "Start", "End", "Strand", "New_id", "Group"])
        df_que_seq = pd.read_table("../../0_monomer_anno/CR_CENtype_repeats_" + str(que) + ".csv", sep=",")

        outfile = str(Dir) + "01_AllMatchList/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_all.matchList"
        o = open(outfile, 'w')
        group_list = ["group1", "group2", "group3", "group4", "group5", "group6", "group7", "group8", "group9", "group10", "group11", "group12", "group13", "group14", "group15",
                      "1+", "2+", "3+", "4+", "5+", "6+", "7+", "8+", "9+", "10+", "11+", "12+", "13+", "14+", "15+"]
        for group in group_list:
            #### Pairwise comparisons are made ONLY within the same group
            df_g_ref = df_ref_monomer[(df_ref_monomer["Group"] == str(group)) & (df_ref_monomer["Chr"] == str(each_chr) + "_" + str(ref))]
            df_g_que = df_que_monomer[(df_que_monomer["Group"] == str(group)) & (df_ref_monomer["Chr"] == str(each_chr) + "_" + str(que))]
            for index_r, row_r in df_g_ref.iterrows():
                start_r = row_r["Start"]
                end_r = row_r["End"]
                chr_r = row_r["Chr"]
                New_id_r = row_r["New_id"]
                if New_id_r.split(".")[-2] != "LTR":
                    df_r_seq = df_ref_seq[(df_ref_seq["start"] == float(start_r)) & (df_ref_seq["end"] == float(end_r))].reset_index(drop=True)
                    r_seq = df_r_seq.loc[0, "seq"]
                    for index_q, row_q in df_g_que.iterrows():
                        start_q = row_q["Start"]
                        end_q = row_q["End"]
                        chr_q = row_q["Chr"]
                        New_id_q = row_q["New_id"]
                        if New_id_q.split(".")[-2] != "LTR":
                            df_q_seq = df_que_seq[(df_que_seq["start"] == float(start_q)) & (df_que_seq["end"] == float(end_q))].reset_index(drop=True)
                            q_seq = df_q_seq.loc[0, "seq"]
                            Edit = Levenshtein.distance(str(r_seq), str(q_seq))
                            prop = 2 * Edit / (len(r_seq) + len(q_seq))
                            o.write(str(chr_r) + "\t" + str(New_id_r) + "\t" + str(start_r) + "\t" + str(end_r) + "\t" +str(chr_q) + "\t" + str(New_id_q) + "\t" + str(start_q) + "\t" + str(end_q) + "\t" + str(prop/10000) + "\t" + str(Edit) + "\n")
        os.system('echo "Do %s_%s_%s_group_pairED_all.matchList finished!" >> %s_%s_%s_test.log' % (each_chr, ref, que, each_chr, ref, que))

def change_format_of_dag(ref, que, each_chr, Dir):
    os.system('echo "############## 2 change_format_of_dag #####################" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    if os.path.exists(str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag") and \
            os.path.getsize(str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"):
        os.system('echo "Do %s_%s_%s_group_pairED_order_all.dag exists!" >> %s_%s_%s_test.log' % (each_chr, ref, que, each_chr, ref, que))
    else:
        LTR_infile = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_LRT_blast_best_hit.txt"
        infile = str(Dir) + "01_AllMatchList/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_all.matchList"

        out_all_dag = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_all.dag"
        ########### add LTR, if exist
        if os.path.exists(LTR_infile):
            os.system("cat %s %s > %s" % (infile, LTR_infile, out_all_dag))
        else:
            os.system("cat %s > %s" % (infile, out_all_dag))

        ########### change the format, and use their orders as positions
        out_re_format = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"
        with open(out_all_dag, 'r') as inf:
            with open(out_re_format, 'w') as outf:
                for line in inf: ## "chr_r", "id_r", "s_r", "e_r", "chr_q", "id_q", "s_q", "e_q", "prop", "ED"
                    text = line.strip().split("\t")
                    if text[1].split(".")[-2] == "LTR":
                        chr_r = text[0]
                        id_r = text[1]
                        index_r = id_r.split(".")[-1]
                        chr_q = text[4]
                        id_q = text[5]
                        index_q = id_q.split(".")[-1]
                        ED = text[9]
                        prop = text[8]
                        outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\t" + str(ED) + "\n")
                    else:
                        chr_r = text[0]
                        id_r = text[1]
                        index_r = id_r.split(".")[-1]
                        chr_q = text[4]
                        id_q = text[5]
                        index_q = id_q.split(".")[-1]
                        ED = text[9]
                        if int(ED) == 0:
                            prop = 1e-50
                        elif int(ED) <= 3 and int(ED) > 0:
                            prop = 1e-40
                        elif int(ED) <= 8 and int(ED) > 3:
                            prop = 1e-30
                        elif int(ED) <= 15 and int(ED) > 8:
                            prop = 1e-20
                        elif int(ED) <= 20 and int(ED) > 15:
                            prop = 1e-10
                        else:
                            prop = text[8]
                        outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\t" + str(ED) + "\n")
        os.system("rm %s" % (out_all_dag))
        os.system('echo "Do %s_%s_%s_group_pairED_order_all.dag finished!" >> %s_%s_%s_test.log' % (each_chr, ref, que, each_chr, ref, que))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-chr', '--chr', help='Chr', dest='chr', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    calculate_ED_in_groups(args.ref, args.que, args.chr, args.Dir)
    change_format_of_dag(args.ref, args.que, args.chr, args.Dir)


if __name__ == '__main__':
    main()
