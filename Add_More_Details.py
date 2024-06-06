import pandas as pd
import argparse
import os

def add_black_region(ref, que, each_chr, ED, Dir):
    output = str(ref) + "_" + str(que) + "_12chr_for_plot_ED" + str(ED) + ".txt"
    o = open(output, 'a')
    ############ blank region
    os.system('echo " 5 blank region" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    ####################### BEGIN to mach ##################
    df_ref_anno = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt", sep="\t",header=None, names=["chr_ref", "s_ref", "e_ref", "strand_ref", "id_ref","Group"])  ##包括所有monomer和LTR
    df_que_anno = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt", sep="\t",header=None, names=["chr_que", "s_que", "e_que", "strand_que", "id_que", "Group"])
    df_ref_newID = df_ref_anno[df_ref_anno["chr_ref"] == str(each_chr) + "_" + str(ref)].reset_index(drop=True)
    df_que_newID = df_que_anno[df_que_anno["chr_que"] == str(each_chr) + "_" + str(que)].reset_index(drop=True)
    line_total = df_ref_newID.shape[0] - 1
    for index, row in df_ref_newID.iterrows():
        End = row["e_ref"]
        if index < line_total:
            nest_index = int(index) + 1
            Start_next = df_ref_newID.loc[nest_index, "s_ref"]
            if int(End) + 1 == int(Start_next):
                pass
            else:
                blank_s_r = End - (df_ref_newID["s_ref"].min())
                blank_e_r = Start_next - (df_ref_newID["s_ref"].min())
                blank_id_r = str(ref) + "_blank_" + str(index)
                o.write(str(each_chr) + "\t" + str(blank_id_r) + "\t" + str(blank_s_r) + "\t" + str(blank_e_r) + "\t5\t4.7\t" + str(ref) + "\tpair_refblank_" + str(index) + "\tblank\n")


    line_total2 = df_que_newID.shape[0] - 1
    for index2, row2 in df_que_newID.iterrows():
        End2 = row2["e_que"]
        if index2 < line_total2:
            nest_index2 = int(index2) + 1
            Start_next2 = df_que_newID.loc[nest_index2, "s_que"]
            if int(End2) + 1 == int(Start_next2):
                pass
            else:
                blank_s_q = End2 - (df_que_newID["s_que"].min())
                blank_e_q = Start_next2 - (df_que_newID["s_que"].min())
                blank_id_q = str(que) + "_blank_" + str(index2)
                o.write(str(each_chr) + "\t" + str(blank_id_q) + "\t" + str(blank_s_q) + "\t" + str(blank_e_q) + "\t2\t2.3\t" + str(que) + "\tpair_queblank_" + str(index2) + "\tblank\n")


def add_strand_track(ref, que, each_chr, ED, Dir):
    output = str(ref) + "_" + str(que) + "_12chr_for_plot_ED" + str(ED) + ".txt"
    o = open(output, 'a')
    ############ blank region
    os.system('echo " 5 blank region" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    ####################### BEGIN to mach ##################
    df_ref_anno = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["chr_ref", "s_ref", "e_ref", "strand_ref", "id_ref", "Group"])  ##包括所有monomer和LTR
    df_que_anno = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["chr_que", "s_que", "e_que", "strand_que", "id_que", "Group"])
    df_ref_newID = df_ref_anno[df_ref_anno["chr_ref"] == str(each_chr) + "_" + str(ref)].reset_index(drop=True)
    df_que_newID = df_que_anno[df_que_anno["chr_que"] == str(each_chr) + "_" + str(que)].reset_index(drop=True)
    for index_r, row_r in df_ref_newID.iterrows():
        strand_id_r = str(ref) + "_strand_" + str(index_r)
        LTR_id_r = str(ref) + "_LTRdetail_" + str(index_r)
        s_r = row_r["s_ref"]
        e_r = row_r["e_ref"]
        s_new_r = s_r - (df_ref_newID["s_ref"].min())
        e_new_r = e_r - (df_ref_newID["s_ref"].min())
        strand_r = row_r["strand_ref"]
        group = row_r["Group"]
        if group != "LTR" and group != "other_repeat" and group != "Gene" and group != "Gap":
            if strand_r == "+":
                o.write(str(each_chr) + "\t" + str(strand_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t5.9\t5.6\t" + str(que) + "\tstrandref_" + str(index_r) + "\tplus\n")
            elif strand_r == "-":
                o.write(str(each_chr) + "\t" + str(strand_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t5.9\t5.6\t" + str(que) + "\tstrandref_" + str(index_r) + "\tminus\n")
            else:
                pass
                #o.write(str(each_chr) + "\t" + str(strand_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t5.9\t5.6\t" + str(que) + "\tstrandref_" + str(index_r) + "\tequal\n")
        elif group == "LTR":
            ltrid = row_r["id_ref"]
            if "SZ-22" in ltrid:
                o.write(str(each_chr) + "\t" + str(LTR_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t6.8\t6.5\t" + str(que) + "\tLTRdetailref_" + str(index_r) + "\tSZ-22\n")
            elif "RIRE" in ltrid:
                o.write(str(each_chr) + "\t" + str(LTR_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t6.8\t6.5\t" + str(que) + "\tLTRdetailref_" + str(index_r) + "\tRIRE\n")
            else:
                o.write(str(each_chr) + "\t" + str(LTR_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t6.8\t6.5\t" + str(que) + "\tLTRdetailref_" + str(index_r) + "\totherLTR\n")
        elif group == "other_repeat":
            o.write(str(each_chr) + "\t" + str(LTR_id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t6.8\t6.5\t" + str(que) + "\tLTRdetailref_" + str(index_r) + "\tDNAtrans\n")
        else:
            pass

    for index_q, row_q in df_que_newID.iterrows():
        strand_id_q = str(que) + "_strand_" + str(index_q)
        LTR_id_q = str(que) + "_LTRdetail_" + str(index_q)
        s_q = row_q["s_que"]
        e_q = row_q["e_que"]
        s_new_q = s_q - (df_que_newID["s_que"].min())
        e_new_q = e_q - (df_que_newID["s_que"].min())
        strand_q = row_q["strand_que"]
        group_q = row_q["Group"]
        if group_q != "LTR" and group_q != "other_repeat" and group_q != "Gene" and group_q != "Gap":
            if strand_q == "+":
                o.write(str(each_chr) + "\t" + str(strand_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t1.1\t1.4\t" + str(que) + "\tstrandque_" + str(index_q) + "\tplus\n")
            elif strand_q == "-":
                o.write(str(each_chr) + "\t" + str(strand_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t1.1\t1.4\t" + str(que) + "\tstrandque_" + str(index_q) + "\tminus\n")
            else:
                pass
                # o.write(str(each_chr) + "\t" + str(strand_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t1.1\t1.4\t" + str(que) + "\tstrandque_" + str(index_q) + "\tequal\n")
        elif group_q == "LTR":
            ltrid_q = row_q["id_que"]
            if "SZ-22" in ltrid_q:
                o.write(str(each_chr) + "\t" + str(LTR_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t0.2\t0.5\t" + str(que) + "\tLTRdetailref_" + str(index_q) + "\tSZ-22\n")
            elif "RIRE" in ltrid_q:
                o.write(str(each_chr) + "\t" + str(LTR_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t0.2\t0.5\t" + str(que) + "\tLTRdetailref_" + str(index_q) + "\tRIRE\n")
            else:
                o.write(str(each_chr) + "\t" + str(LTR_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t0.2\t0.5\t" + str(que) + "\tLTRdetailref_" + str(index_q) + "\totherLTR\n")
        elif group_q == "other_repeat":
            o.write(str(each_chr) + "\t" + str(LTR_id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t0.2\t0.5\t" + str(que) + "\tLTRdetailref_" + str(index_q) + "\tDNAtrans\n")
        else:
            pass


    #### black box
    s_box_r = (df_ref_newID["s_ref"].min()) - (df_ref_newID["s_ref"].min())
    e_box_r = (df_ref_newID["e_ref"].max()) - (df_ref_newID["s_ref"].min())
    o.write(str(each_chr) + "\tBOXid1\t" + str(s_box_r) + "\t" + str(e_box_r) + "\t5\t4.7\t" + str(ref) + "\tBOXref_1\tbox\n")
    o.write(str(each_chr) + "\tBOXid2\t" + str(s_box_r) + "\t" + str(e_box_r) + "\t5.9\t5.6\t" + str(ref) + "\tBOXref_2\tbox\n")
    o.write(str(each_chr) + "\tBOXid3\t" + str(s_box_r) + "\t" + str(e_box_r) + "\t6.8\t6.5\t" + str(ref) + "\tBOXref_3\tbox\n")
    s_box_q = (df_que_newID["s_que"].min()) - (df_que_newID["s_que"].min())
    e_box_q = (df_que_newID["e_que"].max()) - (df_que_newID["s_que"].min())
    o.write(str(each_chr) + "\tBOXid4\t" + str(s_box_q) + "\t" + str(e_box_q) + "\t2\t2.3\t" + str(que) + "\tBOXque_1\tbox\n")
    o.write(str(each_chr) + "\tBOXid5\t" + str(s_box_q) + "\t" + str(e_box_q) + "\t1.1\t1.4\t" + str(que) + "\tBOXque_2\tbox\n")
    o.write(str(each_chr) + "\tBOXid6\t" + str(s_box_q) + "\t" + str(e_box_q) + "\t0.2\t0.5\t" + str(que) + "\tBOXque_3\tbox\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-c', '--chr', help='chr', dest='chr', required=True)
    parser.add_argument('-e', '--ed', help='ED', dest='ED', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    add_black_region(args.ref, args.que, args.chr, args.ED, args.Dir)
    add_strand_track(args.ref, args.que, args.chr, args.ED, args.Dir)

if __name__ == '__main__':
    main()