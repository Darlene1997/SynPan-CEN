import pandas as pd
import argparse
import os

def merge_for_plot(ref, que, each_chr, ED, Dir):
    os.system('echo "Start merge_for_plot" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    ####################### BEGIN to match ##################
    df_ref_newID = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["chr_ref", "s_ref", "e_ref", "strand_ref", "id_ref", "Group"])
    df_que_newID = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt", sep="\t", header=None, names=["chr_que", "s_que", "e_que", "strand_que", "id_que", "Group"])
    df_que_newID_chr = df_que_newID[df_que_newID["chr_que"] == str(each_chr) + "_" + str(que)].reset_index(drop=True)
    df_ref_newID_chr = df_ref_newID[df_ref_newID["chr_ref"] == str(each_chr) + "_" + str(ref)].reset_index(drop=True)

    #### Put all query ids in a list for subsequent screening
    list_que_id = []
    dict_id_s_que = {}
    dict_id_e_que = {}
    for n, rr in df_que_newID_chr.iterrows():
        group = rr["Group"]
        s_que = rr["s_que"]
        e_que = rr["e_que"]
        if group == "LTR" or group == "other_repeat":
            id_que = str(rr["chr_que"]) + ".LTR." + str(n)
        else:
            id_que = rr["id_que"]

        list_que_id.append(id_que)
        dict_id_s_que[id_que] = s_que
        dict_id_e_que[id_que] = e_que

    #### Extract all sg pairs into the corresponding dictionary
    dict_id_match = {}
    ED_all = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg_filter"
    if os.path.exists(ED_all):
        if not os.path.getsize(ED_all):
            print(str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg_filter --- It's empty!")
        else:
            df_merge = pd.read_table(ED_all, sep="\t", header=None)
            df_merge.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2"]
            if df_merge.loc[0, "chr_1"] == str(each_chr) + "_" + str(ref):
                df_merge_select = df_merge[["id_1", "s_1", "e_1", "id_2", "s_2", "e_2"]]
                df_merge_select.columns = ["id_ref", "s_ref", "e_ref", "id_" + str(que), "s_" + str(que), "e_" + str(que)]
                df_merge_uniq = df_merge_select.drop_duplicates(subset=['id_ref'], keep='first', ignore_index=True)
                df_merge_uniq2 = df_merge_uniq.drop_duplicates(subset=["id_" + str(que)], keep='first', ignore_index=True)
                for i, r in df_merge_uniq2.iterrows():
                    id_r = r["id_ref"]
                    id_q = r["id_" + str(que)]
                    dict_id_match[id_r] = id_q
            else:
                df_merge_select = df_merge[["id_2", "s_2", "e_2", "id_1", "s_1", "e_1"]]
                df_merge_select.columns = ["id_ref", "s_ref", "e_ref", "id_" + str(que), "s_" + str(que), "e_" + str(que)]
                #### Only the one-to-one pairs of the merged files are reserved
                df_merge_uniq = df_merge_select.drop_duplicates(subset=['id_ref'], keep='first', ignore_index=True)
                df_merge_uniq2 = df_merge_uniq.drop_duplicates(subset=["id_" + str(que)], keep='first', ignore_index=True)
                for iii, rrr in df_merge_uniq2.iterrows():
                    id_r = rrr["id_ref"]
                    id_q = rrr["id_" + str(que)]
                    dict_id_match[id_r] = id_q
    else:
        print(str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg_filter --- It's not exist!")

    #### Get the desired format for the drawing
    out_for_plot = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_for_plot_ED" + str(ED) + ".txt"
    o = open(out_for_plot, 'w')

    os.system('echo " 1 matched ids" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    ############  matched ids
    for index_r, row_r in df_ref_newID_chr.iterrows():
        chr = row_r["chr_ref"].split("_")[0]  ##Chr01
        s_r = row_r["s_ref"]
        e_r = row_r["e_ref"]
        group = row_r["Group"]
        if group == "LTR" or group == "other_repeat":
            id_r = str(row_r["chr_ref"]) + ".LTR." + str(index_r)
        else:
            id_r = row_r["id_ref"]
        s_new_r = s_r - (df_ref_newID_chr["s_ref"].min())
        e_new_r = e_r - (df_ref_newID_chr["s_ref"].min())
        if str(id_r).find("LTR") == -1:  ### If LTR does not exist in the ref id, -1 is returned
            if str(id_r).find(each_chr) != -1:
                group_r1 = str(id_r).split(".")[0]
                group_r2 = str(group_r1).split("_")[3]
                o.write(str(chr) + "\t" + str(id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t5\t4.7\t" + str(ref) + "\tpair" + str(index_r) + "\t" + str(group_r2) + "\n")
        else:
            if str(id_r).find("LTR") != -1:
                group_r2 = "LTR"
                o.write(str(chr) + "\t" + str(id_r) + "\t" + str(s_new_r) + "\t" + str(e_new_r) + "\t5\t4.7\t" + str(ref) + "\tpair" + str(index_r) + "\t" + str(group_r2) + "\n")

        if id_r in dict_id_match.keys():
            id_q = dict_id_match[id_r]
            s_q = dict_id_s_que[id_q]
            e_q = dict_id_e_que[id_q]
            s_new_q = s_q - (df_que_newID_chr["s_que"].min())
            e_new_q = e_q - (df_que_newID_chr["s_que"].min())
            if str(id_q).find("LTR") == -1:  ### If LTR does not exist in the ref id, -1 is returned
                if str(id_q).find(each_chr) != -1:
                    group_q1 = str(id_q).split(".")[0]
                    group_q2 = str(group_q1).split("_")[3]
                    o.write(str(chr) + "\t" + str(id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t2\t2.3\t" + str(que) + "\tpair" + str(index_r) + "\t" + str(group_q2) + "\n")
                if id_q in list_que_id:
                    list_que_id.remove(id_q)
            else:
                if str(id_q).find("LTR") != -1:
                    group_q2 = "LTR"
                    o.write(str(chr) + "\t" + str(id_q) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t2\t2.3\t" + str(que) + "\tpair" + str(index_r) + "\t" + str(group_q2) + "\n")
                if id_q in list_que_id:
                    list_que_id.remove(id_q)

    ############ no matched ids
    os.system('echo " 2 no matched ids" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    npair = df_ref_newID_chr.shape[0]
    for each in list_que_id:
        chr = str(each).split("_")[0]   ##Chr01
        s_q = dict_id_s_que[each]
        e_q = dict_id_e_que[each]
        s_new_q = s_q - (df_que_newID_chr["s_que"].min())
        e_new_q = e_q - (df_que_newID_chr["s_que"].min())
        if str(each).find("LTR") == -1:  ### If LTR does not exist in the ref id, -1 is returned
            if str(each).find(each_chr) != -1:
                group_q1 = str(each).split(".")[0]
                group_q2 = str(group_q1).split("_")[3]
                o.write(str(chr) + "\t" + str(each) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t2\t2.3\t" + str(que) + "\tpair" + str(npair) + "\t" + str(group_q2) + "\n")
                npair += 1
        else:
            if str(each).find("LTR") != -1:
                group_q2 = "LTR"
                o.write(str(chr) + "\t" + str(each) + "\t" + str(s_new_q) + "\t" + str(e_new_q) + "\t2\t2.3\t" + str(que) + "\tpair" + str(npair) + "\t" + str(group_q2) + "\n")
                npair += 1

    ############ genes
    os.system('echo " 3 gene" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    df_centromere = pd.read_table(str(Dir) + "01_CR_70m.bed", sep="\t")  ## Chr	start	end
    df_cent_r = df_centromere[df_centromere["Chr"] == str(each_chr) + "_" + str(ref)].reset_index(drop=True)
    df_cent_q = df_centromere[df_centromere["Chr"] == str(each_chr) + "_" + str(que)].reset_index(drop=True)
    if df_cent_r.shape[0] != 0:
        cent_s_r = df_cent_r.loc[0, "start"]
        cent_e_r = df_cent_r.loc[0, "end"]
        df_gene_r = pd.read_table("all_asm_70_CR_Msu7_clean_filtered3.gff3", comment="#", sep="\t", header=None, names=["chr", "msu7", "type", "start", "end", "l", "m", "n", "anno"])
        df_ref_centro_gene = df_gene_r[(df_gene_r["chr"] == str(each_chr) + "_" + str(ref) + "_centro") & (df_gene_r["type"] == "gene") & (df_gene_r["start"] >= cent_s_r) & (df_gene_r["end"] <= cent_e_r)].reset_index(drop=True)

        if df_ref_centro_gene.shape[0] != 0:
            for index_gene_r, row_gene_r in df_ref_centro_gene.iterrows():
                start_gene_r = row_gene_r["start"] - (df_ref_newID_chr["s_ref"].min())
                end_gene_r = row_gene_r["end"] - (df_ref_newID_chr["s_ref"].min())
                gene_id_r = str(ref) + "_gene_" + str(index_gene_r)
                o.write(str(each_chr) + "\t" + str(gene_id_r) + "\t" + str(start_gene_r) + "\t" + str(end_gene_r) + "\t5\t4.7\t" + str(ref) + "\tpair_refgene_" + str(index_gene_r) + "\tgene\n")

    if df_cent_q.shape[0] != 0:
        cent_s_q = df_cent_q.loc[0, "start"]
        cent_e_q = df_cent_q.loc[0, "end"]
        df_gene_q = pd.read_table("all_asm_70_CR_Msu7_clean_filtered3.gff3", comment="#", sep="\t", header=None, names=["chr", "msu7", "type", "start", "end", "l", "m", "n", "anno"])
        df_que_centro_gene = df_gene_q[(df_gene_q["chr"] == str(each_chr) + "_" + str(que) + "_centro") & (df_gene_q["type"] == "gene") & (df_gene_q["start"] >= cent_s_q) & (df_gene_q["end"] <= cent_e_q)].reset_index(drop=True)
        if df_que_centro_gene.shape[0] != 0:
            for index_gene_q, row_gene_q in df_que_centro_gene.iterrows():
                start_gene_q = row_gene_q["start"] - (df_que_newID_chr["s_que"].min())
                end_gene_q = row_gene_q["end"] - (df_que_newID_chr["s_que"].min())
                gene_id_q = str(que) + "_gene_" + str(index_gene_q)
                o.write(str(each_chr) + "\t" + str(gene_id_q) + "\t" + str(start_gene_q) + "\t" + str(end_gene_q) + "\t2\t2.3\t" + str(que) + "\tpair_quegene_" + str(index_gene_q) + "\tgene\n")

    ############ gaps
    os.system('echo " 4 gaps" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    df_gap = pd.read_table("gap", sep="\t", header=None, names=["Chr", "start", "end"])
    if df_cent_r.shape[0] != 0:
        cent_s_r = df_cent_r.loc[0, "start"]
        cent_e_r = df_cent_r.loc[0, "end"]
        df_ref_gap = df_gap[(df_gap["Chr"] == str(each_chr) + "_" + str(ref)) & (df_gap["start"] >= cent_s_r) & (df_gap["end"] <= cent_e_r)].reset_index(drop=True)
        if df_ref_gap.shape[0] != 0:
            for index_gap_r, row_gap_r in df_ref_gap.iterrows():
                gap_s_r = row_gap_r["start"] - (df_ref_newID_chr["s_ref"].min())
                gap_e_r = row_gap_r["end"] - (df_ref_newID_chr["s_ref"].min())
                gap_id_r = str(ref) + "_gap_" + str(index_gap_r)
                if gap_s_r > 0:
                    o.write(str(each_chr) + "\t" + str(gap_id_r) + "\t" + str(gap_s_r) + "\t" + str(gap_e_r) + "\t5.3\t4.7\t" + str(ref) + "\tpair_gap_" + str(index_gap_r) + "\tgap\n")

    if df_cent_q.shape[0] != 0:
        cent_s_q = df_cent_q.loc[0, "start"]
        cent_e_q = df_cent_q.loc[0, "end"]
        df_que_gap = df_gap[(df_gap["Chr"] == str(each_chr) + "_" + str(que)) & (df_gap["start"] >= cent_s_q) & (df_gap["end"] <= cent_e_q)].reset_index(drop=True)
        if df_que_gap.shape[0] != 0:
            for index_gap_q, row_gap_q in df_que_gap.iterrows():
                gap_s_q = row_gap_q["start"] - (df_que_newID_chr["s_que"].min())
                gap_e_q = row_gap_q["end"] - (df_que_newID_chr["s_que"].min())
                gap_id_q = str(ref) + "_gap_" + str(index_gap_q)
                if gap_s_q > 0:
                    o.write(str(each_chr) + "\t" + str(gap_id_q) + "\t" + str(gap_s_q) + "\t" + str(gap_e_q) + "\t2.3\t2.3\t" + str(ref) + "\tpair_gap_" + str(index_gap_q) + "\tgap\n")

    ############ blank region
    os.system('echo " 5 blank region" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    line_total = df_ref_newID_chr.shape[0] - 1
    for index, row in df_ref_newID_chr.iterrows():
        End = row["e_ref"]
        if index < line_total:
            nest_index = int(index) + 1
            Start_next = df_ref_newID_chr.loc[nest_index, "s_ref"]
            if int(End) + 1 == int(Start_next):
                pass
            else:
                blank_s_r = End - (df_ref_newID_chr["s_ref"].min())
                blank_e_r = Start_next - (df_ref_newID_chr["s_ref"].min())
                blank_id_r = str(ref) + "_blank_" + str(index)
                o.write(str(each_chr) + "\t" + str(blank_id_r) + "\t" + str(blank_s_r) + "\t" + str(blank_e_r) + "\t5\t4.7\t" + str(ref) + "\tpair_refblank_" + str(index) + "\tblank\n")

    line_total2 = df_que_newID_chr.shape[0] - 1
    for index2, row2 in df_que_newID_chr.iterrows():
        End2 = row2["e_que"]
        if index2 < line_total2:
            nest_index2 = int(index2) + 1
            Start_next2 = df_que_newID_chr.loc[nest_index2, "s_que"]
            if int(End2) + 1 == int(Start_next2):
                pass
            else:
                blank_s_q = End2 - (df_que_newID_chr["s_que"].min())
                blank_e_q = Start_next2 - (df_que_newID_chr["s_que"].min())
                blank_id_q = str(que) + "_blank_" + str(index2)
                o.write(str(each_chr) + "\t" + str(blank_id_q) + "\t" + str(blank_s_q) + "\t" + str(blank_e_q) + "\t2\t2.3\t" + str(que) + "\tpair_queblank_" + str(index2) + "\tblank\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-c', '--chr', help='chr', dest='chr', required=True)
    parser.add_argument('-e', '--ed', help='ED', dest='ED', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    merge_for_plot(args.ref, args.que, args.chr, args.ED, args.Dir)

if __name__ == '__main__':
    main()
