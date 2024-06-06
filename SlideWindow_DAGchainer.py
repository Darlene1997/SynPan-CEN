import pandas as pd
import argparse
import os

def ED_pairnumber_pick(ref, que, each_chr, ED, Dir):
    ### get total elements number
    newID_ref = str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt"
    temp_r = str(each_chr) + "_" + str(ref) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_ref, each_chr, temp_r))
    MAX_order_s_r = count_lines(temp_r)
    os.system("rm %s" % (temp_r))

    ### search for the existing framework
    ED_frame_select_finial = 20
    for frame_range in range(0, 21):
        out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(frame_range) + "_framework.sg"
        if os.path.exists(out_sg) and os.path.getsize(out_sg):
            df_out_sg = pd.read_table(out_sg)
            out_sg_pairnumber = df_out_sg.shape[0]
            if (20 * out_sg_pairnumber > MAX_order_s_r) or (out_sg_pairnumber > 100):
                ED_frame_select_finial = int(frame_range)
                break

    ED_all = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg_filter"
    if os.path.exists(ED_all) and os.path.getsize(ED_all):
        pass
    else:
        os.system('echo "Do get_bigger_ED" >> %s_%s_%s_test.log' % (each_chr, ref, que))
        get_bigger_ED(ref, que, each_chr, ED, ED_frame_select_finial, Dir)

def get_bigger_ED(ref, que, each_chr, ED, ED_frame_select, Dir):
    ### get total elements number
    ## Reference
    newID_ref = str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt"
    temp_r = str(each_chr) + "_" + str(ref) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_ref, each_chr, temp_r))
    MAX_order_s_r = count_lines(temp_r)
    os.system("rm %s" % (temp_r))
    ##Query
    newID_que = str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt"
    temp_q = str(each_chr) + "_" + str(que) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_que, each_chr, temp_q))
    MAX_order_s_q = count_lines(temp_q)
    os.system("rm %s" % (temp_q))

    os.system('echo "MAX_order_s_r: %s , MAX_order_s_q: %s" >> %s_%s_%s_test.log' % (MAX_order_s_r, MAX_order_s_q, each_chr, ref, que))


    out_sg_ED_frame = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework.sg"
    df_sg = pd.read_table(out_sg_ED_frame, sep="\t", header=None)
    df_sg.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue", "score"]
    ### sort by the reference. e.g. NIP
    if df_sg.loc[0, "chr_1"].split("_")[1] == "NIP":
        df_sg_sort = df_sg.sort_values(by="s_1").reset_index(drop=True)
        #### Using slide the window for partition, fix 20 pairs per window, and then use larger ED to run DAGchainer
        pair_num = df_sg_sort.shape[0]
        os.system('echo "chr_1--pair_num: %s" >> %s_%s_%s_test.log' % (pair_num, each_chr, ref, que))
        pair_now = 0
        window = 20
        block_number = 0
        while pair_now <= pair_num:
            os.system('echo "pair_now: %s" >> %s_%s_%s_test.log' % (pair_now, each_chr, ref, que))
            if block_number == 0:
                os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                min_s_1 = 0
                max_s_1 = df_sg_sort.loc[0, "s_1"]
                min_s_2 = 0
                max_s_2 = df_sg_sort.loc[0, "s_2"]
                os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                block_number += 1
            else:
                block_s = pair_now
                block_e = pair_now + window
                if block_e < pair_num:
                    os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    subset_df = df_sg_sort.iloc[block_s:block_e+1]
                    min_s_1 = subset_df['s_1'].min()
                    max_s_1 = subset_df['s_1'].max()
                    min_s_2 = subset_df['s_2'].min()
                    max_s_2 = subset_df['s_2'].max()
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    block_number += 1
                    pair_now += window
                elif pair_now == pair_num:
                    ### from the last block to the last monomer
                    os.system('echo "block_number: %s (same)" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    min_s_1 = df_sg_sort.iloc[-1]["s_1"]
                    max_s_1 = MAX_order_s_r
                    min_s_2 = df_sg_sort.iloc[-1]["s_2"]
                    max_s_2 = MAX_order_s_q
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    pair_now += window
                else:
                    ### the last block
                    os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    subset_df = df_sg_sort.iloc[block_s:pair_num+1]
                    min_s_1 = subset_df['s_1'].min()
                    max_s_1 = subset_df['s_1'].max()
                    min_s_2 = subset_df['s_2'].min()
                    max_s_2 = subset_df['s_2'].max()
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    block_number += 1
                    pair_now += window

                    ### from the last block to the last monomer
                    min_s_1 = df_sg_sort.iloc[-1]["s_1"]
                    max_s_1 = MAX_order_s_r
                    min_s_2 = df_sg_sort.iloc[-1]["s_2"]
                    max_s_2 = MAX_order_s_q
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)


        all_DAG_out = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg"
        os.system("cat %s | cut -f1,2,3,4,5,6,7,8 | sort -k5,5 -k7,7n > %s_%s_%s_pairED%s.sg_filter" % (all_DAG_out, each_chr, ref, que, ED))#### The last step is here !!
        out_count = str(ref) + "_" + str(que) + "_count_ED_pairs.txt"
        with open(out_count, 'a') as outtf:
            df_filter = pd.read_table(all_DAG_out)
            df_without_duplicates = df_filter.drop_duplicates()
            pair_count = df_without_duplicates.shape[0]
            outtf.write(str(ref) + "\t" + str(que) + "\t" + str(each_chr) + "\t" + str(ED) + "\t" + str(pair_count) + "\n")
        os.system('echo "#####################  Finish geting %s_%s_%s_pairED%s.sg_filter  #######################" >> %s_%s_%s_test.log' % (each_chr, ref, que, ED, each_chr, ref, que))
    else:
        df_sg_sort = df_sg.sort_values(by="s_2").reset_index(drop=True)
        #### Using slide the window for partition, fix 20 pairs per window, and then use larger ED to run DAGchainer
        pair_num = df_sg_sort.shape[0]
        os.system('echo "chr_2--pair_num: %s" >> %s_%s_%s_test.log' % (pair_num, each_chr, ref, que))
        pair_now = 0
        window = 20
        block_number = 0
        while pair_now <= pair_num:
            os.system('echo "pair_now: %s" >> %s_%s_%s_test.log' % (pair_now, each_chr, ref, que))
            if block_number == 0:
                os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                min_s_1 = 0
                max_s_1 = df_sg_sort.loc[0, "s_2"]
                min_s_2 = 0
                max_s_2 = df_sg_sort.loc[0, "s_1"]
                os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                block_number += 1
            else:
                block_s = pair_now
                block_e = pair_now + window
                if block_e <= pair_num:
                    os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    subset_df = df_sg_sort.iloc[block_s:block_e+1]
                    min_s_1 = subset_df['s_2'].min()
                    max_s_1 = subset_df['s_2'].max()
                    min_s_2 = subset_df['s_1'].min()
                    max_s_2 = subset_df['s_1'].max()
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    block_number += 1
                    pair_now += window
                elif pair_now == pair_num:
                    ### from the last block to the last monomer
                    os.system('echo "block_number: %s (same)" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    min_s_1 = df_sg_sort.iloc[-1]["s_1"]
                    max_s_1 = MAX_order_s_r
                    min_s_2 = df_sg_sort.iloc[-1]["s_2"]
                    max_s_2 = MAX_order_s_q
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    pair_now += window
                else:
                    ### the last block
                    os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    subset_df = df_sg_sort.iloc[block_s:pair_num+1]
                    min_s_1 = subset_df['s_2'].min()
                    max_s_1 = subset_df['s_2'].max()
                    min_s_2 = subset_df['s_1'].min()
                    max_s_2 = subset_df['s_1'].max()
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
                    block_number += 1
                    pair_now += window

                    ### from the last block to the last monomer
                    os.system('echo "block_number: %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                    min_s_1 = df_sg_sort.iloc[-1]["s_2"]
                    max_s_1 = MAX_order_s_r
                    min_s_2 = df_sg_sort.iloc[-1]["s_1"]
                    max_s_2 = MAX_order_s_q
                    os.system('echo "Do out_bigger_ED_for_panduan" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    os.system('echo "Ref: %s-%s; Que: %s-%s" >> %s_%s_%s_test.log' % (min_s_1, max_s_1, min_s_2, max_s_2, each_chr, ref, que))
                    out_bigger_ED_for_panduan(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)



        all_DAG_out = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + ".sg"
        os.system("cat %s | cut -f1,2,3,4,5,6,7,8 | sort -k5,5 -k7,7n > %s_%s_%s_pairED%s.sg_filter" % (all_DAG_out, each_chr, ref, que, ED))  #### The last step is here !!
        out_count = str(ref) + "_" + str(que) + "_count_ED_pairs.txt"
        with open(out_count, 'a') as outtf:
            df_filter = pd.read_table(all_DAG_out)
            df_without_duplicates = df_filter.drop_duplicates()
            pair_count = df_without_duplicates.shape[0]
            outtf.write(str(ref) + "\t" + str(que) + "\t" + str(each_chr) + "\t" + str(ED) + "\t" + str(pair_count) + "\n")
        os.system('echo "#####################  Finish geting %s_%s_%s_pairED%s.sg_filter  #######################" >> %s_%s_%s_test.log' % (each_chr, ref, que, ED, each_chr, ref, que))

def out_bigger_ED_for_panduan(ref, que, each_chr, i, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir):
    all_DAG_out = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(i) + ".sg"
    in_re_format = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"
    out_bigger_ED = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(i) + "_block" + str(block_number)
    with open(in_re_format, 'r') as inf:
        with open(out_bigger_ED, 'w') as outf:
            for line in inf:  ## "chr_r", "id_r", "s_r", "e_r", "chr_q", "id_q", "s_q", "e_q", "prop", "ED"
                text = line.strip().split("\t")
                chr_r = text[0]
                id_r = text[1]
                index_r = text[2]
                chr_q = text[4]
                id_q = text[5]
                index_q = text[6]
                prop = text[8]
                ED_for_check = text[9]
                if (int(index_r) >= int(min_s_1)) and (int(index_r) <= int(max_s_1)) and (int(index_q) >= int(min_s_2)) and (int(index_q) <= int(max_s_2)) and (int(ED_for_check) <= int(i)):
                    outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\n")

    os.system('echo "Do get_DAGchainer_synteny" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    get_DAGchainer_synteny(ref, que, each_chr, i, block_number, Dir)
    DAG_out = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(i) + "_block" + str(block_number) + ".sg"
    os.system("cat %s >> %s" % (DAG_out, all_DAG_out))
    os.system("rm %s" % (out_bigger_ED))
    os.system("rm %s.sg" % (out_bigger_ED))
    os.system("rm %s.aligncoords" % (out_bigger_ED))

def get_DAGchainer_synteny(ref, que, each_chr, ED, block_number, Dir):
    out_bigger_ED = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number)
    ##### split again for some big block
    threshold = 30 * 1024 * 1024  # 30MB
    file_size = os.path.getsize(out_bigger_ED)
    if file_size > threshold:
        os.system("python ../03_supp_split_block.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
        out_bigger_ED_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1"
        out_bigger_ED_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2"
        DAG_out = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + ".sg"
        out_sg_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1.sg"
        out_sg_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2.sg"
        os.system("cat %s >> %s" % (out_sg_1, DAG_out))
        os.system("cat %s >> %s" % (out_sg_2, DAG_out))
        os.system("rm %s %s" % (out_bigger_ED_1, out_bigger_ED_2))
        os.system("rm %s.sg %s.sg" % (out_bigger_ED_1, out_bigger_ED_2))
        os.system("rm %s.aligncoords %s.aligncoords" % (out_bigger_ED_1, out_bigger_ED_2))
        out_bigger_ED_3 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_3"
        out_bigger_ED_4 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str( block_number) + "_4"
        out_sg_3 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_3.sg"
        out_sg_4 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_4.sg"
        if os.path.exists(out_bigger_ED_3) or os.path.exists(out_bigger_ED_4):
            os.system("cat %s >> %s" % (out_sg_3, DAG_out))
            os.system("cat %s >> %s" % (out_sg_4, DAG_out))
            os.system("rm %s %s" % (out_bigger_ED_3, out_bigger_ED_4))
            os.system("rm %s.sg %s.sg" % (out_bigger_ED_3, out_bigger_ED_4))
            os.system("rm %s.aligncoords %s.aligncoords" % (out_bigger_ED_3, out_bigger_ED_4))
    else:
        DAG_block_file = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + ".aligncoords"
        out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + ".sg"
        if os.path.exists(out_sg) and os.path.getsize(out_sg):
            pass
        else:
            if os.path.exists(DAG_block_file) and os.path.getsize(DAG_block_file):
                pass
            else:
                os.system('echo "DAGchainer for block %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
                os.system("perl ./software/dagchainer/bin/run_DAG_chainer.pl -i %s -g 1 -D 20 -A 6 -Z 50 -e -3f" % (out_bigger_ED))  ## Specify your PATH to DAGCHAINER here

        sz = os.path.getsize(DAG_block_file)
        if not sz:
            print(str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + ".aligncoords --- It's empty!")
        else:
            os.system('echo "Start geting pairED%s_block%s.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))
            add_each_block_without_overlap(ref, que, each_chr, ED, block_number)
            os.system('echo "Finish geting pairED%s_block%s.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))


def count_lines(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file)

def add_each_block_without_overlap(ref, que, each_chr, ED, block_number):
    DAG_block_file = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str( block_number) + ".aligncoords"
    out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + ".sg"
    ###### !! Important: Add non-overlapping results in block order
    block_data = []
    block_num = 0
    position_dict_ref = {}
    position_dict_que = {}
    with open(DAG_block_file, 'r') as inf:
        with open(out_sg, 'w') as outf:
            for lines in inf:
                line = lines.strip()
                if line.startswith('#'):
                    if block_data:
                        block_num += 1
                        if block_num == 1:
                            df_block_data = pd.DataFrame([x.strip().split('\t') for x in block_data])
                            df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2","evalue", "score"]
                            df_block_data['s_1'] = df_block_data['s_1'].astype(int)
                            df_block_data['s_2'] = df_block_data['s_2'].astype(int)
                            min_s_1 = df_block_data['s_1'].min()
                            max_s_1 = df_block_data['s_1'].max()
                            min_s_2 = df_block_data['s_2'].min()
                            max_s_2 = df_block_data['s_2'].max()
                            position_dict_ref[block_num] = [min_s_1, max_s_1]
                            position_dict_que[block_num] = [min_s_2, max_s_2]
                            for lines_block in block_data:
                                line_block = lines_block.strip()
                                outf.write(line_block + "\n")

                        else:
                            df_block_data = pd.DataFrame([x.strip().split('\t') for x in block_data])
                            df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2","evalue", "score"]
                            df_block_data['s_1'] = df_block_data['s_1'].astype(int)
                            df_block_data['s_2'] = df_block_data['s_2'].astype(int)
                            find_min_s_1 = []
                            find_max_s_1 = []
                            find_min_s_2 = []
                            find_max_s_2 = []
                            for index, row in df_block_data.iterrows():
                                s_1 = row["s_1"]
                                s_2 = row["s_2"]
                                n_not_in_count_1 = 0
                                n_not_in_count_2 = 0
                                for block_number1, values1 in position_dict_ref.items():
                                    if s_1 not in range(values1[0], values1[1] + 1):
                                        n_not_in_count_1 += 1
                                for block_number2, values2 in position_dict_que.items():
                                    if s_2 not in range(values2[0], values2[1] + 1):
                                        n_not_in_count_2 += 1
                                if n_not_in_count_1 == len(position_dict_ref) and n_not_in_count_2 == len(
                                        position_dict_que):
                                    row_string = '\t'.join(map(str, row))
                                    outf.write(row_string + '\n')
                                    find_min_s_1.append(s_1)
                                    find_max_s_1.append(s_1)
                                    find_min_s_2.append(s_2)
                                    find_max_s_2.append(s_2)
                            if len(find_min_s_1) > 0:
                                position_dict_ref[block_num] = [min(find_min_s_1), max(find_min_s_1)]
                                position_dict_que[block_num] = [min(find_min_s_2), max(find_min_s_2)]
                        block_data = []
                else:
                    block_data.append(line)

            if block_data:
                df_block_data = pd.DataFrame([x.strip().split('\t') for x in block_data])
                df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue", "score"]
                df_block_data['s_1'] = df_block_data['s_1'].astype(int)
                df_block_data['s_2'] = df_block_data['s_2'].astype(int)
                for index, row in df_block_data.iterrows():
                    s_1 = row["s_1"]
                    s_2 = row["s_2"]
                    n_not_in_count_1 = 0
                    n_not_in_count_2 = 0
                    for block_number1, values1 in position_dict_ref.items():
                        if s_1 not in range(values1[0], values1[1] + 1):
                            n_not_in_count_1 += 1
                    for block_number2, values2 in position_dict_que.items():
                        if s_2 not in range(values2[0], values2[1] + 1):
                            n_not_in_count_2 += 1
                    if n_not_in_count_1 == len(position_dict_ref) and n_not_in_count_2 == len(position_dict_que):
                        row_string = '\t'.join(map(str, row))
                        outf.write(row_string + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-chr', '--chr', help='Chr', dest='chr', required=True)
    parser.add_argument('-e', '--ed', help='ED', dest='ED', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    ED_pairnumber_pick(args.ref, args.que, args.chr, args.ED, args.Dir)


if __name__ == '__main__':
    main()
