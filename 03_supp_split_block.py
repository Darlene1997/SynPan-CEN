import os
import pandas as pd
import argparse

def split_block_add(ref, que, each_chr, ED, block_number, Dir):
    os.system('echo "%s_%s_%s_pairED%s_block%s is too big !!! Split it !!!" >> %s_%s_%s_test.log' % (each_chr, ref, que, ED, block_number, each_chr, ref, que))
    ### get total elements number
    newID_ref = str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt"
    temp_r = str(each_chr) + "_" + str(ref) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_ref, each_chr, temp_r))
    MAX_order_s_r = count_lines(temp_r)
    os.system("rm %s" % (temp_r))
    newID_que = str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt"
    temp_q = str(each_chr) + "_" + str(que) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_que, each_chr, temp_q))
    MAX_order_s_q = count_lines(temp_q)
    os.system("rm %s" % (temp_q))

    ### search for the existing framework
    ED_frame_select_finial = 0
    for ED_frame_select in range(0, 20):
        out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework.sg"
        if os.path.exists(out_sg) and os.path.getsize(out_sg):
            df_out_sg = pd.read_table(out_sg)
            out_sg_pairnumber = df_out_sg.shape[0]
            if 20 * out_sg_pairnumber > MAX_order_s_r:
                ED_frame_select_finial = int(ED_frame_select)
                break

    out_sg_ED_frame = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select_finial) + "_framework.sg"
    df_sg = pd.read_table(out_sg_ED_frame, sep="\t", header=None)
    df_sg.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue", "score"]
    threshold = 20 * 1024 * 1024  # 20MB
    block_s = (int(block_number) - 1) * 20
    block_e = int(block_number) * 20
    out_bigger_ED_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1"
    out_bigger_ED_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2"
    if df_sg.loc[0, "chr_1"].split("_")[1] == "NIP":
        df_sg_sort = df_sg.sort_values(by="s_1").reset_index(drop=True)
        if int(block_number) == 0:
            min_s_1 = 0
            max_s_1 = df_sg_sort.loc[0, "s_1"]
            min_s_2 = 0
            max_s_2 = df_sg_sort.loc[0, "s_2"]
            get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

            file_size_1 = os.path.getsize(out_bigger_ED_1)
            file_size_2 = os.path.getsize(out_bigger_ED_2)
            if (file_size_1 > threshold) or (file_size_2 > threshold):
                os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
            else:
                get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

        else:
            if block_e <= df_sg.shape[0]:
                subset_df = df_sg_sort.iloc[block_s:block_e + 1]
                min_s_1 = subset_df['s_1'].min()
                max_s_1 = subset_df['s_1'].max()
                min_s_2 = subset_df['s_2'].min()
                max_s_2 = subset_df['s_2'].max()
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

            elif (block_s <= df_sg.shape[0]) and (block_e > df_sg.shape[0]):
                subset_df = df_sg_sort.iloc[block_s:df_sg.shape[0] + 1]
                min_s_1 = subset_df['s_1'].min()
                max_s_1 = subset_df['s_1'].max()
                min_s_2 = subset_df['s_2'].min()
                max_s_2 = subset_df['s_2'].max()
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

            else:
                min_s_1 = df_sg_sort.iloc[-1]["s_1"]
                max_s_1 = MAX_order_s_r
                min_s_2 = df_sg_sort.iloc[-1]["s_2"]
                max_s_2 = MAX_order_s_q
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

    else:
        df_sg_sort = df_sg.sort_values(by="s_2").reset_index(drop=True)
        if int(block_number) == 0:
            min_s_1 = 0
            max_s_1 = df_sg_sort.loc[0, "s_2"]
            min_s_2 = 0
            max_s_2 = df_sg_sort.loc[0, "s_1"]

            get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)
            threshold = 10 * 1024 * 1024
            file_size_1 = os.path.getsize(out_bigger_ED_1)
            file_size_2 = os.path.getsize(out_bigger_ED_2)
            if (file_size_1 > threshold) or (file_size_2 > threshold):
                os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
            else:
                get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)
        else:
            if block_e <= df_sg.shape[0]:
                subset_df = df_sg_sort.iloc[block_s:block_e + 1]
                min_s_1 = subset_df['s_2'].min()
                max_s_1 = subset_df['s_2'].max()
                min_s_2 = subset_df['s_1'].min()
                max_s_2 = subset_df['s_1'].max()
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

            elif (block_s <= df_sg.shape[0]) and (block_e > df_sg.shape[0]):
                subset_df = df_sg_sort.iloc[block_s:df_sg.shape[0] + 1]
                min_s_1 = subset_df['s_2'].min()
                max_s_1 = subset_df['s_2'].max()
                min_s_2 = subset_df['s_1'].min()
                max_s_2 = subset_df['s_1'].max()
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)

            else:
                min_s_1 = df_sg_sort.iloc[-1]["s_2"]
                max_s_1 = MAX_order_s_r
                min_s_2 = df_sg_sort.iloc[-1]["s_1"]
                max_s_2 = MAX_order_s_q
                get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1, min_s_2, max_s_2, Dir)

                threshold = 10 * 1024 * 1024
                file_size_1 = os.path.getsize(out_bigger_ED_1)
                file_size_2 = os.path.getsize(out_bigger_ED_2)
                if (file_size_1 > threshold) or (file_size_2 > threshold):
                    os.system("python ../03_supp_split_block2.py -r %s -q %s -chr %s -e %s -b %s -d %s" % (ref, que, each_chr, ED, block_number, Dir))
                else:
                    get_DAGchainer_synteny(ref, que, each_chr, ED, block_number)


def get_block(ref, que, each_chr, ED, block_number, min_s_1, max_s_1,  min_s_2, max_s_2, Dir):
    in_re_format = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"
    out_bigger_ED_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1"
    out_bigger_ED_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2"

    block_middle_1 = int((min_s_1 + max_s_1) / 2)
    block_middle_2 = int((min_s_2 + max_s_2) / 2)
    os.system('echo "block_number: %s; middle 1: %s" >> %s_%s_%s_test.log' % (block_number, block_middle_1, each_chr, ref, que))
    with open(in_re_format, 'r') as inf:
        with open(out_bigger_ED_1, 'w') as outf:
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
                if (int(index_r) >= int(min_s_1)) and (int(index_r) <= int(block_middle_1)) and (int(index_q) >= int(min_s_2)) and (int(index_q) <= int(block_middle_2)) and (int(ED_for_check) <= int(ED)):
                    outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\n")

    os.system('echo "block_number: %s; middle 2: %s" >> %s_%s_%s_test.log' % (block_number, block_middle_2, each_chr, ref, que))
    with open(in_re_format, 'r') as inf:
        with open(out_bigger_ED_2, 'w') as outf:
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
                if (int(index_r) >= int(block_middle_1)) and (int(index_r) <= int(max_s_1)) and (int(index_q) >= int(block_middle_2)) and (int(index_q) <= int(max_s_2)) and (int(ED_for_check) <= int(ED)):
                    outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\n")


def get_DAGchainer_synteny(ref, que, each_chr, ED, block_number):
    out_bigger_ED_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1"
    DAG_block_file_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1.aligncoords"
    out_sg_1 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1.sg"
    out_bigger_ED_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2"
    DAG_block_file_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2.aligncoords"
    out_sg_2 = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2.sg"
    ########################## the first block ########################
    if os.path.exists(out_sg_1) and os.path.getsize(out_sg_1):
        pass
    else:
        if os.path.exists(DAG_block_file_1) and os.path.getsize(DAG_block_file_1):
            pass
        else:
            os.system('echo "DAGchainer for block %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
            os.system("perl ./software/dagchainer/bin/run_DAG_chainer.pl -i %s -g 1 -D 20 -A 6 -Z 50 -e -3f" % (out_bigger_ED_1))  ## Specify your PATH to DAGCHAINER here

    sz = os.path.getsize(DAG_block_file_1)
    if not sz:
        print(str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_1.aligncoords --- It's empty!")
    else:
        os.system('echo "Start geting pairED%s_block%s_1.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))
        add_each_block_without_overlap(ref, que, each_chr, ED, block_number)
        os.system('echo "Finish geting pairED%s_block%s_1.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))

    ########################## the last block ########################
    if os.path.exists(out_sg_2) and os.path.getsize(out_sg_2):
        pass
    else:
        if os.path.exists(DAG_block_file_2) and os.path.getsize(DAG_block_file_2):
            pass
        else:
            os.system('echo "DAGchainer for block %s" >> %s_%s_%s_test.log' % (block_number, each_chr, ref, que))
            os.system("perl ./software/dagchainer/bin/run_DAG_chainer.pl -i %s -g 1 -D 20 -A 6 -Z 50 -e -3f" % (out_bigger_ED_2))

    sz = os.path.getsize(DAG_block_file_2)
    if not sz:
        print(str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(block_number) + "_2.aligncoords --- It's empty!")
    else:
        os.system('echo "Start geting pairED%s_block%s_2.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))
        add_each_block_without_overlap(ref, que, each_chr, ED, block_number)
        os.system('echo "Finish geting pairED%s_block%s_2.sg" >> %s_%s_%s_test.log' % (ED, block_number, each_chr, ref, que))


def count_lines(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file)


def add_each_block_without_overlap(ref, que, each_chr, ED, block_number):
    DAG_block_file = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(
        block_number) + ".aligncoords"
    out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED" + str(ED) + "_block" + str(
        block_number) + ".sg"
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
                            df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2",
                                                     "evalue", "score"]
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
                            df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2",
                                                     "evalue", "score"]
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
                df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue",
                                         "score"]
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
    parser.add_argument('-b', '--blcok', help='block', dest='block', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    split_block_add(args.ref, args.que, args.chr, args.ED, args.block, args.Dir)


if __name__ == '__main__':
    main()