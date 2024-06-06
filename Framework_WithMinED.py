import pandas as pd
import argparse
import os


def get_frame(ref, que, each_chr, Dir):
    ### get total elements number
    os.system('echo "############## 3 get_frame #####################" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    newID_ref = str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt"
    temp_r = str(each_chr) + "_" + str(ref) + "_anno.rmdup.txt"
    os.system("less %s | grep %s > %s" % (newID_ref, each_chr, temp_r))
    MAX_order_s_r = count_lines(temp_r)
    os.system("rm %s" % (temp_r))

    ### search for the existing framework
    for frame_range in range(0, 21):
        out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(frame_range) + "_framework.sg"
        if os.path.exists(out_sg) and os.path.getsize(out_sg):
            df_out_sg = pd.read_table(out_sg)
            out_sg_pairnumber = df_out_sg.shape[0]
            if (20 * out_sg_pairnumber > MAX_order_s_r) or (out_sg_pairnumber > 100):
                os.system('echo "The ED selected for frame is %s " >> %s_%s_%s_test.log' % (frame_range, each_chr, ref, que))
                os.system('echo "THEN Do ED_pairnumber_pick" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                break
            else:
                os.system("rm %s_%s_%s_pairED_frame%s_addLTR.aligncoords" % (each_chr, ref, que, frame_range))
                os.system("rm %s_%s_%s_pairED_frame%s_addLTR" % (each_chr, ref, que, frame_range))
        else:
            os.system('echo "Do use_ED_as_framework_of_partial_DAGchainer" >> %s_%s_%s_test.log' % (each_chr, ref, que))
            use_ED_as_framework_of_partial_DAGchainer(ref, que, each_chr, frame_range, Dir)
            if os.path.exists(out_sg) and os.path.getsize(out_sg):
                df_out_sg = pd.read_table(out_sg)
                out_sg_pairnumber = df_out_sg.shape[0]
                if (20 * out_sg_pairnumber > MAX_order_s_r) or (out_sg_pairnumber > 100):
                    os.system('echo "The ED selected for frame is %s " >> %s_%s_%s_test.log' % (frame_range, each_chr, ref, que))
                    os.system('echo "THEN Do ED_pairnumber_pick" >> %s_%s_%s_test.log' % (each_chr, ref, que))
                    break
                else:
                    os.system("rm %s_%s_%s_pairED_frame%s_addLTR.aligncoords" % (each_chr, ref, que, frame_range))
                    os.system("rm %s_%s_%s_pairED_frame%s_addLTR" % (each_chr, ref, que, frame_range))
            else:
                os.system("rm %s_%s_%s_pairED_frame%s_addLTR.aligncoords" % (each_chr, ref, que, frame_range))
                os.system("rm %s_%s_%s_pairED_frame%s_addLTR" % (each_chr, ref, que, frame_range))



def use_ED_as_framework_of_partial_DAGchainer(ref, que, each_chr, ED_frame_select, Dir):
    out_ED_addLTRdag = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED_frame" + str(ED_frame_select) + "_addLTR"
    if os.path.exists(out_ED_addLTRdag) and os.path.getsize(out_ED_addLTRdag):
        pass
    else:
        in_re_format = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"
        with open(in_re_format, 'r') as inf:
            with open(out_ED_addLTRdag, 'w') as outf:
                for line in inf: ## "chr_r", "id_r", "s_r", "e_r", "chr_q", "id_q", "s_q", "e_q", "prop", "ED"
                    text = line.strip().split("\t")
                    chr_r = text[0]
                    id_r = text[1]
                    index_r = text[2]
                    chr_q = text[4]
                    id_q = text[5]
                    index_q = text[6]
                    prop = text[8]
                    ED = text[9]
                    if int(ED) <= int(ED_frame_select):
                        outf.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(index_r) + "\t" + str(index_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(index_q) + "\t" + str(index_q) + "\t" + str(prop) + "\n")

    ###### get 'sg' format as the result for framework construction
    out_sg_temp = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework_temp.sg"
    DAG_file = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_pairED_frame" + str(ED_frame_select) + "_addLTR.aligncoords"
    if os.path.exists(DAG_file) and os.path.getsize(DAG_file):
        pass
    else:
        os.system('echo "To get %s_%s_%s_pairED_frame%s_addLTR.aligncoords" >> %s_%s_%s_test.log' % (each_chr, ref, que, ED_frame_select, each_chr, ref, que))
        os.system("perl ./software/dagchainer/bin/run_DAG_chainer.pl -i %s -g 1 -D 20 -A 6 -Z 50 -e -3f" % (out_ED_addLTRdag))  ## Specify your PATH to DAGCHAINER here

    ###### 对DAG结果文件进行处理：按照每个block的顺序，逐个添加
    if not os.path.getsize(DAG_file):
        os.system('echo "To get %s_%s_%s_pairED_frame%s_addLTR.aligncoords --- It is empty!" >> %s_%s_%s_test.log' % (each_chr, ref, que, ED_frame_select, each_chr, ref, que))
    else:
        os.system('echo "Start geting temp ED_frame%s_select_framework.sg" >> %s_%s_%s_test.log' % (ED_frame_select, each_chr, ref, que))
        ###### !! Important: Add non-overlapping results in block order
        block_data = []
        block_num = 0
        position_dict_ref = {}
        position_dict_que = {}
        with open(DAG_file, 'r') as inf:
            with open(out_sg_temp, 'w') as outf:
                for lines in inf:
                    line = lines.strip()
                    if line.startswith('#'):
                        if block_data:
                            block_num += 1
                            if block_num == 1:
                                df_block_data = pd.DataFrame([x.strip().split('\t') for x in block_data])
                                df_block_data.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue", "score"]
                                df_block_data['s_1'] = df_block_data['s_1'].astype(int)
                                df_block_data['s_2'] = df_block_data['s_2'].astype(int)
                                min_s_1 = df_block_data['s_1'].min()
                                max_s_1 = df_block_data['s_1'].max()
                                min_s_2 = df_block_data['s_2'].min()
                                max_s_2 = df_block_data['s_2'].max()
                                position_dict_ref[block_num] = [min_s_1, max_s_1]
                                position_dict_que[block_num] = [min_s_2, max_s_2]
                                #### Keep everything in the first block
                                for lines_block in block_data:
                                    line_block = lines_block.strip()
                                    outf.write(line_block + "\n")
                                block_data = []
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
                                    #### 如果这个pair对的位置不在任何一个已有的bolck中，则添加它
                                    if n_not_in_count_1 == len(position_dict_ref) and n_not_in_count_2 == len(position_dict_que):
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
                        ###将非'#'开头的行添加到当前block的数据中
                        block_data.append(line)

        os.system('echo "Finish geting temp ED_frame%s_select_framework.sg" >> %s_%s_%s_test.log' % (ED_frame_select, each_chr, ref, que))

        os.system('echo "Start geting final ED_frame%s_select_framework.sg" >> %s_%s_%s_test.log' % (ED_frame_select, each_chr, ref, que))
        if os.path.exists(out_sg_temp) and os.path.getsize(out_sg_temp):
            get_mostnum_order_pairs(each_chr, ref, que, ED_frame_select)
            out_test_temp = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework_test.sg"
            os.system("rm %s %s" % (out_test_temp, out_sg_temp))

def get_mostnum_order_pairs(each_chr, ref, que, ED_frame_select):
    out_sg_temp = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework_temp.sg"
    out_sg = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework.sg"
    out_test_temp = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED" + str(ED_frame_select) + "_framework_test.sg"

    df_sg = pd.read_table(out_sg_temp, sep="\t", header=None)
    df_sg.columns = ["chr_1", "id_1", "s_1", "e_1", "chr_2", "id_2", "s_2", "e_2", "evalue", "score"]
    block_num = 1
    order_num = 0
    if df_sg.loc[0, "chr_1"].split("_")[1] == "NIP":
        with open(out_test_temp, 'w') as outf:
            df_sg_sort = df_sg.sort_values(by="s_1").reset_index(drop=True)
            for index, row in df_sg_sort.iterrows():
                start_que = row["s_2"]
                if index == 0:
                    pass
                elif index != df_sg_sort.shape[0] - 1:
                    start_que_last = df_sg_sort.loc[index - 1, "s_2"]
                    if start_que > start_que_last and start_que < start_que_last + 20:
                        order_num += 1
                    elif index != df_sg_sort.shape[0] - 1:
                        index_s = index - order_num - 1
                        block_start = df_sg_sort.loc[index_s, "s_2"]
                        index_e = index - 1
                        block_end = df_sg_sort.loc[index_e, "s_2"]
                        outf.write(str(each_chr) + "\t" + str(ref) + "\t" + str(que) + "\t" + str(block_num) + "\t" + str(block_start) + "\t" + str(block_end) + "\t" + str(index_s) + "\t" + str(index_e) + "\t" + str(order_num + 1) + "\n")
                        order_num = 0
                        block_num += 1
                else:
                    index_s = index - order_num - 1
                    index_e = index
                    block_start = df_sg_sort.loc[index_s, "s_2"]
                    block_end = df_sg_sort.loc[index_e, "s_2"]
                    outf.write(str(each_chr) + "\t" + str(ref) + "\t" + str(que) + "\t" + str(block_num) + "\t" + str(block_start) + "\t" + str(block_end) + "\t" + str(index_s) + "\t" + str(index_e) + "\t" + str(order_num + 1) + "\n")

        df_block = pd.read_table(out_test_temp, sep="\t", header=None)
        df_block.columns = ["chr", "ref", "que", "block_num", "block_start", "block_end", "index_start", "index_end","order_num"]
        df_block["order_num"] = df_block["order_num"].astype(int)
        totalline = df_block.shape[0]
        index_dict_count = {}
        for indexget in range(0, totalline):
            df_pandaun = df_block.iloc[indexget:totalline + 1]
            first_start = df_pandaun.loc[indexget, "block_start"]
            first_order_num = df_pandaun.loc[indexget, "order_num"]
            panduan_count = first_order_num
            bigger = [first_start]
            for indexpanduan, rowpanduan in df_pandaun.iterrows():
                block_start_panduan = rowpanduan["block_start"]
                order_num_panduan = rowpanduan["order_num"]
                if block_start_panduan > bigger[-1]:
                    bigger.append(block_start_panduan)
                    panduan_count += order_num_panduan
                else:
                    pass
            index_dict_count[indexget] = panduan_count

        max_indexget = max(index_dict_count, key=index_dict_count.get)
        df_get_sg = df_block.iloc[max_indexget:totalline + 1]
        first_startsg = df_block.loc[max_indexget, "block_start"]
        if_bigger = [first_startsg]
        for indexgetsg, rowgetsg in df_get_sg.iterrows():
            block_start_sg = rowgetsg["block_start"]
            index_start_sg = rowgetsg["index_start"]
            index_end_sg = rowgetsg["index_end"]
            if block_start_sg >= if_bigger[-1]:
                if_bigger.append(block_start_sg)
                subset_df = df_sg_sort.iloc[index_start_sg:index_end_sg + 1]
                subset_df.to_csv(out_sg, sep="\t", header=False, index=False, mode='a')
    else:
        with open(out_test_temp, 'w') as outf:
            df_sg_sort = df_sg.sort_values(by="s_2").reset_index(drop=True)
            for index, row in df_sg_sort.iterrows():
                start_que = row["s_1"]
                if index == 0:
                    pass
                elif index != df_sg_sort.shape[0] - 1:
                    start_que_last = df_sg_sort.loc[index - 1, "s_1"]
                    if start_que > start_que_last and start_que < start_que_last + 20:
                        order_num += 1
                    else:
                        index_s = index - order_num - 1
                        block_start = df_sg_sort.loc[index_s, "s_1"]
                        index_e = index - 1
                        block_end = df_sg_sort.loc[index_e, "s_1"]
                        outf.write(str(each_chr) + "\t" + str(ref) + "\t" + str(que) + "\t" + str(block_num) + "\t" + str(block_start) + "\t" + str(block_end) + "\t" + str(index_s) + "\t" + str(index_e) + "\t" + str(order_num + 1) + "\n")

                        order_num = 0
                        block_num += 1
                else:
                    index_s = index - order_num - 1
                    block_start = df_sg_sort.loc[index_s, "s_1"]
                    index_e = index
                    block_end = df_sg_sort.loc[index_e, "s_1"]
                    outf.write(str(each_chr) + "\t" + str(ref) + "\t" + str(que) + "\t" + str(block_num) + "\t" + str(block_start) + "\t" + str(block_end) + "\t" + str(index_s) + "\t" + str(index_e) + "\t" + str(order_num + 1) + "\n")

        df_block = pd.read_table(out_test_temp, sep="\t", header=None)
        df_block.columns = ["chr", "ref", "que", "block_num", "block_start", "block_end", "index_start", "index_end","order_num"]
        df_block["order_num"] = df_block["order_num"].astype(int)
        totalline = df_block.shape[0]
        index_dict_count = {}
        for indexget in range(0, totalline):
            df_pandaun = df_block.iloc[indexget:totalline + 1]
            first_start = df_pandaun.loc[indexget, "block_start"]
            first_order_num = df_pandaun.loc[indexget, "order_num"]
            panduan_count = first_order_num
            bigger = [first_start]
            for indexpanduan, rowpanduan in df_pandaun.iterrows():
                block_start_panduan = rowpanduan["block_start"]
                order_num_panduan = rowpanduan["order_num"]
                if block_start_panduan > bigger[-1]:
                    bigger.append(block_start_panduan)
                    panduan_count += order_num_panduan
                else:
                    pass
            index_dict_count[indexget] = panduan_count

        max_indexget = max(index_dict_count, key=index_dict_count.get)
        df_get_sg = df_block.iloc[max_indexget:totalline + 1]
        first_startsg = df_block.loc[max_indexget, "block_start"]
        if_bigger = [first_startsg]
        for indexgetsg, rowgetsg in df_get_sg.iterrows():
            block_start_sg = rowgetsg["block_start"]
            index_start_sg = rowgetsg["index_start"]
            index_end_sg = rowgetsg["index_end"]
            if block_start_sg >= if_bigger[-1]:
                if_bigger.append(block_start_sg)
                subset_df = df_sg_sort.iloc[index_start_sg:index_end_sg + 1]
                subset_df.to_csv(out_sg, sep="\t", header=False, index=False, mode='a')

def count_lines(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-chr', '--chr', help='Chr', dest='chr', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    get_frame(args.ref, args.que, args.chr, args.Dir)


if __name__ == '__main__':
    main()
