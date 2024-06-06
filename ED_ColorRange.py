import os
import pandas as pd
import argparse

def filter_in_each_window(ref, que, each_chr, ED, Dir):
    os.system('echo "Start filter_in_each_window" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    infile = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_for_plot_ED" + str(ED) + ".txt"
    df_all = pd.read_table(infile, sep="\t", header=None, names=["chr", "id", "start", "end", "y", "y_pair", "material", "pair", "family"])
    df_pair = df_all[df_all.duplicated(subset=["pair"], keep=False)]
    df_select_ref = df_pair[(df_pair["material"] == str(ref))].reset_index(drop=True)
    df_pair_sorted = df_select_ref.sort_values(by='start')
    df_pair_list = df_pair_sorted.drop_duplicates(subset=['pair'], keep='first', ignore_index=True)
    pair_list = df_pair_list['pair'].tolist()
    out_temp = str(each_chr) + "_" + str(ref) + "_" + str(que) + "pair_info_for_window_ED" + str(ED) + ".txt"
    outfile = str(each_chr) + "_" + str(ref) + "_" + str(que) + "_for_plot_ED" + str(ED) + "_filtered.txt"
    with open(out_temp, 'w') as ouf:
        for index, pair in enumerate(pair_list):
            if index == 0:
                df_select_ref = df_pair[(df_pair["pair"] == str(pair)) & (df_pair["material"] == str(ref))].reset_index(drop=True)
                id_r = df_select_ref.loc[0, "id"]
                index_r = str(id_r).split(".")[-1]
                df_select_que = df_pair[(df_pair["pair"] == str(pair)) & (df_pair["material"] == str(que))].reset_index(drop=True)
                id_q = df_select_que.loc[0, "id"]
                index_q = str(id_q).split(".")[-1]
                index_last = index_q
                index_order = int(index_q) - int(index_last)
                ouf.write(str(pair) + "\t" + str(index_r) + "\t" + str(index_q) + "\t" + str(index_order) + "\n")
            else:
                df_select_ref = df_pair[(df_pair["pair"] == str(pair)) & (df_pair["material"] == str(ref))].reset_index(drop=True)
                id_r = df_select_ref.loc[0, "id"]
                index_r = str(id_r).split(".")[-1]
                df_select_que = df_pair[(df_pair["pair"] == str(pair)) & (df_pair["material"] == str(que))].reset_index(drop=True)
                id_q = df_select_que.loc[0, "id"]
                index_q = str(id_q).split(".")[-1]
                df_select_que_last = df_pair[(df_pair["pair"] == str(pair_list[index - 1])) & (df_pair["material"] == str(que))].reset_index(drop=True)
                id_q_last = df_select_que_last.loc[0, "id"]
                index_last = str(id_q_last).split(".")[-1]
                index_order = int(index_q) - int(index_last)
                ouf.write(str(pair) + "\t" + str(index_r) + "\t" + str(index_q) + "\t" + str(index_order) + "\n")

    ### Remove rows that are not in order
    df_selected = pd.read_table(out_temp, sep="\t", header=None, names=["Pair_name", "index_ref", "index_que", "index_check"])
    pairs_should_change = []
    for index_select, row_select in df_selected.iterrows():
        index_check = row_select["index_check"]
        index_que = row_select["index_que"]
        n = 1
        if index_check < 0:
            while index_que - df_selected.loc[index_select - n - 1, "index_que"] < 0:
                n += 1
                if index_select - n - 1 == 0:
                    n += 1
                    break
        if index_check < 0:
            for rownumber in range(0, n):
                Pair_name = df_selected.loc[index_select - rownumber - 1, "Pair_name"]
                if Pair_name not in pairs_should_change:
                    pairs_should_change.append(Pair_name)

    print("pairs_should_change: " + str(pairs_should_change))
    os.system('echo "pairs_should_change: %s" >> %s_%s_%s_test.log' % (pairs_should_change, each_chr, ref, que))
    for pairs in pairs_should_change:
        df_all.loc[(df_all["pair"] == str(pairs)) & (df_all["material"] == str(que)), "pair"] = pairs + "_1"

    os.system('echo "Finish filter_in_each_window" >> %s_%s_%s_test.log' % (each_chr, ref, que))

    ######### Extract the new drawing format with the filtered dataframe above
    os.system('echo "Geting new format for plot" >> %s_%s_%s_test.log' % (each_chr, ref, que))
    #### ED
    infile = str(Dir) + "03_allOrderDAG/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_group_pairED_order_all.dag"
    df_uniq_matchlist = pd.read_table(infile, sep="\t", header=None, names=["chr_r", "id_r", "s_r", "e_r", "chr_q", "id_q", "s_q", "e_q", "prop", "ED"])
    #### pairs
    df_only_pair = df_all[df_all.duplicated(subset=["pair"], keep=False)]  ### 只取pair对的进行过滤
    df_pair_list_new = df_only_pair.drop_duplicates(subset=['pair'], keep='first', ignore_index=True)
    pair_list_new = df_pair_list_new['pair'].tolist()
    ##### out for plot
    out_ED_colorrange = str(each_chr) + "_" + str(ref) + "_" + str(que) + "pair_colorrange_ED" + str(ED) + ".txt"
    with open(out_ED_colorrange, 'w') as outf:
        for index, pair in enumerate(pair_list_new):
            df_selected_r = df_only_pair[(df_only_pair["pair"] == str(pair)) & (df_only_pair["material"] == str(ref))].reset_index(drop=True)
            id_r = df_selected_r.loc[0, "id"]  ## Chr12_NIP.LTR.71
            if id_r.split(".")[-2] != "LTR":
                s1 = df_selected_r.loc[0, "start"]
                e1 = df_selected_r.loc[0, "end"]
                y1 = float(df_selected_r.loc[0, "y"]) - 0.3
                df_selected_q = df_only_pair[(df_only_pair["pair"] == str(pair)) & (df_only_pair["material"] == str(que))].reset_index(drop=True)
                s2 = df_selected_q.loc[0, "start"]
                e2 = df_selected_q.loc[0, "end"]
                id_q = df_selected_q.loc[0, "id"]
                y2 = float(df_selected_q.loc[0, "y"]) + 0.3
                df_ED = df_uniq_matchlist[(df_uniq_matchlist["id_r"] == str(id_r)) & (df_uniq_matchlist["id_q"] == str(id_q))].reset_index(drop=True)
                if df_ED.shape[0] != 0:
                    ED_pair = df_ED.loc[0, "ED"]
                    outf.write(str(each_chr) + "\t" + str(s1) + "\t" + str(y1) + "\t" + str(index) + "\t1\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(e1) + "\t" + str(y1) + "\t" + str(index) + "\t2\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(e2) + "\t" + str(y2) + "\t" + str(index) + "\t3\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(s2) + "\t" + str(y2) + "\t" + str(index) + "\t4\t" + str(ED_pair) + "\n")
            else:
                s1 = df_selected_r.loc[0, "start"]
                e1 = df_selected_r.loc[0, "end"]
                y1 = float(df_selected_r.loc[0, "y"]) - 0.3
                df_selected_q = df_only_pair[(df_only_pair["pair"] == str(pair)) & (df_only_pair["material"] == str(que))].reset_index(drop=True)
                s2 = df_selected_q.loc[0, "start"]
                e2 = df_selected_q.loc[0, "end"]
                id_q = df_selected_q.loc[0, "id"]
                y2 = float(df_selected_q.loc[0, "y"]) + 0.3
                df_ED = df_uniq_matchlist[(df_uniq_matchlist["id_r"] == str(id_r)) & (df_uniq_matchlist["id_q"] == str(id_q))].reset_index(drop=True)

                if df_ED.shape[0] != 0:
                    ED_pair = "LTR"
                    outf.write(str(each_chr) + "\t" + str(s1) + "\t" + str(y1) + "\t" + str(index) + "\t1\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(e1) + "\t" + str(y1) + "\t" + str(index) + "\t2\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(e2) + "\t" + str(y2) + "\t" + str(index) + "\t3\t" + str(ED_pair) + "\n")
                    outf.write(str(each_chr) + "\t" + str(s2) + "\t" + str(y2) + "\t" + str(index) + "\t4\t" + str(ED_pair) + "\n")



    df_all.to_csv(outfile, sep="\t", header=False, index=False)
    os.system("rm %s" % (out_temp))
    os.system("cat %s >> %s_%s_12chr_for_plot_ED%s.txt" % (outfile, ref, que, ED))
    os.system("rm %s" % (outfile))
    os.system("rm %s_%s_%s_for_plot_ED%s.txt" % (each_chr, ref, que, ED))
    os.system("cat %s >> %s_%s_12chr_ColorRange_ED%s.txt" % (out_ED_colorrange, ref, que, ED))
    os.system("rm %s" % (out_ED_colorrange))
    os.system('echo "Finish geting new format for plot" >> %s_%s_%s_test.log' % (each_chr, ref, que))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-c', '--chr', help='chr', dest='chr', required=True)
    parser.add_argument('-e', '--ed', help='ED', dest='ED', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    filter_in_each_window(args.ref, args.que, args.chr, args.ED, args.Dir)

if __name__ == '__main__':
    main()