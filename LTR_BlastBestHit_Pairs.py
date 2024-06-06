import os
import pandas as pd
import argparse

def get_LTR_temp_bed(ref, que, each_chr, Dir):
    #### Prepare the annotation files in advance
    LTR_r = str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt"
    df_LTR_r = pd.read_table(LTR_r, sep="\t", header=None, names=["Chr", "Start", "End", "Strand", "ID", "Family"])
    LTR_q = str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt"
    df_LTR_q = pd.read_table(LTR_q, sep="\t", header=None, names=["Chr", "Start", "End", "Strand", "ID", "Family"])

    ############################### Get TE sequences from the fasta file  ################################
    ######   reference centromeric region   #########
    if os.path.exists(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_temp_LTR.bed"):
        pass
    else:
        df_chr_LTR_r = df_LTR_r[(df_LTR_r["Chr"] == str(each_chr) + "_" + str(ref)) & ((df_LTR_r["Family"] == "LTR") | (df_LTR_r["Family"] == "other_repeat"))].reset_index(drop=True)
        if df_chr_LTR_r.shape[0] >= 1:
            out_r = str(Dir)  + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_temp_LTR.bed"
            with open(out_r, 'w') as o_r:
                for index_r, row_r in df_chr_LTR_r.iterrows():
                    Start_r = row_r["Start"]
                    End_r = row_r["End"]
                    o_r.write(str(each_chr) + "_" + str(ref) + "\t" + str(Start_r) + "\t" + str(End_r) + "\n")

    ######   query centromeric region   #########
    if os.path.exists(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(que) + "_temp_LTR.bed"):
        pass
    else:
        df_chr_LTR_q = df_LTR_q[(df_LTR_q["Chr"] == str(each_chr) + "_" + str(que)) & ((df_LTR_q["Family"] == "LTR") | (df_LTR_q["Family"] == "other_repeat"))].reset_index(drop=True)
        if df_chr_LTR_q.shape[0] >= 1:
            out_q = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(que) + "_temp_LTR.bed"
            with open(out_q, 'w') as o_q:
                for index_q, row_q in df_chr_LTR_q.iterrows():
                    Start_q = row_q["Start"]
                    End_q = row_q["End"]
                    o_q.write(str(each_chr) + "_" + str(que) + "\t" + str(Start_q) + "\t" + str(End_q) + "\n")


def get_LTR_seq_fasta_blast(ref, que, each_chr, Dir):
    if (os.path.exists(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(que) + "_temp_LTR.bed")) and \
            (os.path.exists(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_temp_LTR.bed")):
        if os.path.exists(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_LRT_blast_best_hit.txt"):
            pass
        else:
            fastafile = str(ref) + ".fasta"
            outfasta_r = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_LTR_for_blast.fasta"
            os.system("bedtools getfasta -fi %s -bed ../02_LTRBlast/%s_%s_temp_LTR.bed -fo %s" % (fastafile, each_chr, ref, outfasta_r))  ## 利用bed文件从fasta文件中提取序列的方法
            fastafile_q = str(que) + ".fasta"
            outfasta_q = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(que) + "_LTR_for_blast.fasta"
            os.system("bedtools getfasta -fi %s -bed ../02_LTRBlast/%s_%s_temp_LTR.bed -fo %s" % (fastafile_q, each_chr, que, outfasta_q))  ## 利用bed文件从fasta文件中提取序列的方法

            ### change LTR id format
            ## Reference
            df_temp_LTR_r = pd.read_table(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_temp_LTR.bed", sep='\t', header=None, names=["Chr", "Start", "End"])
            df_new_id_r = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(ref) + "_full_annotation.rmdup.txt", sep='\t', header=None, names=["Chr", "Start", "End", "Strand", "New_id", "Group"])
            df_new_id_r_chr = df_new_id_r[df_new_id_r["Chr"] == str(each_chr) + "_" + str(ref)].reset_index(drop=True)
            LTR_id_change_r = {}
            for index_temp_LTR_r, row_temp_LTR_r in df_temp_LTR_r.iterrows():
                Start = row_temp_LTR_r["Start"]
                End = row_temp_LTR_r["End"]
                df_select_new_id_r = df_new_id_r_chr[(df_new_id_r_chr["Start"] == Start) & (df_new_id_r_chr["End"] == End)]
                old_id = str(row_temp_LTR_r["Chr"]) + ":" + str(row_temp_LTR_r["Start"]) + "-" + str(row_temp_LTR_r["End"])
                new_id_r = None
                for iidex_r, rrow_r in df_select_new_id_r.iterrows():
                    new_id_r = str(rrow_r["Chr"]) + ".LTR." + str(iidex_r)
                LTR_id_change_r[old_id] = new_id_r

            ## Query
            df_temp_LTR_q = pd.read_table(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(que) + "_temp_LTR.bed", sep='\t', header=None, names=["Chr", "Start", "End"])
            df_new_id_q = pd.read_table(str(Dir) + "00_AllAnnotation/" + str(que) + "_full_annotation.rmdup.txt", sep='\t', header=None, names=["Chr", "Start", "End", "Strand", "New_id", "Group"])
            df_new_id_q_chr = df_new_id_q[df_new_id_q["Chr"] == str(each_chr) + "_" + str(que)].reset_index(drop=True)
            LTR_id_change_q = {}
            for index_temp_LTR_q, row_temp_LTR_q in df_temp_LTR_q.iterrows():
                Start = row_temp_LTR_q["Start"]
                End = row_temp_LTR_q["End"]
                df_select_new_id_q = df_new_id_q_chr[(df_new_id_q_chr["Start"] == Start) & (df_new_id_q_chr["End"] == End)]
                old_id_q = str(row_temp_LTR_q["Chr"]) + ":" + str(row_temp_LTR_q["Start"]) + "-" + str(row_temp_LTR_q["End"])
                new_id_q = None
                for iidex_q, rrow_q in df_select_new_id_q.iterrows():
                    new_id_q = str(rrow_q["Chr"]) + ".LTR." + str(iidex_q)
                LTR_id_change_q[old_id_q] = new_id_q

            ###   BLAST results to matchlist format(used for DAGchainer)
            blastoutdb = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_LTR"
            blastout = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + ".blast"
            os.system("makeblastdb -in %s -dbtype nucl -out %s" % (outfasta_r, blastoutdb))
            os.system("blastn -query %s -db %s -evalue 1e-6 -outfmt 6 -num_threads 6 -out %s" % (outfasta_q, blastoutdb, blastout))
            sz = os.path.getsize(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + ".blast")
            if not sz:
                print(str(each_chr) + "_" + str(ref) + "_" + str(que) + ".blast --- It's empty!")
            else:
                df_blast = pd.read_table(str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + ".blast", sep="\t", header=None)
                df_blast.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gap", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

                outfile = str(Dir) + "02_LTRBlast/" + str(each_chr) + "_" + str(ref) + "_" + str(que) + "_LRT_blast_best_hit.txt"
                o = open(outfile, 'w')
                for index, row in df_blast.iterrows():
                    chr_r = str(row["sseqid"]).split(":")[0]
                    id_r = LTR_id_change_r[str(row["sseqid"])]
                    new_s_r = int(str(row["sseqid"]).split(":")[1].split("-")[0])
                    new_e_r = int(str(row["sseqid"]).split(":")[1].split("-")[1])
                    chr_q = str(row["qseqid"]).split(":")[0]
                    id_q = LTR_id_change_q[str(row["qseqid"])]
                    new_s_q = int(str(row["qseqid"]).split(":")[1].split("-")[0])
                    new_e_q = int(str(row["qseqid"]).split(":")[1].split("-")[1])
                    evalue = row["evalue"]
                    o.write(str(chr_r) + "\t" + str(id_r) + "\t" + str(new_s_r) + "\t" + str(new_e_r) + "\t" + str(chr_q) + "\t" + str(id_q) + "\t" + str(new_s_q) + "\t" + str(new_e_q) + "\t" + str(evalue) + "\t" + str(0) + "\n")
            os.system("rm ../02_LTRBlast/%s_%s_%s.blast" % (each_chr, ref, que))
            os.system("rm ../02_LTRBlast/%s_%s_LTR.n*" % (each_chr, ref))
            os.system("rm ../02_LTRBlast/%s_%s_LTR_for_blast.fasta" % (each_chr, ref))
            os.system("rm ../02_LTRBlast/%s_%s_LTR_for_blast.fasta" % (each_chr, que))
    else:
        print("No LTR !")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='ref', dest='ref', required=True)
    parser.add_argument('-q', '--que', help='que', dest='que', required=True)
    parser.add_argument('-chr', '--chr', help='Chr', dest='chr', required=True)
    parser.add_argument('-d', '--Dir', help='Dir', dest='Dir', required=True)
    args = parser.parse_args()
    get_LTR_temp_bed(args.ref, args.que, args.chr, args.Dir)
    get_LTR_seq_fasta_blast(args.ref, args.que, args.chr, args.Dir)

if __name__ == '__main__':
    main()