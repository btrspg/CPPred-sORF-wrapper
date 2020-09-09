#! /usr/bin/env python3
import argparse as agp
import os
import tempfile
from subprocess import check_call

# from feamodule import ORF_length as len
import Bio.SeqIO as Seq

from cppredsorf import CTD, fickett, FrameKmer, ProtParam as PP, ORF_length_11codon_GCcount as len
from cppredsorf.utils import get_model_range_hexamer


def print_and_run(cmd):
    '''

    '''
    print("CMD:" + cmd)
    check_call(cmd, shell=True)


def get_tempdir():
    return tempfile.TemporaryDirectory()


def path_file(path, f):
    '''

    '''
    return os.path.join(path, f)


def coding_nocoding_potential(input_file):
    coding = {}
    noncoding = {}
    for line in open(input_file).readlines():
        fields = line.split()
        if fields[0] == 'hexamer': continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])
    return coding, noncoding


def output_feature(seq_file, hex_file, species, tmpdir):
    tmp = open(path_file(tmpdir, 'test.f_svm'), 'w')
    feature = open(path_file(tmpdir, 'test.feature'), 'w')
    out_label = 1
    coding, noncoding = coding_nocoding_potential(hex_file)
    if species == "Human":
        feature.write("\t".join(map(str, ["#ID", "ORF-integrity", "ORF-coverage", "Instability", "T2", "C0", "PI",
                                          "ORF-length", "AC", "T0", "G0", "C2", "A4", "G2", "TG", "A0", "TC", "G1",
                                          "C3", "T3", "A1", "GC", "T1", "G4", "C1", "G3", "A3", "Gravy", "Hexamer",
                                          "C4", "AG", "Fickett", "A2", "T4", "C", "G", "A", "T", "mRNN_11codon",
                                          "GCcount"])) + "\n")
    if species == "Integrated":
        #		feature.write("\t".join(map(str,["#ID","ORF-coverage","ORF-integrity","GC","Instability","ORF-length","T0","Fickett","G2","C3","PI","A3","C1","G3","Hexamer","TG","G1","TC","A0","A1","AC","C2","G0","T4","C0","A4","G","A2","T","T3","G4","C4","Grary","T2","AG","AT","T1","A","C","mRNN_11codon","GCcount"]))+"\n")
        feature.write("\t".join(map(str,
                                    ["#ID", "ORF-coverage", "PI", "Hexamer", "ORF_length", "Fickett", "C0", "C4", "T1",
                                     "TG", "GCcount", "T4", "TC", "GC", "A4", "G4", "mRNN_11codon", "AC", "T3", "C3",
                                     "A0", "Instability", "A1", "ORF-integrity", "G1", "AG", "T0", "C1", "G3", "G0",
                                     "Gravy", "AT", "A3", "T2", "C", "C2", "G2", "A", "A2", "G"])) + "\n")
    for seq in Seq.parse(seq_file, 'fasta'):
        seqid = seq.id
        A, T, G, C, AT, AG, AC, TG, TC, GC, A0, A1, A2, A3, A4, T0, T1, T2, T3, T4, G0, G1, G2, G3, G4, C0, C1, C2, C3, C4 = CTD.CTD(
            seq.seq)
        insta_fe, PI_fe, gra_fe = PP.param(seq.seq)
        fickett_fe = fickett.fickett_value(seq.seq)
        hexamer = FrameKmer.kmer_ratio(seq.seq, 6, 3, coding, noncoding)
        Len, Cov, inte_fe, mRNN_11codon, GCcount = len.len_cov(seq.seq)
        if species == "Human":
            tem = [inte_fe, Cov, insta_fe, T2, C0, PI_fe, Len, AC, T0, G0, C2, A4, G2, TG, A0, TC, G1, C3, T3, A1, GC,
                   T1, G4, C1, G3, A3, gra_fe, hexamer, C4, AG, fickett_fe, A2, T4, C, G, A, T, mRNN_11codon, GCcount]
            feature.write("\t".join(map(str,
                                        [seqid, inte_fe, Cov, insta_fe, T2, C0, PI_fe, Len, AC, T0, G0, C2, A4, G2, TG,
                                         A0, TC, G1, C3, T3, A1, GC, T1, G4, C1, G3, A3, gra_fe, hexamer, C4, AG,
                                         fickett_fe, A2, T4, C, G, A, T, mRNN_11codon, GCcount])) + "\n")
        elif species == "Integrated":
            #			tem = [Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C,mRNN_11codon,GCcount]
            #			feature.write("\t".join(map(str,[seqid,Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C,mRNN_11codon,GCcount]))+"\n")
            tem = [Cov, PI_fe, hexamer, Len, fickett_fe, C0, C4, T1, TG, GCcount, T4, TC, GC, A4, G4, mRNN_11codon, AC,
                   T3, C3, A0, insta_fe, A1, inte_fe, G1, AG, T0, C1, G3, G0, gra_fe, AT, A3, T2, C, C2, G2, A, A2, G]
            feature.write("\t".join(map(str,
                                        [seqid, Cov, PI_fe, hexamer, Len, fickett_fe, C0, C4, T1, TG, GCcount, T4, TC,
                                         GC, A4, G4, mRNN_11codon, AC, T3, C3, A0, insta_fe, A1, inte_fe, G1, AG, T0,
                                         C1, G3, G0, gra_fe, AT, A3, T2, C, C2, G2, A, A2, G])) + "\n")
        else:
            raise ValueError('species should be human or integrated')
        tmp.write(str(out_label))
        for label, item in enumerate(tem):
            tmp.write(' ' + str(label + 1) + ':' + str(item))
        tmp.write('\n')
    tmp.close()


def predict(range_file, model_file, libsvm_bin, tmpdir):
    svm_scale = path_file(libsvm_bin, 'svm-scale ') + ' -r ' + range_file + ' ' \
                + path_file(tmpdir, 'test.f_svm ') + ' > ' + path_file(tmpdir, 'test.scaled ')
    svm_predict = path_file(libsvm_bin, 'svm-predict ') + ' -b 1 ' \
                  + path_file(tmpdir, 'test.scaled ') + model_file + ' ' \
                  + path_file(tmpdir, 'tmp.txt ') + ' > ' + path_file(tmpdir, 'tmp2.txt ')
    # os.system('../libsvm-3.22/svm-scale -r '+ range_file + ' test.f_svm  > test.scaled ')
    # os.system('../libsvm-3.22/svm-predict -b 1 test.scaled ' + model_file +' tmp.txt >  tmp2.txt')
    # os.system(
    #     libsvm_bin + '/svm-scale -r ' + range_file + ' ' + seq_fname + '_test.f_svm  >' + seq_fname + '_test.scaled ')
    # print ('../libsvm-3.22/svm-scale -r '+ range_file +' '+seq_fname+'_test.f_svm  >' +seq_fname+'_test.scaled ')
    # return 0
    # os.system('../libsvm-3.22/svm-predict -b 1 test.scaled ' + model_file +' '+seq_fname+'_tmp.txt >' +seq_fname+'_tmp2.txt')
    # os.system(
    #     libsvm_bin + '/svm-predict -b 1' + ' ' + seq_fname + '_test.scaled ' + model_file + ' ' + seq_fname + '_tmp.txt >' + seq_fname + '_tmp2.txt')
    # print ('../libsvm-3.22/svm-predict -b 1'+' '+seq_fname+'_test.scaled ' + model_file +' '+seq_fname+'_tmp.txt >' +seq_fname+'_tmp2.txt')
    # return 0
    print_and_run('ls -R '+ tmpdir)
    print_and_run(svm_scale)
    print_and_run(svm_predict)
    os.system('ls -R '+ tmpdir)
    os.system("head "+path_file(tmpdir, 'tmp.txt '))
    os.system("head " + path_file(tmpdir, 'tmp2.txt '))
    coding_poten = open(path_file(tmpdir, 'coding_potential'), 'w')
    coding_poten.write("\t".join(map(str, ["table", path_file(tmpdir, 'coding_potential')])) + "\n")

    for line in open(path_file(tmpdir, 'tmp.txt'), 'r').readlines():
        if line[0] == "l":
            continue
        coding_potential = line.split(" ")[1]
        if line.split(" ")[0] == "1":
            coding_poten.write("\t".join(map(str, ["coding", coding_potential])) + "\n")
        else:
            coding_poten.write("\t".join(map(str, ["noncoding", coding_potential])) + "\n")


def merge(tmpdir, output_file):
    cmd = "paste " + path_file(tmpdir, 'test.feature ') + path_file(tmpdir, 'coding_potential ') \
          + ' > ' + output_file
    print_and_run(cmd)


def deleted(tmpdir):
    tmpdir.cleanup()
    # os.system("rm" + " " + seq_fname + "_test.*")
    # os.system("rm" + " " + seq_fname + "_coding_potential")
    # os.system("rm" + " " + seq_fname + "_tmp*")


def main():
    parser = agp.ArgumentParser()
    parser.add_argument('-i', '--RNA_file', help="the input FASTA file of RNA sequence")
    parser.add_argument('--species', '-s', help="the species", choices=['Human', 'Integrated'])
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('--libsvm-bin', help='libsvm bin', default='/usr/local/bin')
    args = parser.parse_args()
    tmpdir = get_tempdir()
    m_f, r_f, h_f = get_model_range_hexamer()
    output_feature(args.RNA_file, h_f, args.species, tmpdir.name)
    predict(r_f, m_f, args.libsvm_bin,tmpdir.name)
    merge(tmpdir.name,args.outfile)
    deleted(tmpdir)


if __name__ == '__main__':
    main()
