from Bio import SeqIO
from Bio.Seq import Seq

# This script takes a nextalign aligned set of sequences, and creates
# an output file that includes only gene coding regions.

projectdir = "your_directory"
#outDir is the directory with the output of nextalign
outDir = projectdir + '/output'
codeDir = projectdir + '/CodingRegions'

#fastaFiles 
fastaFiles = ["alignment.fasta"]

#Original genemap had overlapping final gene
#genemap = [[266,13468],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27756,27887],[27894,28259],[28274,29533],[28284,28577]]
genemap = [[266,13467],[13468,21555],[21563,25384],[25393,26220],[26245,26472],[26523,27191],[27202,27387],[27394,27759],[27756,27887],[27894,28259],[28274,29533],[29558,29674]]

masks = [150,153,635,1707,1895,2091,2094,2198,2604,3145,3564,3639,3778,4050,4221,5011,5257,5736,5743,5744,6167,6255,6869,6874,8022,8026,8328,8790,8827,8828,8835,8886,8887,9039,10129,10239,10554,10716,11074,11083,11535,13117,13402,13408,13476,13571,13599,13687,14222,14223,14225,14277,14851,14852,15435,15521,15922,16290,16887,17178,17179,17182,17567,19286,19298,19299,19484,19548,20056,20123,20465,21149,21151,21209,21212,21550,21551,21575,21968,21987,22335,22516,22521,22651,22661,22802,24389,24390,24410,24622,24933,25202,25381,25382,26549,27658,27660,27760,27761,27784,28184,28253,28985,29037,29039,29058,29425,29553,29594,29783,29827,29830]

for f in fastaFiles:
    fastaFile = outDir + '/' + f
    outFile = codeDir + '/' + f
    coding = []

    with open(fastaFile) as fasta:
        for sample in SeqIO.parse(fasta, "fasta"):
            masked = Seq("")
            prevpos = 0
            for m in masks:
                masked += sample.seq[prevpos:m-1]
                masked += '-'
                prevpos = m
            masked += sample.seq[prevpos:len(sample.seq) - 1]

            chopped = Seq("")
            gnum = 0
            for g in genemap:
                chopped += masked[g[0]-1:g[1]]
                if 0 == gnum:
                    chopped += '-' # Addressing "Ribosomal Slippage" between Orf1a and Orf1b
                gnum += 1
                print("Gene [" + str(g[0]) + ":" + str(g[1]) + "]: " + masked[g[0]-1:(g[0]-1)+10] + "..." + masked[g[1]-10:g[1]])
                
            sample.seq = chopped
            coding.append(sample)

    SeqIO.write(coding, outFile, "fasta")

