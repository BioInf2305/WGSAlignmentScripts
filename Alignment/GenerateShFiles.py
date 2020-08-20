import sys
import subprocess
import argparse

class ProduceAlignmentScript:

    def __init__(self,sampleFile,sample):
        self.sampleFile=sampleFile
        self.sample=sample

    def MakeSlurmFiles(self):
        cnt=-2
        dest=open("RunAlignment"+self.sample+".sh",'w')
        dest.write("#!/bin/bash")
        dest.write("\n")
        dest.write("#SBATCH --time=40:00:00")
        dest.write("\n")
        dest.write("#SBATCH --ntasks=16")
        dest.write("\n")
        dest.write("#SBATCH --mem=16000")
        dest.write("\n")
        dest.write("#SBATCH --nodes=1")
        dest.write("\n")
        dest.write("#SBATCH --output=Output_"+self.sample+"_align_%j.txt")
        dest.write("\n")
        dest.write("#SBATCH --error=Error_"+self.sample+"_align_%j.txt")
        dest.write("\n")
        dest.write("#SBATCH --job-name="+self.sample+"_align")
        dest.write("\n")
        dest.write("\n")
        dest.write("export PATH=/home/maulik/software/bwa:$PATH")
        dest.write("\n")
        dest.write("export PATH=/home/maulik/software/FastQC:$PATH")
        dest.write("\n")
        dest.write("export PATH=/home/maulik/software/samtools:$PATH")
        dest.write("\n")
        dest.write("export PATH=/home/maulik/software/sickle:$PATH")
        dest.write("\n")
        with open(self.sampleFile) as source:
            for line in source:
                cnt+=1
                a=line.rstrip().split()
                if cnt==-1:
                    path=a[0]
                    RunDir="Dir_"+self.sample
                    dest.write("mkdir "+RunDir)
                    dest.write("\n")
                    dest.write("echo \"running sickle on sample "+self.sample+"\"")
                    dest.write("\n")
                elif cnt==0:
                    Reference=line.rstrip()
                else:
                    Read1=path+"/"+a[0]
                    Read2=path+"/"+a[1]
                    print(Read1,Read2)
                    Encoding=subprocess.run(['sh', 'CheckEncoding.sh',Read1], stdout=subprocess.PIPE).stdout.decode('utf-8')
                    if "33" in Encoding:
                        Encoding="sanger"
                    elif "ed+64" in Encoding:
                        Encoding="illumina"
                    elif "Solexa" in Encoding:
                        Encoding="solexa"
                    else:
                        print("ERROR:unknown encoding")
                        sys.exit()
                    print("fastq encoding for the sample "+self.sample+" is "+Encoding)
                    dest.write("sickle pe -f "+Read1+" -r "+Read2+" -t "+Encoding+" -l 50 -g -o "+"./"+\
                    RunDir+"/trimmed_"+a[0]+" -p "+"./"+RunDir+"/trimmed_"+a[1]+" -s "+"./"+RunDir+"/trUnpair_"+self.sample+".fq.gz")
                    dest.write("\n")
                    dest.write("echo \"running fastqc for sample "+self.sample+"\"")
                    dest.write("\n")
                    dest.write("mkdir ./"+RunDir+"/fastqc_out_"+self.sample+"_"+str(cnt))
                    dest.write("\n")
                    dest.write("fastqc -o ./"+RunDir+"/fastqc_out_"+self.sample+"_"+str(cnt)+" -t 16 "+"./"+RunDir+"/trimmed_"+a[0])
                    dest.write("\n")
                    dest.write("fastqc -o ./"+RunDir+"/fastqc_out_"+self.sample+"_"+str(cnt)+" -t 16 "+"./"+RunDir+"/trimmed_"+a[1])
                    dest.write("\n")
                    dest.write("echo \"finished running fastqc\"")
                    dest.write("\n")
                    dest.write("cd ./"+RunDir+"/fastqc_out_"+self.sample+"_"+str(cnt))
                    dest.write("\n")
                    dest.write("for i in $(ls *.zip);do unzip $i;done")
                    dest.write("\n")
                    dest.write("find . -type f -name \"summary.txt\"|while read fname;do awk \'$0~/FAIL/ || $0~/WARN/{print $0}\' $fname>>../"+ \
                    "fastqc_summary_fail_warn"+str(cnt)+".txt;done")
                    dest.write("\n")
                    dest.write("cd ../../")
                    dest.write("\n")
                    dest.write("echo \"starting bwa alignment\"")
                    dest.write("\n")
                    dest.write("bwa mem -t 16"+" -M -R '@RG\\tID:sample_"+str(cnt)+"\\tSM:"+self.sample+"\\tPL:ILLUMINA' "+ \
                    Reference+" ./"+RunDir+"/trimmed_"+a[0]+" ./"+RunDir+"/trimmed_"+a[1]+" >./"+RunDir+"/"+self.sample+"_"+str(cnt)+".sam")
                    dest.write("\n")
                    dest.write("samtools view -@ 16 -Shb "+"./"+RunDir+"/"+self.sample+"_"+str(cnt)+".sam >"+"./"+RunDir+"/"+\
                    self.sample+"_"+str(cnt)+".bam")
                    dest.write("\n")
                    dest.write("rm "+"./"+RunDir+"/"+self.sample+"_"+str(cnt)+".sam")
                    dest.write("\n")
                    dest.write("echo start sorting")
                    dest.write("\n")
                    dest.write("samtools sort -@ 16 -m 1500M -O BAM "+"./"+RunDir+"/"+self.sample+"_"+str(cnt)+".bam"+" -o ./"+RunDir+"/"+self.sample+"_"+str(cnt)+".sorted.bam")
                    dest.write("\n")
                    dest.write("rm ./"+RunDir+"/"+self.sample+"_"+str(cnt)+".bam")
                    dest.write("\n")
        dest.write("#total number of bam files "+str(cnt))
        dest.write("\n")
        if cnt>1:
            BamList=[]
            for i in range(1,cnt+1):
                if i==1:
                    dest.write("samtools view -H ./"+RunDir+"/"+self.sample+"_"+str(i)+".sorted.bam "+"|sed 's/SM:unknown/SM:"+self.sample+\
                    "/'"+"|sed 's/PL:sanger/PL:ILLUMINA/' >./"+RunDir+"/newheader.txt")
                    dest.write("\n")
                else:
                    dest.write("samtools view -H ./"+RunDir+"/"+self.sample+"_"+str(i)+".sorted.bam "+"|sed 's/SM:unknown/SM:"+self.sample+\
                    "/'"+"|sed 's/PL:sanger/PL:ILLUMINA/'"+"|sed 's/ID:bwa/ID:bwa_"+str(i)+"/'|grep @RG >>./"+RunDir+"/newheader.txt")
                    dest.write("\n")
                BamList.append("./"+RunDir+"/"+self.sample+"_"+str(i)+".sorted.bam")
            dest.write("samtools merge "+"./"+RunDir+"/"+self.sample+".tmp.bam "+" ".join(BamList))
            dest.write("\n")
            dest.write("samtools reheader "+"./"+RunDir+"/newheader.txt "+"./"+RunDir+"/"+self.sample+".tmp.bam >"+"./"+RunDir+"/"+self.sample+"_rh.bam")
            dest.write("\n")
            dest.write("rm "+"./"+RunDir+"/"+self.sample+".tmp.bam")
            dest.write("\n")
        else:
            dest.write("mv ./"+RunDir+"/"+self.sample+"_"+str(cnt)+".sorted.bam "+"./"+RunDir+"/"+self.sample+"_rh.bam")
            dest.write("\n")
        dest.write("echo \"starting duplicate-removals of reads using picard\"")
        dest.write("\n")
        dest.write("java -Xmx8g -jar ~/software/picard/build/libs/picard.jar MarkDuplicates I="+"./"+RunDir+"/"+self.sample+"_rh.bam O="+\
        "./"+RunDir+"/"+self.sample+"_rh.rmDupli.bam M="+self.sample+ \
        "rmDup.metri.txt AS=true REMOVE_DUPLICATES=true")
        dest.write("\n")
        dest.write("##re-alignment using GATK")
        dest.write("\n")
        dest.write("samtools index "+"./"+RunDir+"/"+self.sample+"_rh.rmDupli.bam")
        dest.write("\n")
        dest.write("java -Xmx8g -jar /home/maulik/software/AlignmentPipeline/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar "+\
        "-nt 16 -T RealignerTargetCreator -R "+Reference+" -I "+\
        "./"+RunDir+"/"+self.sample+"_rh.rmDupli.bam -o "+"./"+RunDir+"/"+self.sample+"_rh.intervals")
        dest.write("\n")
        dest.write("java -Xmx8g -jar /home/maulik/software/AlignmentPipeline/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar "+\
        "-T IndelRealigner -R "+Reference+" -I "+\
        "./"+RunDir+"/"+self.sample+"_rh.rmDupli.bam "+"-targetIntervals "+"./"+RunDir+"/"+self.sample+"_rh.intervals -o "+\
        "./"+RunDir+"/"+self.sample+"_rh.rmDupliIndelRealigned.bam")
        dest.write("\n")
        dest.write("samtools index "+"./"+RunDir+"/"+self.sample+"_rh.rmDupliIndelRealigned.bam")
        dest.write("\n")
        dest.write("samtools flagstat "+"./"+RunDir+"/"+self.sample+"_rh.rmDupliIndelRealigned.bam >"+"./"+RunDir+\
        "/BasicStats.txt")
        dest.write("\n")
        dest.write("samtools depth -aa "+"./"+RunDir+"/"+self.sample+"_rh.rmDupliIndelRealigned.bam"+\
        "|awk '{sum+=$3;count++}END{print sum/count}' >"+"./"+RunDir+"/FinalDepth.txt")
        dest.write("\n")
        dest.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This python script will make the slurm script to run the whole genome alignment pipeline.\
    \n The pipeline starts with filtering of WGS data using sickle, followed by quality assement using FastQC. \n For the alignment purpose \
    bwa mem is used, followed by duplicate reads removal using picard and realignment using GATK3.8 pipeline. \n This script also requires \
    CheckEncoding.sh in the same folder as this script")
    parser.add_argument('-f',"--inFile",metavar="FILE",help="For the format refer to example file Example1.txt",required=True)
    parser.add_argument('-s',"--sampleId",metavar="STRING",help="This will be the name of the output files as well as the directory with prefix Dir_",required=True)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        makeScript=ProduceAlignmentScript(args.inFile,args.sampleId)
        makeScript.MakeSlurmFiles()
