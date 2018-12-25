#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @PROJECT : rna-seq-pipeline
# @Time    : 2018/3/30 15:09
# @Author  : Chen Yuelong
# @Mail    : yuelong.chen@oumeng.com.cn
# @File    : utils.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals
import sys, os
import re

sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))


def get_list(file):
    '''
    通过样本list生成dict
    :param file: 样本list文件
    :return: dict
    '''
    sampledict = {}
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # print(line)
            cells = re.split('[ |\t]+', line.strip('\n'))
            # print(cells)
            sampledict[cells[0]] = cells[1:]
    return sampledict


def generate_whole_makefils(args):
    '''

    :param args:
    :return:
    '''
    sampledict = get_list(args.input_list)
    with open(args.makefile, 'w') as f:
        f.write('include $(CONFIG)\n')
        f.write('export PATH:=/Bioinfo/MDDRD2/PMO/chenyl/venv_python/py_RNACocktail/bin:${PATH} \n')
        # f.write('echo $(PATH)\n')
        targ_n = 0
        diff_depends = []
        control_bam = []
        case_bam = []
        control_sample=[]
        case_sample=[]
        ALL=[]
        for sample, infos in sampledict.items():
            target_sample=[]
            targ_n = targ_n + 1
            fq1 = infos.pop(0)
            fq2 = infos.pop(0)

            qcdir ='{}/{}_{}/QC'.format(args.outdir,args.project,sample)
            aligndir ='{}/{}_{}/ALIGN'.format(args.outdir,args.project,sample)
            reconstructdir ='{}/{}_{}/RECONSTRUCT'.format(args.outdir,args.project,sample)
            variantdir ='{}/{}_{}/VARIANT'.format(args.outdir,args.project,sample)
            annodir='{}/{}_{}/ANNO'.format(args.outdir,args.project,sample)
            workdir ='{}/{}_{}/WORKDIR'.format(args.outdir,args.project,sample)


            cleanfq1 ='{}/{}.clean.R1.fastq.gz'.format(qcdir,sample)
            cleanfq2 ='{}/{}.clean.R2.fastq.gz'.format(qcdir,sample)
            cleanupfq1 ='{}/{}.cleanUp.R1.fastq.gz'.format(qcdir,sample)
            cleanupfq2 ='{}/{}.cleanUp.R2.fastq.gz'.format(qcdir,sample)


            mk_target = generate_mkdir_makfile(f,targ_n,qcdir,
                                              aligndir,reconstructdir,
                                              variantdir,annodir,workdir)
            target_sample.append(mk_target)
            qc_target=generate_fastqc_makefile(f,targ_n,qcdir,workdir,fq1,fq2,
                                                 cleanfq1,cleanupfq1,cleanfq2,cleanupfq2,
                                                 mk_target)
            target_sample.append(qc_target)
            align_target=generate_align_makefile(f,targ_n,aligndir,workdir,cleanfq1,cleanfq2,sample,qc_target)
            target_sample.append(align_target)
            align_bam = '{aligndir}/hisat2/{sample}/alignments.sorted.bam'.format(
                aligndir=aligndir,
                sample=sample
            )
            if infos[0] == 'control':
                control_bam.append(align_bam)
                control_sample.append(sample)
            elif infos[0] == 'case':
                case_bam.append(align_bam)
                case_sample.append(sample)
            else:
                raise ValueError('input_list.txt:\n'
                                 'your type is **{}**!!!\n'
                                 'sample type must be either `case` or `control`!!!'.format(infos[0]))
            diff_depends.append(align_target)
            recon_target = generate_reco_makefile(f,targ_n,aligndir,workdir,reconstructdir,sample,align_target)
            variant_target = generate_vari_makefile(f,targ_n,aligndir,workdir,variantdir,sample,align_target)
            target_sample.append(recon_target)
            target_sample.append(variant_target)
            ALL.append(generate_target(f,sample,target_sample))
        targ_n = targ_n + 1
        diffdir = '{}/{}_DIFF'.format(args.outdir, args.project)
        coexpdir = '{}/{}_COEXP'.format(args.outdir, args.project)
        coworkdir = '{}/{}_WORKDIR'.format(args.outdir, args.project)
        comkdir = generate_mkdir_makfile(f, targ_n, diffdir, coexpdir, coworkdir)
        diff_depends.append(comkdir)
        dif_target = generate_diff_makefile(f,targ_n,diffdir,coworkdir,
                                            case_bam,case_sample,control_bam,control_sample,
                                            diff_depends)
        coexp_target = generate_coexp_makefile(f,targ_n,coworkdir,coexpdir,control_sample,case_sample,dif_target)
        ALL.append(dif_target)
        ALL.append(coexp_target)
        alltarget = generate_target(f,'ALL',ALL)
        print(alltarget)




def generate_target(fbuffer,name,depends):
    '''

    :param fbuffer:
    :param name:
    :param depends:
    :return:
    '''
    deps=' '.join(depends)
    fbuffer.write('{}:{}\n'.format(name,deps))
    return name




def generate_mkdir_makfile(fbuffer,targ_n,*dir):
    '''

    :param fbuffer:
    :param targ_n:
    :param dir:
    :return:
    '''
    name = 'mkdir{}'.format(targ_n)
    fbuffer.write('{}:\n'.format(name))
    # print('{}:'.format(name))
    for i in dir:
        fbuffer.write('\tmkdir -p {}\n'.format(i))
        # print('\tmkdir -p {}'.format(i))
    return name



def generate_fastqc_makefile(fbuffer,targ_n,qcdir,workdir,fq1,fq2,clean1,cleanup1,clean2,cleanup2,*depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param qcdir:
    :param workdir:
    :param fq1:
    :param fq2:
    :param clean1:
    :param cleanup1:
    :param clean2:
    :param cleanup2:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'fastqc{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name,deps))
    fbuffer.write('\t$(FASTQC) {fq1} {fq2} --outdir {qcdir} --threads 1 -extract -q -f fastq --java $(JAVA)\n'.
                  format(fq1=fq1,fq2=fq2,qcdir=qcdir))
    fbuffer.write('\t$(JAVA) -jar -XX:ParallelGCThreads=4 -Xmx10g $(TRIMMOMATIC) PE {fq1} {fq2} '
                  '{clean1} {cleanup1} {clean2} {cleanup2} -trimlog {workdir}/fastqc.trim.log '
                  '-threads 4 -phred33 ILLUMINACLIP:$(TRUSEQPE):2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
                  'MINLEN:36\n'.
                  format(fq1=fq1,fq2=fq2,clean1=clean1,clean2=clean2,cleanup1=cleanup1,cleanup2=cleanup2,
                         workdir=workdir))
    return name



def generate_align_makefile(fbuffer,targ_n,aligndir,workdir,clean1,clean2,sample,*depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param aligndir:
    :param workdir:
    :param clean1:
    :param clean2:
    :param sample:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'align{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name,deps))
    fbuffer.write('\t$(RNACOCKTAIL) align --align_idx $(REFHISAT2) --outdir {aligndir} '
                  '--workdir {workdir} --ref_gtf $(GENEGTF) --1 {clean1} --2 {clean2} '
                  '--hisat2 $(HISAT2) --hisat2_sps $(HISAT2SPS) --samtools $(SAMTOOLS) '
                  '	--threads 10 --sample {sample}\n'.format(
        aligndir=aligndir,
        workdir=workdir,
        clean1=clean1,
        clean2=clean2,
        sample=sample
    ))
    return name


def generate_reco_makefile(fbuffer,targ_n,aligndir,workdir,reconstructdir,sample,*depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param aligndir:
    :param workdir:
    :param reconstructdir:
    :param sample:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'reconstruct{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name, deps))
    fbuffer.write('\t$(RNACOCKTAIL) reconstruct --alignment_bam {aligndir}/hisat2/{sample}/alignments.sorted.bam '
                  '--outdir {reconstructdir} --workdir {workdir} --ref_gtf $(GENEGTF) --stringtie $(STRINGTIE) '
                  '--threads 10 --sample {sample}\n'.format(
        aligndir=aligndir,
        sample=sample,
        workdir=workdir,
        reconstructdir=reconstructdir
    ))
    return name


def generate_vari_makefile(fbuffer,targ_n,aligndir,workdir,variantdir,sample,*depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param aligndir:
    :param workdir:
    :param variantdir:
    :param sample:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'variant{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name, deps))
    fbuffer.write('\t$(RNACOCKTAIL) variant --alignment {aligndir}/hisat2/{sample}/alignments.sorted.bam '
                  '--CleanSam --outdir {variantdir} --workdir {workdir} --picard $(PICARD) --gatk $(GATK) '
                  '--threads 10 --sample {sample} --ref_genome $(REF) --IndelRealignment --knownsites $(DBSNP) '
                  '--java $(JAVA)\n'.format(
        aligndir=aligndir,
        sample=sample,
        workdir=workdir,
        variantdir=variantdir
    ))
    return name


def generate_anno_makefile(fbuffer,targ_n,aligndir,workdir,reconstructdir,sample,*depend):
    deps = ' '.join(depend)
    name = 'anno{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name, deps))


def generate_diff_makefile(fbuffer,targ_n,diffdir,coworkdir,case_bam,case_sample,control_bam,control_sample,depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param diffdir:
    :param coworkdir:
    :param case_bam:
    :param case_sample:
    :param control_bam:
    :param control_sample:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'diff{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name, deps))
    fbuffer.write('\t$(RNACOCKTAIL) diff --alignments {abams} {bbams} '
                  '--sample {anames} {bnames} --ref_gtf $(GENEGTF) --outdir {diffdir} --workdir {workdir} '
                  '--featureCounts $(FEATURECOUNTS) --R $(R) --stringtie $(STRINGTIE)\n'.format(
        abams=','.join(control_bam),
        bbams=','.join(case_bam),
        anames=','.join(control_sample),
        bnames=','.join(case_sample),
        diffdir=diffdir,
        workdir=coworkdir
    ))
    fbuffer.write('\t$(COR) $(PLOTR) {workdir}/deseq2/{samples}/deseq2.rda 0.1 {diffdir}/deseq2/{samples}\n'.format(
        workdir=coworkdir,
        diffdir=diffdir,
        samples='-'.join([','.join(control_sample), ','.join(case_sample)])
    ))
    return name


def generate_coexp_makefile(fbuffer,targ_n,coworkdir,coexpdir,control_sample,case_sample,*depend):
    '''

    :param fbuffer:
    :param targ_n:
    :param coworkdir:
    :param coexpdir:
    :param control_sample:
    :param case_sample:
    :param depend:
    :return:
    '''
    deps = ' '.join(depend)
    name = 'coexp{}'.format(targ_n)
    fbuffer.write('{}:{}\n'.format(name, deps))
    fbuffer.write('\t$(COR) $(COSEQR) {coworkdir}/deseq2/{samples}/deseq2.rda {coexpdir} 0.1\n'.format(
        coworkdir=coworkdir,
        coexpdir=coexpdir,
        samples='-'.join([','.join(control_sample),','.join(case_sample)])
    ))
    return name




def main():
    '''
    测试流程
    '''
    # get_args()
    # test = get_list('/Bioinfo/MDDRD2/PMO/chenyl/PROJECTS/RNA-seq/examples/input_list.txt')
    # for i, j in test.items():
    #     print('{}:{}'.format(i, j))
    # pass
    # a = 1
    # b = 2
    test=generate_mkdir_makfile('', 2, 'test1','test2','test3')
    print(test)


if __name__ == '__main__':
    main()
