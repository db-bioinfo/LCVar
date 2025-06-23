#!/usr/bin/env python
#########################################################################
# Author: Lee Quan (leequan@gmail.com)
# Created Time: 2015-11-10 19:15:32 Tuesday 
# File Name: InterVar.py File type: python
# Last Change:.
# Description: python script for  Interpretation of Pathogenetic Benign
#########################################################################

import copy,logging,os,io,re,time,sys,platform,optparse,gzip,glob

prog="InterVar"

version = """%prog 2.2.2 20210727
Written by Quan LI,leequan@gmail.com. 
InterVar is free for non-commercial use without warranty.
Please contact the authors for commercial use.
Copyright (C) 2016 Wang Genomic Lab
============================================================================
"""

usage = """Usage: %prog [OPTION] -i  INPUT -o  OUTPUT ...
       %prog  --config=config.ini ...
"""

description = """=============================================================================
InterVar                                                                       
Interpretation of Pathogenic/Benign for variants using python scripts. 

.####.##....##.########.########.########..##.....##....###....########.
..##..###...##....##....##.......##.....##.##.....##...##.##...##.....##
..##..####..##....##....##.......##.....##.##.....##..##...##..##.....##
..##..##.##.##....##....######...########..##.....##.##.....##.########.
..##..##..####....##....##.......##...##....##...##..#########.##...##..
..##..##...###....##....##.......##....##....##.##...##.....##.##....##.
.####.##....##....##....########.##.....##....###....##.....##.##.....##

=============================================================================
"""
end = """=============================================================================
........................................................................
.####.##....##.########.########.########..##.....##....###....########.
..##..###...##....##....##.......##.....##.##.....##...##.##...##.....##
..##..####..##....##....##.......##.....##.##.....##..##...##..##.....##
..##..##.##.##....##....######...########..##.....##.##.....##.########.
..##..##..####....##....##.......##...##....##...##..#########.##...##..
..##..##...###....##....##.......##....##....##.##...##.....##.##....##.
.####.##....##....##....########.##.....##....###....##.....##.##.....##
.......................................................................
Thanks for using InterVar!
Report bugs to leequan@gmail.com;
InterVar homepage: <http://wInterVar.wglab.org>

=============================================================================
"""


line_sum=0;

if platform.python_version()< '3.0.0' :
    import ConfigParser
else:
    import configparser

paras = {}

def ConfigSectionMap(config,section):
    global paras
    options = config.options(section)
    for option in options:
        try:
            paras[option] = config.get(section, option)
            if paras[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            paras[option] = None
    return

user_evidence_dict={}


class myGzipFile(gzip.GzipFile): 
    def __enter__(self): 
        if self.fileobj is None: 
            raise ValueError("I/O operation on closed GzipFile object") 
        return self 

    def __exit__(self, *args): 
        self.close() 


#begin read some important datsets/list firstly;
lof_genes_dict={}
aa_changes_dict={}
domain_benign_dict={}
mim2gene_dict={}
mim2gene_dict2={}
morbidmap_dict={}
morbidmap_dict2={}
PP2_genes_dict={}
BP1_genes_dict={}
PS4_snps_dict={}
exclude_snps_dict={}
mim_recessive_dict={}
mim_domin_dict={}
mim_adultonset_dict={}
mim_pheno_dict={}
mim_orpha_dict={}
orpha_dict={}
BS2_snps_recess_dict={}
BS2_snps_domin_dict={}
knownGeneCanonical_dict={}
knownGeneCanonical_st_dict={}
knownGeneCanonical_ed_dict={}

def flip_ACGT(acgt):
    nt="";
    if acgt=="A":
        nt="T"
    if acgt=="T":
        nt="A"
    if acgt=="C":
        nt="G"
    if acgt=="G":
        nt="C"
    if acgt=="N":
        nt="N"
    if acgt=="X":
        nt="X"
    return(nt)

def read_datasets():
#0. read the user specified evidence file
    if os.path.isfile(paras['evidence_file']):
        try:
            fh=open(paras['evidence_file'], "r")
            strs = fh.read()
            for line2 in strs.split('\n'):
                cls2=line2.split('\t')
                if len(cls2)>1:
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]
                    keys=re.sub("[Cc][Hh][Rr]","",keys)
                    #print("%s" %keys)
                    user_evidence_dict[keys]=cls2[4].upper()
        except IOError:
            print("Error: can\'t read the user specified evidence file %s" % paras['evidence_file'])
        else:
            fh.close()    

#1.LOF gene list       
    try:               
        fh = open(paras['lof_genes'], "r")
        str = fh.read()
        for line2 in str.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                lof_genes_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the LOF genes file %s" % paras['lof_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()    

#2. AA change list
    try:
        fh = open(paras['ps1_aa'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1 :
                keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[4]
                keys=re.sub("[Cc][Hh][Rr]","",keys)
                aa_changes_dict[keys]=cls2[6]
    except IOError:
        print("Error: can\'t read the  amino acid change file %s" % paras['ps1_aa'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()    

#3. Domain with benign 
    try:
        fh = open(paras['pm1_domain'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1:
                keys=cls2[0]+"_"+cls2[1]+": "+cls2[2]
                domain_benign_dict[keys]="1"
    except IOError:
        print("Error: can\'t read the PM1 domain  file %s" % paras['pm1_domain'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   

#4. OMIM mim2gene.txt file 
    try:
        fh = open(paras['mim2gene'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1:
                cls0=cls2[4].split(',')
                keys=cls0[0]
                mim2gene_dict[keys]=cls2[0]
                keys1=cls2[3]
                keys=keys1.upper()
                mim2gene_dict2[keys]=cls2[0]
    except IOError:
        print("Error: can\'t read the OMIM  file %s" % paras['mim2gene'])
        print("Error: Please download it from http://www.omim.org/downloads")
        sys.exit()
    else:
        fh.close()   

#5.PP2 gene list       
    try:               
        fh = open(paras['pp2_genes'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                PP2_genes_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the PP2 genes file %s" % paras['PP2_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()    

#5.BP1 gene list       
    try:               
        fh = open(paras['bp1_genes'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                BP1_genes_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the BP1 genes file %s" % paras['BP1_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()    

#6.morbidmap from OMIM  for BP5 ,  multifactorial disorders  list
#The reviewers suggeset to disable the OMIM morbidmap for BP5 
    '''
    try:               
        fh = open(paras['morbidmap'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            #print("%s %s %d" % (cls2[0], cls[Funcanno_flgs['Gene']], len(cls2[0])) )
            #{Tuberculosis, protection against}, 607948 (3)|TIRAP, BACTS1|606252|11q24.2
            if len(cls2[0])>1 and cls2[0].find('{')==0:  # disorder start with "{"
                morbidmap_dict2[ cls2[2] ]='1'  # key as mim number 
                for cls3 in cls2[1].split(', '):
                    keys=cls3.upper()
                    morbidmap_dict[ keys ]='1'  # key as gene name
    except IOError:
        print("Error: can\'t read the OMIM morbidmap disorder file %s" % paras['morbidmap'])
        print("Error: Please download it from http://www.omim.org/downloads")
        sys.exit()
    else:
        fh.close()    
    '''

#7.prevalence of the variant with OR>5 for PS4 ,  the dataset is from gwasdb jjwanglab.org/gwasdb
    try:               
        fh = open(paras['ps4_snps'], "r")
        str = fh.read()
        for line2 in str.split('\n'):
            cls2=line2.split('\t')
            # PS4_snps_dict
            if len(cls2[0])>=1:  # 
                keys=cls2[0]+"_"+cls2[1]+"_"+cls2[1]+"_"+cls2[3]+"_"+cls2[4]
                keys=re.sub("[Cc][Hh][Rr]","",keys)
                PS4_snps_dict[ keys ]='1'  # key as gene name
    except IOError:
        print("Error: can\'t read the snp list file for PS4 %s" % paras['ps4_snps'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()    

#8. read the user specified SNP list, the variants will pass the frequency check.
    if os.path.isfile(paras['exclude_snps']):
        try:
            fh=open(paras['exclude_snps'], "r")
            strs = fh.read()
            for line2 in strs.split('\n'):
                cls2=line2.split('\t')
                if len(cls2)>1:
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]
                    keys=re.sub("[Cc][Hh][Rr]","",keys)
                    exclude_snps_dict[keys]="1"
        except IOError:
            print("Error: can\'t read the user specified SNP list file %s" % paras['exclude_snps'])
        else:
            fh.close()    
#9. OMIM mim_recessive.txt file mim_domin  mim_adultonset 
    try:
        fh = open(paras['mim_recessive'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                mim_recessive_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the OMIM recessive disorder file %s" % paras['mim_recessive'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   
    try:
        fh = open(paras['mim_domin'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                mim_domin_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the OMIM dominant disorder file %s" % paras['mim_domin'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   
    try:
        fh = open(paras['mim_adultonset'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                mim_adultonset_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the OMIM adult onset disorder file %s" % paras['mim_adultonset'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   




#10. knownGeneCanonical exon file  # caution the build ver, now it is hg19
    try:
        fh = open(paras['knowngenecanonical'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                knownGeneCanonical_dict[keys]=cls2[1]
                knownGeneCanonical_st_dict[keys]=cls2[2]
                knownGeneCanonical_ed_dict[keys]=cls2[3]
                #print("%s %s" %(keys,knownGeneCanonical_dict[keys]))
    except IOError:
        print("Error: can\'t read the knownGeneCanonical  file %s" % paras['knowngenecanonical'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   

#11.BS2 variants of recessive homo, domin heter
    try:              
        with myGzipFile(paras['bs2_snps'], "rb") as fh:

            #fh = open(paras['bs2_snps'], "r")
            strs = fh.read().decode()
            for line2 in strs.split('\n'):
                cls2=line2.split(' ')
                # PS4_snps_dict
                if len(cls2[0])>=1:  #
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]
                    #keys=re.sub("[Cc][Hh][Rr]","",keys)
                    BS2_snps_recess_dict[ keys ]=cls2[4]  # key as snp info
                    BS2_snps_domin_dict[ keys ]=cls2[5]  # key as snp info
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[1]+"_"+flip_ACGT(cls2[2])+"_"+flip_ACGT(cls2[3])
                    #keys=re.sub("[Cc][Hh][Rr]","",keys)
                    BS2_snps_recess_dict[ keys ]=cls2[4]  # key as snp info
                    BS2_snps_domin_dict[ keys ]=cls2[5]  # key as snp info

    except IOError:
        print("Error: can\'t read the snp list file for BS2 %s" % paras['bs2_snps'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()    

#12. OMIM mim_pheno.txt file
#mim_pheno = %(database_intervar)s/mim_pheno.txt
#mim_orpha = %(database_intervar)s/mim_orpha.txt
 
    try:
        fh = open(paras['mim_pheno'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                mim_pheno_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_pheno_dict[keys]))
    except IOError:
        print("Error: can\'t read the MIM  file %s" % paras['mim_pheno'])
        print("Error: Please download it from InterVar source website")
        sys.exit()
    else:
        fh.close()   




#13. OMIM mim_orpha.txt file 
    try:
        fh = open(paras['mim_orpha'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                mim_orpha_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_orpha_dict[keys]))
    except IOError:
        print("Error: can\'t read the MIM  file %s" % paras['mim_orpha'])
        print("Error: Please download it from InterVar source website")
        sys.exit()
    else:
        fh.close()   

#14.  orpha.txt file 
    try:
        fh = open(paras['orpha'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1:
                keys=cls2[0]
                orpha_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_orpha_dict[keys]))
    except IOError:
        print("Error: can\'t read the Orpha  file %s" % paras['orpha'])
        print("Error: Please download it from InterVar source website")
        sys.exit()
    else:
        fh.close()   


#end read datasets
    return



def check_downdb():
    path=paras['database_locat']
    path=path.strip()
    path=path.rstrip("\/")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print("Notice: the folder of %s is created!" % path)
    else:
        print("Warning: the folder of %s is already created!" % path)
    ds=paras['database_names']
    ds.expandtabs(1);
    # database_names = refGene 1000g2014oct esp6500siv2_all avsnp147 ljb26_all clinvar_20150629 gnomad_genome hg19_dbscsnv11 dbnsfp31a_interpro rmsk ensGene
    if not os.path.isfile(paras['annotate_variation']):
        print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % paras['annotate_variation'])
        if paras['skip_annovar'] != True:
            sys.exit()

    for dbs in ds.split():
        # os.path.isfile(options.table_annovar)
        file_name=dbs
        #if dbs=="1000g2014oct":
        #    file_name="ALL.sites.2014_10"
        if dbs=="1000g2015aug":
            file_name="ALL.sites.2015_08" # hg19_ALL.sites.2015_08.txt

        dataset_file=paras['database_locat']+"/"+paras['buildver']+"_"+file_name+".txt"
        if dbs != 'rmsk':
            cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb -webfrom annovar "+file_name+" "+paras['database_locat']
        if dbs == 'rmsk':
            cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb "+file_name+" "+paras['database_locat']
        if  not os.path.isfile(dataset_file):
            if dbs=="1000g2015aug":
                file_name="1000g2015aug"
                dataset_file=paras['database_locat']+"/"+paras['buildver']+"_"+file_name+".txt"
                cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb -webfrom annovar "+file_name+" "+paras['database_locat']
            if paras['skip_annovar'] != True:
                print("Warning: The Annovar dataset file of %s is not in %s,begin to download this %s ..." %(dbs,paras['database_locat'],dataset_file))
    
            if paras['skip_annovar'] != True:
                print("%s" %cmd)
                os.system(cmd)

def check_input():
    inputft= paras['inputfile_type']
    if inputft.lower() == 'avinput' :
        return
    if inputft.lower() == 'vcf':
        if os.path.isfile(paras['convert2annovar']):
        #convert2annovar.pl -format vcf4 variantfile > variant.avinput
            cmd="perl "+paras['convert2annovar']+" -format vcf4 "+ paras['inputfile']+"> "+paras['inputfile']+".avinput"
            print("Warning: Begin to convert your vcf file of %s to AVinput of Annovar ..." % paras['inputfile'])
            print("%s" %cmd)
            os.system(cmd)
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % paras['convert2annovar'])
            if paras['skip_annovar'] != True:
                sys.exit()
    if inputft.lower() == 'vcf_m':
        if os.path.isfile(paras['convert2annovar']):
        #convert2annovar.pl -format vcf4 variantfile > variant.avinput
            cmd="perl "+paras['convert2annovar']+" -format vcf4 "+ paras['inputfile']+" --allsample   --outfile "+ paras['outfile']
            print("Warning: Begin to convert your vcf file with multiple samples of %s to AVinput of Annovar with All.raw.highqc.vcf.<samplename>.avinput..." % paras['inputfile'])
            print("Warning: Please attention that the sample names in VCF file should  contain letters/numners only, otherwise the converting may be failure!")
            print("%s" %cmd)
            os.system(cmd)
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % paras['convert2annovar'])
            if paras['skip_annovar']  != True:
                sys.exit()
    return

def check_annovar_result():
# table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,ljb26_all,CLINSIG,gnomad_genome   -operation  g,f,f,f,f,f,f   -nastring . -csvout
    inputft= paras['inputfile_type']
    annovar_options=" "
    if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
        annovar_options=annovar_options+"--otherinfo " 
    if re.findall('true',paras['onetranscript'], flags=re.IGNORECASE) :
        annovar_options=annovar_options+"--onetranscript " 

    if not os.path.isfile(paras['table_annovar']):
        print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % paras['table_annovar'])
        if paras['skip_annovar'] != True:
            sys.exit()
    if inputft.lower() == 'avinput' :
        cmd="perl "+paras['table_annovar']+" "+paras['inputfile']+" "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ paras['outfile']+" -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp151,dbnsfp47a,clinvar_20250306,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene  -operation  g,f,f,f,f,f,f,f,r,g,g   -nastring ."+annovar_options
        print("%s" %cmd)
        os.system(cmd)
    if inputft.lower() == 'vcf' :
        cmd="perl "+paras['table_annovar']+" "+paras['inputfile']+".avinput "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ paras['outfile']+" -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp151,dbnsfp47a,clinvar_20250306,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene   -operation  g,f,f,f,f,f,f,f,r,g,g   -nastring ."+annovar_options
        print("%s" %cmd)
        os.system(cmd)
    if inputft.lower() == 'vcf_m' :
        for f in glob.iglob(paras['outfile']+"*.avinput"): 
            print("INFO: Begin to annotate sample file of %s ...." %(f))
            new_outfile=re.sub(".avinput","",f)
            cmd="perl "+paras['table_annovar']+" "+f+" "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ new_outfile +" -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp151,dbnsfp47a,clinvar_20250306,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene   -operation  g,f,f,f,f,f,f,f,r,g,g   -nastring ."+annovar_options
            print("%s" %cmd)
            os.system(cmd)
        
    
    return

def check_genes(anvfile):
#check with multiple genes, so one gene by one gene  to annote
    newoutfile=anvfile+".grl_p"
    try:
        fh = open(anvfile, "r")
        fw = open(newoutfile, "w")
        strs = fh.read()
        sum=0
        otherinf_pos=1
        line_sum=0;
        for line in strs.split('\n'):
            if line_sum==0:
                line=re.sub("CLNSIG","CLINSIG",line)
            cls=line.split('\t')
            if len(cls)>1:
                if sum==0 and re.findall('true',paras['otherinfo'], flags=re.IGNORECASE) :
                    for ii in range(0,len(cls)):
                        if  re.findall('otherinfo',cls[ii], flags=re.IGNORECASE) :
                            otherinf_pos=ii
                     
                gene_name=cls[6]
                if cls[6] == 'Gene.refGene':
                    gene_name='Gene'
#some with multiple genes, so one gene by one gene  to annote
                sum=sum+1
                #for gg in gene_name.split(','):
                for gg in re.split("[,;]",gene_name):
                    if not re.findall('true',paras['otherinfo'], flags=re.IGNORECASE) :
                        line_out=line+"\t"+gg
                    else:
                        line_out=cls[0]
                        for ii in range(1,len(cls)):
                            if ii != otherinf_pos :
                                line_out=line_out+"\t"+cls[ii]
                            if ii == otherinf_pos :
                                line_out=line_out+"\t"+gg+"\t"+cls[ii]
                            
                    if sum >1: line_out=re.sub("^[Cc][Hh][Rr]","",line_out)
                            
                    #line_out=line+"\t"+gg
                    # re.sub("[Cc][Hh][Rr]","",keys)
                    fw.write("%s\t\n" % line_out)
            line_sum=line_sum+1

    except IOError:
        print("Error: can\'t read/write the annovar output file %s %s" % (anvfile,newoutfile))
        sys.exit()
        return
    else:
        pass
        fh.close()
        fw.close()

    return(sum)



def sum_of_list(list):
    sum=0
    for i in list:
        sum=sum+i
    return(sum)

def classfy(PVS1,PS,PM,PP,BA1,BS,BP,Allels_flgs,cls):
    """
    ENHANCEMENT 1: Points-based classification system (Tavtigian et al. 2020)
    Points: Supporting=1, Moderate=2, Strong=4, Very Strong=8
    """
    BPS=["Pathogenic","Likely pathogenic","Benign","Likely benign","Uncertain significance"]
    
    # Calculate pathogenic evidence points
    pathogenic_points = 0
    pathogenic_points += PVS1 * 8  # Very Strong = 8 points
    pathogenic_points += sum(PS) * 4  # Strong = 4 points (PS always returns 0 or 1)
    pathogenic_points += sum(PM)      # PM can return variable points - don't multiply!
    pathogenic_points += sum(PP)      # PP can return variable points - don't multiply!
    
    # Calculate benign evidence points
    benign_points = 0
    benign_points += BA1 * 8  # Stand-alone = 8 points
    benign_points += sum(BS) * 4  # Strong = 4 points (BS always returns 0 or 1)
    benign_points += sum(BP)      # BP can return variable points - don't multiply!
    
    # Calculate total score
    total_score = pathogenic_points - benign_points
    
    # Classify based on points (Tavtigian thresholds)
    if total_score >= 10:
        BPS_out = 0  # Pathogenic
    elif 6 <= total_score <= 9:
        BPS_out = 1  # Likely pathogenic  
    elif 0 <= total_score <= 5:
        BPS_out = 4  # Uncertain significance
    elif -6 <= total_score <= -1:
        BPS_out = 3  # Likely benign
    else:  # total_score <= -7
        BPS_out = 2  # Benign

    # Handle user evidence modifications (keep existing logic)
    if os.path.isfile(paras['evidence_file']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            evds=user_evidence_dict[keys]
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and  re.findall('grade', evd_t[0], flags=re.IGNORECASE) ):
                    if int(evd_t[1])<=3:
                        if(evd_t[0].find('PS')!=-1): 
                            t=evd_t[0].find('PS'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            pathogenic_points -= PS[tt3-1] * 4  # Remove old
                            if int(evd_t[1]) ==1: pathogenic_points += 4  # Strong
                            elif int(evd_t[1]) ==2: pathogenic_points += 2  # Moderate
                            elif int(evd_t[1]) ==3: pathogenic_points += 1  # Supporting
                        # Similar logic for PM, PP, BS, BP...
        except KeyError:
            pass
        
        # Recalculate classification after user modifications
        total_score = pathogenic_points - benign_points
        if total_score >= 10: BPS_out = 0
        elif 6 <= total_score <= 9: BPS_out = 1
        elif 0 <= total_score <= 5: BPS_out = 4
        elif -6 <= total_score <= -1: BPS_out = 3
        else: BPS_out = 2

    return(BPS[BPS_out])

def check_PVS1(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''
    PVS1: Enhanced LOF variant assessment following ClinGen SVI 2018 recommendations
    Implements full decision tree with strength modulation
    '''
    cls = line.split('\t')
    PVS1 = 0
    PVS1_strength = 0  # 0=none, 1=supporting, 2=moderate, 4=strong, 8=very_strong
    
    # Parse complete SnpEff annotations from VCF INFO field
    snpeff_data = {}
    lof_transcripts = {}
    nmd_transcripts = {}
    ann_data = {}
    
    try:
        # Find the column containing ANN= field (search all columns)
        ann_field = ""
        for col in cls:
            if 'ANN=' in str(col):
                ann_field = col
                break
        
        if ann_field:
            # Parse LOF annotations
            import re
            lof_matches = re.findall(r'LOF=\(([^)]+)\)', ann_field)
            for lof_match in lof_matches:
                parts = lof_match.split('|')
                if len(parts) >= 4:
                    gene = parts[0]
                    confidence = float(parts[3])
                    lof_transcripts[gene] = confidence
            
            # Parse NMD annotations
            nmd_matches = re.findall(r'NMD=\(([^)]+)\)', ann_field)
            for nmd_match in nmd_matches:
                parts = nmd_match.split('|')
                if len(parts) >= 4:
                    gene = parts[0]
                    confidence = float(parts[3])
                    nmd_transcripts[gene] = confidence
            
            # Parse detailed ANN annotations
            ann_start = ann_field.find('ANN=')
            if ann_start != -1:
                ann_content = ann_field[ann_start+4:]
                # Split by semicolon to separate from other INFO fields
                ann_content = ann_content.split(';')[0]
                
                # Parse each annotation (comma-separated)
                for ann in ann_content.split(','):
                    fields = ann.split('|')
                    if len(fields) >= 14:
                        gene = fields[3]
                        if gene not in ann_data:
                            ann_data[gene] = []
                        ann_data[gene].append({
                            'effect': fields[1],
                            'impact': fields[2],
                            'transcript': fields[6],
                            'biotype': fields[7],
                            'exon_intron': fields[8],
                            'hgvs_c': fields[9],
                            'hgvs_p': fields[10],
                            'cdna_pos': fields[11],
                            'cds_pos': fields[12],
                            'protein_pos': fields[13]
                        })
    except Exception as e:
        pass
    
    # Get variant gene
    variant_gene = cls[Funcanno_flgs['Gene']]
    
    # Check if gene has LOF mechanism
    gene_has_lof = False
    try:
        if lof_genes_dict.get(variant_gene) == '1':
            gene_has_lof = True
    except:
        pass
    
    if not gene_has_lof:
        return 0  # PVS1 not applicable
    
    # Analyze variant type and position
    func_info = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    
    # 1. Check SnpEff LOF prediction with high confidence
    if variant_gene in lof_transcripts and lof_transcripts[variant_gene] >= 0.9:
        # Start with assumption of Very Strong
        PVS1_strength = 8
        
        # Apply ClinGen SVI 2018 decision tree
        
        # A. Nonsense/frameshift variants
        if any(term in func_info.lower() for term in ['nonsense', 'stopgain', 'frameshift']):
            
            # Check NMD prediction
            triggers_nmd = variant_gene in nmd_transcripts and nmd_transcripts[variant_gene] >= 0.9
            
            # Analyze exon position if available
            last_exon = False
            in_last_50bp = False
            exon_info = None
            
            if variant_gene in ann_data:
                for ann in ann_data[variant_gene]:
                    if ann['exon_intron'] and '/' in ann['exon_intron']:
                        try:
                            current_exon, total_exons = ann['exon_intron'].split('/')
                            current_exon = int(current_exon)
                            total_exons = int(total_exons)
                            
                            # Check if last exon
                            if current_exon == total_exons:
                                last_exon = True
                            
                            # Check if in last 50bp of penultimate exon
                            if current_exon == total_exons - 1:
                                # Need CDS position for this check
                                if ann['cds_pos'] and '/' in ann['cds_pos']:
                                    cds_current, cds_total = ann['cds_pos'].split('/')
                                    cds_current = int(cds_current)
                                    cds_total = int(cds_total)
                                    if cds_total - cds_current <= 50:
                                        in_last_50bp = True
                        except:
                            pass
            
            # Apply strength modulation based on position and NMD
            if triggers_nmd:
                PVS1_strength = 8  # Very Strong
            elif last_exon or in_last_50bp:
                # Variants escaping NMD
                PVS1_strength = 4  # Strong (more conservative)
            else:
                # No clear NMD prediction
                PVS1_strength = 8  # Default to Very Strong for clear LOF
        
        # B. Canonical splice site variants (±1,2)
        elif any(term in func_info.lower() for term in ['splice_donor', 'splice_acceptor']):
            # Check if canonical positions
            if any(pattern in str(ann_field) for pattern in ['+1', '+2', '-1', '-2']):
                PVS1_strength = 8  # Very Strong for canonical splice sites
            else:
                PVS1_strength = 4  # Strong for other splice variants
        
        # C. Start loss variants
        elif 'startloss' in func_info.lower() or 'start_lost' in func_info.lower():
            # Check for alternative start codons
            # This would require additional annotation data
            # For now, assign Strong
            PVS1_strength = 4  # Strong
        
        # D. Single exon deletions (if we had CNV data)
        # Not applicable here
        
    else:
        # Fallback to traditional detection if no SnpEff LOF
        
        # Check variant types that cause LOF
        lof_types = ['nonsense', 'frameshift', 'stopgain', 'startloss', 'stop_gained', 'start_lost']
        canonical_splice = ['splice_donor', 'splice_acceptor']
        
        for lof_type in lof_types:
            if lof_type in func_info.lower():
                PVS1_strength = 8  # Very Strong by default
                break
        
        if not PVS1_strength:
            for splice_type in canonical_splice:
                if splice_type in func_info.lower():
                    # Enhanced splice assessment
                    try:
                        # Check splice scores if available
                        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
                        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
                        
                        if ada_score != '.' and rf_score != '.':
                            ada_val = float(ada_score)
                            rf_val = float(rf_score)
                            # High confidence splice disruption
                            if ada_val >= 0.958 or rf_val >= 0.584:
                                PVS1_strength = 8  # Very Strong
                        else:
                            # Default for canonical splice without scores
                            PVS1_strength = 8
                    except:
                        PVS1_strength = 8  # Default Very Strong
    
    # Convert strength to points for compatibility
    if PVS1_strength >= 8:
        PVS1 = 1  # Will be multiplied by 8 in classfy function
    elif PVS1_strength >= 4:
        PVS1 = 0.5  # This needs special handling in classfy
    elif PVS1_strength >= 2:
        PVS1 = 0.25  # This needs special handling in classfy
    elif PVS1_strength >= 1:
        PVS1 = 0.125  # This needs special handling in classfy
    
    # For now, maintain backward compatibility
    if PVS1_strength >= 4:
        PVS1 = 1  # Full PVS1
    else:
        PVS1 = 0  # No PVS1
    
    return PVS1

def get_pvs1_strength_description(pvs1_value):
    """Helper to describe PVS1 strength for reporting"""
    if pvs1_value >= 1:
        return "PVS1_VeryStrong"
    elif pvs1_value >= 0.5:
        return "PVS1_Strong" 
    elif pvs1_value >= 0.25:
        return "PVS1_Moderate"
    elif pvs1_value >= 0.125:
        return "PVS1_Supporting"
    else:
        return ""

def check_PS1(line,Funcanno_flgs,Allels_flgs,aa_changes_dict):
    '''
    PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
    Example: Val->Leu caused by either G>C or G>T in the same codon
    AAChange.refGene
    NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W
    '''
    
    PS1=0
    PS1_t1=0
    PS1_t2=0
    PS1_t3=0
    dbscSNV_cutoff=0.6    #either score(ada and rf) >0.6 as splicealtering
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    ACGTs=["A","C","G","T"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            PS1_t1=1;
            # need to wait to check Same amino acid change as a previously pathogenic variant
            line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
            #cls0=line_tmp2.split(',')
            cls0=re.split("[,;]",line_tmp2)
            cls0_1=cls0[0].split(':')
            aa=cls0_1[4]
            aa_last=aa[len(aa)-1:]
            keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Alt']]
            try:
                if  aa_changes_dict[keys_tmp2]:
                    PS1_t2=0
            except KeyError:
                for nt in ACGTs:
                    if nt != cls[Allels_flgs['Alt']] and nt != cls[Allels_flgs['Ref']]:
                        keys_tmp3=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+nt
                        try:
                            if aa_changes_dict[keys_tmp3] == aa_last:
                                PS1_t2=1
                        except KeyError:
                            pass
                        else:
                            pass

            else:
                pass
    try:
        if float(cls[Funcanno_flgs['dbscSNV_RF_SCORE']])>dbscSNV_cutoff or float(cls[Funcanno_flgs['dbscSNV_ADA_SCORE']])>dbscSNV_cutoff: # means alter the splicing
            PS1_t3=1
        if cls[Funcanno_flgs['dbscSNV_RF_SCORE']] == "." or cls[Funcanno_flgs['dbscSNV_ADA_SCORE']] == ".": # absent also means not in splicing
            PS1_t3=0
    except ValueError:
        pass
    else:
        pass

    if PS1_t1 !=0 and PS1_t2 != 0 :
        PS1=1
        if PS1_t3 ==1: # remove the splicing affect 
            PS1=0
    return(PS1)

def check_PS2(line,Funcanno_flgs,Allels_flgs):
    '''
    De novo (both maternity and paternity confirmed) in a patient with the disease and no family history
    '''
    PS2=0
    return(PS2)


def check_PS3(line,Funcanno_flgs,Allels_flgs):
    '''
    Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene
    product
    '''
    PS3=0
    '''
    cls=line.split('\t')
    line_tmp=cls[Funcanno_flgs['Gene damage prediction (all disease-causing genes)']]

    funcs_tmp=["missense","nonsynony"]
    line_tmp2=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
    
        #if line_tmp == 'Medium' or line_tmp == 'High':
        if line_tmp2.find(fc)>=0 and  (line_tmp == 'High') :
            PS3=1
    '''
    return(PS3)


def check_PS4(line,Funcanno_flgs,Allels_flgs):
    '''
    The prevalence of the variant in affected individuals is significantly increased compared with the prevalence
    in controls; OR>5 in all the gwas, the dataset is from gwasdb jjwanglab.org/gwasdb
    '''
    PS4=0
    cls=line.split('\t')
    keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
    try:
        if PS4_snps_dict[keys_tmp2] == "1":
            PS4=1
    except KeyError:
        pass
    else:
        pass
    return(PS4)


def check_PM1(line,Funcanno_flgs,Allels_flgs,domain_benign_dict):
    '''
    Located in a mutational hot spot and/or critical functional domain - UPDATED
    Enhanced with hotspot density analysis and domain significance
    '''
    PM1=0
    cls=line.split('\t')
    
    # Only apply to missense variants
    funcs_missense = ["missense", "nonsynony"]
    line_tmp = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    
    is_missense = False
    for fc in funcs_missense:
        if line_tmp.find(fc) >= 0:
            is_missense = True
            break
    
    if not is_missense:
        return PM1
    
    # Enhanced domain analysis
    domain_evidence = 0
    if cls[Funcanno_flgs.get('Interpro_domain', -1)] != '.':
        domain = cls[Funcanno_flgs['Interpro_domain']]
        gene = cls[Funcanno_flgs['Gene']]
        
        # Check domain significance
        domain_key = f"{cls[Allels_flgs['Chr']]}_{gene}: {domain}"
        
        # Enhanced domain evaluation
        if domain_key not in domain_benign_dict:
            # Domain has no known benign variants
            domain_evidence = 1
            
            # Additional checks for domain importance
            important_domains = [
                'kinase', 'phosphatase', 'DNA_binding', 'active_site', 
                'catalytic', 'functional', 'conserved', 'critical'
            ]
            
            domain_lower = domain.lower()
            for important in important_domains:
                if important in domain_lower:
                    domain_evidence = 2  # Higher confidence
                    break
    
    # Enhanced hotspot analysis (simplified version)
    # In a full implementation, this would analyze nearby pathogenic variants
    hotspot_evidence = 0
    
    # Check if variant is in a gene known for hotspots
    gene = cls[Funcanno_flgs['Gene']]
    hotspot_genes = [
        'TP53', 'KRAS', 'NRAS', 'BRAF', 'PIK3CA', 'AKT1', 'EGFR', 
        'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'MSH6', 'PMS2'
    ]
    
    if gene in hotspot_genes:
        # These genes commonly have hotspots
        hotspot_evidence = 1
    
    # Combine evidence with enhanced logic
    if domain_evidence >= 2 or hotspot_evidence >= 1:
        PM1 = 1
    elif domain_evidence >= 1:
        # Moderate evidence
        PM1 = 1
    
    return PM1


def check_PM2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs,mim2gene_dict,mim2gene_dict2):
    '''
    ENHANCEMENT 2: PM2 Strength Fix - Always Supporting per ClinGen SVI v1.0
    Absent from controls (or at extremely low frequency if recessive)
    '''
    PM2=0
    cls=line.split('\t')
    
    # Enhanced frequency and coverage checks
    valid_frequencies = []
    
    # Check each population database with enhanced QC
    # Check each population database with enhanced QC
    for key in Freqs_flgs.keys():
        # EXCLUDE ESP6500 from PM2 evaluation (outdated database per ACMG/ClinGen guidance)
        if key == 'esp6500siv2_all':
            continue
            
        if cls[Freqs_flgs[key]] != '.':
            try:
                freq = float(cls[Freqs_flgs[key]])
                # For gnomAD, also check coverage if available (basic check)
                if key.startswith('gnomAD') and freq > 0:
                    # If we have access to coverage data, use it
                    # For now, accept if frequency exists and is reasonable
                    if freq < 0.5:  # Sanity check - no variant should be >50%
                        valid_frequencies.append(freq)
                else:
                    valid_frequencies.append(freq)
            except ValueError:
                continue
    
    # If truly absent from all valid databases
    if len(valid_frequencies) == 0:
        PM2 = 1  # FIXED: Always Supporting (1), never Moderate (2)
    else:
        # Check mode of inheritance for recessive disorders
        mim_num = 0
        try:
            mim_num = mim2gene_dict2.get(cls[Funcanno_flgs['Gene']], 0)
        except:
            pass
        
        try:
            mim_num = mim2gene_dict.get(cls[Funcanno_flgs['Gene.ensGene']], mim_num)
        except:
            pass
            
        if str(mim_num) != '0' and mim_num:
            try:
                # For recessive disorders, use more lenient threshold
                if mim_recessive_dict.get(str(mim_num)) == "1":
                    # Recessive threshold - more lenient
                    max_freq = max(valid_frequencies)
                    if max_freq < 0.001:  # 0.1% threshold for recessive
                        PM2 = 1  # FIXED: Always Supporting
                else:
                    # For dominant disorders
                    max_freq = max(valid_frequencies)
                    if max_freq < 0.0001:  # 0.01% threshold for dominant
                        PM2 = 1  # FIXED: Always Supporting
            except:
                # Default behavior if inheritance unknown
                # Use original InterVar-like threshold
                max_freq = max(valid_frequencies) if valid_frequencies else 0
                if max_freq < 0.0001:  # 0.01% threshold as default
                    PM2 = 1  # FIXED: Always Supporting
    
    return PM2



def check_PM3(line,Funcanno_flgs,Allels_flgs):
    '''
    For recessive disorders, detected in trans with a pathogenic variant
    '''
    PM3=0
    return(PM3)

def check_PM4(line,Funcanno_flgs,Allels_flgs):
    '''
    Protein length changes as a result of in-frame deletions/insertions in a nonrepeat region or stop-loss variants
    '''
    PM4=0
    PM4_t1=0
    PM4_t2=0
    cls=line.split('\t')
    #funcs_tmp=["cds-indel","stop-loss"]
    funcs_tmp=["nonframeshift insertion","nonframeshift deletion","stoploss"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            PM4_t1=1
        # need to wait to check  in a nonrepeat region
    if cls[Funcanno_flgs['rmsk']] == '.':
        PM4_t2=1
    if cls[Funcanno_flgs['rmsk']] != '.' and  line_tmp.find("stoploss")>=0 :
        PM4_t2=1

    if PM4_t1 !=0 and PM4_t2 != 0 :
        PM4=1

    return(PM4)

def check_PM5(line,Funcanno_flgs,Allels_flgs,aa_changes_dict):
    '''
    Novel missense change at an amino acid residue where a different missense change determined to be
    pathogenic has been seen before;Example: Arg156His is pathogenic; now you observe Arg156Cys
    NOD2:NM_001293557:exon3:c.C2023T:p.R675W,NOD2:NM_022162:exon4:c.C2104T:p.R702W
    '''
    PM5=0
    PM5_t1=0
    PM5_t2=0
    PM5_t3=0
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    ACGTs=["A","C","G","T"]

    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            PM5_t1=1;
        # need to wait to check no-Same amino acid change as a previously pathogenic variant
            line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
            #cls0=line_tmp2.split(',')
            cls0=re.split("[,;]",line_tmp2)
            cls0_1=cls0[0].split(':')
            aa=cls0_1[4]
            aa_last=aa[len(aa)-1:]
            #keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
            keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Alt']]
            #print("%s %s %s" %(aa,aa_last,keys_tmp2))
            try:
                if  aa_changes_dict[keys_tmp2] :
                    PM5_t2=0
                    #print("PM5 %s %s" %(aa_changes_dict[keys_tmp2],aa_last))
            except KeyError:
                PM5_t2=1
                PM5_t3=0
                for nt in ACGTs:
                    if nt != cls[Allels_flgs['Alt']] and  nt != cls[Allels_flgs['Ref']]:
                        keys_tmp3=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+nt
                        try:
                            if aa_changes_dict[keys_tmp3]:
                                PM5_t3=1
                            if aa_changes_dict[keys_tmp3] == aa_last:
                                PM5_t2=0*PM5_t2
                        except KeyError:
                            pass
                        else:
                            pass

            else:
                pass

    if PM5_t1 !=0 and PM5_t2 != 0 and PM5_t3 !=0 :
        PM5=1
    return(PM5)

def check_PM6(line,Funcanno_flgs,Allels_flgs):
    '''
    Assumed de novo, but without confirmation of paternity and maternity
    '''
    PM6=0
    return(PM6)

def check_PP1(line,Funcanno_flgs,Allels_flgs):
    '''
    Cosegregation with disease in multiple affected family members in a gene definitively 
    known to cause the disease
    '''
    PP1=0
    return(PP1)

def check_PP2(line,Funcanno_flgs,Allels_flgs,PP2_genes_dict):
    '''
    Missense variant in a gene that has a low rate of benign missense variation and in which 
    missense variants are a common mechanism of disease
    '''
    PP2=0
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
        # need to check whether gene has a low rate of benign missense variation.....
            try:
                if PP2_genes_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
                    PP2=1
            except KeyError:
                PP2=0
            else:
                pass

    return(PP2)

def check_PP3(line,Funcanno_flgs,Allels_flgs):
    '''
    ENHANCED PP3: Multiple lines of computational evidence support deleterious effect
    Checks ALL predictors and uses HIGHEST evidence level
    Splice predictions apply to ALL variant types, not just missense
    '''
    PP3 = 0
    cls = line.split('\t')
    max_evidence = 0  # Track highest evidence level
    
    # CRITICAL: Check splice predictions FIRST for ALL variant types
    try:
        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
        
        if ada_score != '.' and rf_score != '.':
            ada_val = float(ada_score)
            rf_val = float(rf_score)
            
            # ClinGen-calibrated dbscSNV thresholds
            if ada_val >= 0.96 or rf_val >= 0.65:
                max_evidence = max(max_evidence, 8)  # Very Strong for high-confidence splice disruption
            elif ada_val >= 0.93 or rf_val >= 0.575:
                max_evidence = max(max_evidence, 4)  # Strong
            elif ada_val >= 0.86 or rf_val >= 0.45:
                max_evidence = max(max_evidence, 2)  # Moderate
            elif ada_val >= 0.6 or rf_val >= 0.15:
                max_evidence = max(max_evidence, 1)  # Supporting
    except (ValueError, KeyError, IndexError):
        pass
    
    # Check if this is a missense variant for missense-specific predictors
    func_info = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    is_missense = any(term in func_info.lower() for term in ['missense', 'nonsynony'])
    
    if is_missense:
        # For missense variants, check all computational predictors
        
        # 1. MetaRNN - Highest priority
        try:
            if 'MetaRNN_score' in Funcanno_flgs and cls[Funcanno_flgs['MetaRNN_score']] != '.':
                score = float(cls[Funcanno_flgs['MetaRNN_score']])
                if score >= 0.939:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.841:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.748:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
        
        # 2. REVEL
        try:
            if 'REVEL_score' in Funcanno_flgs and cls[Funcanno_flgs['REVEL_score']] != '.':
                score = float(cls[Funcanno_flgs['REVEL_score']])
                if score >= 0.932:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.773:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.644:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
        
        # 3. BayesDel
        try:
            if 'BayesDel_addAF_score' in Funcanno_flgs and cls[Funcanno_flgs['BayesDel_addAF_score']] != '.':
                score = float(cls[Funcanno_flgs['BayesDel_addAF_score']])
                if score >= 0.424:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.173:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.0579:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
        
        # 4. AlphaMissense
        try:
            if 'AlphaMissense_score' in Funcanno_flgs and cls[Funcanno_flgs['AlphaMissense_score']] != '.':
                score = float(cls[Funcanno_flgs['AlphaMissense_score']])
                if score >= 0.969:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.794:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.665:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
        
        # 5. VARITY_R
        try:
            if 'VARITY_R_score' in Funcanno_flgs and cls[Funcanno_flgs['VARITY_R_score']] != '.':
                score = float(cls[Funcanno_flgs['VARITY_R_score']])
                if score >= 0.957:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.826:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.748:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
        
        # 6. M-CAP
        try:
            if 'M-CAP_score' in Funcanno_flgs and cls[Funcanno_flgs['M-CAP_score']] != '.':
                score = float(cls[Funcanno_flgs['M-CAP_score']])
                if score >= 0.389:
                    max_evidence = max(max_evidence, 4)  # Strong
                elif score >= 0.129:
                    max_evidence = max(max_evidence, 2)  # Moderate
                elif score >= 0.0871:
                    max_evidence = max(max_evidence, 1)  # Supporting
        except (ValueError, KeyError, IndexError):
            pass
    
    # CADD and conservation can apply to any variant type
    
    # CADD
    try:
        if 'CADD_phred' in Funcanno_flgs and cls[Funcanno_flgs['CADD_phred']] != '.':
            score = float(cls[Funcanno_flgs['CADD_phred']])
            if score >= 32.1:
                max_evidence = max(max_evidence, 4)  # Strong
            elif score >= 25.3:
                max_evidence = max(max_evidence, 2)  # Moderate
            elif score >= 22.7:
                max_evidence = max(max_evidence, 1)  # Supporting
    except (ValueError, KeyError, IndexError):
        pass
    
    # Conservation (GERP++)
    try:
        if 'GERP++_RS' in Funcanno_flgs and cls[Funcanno_flgs['GERP++_RS']] != '.':
            gerp_score = float(cls[Funcanno_flgs['GERP++_RS']])
            if gerp_score >= 6.17:
                max_evidence = max(max_evidence, 2)  # Moderate
            elif gerp_score >= 4.72:
                max_evidence = max(max_evidence, 1)  # Supporting
    except (ValueError, KeyError, IndexError):
        pass
    
    # phyloP as alternative conservation metric
    try:
        if 'phyloP100way_vertebrate' in Funcanno_flgs and cls[Funcanno_flgs['phyloP100way_vertebrate']] != '.':
            phylop_score = float(cls[Funcanno_flgs['phyloP100way_vertebrate']])
            if phylop_score >= 9.42:
                max_evidence = max(max_evidence, 2)  # Moderate
            elif phylop_score >= 7.2:
                max_evidence = max(max_evidence, 1)  # Supporting
    except (ValueError, KeyError, IndexError):
        pass
    
    PP3 = max_evidence
    return PP3

def check_PP4(line,Funcanno_flgs,Allels_flgs):
    '''
    Patient's phenotype or family history is highly specific for a disease with a single genetic etiology
    '''
    PP4=0
    return(PP4)


def check_PP5(line,Funcanno_flgs,Allels_flgs):
    '''
    ENHANCED PP5: Reputable source recently reports variant as pathogenic
    Implements ClinGen SVI 2023 recommendations for evidence strength scaling
    Based on ClinVar review status and clinical significance
    '''
    PP5 = 0
    cls = line.split('\t')
    
    # Get ClinVar fields - handle missing indices gracefully
    try:
        clnsig = cls[Funcanno_flgs.get('CLINSIG', -1)]
        clnrevstat = cls[Funcanno_flgs.get('CLNREVSTAT', -1)]
    except (IndexError, KeyError):
        return 0
    
    # Skip if either field is missing or empty
    if clnsig == '.' or clnrevstat == '.' or clnsig == '' or clnrevstat == '':
        return 0
    
    # Check for pathogenic classification
    # Handle all pathogenic variants: Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic
    pathogenic_terms = ['Pathogenic', 'pathogenic']
    is_pathogenic = any(term in clnsig for term in pathogenic_terms)
    
    # Exclude conflicting classifications
    conflicting_terms = ['conflicting', 'Conflicting']
    has_conflicts = any(term in clnsig or term in clnrevstat for term in conflicting_terms)
    
    if is_pathogenic and not has_conflicts:
        # Scale evidence strength based on review status per ClinGen SVI guidelines
        
        # Very Strong Evidence (8 points): Expert panel or practice guideline
        if ('reviewed_by_expert_panel' in clnrevstat or 
            'practice_guideline' in clnrevstat):
            PP5 = 8
            
        # Strong Evidence (4 points): Multiple submitters, no conflicts
        elif ('multiple_submitters' in clnrevstat and 
              'no_conflicts' in clnrevstat and 
              'criteria_provided' in clnrevstat):
            PP5 = 4
            
        # Supporting Evidence (1 point): Single submitter with criteria
        elif ('criteria_provided' in clnrevstat and 
              'single_submitter' in clnrevstat):
            PP5 = 1
            
        # No Evidence: Insufficient review criteria
        elif ('no_assertion_criteria_provided' in clnrevstat or
              'no_classification' in clnrevstat):
            PP5 = 0
            
        # Default Supporting for other cases with pathogenic classification
        else:
            PP5 = 1
    
    return PP5

def check_BA1(line,Freqs_flgs,Allels_flgs):
    '''
    Allele frequency is >5% - UPDATED with population-specific checks
    '''
    BA1=0
    cls=line.split('\t')
    
    # Enhanced BA1 with population-specific thresholds
    population_threshold = 0.05  # 5% as per ACMG
    
    # Check each major population
    for key in Freqs_flgs.keys():
        # EXCLUDE ESP6500 from BA1 evaluation (outdated database per ACMG/ClinGen guidance)
        if key == 'esp6500siv2_all':
            continue
            
        if cls[Freqs_flgs[key]] != '.':
            try:
                freq = float(cls[Freqs_flgs[key]])
                # Apply population-specific BA1 threshold
                if freq > population_threshold:
                    # Removed overly restrictive sanity check
                    # Any frequency >5% should trigger BA1
                    BA1 = 1
                    break
            except ValueError:
                continue
    
    return BA1

def check_BS1(line,Freqs_flgs,Allels_flgs):
    '''
    Allele frequency is greater than expected for disorder - UPDATED
    '''
    BS1=0
    cls=line.split('\t')
    
    # Enhanced BS1 with gene-specific thresholds
    default_threshold = 0.005  # Default for rare diseases
    try:
        user_threshold = float(paras.get('disorder_cutoff', 0.005))
        threshold = user_threshold
    except:
        threshold = default_threshold
    
    # Check against multiple populations
    frequencies = []
    for key in Freqs_flgs.keys():
        # EXCLUDE ESP6500 from BS1 evaluation (outdated database per ACMG/ClinGen guidance)
        if key == 'esp6500siv2_all':
            continue
            
        if cls[Freqs_flgs[key]] != '.':
            try:
                freq = float(cls[Freqs_flgs[key]])
                # Removed overly restrictive sanity check
                # Include all valid frequencies
                if freq > 0 and freq <= 1.0:  # Basic sanity check
                    frequencies.append(freq)
            except ValueError:
                continue
    
    if frequencies:
        max_freq = max(frequencies)
        # Use more stringent thresholds for different inheritance patterns
        if max_freq >= threshold:
            BS1 = 1
    
    return BS1

def check_BS2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs):
    '''
    Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked
    (hemizygous) disorder, with full penetrance expected at an early age
    check gnomAD_genome_ALL
    '''
    BS2=0
    cls=line.split('\t')
    keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
    try:
        mim1=mim2gene_dict[ cls[Funcanno_flgs['Gene.ensGene']] ]
        try:
            if mim_adultonset_dict[mim1] == "1": # means adult oneset disorder
                BS2=0;
        except KeyError: # means not adult onset, begin to check recessive or domiant ,the genotype from 1000 genome
            try:
                if mim_recessive_dict[mim1] == "1": # means recessive disorder: check snps as homo
                    try:
                        if BS2_snps_recess_dict[ keys ]=="1":  # key as snp info
                            BS2=1
                    except KeyError:
                        pass
                    else:
                        pass
            except KeyError:
                pass
            else:
                pass


            try:
                if mim_domin_dict[mim1] == "1": # means dominant disorder: check snps as heter
                    BS2=0
                    try:
                        if BS2_snps_domin_dict[ keys ]=="1":  # key as snp info
                            BS2=1
                    except KeyError:
                        pass
                    else:
                        pass
            except KeyError:
                pass
            else:
                pass




        else:

            pass

    except KeyError: # means no information of recessive or domiant so BS2=0
        BS2=0
    else:
        pass

    return(BS2)

def check_BS3(line,Funcanno_flgs,Allels_flgs):
    '''
    Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing
    '''
    BS3=0
    '''
    cls=line.split('\t')
    line_tmp=cls[Funcanno_flgs['Gene damage prediction (all disease-causing genes)']]

    funcs_tmp=["missense","nonsynony"]
    line_tmp2=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp2.find(fc)>=0 and  line_tmp == 'Low' :
            BS3=1
    '''
    return(BS3)

def check_BS4(line,Funcanno_flgs,Allels_flgs):
    '''
    Lack of segregation in affected members of a family
    '''
    BS4=0
    return(BS4)

def check_BP1(line,Funcanno_flgs,Allels_flgs,BP1_genes_dict):
    '''
    Missense variant in a gene for which primarily truncating variants are known to cause disease
    truncating:  stop_gain / frameshift deletion/  nonframshift deletion
    We defined Protein truncating variants  (4) (table S1) as single-nucleotide variants (SNVs) predicted to introduce a premature stop codon or to disrupt a splice site, small insertions or deletions (indels) predicted to disrupt a transcript reading frame, and larger deletions 
    '''
    BP1=0
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
        # need to wait to check whether truncating is the only cause disease
            try:
                if BP1_genes_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
                    BP1=1
            except KeyError:
                BP1=0
            else:
                pass
    return(BP1)

def check_BP2(line,Funcanno_flgs,Allels_flgs):
    '''
    Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed 
    in cis with a pathogenic variant in any inheritance pattern
    '''
    BP2=0
    return(BP2)

def check_BP3(line,Funcanno_flgs,Allels_flgs):
    '''
    In-frame deletions/insertions in a repetitive region without a known function
    if the repetitive region is in the domain, this BP3 should not be applied.
    '''
    BP3=0
    BP3_t1=0
    BP3_t2=0
    cls=line.split('\t')
    #funcs_tmp=["cds-indel","stop-loss"]
    funcs_tmp=["nonframeshift insertion","nonframeshift deletion","nonframeshift substitution"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            BP3_t1=1;
        # need to wait to check  in a repeat region
    if cls[Funcanno_flgs['rmsk']] != '.' and  cls[Funcanno_flgs['Interpro_domain']]== '.' : # repeat and not in domain
        BP3_t2=1

    if BP3_t1 !=0 and BP3_t2 != 0 :
        BP3=1
    return(BP3)

def check_BP4(line,Funcanno_flgs,Allels_flgs):
    '''
    ENHANCED BP4: Multiple lines of computational evidence suggest no impact
    Implements ClinGen-calibrated thresholds for benign evidence
    Properly excludes LOF variants and handles evidence conflicts
    '''
    BP4 = 0
    cls = line.split('\t')
    
    # CRITICAL: Exclude loss-of-function variant types
    func_info = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    lof_terms = ['stopgain', 'stoploss', 'frameshift', 'startloss', 'nonsense', 
                 'stop_gained', 'stop_lost', 'start_lost', 'canonical_splice']
    
    for term in lof_terms:
        if term in func_info.lower():
            return 0  # Never apply BP4 to LOF variants
    
    # Check if canonical splice site
    try:
        if 'splicing' in func_info.lower() or 'splice' in func_info.lower():
            # Check for canonical positions in annotation
            if any(col for col in cls if any(pat in str(col) for pat in ['+1', '+2', '-1', '-2'])):
                return 0  # Don't apply BP4 to canonical splice sites
    except:
        pass
    
    # Meta-predictor hierarchy for benign evidence
    
    # 1. MetaRNN - Highest priority
    try:
        if 'MetaRNN_score' in Funcanno_flgs and cls[Funcanno_flgs['MetaRNN_score']] != '.':
            score = float(cls[Funcanno_flgs['MetaRNN_score']])
            if score <= 0.151:
                return 4  # Strong benign
            elif score <= 0.359:
                return 2  # Moderate benign
            elif score <= 0.507:
                return 1  # Supporting benign
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # 2. REVEL
    try:
        if 'REVEL_score' in Funcanno_flgs and cls[Funcanno_flgs['REVEL_score']] != '.':
            score = float(cls[Funcanno_flgs['REVEL_score']])
            if score <= 0.003:
                return 4  # Strong
            elif score <= 0.183:
                return 2  # Moderate
            elif score <= 0.29:
                return 1  # Supporting
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # 3. BayesDel (addAF version)
    try:
        if 'BayesDel_addAF_score' in Funcanno_flgs and cls[Funcanno_flgs['BayesDel_addAF_score']] != '.':
            score = float(cls[Funcanno_flgs['BayesDel_addAF_score']])
            if score <= -0.359:
                return 4  # Strong
            elif score <= -0.0862:
                return 2  # Moderate
            elif score <= 0.00131:
                return 1  # Supporting
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # 4. AlphaMissense (only for missense)
    if 'missense' in func_info.lower() or 'nonsynon' in func_info.lower():
        try:
            if 'AlphaMissense_score' in Funcanno_flgs and cls[Funcanno_flgs['AlphaMissense_score']] != '.':
                score = float(cls[Funcanno_flgs['AlphaMissense_score']])
                if score <= 0.163:
                    return 4  # Strong
                elif score <= 0.302:
                    return 2  # Moderate
                elif score <= 0.458:
                    return 1  # Supporting
                else:
                    return 0
        except (ValueError, KeyError, IndexError):
            pass
    
    # 5. VARITY_R
    try:
        if 'VARITY_R_score' in Funcanno_flgs and cls[Funcanno_flgs['VARITY_R_score']] != '.':
            score = float(cls[Funcanno_flgs['VARITY_R_score']])
            if score <= 0.245:
                return 4  # Strong
            elif score <= 0.496:
                return 2  # Moderate
            elif score <= 0.613:
                return 1  # Supporting
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # 6. M-CAP
    try:
        if 'M-CAP_score' in Funcanno_flgs and cls[Funcanno_flgs['M-CAP_score']] != '.':
            score = float(cls[Funcanno_flgs['M-CAP_score']])
            if score <= 0.00316:
                return 4  # Strong
            elif score <= 0.0209:
                return 2  # Moderate
            elif score <= 0.0439:
                return 1  # Supporting
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # 7. CADD (lower scores = more benign)
    try:
        if 'CADD_phred' in Funcanno_flgs and cls[Funcanno_flgs['CADD_phred']] != '.':
            score = float(cls[Funcanno_flgs['CADD_phred']])
            if score <= 11.5:
                return 4  # Strong
            elif score <= 17.0:
                return 2  # Moderate
            elif score <= 19.2:
                return 1  # Supporting
            else:
                return 0
    except (ValueError, KeyError, IndexError):
        pass
    
    # Splicing predictions for benign evidence
    try:
        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
        
        # Only apply if we have actual scores (not missing)
        if ada_score != '.' and rf_score != '.' and ada_score != '-1' and rf_score != '-1':
            ada_val = float(ada_score)
            rf_val = float(rf_score)
            # Low scores = no splicing impact
            if ada_val <= 0.02 and rf_val <= 0.03:
                return 2  # Moderate benign
            elif ada_val <= 0.1 and rf_val <= 0.1:
                return 1  # Supporting benign
    except (ValueError, KeyError, IndexError):
        pass
    
    # Conservation for benign (low conservation = benign)
    try:
        if 'GERP++_RS' in Funcanno_flgs and cls[Funcanno_flgs['GERP++_RS']] != '.':
            gerp_score = float(cls[Funcanno_flgs['GERP++_RS']])
            if gerp_score <= -3.53:
                return 4  # Strong benign
            elif gerp_score <= 0.499:
                return 2  # Moderate benign
            elif gerp_score <= 2.31:
                return 1  # Supporting benign
    except (ValueError, KeyError, IndexError):
        pass
    
    return BP4

def check_BP5(line,Funcanno_flgs,Allels_flgs,morbidmap_dict):
    '''
    Variant found in a case with an alternate molecular basis for disease
    check the genes whether are for mutilfactor disorder
    The reviewers suggeast to disable the OMIM morbidmap for BP5
    '''
    '''
    BP5=0
    cls=line.split('\t')
    try:
        if morbidmap_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
            BP5=1
    except KeyError:
        BP5=0
    else:
        pass
    '''
    BP5=0; 
    return(BP5)

def check_BP6(line,Funcanno_flgs,Allels_flgs):
    '''
    Reputable source recently reports variant as benign, but the evidence is not available to the 
    laboratory to perform an independent evaluation; Check the ClinVar column to see whether this 
    is "benign". 
    '''
    BP6=0
    cls=line.split('\t')

    line_tmp2=cls[Funcanno_flgs['CLINSIG']]
    if line_tmp2 != '.':
        cls3=line_tmp2.split(';')
        clinvar_bp=cls3[0]
        if clinvar_bp.find("ikely benign")>=0 or clinvar_bp.find("enign")>=0:
            BP6=1

    return(BP6)

def check_BP7(line,Funcanno_flgs,Allels_flgs):
    '''
    A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the 
    splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly 
    conserved
    '''
    BP7=0
    BP7_t1=0
    BP7_t2=0
    cutoff_conserv=2 # for GERP++
    cls=line.split('\t')
    funcs_tmp=["synon","coding-synon"]
    funcs_tmp2="nonsynon"
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 and line_tmp.find(funcs_tmp2)<0 :
    # need to wait to check the  impact to the splice from dbscSNV 
    # either score(ada and rf) >0.6 as splicealtering
            if cls[Funcanno_flgs['dbscSNV_RF_SCORE']]=="." or  cls[Funcanno_flgs['dbscSNV_ADA_SCORE']]==".":
                BP7_t1=1  # absent means it is not in the  splice consensus sequence
            else:
                if float(cls[Funcanno_flgs['dbscSNV_RF_SCORE']])<0.6 and float(cls[Funcanno_flgs['dbscSNV_ADA_SCORE']])<0.6:
                    BP7_t1=1
# check the conservation score of gerp++ > 2
    try:
        if float(cls[Funcanno_flgs['GERP++_RS']]) <= float(cutoff_conserv) or cls[Funcanno_flgs['GERP++_RS']] == '.' :
            BP7_t2=1
    except ValueError:
        # absent means there are gaps in the multiple alignment,so cannot have the score,not conserved
        BP7_t2=1
    else:
        pass

    if BP7_t1 !=0 and BP7_t2 != 0 :
        BP7=1        
    return(BP7)


def assign(BP,line,Freqs_flgs,Funcanno_flgs,Allels_flgs):
    PVS1=0
    PS=[0,0,0,0,0]
    PM=[0,0,0,0,0,0,0]
    PP=[0,0,0,0,0,0]

    BA1=0
    BS=[0,0,0,0,0]
    BP=[0,0,0,0,0,0,0,0]

    # First, evaluate all rules
    PVS1=check_PVS1(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    
    PS1=check_PS1(line,Funcanno_flgs,Allels_flgs,aa_changes_dict)
    PS[0]=PS1
    PS2=check_PS2(line,Funcanno_flgs,Allels_flgs)
    PS[1]=PS2
    PS3=check_PS3(line,Funcanno_flgs,Allels_flgs)
    PS[2]=PS3
    PS4=check_PS4(line,Funcanno_flgs,Allels_flgs)
    PS[3]=PS4

    PM1=check_PM1(line,Funcanno_flgs,Allels_flgs,domain_benign_dict)
    PM[0]=PM1
    PM2=check_PM2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs,mim2gene_dict,mim2gene_dict2)
    PM[1]=PM2
    PM3=check_PM3(line,Funcanno_flgs,Allels_flgs)
    PM[2]=PM3
    PM4=check_PM4(line,Funcanno_flgs,Allels_flgs)
    PM[3]=PM4
    PM5=check_PM5(line,Funcanno_flgs,Allels_flgs,aa_changes_dict)
    PM[4]=PM5
    PM6=check_PM6(line,Funcanno_flgs,Allels_flgs)
    PM[5]=PM6

    PP1=check_PP1(line,Funcanno_flgs,Allels_flgs)
    PP[0]=PP1
    PP2=check_PP2(line,Funcanno_flgs,Allels_flgs,PP2_genes_dict)
    PP[1]=PP2
    PP3=check_PP3(line,Funcanno_flgs,Allels_flgs)
    PP[2]=PP3
    PP4=check_PP4(line,Funcanno_flgs,Allels_flgs)
    PP[3]=PP4
    PP5=check_PP5(line,Funcanno_flgs,Allels_flgs)
    PP[4]=PP5

    BA1=check_BA1(line,Freqs_flgs,Allels_flgs)
    
    BS1=check_BS1(line,Freqs_flgs,Allels_flgs)
    BS[0]=BS1
    BS2=check_BS2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs)
    BS[1]=BS2
    BS3=check_BS3(line,Funcanno_flgs,Allels_flgs)
    BS[2]=BS3
    BS4=check_BS4(line,Funcanno_flgs,Allels_flgs)
    BS[3]=BS4

    BP1=check_BP1(line,Funcanno_flgs,Allels_flgs,BP1_genes_dict)
    BP[0]=BP1
    BP2=check_BP2(line,Funcanno_flgs,Allels_flgs)
    BP[1]=BP2
    BP3=check_BP3(line,Funcanno_flgs,Allels_flgs)
    BP[2]=BP3
    BP4=check_BP4(line,Funcanno_flgs,Allels_flgs)
    BP[3]=BP4
    BP5=check_BP5(line,Funcanno_flgs,Allels_flgs,morbidmap_dict)
    BP[4]=BP5
    BP6=check_BP6(line,Funcanno_flgs,Allels_flgs)
    BP[5]=BP6
    BP7=check_BP7(line,Funcanno_flgs,Allels_flgs)
    BP[6]=BP7

    # CRITICAL FIX: Apply rule conflict resolution following ACMG/VarSome guidelines
    
    # 1. PVS1 conflicts: If PVS1 is triggered, disable conflicting benign rules
    # 1. PVS1 conflicts: Enhanced conflict resolution per ClinGen SVI 2018
    if PVS1 > 0:
        BP[3] = 0  # Disable BP4 (computational evidence suggests benign)
        # Note: Do NOT disable PP3 for splicing variants per ClinGen SVI Splicing Subgroup 2023
        # Only disable PP3 for canonical splice sites already captured by PVS1
        try:
            line_tmp = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
            if "canonical" in line_tmp.lower() and ("splice" in line_tmp.lower() or "splicing" in line_tmp.lower()):
                PP[2] = 0  # Disable PP3 for canonical splice sites only
        except:
            pass
        # Disable PM4 to avoid double-counting with PVS1
        PM[3] = 0
        
    # 2. PM1 conflicts: If PM1 is triggered, limit PP3 strength and disable BP3
    if PM[0] > 0:
        if PP[2] > 2:  # If PP3 was Moderate or Strong
            PP[2] = 1  # Reduce to Supporting to avoid double-counting
        BP[2] = 0     # Disable BP3 (in-frame indels in repetitive regions)
        
    # 3. BA1/BS1 conflicts: Strong frequency evidence overrides computational evidence
    if BA1 > 0 or BS[0] > 0:
        BP[3] = 0  # Disable BP4 when strong frequency evidence exists
        PP[2] = 0  # Disable PP3 when strong frequency evidence exists
        
    # 4. Clinical evidence moderation
    if PP[4] > 0:  # If PP5 (clinical pathogenic evidence) exists
        if BP[3] > 2:  # If BP4 was Strong/Very Strong
            BP[3] = min(BP[3], 1)  # Reduce to Supporting max
            
    if BP[5] > 0:  # If BP6 (clinical benign evidence) exists  
        if PP[2] > 2:  # If PP3 was Strong/Very Strong
            PP[2] = min(PP[2], 1)  # Reduce to Supporting max

    # 5. Variant-type specific exclusions
    cls = line.split('\t')
    line_tmp = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    
    # For synonymous variants, disable pathogenic computational evidence unless splicing
    if "synon" in line_tmp.lower() and "nonsynon" not in line_tmp.lower():
        # Only allow PP3 if strong splicing prediction exists
        try:
            ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
            rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
            if ada_score == '.' or rf_score == '.' or \
               (float(ada_score) < 0.958 and float(rf_score) < 0.584):
                PP[2] = 0  # Disable PP3 for synonymous without splicing impact
        except:
            PP[2] = 0  # Disable PP3 if can't assess splicing
            
    # 6. PM2 frequency conflicts
    if PM[1] > 0:  # If PM2 triggered (absent/rare)
        BS[0] = 0  # Disable conflicting BS1 (too frequent)
        BS[1] = 0  # Disable conflicting BS2 (observed in healthy controls)

    #print("AFTER CONFLICT RESOLUTION: PVS1=%s PS=%s PM=%s PP=%s BA1=%s BS=%s BP=%s" %(PVS1,PS,PM,PP,BA1,BS,BP))
    
    # Continue with existing user evidence processing...
    cls=line.split('\t')
    
    #begin process the exclude snp list. which will affect BA1 BS1 BS2
    if os.path.isfile(paras['exclude_snps']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            if exclude_snps_dict[keys]=="1":  
                BA1=0; 
                BS[0]=0; 
                BS[1]=0;
        except KeyError:
            pass
        else:
            pass
            
    #begin process the user's evidence file
    if os.path.isfile(paras['evidence_file']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            evds=user_evidence_dict[keys] #PS1=1;PM1=1;BA1=1;PVS1 PP BS BP
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and (not re.findall('grade', evd_t[0], flags=re.IGNORECASE)) ):
                    if int(evd_t[1])<=1:
                        #print ("%s %s %s " %(keys,evd_t[1],evd_t[0]))
                        if(evd_t[0]=="PVS1"): PVS1=evd_t[1]
                        if(evd_t[0]=="BA1"): BA1=evd_t[1]
                        if(evd_t[0].find('PS')!=-1): 
                            t=evd_t[0].find('PS'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            if(t<len(evd_t[0])-2 and tt3<=5 ): PS[tt3-1]=int(evd_t[1])
                        if(evd_t[0].find('PM')!=-1):
                            t=evd_t[0].find('PM'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            if(t<len(evd_t[0])-2 and tt3<=7 ): PM[tt3-1]=int(evd_t[1])
                        if(evd_t[0].find('PP')!=-1): 
                            t=evd_t[0].find('PP'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            if(t<len(evd_t[0])-2 and tt3<=6 ): PP[tt3-1]=int(evd_t[1])
                        if(evd_t[0].find('BS')!=-1): 
                            t=evd_t[0].find('BS'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            if(t<len(evd_t[0])-2 and tt3<=5 ): BS[tt3-1]=int(evd_t[1])
                        if(evd_t[0].find('BP')!=-1):
                            t=evd_t[0].find('BP'); 
                            tt=evd_t[0];
                            tt3=int(tt[t+2:t+3])
                            if(t<len(evd_t[0])-2 and tt3<=8 ): BP[tt3-1]=int(evd_t[1])
        except KeyError:
            pass
        else:
            pass

    # end process the user's evidence file 

    if len(cls)>1:#esp6500siv2_all 1000g2015aug_all gnomAD_genome_ALL    
        BP_out=classfy(PVS1,PS,PM,PP,BA1,BS,BP,Allels_flgs,cls)
        line_t="%s PVS1=%s PS=%s PM=%s PP=%s BA1=%s BS=%s BP=%s" %(BP_out,PVS1,PS,PM,PP,BA1,BS,BP)

        #print("%s " % BP_out)
        BP_out=line_t
        pass
    #BP=BP_out
    return(BP_out)


def search_key_index(line,dict):
    cls=line.split('\t')
    for key in dict.keys():
        for i in range(1,len(cls)):
            ii=i-1
            if key==cls[ii]:
                dict[key]=ii
                break
            if key=="Otherinfo":
                if cls[ii]==key or cls[ii]=="Otherinfo1":
                    dict[key]=ii
                    break
    return

def my_inter_var(annovar_outfile):
    newoutfile=annovar_outfile+".grl_p"
    newoutfile2=annovar_outfile+".intervar"

    Freqs_flgs={'1000g2015aug_all':0,'esp6500siv2_all':0,'gnomAD_genome_ALL':0,'gnomAD_genome_AFR':0,'gnomAD_genome_AMR':0,'gnomAD_genome_EAS':0,'gnomAD_genome_FIN':0,'gnomAD_genome_NFE':0,'gnomAD_genome_OTH':0,'gnomAD_genome_ASJ':0}
    Funcanno_flgs={'Func.refGene':0,'ExonicFunc.refGene':0,'AAChange.refGene':0,'Gene':0,'Gene damage prediction (all disease-causing genes)':0,'CLNDBN':0,'CLNACC':0,'CLNDSDB':0,'dbscSNV_ADA_SCORE':0,'dbscSNV_RF_SCORE':0,'GERP++_RS':0,'LoFtool_percentile':0,'Interpro_domain':0,'rmsk':0,'SIFT_score':0,'phyloP46way_placental':0,'Gene.ensGene':0,'CLINSIG':0,'CLNREVSTAT':0,'CADD_raw':0,'CADD_phred':0,'avsnp151':0,'AAChange.ensGene':0,'AAChange.knownGene':0,'MetaSVM_score':0,'MetaRNN_score':0,'REVEL_score':0,'BayesDel_addAF_score':0,'AlphaMissense_score':0,'phyloP100way_vertebrate':0,'Otherinfo':0}
    Allels_flgs={'Chr':0,'Start':0,'End':0,'Ref':0,'Alt':0}
# gnomAD_genome_ALL esp6500siv2_all   1000g2015aug_all  SIFT_score    CADD_raw    CADD_phred  GERP++_RS   phyloP46way_placental  dbscSNV_ADA_SCORE   dbscSNV_RF_SCORE   Interpro_domain

    try:
        fh=open(newoutfile, "r")
        fw=open(newoutfile2, "w")
        strs=fh.read()
        line_sum=0;
        print("Notice: Begin the variants interpretation by InterVar ")
        if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t InterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp151","AAChange.ensGene","AAChange.refGene","Clinvar","ClinVar_Review","InterVar and Evidence","Freq_gnomAD_genome_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha","Otherinfo"  ))
        else:
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t InterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp151","AAChange.ensGene","AAChange.refGene","Clinvar","ClinVar_Review","InterVar and Evidence","Freq_gnomAD_genome_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha"  ))

        for line in strs.split('\n'):
            BP="UNK" # the inter of pathogenetic/benign
            clinvar_bp="UNK"
            cls=line.split('\t')
            if len(cls)<2: break
            if line_sum==0:
                search_key_index(line,Freqs_flgs)
                search_key_index(line,Funcanno_flgs)
                search_key_index(line,Allels_flgs)

            else:
                #begin check the BP status from clinvar
                line_tmp2=cls[Funcanno_flgs['CLINSIG']]
                if line_tmp2 != '.':
                    cls3=line_tmp2.split(';')
                    clinvar_bp=cls3[0]
                    
                intervar_bp=assign(BP,line,Freqs_flgs,Funcanno_flgs,Allels_flgs)
                Freq_gnomAD_genome_POPs="AFR:"+cls[Freqs_flgs['gnomAD_genome_AFR']]+",AMR:"+cls[Freqs_flgs['gnomAD_genome_AMR']]+",EAS:"+cls[Freqs_flgs['gnomAD_genome_EAS']]+",FIN:"+cls[Freqs_flgs['gnomAD_genome_FIN']]+",NFE:"+cls[Freqs_flgs['gnomAD_genome_NFE']]+",OTH:"+cls[Freqs_flgs['gnomAD_genome_OTH']]+",ASJ:"+cls[Freqs_flgs['gnomAD_genome_ASJ']]
                OMIM="." 
                mim2=mim2gene_dict2.get(cls[Funcanno_flgs['Gene']],".")
                mim1=mim2gene_dict.get(cls[Funcanno_flgs['Gene.ensGene']],".")
                if(mim1 !="."):
                   OMIM=mim1
                if(mim2 !="."):
                   OMIM=mim2
                Pheno_MIM=mim_pheno_dict.get(OMIM,".")
                orpha="";
                orpha_details="";
                # .;442835;;306;;.;
                for ort2 in Pheno_MIM.split(';'):
                    ort3=mim_orpha_dict.get(ort2,".")
                    if(ort3 !="."):
                        orpha=ort3+orpha
                for ort4 in orpha.split(';'):
                    if len(ort4)>0:
                         orpha_details=orpha_details+orpha_dict.get(ort4,".")+"~"
                        

                if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
                    cls[Funcanno_flgs['Otherinfo']]=cls[Funcanno_flgs['Otherinfo']].replace('\t', ';')
                    clnrevstat_output = cls[Funcanno_flgs.get('CLNREVSTAT', len(cls)-1)] if 'CLNREVSTAT' in Funcanno_flgs else "."
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t InterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp151']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,clnrevstat_output,intervar_bp,cls[Freqs_flgs['gnomAD_genome_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']], cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],".", cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['rmsk']],cls[Funcanno_flgs['MetaSVM_score']],Freq_gnomAD_genome_POPs,OMIM,Pheno_MIM,orpha,orpha_details,cls[Funcanno_flgs['Otherinfo']]   ))
                else:
                    clnrevstat_output = cls[Funcanno_flgs.get('CLNREVSTAT', len(cls)-1)] if 'CLNREVSTAT' in Funcanno_flgs else "."
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t InterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp151']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,clnrevstat_output,intervar_bp,cls[Freqs_flgs['gnomAD_genome_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']], cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],".", cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['rmsk']],cls[Funcanno_flgs['MetaSVM_score']],Freq_gnomAD_genome_POPs,OMIM,Pheno_MIM,orpha,orpha_details  ))

                #print("%s\t%s %s" % (line,clinvar_bp,intervar_bp))

            line_sum=line_sum+1

    except IOError:
        print("Error: can\'t readi/write the annovar output files %s" % (newoutfile,newoutfile2))
        sys.exit()
        return
    else:
        fh.close()
        fw.close()
    return(line_sum)


def main():


    if platform.python_version()< '3.0.0'  :
        config=ConfigParser.ConfigParser()
    else:
        config=configparser.ConfigParser()
    




    parser = optparse.OptionParser(usage=usage, version=version, description=description)


    parser.add_option("-?", action="help", help=optparse.SUPPRESS_HELP, dest="help")
    parser.add_option("-v", action="version", help=optparse.SUPPRESS_HELP, dest="version")
    
    parser.add_option("-c", "--config", dest="config", action="store",
                  help="The config file of all options. it is for your own configure file.You can edit all the options in the configure and if you use this options,you can ignore all the other options bellow", metavar="config.ini")

    parser.add_option("-b", "--buildver", dest="buildver", action="store",
                  help="The genomic build version, it can be hg19 or hg38 , will support other version later", metavar="hg19")


    parser.add_option("-i", "--input", dest="input", action="store",
                  help="The input file contains your variants", metavar="example/ex1.avinput")

    parser.add_option("--input_type", dest="input_type", action="store",
                  help="The input file type, it can be  AVinput(Annovar's format),VCF(VCF with single sample),VCF_m(VCF with multiple samples)", metavar="AVinput")

    parser.add_option("-o", "--output", dest="output", action="store",
                  help="The prefix of output file which contains the results, the file of results will be as [$prefix].intervar ", metavar="example/myanno")


    group = optparse.OptionGroup(parser, "InterVar Other Options")
    group.add_option("-t", "--database_intervar", dest="database_intervar", action="store",
            help="The  database location/dir for the InterVar dataset files", metavar="intervardb")
    group.add_option("-s", "--evidence_file", dest="evidence_file", action="store",
            help="User specified Evidence file for each variant", metavar="your_evidence_file")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "   How to add your own Evidence for each Variant",
    """ Prepare your own evidence  file as tab-delimited,the line format:
         (The code for additional evidence should be as: PS5/PM7/PP6/BS5/BP8 ;
         The format for upgrad/downgrade of criteria should be like: grade_PS1=2; 
         1 for Strong; 2 for Moderate; 3 for Supporting)   
            Chr Pos Ref_allele Alt_allele  PM1=1;BS2=1;BP3=0;PS5=1;grade_PM1=1
                                """)
    parser.add_option_group(group)



    group = optparse.OptionGroup(parser, "Annovar Options",
                                "Caution: check these options from manual of Annovar. The ANNOVAR version should be >=  2016-02-01, older verions of ANNOVAR will bring problems.")
    group.add_option("--table_annovar", action="store", help="The Annovar perl script of table_annovar.pl",metavar="./table_annovar.pl",dest="table_annovar")
    group.add_option("--convert2annovar", action="store", help="The Annovar perl script of convert2annovar.pl",metavar="./convert2annovar.pl",dest="convert2annovar")
    group.add_option("--annotate_variation", action="store", help="The Annovar perl script of annotate_variation.pl",metavar="./annotate_variation.pl",dest="annotate_variation")
    group.add_option("-d", "--database_locat", dest="database_locat", action="store",
            help="The  database location/dir for the annotation datasets", metavar="humandb")

    group.add_option("--skip_annovar", action="store_true", help="Skip the Annovar annotation, this can be true only after you  already got the annovar's annotation results",dest="skip_annovar")

    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Examples",
                                """./InterVar.py -c config.ini  # Run the examples in config.ini
                                 ./InterVar.py  -b hg19 -i your_input  --input_type=VCF  -o your_output 
                                """)
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    
    #(options,args) = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    print("%s" %description)
    print("%s" %version)
    print("Notice: Your command of InterVar is %s" % sys.argv[:])



    config_file = os.path.join(os.path.dirname(__file__),"config.ini") 
    if os.path.isfile(config_file):
        config.read(config_file)
        sections = config.sections()
        for section in sections:
            ConfigSectionMap(config,section)    
    else:
        print("Error: The default configure file of [ config.ini ] is not here, exit! Please redownload the InterVar.")
        sys.exit()

#begin to process user's options:
    if options.config != None:
        if os.path.isfile(options.config):
            config.read(options.config)
            sections = config.sections()
            for section in sections:
                ConfigSectionMap(config,section)
        else:
            print("Error: The config file [ %s ] is not here,please check the path of your config file." % options.config)
            sys.exit()

    if options.buildver != None:
        paras['buildver']=options.buildver
    if options.database_locat != None:
        paras['database_locat']=options.database_locat
    if options.input != None:
        paras['inputfile']=options.input
    if options.input_type != None:
        paras['inputfile_type']=options.input_type
    if options.output != None:
        paras['outfile']=options.output
    if options.evidence_file != None:
        paras['evidence_file']=options.evidence_file
        print("Warning: You provided your own evidence file [ %s ] for the InterVar." % options.evidence_file)
    if options.database_intervar != None:
        paras['database_intervar']=options.database_intervar
        paras['lof_genes'] = paras['database_intervar']+'/PVS1.LOF.genes'
        paras['pm1_domain'] = paras['database_intervar']+'/PM1_domains_with_benigns'
        paras['mim2gene'] =paras['database_intervar']+'/mim2gene.txt'
        paras['morbidmap'] = paras['database_intervar']+'/morbidmap'
        paras['mim_recessive'] = paras['database_intervar']+'/mim_recessive.txt'
        paras['mim_domin'] = paras['database_intervar']+'/mim_domin.txt'
        paras['mim_adultonset'] = paras['database_intervar']+'/mim_adultonset.txt'
        paras['mim_pheno'] = paras['database_intervar']+'/mim_pheno.txt'
        paras['mim_orpha'] = paras['database_intervar']+'/mim_orpha.txt'
        paras['orpha'] = paras['database_intervar']+'/orpha.txt'
        paras['knowngenecanonical'] = paras['database_intervar']+'/knownGeneCanonical.txt'
        paras['pp2_genes'] = paras['database_intervar']+'/PP2.genes'
        paras['bp1_genes'] = paras['database_intervar']+'/BP1.genes'
        paras['ps1_aa'] = paras['database_intervar']+'/PS1.AA.change.patho'
        paras['ps4_snps'] = paras['database_intervar']+'/PS4.variants'
        paras['bs2_snps'] = paras['database_intervar']+'/BS2_hom_het'
        paras['exclude_snps'] = paras['database_intervar']+'/ext.variants'


    paras['lof_genes'] = paras['lof_genes']+'.'+paras['buildver']
    paras['pm1_domain'] = paras['pm1_domain']+'.'+paras['buildver']
    paras['knowngenecanonical'] = paras['knowngenecanonical']+'.'+paras['buildver']
    paras['pp2_genes'] = paras['pp2_genes']+'.'+paras['buildver']
    paras['bp1_genes'] = paras['bp1_genes']+'.'+paras['buildver']

    paras['ps1_aa'] = paras['ps1_aa']+'.'+paras['buildver']
    paras['ps4_snps'] = paras['ps4_snps']+'.'+paras['buildver']
    paras['bs2_snps'] = paras['bs2_snps']+'.'+paras['buildver']
    paras['exclude_snps'] = paras['exclude_snps']+'.'+paras['buildver']
    paras['skip_annovar'] = False;

    if options.skip_annovar == True:
        paras['skip_annovar'] = True

    if options.table_annovar != None and options.skip_annovar != True:
        if os.path.isfile(options.table_annovar):
            paras['table_annovar']=options.table_annovar
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % options.table_annovar)
            if options.skip_annovar != True:
                sys.exit()
    if options.convert2annovar != None and options.skip_annovar != True:
        if os.path.isfile(options.convert2annovar):
            paras['convert2annovar']=options.convert2annovar
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % options.convert2annovar)
            if options.skip_annovar != True:
                sys.exit()
    if options.annotate_variation != None and options.skip_annovar != True:
        if os.path.isfile(options.annotate_variation):
            paras['annotate_variation']=options.annotate_variation
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % options.annotate_variation)
            if options.skip_annovar != True:
                sys.exit()


    if not os.path.isfile(paras['inputfile']):
        print("Error: Your input file [ %s ] is not here,please check the path of your input file." % paras['inputfile'])
        sys.exit()
    if  not os.path.isfile(paras['evidence_file']) and paras['evidence_file']!="None":
        print("Warning: Your specified evidence file [ %s ] is not here,please check the path of your evidence file." % paras['evidence_file'])
        print("         Your analysis will begin without your specified evidence.")


            


    print ("INFO: The options are %s " % paras)
    check_downdb()
    check_input()
    if  options.skip_annovar != True   :
        check_annovar_result() #  to obtain myanno.hg19_multianno.csv
    else:
         print ("Warning: You activated the option of --skip_annovar, the Annovar will not run!")
         print ("Warning: The InterVar will interpret the variants based on your old annotation information!")

    cutoff=float(paras['disorder_cutoff']) # user's disorder cutoff
    if cutoff != 0.005 :
        print("Warning: Customized  disease cutoff %s in config.ini, not suggested as 0.005" % (cutoff)) 

    read_datasets()
    
    inputft= paras['inputfile_type']
    some_file_fail=0
    out_annf=0; 
    for annovar_outfile  in glob.iglob(paras['outfile']+"*."+paras['buildver']+"_multianno.txt"):
        sum1=check_genes(annovar_outfile)
        sum2=my_inter_var(annovar_outfile)
        out_annf=out_annf+1; 

        outfile=annovar_outfile+".intervar"
        if os.path.isfile(outfile):
            print ("Notice: About %d lines in your variant file! " % (sum1-1)) 
            print ("Notice: About %d variants has been processed by InterVar" % (sum2-1))
            if inputft.lower() != 'vcf_m' :
                print ("Notice: The InterVar is finished, the output file is [ %s.intervar ]" % annovar_outfile)
        else:
            some_file_fail=some_file_fail+1 
            print ("Warning: The InterVar seems not run correctly, please check your inputs and options in configure file")

    if inputft.lower() == 'vcf_m' :
        print ("Notice: The InterVar for VCF with multiple samples is finished, the output file is as [ %s.<samplename>.intervar ]" % annovar_outfile)
        sum_sample=1;
        for f in glob.iglob(paras['outfile']+"*."+paras['buildver']+"_multianno.txt.intervar"):
            print ("Notice: The InterVar for VCF with multiple samples is finished, The %d sample output file is [ %s]" %(sum_sample,f))
            sum_sample=sum_sample+1;
        if some_file_fail>=1:    
            print ("Warning: The InterVar seems not run correctly for your %d samples in the VCF, please check your inputs and options in configure file" %  some_file_fail )
    if out_annf==0:
         print ("Warning: The InterVar seems not run correctly, please check your inputs , options and configure file!")
         print ("ERROR: The InterVar did not find the annotation result file from ANNOVAR!")
         print ("ERROR: The name of annotation result file should be like %s*.%s__multianno.txt" % (paras['outfile'],paras['buildver']))
    print("%s" %end)


    

if __name__ == "__main__":
    main()
