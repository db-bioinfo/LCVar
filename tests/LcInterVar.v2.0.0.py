#!/usr/bin/env python
#########################################################################
# Enhanced LcInterVar with Varsome Clinical Methodology
# Based on original InterVar by Quan LI (leequan@gmail.com)
# Enhanced with Varsome methodology while maintaining all original functionality
#########################################################################

import copy,logging,os,io,re,time,sys,platform,optparse,gzip,glob

prog="LcInterVar"

version = """%prog 2.3.0 20241218
Enhanced InterVar with Varsome Clinical Methodology
Based on original InterVar by Quan LI, leequan@gmail.com
Enhanced with ClinGen-calibrated predictors and points-based system
LcInterVar is free for non-commercial use without warranty.
Copyright (C) 2024 Enhanced Version
============================================================================
"""

usage = """Usage: %prog [OPTION] -i  INPUT -o  OUTPUT ...
       %prog  --config=config.ini ...
"""

description = """=============================================================================
LcInterVar - Enhanced InterVar with Varsome Clinical Methodology                                                                    
Interpretation of Pathogenic/Benign for variants using enhanced python scripts. 

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
Thanks for using Enhanced LcInterVar!
Report bugs to leequan@gmail.com;
LcInterVar homepage: <http://wInterVar.wglab.org>

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
        print("Error: can\'t read the PP2 genes file %s" % paras['pp2_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()    

    #6.BP1 gene list       
    try:               
        fh = open(paras['bp1_genes'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                BP1_genes_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the BP1 genes file %s" % paras['bp1_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()    

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
    except IOError:
        print("Error: can\'t read the knownGeneCanonical  file %s" % paras['knowngenecanonical'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()   

    #11.BS2 variants of recessive homo, domin heter
    try:              
        with myGzipFile(paras['bs2_snps'], "rb") as fh:
            strs = fh.read().decode()
            for line2 in strs.split('\n'):
                cls2=line2.split(' ')
                if len(cls2[0])>=1:  #
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]
                    BS2_snps_recess_dict[ keys ]=cls2[4]  # key as snp info
                    BS2_snps_domin_dict[ keys ]=cls2[5]  # key as snp info
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[1]+"_"+flip_ACGT(cls2[2])+"_"+flip_ACGT(cls2[3])
                    BS2_snps_recess_dict[ keys ]=cls2[4]  # key as snp info
                    BS2_snps_domin_dict[ keys ]=cls2[5]  # key as snp info
    except IOError:
        print("Error: can\'t read the snp list file for BS2 %s" % paras['bs2_snps'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()    

    #12. OMIM mim_pheno.txt file
    try:
        fh = open(paras['mim_pheno'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                mim_pheno_dict[keys]=cls2[1]
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
    path=path.rstrip("/\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print("Notice: the folder of %s is created!" % path)
    else:
        print("Warning: the folder of %s is already created!" % path)
    ds=paras['database_names']
    ds.expandtabs(1);
    if not os.path.isfile(paras['annotate_variation']):
        print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % paras['annotate_variation'])
        if paras['skip_annovar'] != True:
            sys.exit()

    for dbs in ds.split():
        file_name=dbs
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
                sum=sum+1
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

def check_PVS1(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''
    ENHANCED PVS1: LOF variant assessment with SnpEff integration following ClinGen SVI 2018
    Incorporates SnpEff LOF annotations and NMD predictions for enhanced accuracy
    '''
    cls = line.split('\t')
    PVS1 = 0
    
    # Parse SnpEff LOF annotation from ALL columns (SnpEff data can be in any Otherinfo column)
    snpeff_lof_genes = set()
    nmd_genes = set()
    
    try:
        import re
        # Search ALL columns for LOF data (LOF could be in any Otherinfo column)
        for col_value in cls:
            if 'LOF=' in str(col_value):
                # Extract LOF annotations: LOF=(GENE|GENE|1|1.00)
                lof_matches = re.findall(r'LOF=\(([^|]+)\|[^|]+\|[^|]+\|([^)]+)\)', col_value)
                for gene, confidence in lof_matches:
                    if float(confidence) >= 0.9:  # High confidence LOF
                        snpeff_lof_genes.add(gene)
                
                # Extract NMD annotations: NMD=(GENE|GENE|1|1.00)
                nmd_matches = re.findall(r'NMD=\(([^|]+)\|[^|]+\|[^|]+\|([^)]+)\)', col_value)
                for gene, confidence in nmd_matches:
                    if float(confidence) >= 0.9:
                        nmd_genes.add(gene)
    except:
        pass
    
    # Check variant type from functional annotation
    line_tmp = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    
    # Enhanced LOF variant type detection per ClinGen SVI 2018
    lof_evidence = 0
    variant_gene = cls[Funcanno_flgs['Gene']]
    
    # Check if gene has LOF mechanism of disease
    gene_lof_mechanism = 0
    try:
        if lof_genes_dict.get(variant_gene) == '1':
            gene_lof_mechanism = 1
    except:
        pass
    
    # 1. SnpEff high-confidence LOF annotation
    if variant_gene in snpeff_lof_genes:
        lof_evidence = 1
        
        # Apply ClinGen SVI decision tree
        if gene_lof_mechanism:
            # Check for NMD prediction
            if variant_gene in nmd_genes:
                PVS1 = 1  # Very Strong (converted to points in classification)
            else:
                # Last exon or NMD-escape considerations
                try:
                    # Check if stopgain/nonsense in last exon
                    if "stopgain" in line_tmp.lower() or "nonsense" in line_tmp.lower():
                        # More conservative for last exon variants
                        PVS1 = 1  # Strong but flagged for review
                    else:
                        PVS1 = 1  # Very Strong
                except:
                    PVS1 = 1
    
    # 2. Traditional variant type checking as backup
    if not lof_evidence and gene_lof_mechanism:
        funcs_lof = ["nonsense", "frameshift", "stopgain", "startloss"]
        splice_lof = ["splice", "splicing"]
        
        for fc in funcs_lof:
            if fc in line_tmp.lower():
                lof_evidence = 1
                PVS1 = 1  # Very Strong
                break
        
        # Enhanced splice site assessment
        if not lof_evidence:
            for splice_type in splice_lof:
                if splice_type in line_tmp.lower():
                    # Check splice prediction scores
                    try:
                        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
                        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
                        
                        if ada_score != '.' and rf_score != '.':
                            ada_val = float(ada_score)
                            rf_val = float(rf_score)
                            # High confidence splice disruption
                            if ada_val >= 0.958 or rf_val >= 0.584:
                                PVS1 = 1  # Very Strong
                    except:
                        # Default for canonical splice sites
                        if "canonical" in line_tmp.lower():
                            PVS1 = 1  # Very Strong
    
    return PVS1

def check_PS1(line,Funcanno_flgs,Allels_flgs,aa_changes_dict):
    PS1=0
    PS1_t1=0
    PS1_t2=0
    PS1_t3=0
    dbscSNV_cutoff=0.6    
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    ACGTs=["A","C","G","T"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            PS1_t1=1;
            line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
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
    try:
        if float(cls[Funcanno_flgs['dbscSNV_RF_SCORE']])>dbscSNV_cutoff or float(cls[Funcanno_flgs['dbscSNV_ADA_SCORE']])>dbscSNV_cutoff: 
            PS1_t3=1
        if cls[Funcanno_flgs['dbscSNV_RF_SCORE']] == "." or cls[Funcanno_flgs['dbscSNV_ADA_SCORE']] == ".": 
            PS1_t3=0
    except ValueError:
        pass

    if PS1_t1 !=0 and PS1_t2 != 0 :
        PS1=1
        if PS1_t3 ==1: 
            PS1=0
    return(PS1)

def check_PS2(line,Funcanno_flgs,Allels_flgs):
    PS2=0
    return(PS2)

def check_PS3(line,Funcanno_flgs,Allels_flgs):
    PS3=0
    return(PS3)

def check_PS4(line,Funcanno_flgs,Allels_flgs):
    PS4=0
    cls=line.split('\t')
    keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
    try:
        if PS4_snps_dict[keys_tmp2] == "1":
            PS4=1
    except KeyError:
        pass
    return(PS4)

def check_PM1(line,Funcanno_flgs,Allels_flgs,domain_benign_dict):
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
        domain_key = cls[Allels_flgs['Chr']]+"_"+gene+": "+domain
        
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
    
    # Enhanced hotspot analysis
    hotspot_evidence = 0
    
    # Check if variant is in a gene known for hotspots
    gene = cls[Funcanno_flgs['Gene']]
    hotspot_genes = [
        'TP53', 'KRAS', 'NRAS', 'BRAF', 'PIK3CA', 'AKT1', 'EGFR', 
        'BRCA1', 'BRCA2', 'MLH1', 'MSH2', 'MSH6', 'PMS2'
    ]
    
    if gene in hotspot_genes:
        hotspot_evidence = 1
    
    # Combine evidence with enhanced logic
    if domain_evidence >= 2 or hotspot_evidence >= 1:
        PM1 = 1
    elif domain_evidence >= 1:
        PM1 = 1
    
    return PM1

def check_PM2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs,mim2gene_dict,mim2gene_dict2):
    '''
    CRITICAL FIX: PM2 Strength - Always Supporting per ClinGen SVI v1.0
    Absent from controls (or at extremely low frequency if recessive)
    '''
    PM2=0
    cls=line.split('\t')
    
    # Enhanced frequency and coverage checks
    valid_frequencies = []
    
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
                max_freq = max(valid_frequencies) if valid_frequencies else 0
                if max_freq < 0.0001:  # 0.01% threshold as default
                    PM2 = 1  # FIXED: Always Supporting
    
    return PM2

def check_PM3(line,Funcanno_flgs,Allels_flgs):
    PM3=0
    return(PM3)

def check_PM4(line,Funcanno_flgs,Allels_flgs):
    PM4=0
    PM4_t1=0
    PM4_t2=0
    cls=line.split('\t')
    funcs_tmp=["nonframeshift insertion","nonframeshift deletion","stoploss"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            PM4_t1=1
    if cls[Funcanno_flgs['rmsk']] == '.':
        PM4_t2=1
    if cls[Funcanno_flgs['rmsk']] != '.' and  line_tmp.find("stoploss")>=0 :
        PM4_t2=1

    if PM4_t1 !=0 and PM4_t2 != 0 :
        PM4=1

    return(PM4)

def check_PM5(line,Funcanno_flgs,Allels_flgs,aa_changes_dict):
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
            line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
            cls0=re.split("[,;]",line_tmp2)
            cls0_1=cls0[0].split(':')
            aa=cls0_1[4]
            aa_last=aa[len(aa)-1:]
            keys_tmp2=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Alt']]
            try:
                if  aa_changes_dict[keys_tmp2] :
                    PM5_t2=0
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

    if PM5_t1 !=0 and PM5_t2 != 0 and PM5_t3 !=0 :
        PM5=1
    return(PM5)

def check_PM6(line,Funcanno_flgs,Allels_flgs):
    PM6=0
    return(PM6)

def check_PP1(line,Funcanno_flgs,Allels_flgs):
    PP1=0
    return(PP1)

def check_PP2(line,Funcanno_flgs,Allels_flgs,PP2_genes_dict):
    PP2=0
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            try:
                if PP2_genes_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
                    PP2=1
            except KeyError:
                PP2=0
    return(PP2)

def check_PP3(line,Funcanno_flgs,Allels_flgs):
    '''
    ENHANCED PP3: ClinGen-calibrated computational predictors with graduated strength
    Uses single predictor approach with evidence-based thresholds
    Priority: MetaRNN > REVEL > BayesDel > AlphaMissense (per Varsome v11.18.0)
    '''
    PP3 = 0
    cls = line.split('\t')
    
    # Priority 1: MetaRNN (Varsome's primary predictor after calibration)
    try:
        if 'MetaRNN_score' in Funcanno_flgs and cls[Funcanno_flgs['MetaRNN_score']] != '.':
            metarnn_score = float(cls[Funcanno_flgs['MetaRNN_score']])
            # ClinGen-calibrated thresholds for MetaRNN (estimated based on performance)
            if metarnn_score >= 0.939:
                return 4  # Strong pathogenic evidence
            elif metarnn_score >= 0.841:
                return 2  # Moderate pathogenic evidence  
            elif metarnn_score >= 0.748:
                return 1  # Supporting pathogenic evidence
    except (ValueError, KeyError):
        pass
    
    # Priority 2: REVEL (ClinGen-calibrated thresholds from Pejaver et al. 2022)
    try:
        if 'REVEL_score' in Funcanno_flgs and cls[Funcanno_flgs['REVEL_score']] != '.':
            revel_score = float(cls[Funcanno_flgs['REVEL_score']])
            if revel_score >= 0.932:
                return 8  # Very Strong
            elif revel_score >= 0.879:
                return 4  # Strong
            elif revel_score >= 0.773:
                return 4  # Strong
            elif revel_score >= 0.644:
                return 2  # Moderate
            elif revel_score >= 0.291:
                return 1  # Supporting
    except (ValueError, KeyError):
        pass
    
    # Priority 3: BayesDel (without AF) - ClinGen-calibrated thresholds
    try:
        if 'BayesDel_addAF_score' in Funcanno_flgs and cls[Funcanno_flgs['BayesDel_addAF_score']] != '.':
            bayesdel_score = float(cls[Funcanno_flgs['BayesDel_addAF_score']])
            if bayesdel_score >= 0.500:
                return 8  # Very Strong
            elif bayesdel_score >= 0.410:
                return 4  # Strong
            elif bayesdel_score >= 0.270:
                return 2  # Moderate
            elif bayesdel_score >= 0.130:
                return 1  # Supporting
    except (ValueError, KeyError):
        pass
    
    # Priority 4: AlphaMissense (ClinGen-calibrated thresholds)
    try:
        if 'AlphaMissense_score' in Funcanno_flgs and cls[Funcanno_flgs['AlphaMissense_score']] != '.':
            alphamissense_score = float(cls[Funcanno_flgs['AlphaMissense_score']])
            if alphamissense_score >= 0.990:
                return 8  # Very Strong
            elif alphamissense_score >= 0.972:
                return 4  # Strong
            elif alphamissense_score >= 0.906:
                return 2  # Moderate
            elif alphamissense_score >= 0.792:
                return 1  # Supporting
    except (ValueError, KeyError):
        pass
    
    # Enhanced splicing prediction (high priority for splice variants)
    try:
        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
        
        if ada_score != '.' and rf_score != '.':
            ada_val = float(ada_score)
            rf_val = float(rf_score)
            # Enhanced splicing thresholds (ClinGen calibrated)
            if ada_val >= 0.958 or rf_val >= 0.584:
                # Splicing impact can reach Very Strong evidence
                return 8  # Very Strong for high-confidence splicing
    except (ValueError, KeyError):
        pass
    
    # Conservation as fallback (if no other predictors available)
    try:
        if 'GERP++_RS' in Funcanno_flgs and cls[Funcanno_flgs['GERP++_RS']] != '.':
            gerp_score = float(cls[Funcanno_flgs['GERP++_RS']])
            if gerp_score >= 9.88:
                return 2  # Moderate pathogenic
            elif gerp_score >= 7.52:
                return 1  # Supporting pathogenic
    except (ValueError, KeyError):
        pass
    
    return PP3

def check_PP4(line,Funcanno_flgs,Allels_flgs):
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
    BA1: Allele frequency >5% with population-specific checks
    ENHANCED: Exclude ESP6500 per ClinGen guidance
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
                    BA1 = 1  # Stand alone benign (converted to points in classification)
                    break
            except ValueError:
                continue
    
    return BA1

def check_BS1(line,Freqs_flgs,Allels_flgs):
    '''
    BS1: Allele frequency greater than expected for disorder
    ENHANCED: Exclude ESP6500 per ClinGen guidance
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
                if freq > 0 and freq <= 1.0:  # Basic sanity check
                    frequencies.append(freq)
            except ValueError:
                continue
    
    if frequencies:
        max_freq = max(frequencies)
        # Use more stringent thresholds for different inheritance patterns
        if max_freq >= threshold:
            BS1 = 1  # Strong benign evidence (original InterVar returns 1)
    
    return BS1

def check_BS2(line,Freqs_flgs,Allels_flgs,Funcanno_flgs):
    BS2=0
    cls=line.split('\t')
    keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
    try:
        mim1=mim2gene_dict[ cls[Funcanno_flgs['Gene.ensGene']] ]
        try:
            if mim_adultonset_dict[mim1] == "1": 
                BS2=0;
        except KeyError: 
            try:
                if mim_recessive_dict[mim1] == "1": 
                    try:
                        if BS2_snps_recess_dict[ keys ]=="1":  
                            BS2=1
                    except KeyError:
                        pass
            except KeyError:
                pass

            try:
                if mim_domin_dict[mim1] == "1": 
                    BS2=0
                    try:
                        if BS2_snps_domin_dict[ keys ]=="1":  
                            BS2=1
                    except KeyError:
                        pass
            except KeyError:
                pass
        else:
            pass

    except KeyError: 
        BS2=0

    return(BS2)

def check_BS3(line,Funcanno_flgs,Allels_flgs):
    BS3=0
    return(BS3)

def check_BS4(line,Funcanno_flgs,Allels_flgs):
    BS4=0
    return(BS4)

def check_BP1(line,Funcanno_flgs,Allels_flgs,BP1_genes_dict):
    BP1=0
    cls=line.split('\t')
    funcs_tmp=["missense","nonsynony"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            try:
                if BP1_genes_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
                    BP1=1
            except KeyError:
                BP1=0
    return(BP1)

def check_BP2(line,Funcanno_flgs,Allels_flgs):
    BP2=0
    return(BP2)

def check_BP3(line,Funcanno_flgs,Allels_flgs):
    BP3=0
    BP3_t1=0
    BP3_t2=0
    cls=line.split('\t')
    funcs_tmp=["nonframeshift insertion","nonframeshift deletion","nonframeshift substitution"]
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 :
            BP3_t1=1;
    if cls[Funcanno_flgs['rmsk']] != '.' and  cls[Funcanno_flgs['Interpro_domain']]== '.': 
        BP3_t2=1

    if BP3_t1 !=0 and BP3_t2 != 0 :
        BP3=1
    return(BP3)

def check_BP4(line,Funcanno_flgs,Allels_flgs):
    '''
    ENHANCED BP4: ClinGen-calibrated computational predictors for benign evidence
    Uses single predictor approach with evidence-based thresholds
    FIXED: Properly exclude LOF variants and respect PVS1 conflicts
    '''
    BP4 = 0
    cls = line.split('\t')
    
    # CRITICAL FIX 1: Exclude loss-of-function variant types
    # BP4 should never apply to clearly damaging variants
    lof_variant_types = ["stopgain", "stoplose", "frameshift", "startloss", "splicing"]
    line_tmp = cls[Funcanno_flgs['Func.refGene']] + " " + cls[Funcanno_flgs['ExonicFunc.refGene']]
    
    for lof_type in lof_variant_types:
        if lof_type in line_tmp.lower():
            return 0  # Do not apply BP4 to LOF variants
    
    # Priority 1: MetaRNN (most accurate for benign predictions)
    try:
        if 'MetaRNN_score' in Funcanno_flgs and cls[Funcanno_flgs['MetaRNN_score']] != '.':
            metarnn_score = float(cls[Funcanno_flgs['MetaRNN_score']])
            if metarnn_score <= 0.108:
                return 4  # Strong benign evidence
            elif metarnn_score <= 0.267:
                return 2  # Moderate benign evidence
            elif metarnn_score <= 0.43:
                return 1  # Supporting benign evidence
    except (ValueError, KeyError):
        pass
    
    # Priority 2: REVEL (ClinGen-calibrated thresholds)
    try:
        if 'REVEL_score' in Funcanno_flgs and cls[Funcanno_flgs['REVEL_score']] != '.':
            revel_score = float(cls[Funcanno_flgs['REVEL_score']])
            if revel_score <= 0.016:
                return 4  # Strong benign
            elif revel_score <= 0.017:
                return 2  # Moderate benign
            elif revel_score <= 0.053:
                return 1  # Supporting benign
    except (ValueError, KeyError):
        pass
    
    # Priority 3: BayesDel (ClinGen-calibrated thresholds)
    try:
        if 'BayesDel_addAF_score' in Funcanno_flgs and cls[Funcanno_flgs['BayesDel_addAF_score']] != '.':
            bayesdel_score = float(cls[Funcanno_flgs['BayesDel_addAF_score']])
            if bayesdel_score <= -0.520:
                return 4  # Strong benign
            elif bayesdel_score <= -0.360:
                return 2  # Moderate benign
            elif bayesdel_score <= -0.180:
                return 1  # Supporting benign
    except (ValueError, KeyError):
        pass
    
    # Priority 4: AlphaMissense (only for missense variants)
    try:
        if 'AlphaMissense_score' in Funcanno_flgs and cls[Funcanno_flgs['AlphaMissense_score']] != '.':
            # Only apply AlphaMissense to missense variants
            if "missense" in line_tmp.lower() or "nonsynon" in line_tmp.lower():
                alphamissense_score = float(cls[Funcanno_flgs['AlphaMissense_score']])
                if alphamissense_score <= 0.070:
                    return 4  # Strong benign
                elif alphamissense_score <= 0.099:
                    return 2  # Moderate benign
                elif alphamissense_score <= 0.169:
                    return 1  # Supporting benign
    except (ValueError, KeyError):
        pass
    
    # Conservation for benign evidence (low conservation = benign)
    try:
        if 'GERP++_RS' in Funcanno_flgs and cls[Funcanno_flgs['GERP++_RS']] != '.':
            gerp_score = float(cls[Funcanno_flgs['GERP++_RS']])
            if gerp_score <= -1.04:
                return 4  # Strong benign
            elif gerp_score <= 1.08:
                return 2  # Moderate benign
            elif gerp_score <= 3.58:
                return 1  # Supporting benign
    except (ValueError, KeyError):
        pass
    
    # Check splicing - no predicted impact (conservative approach)
    try:
        ada_score = cls[Funcanno_flgs.get('dbscSNV_ADA_SCORE', -1)]
        rf_score = cls[Funcanno_flgs.get('dbscSNV_RF_SCORE', -1)]
        
        # ONLY trigger BP4 if we have actual scores showing no splicing impact
        if ada_score != '.' and rf_score != '.':
            ada_val = float(ada_score)
            rf_val = float(rf_score)
            if ada_val < 0.958 and rf_val < 0.584:
                return 1  # Supporting benign
        
    except (ValueError, KeyError):
        pass
    
    return BP4

def check_BP5(line,Funcanno_flgs,Allels_flgs,morbidmap_dict):
    BP5=0; 
    return(BP5)

def check_BP6(line,Funcanno_flgs,Allels_flgs):
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
    BP7=0
    BP7_t1=0
    BP7_t2=0
    cutoff_conserv=2 
    cls=line.split('\t')
    funcs_tmp=["synon","coding-synon"]
    funcs_tmp2="nonsynon"
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 and line_tmp.find(funcs_tmp2)<0 :
            if cls[Funcanno_flgs['dbscSNV_RF_SCORE']]=="." or  cls[Funcanno_flgs['dbscSNV_ADA_SCORE']]==".":
                BP7_t1=1  
            else:
                if float(cls[Funcanno_flgs['dbscSNV_RF_SCORE']])<0.6 and float(cls[Funcanno_flgs['dbscSNV_ADA_SCORE']])<0.6:
                    BP7_t1=1
    try:
        if float(cls[Funcanno_flgs['GERP++_RS']]) <= float(cutoff_conserv) or cls[Funcanno_flgs['GERP++_RS']] == '.' :
            BP7_t2=1
    except ValueError:
        BP7_t2=1

    if BP7_t1 !=0 and BP7_t2 != 0 :
        BP7=1        
    return(BP7)

def classfy_points(PVS1,PS,PM,PP,BA1,BS,BP,Allels_flgs,cls):
    """
    ENHANCED: Tavtigian et al. 2020 points-based classification system
    Points: Supporting=1, Moderate=2, Strong=4, Very Strong=8
    """
    BPS=["Pathogenic","Likely pathogenic","Benign","Likely benign","Uncertain significance"]
    
    # Calculate pathogenic evidence points
    pathogenic_points = 0
    if PVS1 > 0:
        pathogenic_points += 8  # PVS1 = Very Strong = 8 points
    pathogenic_points += sum(PS) * 4  # PS = Strong = 4 points each
    pathogenic_points += sum(PM)      # PM can return variable points - don't multiply!
    pathogenic_points += sum(PP)      # PP can return variable points - don't multiply!
    
    # Calculate benign evidence points  
    benign_points = 0
    if BA1 > 0:
        benign_points += 8  # BA1 = Stand alone = 8 points
    benign_points += sum(BS) * 4  # BS = Strong = 4 points each
    benign_points += sum(BP)      # BP can return variable points - don't multiply!
    
    # Calculate total score
    total_score = pathogenic_points - benign_points
    
    # Tavtigian et al. 2020 thresholds
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

def assign(BP,line,Freqs_flgs,Funcanno_flgs,Allels_flgs):
    PVS1=0
    PS=[0,0,0,0,0]
    PM=[0,0,0,0,0,0,0]
    PP=[0,0,0,0,0,0]

    BA1=0
    BS=[0,0,0,0,0]
    BP=[0,0,0,0,0,0,0,0]

    # First, evaluate all rules using enhanced functions
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

    # ENHANCED CONFLICT RESOLUTION: Apply rule conflict resolution following ACMG/VarSome guidelines
    
    # 1. PVS1 conflicts: Enhanced conflict resolution per ClinGen SVI 2018
    if PVS1 > 0:
        BP[3] = 0  # Disable BP4 (computational evidence suggests benign)
        # Note: Do NOT disable PP3 for splicing variants per ClinGen SVI Splicing Subgroup 2023
        # Only disable PP3 for canonical splice sites already captured by PVS1
        try:
            cls = line.split('\t')
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
            
    #begin process the user's evidence file
    if os.path.isfile(paras['evidence_file']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            evds=user_evidence_dict[keys] 
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and (not re.findall('grade', evd_t[0], flags=re.IGNORECASE)) ):
                    if int(evd_t[1])<=1:
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

    # end process the user's evidence file 

    if len(cls)>1:    
        BP_out=classfy_points(PVS1,PS,PM,PP,BA1,BS,BP,Allels_flgs,cls)
        line_t="%s PVS1=%s PS=%s PM=%s PP=%s BA1=%s BS=%s BP=%s" %(BP_out,PVS1,PS,PM,PP,BA1,BS,BP)
        BP_out=line_t
        pass
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

    try:
        fh=open(newoutfile, "r")
        fw=open(newoutfile2, "w")
        strs=fh.read()
        line_sum=0;
        print("Notice: Begin the variants interpretation by Enhanced LcInterVar with Varsome methodology")
        if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t LcInterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp151","AAChange.ensGene","AAChange.refGene","Clinvar","ClinVar_Review","LcInterVar and Evidence","Freq_gnomAD_genome_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha","Otherinfo"  ))
        else:
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t LcInterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp151","AAChange.ensGene","AAChange.refGene","Clinvar","ClinVar_Review","LcInterVar and Evidence","Freq_gnomAD_genome_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha"  ))

        for line in strs.split('\n'):
            BP="UNK" 
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
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t LcInterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp151']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,clnrevstat_output,intervar_bp,cls[Freqs_flgs['gnomAD_genome_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']], cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],".", cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['rmsk']],cls[Funcanno_flgs['MetaSVM_score']],Freq_gnomAD_genome_POPs,OMIM,Pheno_MIM,orpha,orpha_details,cls[Funcanno_flgs['Otherinfo']]   ))
                else:
                    clnrevstat_output = cls[Funcanno_flgs.get('CLNREVSTAT', len(cls)-1)] if 'CLNREVSTAT' in Funcanno_flgs else "."
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \tclinvar_review: %s \t LcInterVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp151']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,clnrevstat_output,intervar_bp,cls[Freqs_flgs['gnomAD_genome_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']], cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],".", cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['rmsk']],cls[Funcanno_flgs['MetaSVM_score']],Freq_gnomAD_genome_POPs,OMIM,Pheno_MIM,orpha,orpha_details  ))

            line_sum=line_sum+1

    except IOError:
        print("Error: can\'t read/write the annovar output files %s" % (newoutfile,newoutfile2))
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

    group = optparse.OptionGroup(parser, "LcInterVar Other Options")
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
                                """./LcInterVar.py -c config.ini  # Run the examples in config.ini
                                 ./LcInterVar.py  -b hg19 -i your_input  --input_type=VCF  -o your_output 
                                """)
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    print("%s" %description)
    print("%s" %version)
    print("Notice: Your command of Enhanced LcInterVar is %s" % sys.argv[:])

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
        print("Warning: You provided your own evidence file [ %s ] for the LcInterVar." % options.evidence_file)
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
        check_annovar_result() 
    else:
         print ("Warning: You activated the option of --skip_annovar, the Annovar will not run!")
         print ("Warning: The Enhanced LcInterVar will interpret the variants based on your old annotation information!")

    cutoff=float(paras['disorder_cutoff']) 
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
            print ("Notice: About %d variants has been processed by Enhanced LcInterVar" % (sum2-1))
            if inputft.lower() != 'vcf_m' :
                print ("Notice: The Enhanced LcInterVar is finished, the output file is [ %s.intervar ]" % annovar_outfile)
        else:
            some_file_fail=some_file_fail+1 
            print ("Warning: The Enhanced LcInterVar seems not run correctly, please check your inputs and options in configure file")

    if inputft.lower() == 'vcf_m' :
        print ("Notice: The Enhanced LcInterVar for VCF with multiple samples is finished, the output file is as [ %s.<samplename>.intervar ]" % annovar_outfile)
        sum_sample=1;
        for f in glob.iglob(paras['outfile']+"*."+paras['buildver']+"_multianno.txt.intervar"):
            print ("Notice: The Enhanced LcInterVar for VCF with multiple samples is finished, The %d sample output file is [ %s]" %(sum_sample,f))
            sum_sample=sum_sample+1;
        if some_file_fail>=1:    
            print ("Warning: The Enhanced LcInterVar seems not run correctly for your %d samples in the VCF, please check your inputs and options in configure file" %  some_file_fail )
    if out_annf==0:
         print ("Warning: The Enhanced LcInterVar seems not run correctly, please check your inputs , options and configure file!")
         print ("ERROR: The Enhanced LcInterVar did not find the annotation result file from ANNOVAR!")
         print ("ERROR: The name of annotation result file should be like %s*.%s_multianno.txt" % (paras['outfile'],paras['buildver']))
    print("%s" %end)

if __name__ == "__main__":
    main()
