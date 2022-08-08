import os

input_file = open("xlms_data.csv", 'r')
output_file = open("master.csv", 'w')

################################################################################
############################## To make master.csv ##############################
################################################################################

for line_num in input_file.readlines():
    #print(line_num.strip())
    fields = line_num.strip().split(',')
    #print(fields)
    if fields[0].startswith("Link"):
#        output_file.write(line_num)
        continue
    if fields[1]=="MTA1" or fields[1]=="HDAC1" or fields[1]=="MBD3" or fields[1]=="P66A" or fields[1]=="RBBP4":
        if fields[3]=="MTA1" or fields[3]=="HDAC1" or fields[3]=="MBD3" or fields[3]=="P66A" or fields[3]=="RBBP4":
            output_file.write(line_num)
            #print(line_num.strip())

input_file.close()
output_file.close()


################################################################################
############################ To make XL master files ###########################
################################################################################

output_file = open("master.csv", 'r')
adh_master_file = open("adh_master.dat", 'w')
bs3dss_master_file = open("bs3dss_master.dat", 'w')
dmtmm_master_file = open("dmtmm_master.dat", 'w')

adh_total=0
bs3dss_total=0
dmtmm_total=0
i=0
for ln in output_file.readlines():
    entries = ln.strip().split(",")
    i=i+1
    #print(entries)
    if entries[0].startswith("AD"):
        adh_master_file.write(ln)
        adh_total=adh_total+1
    if entries[0].startswith("BS3"):
        bs3dss_master_file.write(ln)
        bs3dss_total=bs3dss_total+1
    if entries[0].startswith("DMTM"):
        dmtmm_master_file.write(ln)
        dmtmm_total=dmtmm_total+1

print("Total XLs: ", adh_total+bs3dss_total+dmtmm_total)
output_file.close()
adh_master_file.close()
bs3dss_master_file.close()
dmtmm_master_file.close()


# ################################################################################
# ##################### To separate inter and intra XL files #####################
# ################################################################################
#
# adh_file = open("adh_master.csv", 'r')
# bs3dss_file = open("bs3dss_master.csv", 'r')
# dmtmm_file = open("dmtmm_master.csv", 'r')
#
# adh_intra_file = open("adh_intra.dat", 'a')
# bs3dss_intra_file = open("bs3dss_intra.dat", 'a')
# dmtmm_intra_file = open("dmtmm_intra.dat", 'a')
#
# adh_inter_file = open("adh_inter.dat", 'a')
# bs3dss_inter_file = open("bs3dss_inter.dat", 'a')
# dmtmm_inter_file = open("dmtmm_inter.dat", 'a')
#
# a=0
# b=0
# c=0
# intra_adh=0
# inter_adh=0
# intra_bs3dss=0
# inter_bs3dss=0
# intra_dmtmm=0
# inter_dmtmm=0
#
# ############################## For ADH ##############################
# for ln in adh_file.readlines():
#     adh=ln.strip().split(',')
#     a=a+1
#
#     if adh[1]==adh[3]:
#         intra_adh=intra_adh+1
#         adh_intra_file.write(ln)
#     if adh[1]!=adh[3]:
#         inter_adh=inter_adh+1
#         adh_inter_file.write(ln)
# print("ADH Total : ", a)
# print("Intra-ADH: ", intra_adh, "Inter-ADH: ", inter_adh, "ADH Total=(inter+intra): ", intra_adh+inter_adh)
#
# ############################## For BS3DSS ##############################
# for ln in bs3dss_file.readlines():
#     bs3dss=ln.strip().split(',')
#     b=b+1
#
#     if bs3dss[1]==bs3dss[3]:
#         intra_bs3dss=intra_bs3dss+1
#         bs3dss_intra_file.write(ln)
#     if bs3dss[1]!=bs3dss[3]:
#         inter_bs3dss=inter_bs3dss+1
#         bs3dss_inter_file.write(ln)
# print("\nBS3DSS Total: ", b)
# print("Intra-BS3DSS: ", intra_bs3dss, "Inter-BS3DSS: ", inter_bs3dss, "BS3DSS-Total=(inter+intra): ", intra_bs3dss+inter_bs3dss)
#
# ############################## For DMTMM ##############################
# for ln in dmtmm_file.readlines():
#     dmtmm=ln.strip().split(',')
#     c=c+1
#
#     if dmtmm[1]==dmtmm[3]:
#         intra_dmtmm=intra_dmtmm+1
#         dmtmm_intra_file.write(ln)
#     if dmtmm[1]!=dmtmm[3]:
#         inter_dmtmm=inter_dmtmm+1
#         dmtmm_inter_file.write(ln)
# print("\nDMTMM-Total: ", c)
# print("Intra-DMTMM: ", intra_dmtmm, "Inter-DMTMM: ", inter_dmtmm, "DMTMM-Total=(inter+intra): ", intra_dmtmm+inter_dmtmm)
#
# adh_file.close()
# bs3dss_file.close()
# dmtmm_file.close()
#
# adh_intra_file.close()
# bs3dss_intra_file.close()
# dmtmm_intra_file.close()
#
# adh_inter_file.close()
# bs3dss_inter_file.close()
# dmtmm_inter_file.close()


# ################################################################################
# ########################### To remove intra-PDB XLs ############################
# ################################################################################
#
# adh_master_file = open("adh_master.dat", 'r')
# bs3dss_master_file = open("bs3dss_master.dat", 'r')
# dmtmm_master_file = open("dmtmm_master.dat", 'r')
#
# adh_file = open("adh.dat", 'w')
# bs3dss_file = open("bs3dss.dat", 'w')
# dmtmm_file = open("dmtmm.dat", 'w')
#
# i=0
# j=0
# k=0
#
# ################################### For ADH ####################################
# for ln in adh_master_file.readlines():
#     adh = ln.strip().split(',')
#     if adh[1]==adh[3] and adh[1]=="MBD3":
#         if adh[2]>
