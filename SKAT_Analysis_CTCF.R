#This is the full code run for the analysis on the transcription factor marker CTCF using the SKAT package in R. This format was followed for the other three transcription factor markers.
install.packages("SKAT")

library(SKAT)

setwd("Research_Project_Data/Plink_Output/CTCF")

Generate_SSD_SetID(File.Bed = "1_CTCF.bed", File.Bim = "1_CTCF.bim", File.Fam = "1_updated_CTCF.fam", File.SetID = "1_CTCF.SetID", File.SSD = "1_CTCF.SSD", File.Info = "1_CTCF.info")

FAM<-Read_Plink_FAM(Filename = "1_updated_CTCF.fam", Is.binary=FALSE)
y<-FAM$Phenotype

SSD.INFO<-Open_SSD("1_CTCF.SSD", "1_CTCF.info")

obj<-SKAT_Null_Model(y ~ 1, out_type="D")

out<-SKAT.SSD.All(SSD.INFO, obj)
out

Close_SSD()


#####Chr2################
Generate_SSD_SetID(File.Bed = "2_CTCF.bed", File.Bim = "2_CTCF.bim", File.Fam = "2_updated_CTCF.fam", File.SetID = "2_CTCF.SetID", File.SSD = "2_CTCF.SSD", File.Info = "2_CTCF.info")

FAM_2<-Read_Plink_FAM(Filename = "2_updated_CTCF.fam", Is.binary=FALSE)
y_2<-FAM_2$Phenotype

SSD.INFO_2<-Open_SSD("2_CTCF.SSD", "2_CTCF.info")

obj_2<-SKAT_Null_Model(y_2 ~ 1, out_type="D")

out_2<-SKAT.SSD.All(SSD.INFO_2, obj_2)
out_2

Close_SSD()


#####Chr3################
Generate_SSD_SetID(File.Bed = "3_CTCF.bed", File.Bim = "3_CTCF.bim", File.Fam = "3_updated_CTCF.fam", File.SetID = "3_CTCF.SetID", File.SSD = "3_CTCF.SSD", File.Info = "3_CTCF.info")

FAM_3<-Read_Plink_FAM(Filename = "3_updated_CTCF.fam", Is.binary=FALSE)
y_3<-FAM_3$Phenotype

SSD.INFO_3<-Open_SSD("3_CTCF.SSD", "3_CTCF.info")

obj_3<-SKAT_Null_Model(y_3 ~ 1, out_type="D")

out_3<-SKAT.SSD.All(SSD.INFO_3, obj_3)
out_3

Close_SSD()


#####Chr4################
Generate_SSD_SetID(File.Bed = "4_CTCF.bed", File.Bim = "4_CTCF.bim", File.Fam = "4_updated_CTCF.fam", File.SetID = "4_CTCF.SetID", File.SSD = "4_CTCF.SSD", File.Info = "4_CTCF.info")

FAM_4<-Read_Plink_FAM(Filename = "4_updated_CTCF.fam", Is.binary=FALSE)
y_4<-FAM_4$Phenotype

SSD.INFO_4<-Open_SSD("4_CTCF.SSD", "4_CTCF.info")

obj_4<-SKAT_Null_Model(y_4 ~ 1, out_type="D")

out_4<-SKAT.SSD.All(SSD.INFO_4, obj_4)
out_4

Close_SSD()


#####Chr5################
Generate_SSD_SetID(File.Bed = "5_CTCF.bed", File.Bim = "5_CTCF.bim", File.Fam = "5_updated_CTCF.fam", File.SetID = "5_CTCF.SetID", File.SSD = "5_CTCF.SSD", File.Info = "5_CTCF.info")

FAM_5<-Read_Plink_FAM(Filename = "5_updated_CTCF.fam", Is.binary=FALSE)
y_5<-FAM_5$Phenotype

SSD.INFO_5<-Open_SSD("5_CTCF.SSD", "5_CTCF.info")

obj_5<-SKAT_Null_Model(y_5 ~ 1, out_type="D")

out_5<-SKAT.SSD.All(SSD.INFO_5, obj_5)
out_5

Close_SSD()


#####Chr6################
Generate_SSD_SetID(File.Bed = "6_CTCF.bed", File.Bim = "6_CTCF.bim", File.Fam = "6_updated_CTCF.fam", File.SetID = "6_CTCF.SetID", File.SSD = "6_CTCF.SSD", File.Info = "6_CTCF.info")

FAM_6<-Read_Plink_FAM(Filename = "6_updated_CTCF.fam", Is.binary=FALSE)
y_6<-FAM_6$Phenotype

SSD.INFO_6<-Open_SSD("6_CTCF.SSD", "6_CTCF.info")

obj_6<-SKAT_Null_Model(y_6 ~ 1, out_type="D")

out_6<-SKAT.SSD.All(SSD.INFO_6, obj_6)
out_6

Close_SSD()


#####Chr7################
Generate_SSD_SetID(File.Bed = "7_CTCF.bed", File.Bim = "7_CTCF.bim", File.Fam = "7_updated_CTCF.fam", File.SetID = "7_CTCF.SetID", File.SSD = "7_CTCF.SSD", File.Info = "7_CTCF.info")

FAM_7<-Read_Plink_FAM(Filename = "7_updated_CTCF.fam", Is.binary=FALSE)
y_7<-FAM_7$Phenotype

SSD.INFO_7<-Open_SSD("7_CTCF.SSD", "7_CTCF.info")

obj_7<-SKAT_Null_Model(y_7 ~ 1, out_type="D")

out_7<-SKAT.SSD.All(SSD.INFO_7, obj_7)
out_7

Close_SSD()


#####Chr8################
Generate_SSD_SetID(File.Bed = "8_CTCF.bed", File.Bim = "8_CTCF.bim", File.Fam = "8_updated_CTCF.fam", File.SetID = "8_CTCF.SetID", File.SSD = "8_CTCF.SSD", File.Info = "8_CTCF.info")

FAM_8<-Read_Plink_FAM(Filename = "8_updated_CTCF.fam", Is.binary=FALSE)
y_8<-FAM_8$Phenotype

SSD.INFO_8<-Open_SSD("8_CTCF.SSD", "8_CTCF.info")

obj_8<-SKAT_Null_Model(y_8 ~ 1, out_type="D")

out_8<-SKAT.SSD.All(SSD.INFO_8, obj_8)
out_8

Close_SSD()


#####Chr9################
Generate_SSD_SetID(File.Bed = "9_CTCF.bed", File.Bim = "9_CTCF.bim", File.Fam = "9_updated_CTCF.fam", File.SetID = "9_CTCF.SetID", File.SSD = "9_CTCF.SSD", File.Info = "9_CTCF.info")

FAM_9<-Read_Plink_FAM(Filename = "9_updated_CTCF.fam", Is.binary=FALSE)
y_9<-FAM_9$Phenotype

SSD.INFO_9<-Open_SSD("9_CTCF.SSD", "9_CTCF.info")

obj_9<-SKAT_Null_Model(y_9 ~ 1, out_type="D")

out_9<-SKAT.SSD.All(SSD.INFO_9, obj_9)
out_9

Close_SSD()


#####Chr10################
Generate_SSD_SetID(File.Bed = "10_CTCF.bed", File.Bim = "10_CTCF.bim", File.Fam = "10_updated_CTCF.fam", File.SetID = "10_CTCF.SetID", File.SSD = "10_CTCF.SSD", File.Info = "10_CTCF.info")

FAM_10<-Read_Plink_FAM(Filename = "10_updated_CTCF.fam", Is.binary=FALSE)
y_10<-FAM_10$Phenotype

SSD.INFO_10<-Open_SSD("10_CTCF.SSD", "10_CTCF.info")

obj_10<-SKAT_Null_Model(y_10 ~ 1, out_type="D")

out_10<-SKAT.SSD.All(SSD.INFO_10, obj_10)
out_10

Close_SSD()


#####Chr11################
Generate_SSD_SetID(File.Bed = "11_CTCF.bed", File.Bim = "11_CTCF.bim", File.Fam = "11_updated_CTCF.fam", File.SetID = "11_CTCF.SetID", File.SSD = "11_CTCF.SSD", File.Info = "11_CTCF.info")

FAM_11<-Read_Plink_FAM(Filename = "11_updated_CTCF.fam", Is.binary=FALSE)
y_11<-FAM_11$Phenotype

SSD.INFO_11<-Open_SSD("11_CTCF.SSD", "11_CTCF.info")

obj_11<-SKAT_Null_Model(y_11 ~ 1, out_type="D")

out_11<-SKAT.SSD.All(SSD.INFO_11, obj_11)
out_11

Close_SSD()


#####Chr12################
Generate_SSD_SetID(File.Bed = "12_CTCF.bed", File.Bim = "12_CTCF.bim", File.Fam = "12_updated_CTCF.fam", File.SetID = "12_CTCF.SetID", File.SSD = "12_CTCF.SSD", File.Info = "12_CTCF.info")

FAM_12<-Read_Plink_FAM(Filename = "12_updated_CTCF.fam", Is.binary=FALSE)
y_12<-FAM_12$Phenotype

SSD.INFO_12<-Open_SSD("12_CTCF.SSD", "12_CTCF.info")

obj_12<-SKAT_Null_Model(y_12 ~ 1, out_type="D")

out_12<-SKAT.SSD.All(SSD.INFO_12, obj_12)
out_12

Close_SSD()


#####Chr13################
Generate_SSD_SetID(File.Bed = "13_CTCF.bed", File.Bim = "13_CTCF.bim", File.Fam = "13_updated_CTCF.fam", File.SetID = "13_CTCF.SetID", File.SSD = "13_CTCF.SSD", File.Info = "13_CTCF.info")

FAM_13<-Read_Plink_FAM(Filename = "13_updated_CTCF.fam", Is.binary=FALSE)
y_13<-FAM_13$Phenotype

SSD.INFO_13<-Open_SSD("13_CTCF.SSD", "13_CTCF.info")

obj_13<-SKAT_Null_Model(y_13 ~ 1, out_type="D")

out_13<-SKAT.SSD.All(SSD.INFO_13, obj_13)
out_13

Close_SSD()


#####Chr14################
Generate_SSD_SetID(File.Bed = "14_CTCF.bed", File.Bim = "14_CTCF.bim", File.Fam = "14_updated_CTCF.fam", File.SetID = "14_CTCF.SetID", File.SSD = "14_CTCF.SSD", File.Info = "14_CTCF.info")

FAM_14<-Read_Plink_FAM(Filename = "14_updated_CTCF.fam", Is.binary=FALSE)
y_14<-FAM_14$Phenotype

SSD.INFO_14<-Open_SSD("14_CTCF.SSD", "14_CTCF.info")

obj_14<-SKAT_Null_Model(y_14 ~ 1, out_type="D")

out_14<-SKAT.SSD.All(SSD.INFO_14, obj_14)
out_14

Close_SSD()


#####Chr15################
Generate_SSD_SetID(File.Bed = "15_CTCF.bed", File.Bim = "15_CTCF.bim", File.Fam = "15_updated_CTCF.fam", File.SetID = "15_CTCF.SetID", File.SSD = "15_CTCF.SSD", File.Info = "15_CTCF.info")

FAM_15<-Read_Plink_FAM(Filename = "15_updated_CTCF.fam", Is.binary=FALSE)
y_15<-FAM_15$Phenotype

SSD.INFO_15<-Open_SSD("15_CTCF.SSD", "15_CTCF.info")

obj_15<-SKAT_Null_Model(y_15 ~ 1, out_type="D")

out_15<-SKAT.SSD.All(SSD.INFO_15, obj_15)
out_15

Close_SSD()


#####Chr16################
Generate_SSD_SetID(File.Bed = "16_CTCF.bed", File.Bim = "16_CTCF.bim", File.Fam = "16_updated_CTCF.fam", File.SetID = "16_CTCF.SetID", File.SSD = "16_CTCF.SSD", File.Info = "16_CTCF.info")

FAM_16<-Read_Plink_FAM(Filename = "16_updated_CTCF.fam", Is.binary=FALSE)
y_16<-FAM_16$Phenotype

SSD.INFO_16<-Open_SSD("16_CTCF.SSD", "16_CTCF.info")

obj_16<-SKAT_Null_Model(y_16 ~ 1, out_type="D")

out_16<-SKAT.SSD.All(SSD.INFO_16, obj_16)
out_16

Close_SSD()


#####Chr17################
Generate_SSD_SetID(File.Bed = "17_CTCF.bed", File.Bim = "17_CTCF.bim", File.Fam = "17_updated_CTCF.fam", File.SetID = "17_CTCF.SetID", File.SSD = "17_CTCF.SSD", File.Info = "17_CTCF.info")

FAM_17<-Read_Plink_FAM(Filename = "17_updated_CTCF.fam", Is.binary=FALSE)
y_17<-FAM_17$Phenotype

SSD.INFO_17<-Open_SSD("17_CTCF.SSD", "17_CTCF.info")

obj_17<-SKAT_Null_Model(y_17 ~ 1, out_type="D")

out_17<-SKAT.SSD.All(SSD.INFO_17, obj_17)
out_17

Close_SSD()


#####Chr18################
Generate_SSD_SetID(File.Bed = "18_CTCF.bed", File.Bim = "18_CTCF.bim", File.Fam = "18_updated_CTCF.fam", File.SetID = "18_CTCF.SetID", File.SSD = "18_CTCF.SSD", File.Info = "18_CTCF.info")

FAM_18<-Read_Plink_FAM(Filename = "18_updated_CTCF.fam", Is.binary=FALSE)
y_18<-FAM_18$Phenotype

SSD.INFO_18<-Open_SSD("18_CTCF.SSD", "18_CTCF.info")

obj_18<-SKAT_Null_Model(y_18 ~ 1, out_type="D")

out_18<-SKAT.SSD.All(SSD.INFO_18, obj_18)
out_18

Close_SSD()


#####Chr19################
Generate_SSD_SetID(File.Bed = "19_CTCF.bed", File.Bim = "19_CTCF.bim", File.Fam = "19_updated_CTCF.fam", File.SetID = "19_CTCF.SetID", File.SSD = "19_CTCF.SSD", File.Info = "19_CTCF.info")

FAM_19<-Read_Plink_FAM(Filename = "19_updated_CTCF.fam", Is.binary=FALSE)
y_19<-FAM_19$Phenotype

SSD.INFO_19<-Open_SSD("19_CTCF.SSD", "19_CTCF.info")

obj_19<-SKAT_Null_Model(y_19 ~ 1, out_type="D")

out_19<-SKAT.SSD.All(SSD.INFO_19, obj_19)
out_19

Close_SSD()


#####Chr20################
Generate_SSD_SetID(File.Bed = "20_CTCF.bed", File.Bim = "20_CTCF.bim", File.Fam = "20_updated_CTCF.fam", File.SetID = "20_CTCF.SetID", File.SSD = "20_CTCF.SSD", File.Info = "20_CTCF.info")

FAM_20<-Read_Plink_FAM(Filename = "20_updated_CTCF.fam", Is.binary=FALSE)
y_20<-FAM_20$Phenotype

SSD.INFO_20<-Open_SSD("20_CTCF.SSD", "20_CTCF.info")

obj_20<-SKAT_Null_Model(y_20 ~ 1, out_type="D")

out_20<-SKAT.SSD.All(SSD.INFO_20, obj_20)
out_20

Close_SSD()


#####Chr21################
Generate_SSD_SetID(File.Bed = "21_CTCF.bed", File.Bim = "21_CTCF.bim", File.Fam = "21_updated_CTCF.fam", File.SetID = "21_CTCF.SetID", File.SSD = "21_CTCF.SSD", File.Info = "21_CTCF.info")

FAM_21<-Read_Plink_FAM(Filename = "21_updated_CTCF.fam", Is.binary=FALSE)
y_21<-FAM_21$Phenotype

SSD.INFO_21<-Open_SSD("21_CTCF.SSD", "21_CTCF.info")

obj_21<-SKAT_Null_Model(y_21 ~ 1, out_type="D")

out_21<-SKAT.SSD.All(SSD.INFO_21, obj_21)
out_21

Close_SSD()


#####Chr22###############
Generate_SSD_SetID(File.Bed = "22_CTCF.bed", File.Bim = "22_CTCF.bim", File.Fam = "22_updated_CTCF.fam", File.SetID = "22_CTCF.SetID", File.SSD = "22_CTCF.SSD", File.Info = "22_CTCF.info")

FAM_22<-Read_Plink_FAM(Filename = "22_updated_CTCF.fam", Is.binary=FALSE)
y_22<-FAM_22$Phenotype

SSD.INFO_22<-Open_SSD("22_CTCF.SSD", "22_CTCF.info")

obj_22<-SKAT_Null_Model(y_22 ~ 1, out_type="D")

out_22<-SKAT.SSD.All(SSD.INFO_22, obj_22)
out_22

Close_SSD()


#####ChrX###############
Generate_SSD_SetID(File.Bed = "X_CTCF.bed", File.Bim = "X_CTCF.bim", File.Fam = "X_updated_CTCF.fam", File.SetID = "X_CTCF.SetID", File.SSD = "X_CTCF.SSD", File.Info = "X_CTCF.info")

FAM_X<-Read_Plink_FAM(Filename = "X_updated_CTCF.fam", Is.binary=FALSE)
y_X<-FAM_X$Phenotype

SSD.INFO_X<-Open_SSD("X_CTCF.SSD", "X_CTCF.info")

obj_X<-SKAT_Null_Model(y_X ~ 1, out_type="D")

out_X<-SKAT.SSD.All(SSD.INFO_X, obj_X)
out_X

Close_SSD()
