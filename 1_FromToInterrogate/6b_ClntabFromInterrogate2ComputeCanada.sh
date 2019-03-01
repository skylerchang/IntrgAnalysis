#!/bin/bash

#cd to equivalent of '/intrg_storage/users/smkeller/analyses/Run02_feTR_normal_Q30L50/last_run_res/Run02_feTR_normal--results'

#create folder in which to copy all clntab files
mkdir ../Clntab
#move clntab files to folder
for d in *001/; do echo $d; echo $d"/"*".clntab"; cp $d"/"*".clntab" ../Clntab/; done
#compress folder
tar -czvf ../Clntab.tar.gz ../Clntab

#********** transfer tar ball to compute canada (from compute canada or local) ***************
rsync -ave ssh keeper@s0.arrest.tools:/intrg_storage/users/smkeller/analyses/Run02_feTR_normal_Q30L50/last_run_res/Clntab.tar.gz ./


#clean-up
rm -rf ../Clntab


#============== option 2: transfer the entire "--results" folder 

cd /intrg_storage/users/smkeller/analyses/Run02_feTR_normal_Q30L50/last_run_res
tar -czvf Clntab.tar.gz Run02_feTR_normal_Q30L50--results/

#********** transfer tar ball to compute canada (from compute canada or local) ***************
rsync -ave ssh keeper@s0.arrest.tools:/intrg_storage/users/smkeller/analyses/Run02_feTR_normal_Q30L50/last_run_res/Clntab.tar.gz ./

tar -xzvf Clntab.tar.gz 
mv Run02_feTR_normal_Q30L50--results/ ClntabPlus/

#create folder in which to copy all clntab files
mkdir ../Clntab
#move clntab files to folder
for d in *001/; do echo $d; echo $d"/"*".clntab"; cp $d"/"*".clntab" ../Clntab/; done