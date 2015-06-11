source ~/.tcshrc
rehash
echo $PYTHONPATH

set dir=$1
set beg=$2

if ($#argv < 1) then
   echo "Usage: csh analyze_DTI.csh <dir>"
   exit
endif

#cd $dir
#if ($#argv == 2) then
#   python ~/bin/DTI_analysis.py lambda_ -b $beg
#   python ~/bin/DTI_analysis2.py lambda_
#   python ~/bin/rename_dti.py lambda_ l$beg
#else 
#   python ~/bin/DTI_analysis.py lambda_
#   python ~/bin/DTI_analysis2.py lambda_
#   python ~/bin/rename_dti.py lambda_ all
#endif
#
#cd ..


#first 2ns
#python ~/bin/DTI_analysis.py lambda_ -e 2000
#python ~/bin/DTI_analysis2.py lambda_
#python ~/bin/rename_dti.py lambda_ f2ns

#last 8ns
#python ~/bin/DTI_analysis.py lambda_ -b 2000
#python ~/bin/DTI_analysis2.py lambda_      
#python ~/bin/rename_dti.py lambda_ l8ns

#last 2ns
python ~/bin/DTI_analysis.py lambda_ -b 8000
python ~/bin/DTI_analysis2.py lambda_       
python ~/bin/rename_dti.py lambda_ l2ns      

#all
#python ~/bin/DTI_analysis.py lambda_
#python ~/bin/DTI_analysis2.py lambda_
#python ~/bin/rename_dti.py lambda_ all

