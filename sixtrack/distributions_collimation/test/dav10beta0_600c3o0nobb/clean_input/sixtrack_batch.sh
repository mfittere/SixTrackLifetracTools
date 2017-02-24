#!/bin/bash

# ALWAYS check:
# beam 1 or 2
# batch queue


# LocalPWD : where on the mac should the files be copied. Does the directory exist?
# post-processing programs to run
# which output files should be copied back
# should previous dirs run* be deleted


PWD=`pwd`
LocalPWD="/media/Daniele/RunII/Debris/2016/B1_IP1/"

beam=b1
LIMIT=354

echo $PWD
scp -r $PWD/clean_input dmirarch@pcen33066:"${LocalPWD}"
scp $PWD/sixtrack_batch.sh dmirarch@pcen33066:"${LocalPWD}clean_input"






rm -r run*

for ((a=1; a <= LIMIT ; a++))
  do
  index=$a
  zero=0
  while [ "${#index}" -lt "4" ]
    do
    index=$zero$index
  done 
  mkdir run$index

  cat > run$index/SixTr$index.job << EOF
  #!/bin/bash
  cp $PWD/clean_input/* .
  cp /afs/cern.ch/work/d/dmirarch/RunII/inputs/input_inelDeb.dat .

# a = pack number 100 = number of packs
  ./make_distr input_inelDeb.dat $a 100
# sixtrack executable
  ./SixTrack_SVN4515_debug_final > screenout 
  ./BeamLossPattern_2005-04-30_gcc2.9 lowb tracks2.dat BLP_out allapert.b1
  perl -pi -e 's/\0/ /g' LPI_BLP_out.s
  ./CleanInelastic_2013-08-19 FLUKA_impacts.dat LPI_BLP_out.s CollPositionsRunII.b1.dat coll_summary.dat



  cp coll_summary.dat FirstImpacts.dat FLUKA_impacts.dat impacts_fake.dat LPI_BLP_out.s $PWD/run$index/
#  cp coll_summary.dat efficiency.dat $PWD/run$index/


  scp -r $PWD/run${index} dmirarch@pcen33066:"${LocalPWD}"
  rm -rf $PWD/run${index}
  exit
EOF
if [ -d "run$index" ]; then
    cd run$index
    chmod 777 SixTr$index.job
    bsub -q 8nh -R "rusage[pool=50000]" SixTr$index.job
    cd ..
fi
  
done
