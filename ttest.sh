for i in {0..12}; do OMP_NUM_THREADS=$i; sed "s/XXX/$i/" noes.ROC.0.ini > temp.ini; echo $i >> times.out; time Packages/User/noes-corr --config temp.ini >> times.out; done
