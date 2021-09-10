cd /home/ignatius/PostDoc/2021/ProteomeRiver

Rscript -e 'devtools::document()'

cd ..

tar -czvf ProteomeRiver.tar.gz ProteomeRiver/

R CMD INSTALL ProteomeRiver.tar.gz
