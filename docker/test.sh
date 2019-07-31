sudo docker build -t pipeline .

echo 'docker image create successfully'

sudo docker run -it \
-v $(pwd)/Pipeline:/home/pipeline/Pipeline \
-w /home/pipeline/Pipeline \
pipeline \
python pipeline.py \
--run-ccube \
--run-pyclone \
--run-phylowgs \
--phylowgs-burnin-samples 2 \
--phylowgs-mcmc-samples 5 \
--phylowgs-mh-iterations 2 \
--pyclone-mcmc-iterations 5 \
--pyclone-burnin-samples 2 \
--run-nmf --run-lda --lda-nmf-input data/Data/Phase_two/Subclonal_Structure/ \
--samples-to-run data/selectfile.txt
