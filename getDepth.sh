#PBS -l walltime=400:00:00
#PBS -l mem=2500gb
#PBS -l nodes=1:ppn=1
#PBS -M ian.beddows@vai.org
#PBS -m abe
#PBS -N GROseq

cd /home/ian.beddows/go/GROP_20190515_GroSeq/

perl src/get_region_from_depthfile.pl -region_file induced_repressed_data_for_file_creation_updated.txt -depth_file deliverables/depth_of_coverage.filtered.txt -outfile induced_repressed_depth_updated.txt
