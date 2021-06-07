
PROTECTED_FOLDERNAME=protected_maf_hg38
PUBLIC_FOLDERNAME=public_somatic_maf_hg38
CLINICAL_FOLDERNAME=clinical
DATA_FOLDER=data_tcga_wes
gdc_path=~/dl_tools_centos/gdc-client/bin



mkdir $DATA_FOLDER

cd TCGA_WES_raw_dl
$gdc_path/gdc-client download -m utils/gdc_manifest.2019-05-19.txt -t utils/gdc-user-token.2019-06-27T14_44_52.487Z.txt

for folder in `ls .`; do 
if [[ $folder = "utils" ]]; then
echo "just the manifest"
else
mutation_file=`ls $folder/*.gz`
#mutation_file=`ls e$folder/*.maf`
cancer_loc=`echo $mutation_file | cut -d "/" -f 2 | cut -d "." -f 2`;
gunzip $mutation_file
full_filename=`ls $folder/*.maf`
filename=`echo $full_filename | cut -d "/" -f 2`;
protection_status=`echo $mutation_file | cut -d "/" -f 2 | cut -d "." -f 7`;
mkdir -p ../$DATA_FOLDER/$cancer_loc/$PROTECTED_FOLDERNAME
mkdir -p ../$DATA_FOLDER/$cancer_loc/$PUBLIC_FOLDERNAME
echo $filename
if [ "$protection_status" = "protected" ]; then 
    ln -s $PWD/$full_filename ../$DATA_FOLDER/$cancer_loc/$PROTECTED_FOLDERNAME/$filename
elif [ "$protection_status" = "somatic" ]; then 
    ln -s $PWD/$full_filename ../$DATA_FOLDER/$cancer_loc/$PUBLIC_FOLDERNAME/$filename
fi
if [ ! -f "../$DATA_FOLDER/$cancer_loc/annotations.txt" ] ; then
    ln -s $PWD/$folder/annotations.txt ../$DATA_FOLDER/$cancer_loc/annotations.txt
fi
fi
done

mkdir clinical_data
cd clinical_data
for cancer_loc in ACC BRCA CHOL DLBC GBM KICH KIRP LGG LUAD MESO PAAD PRAD SARC STAD THCA UCEC UVM BLCA CESC COAD ESCA HNSC KIRC LAML LIHC LUSC OV PCPG READ SKCM TGCT THYM UCS ; do
lower_cancer_loc=`echo "${cancer_loc,,}"`
echo $cancer_loc
echo $lower_cancer_loc
done
wget http://download.cbioportal.org/${lower_cancer_loc}\_tcga.tar.gz
mkdir $cancer_loc
cd $cancer_loc
tar -xvf ../${lower_cancer_loc}\_tcga.tar.gz
cd ..
mkdir -p ../../$DATA_FOLDER/$cancer_loc/$CLINICAL_FOLDERNAME
ln -s $PWD/$cancer_loc/data_bcr_clinical_data_patient.txt ../../$DATA_FOLDER/$cancer_loc/$CLINICAL_FOLDERNAME/data_bcr_clinical_data_patient.txt ;
ln -s $PWD/$cancer_loc/data_bcr_clinical_data_sample.txt ../../$DATA_FOLDER/$cancer_loc/$CLINICAL_FOLDERNAME/data_bcr_clinical_data_sample.txt ;
done

