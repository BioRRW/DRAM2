#! /bin/bash

echo -n Enter KEGG user name: 
read UNAME
echo -n Enter KEGG password: 
read -s PASS
echo pulling from ftp://$UNAME:$PASS@ftp.kegg.net/kegg/README.kegg
# make directory structure
mkdir -p kegg/genes/
wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/README.kegg
wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/RELEASE
wget -r ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/MD5.genes
wget -r  -A "*\.pep\.gz" ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/organisms
wget -r  -A "*\.kff\.gz" ftp://$UNAME:$PASS@ftp.kegg.net/kegg/genes/organisms

gzip -d ftp.kegg.net/kegg/genes/organisms/*/*pep.gz
gzip -d ftp.kegg.net/kegg/genes/organisms/*/*\.kff\.gz

#cat ftp.kegg.net/kegg/genes/organisms/*/*pep >> kegg-all-orgs_20220129.pep
#cat ftp.kegg.net/kegg/genes/organisms/*/*kff >> kegg-all-orgs_20220129.kff

