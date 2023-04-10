Updating Databases
===================

Updating the database is supper simple and for most you will just run::

    dram2 -d <where to put the db> build_db

To put the config in a specific location::

    dram2 -d <where to put the db> build_db -y <yaml_loc>

To update::

    dram2 -d <where to put the db> build_db  -u

To make only a few dbs::

    dram2 -d <where to put the db> build_db  --make_db <db_name>


What about kegg and camper
------------------------

KEGG will allwase need to be independently downloaded and CAMPER to wile its repository is private. The input argument lets you do this. 

Download
^^^^^^^

For camper just go to the release page on the github page and download the tag_gz release.

for kegg use the scripts in this repository.

