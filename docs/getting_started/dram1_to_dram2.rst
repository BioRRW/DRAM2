.. _dram1_to_dram2:

Using DRAM2: for DRAM1 Users
=====================

Getting started with DRAM2 for those that are most familiar with DRAM1.


Welcome to Members of the Wrighton Lab!
------------------------------------

This document is meant for you only. You as a lab member are privy to the information on DRAM2 before anyone else. Please take care to note any concerns you have with this document, as it will be used to make a generalized version for the public.

The world will be very grateful for the time you take to look over this document!

The point of this document
--------------------------

As you know we are migrating to DRAM2 and that will work differently to DRAM1, I hope better. This document is intended to get new users in the Wrighton lab ready to use DRAM2 having already experienced DRAM1. The DRAM2 Documentation will be different to new users as compared to those who have already used DRAM2, but you still may get something out of it, so feel free to look at it [HERE](link). This document will focus mostly on the differences and new features in DRAM2, and less on how DRAM2 works fundamentally.

As you look over this document, notice that DRAM2 differs from DRAM1 in these key ways:

 * It has a project config that can store run data. This is used for checking, and lets us skip most arguments.
 * It has more steps than DRAM1, but each step is simpler.
 * DRAM-v is not part of DRAM2, it may be in the future, but not yet.
 * There is only one command `dram2` which lets you use all sub commands.
 * Only DRAM2 has Adjectives and Phylogenetic Trees

Preview
-------

Before we start I want to just lay out the simplest form of the pipeline that we will use as an example to compare DRAM1 and DRAM2. Like the points above, it would be good to keep these at the back of your mind while you read the rest of this document.

.. code-block:: bash

   DRAM.py annotate \
       	-i "" \
       	-o soil/dram1_out \
       	--use_camper --use_fegenie --use_sulfur



Grab an Example Directory
------------------------


Activate the Environment
------------------------

First let's activate the DRAM2 environment, this is just like with DRAM1 except it is DRAM2. The first version of DRAM2 is DRAM2.0.rc1 this name indicates that it is DRAM major release 2 minor release 0 or the first, and then it is a release candidate 1.
So run the command below, and we can get started testing this release candidate.

.. code-block:: bash

   source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
   source /opt/Miniconda2/miniconda2/bin/activate DRAM2BETA


Just a quick not on the setup. The DRAM2 config dose not live with DRAM2. The
global config on zenith is at `/ect/dram_config.yaml`. You as a user can make a
file in you home directories config file ~/.config/dram_config.yaml and that
will replace the global one for you only. We are not going to get into this much
but you should know, because that is a big change.

Explore The help
----------------

With any new program, it is good to explore the help. The DRAM2 help is a lot larger compared to the dram1 and really needs a good checking in order to validate. The DRAM2 command line structure is very hierarchical, much more than DRAM2, so here are some things I think will explain that best.

 * DRAM2 has one main command aka dram2 so no matter what you want to do with DRAM2 it starts with dram2 as apposed to DRAM.py and DRAM-setup.py
 * You provide general/universal arguments to dram2 and then specific arguments to the dram2 sub commands. That should be clear as we go through the help.
 * This will take just a moment to explore, just notice the différance between options and the

First, let's look at the overview:

.. code-block:: bash

	dram2 --help

At the end of this process, we want to be able to have adjectives generated, but the adjectives in DRAM2 have more requirements that need to be met first. To learn more about adjectives and see what specifically is required, you can run:

.. code-block:: bash

   dram2 adjectives --help
   dram2 adjectives eval --help

We will revisit adjectives regularly as we go through this example. Adjectives use several new features of DRAM2, such as database checking and phylogenetic trees.

The first step to any dram project however is probably calling genes.

Call Genes, and start a project
-------------------------------

In DRAM1 calling genes was part of the annotation process but now it is done with the call command and the annotation process only works on already called genes. This adds a step but makes the process a lot simpler for a lay person to follow

First please read the help and make sure you understand it. We will reiterate some of what it says in the next section however.


.. code-block:: bash

   dram2 call --help

The 15 soil genomes are a good place to start for dram. You could call all of the soil genomes with the command below. However, **I suggest you don't run this command**. It would take too long and too much memory even though calling genes is now multi-threaded.vb so let's just select 2 like in the next command.

Before we move on I want to talk about this command, which we are not running.  There are so many things to cover here.

  #. The commands that get passed to dram2 are universal and work with all dram2 sub-commands, but you don't pass those to the sub commands. So `dram2 call -0 soil/test1` would not work. The reverse is also true, you don't pass option to dram2 that go to the sub commands so `dram2 --prodigal_mode train call -0 soil/test1` would not work.  
  #. Additionally, `dram2 call` and `dram2 annotate` allow for a list of arguments after all the options. In both cases these are lists of fastas only one is for called fastas and the other is for uncalled fastas. Arguments are different from options in that they have no flags like no `-f` or `--flag`, they can't be fallowed by options and,  in this case, there can be as many arguments as you want. Before, the path to FASTA files had to be a string. That was ok, but it was confusing at times use a normal file path instead.   

Notice the output is specified by a -o. and it is passed to the dram2 command before the call command runs, the same with the -c command that tells dram the most cores it needs are 2.
In the past, DRAM confused people by having them pass a string to call genes with a python command. So now we let bash handle this. This should be safer and result in less errors.


.. code-block:: bash

   dram2 -o soil/test1 -c 15 call \
      	/home/projects-wrighton-2/DRAM/input_datasets/15_soil_genomes/all_data/*.fasta

I suggest a simple dataset with just two of the 15 soil genomes
rm -r soil/test1

.. code-block:: bash

   dram2 -vv -o soil/test1 -c 2 call \
        	/home/projects-wrighton-2/DRAM/input_datasets/15_soil_genomes/all_data/Cytophaga_hutchinsonii_ATCC_33406.fasta \
       	/home/projects-wrighton-2/DRAM/input_datasets/15_soil_genomes/all_data/Dechloromonas_aromatica_RCB.fasta

calling annotations
-------------------


.. code-block:: bash

  dram2 annotate --help


Calling annotations can be done with a db_set, as seen above, but it can also be done with the use_db flag individually. Seeing as the one above would take a long time you can use these to get a taste of annotations.


.. code-block:: bash

  dram2 -vvvv -o soil/test1 -c 30 annotate --use_db fegenie --use_db camper --use_db methyl

There are also some databases that you may not think of as databases like Heme Motif count and even the stars.


.. code-block:: bash

  dram2 -vvvv -o soil/test1 -c 30 annotate --use_db heme --use_db stats


.. code-block:: bash

  dram2 -vvvv -o soil/test1 -c 30 annotate --use_db dbcan --use_db peptidase --use_db kegg --use_db pfam  -f


.. code-block:: bash

  dram2 -vvvv -o soil/test1 -c 30 distill


.. code-block:: bash

   dram2 -vvvv -o soil/test1 -c 30 phylotree


.. code-block:: bash

   dram2 -vvvv -o soil/distill_this/ -c 30 adjectives eval


S
