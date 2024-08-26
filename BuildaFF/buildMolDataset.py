# use env build-a-force-field

from openff.toolkit.topology import Molecule

# QC Submit packages
from openff.qcsubmit.common_structures import QCSpec
from openff.qcsubmit.factories import TorsiondriveDatasetFactory
from openff.qcsubmit.workflow_components import (
    Scan1D, ScanEnumerator, StandardConformerGenerator
)
# Compute the training set, running calculations
from qcfractal import FractalSnowflakeHandler
from qcportal.client import FractalClient

# first build your moecules that you want to run torsion scans on
training_smiles= ['CC', 'CCCC']

training_moleucles = [Molecule.from_smiles(smile) for smile in training_smiles]

## we now have 3 molecules, how do we turn these into a set of QC calucaltions?

torsion_drive_factory = TorsiondriveDatasetFactory(
    qc_specifications={
        "default": QCSpec(
            method="B3LYP-D3BJ", basis="DZVP", program="psi4"
        )
    },
    workflow=[
        ScanEnumerator(
            torsion_scans=[
                Scan1D(
                    smarts1="[*:1]-[#6:2]-[#6:3]-[*:4]",
                    scan_range1=(-150, 180),
                    scan_increment=[30]
                )
            ]
        ),
        StandardConformerGenerator(toolkit='rdkit',max_conformers=10)
    ],
)

## Create a Torsiondrive Dataset Factor that will produce a dataset of QC torsion drives

# supply: 
#       1. a set of QC Specifications (define the level of theory as well as the program with which we want to run out calcs.
#       we have used the openFF defaults which were carefully chosen to provide a balance of speed and accuracy at reproducing conformational energies)
#       2. a workflow (set of processing steps that will be applied to a set of moleucles. Gnerate a set of conformers to start torsion drives from, 
#       enumerating all tautomers and protomers, fragmenting input molecules to speed up QC calcs)


# conformer generator to make sure we have several conformers to start the torsion drive from 
# smarts pattern to define exactly what torsions should be scanned. Match all unique torsions around a C-C single bond



torsion_drive_dataset = torsion_drive_factory.create_dataset(
    molecules=training_moleucles, 
    dataset_name='Alkane torsion Drives',
    description='contains two linear alkanes, ethane and butane',
    tagline='basic alkane dataset',
)

#make a temporary local QCFractal to better showcase how generating new QC data is:
local_fractal_instance = FractalSnowflakeHandler(ncores=16)
local_fractal_client = FractalClient(local_fractal_instance)

#submit, compute, and storing a QC dataset is as easy as calling the submit command and providing the address of the server taht should handle the submittion
submission = torsion_drive_dataset.submit(local_fractal_client)
local_fractal_client.query_procedures(procedure='torsiondrive', status=None)