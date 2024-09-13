# The Living England Project

## About the Project
The Living England (LE) project, led by Natural England, is a multi-year programme delivering a satellite-derived national habitat map showing the extent and distribution of England’s diverse broad habitats. This is in support of the Defra Environmental Land Management (ELM) schemes and the National Capital and Ecosystem Assessment (NCEA) Programmes. The project uses European Space Agency Sentinel-1 and Sentinel-2 imagery, alongside additional open source datasets, machine learning modelling and is supported by a national targetted field data collection programme.

The current Living England habitat probability map for England is openly available to view on Defra's MAGIC platform, Defra's Data Services Platform, data.gov.uk and NE's Open data Platform, under an Open Government Licence v3.0.

## The Living_England repository

Living England is an ongoing long running project,  the code in this repository reflects the standardised method used to produce the Living England habitat maps from 2024 onwards (Living England 2022-23 map publication). The scripts are subject to change as improvements and optimisations are made to the reproducible analytical workflows. If you would like to be informed of any changes as they happen, please contact the team using the email address provided below.

The code has been specifically written for LE, and has intricacies that are unique to the project that will likely need to be updated or rewritten before being used on other projects. An example of this is the sub-division of England into BioGeographic Zones (bgzs). These are regional divisions to allow for the spatial variations in vegetation composition to be captured and to support acquisition of cloud-free imagery. The boundaries of these zones are based on the National Character Area boundaries and the Sentinel-2 satellite orbits, are are available from: <a href="https://www.data.gov.uk/dataset/8c36c913-eae3-411e-af2f-4d375d40a074/biogeographic-zones-living-england-2021)">Biogeographic Zones</a>.

The repository is split up into 'Workflow' scripts which call the filepaths and functions used in each step of the workflow, and 'Function' scripts where discrete functions are reused for particular tasks. These follow the processing workflow detailed in the diagram below.

<img src="LE_workflow_diagram.png" alt="Living England workflow diagram" width="700" height="600">

### Copyright

© Natural England 2024

### Contact
The Living England team are happy to answer any queries regarding the project, and welcome any suggested improvements and/or bugs reports for the code. The LE team are also keen to hear from anyone who is using the Living England code and applying it in their own work. Please contact the LE team at livingenglandenquiries@naturalengland.org.uk

### Acknowledgements
The Natural England Living England project is funded through the Defra National Capital and Ecosystem Assessment (NCEA), 
and Environmental Land Management schemes (ELMs) programmes. Thanks go to the many colleagues and partners within Natural England, Defra group and JNCC and to the landowners and surveyors engaged witht the project for their collaboration on the project.




