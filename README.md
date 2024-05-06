# SyntheticElongation
CC3D Simulation Files for Synthetic Elongation Project @ Morsut Lab (Curated by Christian Chung)

The files provided in this repository provide the coding files sed to produce the simulations seen in Courte et al., 2024. 

**Files Provided in this Repository**
1) Unparameterized Elongation Model - An unparameterized version of our synthetic elongation system, showcasing ideal elongation based on arbitrary values implemented into CompuCell3D. (Only the base file is provided for space purposes, the metrics altered for the parameter screen can be seen in the Supplementary Excel file showcasing the genomes). 
2) Phase Separation Models - Models highlighting changes in adhesion metrics that contribute to our system's phase separation between grey and red cells (attests to the signaling speed control).
3) Proliferation and Viscocapillary Speed Models - Models testing various _in silico_ proliferation metrics to parameterize their _in vitro_ equivalences. (Only the base file is provided for space purposes, the metrics altered can be seen in the Supplementary Excel file showcasing the genomes).
4) Bounded Morphospace Models - Models demonstrating changes in proliferation control, phase separation, and viscocapillary speed utilized for the morphospace.
5) _in silico_ Recapitulations of Current Elongation Models - _in silico_ Recapitulation of _in vitro_ datasets placed in the morphospace. 

**Getting Started**
1) First, install CompuCell3D v4.1.1 (CC3D) from the CompuCell3D website (https://compucell3d.org/). The version may change at the time of accessing this repository but for our files, v4.1.1 was used.
2) Install CompuCell3D as per the installation file.

**Running the Code**
1) Download the codes in this repository. You can download them individually or use the code tab to download everything as a zip file (as per your interests).

  ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/d43c7036-b13c-4ae5-9e9e-ca85034339fb)

2) Once the code and CompuCell3D have been downloaded and installed, run the Twedit++ file. Below is an image of these files that are contained within the CompuCell3D folder that is downloaded.

  ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/3f7227e7-2de7-44c1-9298-577826301852)

**Running the Files Using Twedit++**
1) Open the Twedit++ file.
2) Open your file by going to "CC3D Project" and selecting "Open CC3D Project."
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/d209932c-054a-4b96-a25f-ece8ea8cffb1)
3) Double click the CC3D file provided in one of the folders for the code in this repository.
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/c405ec57-c3dc-429e-afb7-566efbb16a3f)
4) Double click the simulation name to open a drop down of all the coding files in the folder.
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/a45182be-7956-4a52-9320-408ed1955752)
5) Right click the simulation name and press "Open In Player" to run the programming file.
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/8d47e574-4640-4d76-84bd-f6f8df8481f1)

**Troubleshooting**

If your simulation is having difficulty running or crashes upon initiation of the simulation, this may be due to your system requirements. Simulations are set to 8 cores for parallelization but can be altered to accommodate your system's needs. 
1) After open the CC3D file, navigate to the ELUGMSteppables.py file.
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/7124480d-a585-4479-9914-e4ddea9643ec)
2) Scroll to Line 45 and set the USEDNODES = 1 (or whatever is appropriate for your system).
   ![image](https://github.com/TDEL-SyntheticBio/SyntheticElongation/assets/168860358/1376ce9e-aac4-40d5-81b1-78f54f351d78)


If you have any other questions, please contact me at christianchung5@gmail.com


