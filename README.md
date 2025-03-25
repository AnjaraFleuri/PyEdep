
## Description
Python class to Read/Write/Plot the [edep-sim](https://github.com/ClarkMcGrew/edep-sim/tree/master/io) files.

## Example plots

```
python test.py
```

Can also run in Jupyter Notebook.

## Writer
Summarize list of edep-sim root files into a customized output tree
```py
python3 writer.py 'input_file*.root' 'Marley' output_file.root
```

The branches are defined as follows:
| Branch        | unit           | description  |
| :------------ |:-------------| :-----|
| nu_pdg        | int | neutrino pdg code |
| nu_xs      | cm^2      | interaction cross section |
| nu_proc | int | 10+X: CC; 20+X: NC; (X: 1: QES; 2: RES; 3: DIS; 4: COH, 5: MEC)  |
| nu_nucl | int | interact with proton or neutron |
| E_nu | MeV | true neutrino energy |
| E_avail | MeV | available energy from initial state particles, including mass if meson or lepton |
| E_availList | array, MeV | 0: mu/e; 1: proton; 2: neutron; 3: pi+/-; 4: pi0; 5: gamma; 6: alpha; 7: others |
| E_depoTotal | MeV | total deposited energy  |
| E_depoTotal_track | MeV | total deposited energy for all tracks longer than 2cm  |
| Q_depoTotal | MeV | Deposited energy in charge with Birks Model  |
| Q_depoTotal_th_75keV | MeV | Deposited energy in charge with threshold dQ>75keV with Birks Model |
| Q_depoTotal_th_500keV | MeV | Deposited energy in charge with threshold dQ>500keV with Birks Model |
| Q_depoTotal_dots_th_75keV | MeV | Deposited energy in charge for all blips shorter than 2cm with threshold dQ>75keV with Birks Model |
| Q_depoTotal_MBox | MeV | Deposited energy in charge the Modified Box Model |
| Q_depoTotal_MBox_th_75keV | MeV | Deposited energy in charge with threshold dQ>75keV with the Modified Box Model |
| Q_depoTotal_MBox_th_500keV | MeV | Deposited energy in charge with threshold dQ>500keV with the Modified Box Model |
| L_depoTotal_avg_220PEpMeV | MeV | Deposited energy in light with mean LY 220PE/MeV with Birks Model |
| L_depoTotal_avg_180PEpMeV | MeV | Deposited energy in light with mean LY 180PE/MeV with Birks Model |
| L_depoTotal_avg_140PEpMeV | MeV | Deposited energy in light with mean LY 140PE/MeV with Birks Model |
| L_depoTotal_avg_100PEpMeV | MeV | Deposited energy in light with mean LY 100PE/MeV with Birks Model |
| L_depoTotal_avg_35PEpMeV | MeV | Deposited energy in light with mean LY 35PE/MeV with Birks Model |
| L_depoTotal_MBox_avg_220PEpMeV | MeV | Deposited energy in light with mean LY 220PE/MeV with the Modified Box Model  |
| L_depoTotal_MBox_avg_180PEpMeV | MeV | Deposited energy in light with mean LY 180PE/MeV with the Modified Box Model  |
| L_depoTotal_MBox_avg_140PEpMeV | MeV | Deposited energy in light with mean LY 140PE/MeV with the Modified Box Model  |
| L_depoTotal_MBox_avg_100PEpMeV | MeV | Deposited energy in light with mean LY 100PE/MeV with the Modified Box Model  |
| L_depoTotal_MBox_avg_35PEpMeV | MeV | Deposited energy in light with mean LY 35PE/MeV with the Modified Box Model  |
| E_depoList | array, MeV | similar to E_availList but for deposited energy including all children |
| E_depoList_track | array, MeV | similar to E_availList but for deposited energy for tracks longer than 2cm |
| Q_depoList | array, MeV | similar to E_depoList but for deposited charge with Birks Model |
| Q_depoList_th_75keV | array, MeV | similar to Q_depoList but with a 75keV cut on deposited charge with Birks Model |
| Q_depoList_th_500keV | array, MeV | similar to Q_depoList but with a 500keV cut on deposited charge with Birks Model |
| Q_depoList_dots_th_75keV | array, MeV | similar to Q_depoList but for blips shorter than 2cm with a 75keV cut on deposited charge with Birks Model |
| Q_depoList_MBox | array, MeV | similar to E_depoList but for deposited charge with the Modified Box Model |
| Q_depoList_MBox_th_75keV | array, MeV | similar to Q_depoList but with a 75keV cut on deposited charge with the Modified Box Model |
| Q_depoList_MBox_th_500keV | array, MeV | similar to Q_depoList but with a 500keV cut on deposited charge with the Modified Box Model |
| L_depoList_avg_220PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 220PE/MeV with Birks Model |
| L_depoList_avg_180PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 180PE/MeV with Birks Model |
| L_depoList_avg_140PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 140PE/MeV with Birks Model |
| L_depoList_avg_100PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 100PE/MeV with Birks Model |
| L_depoList_avg_35PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 35PE/MeV with Birks Model |
| L_depoList_MBox_avg_220PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 220PE/MeV with the Modified Box Model |
| L_depoList_MBox_avg_180PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 180PE/MeV with the Modified Box Model |
| L_depoList_MBox_avg_140PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 140PE/MeV with the Modified Box Model |
| L_depoList_MBox_avg_100PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 100PE/MeV with the Modified Box Model |
| L_depoList_MBox_avg_35PEpMeV | array, MeV | similar to E_depoList but for detected light with detector mean light yield 35PE/MeV with the Modified Box Model |
| N_parList | array, int | 0: mu/e; 1: proton; 2: neutron; 3: pi+/-; 4: pi0; 5: gamma; 6: alpha; 7: others |
