import ROOT
from ROOT import TG4Event, TFile, TChain

import sys, os
import numpy as np
import random
import math

class Event:

    def __init__(self, fileName, evgen='Genie'):
        print("event: initilization")
        self.fileName = fileName
        self.evgen = evgen
        self.ReadTree()

        self.currentEntry = 0
        # create a folder to store plots
        self.plotpath = "./plots"
        if not os.path.exists( self.plotpath):
            os.makedirs( self.plotpath)
            print("plotpath '" + self.plotpath + "' did not exist. It has been created!")

    # ------------------------
    def ReadTree(self):
        # self.rootFile = TFile(self.fileName)
        self.simTree = TChain("EDepSimEvents")
        self.simTree.Add(self.fileName)
        self.nEntry = self.simTree.GetEntries()
        # Only Genie has gRooTracker, Marley doesn't
        if self.evgen == 'Genie':
            self.genieTree = TChain("DetSimPassThru/gRooTracker")
            self.genieTree.Add(self.fileName)
            # self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")
            if self.nEntry != self.genieTree.GetEntries():
                print("Edep-sim tree and GENIE tree number of entries do not match!")
                sys.exit()
        elif self.evgen == 'Marley':
            print("Marley events: Skip looking for gRooTracker directory.")
        else:
            print("Unknown event generator!")
            sys.exit()

        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)

    # ------------------------
    def Jump(self, entryNo):
        if entryNo%100 == 0:
            print(f'reading event {entryNo}/{self.nEntry}')

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)

        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')
        
        self.info = {}
        self.info['Event_ID'] = self.currentEntry
        self.info['E_nu'] = 0
        self.info['E_avail'] = 0
        self.info['E_availList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoTotal'] = 0
        self.info['E_depoTotal_track'] = 0
        self.info['Q_depoTotal_dots_th_75keV'] = 0
        self.info['Q_depoTotal'] = 0
        self.info['Q_depoTotal_th_75keV'] = 0
        self.info['Q_depoTotal_th_500keV'] = 0
        self.info['L_depoTotal_avg_APEX_WP'] = 0
        self.info['L_depoTotal_avg_APEX_WP2'] = 0
        self.info['E_depoList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoList_track'] = np.zeros(8)
        self.info['Q_depoList']               = np.zeros(8)
        self.info['Q_depoList_th_75keV']      = np.zeros(8)
        self.info['Q_depoList_th_500keV']     = np.zeros(8)
        self.info['Q_depoList_dots_th_75keV'] = np.zeros(8)
        self.info['L_depoList_avg_APEX_WP'] = np.zeros(8)
        self.info['L_depoList_avg_APEX_WP2'] = np.zeros(8)
        self.info['N_parList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['nu_pdg'] = 0
        self.info['nu_xs'] = self.vertex.GetCrossSection()
        self.info['nu_proc'], self.info['nu_nucl'] = self.GetReaction()
        self.FillEnergyInfo()
        if self.evgen == 'Genie':
            self.ReadGenie()
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

    # ------------------------
    def ReadGenie(self):
        # gRooTracker info
        # https://github.com/GENIE-MC/Generator/blob/master/src/Apps/gNtpConv.cxx#L1837
        # StdHepStatus: 0: initial state; 1: final state particles; others: intermediate transport
        # following assumes the first particle is always the neutrino.
        self.genieTree.GetEntry(self.currentEntry)
        self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
        self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000

    # ------------------------
    def ReadMarley(self):
        # we are mostly looking at nue anyway
        self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'nue': self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'numu': self.info['nu_pdg'] = 14
        if self.GetnuPDGFromFileName() == 'anue': self.info['nu_pdg'] = -12
        if self.GetnuPDGFromFileName() == 'anumu': self.info['nu_pdg'] = -14
        self.info['E_nu'] = self.GetEnuFromFileName()
        if self.currentEntry == 0: # only print once
            print("Marley events: assert info from file name")

    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        if (primaries.size != 1 and self.evgen == 'Genie'):
            print("Number of primaries not equal to 1 (not neutrino vertex)!")
            return

        self.vertex = primaries[0]

    #--------------------------
    def GetReaction(self):
        txt_list = self.vertex.GetReaction().split(';')
        proc = ''
        nucl = 0
        for x in txt_list[2:]:
            if 'proc:' in x:
                proc = x.replace('proc:', '').replace('Weak[','').replace('],', '')
            elif 'N:' in x:
                nucl = int(x.replace('N:', ''))
        # print(proc, nucl)
        proc_num = 0
        if (proc.startswith('CC')):
            proc_num = 10
            proc = proc.replace('CC', '')
        elif (proc.startswith('NC')):
            proc_num = 20
            proc = proc.replace('NC', '')
        else:
            proc_num = 30
        proc_dict = {
            'QES' : 1,
            'RES' : 2,
            'DIS' : 3,
            'COH' : 4,
            'MEC' : 5,
        }
        proc_num += proc_dict.get(proc, 0)
        return proc_num, nucl

    # ------------------------
    def ReadTracks(self):
        self.tracks = np.array(self.event.Trajectories)
        for i in range(self.tracks.size):
            self.tracks[i].energy = {}
            self.tracks[i].association = {}
            self.tracks[i].length = {}
            self.tracks[i].energy['depoTotal'] = 0
            self.tracks[i].energy['depoTotal_charge'] = 0
            self.tracks[i].energy['depoTotal_charge_th_75keV'] = 0
            self.tracks[i].energy['depoTotal_charge_th_500keV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_APEX_WP'] = 0
            self.tracks[i].energy['depoTotal_light_avg_APEX_WP2'] = 0
            self.tracks[i].length['selfDepo'] = 0
            self.tracks[i].association['depoList'] = []
            self.tracks[i].association['children'] = []
            self.tracks[i].association['ancestor'] = i

        for i in range(self.tracks.size):
            track = self.tracks[i]
            parId = track.GetParentId()
            if parId == -1: continue
            self.tracks[parId].association['children'].append(i)
            while parId != -1:
                track.association['ancestor'] = parId
                parId = self.tracks[parId].GetParentId()

    #Section For Light Yield sampling ####################################################
    def get_light_yield_at_sampled_position(self):
        # Open the 2D light yield map
        file = TFile('apex_ly_map.root')
        lymap = file.Get('ly_map')
        
        if not lymap:
            raise ValueError("2D LY map not found in apex_ly_map.root")

        # Get number of bins
        x_bins = lymap.GetNbinsX()
        y_bins = lymap.GetNbinsY()

        # Randomly select a bin within the LY map
        bin_x = np.random.randint(1, x_bins + 1)
        bin_y = np.random.randint(1, y_bins + 1)

        # Convert bin index to real-world coordinates
        true_x_pos = lymap.GetXaxis().GetBinCenter(bin_x)
        true_y_pos = lymap.GetYaxis().GetBinCenter(bin_y)

        # Ensure positions stay within bounds
        true_x_pos = np.clip(true_x_pos, -6.0, 6.0)
        true_y_pos = np.clip(true_y_pos, -6.5, 6.5)

        # Uncertainties used (m)
        sigma_x = 0.28
        sigma_y = 0.47

        # Sample position using Gaussian
        sampled_x = np.random.normal(loc=true_x_pos, scale=sigma_x)                    
        sampled_y = np.random.normal(loc=true_y_pos, scale=sigma_y)

        # Get light yield at the sampled position
        bin_x_truePos = lymap.GetXaxis().FindBin(true_x_pos)
        bin_y_truePos = lymap.GetYaxis().FindBin(true_y_pos)
        
        # Get light yield at the sampled position with uncertainties
        bin_x_truePosUncert = lymap.GetXaxis().FindBin(sampled_x)
        bin_y_truePosUncert = lymap.GetYaxis().FindBin(sampled_y)

        # Get light yield values
        ly_true_pos = lymap.GetBinContent(bin_x_truePos, bin_y_truePos)
        ly_true_pos_uncert = lymap.GetBinContent(bin_x_truePosUncert, bin_y_truePosUncert)

        return ly_true_pos, ly_true_pos_uncert
    
    # ------------------------
    def ReadEnergyDepo(self, detName):
        self.depos = np.array(self.event.SegmentDetectors[detName])

        # add depo info to tracks
        # depoList = self.FindDepoListFromTrack(1)
        # depoEnergy = np.sum([depo.GetEnergyDeposit() for depo in self.depos[depoList]])
        # print('debug: muon deposit energy:', depoEnergy)
        
        # Get light yields
        ly_sampled1, ly_sampled2 = self.get_light_yield_at_sampled_position()
        # Ensure neither ly_sampled1 nor ly_sampled2 is zero
        while ly_sampled1 == 0 or ly_sampled2 == 0:
            ly_sampled1, ly_sampled2 = self.get_light_yield_at_sampled_position()

        mm2cm = 0.1
        # Proceed with calculations
        for i, depo in enumerate(self.depos):
            trkId = depo.Contrib[0]
            edep = depo.GetEnergyDeposit()
            trkLength = depo.GetTrackLength() * mm2cm
            Qdep = self.ChargeBirksLaw(edep, trkLength)
            track = self.tracks[trkId]
            track.association['depoList'].append(i)
            track.energy['depoTotal'] += edep
            track.energy['depoTotal_charge'] += Qdep
            track.length['selfDepo'] += trkLength
            if Qdep > 0.075: track.energy['depoTotal_charge_th_75keV'] += Qdep
            if Qdep > 0.5: track.energy['depoTotal_charge_th_500keV'] += Qdep

            # Light deposit
            # edepSecond = depo.GetSecondaryDeposit()
            # Ideally this should be the light dep
            # how the light yield is actually derived: 25k photons / MeV are simualted and 220 on average are detected
            # the overall photon collection efficiency (PCE) is  PCE = light_yield_number / 21622
            # (dL/W_ph)*PCE = (dL/19.5eV)*PCE = number of photons (taking into account Birks model)
            Ldep = edep - Qdep
            PCE1 = ly_sampled1 / 21622.0
            PCE2 = ly_sampled2 / 21622.0   
            nPE_sampled = (Ldep * 1000000 / 19.5) * PCE1
            nPE_detected = random.gauss(nPE_sampled, math.sqrt(nPE_sampled))
            Ldep_avg_detected = nPE_detected * 19.5 / 1000000 / PCE1
            Ldep_avg_detected2 = nPE_detected * 19.5 / 1000000 / PCE2
            track.energy['depoTotal_light_avg_APEX_WP'] += Ldep_avg_detected
            track.energy['depoTotal_light_avg_APEX_WP2'] += Ldep_avg_detected2

 

    # ------------------------
    def ChargeBirksLaw(self, edep, trkLength):
        # PHYS. REV. D 99, 036009 (2019)
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = A/(1+kQ*dE/dx)
        # A = 0.8, kQ = 0.0972 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83*edep*0.8/(1+0.0972*edep/trkLength/1.4)

    # ------------------------
    def ChargeModifiedBoxModel(self, edep, trkLength):
        # https://lar.bnl.gov/properties/pass.html
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = ln(A+B*dE/dx)/(B*dE/dx)
        # A = 0.930, B = 0.424 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83*edep*np.log(0.930+0.424*edep/trkLength/1.4)/(0.424*edep/trkLength/1.4)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo)
             in enumerate(self.depos)
             if (trkId in depo.Contrib)]
            # if (trkId == depo.GetPrimaryId())]
        return x

    # ------------------------
    def PrintVertex(self):
        if self.evgen == 'Genie':
            self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
            self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000
            print(f"neutrino {self.info['nu_pdg']}: {self.info['E_nu']} MeV")
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

        posx = self.vertex.GetPosition().X()
        posy = self.vertex.GetPosition().Y()
        posz = self.vertex.GetPosition().Z()
        print(f"vertex @: ({posx}, {posy}, {posz}) [mm]")
        print(f"reaction: {self.vertex.GetReaction()}")
        # print(f"interaction #: {self.vertex.GetInteractionNumber()}")

        # print(f"{self.vertex.Particles.size()} particles at the vertex", )
        print(f'{"pdg":>8}{"name":>8}{"trkId":>6}{"mass":>10}{"KE":>10}')
        print(f'{"":>8}{"":>8}{"":>6}{"[MeV]":>10}{"[MeV]":>10}')
        print('-'*(8+8+6+10+10))
        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            print(f'{pdg:>8d}{name:>8s}{trkId:>6d}{mass:>10.2f}{KE:>10.2f}')
        print('-'*(8+8+6+10+10))

        print(f'{self.info}')

    # ------------------------
    def FillEnergyInfo(self):

        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            depoE = self.GetEnergyDepoWithDesendents(trkId)
            track = self.tracks[trkId]
            # if longer than 2cm, assume can reconstruct dE
            if track.length['selfDepo'] > 2:
                depoE_track = track.energy['depoTotal']
                depoQ_dots = 0
            else:
                # dots/blips
                depoE_track = 0
                depoQ_dots = track.energy['depoTotal_charge_th_75keV']
            depoQ        = self.GetChargeDepoWithDesendents(trkId)[0]
            depoQ_75keV  = self.GetChargeDepoWithDesendents(trkId)[1]
            depoQ_500keV = self.GetChargeDepoWithDesendents(trkId)[2]
            depoLY = self.GetLightDepoWithDesendents(trkId)
            depoLY2 = self.GetLightDepoWithDesendents2(trkId)
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            self.info['E_depoTotal'] += depoE
            self.info['E_depoTotal_track'] += depoE_track
            self.info['Q_depoTotal'] += depoQ
            self.info['Q_depoTotal_th_75keV'] += depoQ_75keV
            self.info['Q_depoTotal_th_500keV'] += depoQ_500keV
            self.info['Q_depoTotal_dots_th_75keV'] += depoQ_dots
            self.info['L_depoTotal_avg_APEX_WP'] += depoLY
            self.info['L_depoTotal_avg_APEX_WP2'] += depoLY2
            # fill E_availList: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
            if (pdg in [13, -13, 11, -11]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][0] += (KE + mass)
                self.info['E_depoList'][0]       += depoE
                self.info['E_depoList_track'][0] += depoE_track
                self.info['Q_depoList'][0]                += depoQ
                self.info['Q_depoList_th_75keV'][0]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][0]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][0]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][0] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][0] += depoLY2
                self.info['N_parList'][0] += 1
            elif (pdg == 2212):
                self.info['E_avail'] += KE
                self.info['E_availList'][1] += KE
                self.info['E_depoList'][1]       += depoE
                self.info['E_depoList_track'][1] += depoE_track
                self.info['Q_depoList'][1]                += depoQ
                self.info['Q_depoList_th_75keV'][1]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][1]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][1]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][1] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][1] += depoLY2
                self.info['N_parList'][1] += 1
            elif (pdg == 2112):
                self.info['E_avail'] += KE
                self.info['E_availList'][2] += KE
                self.info['E_depoList'][2]       += depoE
                self.info['E_depoList_track'][2] += depoE_track
                self.info['Q_depoList'][2]                += depoQ
                self.info['Q_depoList_th_75keV'][2]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][2]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][2]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][2] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][2] += depoLY2
                self.info['N_parList'][2] += 1
            elif (pdg in [211, -211]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][3] += (KE + mass)
                self.info['E_depoList'][3]       += depoE
                self.info['E_depoList_track'][3] += depoE_track
                self.info['Q_depoList'][3]                += depoQ
                self.info['Q_depoList_th_75keV'][3]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][3]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][3]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][3] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][3] += depoLY2
                self.info['N_parList'][3] += 1
            elif (pdg == 111):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][4] += (KE + mass)
                self.info['E_depoList'][4]       += depoE
                self.info['E_depoList_track'][4] += depoE_track
                self.info['Q_depoList'][4]                += depoQ
                self.info['Q_depoList_th_75keV'][4]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][4]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][4]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][4] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][4] += depoLY2
                self.info['N_parList'][4] += 1
            elif (pdg == 22):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][5] += (KE + mass)
                self.info['E_depoList'][5]       += depoE
                self.info['E_depoList_track'][5] += depoE_track
                self.info['Q_depoList'][5]                += depoQ
                self.info['Q_depoList_th_75keV'][5]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][5]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][5]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][5] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][5] += depoLY2
                self.info['N_parList'][5] += 1
            elif (pdg == 1000020040):
                self.info['E_avail'] += KE
                self.info['E_availList'][6] += KE
                self.info['E_depoList'][6]       += depoE
                self.info['E_depoList_track'][6] += depoE_track
                self.info['Q_depoList'][6]                += depoQ
                self.info['Q_depoList_th_75keV'][6]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][6]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][6]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][6] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][6] += depoLY2
                self.info['N_parList'][6] += 1
            else:
                self.info['E_avail'] += KE
                self.info['E_availList'][7] += KE
                self.info['E_depoList'][7]       += depoE
                self.info['E_depoList_track'][7] += depoE_track
                self.info['Q_depoList'][7]                += depoQ
                self.info['Q_depoList_th_75keV'][7]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][7]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][7]  += depoQ_dots
                self.info['L_depoList_avg_APEX_WP'][7] += depoLY
                self.info['L_depoList_avg_APEX_WP2'][7] += depoLY2
                self.info['N_parList'][7] += 1

    # ------------------------
    def PrintDepo(self, i):
        # print(f"{self.depos.size} depo points stored in total")
        depo = self.depos[i]
        contrib = depo.GetContributors()
        primaryId = depo.GetPrimaryId()
        edep = depo.GetEnergyDeposit()
        # edepSecond = depo.GetSecondaryDeposit()
        trkLength = depo.GetTrackLength()
        print(contrib, "%6d %.2f %.2f" %(primaryId, edep, trkLength))

    # ------------------------
    def PrintTrack(self, trkId):
        self.PrintTracks(trkId, trkId+1)
        track = self.tracks[trkId]
        children = track.association['children']
        print(f'children: {children}')
        for childId in children:
            child = self.tracks[childId]
            name = child.GetName()
            mom = child.GetInitialMomentum()
            KE = mom.E() - mom.M()
            print(f"{childId} {name}: {KE:.2f} MeV, {child.energy['depoTotal']:.2f} MeV, {child.association['children']}")

        print(f"{track.Points.size()} points stored in track {trkId}")
        for point in track.Points:
            x = point.GetPosition().X()
            y = point.GetPosition().Y()
            z = point.GetPosition().Z()
            t = point.GetPosition().T()
            print(f"{point.GetProcess()}, {point.GetSubprocess()}, {x}, {y}, {z}, {t}")

    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        # print(f"{self.tracks.size} trajectories stored", )
        print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        print('-'*(8+8+6+6+6+10+10+10))

        for track in self.tracks[start:stop]:
            pdg = track.GetPDGCode()
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            ancestor = track.association['ancestor']
            selfDepo = track.energy['depoTotal']
            allDepo = self.GetEnergyDepoWithDesendents(trkId)
            print(f"{pdg:>8d}{name:>8s}{trkId:>6d}{parId:>6d}{ancestor:>6d}{KE:>10.2f}{selfDepo:>10.2f}{allDepo:>10.2f}")

        print('-'*(8+8+6+6+6+10+10+10))


    # ------------------------
    def Next(self):
        if self.currentEntry != self.nEntry -1:
            self.currentEntry += 1
        else:
            self.currentEntry = 0
        self.Jump(self.currentEntry)

    # ------------------------
    def Prev(self):
        if self.currentEntry != 0:
            self.currentEntry -= 1
        else:
            self.currentEntry = self.nEntry -1
        self.Jump(self.currentEntry)

    #-------------------------
    def GetEnergyDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents(childId)
        return energy

    #-------------------------
    def GetChargeDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        charge = track.energy['depoTotal_charge']
        charge_75keV = track.energy['depoTotal_charge_th_75keV']
        charge_500keV = track.energy['depoTotal_charge_th_500keV']
        children = track.association['children']
        for childId in children:
            charge += self.GetChargeDepoWithDesendents(childId)[0]
            charge_75keV += self.GetChargeDepoWithDesendents(childId)[1]
            charge_500keV += self.GetChargeDepoWithDesendents(childId)[2]
        return [charge, charge_75keV, charge_500keV]

    #-------------------------
    def GetLightDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        light_avg_APEX_WP = track.energy['depoTotal_light_avg_APEX_WP']
        children = track.association['children']
        for childId in children:
            light_avg_APEX_WP += self.GetLightDepoWithDesendents(childId)
        
        return light_avg_APEX_WP

    def GetLightDepoWithDesendents2(self, trkId):
        track = self.tracks[trkId]
        light_avg_APEX_WP2 = track.energy['depoTotal_light_avg_APEX_WP2']
        children = track.association['children']
        for childId in children:
            light_avg_APEX_WP2 += self.GetLightDepoWithDesendents2(childId)
        
        return light_avg_APEX_WP2
    #-----------------------
    def GetEnuFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('MeV', '')
        return float(filename)  # Marley file names already in MeV

    #-----------------------
    def GetnuPDGFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-3]
        return filename

    #-----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()
# ------------------------
if __name__ == "__main__":
    event = Event(sys.argv[1])
    #event.Jump(0)
    #event.PrintVertex()
    #event.Next()
    #event.PrintVertex()
    #event.PrintTracks(0,8)
    # event.PrintTrack(1)
