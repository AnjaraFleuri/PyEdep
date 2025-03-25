from event import Event
import sys
import numpy as np
from array import array
from ROOT import TFile, TTree

class Writer:

    def __init__(self, event, outfile='output.root'):
        self.event = event

        self.f_out = TFile(outfile, 'RECREATE')
        self.initOutputTree()

    def initOutputTree(self):
        self.T_out = TTree('Sim', 'Sim') # output tree

        self.Event_ID = array('i', [0]) #
        self.T_out.Branch('Event_ID', self.Event_ID, 'Event_ID/I')

        self.nu_pdg = array('i', [0]) #
        self.T_out.Branch('nu_pdg', self.nu_pdg, 'nu_pdg/I')

        self.nu_xs = array('f', [0])
        self.T_out.Branch('nu_xs', self.nu_xs, 'nu_xs/F')

        self.nu_proc = array('i', [0]) # genie process
        self.T_out.Branch('nu_proc', self.nu_proc, 'nu_proc/I')

        self.nu_nucl = array('i', [0]) # which nucleon
        self.T_out.Branch('nu_nucl', self.nu_nucl, 'nu_nucl/I')

        self.E_nu = array('f', [0]) # true neutrino energy
        self.T_out.Branch('E_nu', self.E_nu, 'E_nu/F')

        self.E_avail = array('f', [0]) # energy availabe from vertex interaction (excluding energy lost inside nuclei)
        self.T_out.Branch('E_avail', self.E_avail, 'E_avail/F')

        self.E_availList = np.zeros((8,), dtype=np.float32) # E avail for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('E_availList', self.E_availList, 'E_availList[8]/F')

        self.E_depoTotal = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal', self.E_depoTotal, 'E_depoTotal/F')

        self.E_depoList = np.zeros((8,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('E_depoList', self.E_depoList, 'E_depoList[8]/F')

        self.E_depoTotal_track = array('f', [0]) # total energy deposit from long tracks
        self.T_out.Branch('E_depoTotal_track', self.E_depoTotal_track, 'E_depoTotal_track/F')

        self.E_depoList_track = np.zeros((8,), dtype=np.float32) # depo for long tracks: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('E_depoList_track', self.E_depoList_track, 'E_depoList_track[8]/F')

        self.Q_depoTotal = array('f', [0]) # total charge deposit from all (charged) tracks
        self.T_out.Branch('Q_depoTotal', self.Q_depoTotal, 'Q_depoTotal/F')

        self.Q_depoList = np.zeros((8,), dtype=np.float32) # charge depo for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('Q_depoList', self.Q_depoList, 'Q_depoList[8]/F')

        self.Q_depoTotal_th_75keV = array('f', [0]) # total charge deposit from all (charged) tracks with dQ cut of 75keV
        self.T_out.Branch('Q_depoTotal_th_75keV', self.Q_depoTotal_th_75keV, 'Q_depoTotal_th_75keV/F')

        self.Q_depoList_th_75keV = np.zeros((8,), dtype=np.float32) # charge depo with 75keV cut for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('Q_depoList_th_75keV', self.Q_depoList_th_75keV, 'Q_depoList_th_75keV[8]/F')

        self.Q_depoTotal_th_500keV = array('f', [0]) # total charge deposit from all (charged) tracks with dQ cut of 500keV
        self.T_out.Branch('Q_depoTotal_th_500keV', self.Q_depoTotal_th_500keV, 'Q_depoTotal_th_500keV/F')

        self.Q_depoList_th_500keV = np.zeros((8,), dtype=np.float32) # charge depo with 500keV cut for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('Q_depoList_th_500keV', self.Q_depoList_th_500keV, 'Q_depoList_th_500keV[8]/F')

        self.Q_depoTotal_dots_th_75keV = array('f', [0]) # total charge deposit from all short dots/blips with dQ cut of 75keV
        self.T_out.Branch('Q_depoTotal_dots_th_75keV', self.Q_depoTotal_dots_th_75keV, 'Q_depoTotal_dots_th_75keV/F')

        self.Q_depoList_dots_th_75keV = np.zeros((8,), dtype=np.float32) # charge depo with 75keV cut for short dots/blips from: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('Q_depoList_dots_th_75keV', self.Q_depoList_dots_th_75keV, 'Q_depoList_dots_th_75keV[8]/F')

        self.L_depoTotal_avg_APEX_WP = array('f', [0]) 
        self.T_out.Branch('L_depoTotal_avg_APEX_WP', self.L_depoTotal_avg_APEX_WP, 'L_depoTotal_avg_APEX_WP/F')

        self.L_depoTotal_avg_APEX_WP2 = array('f', [0]) 
        self.T_out.Branch('L_depoTotal_avg_APEX_WP2', self.L_depoTotal_avg_APEX_WP2, 'L_depoTotal_avg_APEX_WP2/F')

        self.L_depoList_avg_APEX_WP = np.zeros((8,), dtype=np.float32) 
        self.T_out.Branch('L_depoList_avg_APEX_WP', self.L_depoList_avg_APEX_WP, 'L_depoList_avg_APEX_WP[8]]/F')

        self.L_depoList_avg_APEX_WP2 = np.zeros((8,), dtype=np.float32) 
        self.T_out.Branch('L_depoList_avg_APEX_WP2', self.L_depoList_avg_APEX_WP2, 'L_depoList_avg_APEX_WP2[8]]/F')

        self.N_parList = np.zeros((8,), dtype=np.int32) # number of particles for: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.T_out.Branch('N_parList', self.N_parList, 'N_parList[8]/I')

    def Write(self):
        
        # self.stat = {
        # }
        self.f_out.cd()

        for i in range(self.event.nEntry):
        # for i in range(100):
            self.event.Jump(i)

            # proc = str(self.event.info['nu_proc']) + '-' + str(self.event.info['nu_nucl'])
            # v = self.stat.setdefault(proc, 0)
            # self.stat[proc] = v + 
            
            self.Event_ID[0] = self.event.info['Event_ID']
            self.nu_pdg[0] = self.event.info['nu_pdg']
            self.nu_xs[0] = self.event.info['nu_xs']
            self.nu_proc[0] = self.event.info['nu_proc']
            self.nu_nucl[0] = self.event.info['nu_nucl']
            self.E_nu[0] = self.event.info['E_nu']

            self.E_avail[0] = self.event.info['E_avail']
            self.E_depoTotal[0] = self.event.info['E_depoTotal']
            self.E_depoTotal_track[0] = self.event.info['E_depoTotal_track']
            self.Q_depoTotal[0] = self.event.info['Q_depoTotal']
            self.Q_depoTotal_th_75keV[0] = self.event.info['Q_depoTotal_th_75keV']
            self.Q_depoTotal_th_500keV[0] = self.event.info['Q_depoTotal_th_500keV']
            self.Q_depoTotal_dots_th_75keV[0] = self.event.info['Q_depoTotal_dots_th_75keV']
            self.L_depoTotal_avg_APEX_WP[0] = self.event.info['L_depoTotal_avg_APEX_WP']
            self.L_depoTotal_avg_APEX_WP2[0] = self.event.info['L_depoTotal_avg_APEX_WP2']

            self.E_availList[:] = self.event.info['E_availList']
            self.E_depoList[:] = self.event.info['E_depoList']
            self.E_depoList_track[:] = self.event.info['E_depoList_track']
            self.Q_depoList[:] = self.event.info['Q_depoList']
            self.Q_depoList_th_75keV[:] = self.event.info['Q_depoList_th_75keV']
            self.Q_depoList_th_500keV[:] = self.event.info['Q_depoList_th_500keV']
            self.Q_depoList_dots_th_75keV[:] = self.event.info['Q_depoList_dots_th_75keV']
            self.L_depoList_avg_APEX_WP[:] = self.event.info['L_depoList_avg_APEX_WP']
            self.L_depoList_avg_APEX_WP2[:] = self.event.info['L_depoList_avg_APEX_WP2']

            self.N_parList[:] = self.event.info['N_parList']

            self.T_out.Fill()

        self.f_out.cd()
        self.T_out.Write()
        # print(self.stat)

if __name__ == "__main__":
    if (len(sys.argv)>3):
        outfile = sys.argv[3]
    else:
        outfile = 'output.root'
    event = Event(sys.argv[1], sys.argv[2])
    w = Writer(event, outfile)
    w.Write()
