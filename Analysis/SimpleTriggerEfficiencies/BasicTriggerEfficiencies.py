'''
Created on 08 May 2015

@author: Douglas Burns

Calculate Trigger Efficiency and also Added Trigger Efficiency

'''
from __future__ import division
import ROOT 
from ROOT import gROOT, gPad, gStyle, TChain, TFile, TTree, TMath, TH1, TH1F, TH2F, TCanvas, TPad, TAxis, TLegend, TLatex, kRed, kBlue, kGreen
import math

if __name__ == '__main__':


	########## 			SETUP 			##########
	gStyle.SetOptStat("")
	input_file = "../../data/TestTrigger.root"
	inputTree = "storeTrigger/Triggers"
	Chain = TChain(inputTree)
	Chain.Add(input_file)
	Chain.SetBranchStatus("*",1)

	NumberEventsTotal = 0
	Ele27 = Ele32 = Mu20 = Mu24 = 0
	Ele27BTrig = Ele32BTrig = Mu20BTrig = Mu24BTrig = 0
	Ele27TriJet30Trig = Ele32TriJet30Trig = Mu20TriJet30Trig = Mu24TriJet30Trig = 0
	Ele27TriJet304050Trig = Ele32TriJet304050Trig = Mu20TriJet304050Trig = Mu24TriJet304050Trig = 0

	AddedEle27BTrig = AddedEle32BTrig = AddedMu20BTrig = AddedMu24BTrig = 0
	AddedEle27TriJet30Trig = AddedEle32TriJet30Trig = AddedMu20TriJet30Trig = AddedMu24TriJet30Trig = 0
	AddedEle27TriJet304050Trig = AddedEle32TriJet304050Trig = AddedMu20TriJet304050Trig = AddedMu24TriJet304050Trig = 0

	ODDMu20BTrig = ODDAddedMu20BTrig = ODDMu20TriJet30Trig = ODDAddedMu20TriJet30Trig = ODDMu20TriJet304050Trig = ODDAddedMu20TriJet304050Trig = 0

	for event in Chain:

		NumberEventsTotal += 1
		Trigger_Decision = event.__getattr__("Trigger_Decision")
		# Number_of_Triggers = len(Trigger_Decision)
		# print ("Number of Triggers : "), len(Trigger_Decision)

		# for index in range (0, Number_of_Triggers):
			# PassTrigger = Trigger_Decision.at(index)
			# print ("Pass Trigger : "), PassTrigger

		if (Trigger_Decision[0] == 1): Ele27BTrig += 1
		if (Trigger_Decision[1] == 1): Ele27TriJet30Trig += 1
		if (Trigger_Decision[2] == 1): Ele27TriJet304050Trig += 1
		
		if (Trigger_Decision[3] == 1): 
			Ele27 += 1
			if (Trigger_Decision[0] == 1): AddedEle27BTrig += 1 		
			if (Trigger_Decision[1] == 1): AddedEle27TriJet30Trig += 1 
			if (Trigger_Decision[2] == 1): AddedEle27TriJet304050Trig += 1

		if (Trigger_Decision[4] == 1): Ele32BTrig += 1
		if (Trigger_Decision[5] == 1): Ele32TriJet30Trig += 1
		if (Trigger_Decision[6] == 1): Ele32TriJet304050Trig += 1
		
		if (Trigger_Decision[7] == 1): 
			Ele32 += 1
			if (Trigger_Decision[4] == 1): AddedEle32BTrig += 1 		
			if (Trigger_Decision[5] == 1): AddedEle32TriJet30Trig += 1 
			if (Trigger_Decision[6] == 1): AddedEle32TriJet304050Trig += 1 	
		
		if (Trigger_Decision[8] == 1): Mu20BTrig += 1
		if (Trigger_Decision[9] == 1): Mu20TriJet30Trig += 1
		if (Trigger_Decision[10] == 1): Mu20TriJet304050Trig += 1
 		
		if (Trigger_Decision[11] == 1):  
			Mu20 += 1		
	 		if (Trigger_Decision[8] == 1): AddedMu20BTrig += 1 		
			if (Trigger_Decision[9] == 1): AddedMu20TriJet30Trig += 1 
			if (Trigger_Decision[10] == 1): AddedMu20TriJet304050Trig += 1 	

		if (Trigger_Decision[12] == 1): Mu24BTrig += 1
		if (Trigger_Decision[13] == 1): Mu24TriJet30Trig += 1
		if (Trigger_Decision[14] == 1): Mu24TriJet304050Trig += 1
		
		if (Trigger_Decision[15] == 1): 
			Mu24 +=1
			if (Trigger_Decision[12] == 1): AddedMu24BTrig += 1 		
			if (Trigger_Decision[13] == 1): AddedMu24TriJet30Trig += 1 
			if (Trigger_Decision[14] == 1): AddedMu24TriJet304050Trig += 1 

# Weird Happenings #
		if (Trigger_Decision[8] == 1):
			ODDMu20BTrig += 1
			if (Trigger_Decision[11] == 0):
				ODDAddedMu20BTrig += 1 		
		if (Trigger_Decision[9] == 1):
			ODDMu20TriJet30Trig += 1
			if (Trigger_Decision[11] == 0):
				ODDAddedMu20TriJet30Trig += 1 
		if (Trigger_Decision[10] == 1):
			ODDMu20TriJet304050Trig += 1
			if (Trigger_Decision[11] == 0):
				ODDAddedMu20TriJet304050Trig += 1 	

	print ("---------------------------------------------------------------------------")
	print ("Number of Event Total : "), NumberEventsTotal
	print ("---------------------------------------------------------------------------")

	print ("Odd Number Mu20B : "), ODDMu20BTrig
	print ("Odd Number Mu20J : "), ODDMu20TriJet30Trig
	print ("Odd Number Mu20J3 : "), ODDMu20TriJet304050Trig
	print ("Odd Number AMu20B : "), ODDAddedMu20BTrig
	print ("Odd Number AMu20J : "), ODDAddedMu20TriJet30Trig
	print ("Odd Number AMu20J3 : "), ODDAddedMu20TriJet304050Trig

	print ("HLT_Ele27_eta2p1_WP75_Gsf__v1 Total : "), Ele27
	print ("HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), Ele27BTrig
	print ("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 Total : "), Ele27TriJet30Trig
	print ("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 Total : "), Ele27TriJet304050Trig
	print ("HLT_Ele32_eta2p1_WP75_Gsf__v1 Total : "), Ele32
	print ("HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), Ele32BTrig
	print ("HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 Total : "), Ele32TriJet30Trig
	print ("HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 Total : "), Ele32TriJet304050Trig

	print ("HLT_IsoMu20_eta2p1_WP75_Gsf__v1 Total : "), Mu20
	print ("HLT_IsoMu20_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), Mu20BTrig
	print ("HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 : "), Mu20TriJet30Trig
	print ("HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 : "), Mu20TriJet304050Trig
	print ("HLT_IsoMu24_eta2p1_WP75_Gsf__v1 Total : "), Mu24
	print ("HLT_IsoMu24_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), Mu24BTrig
	print ("HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 : "), Mu24TriJet30Trig
	print ("HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 : "), Mu24TriJet304050Trig
	print ("---------------------------------------------------------------------------")
	print ("Added HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), AddedEle27BTrig
	print ("Added HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 Total : "), AddedEle27TriJet30Trig
	print ("Added HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 Total : "), AddedEle27TriJet304050Trig
	print ("Added HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), AddedEle32BTrig
	print ("Added HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 Total : "), AddedEle32TriJet30Trig
	print ("Added HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 Total : "), AddedEle32TriJet304050Trig

	print ("Added HLT_IsoMu20_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), AddedMu20BTrig
	print ("Added HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 : "), AddedMu20TriJet30Trig
	print ("Added HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 : "), AddedMu20TriJet304050Trig
	print ("Added HLT_IsoMu24_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 Total : "), AddedMu24BTrig
	print ("Added HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 : "), AddedMu24TriJet30Trig
	print ("Added HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 : "), AddedMu24TriJet304050Trig
	print ("---------------------------------------------------------------------------")
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 wrt HLT_Ele27_eta2p1_WP75_Gsf_v1"), AddedEle27BTrig/Ele27
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"), Ele27BTrig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 wrt HLT_Ele27_eta2p1_WP75_Gsf_v1"), AddedEle27TriJet30Trig/Ele27
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"), Ele27TriJet30Trig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 wrt HLT_Ele27_eta2p1_WP75_Gsf_v1"), AddedEle27TriJet304050Trig/Ele27
	print ("Efficiency of HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1"), Ele27TriJet30Trig/NumberEventsTotal

	print ("*")
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 wrt HLT_Ele32_eta2p1_WP75_Gsf_v1"), AddedEle32BTrig/Ele32
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"), Ele32BTrig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 wrt HLT_Ele32_eta2p1_WP75_Gsf_v1"), AddedEle32TriJet30Trig/Ele32
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"), Ele32TriJet30Trig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 wrt HLT_Ele32_eta2p1_WP75_Gsf_v1"), AddedEle32TriJet304050Trig/Ele32
	print ("Efficiency of HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1"), Ele32TriJet304050Trig/NumberEventsTotal

	print ("*")
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 wrt HLT_IsoMu20_eta2p1_WP75_Gsf_v1"), AddedMu20BTrig/Mu20
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"), Mu20BTrig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 wrt HLT_IsoMu20_eta2p1_WP75_Gsf_v1"), AddedMu20TriJet30Trig/Mu20
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"), Mu20TriJet30Trig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 wrt HLT_IsoMu20_eta2p1_WP75_Gsf_v1"), AddedMu20TriJet304050Trig/Mu20
	print ("Efficiency of HLT_IsoMu20_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1"), Mu20TriJet304050Trig/NumberEventsTotal

	print ("*")
	print ("Efficiency of HLT_IsoMu24_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1 wrt HLT_IsoMu24_eta2p1_WP75_Gsf_v1"), AddedMu24BTrig/Mu24
	print ("Efficiency of HLT_IsoMu24_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"), Mu24BTrig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet30_v1 wrt HLT_IsoMu24_eta2p1_WP75_Gsf_v1"), AddedMu24TriJet30Trig/Mu24
	print ("Efficiency of HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"), Mu24TriJet30Trig/NumberEventsTotal
	print ("*")
	print ("Efficiency of HLT_IsoMu24_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1 wrt HLT_IsoMu24_eta2p1_WP75_Gsf_v1"), AddedMu24TriJet304050Trig/Mu24















