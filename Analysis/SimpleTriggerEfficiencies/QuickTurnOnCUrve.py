'''
Created on 08 May 2015

@author: Douglas Burns

Calculate Trigger Efficiency and also Added Trigger Efficiency

'''
from __future__ import division
import ROOT 
from ROOT import gROOT, gPad, gStyle, TFile, TMath, TEfficiency, TGraph, TH1F, TCanvas
import math

if __name__ == '__main__':


	########## 			SETUP 			##########
	gStyle.SetOptStat("")

	VariableType = {
	0 : 'Ele27',
	1 : 'Ele32',
	2 : 'Mu20',
	3 : 'Mu24'
	};

	TriggerType = {
	0 : 'SingleTop',
	1 : 'TTBarJet30',
	2 : 'TTBarJet304050'
	};

	# input file
	inputFile = TFile('../../data/TestTrdscigger.root')
	file = TFile("TurnON.root", "RECREATE")

	for TriggerIndex in range (0,len(TriggerType)):
		for VariableIndex in range (0, len(VariableType)):


			if (TriggerType[TriggerIndex] == 'SingleTop'):
				# if (VariableType[VariableIndex] == 'Ele27' or VariableType[VariableIndex] == 'Ele32'):
				if (VariableType[VariableIndex] == 'Ele27'):

					print 'Getting Filtered Pt Hist: ' + VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + VariableType[VariableIndex] + '/matched RECO Jet Observables/Filtered matched Jet Pt'
					print 'Getting Pt Hist: ' + VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + VariableType[VariableIndex] + '/matched RECO Jet Observables/matched Jet Pt'
					# print 'Getting Combined Hist: ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Added ' + TriggerType[TriggerIndex] + ' Trigger Decision'
					filteredpTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + VariableType[VariableIndex] + '/matched RECO Jet Observables/Filtered matched Jet Pt' )
					pTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + VariableType[VariableIndex] + '/matched RECO Jet Observables/matched Jet Pt')
					# combinedHist = inputFile.Get( VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Added ' + TriggerType[TriggerIndex] + ' Trigger Decision')
					print "Got Hists!"

					Filter1 = TEfficiency(filteredpTHist,pTHist)

					filteredpTHist.Draw()
					pTHist.Draw()
					filteredpTHist.Write()
					pTHist.Write()
					########## ETT Diff Canvases ##########

					Canvas1 = TCanvas("ETTCuts","ETTCuts", 0, 0, 800, 600)

					#TEfficiencyDictionary_ETT[cut].SetLineColor(int(cut/50))
					Filter1.SetTitle("ETT Cut at;Pt;Efficiency" )
					Filter1.Draw()
					Filter1.Write()
					gPad.Update()

					Canvas1.Update()
					Canvas1.SaveAs("TurnOn.png")
					# print "***********************************************************************"
					# print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Single Lepton :' ), leptonicHist_Eff
					# print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Hadronic :' ), hadronicHist_Eff
					# print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Hadronic | Single Lepton :' ), combinedHist_Eff

