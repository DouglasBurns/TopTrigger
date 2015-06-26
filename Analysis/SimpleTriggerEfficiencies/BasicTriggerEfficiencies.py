'''
Created on 08 May 2015

@author: Douglas Burns

Calculate Trigger Efficiency and also Added Trigger Efficiency

'''
from __future__ import division
import ROOT 
from ROOT import gROOT, gPad, gStyle, TFile, TMath, TH1F
import math

if __name__ == '__main__':


	########## 			SETUP 			##########
	gStyle.SetOptStat("")

	VariableType = {
	0 : 'Ele27_',
	1 : 'Ele32_',
	2 : 'Mu20_',
	3 : 'Mu24_'
	};

	TriggerType = {
	0 : 'SingleTop',
	1 : 'TTBarJet30',
	2 : 'TTBarJet304050'
	};

	# input file
	inputFile = TFile('../../data/TestTrigger.root')
	for TriggerIndex in range (0,len(TriggerType)):
		for VariableIndex in range (0, len(VariableType)):
			# print 'Getting Leptonic Hist: ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/' + TriggerType[TriggerIndex] + ' Trigger Decision'
			# print 'Getting Hadronic Hist: ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Single Lepton Trigger Decision'
			# print 'Getting Combined Hist: ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Added ' + TriggerType[TriggerIndex] + ' Trigger Decision'
			hadronicHist = inputFile.Get( VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/' + TriggerType[TriggerIndex] + ' Trigger Decision')
			leptonicHist = inputFile.Get( VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Single Lepton Trigger Decision')
			combinedHist = inputFile.Get( VariableType[VariableIndex] + TriggerType[TriggerIndex] + '/Trigger Decision/Added ' + TriggerType[TriggerIndex] + ' Trigger Decision')

			hadronicHist_Eff = hadronicHist.GetBinContent(2)/(hadronicHist.GetBinContent(1)+hadronicHist.GetBinContent(2))
			leptonicHist_Eff = leptonicHist.GetBinContent(2)/(leptonicHist.GetBinContent(1)+leptonicHist.GetBinContent(2))
			combinedHist_Eff = combinedHist.GetBinContent(2)/(combinedHist.GetBinContent(1)+combinedHist.GetBinContent(2))

			print "***********************************************************************"
			print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Single Lepton :' ), leptonicHist_Eff
			print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Hadronic :' ), hadronicHist_Eff
			print ('Efficiency of ' + VariableType[VariableIndex] + TriggerType[TriggerIndex] + ', Hadronic | Single Lepton :' ), combinedHist_Eff

