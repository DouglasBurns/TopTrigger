'''
Created on 08 May 2015

@author: Douglas Burns

Calculate Trigger Efficiency and also Added Trigger Efficiency

'''
from __future__ import division
import ROOT 
from ROOT import gROOT, gPad, gStyle, TFile, TMath, TEfficiency, TGraph, TCollection, TH1F, TF1, TCanvas, TLegend, TGraphPainter, TPaveStats, kRed, kBlue, kGreen
import math

if __name__ == '__main__':


	########## 			SETUP 			##########
	gStyle.SetOptStat("")

	TriggerType = {
	0 : 'SingleTop',
	1 : 'TTBarJet30',
	2 : 'TTBarJet304050',
	};

	# inputFile = TFile('../../StoreTrigger/Trigger_MC_50ns.root')
	# Type = 'TTBar_MC'

	inputFile = TFile('data/Electron.root')
	Type = 'Electron_Data'

	# inputFile = TFile('data/Muon.root')
	# Type = 'Muon_Data'

	# inputFile = TFile('data/Data.root')
	# inputFile = TFile('../../StoreTrigger/Trigger_Data.root')
	# Type = 'E+Mu_Data'

	file = TFile("Turn-on.root", "RECREATE")

	if (Type == 'Electron_Data'):
		VariableType = {
		0 : 'Ele27',
		1 : 'Ele32',
		}
	if (Type == 'Muon_Data'):
		VariableType = {
		0 : 'Mu20', 
		1 : 'Mu24',
		}

	if (Type == 'E+Mu_Data'): 
		VariableType = {
		0 : 'Ele27',
		1 : 'Ele32',
		2 : 'Mu20',
		3 : 'Mu24',
		};

	if (Type == 'TTBar_MC'): 
		VariableType = {
		0 : 'Ele27',
		1 : 'Ele32',
		2 : 'Mu20',
		3 : 'Mu24',
		};

	print Type

	# gStyle.SetOptFit(0111)
	# f1 = TF1("f1","(0.5*[0]*(1+TMath::Erf((x-[1])/(sqrt(x)*[2]))))",0,200)
	# f2 = TF1("f2","(0.5*[0]*(1+TMath::Erf((x-[1])/(sqrt(x)*[2]))))",0,10)
	# f3 = TF1("f2","(0.5*[0]*(1+TMath::Erf((x-[1])/(sqrt(x)*[2]))))",0,1)

	for TriggerIndex in range (0,len(TriggerType)):
		for VariableIndex in range (0, len(VariableType)):

			diffEffJetMult20Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/CrossTrigger_Pass_JetMultiplicity_20')
			JetMult20Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/CrossTrigger_Total_JetMultiplicity_20')

			diffEffHTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/CrossTrigger_Pass_HT')
			HTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/CrossTrigger_Total_HT')
			
			diffEffVertexMultHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Vertices/CrossTrigger_Pass_VertexMultiplicityHist')
			VertexMultHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Vertices/CrossTrigger_Total_VertexMultiplicityHist')
			
			diffEffMETHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/MET/CrossTrigger_Pass_MET')
			METHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/MET/CrossTrigger_Total_MET')
			
			diffEffLeadingLepHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Leading_Lepton/CrossTrigger_Pass_LeptonPt')
			LeadingLepHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Leading_Lepton/CrossTrigger_Total_LeptonPt')


			diffEffJetMult20Hist.Draw()
			JetMult20Hist.Draw()
			diffEffHTHist.Draw()
			HTHist.Draw()
			diffEffVertexMultHist.Draw()
			VertexMultHist.Draw()
			diffEffMETHist.Draw()
			METHist.Draw()
			diffEffLeadingLepHist.Draw()
			LeadingLepHist.Draw()

			JetMult20 = TEfficiency(diffEffJetMult20Hist,JetMult20Hist)
			HT = TEfficiency(diffEffHTHist,HTHist)
			VertexMult = TEfficiency(diffEffVertexMultHist,VertexMultHist)
			MET = TEfficiency(diffEffMETHist,METHist)
			LeadingLep = TEfficiency(diffEffLeadingLepHist,LeadingLepHist)


			CanvasJetMult20 = TCanvas(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_" + "JetMult20","Title", 0, 0, 800, 600)
			CanvasJetMult20.SetGrid()
			JetMult20.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Differential_Efficiency_JetMult20; RECO Jet Pt (GeV); Efficiency" )
			JetMult20.SetMarkerColor(4)
			JetMult20.SetMarkerStyle(21)
			JetMult20.Draw()
			JetMult20.Write()
			gPad.Update()
			JetMult20.Draw("ALP")
			gPad.Update()
			CanvasJetMult20.Update()
			CanvasJetMult20.SaveAs("plots/" + Type + '_' + VariableType[VariableIndex] + '_' +  TriggerType[TriggerIndex] + "_" + "Differential_Efficiency_JetMult20.png" )

			CanvasHT = TCanvas(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_" + "HT","Title", 0, 0, 800, 600)
			CanvasHT.SetGrid()
			HT.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Differential_Efficiency_HT; Number of Jets > 20 GeV; Efficiency" )
			HT.SetMarkerColor(4)
			HT.SetMarkerStyle(21)
			HT.Draw()
			HT.Write()
			gPad.Update()
			HT.Draw("ALP")
			gPad.Update()
			CanvasHT.Update()
			CanvasHT.SaveAs("plots/" + Type + '_' + VariableType[VariableIndex] + '_' +  TriggerType[TriggerIndex] + "_" + "Differential_Efficiency_HT.png" )

			CanvasVertexMult = TCanvas(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_" + "VertexMult","Title", 0, 0, 800, 600)
			CanvasVertexMult.SetGrid()
			VertexMult.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Differential_Efficiency_VertexMult; Number of Vertices; Efficiency" )
			VertexMult.SetMarkerColor(4)
			VertexMult.SetMarkerStyle(21)
			VertexMult.Draw()
			VertexMult.Write()
			gPad.Update()
			VertexMult.Draw("ALP")
			gPad.Update()
			CanvasVertexMult.Update()
			CanvasVertexMult.SaveAs("plots/" + Type + '_' + VariableType[VariableIndex] + '_' +  TriggerType[TriggerIndex] + "_" + "Differential_Efficiency_VertexMult.png" )

			CanvasMET = TCanvas(Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_" + "MET","Title", 0, 0, 800, 600)
			CanvasMET.SetGrid()
			MET.SetTitle(VariableType[VariableIndex]+" Differential_Efficiency_MET; MET (GeV); Efficiency" )
			MET.SetMarkerColor(4)
			MET.SetMarkerStyle(21)
			MET.Draw()
			MET.Write()
			gPad.Update()
			MET.Draw("ALP")
			gPad.Update()
			CanvasMET.Update()
			CanvasMET.SaveAs("plots/" + Type + '_' + VariableType[VariableIndex] + '_' +  TriggerType[TriggerIndex] + "_" + "Differential_Efficiency_MET.png" )

			CanvasLeadingLeptonPt = TCanvas(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_" + "LeadLepPt","Title", 0, 0, 800, 600)
			CanvasLeadingLeptonPt.SetGrid()
			LeadingLep.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Differential_Efficiency_LeadingLeptonPt; RECO Leading Lepton Pt (GeV); Efficiency" )
			LeadingLep.SetMarkerColor(4)
			LeadingLep.SetMarkerStyle(21)
			LeadingLep.Draw()
			LeadingLep.Write()
			gPad.Update()
			LeadingLep.Draw("ALP")
			gPad.Update()
			CanvasLeadingLeptonPt.Update()
			CanvasLeadingLeptonPt.SaveAs("plots/" + Type + '_' + VariableType[VariableIndex] + '_' +  TriggerType[TriggerIndex] + "_" + "Differential_Efficiency_LeadingLeptonPt.png" )



			if (TriggerType[TriggerIndex] == 'SingleTop'):

				if (VariableType[VariableIndex] == 'Ele27' or VariableType[VariableIndex] == 'Ele32'): filteredpTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + VariableType[VariableIndex] + '/matched RECO Jet Observables/Filter 1 matched Jet Pt' )
				if (VariableType[VariableIndex] == 'Mu20' or VariableType[VariableIndex] == 'Mu24'): filteredpTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltJetFilter' + TriggerType[TriggerIndex] + 'Iso' + VariableType[VariableIndex] + 'Eta2p1/matched RECO Jet Observables/Filter 1 matched Jet Pt' )



				pTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/Total_RECO_Jet_Pt')
				filteredCSVHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltCSVFilterSingleTop/matched RECO Jet Observables/Filter 2 matched Jet CSV -log10(1-CSV)' )
				CSVHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/Total_RECO_Jet_BTag')
				filteredBTagHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltCSVFilterSingleTop/matched RECO Jet Observables/Filter 2 matched Jet CSV')
				BTagHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/Total_RECO_Jet_CSV')
				# filteredCSVHist.SetAxisRange(0,4,"X")
				# CSVHist.SetAxisRange(0,4,"X")
				# filteredCSVHist.GetXaxis().SetLimits(0,3.5)
				# CSVHist.GetXaxis().SetLimits(0,3.5)

				filteredpTHist.Draw()
				pTHist.Draw()
				filteredCSVHist.Draw()
				CSVHist.Draw()
				filteredBTagHist.Draw()
				BTagHist.Draw()				

				Filter1 = TEfficiency(filteredpTHist,pTHist)
				Filter2 = TEfficiency(filteredCSVHist,CSVHist)
				Filter3 = TEfficiency(filteredBTagHist,BTagHist)


				CanvasCSV = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex]+ '_CSV',"Title", 0, 0, 800, 600)
				CanvasCSV.SetLogy()
				filteredCSVHist.SetLineColor(kRed)
				CSVHist.Draw()
				filteredCSVHist.Draw("same")
				Leg_CSV = TLegend(0.2,0.8,0.4,0.88)#xmin,ymin,xmax,ymax (in % of canvas)
				Leg_CSV.SetFillColor(0)
				Leg_CSV.SetLineColor(0)
				Leg_CSV.AddEntry(filteredCSVHist, "Filtered CSV" ,"f")
				Leg_CSV.AddEntry(CSVHist, "Inclusive CSV" ,"f")
				Leg_CSV.Draw()
				CanvasCSV.Update()
				CanvasCSV.SaveAs('plots/' + Type + '_'  + VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter1_CSV_Histos.png" )

				CanvasBTag = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex]+ '_BTag',"Title", 0, 0, 800, 600)
				CanvasBTag.SetLogy()
				filteredBTagHist.SetLineColor(kRed)
				BTagHist.Draw()
				filteredBTagHist.Draw("same")
				Leg_BTag = TLegend(0.2,0.8,0.4,0.88)#xmin,ymin,xmax,ymax (in % of canvas)
				Leg_BTag.SetFillColor(0)
				Leg_BTag.SetLineColor(0)
				Leg_BTag.AddEntry(filteredBTagHist, "Filtered BTag" ,"f")
				Leg_BTag.AddEntry(BTagHist, "Inclusive BTag" ,"f")
				Leg_BTag.Draw()
				CanvasBTag.Update()
				CanvasBTag.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter1_BTag_Histos.png" )


				filteredpTHist.Write()
				pTHist.Write()
				filteredCSVHist.Write()
				CSVHist.Write()
				filteredBTagHist.Write()
				BTagHist.Write()

				########## Diff Canvases ##########

				Canvas = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex],"Title", 0, 0, 800, 600)
				Canvas.SetGrid()

				Filter1.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Single Top Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency" )
				Filter1.SetMarkerColor(4)
				Filter1.SetMarkerStyle(21)
				Filter1.Draw()
				Filter1.Write()
				gPad.Update()

				# Filter1.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f1.SetParameters(1.0,30,0.6)
				# Filter1.Fit(f1)
				Filter1.Draw("ALP")
				gPad.Update()

				# stats1 = Filter1.GetPaintedGraph().FindObject("stats")
				# stats1.SetX1NDC(0.5);
				# stats1.SetX2NDC(0.9);
				# stats1.SetY1NDC(0.1);
				# stats1.SetY2NDC(0.3);
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter1_TurnOnPt.png" )



				Filter2.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Single Top Turn On Jet CSV; RECO Jet CSV (pfCSVv2); Efficiency" )
				Filter2.SetMarkerColor(4)
				Filter2.SetMarkerStyle(21)
				Filter2.Draw()
				Filter2.Write()
				gPad.Update()

				# Filter2.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f2.SetParameters(1.0,2,0.3)
				# Filter2.Fit(f2)
				Filter2.Draw("ALP")
				gPad.Update()
				
				# stats2 = Filter2.GetPaintedGraph().FindObject("stats")
				# stats2.SetX1NDC(0.5);
				# stats2.SetX2NDC(0.9);
				# stats2.SetY1NDC(0.1);
				# stats2.SetY2NDC(0.3);
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter2_TurnOnCSV.png" )	

				Filter3.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Single Top Turn On Jet BTag; RECO Jet BTag (pfBTagv2); Efficiency" )
				Filter3.SetMarkerColor(4)
				Filter3.SetMarkerStyle(21)
				Filter3.Draw()
				Filter3.Write()
				gPad.Update()

				# Filter3.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f3.SetParameters(1.0,0.9,0.3)
				# Filter3.Fit(f3)
				Filter3.Draw("ALP")
				gPad.Update()
				
				# stats3 = Filter3.GetPaintedGraph().FindObject("stats")
				# stats3.SetX1NDC(0.5);
				# stats3.SetX2NDC(0.9);
				# stats3.SetY1NDC(0.1);
				# stats3.SetY2NDC(0.3);
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter2_TurnOnBTag.png" )	

			if (TriggerType[TriggerIndex] == 'TTBarJet30'):

				if (VariableType[VariableIndex] == 'Ele27' or VariableType[VariableIndex] == 'Ele32'): filteredpTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'TriCentralPFJet30EleCleaned/matched RECO Jet Observables/Filter 1 matched Jet Pt' )
				if (VariableType[VariableIndex] == 'Mu20' or VariableType[VariableIndex] == 'Mu24'):  filteredpTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02TriCentralPFJet30MuCleaned/matched RECO Jet Observables/Filter 1 matched Jet Pt' )
				pTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/Total_RECO_Jet_Pt')

				Filter1 = TEfficiency(filteredpTHist,pTHist)

				filteredpTHist.Draw()
				pTHist.Draw()

				filteredpTHist.Write()
				pTHist.Write()

				######### Diff Canvases ##########

				Canvas = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex],"Title", 0, 0, 800, 600)
				Canvas.SetGrid()

				Filter1.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Symmetric TTBar Turn On Jet Pt (30 GeV); RECO Jet Pt (GeV); Efficiency" )
				Filter1.SetMarkerColor(4)
				Filter1.SetMarkerStyle(21)
				Filter1.Draw()
				Filter1.Write()
				gPad.Update()

				# Filter1.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f1.SetParameters(1.0,30,0.6)
				# Filter1.Fit(f1)
				Filter1.Draw("ALP")
				gPad.Update()
				
				# stats1 = Filter1.GetPaintedGraph().FindObject("stats")
				# stats1.SetX1NDC(0.5)
				# stats1.SetX2NDC(0.9)
				# stats1.SetY1NDC(0.1)
				# stats1.SetY2NDC(0.3)				
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter1_TurnOnPt.png" )



			if (TriggerType[TriggerIndex] == 'TTBarJet304050'):
				if (VariableType[VariableIndex] == 'Ele27' or VariableType[VariableIndex] == 'Ele32'): 
					filteredpT30Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'TriCentralPFJet30EleCleaned/matched RECO Jet Observables/Filter 1 matched Jet Pt' )
					filteredpT40Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'DiCentralPFJet40EleCleaned/matched RECO Jet Observables/Filter 2 matched Jet Pt' )
					filteredpT50Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'CentralPFJet50EleCleaned/matched RECO Jet Observables/Filter 3 matched Jet Pt' )
					pT30DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'TriCentralPFJet30EleCleaned/Turn On Curves/Filter1_dR' )
					pT40DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'DiCentralPFJet40EleCleaned/Turn On Curves/Filter2_dR' )
					pT50DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hlt' + VariableType[VariableIndex] + 'CentralPFJet50EleCleaned/Turn On Curves/Filter3_dR' )

				if (VariableType[VariableIndex] == 'Mu20' or VariableType[VariableIndex] == 'Mu24'):  
					filteredpT30Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02TriCentralPFJet30MuCleaned/matched RECO Jet Observables/Filter 1 matched Jet Pt' )
					filteredpT40Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02DiCentralPFJet40MuCleaned/matched RECO Jet Observables/Filter 2 matched Jet Pt' )
					filteredpT50Hist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02CentralPFJet50MuCleaned/matched RECO Jet Observables/Filter 3 matched Jet Pt' )
					pT30DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02TriCentralPFJet30MuCleaned/Turn On Curves/Filter1_dR' )
					pT40DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02DiCentralPFJet40MuCleaned/Turn On Curves/Filter2_dR' )
					pT50DeltaRHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/hltIso' + VariableType[VariableIndex] + 'Eta2p1Trk02CentralPFJet50MuCleaned/Turn On Curves/Filter3_dR' )




				pTHist = inputFile.Get( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '/Trigger Observables/Jets/Total_RECO_Jet_Pt')

				# Filter1 = TEfficiency(filteredpT30Hist,pTHist)
				# Filter2 = TEfficiency(filteredpT40Hist,pTHist)
				# Filter3 = TEfficiency(filteredpT50Hist,pTHist)

				Filter1 = TEfficiency(filteredpT30Hist,pTHist)
				Filter2 = TEfficiency(filteredpT40Hist,filteredpT30Hist)
				Filter3 = TEfficiency(filteredpT50Hist,filteredpT40Hist)

				filteredpT30Hist.Draw()
				filteredpT40Hist.Draw()
				filteredpT50Hist.Draw()
				pTHist.Draw()

				filteredpT30Hist.Write()
				filteredpT40Hist.Write()
				filteredpT50Hist.Write()
				pTHist.Write()

				########## Plotting Scripts. ##########

				dRCanvas = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '_dR' ,"Title", 0, 0, 800, 600)
				dRCanvas.SetGrid()	
				dRCanvas.SetLogy()

				pT30DeltaRHist.SetTitle(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_DeltaR Distributions; Delta R; Number of Events")
				pT30DeltaRHist.SetLineColor(kRed)
				pT30DeltaRHist.SetMarkerColor(4)
				pT30DeltaRHist.SetMarkerStyle(21)
				pT30DeltaRHist.Draw()

				pT40DeltaRHist.SetLineColor(kBlue)
				pT40DeltaRHist.SetMarkerColor(4)
				pT40DeltaRHist.SetMarkerStyle(21)
				pT40DeltaRHist.Draw("same")

				pT50DeltaRHist.SetLineColor(kGreen)
				pT50DeltaRHist.SetMarkerColor(4)
				pT50DeltaRHist.SetMarkerStyle(21)
				pT50DeltaRHist.Draw("same")

				dRLeg = TLegend(0.5,0.5,0.88,0.88)		
				dRLeg.SetFillColor(0)
				dRLeg.SetLineColor(0)
				dRLeg.AddEntry(pT30DeltaRHist, ">30" ,"l")
				dRLeg.AddEntry(pT40DeltaRHist, ">40" ,"l")
				dRLeg.AddEntry(pT50DeltaRHist, ">50" ,"l")
				dRLeg.Draw()

				gPad.Update()
				dRCanvas.Update()
				dRCanvas.Write()
				dRCanvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_dR.png" )



				pTCanvas = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + '_pT' ,"Title", 0, 0, 800, 600)
				pTCanvas.SetGrid()	
				pTCanvas.SetLogy()

				pTHist.SetTitle(VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Pt Distributions; Pt (GeV); Number of Events")
				pTHist.SetMarkerColor(4)
				pTHist.SetMarkerStyle(21)
				pTHist.Draw()

				filteredpT30Hist.SetLineColor(kRed)
				filteredpT30Hist.SetMarkerColor(4)
				filteredpT30Hist.SetMarkerStyle(21)
				filteredpT30Hist.Draw("same")

				filteredpT40Hist.SetLineColor(kBlue)
				filteredpT40Hist.SetMarkerColor(4)
				filteredpT40Hist.SetMarkerStyle(21)
				filteredpT40Hist.Draw("same")

				filteredpT50Hist.SetLineColor(kGreen)
				filteredpT50Hist.SetMarkerColor(4)
				filteredpT50Hist.SetMarkerStyle(21)
				filteredpT50Hist.Draw("same")

				pTLeg = TLegend(0.5,0.5,0.88,0.88)		
				pTLeg.SetFillColor(0)
				pTLeg.SetLineColor(0)
				pTLeg.AddEntry(pTHist, "Central Cleaned Jets" ,"l")
				pTLeg.AddEntry(filteredpT30Hist, ">30" ,"l")
				pTLeg.AddEntry(filteredpT40Hist, ">40" ,"l")
				pTLeg.AddEntry(filteredpT50Hist, ">50" ,"l")
				pTLeg.Draw()

				gPad.Update()
				pTCanvas.Update()
				pTCanvas.Write()
				pTCanvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_pT.png" )



				Canvas = TCanvas( VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex],"Title", 0, 0, 800, 600)
				Canvas.SetGrid()

				Filter1.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Assymmetric TTBar Turn On Jet Pt (30 GeV); RECO Jet Pt (GeV); Efficiency" )
				Filter1.SetMarkerColor(4)
				Filter1.SetMarkerStyle(21)
				Filter1.Draw()
				Filter1.Write()
				gPad.Update()

				# Filter1.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f1.SetParameters(1.0,30,0.6)
				# Filter1.Fit(f1)
				Filter1.Draw("ALP")
				gPad.Update()
				
				# stats1 = Filter1.GetPaintedGraph().FindObject("stats")
				# stats1.SetX1NDC(0.5)
				# stats1.SetX2NDC(0.9)
				# stats1.SetY1NDC(0.1)
				# stats1.SetY2NDC(0.3)				
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter1_TurnOnPt.png" )


				Filter2.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Assymmetric TTBar Turn On Jet Pt (40 GeV); RECO Jet Pt (GeV); Efficiency" )
				Filter2.SetMarkerColor(4)
				Filter2.SetMarkerStyle(21)
				Filter2.Draw()
				Filter2.Write()
				gPad.Update()

				# Filter2.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f1.SetParameters(9.9,40,0.6)
				# Filter2.Fit(f1)
				Filter2.Draw("ALP")
				gPad.Update()
				
				# stats2 = Filter2.GetPaintedGraph().FindObject("stats")
				# stats2.SetX1NDC(0.5)
				# stats2.SetX2NDC(0.9)
				# stats2.SetY1NDC(0.1)
				# stats2.SetY2NDC(0.3)				
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter2_TurnOnPt.png" )
				
				Filter3.SetTitle(Type + '_'  +  VariableType[VariableIndex]+" Assymmetric TTBar Turn On Jet Pt (50 GeV); RECO Jet Pt (GeV); Efficiency" )
				Filter3.SetMarkerColor(4)
				Filter3.SetMarkerStyle(21)
				Filter3.Draw()
				Filter3.Write()
				gPad.Update()

				# Filter3.GetPaintedGraph().SetMaximum(1.1)
				gPad.Update()

				# f1.SetParameters(0.9,50,0.6)
				# Filter3.Fit(f1)
				Filter3.Draw("ALP")
				gPad.Update()
				
				# stats3 = Filter3.GetPaintedGraph().FindObject("stats")
				# stats3.SetX1NDC(0.5)
				# stats3.SetX2NDC(0.9)
				# stats3.SetY1NDC(0.1)
				# stats3.SetY2NDC(0.3)				
				Canvas.Update()
				Canvas.SaveAs('plots/' + Type + '_'  +  VariableType[VariableIndex] + '_' + TriggerType[TriggerIndex] + "_Filter3_TurnOnPt.png" )

