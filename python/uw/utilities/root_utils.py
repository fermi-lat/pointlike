"""
A suite of tools for making easy things easy with ROOT.

Author: Damien Parent <dmnparent@gmail.com>
"""

import numpy as np
from ROOT import gROOT, gStyle, TH1F, TH2F, TGraph, gPad, TGaxis, Double
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan

gROOT.SetBatch(True) # force the batch mode

# default font for text, label, ...
default_font = 83

def initialization(batch=True, font=default_font):
    '''-- ROOT initialization --'''
    
    print "Initializing ROOT ..."
    
    # general
    gROOT.Reset()
    gROOT.SetBatch(batch)
    gROOT.SetStyle("Plain")
    
    # gStyle
    gStyle.SetFillColor(0)
    gStyle.SetCanvasColor(10)
    gStyle.SetLineWidth(1)
    gStyle.SetPalette(8)
    gStyle.SetTextFont(font)
    #gStyle.SetTextSize(30)
    
    # Frame
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameFillColor(0)

    # Pad
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetPadBottomMargin(0.1)
    gStyle.SetPadTopMargin(0.01)
    gStyle.SetPadLeftMargin(0.1)
    gStyle.SetPadRightMargin(0.01)
    gStyle.SetPadTickX(1)    # make ticks be on all 4 sides.
    gStyle.SetPadTickY(1)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)

    # histogram
    gStyle.SetHistFillStyle(0)
    gStyle.SetOptTitle(0)
    gStyle.SetTitleSize(0.22)
    gStyle.SetTitleFontSize(10)
    gStyle.SetTitleFont(font)
    gStyle.SetTitleFont(font,"xyz")
    gStyle.SetTitleYOffset(1.0)
    gStyle.SetTitleXOffset(1.0)
    gStyle.SetTitleXSize(0.04)
    gStyle.SetTitleYSize(0.04)
    gStyle.SetTitleX(.15)
    gStyle.SetTitleY(.98)
    gStyle.SetTitleW(.70)
    gStyle.SetTitleH(.05)
    
    # statistics box
    gStyle.SetOptStat(0)
    gStyle.SetStatFont(font)
    gStyle.SetStatFontSize(10)
    gStyle.SetStatX(.91)
    gStyle.SetStatY(.90)
    gStyle.SetStatW(.15)
    gStyle.SetStatH(.15)
        
    # axis labels
    gStyle.SetLabelFont(font)
    gStyle.SetLabelFont(font,"xyz")
    gStyle.SetLabelSize(10,"xyz")
    # gStyle.SetGridColor(1)
    gStyle.SetLegendBorderSize(1);

    
def SetHistoAxis( hist, x_title="", x_title_size=0, x_title_offset=0, x_label_size=0, 
                  y_title="", y_title_size = 0, y_title_offset=0, y_label_size=0,
                  font=default_font, color='black' ):
    
    if color is 'black': kcolor = kBlack
    elif color is 'white': kcolor = kWhite
    elif color is 'red': kcolor = kRed+1
    elif color is 'blue': kcolor = kBlue+1
    elif color is 'green': kcolor = kGreen+2
    elif color is 'gray': kcolor = kGray+2
    elif color is 'orange': kcolor = kOrange+3
    elif color is 'yellow': kcolor = kYellow
    else: print "Warning: color %s is not implemented!"; kcolor = kBlack
        
    # hist.SetLineWidth(1)
    hist.SetLineColor(kcolor)

    hist.GetXaxis().SetTitle(x_title)
    hist.GetXaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(x_title_size)
    hist.GetXaxis().SetTitleOffset(x_title_offset)
    hist.GetXaxis().SetTitleFont(font)
    hist.GetXaxis().SetLabelSize(x_label_size)
    hist.GetXaxis().SetLabelOffset(0.)
    hist.GetXaxis().SetLabelFont(font)
    hist.GetXaxis().SetNdivisions(510);
    
    hist.GetYaxis().SetTitle(y_title)
    hist.GetYaxis().CenterTitle()
    hist.GetYaxis().SetTitleSize(y_title_size)
    hist.GetYaxis().SetTitleOffset(y_title_offset)
    hist.GetYaxis().SetTitleFont(font)
    hist.GetYaxis().SetLabelSize(y_label_size)
    hist.GetYaxis().SetLabelOffset(0.01)
    hist.GetYaxis().SetLabelFont(font)
    hist.GetYaxis().SetNdivisions(505);

def DrawAxis( axis, title="", title_size=0, title_offset=0, label_size=0, font=default_font, div=510 ):
    axis.SetTitle(title)
    axis.CenterTitle()
    axis.SetTitleSize(title_size)
    axis.SetTitleOffset(title_offset)
    axis.SetTitleFont(font)
    axis.SetLabelSize(label_size)
    axis.SetLabelOffset(0.01)
    axis.SetLabelFont(font)
    axis.SetNdivisions(div);
    axis.Draw()

def fit_lorentzian(x, par):
    return (0.5*par[0]*par[1]/np.pi) / np.max( 1.e-10,(x[0]-par[2]) * (x[0]-par[2]) + .25*par[1]*par[1])

def ScaleGraphY(graph, factor):
    '''Scales the y values of a graph with the factor in the argument'''
    for i in range(graph.GetN()):
        x, y = Double(0), Double(0)
        graph.GetPoint(i,x,y)
        y *= factor
        graph.SetPoint(i,x,y)

def zero_suppress(histo):
    '''Zoom on the histogram'''
    min_histo = histo.GetBinContent(histo.GetMinimumBin())
    if min_histo != 0:
        factor = 0.8 - np.sqrt(min_histo)/min_histo
        min_level = int(min_histo*factor)
        if min_level%5 == 0 or min_level%5 == 4: min_level += 3
    else: min_level = 0.01
    histo.SetMinimum(min_level)
    
def get_txtlevel(histo, factor):
    min = histo.GetMinimum()
    max = histo.GetMaximum()
    return (max-min) * factor + min
                               
