"""
A suite of tools for making easy things easy with ROOT.

Author: Damien Parent <dmnparent@gmail.com>
"""

import numpy as np
from ROOT import gROOT
# force the batch mode
gROOT.SetBatch(True)
from ROOT import gStyle    

FONT = 42

def initialization(batch=True):
    '''ROOT initialization.

    parameters:
    -----------
    batch : Run ROOT in batch mode [True/False]
    '''
    print "Initializing ROOT ..."

    # general
    gROOT.Reset()
    gROOT.SetBatch(batch)
    gROOT.SetStyle("Plain")

    # gStyle
    gStyle.SetFillColor(0)
    gStyle.SetCanvasColor(10)
    gStyle.SetLineWidth(1)
    gStyle.SetPalette(1)
    gStyle.SetTextFont(FONT)

    # Frame
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameFillColor(0)

    # Pad
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetPadTopMargin(0.)
    gStyle.SetPadLeftMargin(0.10)
    gStyle.SetPadRightMargin(0.02)
    gStyle.SetPadBottomMargin(0.)
    gStyle.SetPadTickX(1)    # make ticks be on all 4 sides.
    gStyle.SetPadTickY(1)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    
    # histogram
    gStyle.SetHistFillStyle(0)
    gStyle.SetOptTitle(0)
    
    # histogram title
    gStyle.SetTitleSize(0.22)
    gStyle.SetTitleFontSize(2)
    gStyle.SetTitleFont(FONT)
    gStyle.SetTitleFont(FONT,"xyz")
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
    gStyle.SetStatFont(FONT)
    gStyle.SetStatX(.91)
    gStyle.SetStatY(.90)
    gStyle.SetStatW(.15)
    gStyle.SetStatH(.15)
        
    # axis labels
    gStyle.SetLabelFont(FONT)
    gStyle.SetLabelFont(FONT,"xyz")
    gStyle.SetLabelSize(0.035,"xyz")
    # gStyle.SetGridColor(1)
    gStyle.SetLegendBorderSize(1);


def set_axis_thx( hist, x_title = "", x_title_size = 0., x_title_offset = 0., x_label_size = 0., 
                  y_title = "", y_title_size = 0., y_title_offset = 0., y_label_size = 0. ):
    
    hist.GetXaxis().SetTitle(x_title)
    hist.GetXaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(x_title_size)
    hist.GetXaxis().SetTitleOffset(x_title_offset)
    hist.GetXaxis().SetTitleFont(FONT)
    hist.GetXaxis().SetLabelSize(x_label_size)
    hist.GetXaxis().SetLabelOffset(0.01)
    hist.GetXaxis().SetLabelFont(FONT)
    
    hist.GetYaxis().SetTitle(y_title)
    hist.GetYaxis().CenterTitle()
    hist.GetYaxis().SetTitleSize(y_title_size)
    hist.GetYaxis().SetTitleOffset(y_title_offset)
    hist.GetYaxis().SetTitleFont(FONT)
    hist.GetYaxis().SetLabelSize(y_label_size)
    hist.GetYaxis().SetLabelOffset(0.01)
    hist.GetYaxis().SetLabelFont(FONT)

def fit_lorentzian(x, par):
    return (0.5*par[0]*par[1]/np.pi) / np.max( 1.e-10,(x[0]-par[2]) * (x[0]-par[2]) + .25*par[1]*par[1])

