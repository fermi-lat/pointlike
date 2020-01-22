"""
A suite of tools for making easy things easy with ROOT.

Author: Damien Parent <dmnparent@gmail.com>
"""
from ROOT import gROOT, gStyle, TH1F, TH2F, TGraph, gPad, TGaxis, Double, TPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
import numpy as np
from numpy import array as sa

gROOT.SetBatch(True) # force the batch mode

# default font for text, label, ...
default_font = 83

def initialization(batch=True, font=default_font):
    '''-- ROOT initialization --'''
    
    print ("Initializing ROOT ...")
    
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

def map_color(color): 
    if color is 'black': kcolor = kBlack
    elif color is 'white': kcolor = kWhite
    elif color is 'red': kcolor = kRed+1
    #elif color is 'blue': kcolor = kBlue+1
    elif color is 'blue': kcolor = kBlue
    elif color is 'green': kcolor = kGreen+2
    elif color is 'gray': kcolor = kGray+2
    elif color is 'orange': kcolor = kOrange+3
    elif color is 'yellow': kcolor = kYellow
    else: print ("Warning: color %s is not implemented!" % color ); kcolor = kBlack
    return kcolor

def SetHistoAxis( hist, x_title="", x_title_size=0, x_title_offset=0, x_label_size=0, 
                  y_title="", y_title_size = 0, y_title_offset=0, y_label_size=0,
                  font=default_font, color='black', line_width=1 ):
    
    kcolor = map_color(color) 

    hist.SetLineWidth(line_width)
    hist.SetLineColor(kcolor)

    hist.GetXaxis().SetTitle(x_title)
    hist.GetXaxis().CenterTitle()
    hist.GetXaxis().SetTitleSize(x_title_size)
    hist.GetXaxis().SetTitleOffset(x_title_offset)
    hist.GetXaxis().SetTitleFont(font)
    hist.GetXaxis().SetLabelSize(x_label_size)
    hist.GetXaxis().SetLabelOffset(0.)
    hist.GetXaxis().SetLabelFont(font)
    hist.GetXaxis().SetNdivisions(510)
    
    hist.GetYaxis().SetTitle(y_title)
    hist.GetYaxis().CenterTitle()
    hist.GetYaxis().SetTitleSize(y_title_size)
    hist.GetYaxis().SetTitleOffset(y_title_offset)
    hist.GetYaxis().SetTitleFont(font)
    hist.GetYaxis().SetLabelSize(y_label_size)
    hist.GetYaxis().SetLabelOffset(3e-3)
    hist.GetYaxis().SetLabelFont(font)
    hist.GetYaxis().SetNdivisions(505)

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

def histo_extrema(histo):
    """ Find minimum of histogram, including error."""
    y,ye = get_histo_content(histo)
    #if len(y)==0: return 0,0
    m = y > 0
    if m.sum()==0: return 0,0
    minim = ((y-ye)[m]).min()
    maxim = ((y+ye)[m]).max()
    return minim,maxim

def zero_suppress(histo,force_ymin=None,background=None):
    """Zero suppression for a histogram.
    
        background [None] -- if not none, specify the background level,
        and make sure the minimum is below this level"""
    if force_ymin is None:
        ymin = histo_extrema(histo)[0]
        if (background is not None):
            ymin = min(background,ymin)
        if ymin>0:
            #ymin = int(ymin*0.9 - np.sqrt(ymin))
            #ymin = int(ymin*0.9)
            ymin = int(ymin)
            #if ymin%5 == 0: ymin *= 0.9
            if ymin%5 == 0: ymin -= 1
            if ymin <= 0: ymin = 0.01 
        else: ymin = 0.01
    else: ymin = force_ymin
    histo.SetMinimum(ymin)

def scale_radio_profile(tgraph,tmin,tmax,bkg_bin=0.1):
    bkg_bin = 0.15*tgraph.GetN()
    x = np.empty(tgraph.GetN());
    y = np.empty_like(x)
    for i in xrange(len(x)):
        xt=Double();yt=Double()
        tgraph.GetPoint(i,xt,yt)
        x[i] = float(xt)
        y[i] = float(yt)
    y -= np.sort(y)[bkg_bin]
    y *= (tmax-tmin)/y.max()
    y += tmin
    for i in xrange(len(x)):
        tgraph.SetPoint(i,x[i],y[i])

def get_tgraph_content(tgraph):
    y = np.empty(tgraph.GetN())
    x = np.empty(tgraph.GetN())
    for i in xrange(len(y)):
        xt=Double();yt=Double()
        tgraph.GetPoint(i,xt,yt)
        y[i] = float(yt)
        x[i] = float(xt)
    return x,y

def tgraph_min_max(tgraph):
    x,y = get_tgraph_content(tgraph)
    return y.min(),y.max()

def get_histo_content(histo,first_half=False):
    n = histo.GetNbinsX()
    if first_half: n /= 2
    y = np.empty(n)
    yerr = np.empty(n)
    for i in xrange(n):
        y[i] = histo.GetBinContent(i)
        yerr[i] = histo.GetBinError(i)
    return y,yerr
    
def get_txtlevel(histo, factor):
    min = histo.GetMinimum()
    max = histo.GetMaximum()
    return (max-min) * factor + min
                               
def eraselabel(pad, h):
    pad.cd()
    pe = TPad("pe","pe",0,0,pad.GetLeftMargin(),h)
    pe.Draw()
    pe.SetFillColor(pad.GetFillColor())
    pe.SetBorderMode(0)

def BPalette():
    r=sa([0.,0.0,1.0,1.0,1.0])
    b=sa([0., 1.0, 0.0, 0.0, 1.0])
    g=sa([0., 0.0, 0.0, 1.0, 1.0])
    stop=sa([0.,.25,.50,.75,1.0])
    TColor.CreateGradientColorTable(5,stop,r,g,b,100)
    return

def GrayPalette():
    R=sa([0.,1.])
    G=sa([0.,1.])
    B=sa([0.,1.])
    Stop=sa([0.,1.])
    TColor.CreateGradientColorTable(2,Stop,R,G,B,100)
    return
