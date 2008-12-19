import numpy as N
import ROOT as R
import pointlike as pl 
import pyfits
from math import floor,exp,log,sqrt,pow

def plot(filename,source,model):

    print filename
    print type(model)

    global g1
    g1=R.TGraphAsymmErrors()
    g1.SetMarkerStyle(7)
    g1.SetMarkerColor(2)

    global g2
    g2=R.TGraphAsymmErrors()
    g2.SetMarkerStyle(7)
    g2.SetMarkerColor(2)

    file=pyfits.open(filename)

    name=""

    for src in range(0,len(file["SOURCES"].data.field('NAME'))):
        name=file["SOURCES"].data.field('NAME')[src]
        if name==source:
            print "Found "+name 
            emin=file["SOURCES"].data.field('EMIN')[src]
            emax=file["SOURCES"].data.field('EMAX')[src]
            alpha=file["SOURCES"].data.field('ALPHA')[src]
            err_alpha=file["SOURCES"].data.field('ERR_ALPHA')[src]
            energy=file["SOURCES"].data.field(model.get_spec_type()+"_ENERGY")[src]
            model_exposure=file["SOURCES"].data.field(model.get_spec_type()+"_EXPOSURE")[src]
            exposure=file["SOURCES"].data.field("EXPOSURE")[src]
            nphoton=file["SOURCES"].data.field('NPHOTON')[src]
            params=file["SOURCES"].data.field(model.get_spec_type()+"_PAR")[src]

    file.close()

    print params

    global func
    func=set_function(model,params)

    exposure_error_fraction=0.01
    livebin=0
    for bin in range(0,int(len(emin)/2)):

        front=2*bin
        back=2*bin+1

        E_min=emin[front]
        E_max=emax[front]
        E_model=energy[livebin]
        exposure_error_fraction=model.get_exposure_uncertainty(R.Double(E_min),R.Double(E_max))
        delta_E=E_max-E_min

        if nphoton[front]>0:
            flux_front=pow(E_model,2)*nphoton[front]*alpha[front]/(model_exposure[livebin]*delta_E)
            err_flux_front=flux_front*sqrt(pow(err_alpha[front]/alpha[front],2)+pow(exposure_error_fraction,2))
        else:
            flux_front=0
            err_flux_front=0
        
        E_model=energy[livebin]

        if nphoton[back]>0:
            flux_back=pow(E_model,2)*nphoton[back]*alpha[back]/(model_exposure[livebin]*delta_E)
            err_flux_back=flux_back*sqrt(pow(err_alpha[back]/alpha[back],2)+pow(exposure_error_fraction,2))
        else:
            flux_back=0
            err_flux_back=0
                
        flux=(flux_front+flux_back)/2
        err_flux=sqrt(pow(err_flux_front,2)+pow(err_flux_back,2))

        # Make sure the energy range was actually fitted
        if exposure[front]!=0 and exposure[back]!=0:
            livebin=livebin+2

        if flux==0:
            continue

        g1.SetPoint(bin,E_model,flux)
        g1.SetPointError(bin,E_model-E_min,E_max-E_model,err_flux,err_flux)

        model_flux=func.Eval(E_model)

        g2.SetPoint(bin,E_model,(flux-model_flux)/flux)
        g2.SetPointError(bin,E_model-E_min,E_max-E_model,err_flux/flux,err_flux/flux)

    g1.Print()
    g2.Print()

    global canvas
    global pad1
    global pad2
    canvas,pad1,pad2=set_canvas(name)
    
    global box
    box=set_box(g1)

    global res_box
    res_box=set_res_box(g2)

    global legend
    legend=set_legend(name,model,g1,func)

    canvas.Draw()
    box.Draw()
    g1.Draw("P")
    func.Draw("SAME")
    legend.Draw()

    pad2.cd()
    res_box.Draw()
    g2.Draw("P")
    global ref
    ref = R.TF1("ref","0",10,5e5)
    ref.SetLineWidth(1)
    ref.SetLineStyle(7)
    ref.SetLineColor(4)
    ref.Draw("SAME")
    
    canvas.Update()

def set_canvas(name):

    R.gStyle.SetOptStat(0)
    R.gStyle.SetOptTitle(0)
    
    c1 = R.TCanvas("SED","SED: "+name,200,10,500,700)
    c1.Draw()

    p1 = R.TPad("SpecPad","SpecPad",0.0,0.3,1.0,1.0);

    p1.SetLogx()
    p1.SetLogy()

    p1.SetLeftMargin(0.14)
    p1.SetBottomMargin(0.1)
    p1.SetTopMargin(0.04)
    p1.SetRightMargin(0.04)
    p1.SetFillColor(10)
    p1.Draw();

    c1.cd()
    p2 = R.TPad("ResPad2","ResPad",0.0,0.0,1.0,0.3)
    p2.Draw()
    p2.cd()

    p2.SetLogx()

    p2.SetLeftMargin(0.14)
    p2.SetBottomMargin(0.14)
    p2.SetTopMargin(0.04)
    p2.SetRightMargin(0.04)
    p2.SetFillColor(10)

    p1.cd();

    return c1,p1,p2

def set_box(graph):
    xmin=R.Double()
    xmax=R.Double()
    ymin=R.Double()
    ymax=R.Double()
    graph.ComputeRange(xmin,ymin,xmax,ymax)
   
    e_min=100
    e_max=3e5
    flux_min=max(1e-8,ymin*R.Double(0.1))
    flux_max=ymax*R.Double(10.)

    h = R.TH2F("Spectral","Spectrum",10,e_min,e_max,10,flux_min,flux_max)
    h.SetTitle("Spectral Energy Distribution")
    h.GetXaxis().SetTitle("Energy [MeV]")
    h.GetXaxis().SetTitleFont(42)
    h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitle("E^{2} dN/dE [MeV cm^{-2} s^{-1}]")
    h.GetYaxis().SetTitleFont(42)
    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitleOffset(1.5)

    return h

def set_res_box(graph):
    xmin=R.Double()
    xmax=R.Double()
    ymin=R.Double()
    ymax=R.Double()
    graph.ComputeRange(xmin,ymin,xmax,ymax)
   
    e_min=100
    e_max=3e5
    res_min=max(min(0,ymin*R.Double(1.2)),-6)
    res_max=min(max(0,ymax*R.Double(1.2)),6)

    h = R.TH2F("RESIDUAL","RESIDUAL",10,e_min,e_max,10,res_min,res_max)
    h.SetTitle("Spectral Energy Distribution Residual")
    h.GetXaxis().SetTitle("Energy [MeV]")
    h.GetXaxis().SetTitleFont(42)
    h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitle("Fit Residuals - \DeltaF/F")
    h.GetYaxis().SetTitleFont(42)
    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitleOffset(1.5)

    return h

# This function is a kludge
# Check consistency of models with SpectralModel.h
def set_function(model,params):
    print type(model)
        
    scale=100. # [MeV]

    if model.get_spec_type()=="POWER_LAW":
        f1 = R.TF1("f1","pow(x,2)*[0]*pow(x/[2],-1*[1])",10,5e5)
        f1.SetParameters(params[0],params[1],scale)
    elif model.get_spec_type()=="BROKEN_POWER_LAW":
        f1 = R.TF1("f1","(x<[3])*pow(x,2)*[0]*pow(x/[3],-1*[1])+(x>[3])*pow(x,2)*[0]*pow(x/[3],-1*[2])",10,5e5);
        f1.SetParameters(params[0],params[1],params[2],params[3]);
    elif model.get_spec_type()=="EXP_CUTOFF":
        f1 = R.TF1("f1","pow(x,2)*[0]*exp(-1.*x/[2])*pow(x/[3],-1*[1])",10,5e5);
        f1.SetParameters(params[0],params[1],params[2],scale)
    else:
        print "Invalid spectral model type"

    f1.SetLineWidth(1)
    f1.SetLineStyle(7)
    f1.SetLineColor(4)

    return f1

def set_legend(name,model,graph,func):
    legend = R.TLegend(0.7,0.8,0.85,0.9)
    legend.SetHeader(name)
    legend.AddEntry(g1,"Data","P")
    legend.AddEntry(func,model.get_spec_type().replace('_',' '),"L")
    legend.SetFillColor(0)
    legend.SetShadowColor(0)
    legend.SetLineColor(0)
    return legend
