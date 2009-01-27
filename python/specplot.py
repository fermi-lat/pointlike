import numpy as N
import ROOT as R
import pointlike as pl 
import pyfits
import sys
from math import floor,exp,log,sqrt,pow

def plot(filename,source=""):

    print filename
    
    global g1
    g1=R.TGraphAsymmErrors()
    g1.SetMarkerStyle(7)
    g1.SetMarkerColor(1)

    global g2
    g2=R.TGraphAsymmErrors()
    g2.SetMarkerStyle(7)
    g2.SetMarkerColor(1)

    global g3
    g3=R.TGraphAsymmErrors()
    g3.SetMarkerStyle(1)
    g3.SetMarkerColor(2)#2
    g3.SetLineColor(2)

    file=pyfits.open(filename)

    name=""
    
    for src in range(0,len(file["SOURCES"].data.field('NAME'))):
        name=file["SOURCES"].data.field('NAME')[src]
        if name==source or len(file["SOURCES"].data.field('NAME'))==1:
            print "Found "+name 

            type=str(file["SOURCES"].data.field('MODEL')[src])
            print type
            if type=="POWER_LAW":
                model=pl.PowerLaw()
            elif type=="BROKEN_POWER_LAW":
                model=pl.BrokenPowerLaw()
            elif type=="EXP_CUTOFF":
                model=pl.ExpCutoff()
            else:
                sys.exit("%s spectral model not recognized."%(type))
            
            emin=file["SOURCES"].data.field('EMIN')[src]
            emax=file["SOURCES"].data.field('EMAX')[src]
            alpha=file["SOURCES"].data.field('ALPHA')[src]
            err_alpha=file["SOURCES"].data.field('ERR_ALPHA')[src]
            energy=file["SOURCES"].data.field('MODEL_ENERGY')[src]
            model_exposure=file["SOURCES"].data.field('MODEL_EXPOSURE')[src]
            nphoton=file["SOURCES"].data.field('NPHOTON')[src]
            params=file["SOURCES"].data.field('PAR')[src]
            param_errors=file["SOURCES"].data.field('PAR_ERR')[src]

            # Spectral fitting energy range
            fit_range=file["SOURCES"].data.field('FIT_RANGE')[src]
            model.set_energy_range(R.Double(fit_range[0]),R.Double(fit_range[1]))
            lower_range=R.Double(fit_range[0])
            upper_range=R.Double(fit_range[1])

            # Band upper limits
            band_upper_limits=file["SOURCES"].data.field('BAND_UPPER_LIMIT')[src]
            energy_upper_limits=file["SOURCES"].data.field('ENERGY_UPPER_LIMIT')[src]
            exposure_upper_limits=file["SOURCES"].data.field('EXPOSURE_UPPER_LIMIT')[src]

            if model.get_spec_type()=="POWER_LAW":
                scale=R.Double(file["SOURCES"].data.field('DECORRELATION_ENERGY')[src])
                print "Pivot energy = %.1f MeV"%(scale)
                model.set_scale(scale)

    file.close()

    print params

    global func
    func=set_function(model,params)

    global butterfly
    butterfly=set_butterfly(model,params,param_errors)

    livebin=0
    for bin in range(0,int(len(emin)/2)):

        # Do not exceed spectral fitting range
        if livebin >= len(energy):
            break

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

        g1.SetPoint(bin,0,0)
        g1.SetPointError(bin,0,0,0,0)
        g2.SetPoint(bin,0,0)
        g2.SetPointError(bin,0,0,0,0)

        if flux==0:
            continue

        # Ensure energy bin was fitted
        livebin=livebin+2

        g1.SetPoint(bin,E_model,flux)
        g1.SetPointError(bin,E_model-E_min,E_max-E_model,err_flux,err_flux)

        model_flux=func.Eval(E_model)

        g2.SetPoint(bin,E_model,(flux-model_flux)/flux)
        g2.SetPointError(bin,E_model-E_min,E_max-E_model,err_flux/flux,err_flux/flux)

    for band in range(0,len(band_upper_limits)/2):
        front=2*band
        back=2*band+1
        
        E=energy_upper_limits[front]
        E_min=get_edge(emin,E,0)
        E_max=get_edge(emax,E,1)
        delta_E=E_max-E_min
        
        flux_front=pow(E,2)*band_upper_limits[front]/(exposure_upper_limits[front]*delta_E)
        flux_back=pow(E,2)*band_upper_limits[back]/(exposure_upper_limits[back]*delta_E)

        g3.SetPoint(band,E,flux_front+flux_back)
        g3.SetPointError(band,E-E_min,E_max-E,0,0)

    g1.Print()
    g2.Print()
    g3.Print()

    global canvas
    global pad1
    global pad2
    canvas,pad1,pad2=set_canvas(name)
    
    global box
    box=set_box(g1)

    global res_box
    res_box=set_res_box(g2)

    global legend
    legend=set_legend(name,model,g1,func,g3)

    canvas.Draw()
    box.Draw()
    butterfly.Draw("F")
    g1.Draw("P")
    g3.Draw("P")
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

    lower_bound=model.get_lower_bound()
    upper_bound=model.get_upper_bound()
    lower_bound=R.Double(10.)
    upper_bound=R.Double(5.e5)

    if model.get_spec_type()=="POWER_LAW":
        f1 = R.TF1("f1","pow(x,2)*[0]*pow(x/[2],-1*[1])",lower_bound,upper_bound)
        f1.SetParameters(params[0],params[1],model.get_scale())
    elif model.get_spec_type()=="BROKEN_POWER_LAW":
        f1 = R.TF1("f1","(x<[3])*pow(x,2)*[0]*pow(x/[3],-1*[1])+(x>[3])*pow(x,2)*[0]*pow(x/[3],-1*[2])",lower_bound,upper_bound);
        f1.SetParameters(params[0],params[1],params[2],params[3]);
    elif model.get_spec_type()=="EXP_CUTOFF":
        scale=100. # [ MeV ]
        f1 = R.TF1("f1","pow(x,2)*[0]*exp(-1.*x/[2])*pow(x/[3],-1*[1])",lower_bound,upper_bound);
        f1.SetParameters(params[0],params[1],params[2],scale)
    else:
        print "Invalid spectral model type"

    f1.SetLineWidth(1)
    f1.SetLineStyle(7)
    f1.SetLineColor(4)

    return f1

def set_butterfly(model,params,param_errors):
    print type(model)

    def CRD(obj):
        for i in range(0,len(obj)):
            obj[i]=R.Double(obj[i])
        return obj

    scale = R.Double(model.get_scale())
    lower_bound  = R.Double(max(model.get_lower_bound(),100))
    upper_bound  = R.Double(min(model.get_upper_bound(),3e5))

    b1=R.TPolyLine()
    
    if model.get_spec_type()=="POWER_LAW":
    
        model.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1]]))
        b1.SetPoint(0,lower_bound,R.Double(pow(lower_bound,2)*model.get_dNdE(lower_bound)))

        b1.SetPoint(1,scale,R.Double(pow(scale,2)*model.get_dNdE(scale)))

        model.set_params(CRD([params[0]+param_errors[0],params[1]-param_errors[1]]))
        b1.SetPoint(2,upper_bound,R.Double(pow(upper_bound,2)*model.get_dNdE(upper_bound)))

        model.set_params(CRD([params[0]-param_errors[0],params[1]+param_errors[1]]))
        b1.SetPoint(3,upper_bound,R.Double(pow(upper_bound,2)*model.get_dNdE(upper_bound)))

        b1.SetPoint(4,scale,R.Double(pow(scale,2)*model.get_dNdE(scale)))

        model.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1]]))
        b1.SetPoint(5,lower_bound,R.Double(pow(lower_bound,2)*model.get_dNdE(lower_bound)))

        model.set_params([float(params[0]),float(params[1])])

    if model.get_spec_type()=="BROKEN_POWER_LAW":

        ebreak=R.Double(params[3])

        model.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2],params[3]]))
        b1.SetPoint(0,lower_bound,R.Double(pow(lower_bound,2)*model.get_dNdE(lower_bound)))

        b1.SetPoint(1,ebreak,R.Double(pow(ebreak,2)*model.get_dNdE(ebreak)))

        model.set_params(CRD([params[0]+param_errors[0],params[1],params[2]-param_errors[2],params[3]]))
        b1.SetPoint(2,upper_bound,R.Double(pow(upper_bound,2)*model.get_dNdE(upper_bound)))

        model.set_params(CRD([params[0]-param_errors[0],params[1],params[2]+param_errors[2],params[3]]))
        b1.SetPoint(3,upper_bound,R.Double(pow(upper_bound,2)*model.get_dNdE(upper_bound)))

        b1.SetPoint(4,ebreak,R.Double(pow(ebreak,2)*model.get_dNdE(ebreak)))

        model.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2],params[3]]))
        b1.SetPoint(5,lower_bound,R.Double(pow(lower_bound,2)*model.get_dNdE(lower_bound)))

        model.set_params([float(params[0]),float(params[1]),float(params[2]),float(params[3])])

    if model.get_spec_type()=="EXP_CUTOFF":

        m_lo_lo=pl.ExpCutoff()
        m_lo_hi=pl.ExpCutoff()
        m_hi_lo=pl.ExpCutoff()
        m_hi_hi=pl.ExpCutoff()

        m_lo_lo.set_params(CRD([params[0]+param_errors[0],params[1]-param_errors[1],params[2]-param_errors[2]]))
        m_lo_hi.set_params(CRD([params[0]+param_errors[0],params[1]-param_errors[1],params[2]+param_errors[2]]))
        m_hi_lo.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2]-param_errors[2]]))
        m_hi_hi.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2]+param_errors[2]]))

        #m_lo.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2]-param_errors[2]]))
        #m_hi.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2]+param_errors[2]]))

        model.set_params(CRD([params[0]+param_errors[0],params[1]+param_errors[1],params[2]]))

        itr=100

        for i in range(0,itr+1):
            frac=float(i)/float(itr)
            E=exp((1-frac)*log(lower_bound)+frac*log(upper_bound))
            E2dNdE=R.Double(max(pow(E,2)*m_lo_lo.get_dNdE(E),
                                pow(E,2)*m_lo_hi.get_dNdE(E),
                                pow(E,2)*m_hi_lo.get_dNdE(E),
                                pow(E,2)*m_hi_hi.get_dNdE(E)))
            b1.SetPoint(i,E,E2dNdE)

        m_lo_lo.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2]-param_errors[2]]))
        m_lo_hi.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2]+param_errors[2]]))
        m_hi_lo.set_params(CRD([params[0]-param_errors[0],params[1]+param_errors[1],params[2]-param_errors[2]]))
        m_hi_hi.set_params(CRD([params[0]-param_errors[0],params[1]+param_errors[1],params[2]+param_errors[2]]))

        #m_lo.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2]-param_errors[2]]))
        #m_hi.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2]+param_errors[2]]))

        model.set_params(CRD([params[0]-param_errors[0],params[1]-param_errors[1],params[2]]))

        for i in range(itr+1,2*itr+2):
            frac=float(i-itr-1)/float(itr)
            E=exp(frac*log(lower_bound)+(1-frac)*log(upper_bound))
            E2dNdE=R.Double(min(pow(E,2)*m_lo_lo.get_dNdE(E),
                                pow(E,2)*m_lo_hi.get_dNdE(E),
                                pow(E,2)*m_hi_lo.get_dNdE(E),
                                pow(E,2)*m_hi_hi.get_dNdE(E)))
            b1.SetPoint(i,E,E2dNdE)

        model.set_params([float(params[0]),float(params[1]),float(params[2])])   

    b1.SetLineWidth(1)
    b1.SetLineStyle(1)#7
    b1.SetLineColor(2)#4
    b1.SetFillColor(R.kBlue-10)
    
    return b1

def set_legend(name,model,graph,func,upper_limit):
    legend = R.TLegend(0.7,0.8,0.85,0.9)
    legend.SetHeader(name)
    legend.AddEntry(g1,"Data","P")
    legend.AddEntry(func,model.get_spec_type().replace('_',' '),"L")
    legend.AddEntry(g3,"Upper Limit","L")
    legend.SetFillColor(0)
    legend.SetShadowColor(0)
    legend.SetLineColor(0)
    return legend

def get_edge(array,E,num):
    edge=-1
    for i in range(0,len(array)-1):
        if E>array[i] and E<array[i+1]:
                edge=array[i+num]
    return edge
