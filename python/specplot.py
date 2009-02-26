import numpy as N
import ROOT as R
import pointlike as pl 
import pyfits
import sys
from math import floor,exp,log,sqrt,pow,fabs

def get_edge(array,x,which):
    edge=-1
    for i in range(0,len(array)-1):
        if x>array[i] and x<array[i+1]:
                edge=array[i+which]
    return edge

def unique(array):
    unique=1
    for i in range(0,len(array)):
        for j in range(0,i):
            if array[i]==array[j]:
                unique=0
    return unique

class SpectralPlotter:

    def __init__(self,filename,source=""):
        self.filename=filename
        self.source=source
        self.func_form=""
        self.mwfile=""
        self.sed_lower_range=100.
        self.sed_upper_range=3.e5
        self.read_fits()
        self.set_graphs()
        
    def set_mwfile(self,mwfile,mwscale=1.):
        self.mwfile=mwfile
        MWData = pl.MWData(mwfile,mwscale)
        self.mw_E           = MWData.get_E()
        self.mw_dNdE        = MWData.get_dNdE()
        self.mw_dNdE_err_lo = MWData.get_dNdE_err_lo()
        self.mw_dNdE_err_hi = MWData.get_dNdE_err_hi()

    def add_func(self,func_form):
        self.func_form=func_form
        self.comp_func = R.TF1("comp_func",func_form,1,1e9)
        self.comp_func.SetLineStyle(1)
        self.comp_func.SetLineColor(1)
        self.comp_func.SetLineWidth(1)
            
    def SED(self):
        self.set_func()
        self.set_lat_points()
        self.set_mw_points()

        try:
            self.UL()
        except:
            "No upper limit information available - displaying SED only"

        self.set_canvas()
        self.set_sed_box()
        self.set_res_box()

        self.sed_canvas.Draw()
        self.sed_box.Draw()

        self.set_butterfly()
    
        self.butterfly.Draw("F")
        self.gLAT.Draw("P")
        self.func.Draw("SAME")
        if self.func_form != "":
            self.comp_func.Draw("SAME")

        self.set_selected_UL()
        self.gULselected.Draw("P")

        for i in range(0,len(self.selected_arrows)):
            self.selected_arrows[i].Draw()

        if self.gMW.GetN()>0:
            self.gMW.Draw("P")

#        Overlay all upper limits
#        self.gUL.Draw("P")
#        for i in range(0,len(self.arrows)):
#            self.arrows[i].Draw()

        self.set_legend()
        self.legend.Draw()

        self.set_stats()
        self.stats.Draw("SAME")

        self.res_pad.cd()
        self.res_box.Draw()
        self.gLATres.Draw("P")
        if self.gMWres.GetN()>0:
            self.gMWres.Draw("P")
        self.ref = R.TF1("ref","0",1,1e9)
        self.ref.SetLineWidth(1)
        self.ref.SetLineStyle(7)
        self.ref.SetLineColor(4)
        self.ref.Draw("SAME")
    
        self.sed_canvas.Update()
        
    def read_fits(self):
        file=pyfits.open(self.filename)

        self.name=""
        if self.source=="" : self.source=file["SOURCES"].data.field('NAME')[0]
    
        for src in range(0,len(file["SOURCES"].data.field('NAME'))):
            self.name=file["SOURCES"].data.field('NAME')[src]
            if self.name==self.source:
                print "Found "+self.name 

                self.emin=file["SOURCES"].data.field('EMIN')[src]
                self.emax=file["SOURCES"].data.field('EMAX')[src]
                self.nphoton=file["SOURCES"].data.field('NPHOTON')[src]
                self.alpha=file["SOURCES"].data.field('ALPHA')[src]
                self.err_alpha=file["SOURCES"].data.field('ERR_ALPHA')[src]

                try:
                    type=str(file["SOURCES"].data.field('MODEL')[src])
                    print "spectral model = %s"%(type)
                    if type=="POWER_LAW":
                        self.model=pl.PowerLaw()
                    elif type=="BROKEN_POWER_LAW":
                        self.model=pl.BrokenPowerLaw()
                    elif type=="EXP_CUTOFF":
                        self.model=pl.ExpCutoff()
                    else:
                        sys.exit("%s spectral model not recognized."%(type))

                    if self.model.get_spec_type()=="POWER_LAW":
                        self.scale=R.Double(file["SOURCES"].data.field('DECORRELATION_ENERGY')[src])
                        print "Pivot energy = %.1f MeV"%(self.scale)
                        self.model.set_scale(self.scale)

                    self.energy=file["SOURCES"].data.field('MODEL_ENERGY')[src]
                    self.model_exposure=file["SOURCES"].data.field('MODEL_EXPOSURE')[src]
                    self.params=file["SOURCES"].data.field('PAR')[src]
                    self.param_errors=file["SOURCES"].data.field('PAR_ERR')[src]

                    self.flux=file["SOURCES"].data.field('INTEGRAL_FLUX')[src]
                    self.flux_err_hi=file["SOURCES"].data.field('INTEGRAL_FLUX_ERR')[src][0]
                    self.flux_err_lo=file["SOURCES"].data.field('INTEGRAL_FLUX_ERR')[src][1]
                except:
                    print "Did not find spectral fit information in results file"

                # Spectral fitting energy range
                self.fit_range=file["SOURCES"].data.field('FIT_RANGE')[src]
                self.model.set_energy_range(R.Double(self.fit_range[0]),R.Double(self.fit_range[1]))
                self.lower_range=R.Double(self.fit_range[0])
                self.upper_range=R.Double(self.fit_range[1])
                
                # Band upper limits
                try:
                    self.band_upper_limits=file["SOURCES"].data.field('BAND_UPPER_LIMIT')[src]
                    self.energy_upper_limits=file["SOURCES"].data.field('ENERGY_UPPER_LIMIT')[src]
                    self.exposure_upper_limits=file["SOURCES"].data.field('EXPOSURE_UPPER_LIMIT')[src]
                except:
                    print "Did not find any upper limit data in results file"

                # Multiwavelength data
                try:
                    self.mwfile=file["SOURCES"].data.field('MWFILE')[src]
                    self.mw_E=file["SOURCES"].data.field('MW_E')[src]
                    self.mw_dNdE=file["SOURCES"].data.field('MW_DNDE')[src]
                    self.mw_dNdE_err_lo=file["SOURCES"].data.field('MW_DNDE_ERR_LO')[src]
                    self.mw_dNdE_err_hi=file["SOURCES"].data.field('MW_DNDE_ERR_HI')[src]
                except:
                    print "Did not find any multiwavelength data in results file"
                
        file.close()
        self.combined=unique(self.emin)

    def set_graphs(self):
        # LAT data
        self.gLAT=R.TGraphAsymmErrors()
        self.gLAT.SetMarkerStyle(7)
        self.gLAT.SetMarkerColor(1)

        # LAT residuals
        self.gLATres=R.TGraphAsymmErrors()
        self.gLATres.SetMarkerStyle(7)
        self.gLATres.SetMarkerColor(1)

        # Multiwavelength data
        self.gMW=R.TGraphAsymmErrors()
        self.gMW.SetMarkerStyle(25)
        self.gMW.SetMarkerColor(1)
        self.gMW.SetMarkerSize(0.5)

        # Multiwavelength residuals
        self.gMWres=R.TGraphAsymmErrors()
        self.gMWres.SetMarkerStyle(25)
        self.gMWres.SetMarkerColor(1)
        self.gMWres.SetMarkerSize(0.5)

        # Upper limits
        self.gUL=R.TGraphAsymmErrors()
        self.gUL.SetMarkerStyle(1)
        self.gUL.SetMarkerColor(2)
        self.gUL.SetLineColor(2)
        self.arrows=[]

        # Selected upper limits
        self.gULselected=R.TGraphAsymmErrors()
        self.gULselected.SetMarkerStyle(1)
        self.gULselected.SetMarkerColor(2)
        self.gULselected.SetLineColor(2)
        self.selected_arrows=[]

    def set_func(self):
        lower_bound=R.Double(1.)
        upper_bound=R.Double(1.e9)

        if self.model.get_spec_type()=="POWER_LAW":
            self.func = R.TF1("func","pow(x,2)*[0]*pow(x/[2],-1*[1])",lower_bound,upper_bound)
            self.func.SetParameters(self.params[0],self.params[1],self.model.get_scale())
        elif self.model.get_spec_type()=="BROKEN_POWER_LAW":
            self.func = R.TF1("func","(x<[3])*pow(x,2)*[0]*pow(x/[3],-1*[1])+(x>[3])*pow(x,2)*[0]*pow(x/[3],-1*[2])",lower_bound,upper_bound);
            self.func.SetParameters(self.params[0],self.params[1],self.params[2],self.params[3]);
        elif self.model.get_spec_type()=="EXP_CUTOFF":
            scale=100. # [ MeV ]
            self.func = R.TF1("func","pow(x,2)*[0]*exp(-1.*x/[2])*pow(x/[3],-1*[1])",lower_bound,upper_bound);
            self.func.SetParameters(self.params[0],self.params[1],self.params[2],scale)
        else:
            print "Invalid spectral model type"

        self.func.SetLineWidth(1)
        self.func.SetLineStyle(7)
        self.func.SetLineColor(4)

    def set_lat_points(self):
        print "Setting LAT points..."
        npoints=0        
        if self.combined==0:
            print "Separate front and back energy bins"
            livebin=0
            for bin in range(0,len(self.emin)-1):

                if(self.emin[bin+1]!=self.emin[bin] or self.energy[-1]<self.emin[bin] or self.energy[0]>self.emax[bin]):continue
            
                front,back=bin,bin+1

                while (livebin < len(self.energy)-2 and self.energy[livebin]<self.emin[bin]): livebin+=1
                
                if(abs(self.energy[livebin+1]-self.energy[livebin])>1e-5): raise RuntimeError,"Corrupt front/back structure in file." 

                E_min=self.emin[front]
                E_max=self.emax[front]
                E_model=self.energy[livebin]

                if E_model<E_min or E_model>E_max: raise RuntimeError,"Corrupt binning in file." 
        
                exposure_error_fraction=self.model.get_exposure_uncertainty(R.Double(E_min),R.Double(E_max))
                delta_E=E_max-E_min

                if self.nphoton[front]>0:
                    flux_front=pow(E_model,2)*self.nphoton[front]*self.alpha[front]/(self.model_exposure[livebin]*delta_E)
                    err_flux_front=flux_front*sqrt(pow(self.err_alpha[front]/self.alpha[front],2)+pow(exposure_error_fraction,2))
                else:
                    flux_front=0
                    err_flux_front=0

                if self.nphoton[back]>0:
                    flux_back=pow(E_model,2)*self.nphoton[back]*self.alpha[back]/(self.model_exposure[livebin+1]*delta_E)
                    err_flux_back=flux_back*sqrt(pow(self.err_alpha[back]/self.alpha[back],2)+pow(exposure_error_fraction,2))
                else:
                    flux_back=0
                    err_flux_back=0
                
                flux=(flux_front+flux_back)/2
                err_flux=sqrt(pow(err_flux_front,2)+pow(err_flux_back,2))

                self.gLAT.SetPoint(npoints,E_model,flux)
                self.gLAT.SetPointError(npoints,E_model-E_min,E_max-E_model,err_flux,err_flux)

                model_flux=self.func.Eval(E_model)
            
                if flux!=0:
                    self.gLATres.SetPoint(npoints,E_model,(flux-model_flux)/flux)
                    self.gLATres.SetPointError(npoints,E_model-E_min,E_max-E_model,err_flux/flux,err_flux/flux)
                else:   
                    self.gLATres.SetPoint(npoints,E_model,0)
                    self.gLATres.SetPointError(npoints,0,0,0,0)

                npoints+=1

        if self.combined==1:
            print "Combined front and back energy bins"
            livebin=0
            for bin in range(0,len(self.emin)):

                # Do not exceed spectral fitting range
                if livebin >= len(self.energy):
                    break

                E_min=self.emin[bin]
                E_max=self.emax[bin]
                E_model=self.energy[livebin]

                if E_model<E_min or E_model>E_max:
                    continue
        
                exposure_error_fraction=self.model.get_exposure_uncertainty(R.Double(E_min),R.Double(E_max))
                delta_E=E_max-E_min

                if self.nphoton[bin]>0:
                    flux=pow(E_model,2)*self.nphoton[bin]*self.alpha[bin]/(self.model_exposure[livebin]*delta_E)
                    err_flux=flux*sqrt(pow(self.err_alpha[bin]/self.alpha[bin],2)+pow(exposure_error_fraction,2))
                else:
                    flux=0
                    err_flux=0
        
                E_model=self.energy[livebin]

                self.gLAT.SetPoint(npoints,E_model,flux)
                self.gLAT.SetPointError(npoints,E_model-E_min,E_max-E_model,err_flux,err_flux)

                model_flux=self.func.Eval(E_model)
            
                if flux!=0: 
                    self.gLATres.SetPoint(npoints,E_model,(flux-model_flux)/flux)
                    self.gLATres.SetPointError(npoints,E_model-E_min,E_max-E_model,err_flux/flux,err_flux/flux)
                else:
                    self.gLATres.SetPoint(npoints,E_model,0)
                    self.gLATres.SetPointError(npoints,0,0,0,0)

                # Ensure energy bin was fitted
                npoints+=1
                livebin=livebin+1

    def set_mw_points(self):
        if self.mwfile!="":
            # Adjust energy range for plotting
            if(0.8*min(self.mw_E)<100): self.sed_lower_range=0.8*min(self.mw_E)
            if(1.2*max(self.mw_E)>3e5): self.sed_upper_range=1.2*max(self.mw_E)
            
            for i in range(0,len(self.mw_E)):
                E=self.mw_E[i]
                dNdE=self.mw_dNdE[i]
                dNdE_lo=self.mw_dNdE[i]-self.mw_dNdE_err_lo[i]
                dNdE_hi=self.mw_dNdE[i]+self.mw_dNdE_err_hi[i]

                flux=pow(E,2)*dNdE
                flux_hi=pow(E,2)*dNdE_hi
                flux_lo=pow(E,2)*dNdE_lo
            
                self.gMW.SetPoint(i,E,flux)
                self.gMW.SetPointError(i,0,0,flux-flux_lo,flux_hi-flux)

                model_flux=self.func.Eval(E)

                self.gMWres.SetPoint(i,E,(flux-model_flux)/flux)
                self.gMWres.SetPointError(i,0,0,(flux-flux_lo)/flux,(flux_hi-flux)/flux)

    def set_canvas(self):

        R.gStyle.SetOptStat(0)
        R.gStyle.SetOptTitle(0)
    
        self.sed_canvas = R.TCanvas("SED","SED: "+self.name,200,10,500,700)
        self.sed_canvas.Draw()

        self.sed_pad = R.TPad("SpecPad","SpecPad",0.0,0.3,1.0,1.0);

        self.sed_pad.SetLogx()
        self.sed_pad.SetLogy()

        self.sed_pad.SetLeftMargin(0.14)
        self.sed_pad.SetBottomMargin(0.1)
        self.sed_pad.SetTopMargin(0.04)
        self.sed_pad.SetRightMargin(0.04)
        self.sed_pad.SetFillColor(10)
        self.sed_pad.Draw();

        self.sed_canvas.cd()
        self.res_pad = R.TPad("ResPad2","ResPad",0.0,0.0,1.0,0.3)
        self.res_pad.Draw()
        self.res_pad.cd()

        self.res_pad.SetLogx()
        
        self.res_pad.SetLeftMargin(0.14)
        self.res_pad.SetBottomMargin(0.14)
        self.res_pad.SetTopMargin(0.04)
        self.res_pad.SetRightMargin(0.04)
        self.res_pad.SetFillColor(10)

        self.sed_pad.cd();
        
    def set_sed_box(self):
        xmin=R.Double()
        xmax=R.Double()
        ymin=R.Double()
        ymax=R.Double()
        self.gLAT.ComputeRange(xmin,ymin,xmax,ymax)
        
        flux_min=max(1e-8,ymin*R.Double(0.1))
        flux_max=ymax*R.Double(10.)
        
        self.sed_box = R.TH2F("Spectral","Spectrum",10,self.sed_lower_range,self.sed_upper_range,10,flux_min,flux_max)
        self.sed_box.SetTitle("Spectral Energy Distribution")
        self.sed_box.GetXaxis().SetTitle("Energy [MeV]")
        self.sed_box.GetXaxis().SetTitleFont(42)
        self.sed_box.GetXaxis().CenterTitle()
        self.sed_box.GetXaxis().SetTitleOffset(1.1)
        self.sed_box.GetYaxis().SetTitle("E^{2} dN/dE [MeV cm^{-2} s^{-1}]")
        self.sed_box.GetYaxis().SetTitleFont(42)
        self.sed_box.GetYaxis().CenterTitle()
        self.sed_box.GetYaxis().SetTitleOffset(1.5)
        
    def set_res_box(self):
        xmin=R.Double()
        xmax=R.Double()
        ymin=R.Double()
        ymax=R.Double()
        self.gLATres.ComputeRange(xmin,ymin,xmax,ymax)
        
        res_min=max(min(0,ymin*R.Double(1.2)),-6)
        res_max=min(max(0,ymax*R.Double(1.2)),6)
        
        self.res_box = R.TH2F("RESIDUAL","RESIDUAL",10,self.sed_lower_range,self.sed_upper_range,10,res_min,res_max)
        self.res_box.SetTitle("Spectral Energy Distribution Residual")
        self.res_box.GetXaxis().SetTitle("Energy [MeV]")
        self.res_box.GetXaxis().SetTitleFont(42)
        self.res_box.GetXaxis().CenterTitle()
        self.res_box.GetXaxis().SetTitleOffset(1.1)
        self.res_box.GetYaxis().SetTitle("Fit Residuals - \DeltaF/F")
        self.res_box.GetYaxis().SetTitleFont(42)
        self.res_box.GetYaxis().CenterTitle()
        self.res_box.GetYaxis().SetTitleOffset(1.5)
        
    def set_legend(self):
        self.legend = R.TLegend(0.7,0.75,0.9,0.9)
        self.legend.SetHeader(self.name)
        self.legend.AddEntry(self.gLAT,"LAT Data","P")
        self.legend.AddEntry(self.func,self.model.get_spec_type().replace('_',' '),"L")
        if self.gULselected.GetN()>0:
            self.legend.AddEntry(self.gULselected,"1\sigma Upper Limit","L")
        if self.func_form!="":
            self.legend.AddEntry(self.comp_func,"Comparison","L")
        if self.gMW.GetN()>0:
            mwfilename="MW Data"
            if self.mwfile!="": mwfilename=self.mwfile
            self.legend.AddEntry(self.gMW,mwfilename,"P")
        self.legend.SetFillColor(0)
        self.legend.SetShadowColor(0)
        self.legend.SetLineColor(0)

    def UL(self):
        print "Setting upper limits..."
        npoints=0
        if self.combined==0:
            print "Separate front and back energy bins"	
            for band in range(0,len(self.band_upper_limits)/2):
                front=2*band
                back=2*band+1
            
                E=self.energy_upper_limits[front]
                E_min=get_edge(self.emin,E,0)
                E_max=get_edge(self.emax,E,1)
                delta_E=E_max-E_min
            
                flux_front=pow(E,2)*self.band_upper_limits[front]/(self.exposure_upper_limits[front]*delta_E)
                flux_back=pow(E,2)*self.band_upper_limits[back]/(self.exposure_upper_limits[back]*delta_E)

                if(E<=0): continue  
	      
                self.gUL.SetPoint(npoints,E,flux_front+flux_back)
                self.gUL.SetPointError(npoints,E-E_min,E_max-E,0,0)
                npoints+=1

        if self.combined==1:
            print "Combined front and back energy bins"
            for band in range(0,len(self.band_upper_limits)):
                E=self.energy_upper_limits[band]
                E_min=get_edge(self.emin,E,0)
                E_max=get_edge(self.emax,E,1)
                delta_E=E_max-E_min
        
                flux=pow(E,2)*self.band_upper_limits[band]/(self.exposure_upper_limits[band]*delta_E)
                
                if(E<=0): continue  

                self.gUL.SetPoint(npoints,E,flux)
                self.gUL.SetPointError(npoints,E-E_min,E_max-E,0,0)
                npoints+=1

        R.gStyle.SetOptStat(0)
        R.gStyle.SetOptTitle(0)

        self.ul_canvas = R.TCanvas("UPPER LIMIT SED","UPPER LIMIT SED: "+self.name,200,10,500,500)
        self.ul_canvas.Draw()

        self.ul_pad = R.TPad("UL Pad","UL Pad",0.0,0.0,1.0,1.0);
        self.ul_pad.SetLogx()
        self.ul_pad.SetLogy()
        self.ul_pad.SetLeftMargin(0.14)
        self.ul_pad.SetBottomMargin(0.1)
        self.ul_pad.SetTopMargin(0.04)
        self.ul_pad.SetRightMargin(0.04)
        self.ul_pad.SetFillColor(10)
        self.ul_pad.Draw();
        self.ul_pad.cd();
    
        xmin=R.Double()
        xmax=R.Double()
        ymin=R.Double()
        ymax=R.Double()
        self.gUL.ComputeRange(xmin,ymin,xmax,ymax)
   
        flux_min=max(1e-8,ymin*R.Double(0.1))
        flux_max=ymax*R.Double(10.)

        self.ul_box = R.TH2F("UL","UL",10,self.sed_lower_range,self.sed_upper_range,10,flux_min,flux_max)
        self.ul_box.SetTitle("Upper Limit Spectral Energy Distribution")
        self.ul_box.GetXaxis().SetTitle("Energy [MeV]")
        self.ul_box.GetXaxis().SetTitleFont(42)
        self.ul_box.GetXaxis().CenterTitle()
        self.ul_box.GetXaxis().SetTitleOffset(1.1)
        self.ul_box.GetYaxis().SetTitle("E^{2} dN/dE [MeV cm^{-2} s^{-1}]")
        self.ul_box.GetYaxis().SetTitleFont(42)
        self.ul_box.GetYaxis().CenterTitle()
        self.ul_box.GetYaxis().SetTitleOffset(1.5)

        self.ul_legend = R.TLegend(0.7,0.75,0.9,0.9)
        self.ul_legend.SetHeader(self.name)
        self.ul_legend.AddEntry(self.gUL,"1\sigma Upper Limit","L")
        self.ul_legend.SetFillColor(0)
        self.ul_legend.SetShadowColor(0)
        self.ul_legend.SetLineColor(0)

        self.ul_canvas.Draw()
        self.ul_box.Draw()
        self.gUL.Draw("P")

        for i in range(0,self.gUL.GetN()):
            x=R.Double(0)
            y=R.Double(0)
            self.gUL.GetPoint(i,x,y)
            #print '%.3e   %.3e'%(x,y)
            E_min=get_edge(self.emin,x,0)
            E_max=get_edge(self.emax,x,1)
            x=R.Double(exp((log(E_min)+log(E_max))/2.))
            #print '%5i  %5i   %.3e'%(E_min,E_max,y)
            x1=x
            x2=x
            y1=y
            y2=R.Double(0.5*y1)
            self.arrows.append(R.TArrow(x1,y1,x2,y2,0.02,">"))
            self.arrows[i].SetAngle(60)
            self.arrows[i].SetOption(">")
            self.arrows[i].SetLineColor(2)

        for i in range(0,len(self.arrows)):
            self.arrows[i].Draw()
        self.ul_legend.Draw()
        self.ul_canvas.Update()

    def set_selected_UL(self):
        points_kept=0

        x=R.Double(0)
        y=R.Double(0)

        SED_E=R.Double(0)
        SED_flux=R.Double(0)

        # seprate front/back case
        if self.combined==0:
            for i in range(0,self.gUL.GetN()):
                keep_point=1
                self.gUL.GetPoint(i,x,y)

                for j in range(0,len(self.emin)-1):
                    if self.nphoton[j]==0 or self.emin[j]!=self.emin[j+1]: continue
                    E_min=self.emin[j]
                    E_max=self.emax[j]
                    if x<E_min or x>E_max: continue

                    for k in range(0,self.gLAT.GetN()):
                        self.gLAT.GetPoint(k,SED_E,SED_flux)
                        if SED_E<E_min or SED_E>E_max: continue
                    
                        if SED_flux>1.e-9:
                            keep_point=0

                if keep_point==1:
                    self.gULselected.SetPoint(points_kept,x,y)
                    self.gULselected.SetPointError(points_kept,self.gUL.GetErrorXlow(i),self.gUL.GetErrorXhigh(i),0,0)
                    points_kept+=keep_point

        #combined front/back binning case
        else:
            for i in range(0,self.gUL.GetN()):
                keep_point=1
                self.gUL.GetPoint(i,x,y)

                for j in range(0,len(self.emin)):
                    if self.nphoton[j]==0: continue
                    E_min=self.emin[j]
                    E_max=self.emax[j]
                    if x<E_min or x>E_max: continue

                    for k in range(0,self.gLAT.GetN()):
                        self.gLAT.GetPoint(k,SED_E,SED_flux)
                        if SED_E<E_min or SED_E>E_max: continue

                        if SED_flux>2.e-8:
                            keep_point=0

                if keep_point==1:
                    self.gULselected.SetPoint(points_kept,x,y)
                    self.gULselected.SetPointError(points_kept,self.gUL.GetErrorXlow(i),self.gUL.GetErrorXhigh(i),0,0)
                    points_kept+=keep_point

        for i in range(0,self.gULselected.GetN()):
            x=R.Double(0)
            y=R.Double(0)
            self.gULselected.GetPoint(i,x,y)
            E_min=get_edge(self.emin,x,0)
            E_max=get_edge(self.emax,x,1)
            x=R.Double(exp((log(E_min)+log(E_max))/2.))
            x1=x
            x2=x
            y1=y
            y2=R.Double(0.5*y1)
            self.selected_arrows.append(R.TArrow(x1,y1,x2,y2,0.02,">"))
            self.selected_arrows[i].SetAngle(60)
            self.selected_arrows[i].SetOption(">")
            self.selected_arrows[i].SetLineColor(2)

    def set_butterfly(self):
        def CRD(obj):
            for i in range(0,len(obj)):
                obj[i]=R.Double(obj[i])
            return obj

        scale = R.Double(self.model.get_scale())
        lower_bound  = R.Double(max(self.model.get_lower_bound(),self.sed_lower_range))
        upper_bound  = R.Double(min(self.model.get_upper_bound(),self.sed_upper_range))

        self.butterfly=R.TPolyLine()
    
        if self.model.get_spec_type()=="POWER_LAW":
    
            self.model.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]+self.param_errors[1]]))
            self.butterfly.SetPoint(0,lower_bound,R.Double(pow(lower_bound,2)*self.model.get_dNdE(lower_bound)))

            self.butterfly.SetPoint(1,scale,R.Double(pow(scale,2)*self.model.get_dNdE(scale)))

            self.model.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]-self.param_errors[1]]))
            self.butterfly.SetPoint(2,upper_bound,R.Double(pow(upper_bound,2)*self.model.get_dNdE(upper_bound)))

            self.model.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]+self.param_errors[1]]))
            self.butterfly.SetPoint(3,upper_bound,R.Double(pow(upper_bound,2)*self.model.get_dNdE(upper_bound)))

            self.butterfly.SetPoint(4,scale,R.Double(pow(scale,2)*self.model.get_dNdE(scale)))

            self.model.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]-self.param_errors[1]]))
            self.butterfly.SetPoint(5,lower_bound,R.Double(pow(lower_bound,2)*self.model.get_dNdE(lower_bound)))

            self.model.set_params([float(self.params[0]),float(self.params[1])])

        if self.model.get_spec_type()=="BROKEN_POWER_LAW":

            ebreak=R.Double(self.params[3])

            self.model.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2],self.params[3]]))
            self.butterfly.SetPoint(0,lower_bound,R.Double(pow(lower_bound,2)*self.model.get_dNdE(lower_bound)))

            self.butterfly.SetPoint(1,ebreak,R.Double(pow(ebreak,2)*self.model.get_dNdE(ebreak)))

            self.model.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1],self.params[2]-self.param_errors[2],self.params[3]]))
            self.butterfly.SetPoint(2,upper_bound,R.Double(pow(upper_bound,2)*self.model.get_dNdE(upper_bound)))

            self.model.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1],self.params[2]+self.param_errors[2],self.params[3]]))
            self.butterfly.SetPoint(3,upper_bound,R.Double(pow(upper_bound,2)*self.model.get_dNdE(upper_bound)))

            self.butterfly.SetPoint(4,ebreak,R.Double(pow(ebreak,2)*self.model.get_dNdE(ebreak)))

            self.model.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2],self.params[3]]))
            self.butterfly.SetPoint(5,lower_bound,R.Double(pow(lower_bound,2)*self.model.get_dNdE(lower_bound)))

            self.model.set_params([float(self.params[0]),float(self.params[1]),float(self.params[2]),float(self.params[3])])

        if self.model.get_spec_type()=="EXP_CUTOFF":

            m_lo_lo=pl.ExpCutoff()
            m_lo_hi=pl.ExpCutoff()
            m_hi_lo=pl.ExpCutoff()
            m_hi_hi=pl.ExpCutoff()

            m_lo_lo.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2]-self.param_errors[2]]))
            m_lo_hi.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2]+self.param_errors[2]]))
            m_hi_lo.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2]-self.param_errors[2]]))
            m_hi_hi.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2]+self.param_errors[2]]))
            
            self.model.set_params(CRD([self.params[0]+self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2]]))

            itr=100

            for i in range(0,itr+1):
                frac=float(i)/float(itr)
                E=exp((1-frac)*log(lower_bound)+frac*log(upper_bound))
                E2dNdE=R.Double(max(pow(E,2)*m_lo_lo.get_dNdE(E),
                                    pow(E,2)*m_lo_hi.get_dNdE(E),
                                    pow(E,2)*m_hi_lo.get_dNdE(E),
                                    pow(E,2)*m_hi_hi.get_dNdE(E)))
                self.butterfly.SetPoint(i,E,E2dNdE)

            m_lo_lo.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2]-self.param_errors[2]]))
            m_lo_hi.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2]+self.param_errors[2]]))
            m_hi_lo.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2]-self.param_errors[2]]))
            m_hi_hi.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]+self.param_errors[1],self.params[2]+self.param_errors[2]]))

            self.model.set_params(CRD([self.params[0]-self.param_errors[0],self.params[1]-self.param_errors[1],self.params[2]]))

            for i in range(itr+1,2*itr+2):
                frac=float(i-itr-1)/float(itr)
                E=exp(frac*log(lower_bound)+(1-frac)*log(upper_bound))
                E2dNdE=R.Double(min(pow(E,2)*m_lo_lo.get_dNdE(E),
                                    pow(E,2)*m_lo_hi.get_dNdE(E),
                                    pow(E,2)*m_hi_lo.get_dNdE(E),
                                    pow(E,2)*m_hi_hi.get_dNdE(E)))
                self.butterfly.SetPoint(i,E,E2dNdE)

            self.model.set_params([float(self.params[0]),float(self.params[1]),float(self.params[2])])   

        self.butterfly.SetLineWidth(1)
        self.butterfly.SetLineStyle(1)
        self.butterfly.SetLineColor(2)
        self.butterfly.SetFillColor(R.kBlue-10)

    def set_stats(self):
        self.stats=R.TPaveText(2.e2,1.e-7,1.e4,1.e-6)
        self.stats.SetFillColor(0)
        self.stats.SetShadowColor(0)
        self.stats.SetLineColor(0)

        par_dict={"gamma":"#Gamma",
                  "gamma_1":"#Gamma_{1}",
                  "gamma_2":"#Gamma_{2}",
                  "E_b":"E_{b}",
                  "E_c":"E_{c}"}
        unit_dict={"gamma":"",
                  "gamma_1":"",
                  "gamma_2":"",
                  "E_b":"MeV",
                  "E_c":"MeV"}

        for i in range(1,len(self.params)):
            par_name=self.model.get_parName()[i]
            par_str="%s = %.2e +/- %.2e %s"%(par_dict[par_name],self.params[i],self.param_errors[i],unit_dict[par_name])
            self.stats.AddText(par_str)

        if self.model.get_spec_type()=="POWER_LAW":
            self.stats.AddText("Pivot Energy = %.2e MeV"%(self.scale))
        
        self.stats.AddText("F_{100} = %.2e +/- %.2e/%.2e ph cm^{-2} s^{-1}"%(self.flux,self.flux_err_hi,self.flux_err_lo))

if __name__ == "__main__":
    
    import sys
    
    if(len(sys.argv)==2):
        sp=SpectralPlotter(sys.argv[1])
    elif(len(sys.argv)==3):
        sp=SpectralPlotter(sys.argv[1],sys.argv[2])
    elif(len(sys.argv)==4):
        sp=SpectralPlotter(sys.argv[1],sys.argv[2])
        sp.add_func(argv[3])
    else:
       print "Usage: %s result-file <sourcename>"%(sys.argv[0])
       sys.exit(1)

    sp.SED()
    
    raw_input()

  

