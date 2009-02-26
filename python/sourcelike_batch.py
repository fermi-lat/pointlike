import sys
import os
from pointlike import SkyDir
import pyfits

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-m", "--mode", dest="mode", default='catalog_fit',
                  help="Mode to run sourcelike_batch in. default:catalog_fit")
parser.add_option("-c", "--catalog", dest="catalog", default='',
                  help="Name of the source catalog")
parser.add_option("-o", "--outdir", dest="outdir", default='',
                  help="Name of output directory for jobs")
parser.add_option("-r", "--range", dest="range", default="0:9999999",
                  help="source range in format <first>:<last>")
parser.add_option("-s", "--submit", action="store_true", dest="submit",
                  help="Do it. Don't just print submit commands.")
parser.add_option("-l", "--run-local", action="store_true", dest="local",
                  help="Run commands locally.")

(options, args) = parser.parse_args()

if(options.catalog=='' or options.outdir=='' or not os.path.exists(options.catalog)):
   parser.print_help()
   sys.exit(1)

#if( os.path.exists(options.outdir)): raise RuntimeError,"Output directory already exists."

cfile=pyfits.open(options.catalog)
outdir=options.outdir

rtok=options.range.split(':')
if len(rtok)!=2:  raise RuntimeError,"Invalid range. Format is first:last"
imin,imax=int(rtok[0]),int(rtok[1])

if not os.path.exists(outdir): os.mkdir(outdir)
if not os.path.exists(outdir+'/log'): os.mkdir(outdir+'/log')
if not os.path.exists(outdir+'/images'): os.mkdir(outdir+'/images')

batchcmd= 'bsub -q xlong -R rhel40'

nametab=cfile["LAT_Point_Source_Catalog"].data.field("NickName")
ratab=cfile["LAT_Point_Source_Catalog"].data.field("RA")
dectab=cfile["LAT_Point_Source_Catalog"].data.field("DEC")
tstab=cfile["LAT_Point_Source_Catalog"].data.field("Test_Statistic")
cfile.close()

source=[]
for name,ra,dec,ts in zip(nametab,ratab,dectab,tstab):
   outfile=outdir+'/'+name+'.fits'
   imagedir=outdir+'/images'
   logfile=outdir+'/log/'+name+'.log'
   elogfile=outdir+'/log/'+name+'.errlog'   
   source.append({'id':name,'ra':float(ra),'dec':float(dec),'ts':float(ts),
                  'out':outfile,'img':imagedir,'log':logfile,'elog':elogfile})


for i in range(len(source)):
   if (i<imin or i>imax): continue
   print source[i]['id'],source[i]['ra'],source[i]['dec'],"TS=",source[i]['ts']
   batchcmd=batchcmd+' -o %s -e %s'%(source[i]['log'],source[i]['elog'])
   
   if(options.mode=="catalog_fit"):
       imgdir=source[i]['img']+'/'+source[i]['id']
       if not os.path.exists(imgdir): os.mkdir(imgdir)
       cmd='python sourcelike.py -s %s -o %s -i %s -r %.3f -d %.3f -e 500. -c -T -S'\
	  %(source[i]['id'],source[i]['out'],imgdir,source[i]['ra'],source[i]['dec'])
       if (source[i]['ts']>200): cmd+=' -f plaw,bplaw,expcut'    
       if (source[i]['ts']>1000): cmd+=' -m'    
       sdir=SkyDir(source[i]['ra'],source[i]['dec'],SkyDir.EQUATORIAL)
       for j in range(len(source)):
	 sdir2=SkyDir(source[j]['ra'],source[j]['dec'],SkyDir.EQUATORIAL)
	 dpsi=sdir.difference(sdir2)
	 if(dpsi*57.296<10 and source[j]['id']!=source[i]['id']): cmd+=' -B %s:%.4f,%.4f'%(source[j]['id'],source[j]['ra'],source[j]['dec'])
   
   elif(options.mode=="tsmap_fit"):
       imgdir=source[i]['img']+'/'+source[i]['id']
       tsmap=imgdir+'/tsmap_'+source[i]['id']+'.fits'
       if not os.path.exists(imgdir): os.mkdir(imgdir)
       cmd='python sourcelike.py -s %s -F %s -o %s -i %s -e 500. -c -T'\
	  %(source[i]['id'],tsmap,source[i]['out'],imgdir)
   
   elif(options.mode=="tsmap_gen"):
       imgdir=source[i]['img']+'/'+source[i]['id']
       if not os.path.exists(imgdir): os.mkdir(imgdir)
       cmd='python tsmap.py -t %s -r %.3f -d %.3f -e 500. -f 10. -s 0.125'\
	  %(imgdir+'/tsmap_'+source[i]['id']+'.fits',source[i]['ra'],source[i]['dec'])
   
   else: raise RuntimeError,"Mode %s unknown."%(options.mode)
   
   if(options.submit): os.system("%s %s"%(batchcmd,cmd))
   elif(options.local): os.system("%s"%(cmd))
   else: os.system("echo %s %s"%(batchcmd,cmd)) 
   
