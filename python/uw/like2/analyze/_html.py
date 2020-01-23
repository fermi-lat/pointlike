"""
Manage the Web page generation
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/_html.py,v 1.20 2015/08/16 01:11:36 burnett Exp $
"""
import os, glob, re, yaml
import pandas as pd
import numpy as np
#from uw.like2.analyze import app # for the menu

# for tool tips
tooltips_css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #ffffff;
  background-color: #000000;
}
td
{
  background-color: #cccccc;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
"""

style="""
<style type="text/css">
body, th, td {	font-family:verdana,arial,sans-serif;
	font-size:10pt;
	margin:10px;
	background-color:white;
	}
p   { font-size:10pt; margin-left:25pt; }
pre { font-size:10pt; margin-left:25pt; 
    border-style:solid;    border-width:thin;}
h3 { -webkit-margin-after: 0px; -webkit-margin-before: 2em; }
h4, h5, ol, ul, dl {margin-left:25pt;}
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table.topmenu {border-style:solid; border-width:0px}
table, th, td { padding: 3px; }
td {text-align:center;}
td.index {text-align:left;}
th.index {text-align:left;}
td.integer {text-align:right;}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>
"""

menu_header="""<!DOCTYPE html>
<html> 
<head> <title>%(name)s</title>
<style type="text/css">
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table.topmenu {border-style:solid; border-width:0px}
table, th, td { padding: 3px; }
td {text-align:center;}
td.index {text-align:left;}
th.index {text-align:left;}

td.integer {text-align:right;}
body{	font-family:verdana,arial,sans-serif; font-size:10pt;	margin:10px;
	background-color:white;	}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>
</head>
""" 
   
dd_menu_header="""<!DOCTYPE html>
<html> 
<head> <title>%(name)s</title>
<link rel="stylesheet" type="text/css" href="%(include)s/flexdropdown.css" />
<style type="text/css">
body{	font-family:verdana,arial,sans-serif; font-size:10pt;	margin:10px;
	background-color:white;	}
h4 {margin-left:10pt; -webkit-margin-after: 0px; -webkit-margin-before: 2em; }
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
<script type="text/javascript" src="%(include)s/flexdropdown.js">
/***********************************************
* Flex Level Drop Down Menu- (c) Dynamic Drive DHTML code library (www.dynamicdrive.com)
* This notice MUST stay intact for legal use
* Visit Dynamic Drive at http://www.dynamicdrive.com/ for this script and 100s more
***********************************************/
</script>
</head>
"""

model_menu_header="""<!DOCTYPE html>\n<html> 
<head> <title>%(model)s index</title>
<style type="text/css">
body{	font-family:verdana,arial,sans-serif;
	font-size:10pt;	margin:10px;
	background-color:white;
	}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>
</head>
<body>
<h2><a href="%(upper_link)s?skipDecoration">%(upper)s</a>%(model)s</h2>"""

top_nav= """<html> 
<head> <title>%(series)s menu</title>
<style type="text/css">
body{	font-family:verdana,arial,sans-serif; font-size:10pt;	margin:10px;
	background-color:white;	}
h4 {margin-left:15pt;}
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table.topmenu {border-style:solid; border-width:0px}
table, th, td { padding: 3px; }
td {text-align:center;}
td.index {text-align:left;}
th.index {text-align:left;}
td.integer {text-align:right;}

a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style> 
</head>
<body>
<h3><a %(skymodels)s/%(upper)s</h3>""" 

mathjax=r"""<script type="text/x-mathjax-config">
     MathJax.Hub.Config({tex2jax: {
      inlineMath: [['$','$'], ["\\(","\\)"]], 
      displayMath: [ ['$$','$$'],["\\[", "\\]"]],
      processEscapes: true},
      TeX: { equationNumbers: {autoNumber: "AMS"}},});
</script>
<script type="text/javascript"
       src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
"""

def table_menu(max_cols=20):
    """Make a table of analyses vx. models
    """
    models = sorted(glob.glob('../*/plots/index.html'), reverse=True)
    mdict={}
    ddict={}
    anames = set()
    for m in models:
        t = os.path.split(m)[0]
        mname = t.split('/')[-2]
        #dirs = filter(lambda f: os.path.isdir(f), glob.glob(t+'/*'))
        #dnames = map(lambda f: f.split('/')[-1], dirs)
        dnames = [x.split('/')[-2] for x in glob.glob(t+'/*/index.html')]
        map( lambda n: anames.add(n), dnames)
        mdict[mname] = dnames
        for d in dnames:
            if d in ddict.keys():
                ddict[d].append(mname)
            else: ddict[d]=[mname]
    idents = sorted(np.array(list(anames ))) 
    cols = sorted(mdict.keys())
    href = lambda m, a: '-A-%s-B-%s-C-' %(m ,a)
    x = pd.DataFrame(np.array([[href(id,a) if id in ddict[a] else '' for id in cols] for a in idents]),
       columns=cols, index=idents)
    if len(cols)>max_cols: 
        #transpose if too many columns
        h=x.T.to_html()
    else:
        h = x.to_html()
    # <a href="%s/plots/%s/index.html?skipDecoration">X</a>'
    h1 =  h.replace('-A-','<a href="').replace('-B-','/plots/')\
           .replace('-C-', '/index.html?skipDecoration">&#10004;</a>')
           
    def replace_model(instring):
        p = re.compile(r'<th>(.+)</th>')
        replacer = lambda m: '<th class="index"><a href="*/plots/index.html?skipDecoration"><strong>*</strong></a></th>'.replace('*',m.group(1)) 
        return p.sub(replacer, instring)
        
    def replace_analysis(instring):
        p = re.compile('<tr>\s*<td>(.+)<')
        replacer = lambda m: '<tr>\n\t<td class="index">'+m.group(1)+'<'
        return p.sub(replacer, instring)

    h2 = replace_model(h1) ##### oops, now they use <th>
    return h2
    # No links to describe analysis yet?
    #h3 = replace_analysis(h2)
    #return h3
    
def header(title=''):
    """ return HTML for start of a document """
    return '<!DOCTYPE html>\n<head><title>%s</title>\n' %title + style + mathjax +'</head>\n'

class DDmenu():
    """ manage a menu document using the Dynamic Drive DHTML code library (www.dynamicdrive.com)
    """
    def __init__(self, name, depth=2):
        self.menuname = name
        self.doc = dd_menu_header%dict(name=name, include='../'*depth+'plot_browser/includes')\
            + '\n<body>'

    def folder(self, name=None, **kw):
        """ must set id, href, text in kw
        """
        self.doc += '\n<h4> <a href="%(href)s" data-flexmenu="%(id)s"> %(text)s</a> </h4>' % kw
        self.doc += '\n <ul id="%(id)s" class="flexdropdownmenu">' % kw

    def item(self, name):
        """ name is full anchor """
        self.doc += '\n  <li>%s</li>' % name
        #print ('added name %s to menu' % name)
        
    def add_menu(self, menu_html, folder):
        """ add a menu from a file created by this class """
        t1 = open(menu_html).read()
        n,m = t1.find('<body'), t1.find('</body>')
        assert n>0 and m>0, 'Parsing problem: n,m=%d %d\n %s' % (n,m, t1)
        t2 = t1[n:m].split('\n')
        t3='\n'.join(t2[1:-1] ) #
        self.doc += '\n'+t3.replace('index.html', os.path.join(folder, 'index.html'))

    def save(self, filename):
        self.doc += '\n </ul>\n</body>'
        open(filename,'w').write(self.doc)

class HTMLindex():
    """ Manage the web browser pages
    """
    def __init__(self):
        w = os.getcwd().split(os.path.sep)
        self.model = w[-1] #'/'.join(w[-2:])
        self.upper = w[-2]+'/'
        menu = DDmenu('%s index'%self.model, depth=3)
        menu.doc +='\n<h2><a href="../../plot_index.html?skipDecoration">%(upper)s</a>%(model)s</h2>'%\
            dict(upper=self.upper, model=self.model)
        plot_folders = [x.split('/')[1] for x in sorted(glob.glob('plots/*/index.html'))]
        for folder in plot_folders:
            menu_html = os.path.join('plots',folder, 'menu.html')
            if os.path.exists(menu_html):
                menu.add_menu(menu_html, folder)
            else:
                menu.doc += '\n <h4><a href="%s?skipDecoration">%s</a></h4>' % (os.path.join(folder,'index.html'), folder)
        self.menu = menu
   
    def _repr_html_(self):    
        return self.menu.doc
    
    def make_config_link(self):
        print ('***Warning: make_config_link should not be called: use config analysis module')
        
    def create_menu(self):
        self.menu.save(os.path.join(os.getcwd(), 'plots/index.html'))
        print ('wrote menu %s' %os.path.join(os.getcwd(), 'plots/index.html'))
        
    def update_top(self, filename='../plot_index.html'):
        def parse_path(x): 
            'return relative path, model name'
            t = x.split(os.path.sep)
            return  '/'.join(t[1:]) , t[1]
        def parse_model(x):
            return '<a href="%s?skipDecoration"> %s </a>' %(parse_path(x) )
        def model_comment(x):
            a,b=parse_path(x)
            try:
                y = '../'+b+'/config.yaml'
                if os.path.exists(y):
                    return yaml.load(open(y).read()).get('comment', 'no comment')
                return eval(open('../'+b+'/config.txt').read()).get('comment', 'no comment')
            except IOError:
                return  'no comment found'
        
        models = sorted(glob.glob('../*/plots/index.html'), reverse=True)
        assert len(models)>0, 'No models found?'
        self.last_model = parse_path(models[0])[0]
        self.skymodels='<a href="../plot_index.html?skipDecoration">skymodels</a>'
        self.series = os.getcwd().split('/')[-2]
        s = top_nav % self.__dict__
        s += '\n<table class="topmenu">'
        for m in models:
            s += '\n  <tr><td valign="top" class="index"><strong>%s</strong></td>'% parse_model(m)
            s += '\n      <td class="index"> %s </td></tr>' % model_comment(m)
        s += '\n</table>\n</body></html>\n'
        
        s += '\n<hr>'
        s += '\n<h3>Table of analyses and models</h3>\n'
        s +=  table_menu()
        open(filename, 'w').write(s)
        print ('wrote top menu %s' % os.path.join(os.getcwd(),filename))
    
    @staticmethod
    def head(title=''):
        return '<head><title>%s</title>\n'+style+'</head>\n'
