"""
Manage the Web page generation
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/_html.py,v 1.7 2013/08/04 14:55:11 burnett Exp $
"""
import os, glob
import pandas as pd

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
h4, h5, ol, ul {margin-left:25pt;}
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table.topmenu {border-style:solid; border-width:0px}
table, th, td { padding: 3px; }
td {text-align:center;}
td.index {text-align:left;}
td.integer {text-align:right;}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>
"""

menu_header="""<!DOCTYPE html>
<html> 
<head> <title>%(name)s</title>
<style type="text/css">
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
h4 {margin-left:15pt;}
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
<head> <title>Top Nav</title>
<style type="text/css">
body{	font-family:verdana,arial,sans-serif; font-size:10pt;	margin:10px;
	background-color:white;	}
h4 {margin-left:15pt;}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style> 
</head>
<body>
<h3>skymodels/%(upper)s</h3>""" 

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
                menu.doc += '\n <li><a href="%s">%s</a></li>' % (os.path.join(folder,'index.html'), folder)
        self.menu = menu
   
    def _repr_html_(self):    
        return self.menu.doc
    
    def make_config_link(self):
        print '***Warning: make_config_link should not be called: use config analysis module'
        
    def create_menu(self):
        self.menu.save(os.path.join(os.getcwd(), 'plots/index.html'))
        print 'wrote menu %s' %os.path.join(os.getcwd(), 'plots/index.html')
        
    def update_top(self, filename='../plot_index.html'):
        def parse_path(x): 
            'return relative path, model name'
            t = x.split(os.path.sep)
            return  '/'.join(t[1:]) , t[1]
        def parse_model(x):
            return '<a href="%s?skipDecoration"> %s </a>' %(parse_path(x) )
        def model_comment(x):
            a,b=parse_path(x)
            return eval(open('../'+b+'/config.txt').read()).get('comment', 'no comment')
        
        models = sorted(glob.glob('../*/plots/index.html'), reverse=True)
        assert len(models)>0, 'No models found?'
        self.last_model = parse_path(models[0])[0]
        s = top_nav % self.__dict__
        s += '\n<table class="topmenu">'
        for m in models:
            s += '\n  <tr><td valign="top" class="index">%s</td>'% parse_model(m)
            s += '\n      <td class="index"> %s </td></tr>' % model_comment(m)
        s += '\n</table>\n</body></html>\n'
        open(filename, 'w').write(s)
        print 'wrote top menu %s' % os.path.join(os.getcwd(),filename)
    
    @staticmethod
    def head(title=''):
        return '<head><title>%s</title>\n'+style+'</head>\n'