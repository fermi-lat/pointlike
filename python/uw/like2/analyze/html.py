"""
Manage the Web page generation
$Header$
"""
import os, glob
import pandas as pd

class HTMLindex():
    """ Manage the web browser pages
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
    border-style:solid;
    border-width:thin;}
h5 {margin-left:25pt;}
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table.topmenu {border-style:solid; border-width:0px}
table, th, td { padding: 3px; }
td {text-align:center;}
td.index {text-align:left;}
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>"""

    menu_header="""<html> <head> <title>%(model)s index</title> %(style)s 
    <script> function load(){ parent.content.location.href='%(model_summary)s';} </script>
    </head>
<body onload="load()">
<h2><a href="%(upper_link)s?skipDecoration">%(upper)s</a>%(model)s</h2>"""

    top_nav="""<html> <head> <title>Top Nav</title> %(style)s 
    <script> function load(){ parent.menu.location.href = '%(last_model)s';} </script>
    </head>
<body onload="load()">
<h3>skymodels/%(upper)s</h3>""" 

    def __init__(self, folder='plots/*'):
        self.style = HTMLindex.style
    
        w= glob.glob(folder)
        if len(w)==0: 
            print 'Did not find any plot folders under %s' % folder
        z = dict( zip(w, [glob.glob(a+'/*.htm*') for a in w] ) )
        w = os.getcwd().split('/')
        self.upper = w[-2]+'/'
        self.upper_link = '../../plot_index.html'
        self.model = w[-1] #'/'.join(w[-2:])
        self.model_summary='plots/config/index.html'
        s= HTMLindex.menu_header % self.__dict__
        
        def parse_item(x):
            head, tail =os.path.split(x)
            name = os.path.splitext(tail)[0]
            n = name.find('_uw')
            #note the special qualifier for use with the SLAC decorator
            return '<a href="%s?skipDecoration" target="content">%s</a><br>' % (x,name[:n])

        for k in sorted(z.keys()):
            v = z[k]
            if len(v)==0: continue
            index = '%s/index.html'%k
            if index in v:
                v.remove(index)
                s += '\n<h4><a href="%s?skipDecoration" target="content">%s</a></h4>'% (index, k.split('/')[-1])
            else:
                s += '\n<h4>%s</h4>'% k.split('/')[-1]
            s += '\n\t<p>' + '\n\t'.join(map(parse_item, v)) 
        self.ul = s + '</p>\n</body>'
        self.make_config_link()
        
    def _repr_html_(self):    return self.ul
    
    def make_config_link(self):
        html = """<head>%s</head><body><h2><a href="../index.html?skipDecoration">%s</a> - configuration and analysis history files</h2>
        """ %( self.style,self.model)
        for filename in ('config.txt', 'dataset.txt', 'converge.txt', 'summary_log.txt'):
            if not os.path.exists(filename): continue
            html += '<h4>%s</h4>\n<pre>%s</pre>' % (filename, open(filename).read())
        html += '\n</body>'
        if not os.path.exists('plots/config'): os.makedirs('plots/config')
        open('plots/config/index.html','w').write(html)
        print 'wrote plots/config/index.html'
        
    def create_menu(self, filename='plot_index.html'):
        ###summary = open(filename, 'w')
        open(filename, 'w').write(self.ul)
        print 'wrote menu %s' % os.path.join(os.getcwd(),filename)
        # make separate menu for the Decorator browser
        t = self.ul.replace('plots/', '')
        open('plots/index.html', 'w').write(t)
        print 'wrote menu %s' %os.path.join(os.getcwd(), 'plots/index.html')

    def update_top(self, filename='../plot_index.html'):
        def parse_path(x): 
            'return relative path, model name'
            t = x.split('/')
            return  '/'.join(t[1:]) , t[1]
        def parse_model(x):
            return '<a href="%s?skipDecoration"> %s </a>' %(parse_path(x) )
        def model_comment(x):
            a,b=parse_path(x)
            return eval(open('../'+b+'/config.txt').read()).get('comment', 'no comment')
        
        models = sorted(glob.glob('../*/plots/index.html'), reverse=True)
        assert len(models)>0, 'No models found?'
        self.last_model = parse_path(models[0])[0]
        s = HTMLindex.top_nav % self.__dict__
        s += '\n<table class="topmenu">'
        for m in models:
            s += '\n  <tr><td valign="top">%s</td>'% parse_model(m)
            s += '\n      <td> %s </td></tr>' % model_comment(m)
        s += '\n</table>\n</body></html>\n'
        open(filename, 'w').write(s)
        print 'wrote top menu %s' % os.path.join(os.getcwd(),filename)
    @staticmethod
    def head(title=''):
        return '<head><title>%s</title>\n'+HTMLindex.style+'</head>\n'