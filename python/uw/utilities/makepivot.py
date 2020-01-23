"""
Code to make a pivot collection, using the PivotWeb service at deeptalk.phys.washington.edu

author: Toby Burnett <tburnett@uw.edu>
PivotWeb support: Gordon Watts <gwatts@uw.edu>
"""        
       
import os, glob, argparse, zipfile, json
import urllib, httplib, json
import pandas as pd    
import numpy as np    

http_host, urlbase = "deeptalk.phys.washington.edu", "/PivotWeb/api"
#
class Connection(object):
    """ manage a connection to the deeptalk server 
     get and put methods encapsulate JSON conversion
    """
    def __init__(self):
        self.connection = httplib.HTTPConnection(http_host)

    def request(self, api, verb, params, data = None, dataFormat = None):
        """ function to run the request to the web service. No data sent if there is no data.
        api : string
            one of 'Images', 'Collections', 'Facets', 'FacetCategories'
        verb : string
            'PUT', or 'GET' or 'DELETE'
        params : dict
        data : 
        dataFormat :
        """
        # this code adapted from the module PivotViewCollection by G. Watts
        
        apis=('Images', 'Collections', 'Facets', 'FacetCategories','FacetsDB')
        assert api in apis , 'api "%s" not one of %s'% (api, apis)
        arglist = urllib.urlencode(params)
        if len(arglist) > 0:
            arglist = "?" + arglist
        url = "%s/%s%s" % (urlbase, api, arglist)
        dlen = len(data) if data is not None else 0
        headers = {"Content-Length" : dlen}
        if dataFormat is not None:
            headers["Content-Type"] = dataFormat
        self.connection.request(verb, url, data, headers)
        response = self.connection.getresponse()
        if response.status<200 or response.status >= 300:
            print ('api, headers:', api, headers)
            print ('params:', params)
            print ('response:', response)
            print (response.status, response.reason)
            raise Exception(response.reason)
        textBack = response.read()
        return textBack
        
    def put(self, api, params, data=None, dataFormat=None):
        """ send a PUT request
            unless dataFormat is specified, assume data is an object, convert to json.
            return result, attempt to convert from json
        """
        if data is not None and dataFormat is None:
            data = json.dumps(data)
            dataFormat ='application/json; charset=utf-8'
        req = self.request(api, 'PUT', params, data=data, dataFormat=dataFormat)
        try: 
            return json.loads(req)
        except:
            return req
        
    def get(self, api, params={}):
        """ set a GET request, convert response from json on return
        """
        return json.loads(self.request(api, 'GET', params))


class FacetList(object):
    """ manage a table, extract Facet-style entries for Pivot
    Note that they have to be serialized with JSON, which has limited number of types
    """
    def __init__(self, table) :
        """ 
        table : string
            Expect to be a csv file with the first column the name of the item
        """
        ext = os.path.splitext(table)[-1]
        if ext=='.csv':
            self.df = pd.read_csv(table, index_col=0) 
        elif ext=='.pickle':
            self.df = pd.load(table)
        else:
            raise Exception('table extension not .csv or .pickle')
        self.cols = self.df.columns
        self.index = self.df.index
        self.df['name']=self.index
        # this is the list of types that I've found so far that  pandas.read_csv generates
        self.types = [ dict(float64='Number', float='Number', int64='Number',object='Value', bool='bool')[str(t)] for t in self.df.dtypes]
        print ('Loaded CSV file %s: found %d items, %d columns' % (table, len(self.df), len(self.cols)))

    def __getitem__(self, name):
        """
        name : integer or string
            used to select the named row
            
        returns a list of dictionary entries for the values; raises KeyError or IndexError if not found
        """
        item = self.df.ix[name]
        def interpret(type, v): # convert python object to type that JSON likes
            if type=='bool': return ('Value', str(v))    # convert bool to 'True' or 'False'
            if type=='Number': return (type, np.float(v)) # only pass float type to json serialization
            return (type, v)
        t = [dict( [('Name',x), interpret(type,item[x])]) for x,type in zip(self.cols, self.types)]
        return t
        
    def send(self, collection_name, keyname='name'):
        """Send the table directly, using the FacetsDB API
        Assume that the collection already exists, as a new-style list of images without attached facets.
        
        collection_name : string
            set
        keyname : string
        """
        self.connection = Connection()
        cid= self.connection.put('Collections', dict(cname=collection_name ) )['Id']
        print ('Created or attached to existing collection, name="%s", Id=%d' % (collection_name, cid))
        columns = [ dict(Name=name, Values=list(np.array(values.values)), isKey=name==keyname) for name,values in self.df.T.iterrows()]
        self.connection.put('FacetsDB', dict(cID=cid), data=dict(Name=collection_name, Columns=columns))

    def keylist(self, keyname='name'):
        return self.df[keyname].values
        

class CompositeImage(object):
    def __init__(self, folders, combine):
        """ 
        Used by ImageList to look up an image
              (under development)
        """
        
class ImageList(object):
    """ manage a list of images, either in a zip file or a folder
    """

    def __init__(self, image_folder, ext='jpg'):
        if os.path.exists(image_folder+'.zip'):
            unzipper = zipfile.ZipFile(image_folder+'.zip')
            self.filenames =  unzipper.namelist()
            self.opener = unzipper.open
        elif os.path.exists(image_folder):
            # should check that it is a folder
            self.filenames = glob.glob(image_folder+'/*.%s'%ext)
            self.opener = lambda f: open(f, 'rb')
        else:
            raise Exception('image_folder "%s" is not a folder or a zip file' %image_folder)
        assert len(self.filenames)>0, 'No image files found in %s' %image_folder
            
    def match(self, prefix, filename):
        # adjusted name must be followed by '_' or '.'
        f = os.path.split(filename)[-1]
        if not f.startswith(prefix):
            return False
        return True #f[len(prefix)] in ('_', '.', 't')  # for tsmap now
        
    def __call__(self, itemname):
        """ return the binary image, and the filename
        """
        prefix = itemname.strip().replace(' ', '_').replace('+','p')
        for f in self.filenames:
            if self.match(prefix, f):
                #if f.find(prefix)>=0: 
                return self.opener(f).read(), f
        raise Exception('File "%s" derived from image name "%s" not found in image list' % (prefix,itemname))
    
    def add_image(self,  name): 
        """ 
        Add an image 
        
        name : string
            name of the image: expect to find equivalent filename in the list of images
        """
        data, filename= self(name)
        image_type = os.path.splitext(filename)[-1][1:]
        assert image_type in ('png', 'jpg', 'jpeg'), 'Image type "%s" not accepted'

        info = self.connection.put('Images', dict(cId=self.cId, iName=name ), 
                data=data, dataFormat="image/%s" % image_type)
 
    def send(self, collection_name, keynames=None):
        """Send the list of images to a collection.
        
        collection_name : string
        keynames : string, optional
            if set, expect to be a list of names corresponding to the image names
            otherwise, send them all
        """   
        self.connection = Connection()
        self.cId = self.connection.put('Collections', dict(cname=collection_name ) )['Id']
        print ('Created or attached to existing collection, name="%s", Id=%d' % (collection_name, self.cId))
        
        nsent=0
        for filename in self.filenames:
            t = os.path.split(filename)[-1].split('_')
            name = (t[0] if len(t)==2 else ' '.join(t[:2])).replace('p','+')
            if keynames is not None and name not in keynames: continue
            try:
                self.add_image( name)
                nsent +=1
            except Exception as msg:
                print (msg)
        print ('Sent %d/%d images' % (nsent, len(self.filenames)))
        
class MakeCollection(object):
    """ Use the PivotWeb service to create a Pivot collection from a csv file, and associated folder containing images
    
    """
    def __init__(self, cname, 
                    image_folder,  table,
                    folder='.',
                    refresh=False,
                    image_type='jpg',):
        """
        cname : string
            Name to give to the collection
        image_folder : string
            name of sub folder containing the images, or a file with extension '.zip'.
            If the latter, it will take precedence
        table : string
            Expect to be a csv file (exension '.csv') with the first column the name of the item
            Each such name must be part of the filename for an image file. (No check for duplicates)
            Future: allow other format
        folder : string
            base folder to chdir to; expect the others to be specied relative to this
        refresh : bool
            set True to delete the collection if it already exists. 
        """
        self.cname = cname
        if refresh: delete_collection(name=cname)
        self.connection = Connection()
        res = self.connection.put('Collections', dict(cname=self.cname ) )
        self.cId = res['Id']
        print ('Created or attached to existing collection, name="%s", Id=%d' % (cname, self.cId))
        
        os.chdir(os.path.expanduser(folder))
        self.facets = FacetList(table)
        self.images = ImageList(image_folder, ext=image_type)
        self.fill_from_table()
    
    def image_name(self, filename):
        " parse filename to a source name "
        return os.path.split(filename)[-1].split('_')[0].replace('_',' ').replace('p','+')
            
    def fill_from_table(self):
        """ Fill the collection using all entries in the csv file with corresponding images
        """
        for id,item in enumerate(self.facets):
            try:
                self.add_image( self.facets.index[id] , item)
            except Exception as msg:
                print (msg)
                #raise

    def add_image(self,  name, item): 
        """ 
        Add an image and set or replace the associated Facets
        
        name : string
            name of the image: expect to find equivalent filename in the list of images
        item : dict
            properties of the facets
        """
        data, filename= self.images(name)
        image_type = os.path.splitext(filename)[-1][1:]
        assert image_type in ('png', 'jpg', 'jpeg'), 'Image type "%s" not accepted'

        info = self.connection.put('Images', dict(cId=self.cId, iName=name ), 
                data=data, dataFormat="image/%s" % image_type)
        self.connection.put('Facets', dict(cId=self.cId, iId=info['Id'], replace=True ),
            data = item)

def set_format(cId, format='F3', facets=None, quiet=True):
    """ Set numerical format for all or a subset of the Facets in a collection
    """
    c = Connection()
    fc = c.get('FacetCategories',  dict(cID=cId))
    for f in fc:
        if facets is None or f['Name'] in facets:
            before = f.get('NumberFormat', 'F1') # default
            f['NumberFormat']=format
            if not quiet: print ('change format of Facet "%s": %s->%s' % (f['Name'], before, format))
    c.put('FacetCategories', dict(cId=cId), data=fc)

def list_collections():
    c = Connection()
    t=c.get('Collections')
    print ('  Id  Name')
    for x in t:
        print ('%4d  %s' % (x['Id'], x['Name']))

def delete_collection(id=None, name=None, ):
    if name is not None:
        par = dict(cName=name)
    elif id is not None:
        par = dict(cId=id)
    else:
        raise Exception('Neither name nor id specified')
    c = Connection()
    c.request("Collections", "DELETE", par)	
    print ('Collection %s deleted' %par)
    
def list_facets(id):
    c = Connection()
    t = c.get('FacetCategories', dict(cID=id))
    return pd.DataFrame( t )

def main(args):
    MakeCollection(name=args.args[0], table=args.args[1], folder='.', image_folder='counts_dir')
      
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='make a Pivot table from a csv file. Can be viewed at http://deeptalk.phys.washington.edu/PivotWeb')
    parser.add_argument('args', nargs=2, help='name csv file' )
    args = parser.parse_args()
    main(args)
