# Copyright (C) 2011-2013 WikiTeam
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see &lt;http://www.gnu.org/licenses/&gt;.

#######################################################################
# dumpgenerator.py is a script to generate backups of MediaWiki wikis #
# To learn more, read the documentation:                              #
#        http://code.google.com/p/wikiteam/wiki/NewTutorial           #
#######################################################################

import cookielib
import cPickle
import datetime
import getopt
try:
    from hashlib import md5
except ImportError:             # Python 2.4 compatibility
    from md5 import new as md5
import os
import re
import subprocess
import sys
import time
import urllib
import urllib2

def truncateFilename(other={}, filename=''):
    """ Truncate filenames when downloading images with large filenames """
    return filename[:other['filenamelimit']] + md5(filename).hexdigest() + '.' + filename.split('.')[-1]

def delay(config={}):
    """ Add a delay if configured for that """
    if config['delay'] &gt; 0:
        print 'Sleeping... %d seconds...' % (config['delay'])
        time.sleep(config['delay'])

def cleanHTML(raw=''):
    """ Extract only the real wiki content and remove rubbish """
    """ This function is ONLY used to retrieve page titles and file names when no API is available """
    """ DO NOT use this function to extract page content """
    #different "tags" used by different MediaWiki versions to mark where starts and ends content
    if re.search('&lt;!-- bodytext --&gt;', raw):
        raw = raw.split('&lt;!-- bodytext --&gt;')[1].split('&lt;!-- /bodytext --&gt;')[0]
    elif re.search('&lt;!-- start content --&gt;', raw):
        raw = raw.split('&lt;!-- start content --&gt;')[1].split('&lt;!-- end content --&gt;')[0]
    elif re.search('&lt;!-- Begin Content Area --&gt;', raw):
        raw = raw.split('&lt;!-- Begin Content Area --&gt;')[1].split('&lt;!-- End Content Area --&gt;')[0]
    elif re.search('&lt;!-- content --&gt;', raw):
        raw = raw.split('&lt;!-- content --&gt;')[1].split('&lt;!-- mw_content --&gt;')[0]
    elif re.search('&lt;article id="WikiaMainContent" class="WikiaMainContent"&gt;', raw):
        raw = raw.split('&lt;article id="WikiaMainContent" class="WikiaMainContent"&gt;')[1].split('&lt;/article&gt;')[0]
    else:
        print raw[:250]
        print 'This wiki doesn\'t use marks to split content'
        sys.exit()
    return raw

def getNamespacesScraper(config={}):
    """ Hackishly gets the list of namespaces names and ids from the dropdown in the HTML of Special:AllPages """
    """ Function called if no API is available """
    namespaces = config['namespaces']
    namespacenames = {0:''} # main is 0, no prefix
    if namespaces:
        req = urllib2.Request(url=config['index'], data=urllib.urlencode({'title': 'Special:Allpages', }), headers={'User-Agent': getUserAgent()})
        f = urllib2.urlopen(req)
        raw = f.read()
        f.close()
        
        m = re.compile(r'&lt;option [^&gt;]*?value="(?P&lt;namespaceid&gt;\d+)"[^&gt;]*?&gt;(?P&lt;namespacename&gt;[^&lt;]+)&lt;/option&gt;').finditer(raw) # [^&gt;]*? to include selected="selected"
        if 'all' in namespaces:
            namespaces = []
            for i in m:
                namespaces.append(int(i.group("namespaceid")))
                namespacenames[int(i.group("namespaceid"))] = i.group("namespacename")
        else:
            #check if those namespaces really exist in this wiki
            namespaces2 = []
            for i in m:
                if int(i.group("namespaceid")) in namespaces:
                    namespaces2.append(int(i.group("namespaceid")))
                    namespacenames[int(i.group("namespaceid"))] = i.group("namespacename")
            namespaces = namespaces2
    else:
        namespaces = [0]
    
    namespaces = list(set(namespaces)) #uniques
    print '%d namespaces found' % (len(namespaces))
    return namespaces, namespacenames
    
def getNamespacesAPI(config={}):
    """ Uses the API to get the list of namespaces names and ids """
    namespaces = config['namespaces']
    namespacenames = {0:''} # main is 0, no prefix
    if namespaces:
        req = urllib2.Request(url=config['api'], data=urllib.urlencode({'action': 'query', 'meta': 'siteinfo', 'siprop': 'namespaces', 'format': 'xml'}), headers={'User-Agent': getUserAgent()})
        f = urllib2.urlopen(req)
        raw = f.read()
        f.close()
            
        m = re.compile(r'&lt;ns id="(?P&lt;namespaceid&gt;\d+)"[^&gt;]*?/?&gt;(?P&lt;namespacename&gt;[^&lt;]+)?(&lt;/ns&gt;)?').finditer(raw) # [^&gt;]*? to include case="first-letter" canonical= etc.
        if 'all' in namespaces:
            namespaces = []
            for i in m:
                namespaces.append(int(i.group("namespaceid")))
                namespacenames[int(i.group("namespaceid"))] = i.group("namespacename")
        else:
            #check if those namespaces really exist in this wiki
            namespaces2 = []
            for i in m:
                if int(i.group("namespaceid")) in namespaces:
                    namespaces2.append(int(i.group("namespaceid")))
                    namespacenames[int(i.group("namespaceid"))] = i.group("namespacename")
            namespaces = namespaces2
    else:
        namespaces = [0]
    
    namespaces = list(set(namespaces)) #uniques
    print '%d namespaces found' % (len(namespaces))
    return namespaces, namespacenames

def getPageTitlesAPI(config={}):
    """ Uses the API to get the list of page titles """
    titles = []
    namespaces, namespacenames = getNamespacesAPI(config=config)
    for namespace in namespaces:
        if namespace in config['exnamespaces']:
            print '    Skipping namespace =', namespace
            continue
        
        c = 0
        print '    Retrieving titles in the namespace %d' % (namespace)
        headers = {'User-Agent': getUserAgent()}
        apfrom = '!'
        while apfrom:
            sys.stderr.write('.') #progress
            params = {'action': 'query', 'list': 'allpages', 'apnamespace': namespace, 'apfrom': apfrom, 'format': 'xml', 'aplimit': 500}
            data = urllib.urlencode(params)
            req = urllib2.Request(url=config['api'], data=data, headers=headers)
            try:
                f = urllib2.urlopen(req)
            except:
                try:
                    print 'Server is slow... Waiting some seconds and retrying...'
                    time.sleep(10)
                    f = urllib2.urlopen(req)
                except:
                    print 'An error has occurred while retrieving page titles with API'
                    print 'Please, resume the dump, --resume'
                    sys.exit()
            xml = f.read()
            f.close()
            m = re.findall(r'&lt;allpages (?:apfrom|apcontinue)="([^&gt;]+)" /&gt;', xml)
            if m:
                apfrom = undoHTMLEntities(text=m[0]) #&amp;quot; = ", etc
            else:
                apfrom = ''
            m = re.findall(r'title="([^&gt;]+)" /&gt;', xml)
            titles += [undoHTMLEntities(title) for title in m]
            c += len(m)
        print '    %d titles retrieved in the namespace %d' % (c, namespace)
    return titles

def getPageTitlesScraper(config={}):
    """  """
    titles = []
    namespaces, namespacenames = getNamespacesScraper(config=config)
    for namespace in namespaces:
        print '    Retrieving titles in the namespace', namespace
        url = '%s?title=Special:Allpages&amp;namespace=%s' % (config['index'], namespace)
        raw = urllib.urlopen(url).read()
        raw = cleanHTML(raw)
        
        r_title = r'title="(?P&lt;title&gt;[^&gt;]+)"&gt;'
        r_suballpages = ''
        r_suballpages1 = r'&amp;amp;from=(?P&lt;from&gt;[^&gt;]+)&amp;amp;to=(?P&lt;to&gt;[^&gt;]+)"&gt;'
        r_suballpages2 = r'Special:Allpages/(?P&lt;from&gt;[^&gt;]+)"&gt;'
        if re.search(r_suballpages1, raw):
            r_suballpages = r_suballpages1
        elif re.search(r_suballpages2, raw):
            r_suballpages = r_suballpages2
        else:
            pass #perhaps no subpages
        
        deep = 3 # 3 is the current deep of English Wikipedia for Special:Allpages, 3 levels
        c = 0
        checked_suballpages = []
        rawacum = raw
        while r_suballpages and re.search(r_suballpages, raw) and c &lt; deep:
            #load sub-Allpages
            m = re.compile(r_suballpages).finditer(raw)
            for i in m:
                fr = i.group('from')
                
                if r_suballpages == r_suballpages1:
                    to = i.group('to')
                    name = '%s-%s' % (fr, to)
                    url = '%s?title=Special:Allpages&amp;namespace=%s&amp;from=%s&amp;to=%s' % (config['index'], namespace, fr, to) #do not put urllib.quote in fr or to
                elif r_suballpages == r_suballpages2: #fix, esta regexp no carga bien todas? o falla el r_title en este tipo de subpag? (wikiindex)
                    fr = fr.split('&amp;amp;namespace=')[0] #clean &amp;amp;namespace=\d, sometimes happens
                    name = fr
                    url = '%s?title=Special:Allpages/%s&amp;namespace=%s' % (config['index'], name, namespace)
                
                if not name in checked_suballpages:
                    checked_suballpages.append(name) #to avoid reload dupe subpages links
                    raw2 = urllib.urlopen(url).read()
                    raw2 = cleanHTML(raw2)
                    rawacum += raw2 #merge it after removed junk
                    print '    Reading', name, len(raw2), 'bytes', len(re.findall(r_suballpages, raw2)), 'subpages', len(re.findall(r_title, raw2)), 'pages'
            c += 1
        
        c = 0
        m = re.compile(r_title).finditer(rawacum)
        for i in m:
            if not i.group('title').startswith('Special:'):
                if not i.group('title') in titles:
                    titles.append(undoHTMLEntities(text=i.group('title')))
                    c += 1
        print '    %d titles retrieved in the namespace %d' % (c, namespace)
    return titles

def getPageTitles(config={}):
    """  """
    #http://en.wikipedia.org/wiki/Special:AllPages
    #http://archiveteam.org/index.php?title=Special:AllPages
    #http://www.wikanda.es/wiki/Especial:Todas
    print 'Loading page titles from namespaces = %s' % (config['namespaces'] and ','.join([str(i) for i in config['namespaces']]) or 'None')
    print 'Excluding titles from namespaces = %s' % (config['exnamespaces'] and ','.join([str(i) for i in config['exnamespaces']]) or 'None')
    
    titles = []
    if config['api']:
        titles = getPageTitlesAPI(config=config)
    elif config['index']:
        titles = getPageTitlesScraper(config=config)
    
    titles = list(set(titles)) #removing dupes (e.g. in CZ appears Widget:AddThis two times (main namespace and widget namespace))
    titles.sort() #sorting
    
    print '%d page titles loaded' % (len(titles))
    return titles

def getXMLHeader(config={}):
    """ Retrieve a random page to extract XML headers (namespace info, etc) """
    #get the header of a random page, to attach it in the complete XML backup
    #similar to: &lt;mediawiki xmlns="http://www.mediawiki.org/xml/export-0.3/" xmlns:x....
    randomtitle = 'Main_Page' #previously AMF5LKE43MNFGHKSDMRTJ
    xml = getXMLPage(config=config, title=randomtitle, verbose=False)
    header = xml.split('&lt;/mediawiki&gt;')[0]
    if not xml:
        print 'XML export on this wiki is broken, quitting.'
        sys.exit()
    return header

def getXMLFileDesc(config={}, title=''):
    """ Get XML for image description page """
    config['curonly'] = 1 #tricky to get only the most recent desc
    return getXMLPage(config=config, title=title, verbose=False)

def getUserAgent():
    """ Return a cool user-agent to hide Python user-agent """
    useragents = ['Mozilla/5.0 (X11; Linux i686; rv:24.0) Gecko/20100101 Firefox/24.0']
    return useragents[0]

def logerror(config={}, text=''):
    """ Log error in file """
    if text:
        f = open('%s/errors.log' % (config['path']), 'a')
        f.write('%s: %s\n' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), text))
        f.close()

def getXMLPageCore(headers={}, params={}, config={}):
    """  """
    #returns a XML containing params['limit'] revisions (or current only), ending in &lt;/mediawiki&gt;
    #if retrieving params['limit'] revisions fails, returns a current only version
    #if all fail, returns the empty string
    xml = ''
    c = 0
    maxseconds = 100 #max seconds to wait in a single sleeping
    maxretries = 5 # x retries and skip
    increment = 20 #increment every retry
    while not re.search(r'&lt;/mediawiki&gt;', xml):
        if c &gt; 0 and c &lt; maxretries:
            wait = increment * c &lt; maxseconds and increment * c or maxseconds # incremental until maxseconds
            print '    XML for "%s" is wrong. Waiting %d seconds and reloading...' % (params['pages'], wait)
            time.sleep(wait)
            if params['limit'] &gt; 1: # reducing server load requesting smallest chunks (if curonly then limit = 1 from mother function)
                params['limit'] = params['limit'] / 2 # half
        if c &gt;= maxretries:
            print '    We have retried %d times' % (c)
            print '    MediaWiki error for "%s", network error or whatever...' % (params['pages'])
            # If it's not already what we tried: our last chance, preserve only the last revision...
            # config['curonly'] means that the whole dump is configured to save nonly the last
            # params['curonly'] should mean that we've already tried this fallback, because it's set by the following if and passed to getXMLPageCore
            if not config['curonly']: 
                print '    Trying to save only the last revision for this page...'
                params['curonly'] = 1
                logerror(config=config, text='Error while retrieving the full history of "%s". Trying to save only the last revision for this page' % (params['pages']))
                return getXMLPageCore(headers=headers, params=params, config=config)
            else:
                print '    Saving in the errors log, and skipping...'
                logerror(config=config, text='Error while retrieving the last revision of "%s". Skipping.' % (params['pages']))
                return '' # empty xml
        
        data = urllib.urlencode(params)
        req = urllib2.Request(url=config['index'], data=data, headers=headers)
        try:
            f = urllib2.urlopen(req)
        except:
            try:
                print 'Server is slow... Waiting some seconds and retrying...'
                time.sleep(15)
                f = urllib2.urlopen(req)
            except:
                print 'An error have occurred while retrieving "%s"' % (params['pages'])
                print 'Please, resume the dump, --resume'
                sys.exit()
                # The error is usually temporary, but we exit the dump altogether.
        xml = f.read()
        c += 1
    
    return xml

def getXMLPage(config={}, title='', verbose=True):
    """  """
    #return the full history (or current only) of a page
    #if server errors occurs while retrieving the full page history, it may return [oldest OK versions] + last version, excluding middle revisions, so it would be partialy truncated
    #http://www.mediawiki.org/wiki/Manual_talk:Parameters_to_Special:Export#Parameters_no_longer_in_use.3F
    
    limit = 1000
    truncated = False
    title_ = title
    title_ = re.sub(' ', '_', title_)
    #do not convert &amp; into %26, title_ = re.sub('&amp;', '%26', title_)
    headers = {'User-Agent': getUserAgent()}
    params = {'title': 'Special:Export', 'pages': title_, 'action': 'submit', }
    if config['curonly']:
        params['curonly'] = 1
        params['limit'] = 1
    else:
        params['offset'] = '1' # 1 always &lt; 2000s
        params['limit'] = limit
    if config.has_key('templates') and config['templates']: #in other case, do not set params['templates']
        params['templates'] = 1
    
    xml = getXMLPageCore(headers=headers, params=params, config=config)

    #if complete history, check if this page history has &gt; limit edits, if so, retrieve all using offset if available
    #else, warning about Special:Export truncating large page histories
    r_timestamp = r'&lt;timestamp&gt;([^&lt;]+)&lt;/timestamp&gt;'
    if not config['curonly'] and re.search(r_timestamp, xml): # search for timestamps in xml to avoid analysing empty pages like Special:Allpages and the random one
        while not truncated and params['offset']: #next chunk
            params['offset'] = re.findall(r_timestamp, xml)[-1] #get the last timestamp from the acum XML
            xml2 = getXMLPageCore(headers=headers, params=params, config=config)
            
            if re.findall(r_timestamp, xml2): #are there more edits in this next XML chunk or no &lt;page&gt;&lt;/page&gt;?
                if re.findall(r_timestamp, xml2)[-1] == params['offset']:
                    #again the same XML, this wiki does not support params in Special:Export, offer complete XML up to X edits (usually 1000)
                    print 'ATTENTION: This wiki does not allow some parameters in Special:Export, therefore pages with large histories may be truncated'
                    truncated = True
                    break
                else:
                    """    &lt;/namespaces&gt;
                      &lt;/siteinfo&gt;
                      &lt;page&gt;
                        &lt;title&gt;Main Page&lt;/title&gt;
                        &lt;id&gt;15580374&lt;/id&gt;
                        &lt;restrictions&gt;edit=sysop:move=sysop&lt;/restrictions&gt; (?)
                        &lt;revision&gt;
                          &lt;id&gt;418009832&lt;/id&gt;
                          &lt;timestamp&gt;2011-03-09T19:57:06Z&lt;/timestamp&gt;
                          &lt;contributor&gt;
                    """
                    #offset is OK in this wiki, merge with the previous chunk of this page history and continue
                    xml = xml.split('&lt;/page&gt;')[0] + '    &lt;revision&gt;' + ('&lt;revision&gt;'.join(xml2.split('&lt;revision&gt;')[1:]))
            else:
                params['offset'] = '' #no more edits in this page history
    
    if verbose:
    	numberofedits = len(re.findall(r_timestamp, xml))
    	if (numberofedits == 1):
    		print '    %s, %s edit' % (title, numberofedits)
    	else:
	        print '    %s, %s edits' % (title, numberofedits)
    
    return xml

def cleanXML(xml=''):
    """ Trim redundant info """
    #do not touch XML codification, leave AS IS
    if re.search(r'&lt;/siteinfo&gt;\n', xml) and re.search(r'&lt;/mediawiki&gt;', xml):
        xml = xml.split('&lt;/siteinfo&gt;\n')[1]
        xml = xml.split('&lt;/mediawiki&gt;')[0]
    return xml

def generateXMLDump(config={}, titles=[], start=''):
    """  """
    print 'Retrieving the XML for every page from "%s"' % (start and start or 'start')
    header = getXMLHeader(config=config)
    footer = '&lt;/mediawiki&gt;\n' #new line at the end
    xmlfilename = '%s-%s-%s.xml' % (domain2prefix(config=config), config['date'], config['curonly'] and 'current' or 'history')
    xmlfile = ''
    lock = True
    if start:
        #remove the last chunk of xml dump (it is probably incomplete)
        xmlfile = open('%s/%s' % (config['path'], xmlfilename), 'r')
        xmlfile2 = open('%s/%s2' % (config['path'], xmlfilename), 'w')
        prev = ''
        c = 0
        for l in xmlfile:
            #removing &lt;page&gt;\n until end of file
            if c != 0: #lock to avoid write an empty line at the begining of file
                if not re.search(r'&lt;title&gt;%s&lt;/title&gt;' % (start), l): 
                    xmlfile2.write(prev)
                else:
                    break
            c += 1
            prev = l
        xmlfile.close()
        xmlfile2.close()
        #subst xml with xml2
        os.remove('%s/%s' % (config['path'], xmlfilename)) #remove previous xml dump
        os.rename('%s/%s2' % (config['path'], xmlfilename), '%s/%s' % (config['path'], xmlfilename)) #move correctly truncated dump to its real name
    else:
        #requested complete xml dump
        lock = False
        xmlfile = open('%s/%s' % (config['path'], xmlfilename), 'w')
        xmlfile.write(header)
        xmlfile.close()
    
    xmlfile = open('%s/%s' % (config['path'], xmlfilename), 'a')
    c = 1
    for title in titles:
        if not title.strip():
            continue
        if title == start: #start downloading from start, included
            lock = False
        if lock:
            continue
        delay(config=config)
        if c % 10 == 0:
            print 'Downloaded %d pages' % (c)
        xml = getXMLPage(config=config, title=title)
        xml = cleanXML(xml=xml)
        if not xml:
            logerror(config=config, text='The page "%s" was missing in the wiki (probably deleted)' % (title))
        #here, XML is a correct &lt;page&gt; &lt;/page&gt; chunk or 
        #an empty string due to a deleted page (logged in errors log) or
        #an empty string due to an error while retrieving the page from server (logged in errors log)
        xmlfile.write(xml)
        c += 1
    xmlfile.write(footer)
    xmlfile.close()
    print 'XML dump saved at...', xmlfilename

def saveTitles(config={}, titles=[]):
    """ Save title list in a file """
    #save titles in a txt for resume if needed
    titlesfilename = '%s-%s-titles.txt' % (domain2prefix(config=config), config['date'])
    titlesfile = open('%s/%s' % (config['path'], titlesfilename), 'w')
    titlesfile.write('\n'.join(titles))
    titlesfile.write('\n--END--')
    titlesfile.close()
    print 'Titles saved at...', titlesfilename

def saveImageFilenamesURL(config={}, images=[]):
    """ Save image list in a file """
    #save list of images and their urls
    imagesfilename = '%s-%s-images.txt' % (domain2prefix(config=config), config['date'])
    imagesfile = open('%s/%s' % (config['path'], imagesfilename), 'w')
    imagesfile.write('\n'.join(['%s\t%s\t%s' % (filename, url, uploader) for filename, url, uploader in images]))
    imagesfile.write('\n--END--')
    imagesfile.close()
    print 'Image filenames and URLs saved at...', imagesfilename

def getImageFilenamesURL(config={}):
    """ Retrieve file list: filename, url, uploader """
    print 'Retrieving image filenames'
    r_next = r'(?&lt;!&amp;amp;dir=prev)&amp;amp;offset=(?P&lt;offset&gt;\d+)&amp;amp;' # (?&lt;! http://docs.python.org/library/re.html
    images = []
    offset = '29990101000000' #january 1, 2999
    limit = 5000
    retries = 5
    while offset:
        #5000 overload some servers, but it is needed for sites like this with no next links http://www.memoryarchive.org/en/index.php?title=Special:Imagelist&amp;sort=byname&amp;limit=50&amp;wpIlMatch=
        req = urllib2.Request(url=config['index'], data=urllib.urlencode({'title': 'Special:Imagelist', 'limit': limit, 'offset': offset, }), headers={'User-Agent': getUserAgent()})
        f = urllib2.urlopen(req)
        raw = f.read()
        f.close()
        if re.search(ur'(?i)(allowed memory size of \d+ bytes exhausted|Call to a member function getURL)', raw): # delicate wiki
            if limit &gt; 10:
                print 'Error: listing %d images in a chunk is not possible, trying tiny chunks' % (limit)
                limit = limit/10
                continue
            elif retries &gt; 0: # waste retries, then exit
                retries -= 1
                print 'Retrying...'
                continue
            else:
                print 'No more retries, exit...'
                break
        
        raw = cleanHTML(raw)
        #archiveteam 1.15.1 &lt;td class="TablePager_col_img_name"&gt;&lt;a href="/index.php?title=File:Yahoovideo.jpg" title="File:Yahoovideo.jpg"&gt;Yahoovideo.jpg&lt;/a&gt; (&lt;a href="/images/2/2b/Yahoovideo.jpg"&gt;file&lt;/a&gt;)&lt;/td&gt;
        #wikanda 1.15.5 &lt;td class="TablePager_col_img_user_text"&gt;&lt;a href="/w/index.php?title=Usuario:Fernandocg&amp;amp;action=edit&amp;amp;redlink=1" class="new" title="Usuario:Fernandocg (pÃ¡gina no existe)"&gt;Fernandocg&lt;/a&gt;&lt;/td&gt;
        r_images1 = r'(?im)&lt;td class="TablePager_col_img_name"&gt;&lt;a href[^&gt;]+title="[^:&gt;]+:(?P&lt;filename&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;a href="(?P&lt;url&gt;[^&gt;]+/[^&gt;/]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;/td&gt;\s*&lt;td class="TablePager_col_img_user_text"&gt;&lt;a[^&gt;]+&gt;(?P&lt;uploader&gt;[^&lt;]+)&lt;/a&gt;&lt;/td&gt;'
        #wikijuegos 1.9.5 http://softwarelibre.uca.es/wikijuegos/Especial:Imagelist old mediawiki version
        r_images2 = r'(?im)&lt;td class="TablePager_col_links"&gt;&lt;a href[^&gt;]+title="[^:&gt;]+:(?P&lt;filename&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;a href="(?P&lt;url&gt;[^&gt;]+/[^&gt;/]+)"&gt;[^&lt;]+&lt;/a&gt;&lt;/td&gt;\s*&lt;td class="TablePager_col_img_timestamp"&gt;[^&lt;]+&lt;/td&gt;\s*&lt;td class="TablePager_col_img_name"&gt;[^&lt;]+&lt;/td&gt;\s*&lt;td class="TablePager_col_img_user_text"&gt;&lt;a[^&gt;]+&gt;(?P&lt;uploader&gt;[^&lt;]+)&lt;/a&gt;&lt;/td&gt;'
        #gentoowiki 1.18 &lt;tr&gt;&lt;td class="TablePager_col_img_timestamp"&gt;18:15, 3 April 2011&lt;/td&gt;&lt;td class="TablePager_col_img_name"&gt;&lt;a href="/wiki/File:Asus_eeepc-1201nl.png" title="File:Asus eeepc-1201nl.png"&gt;Asus eeepc-1201nl.png&lt;/a&gt; (&lt;a href="/w/images/2/2b/Asus_eeepc-1201nl.png"&gt;file&lt;/a&gt;)&lt;/td&gt;&lt;td class="TablePager_col_thumb"&gt;&lt;a href="/wiki/File:Asus_eeepc-1201nl.png" class="image"&gt;&lt;img alt="" src="/w/images/thumb/2/2b/Asus_eeepc-1201nl.png/180px-Asus_eeepc-1201nl.png" width="180" height="225" /&gt;&lt;/a&gt;&lt;/td&gt;&lt;td class="TablePager_col_img_size"&gt;37 KB&lt;/td&gt;&lt;td class="TablePager_col_img_user_text"&gt;&lt;a href="/w/index.php?title=User:Yannails&amp;amp;action=edit&amp;amp;redlink=1" class="new" title="User:Yannails (page does not exist)"&gt;Yannails&lt;/a&gt;&lt;/td&gt;&lt;td class="TablePager_col_img_description"&gt;&amp;#160;&lt;/td&gt;&lt;td class="TablePager_col_count"&gt;1&lt;/td&gt;&lt;/tr&gt;
        r_images3 = r'(?im)&lt;td class="TablePager_col_img_name"&gt;&lt;a[^&gt;]+title="[^:&gt;]+:(?P&lt;filename&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;a href="(?P&lt;url&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;/td&gt;&lt;td class="TablePager_col_thumb"&gt;&lt;a[^&gt;]+&gt;&lt;img[^&gt;]+&gt;&lt;/a&gt;&lt;/td&gt;&lt;td class="TablePager_col_img_size"&gt;[^&lt;]+&lt;/td&gt;&lt;td class="TablePager_col_img_user_text"&gt;&lt;a[^&gt;]+&gt;(?P&lt;uploader&gt;[^&lt;]+)&lt;/a&gt;&lt;/td&gt;'
        #http://www.memoryarchive.org/en/index.php?title=Special:Imagelist&amp;sort=byname&amp;limit=50&amp;wpIlMatch=
        #(&lt;a href="/en/Image:109_0923.JPG" title="Image:109 0923.JPG"&gt;desc&lt;/a&gt;) &lt;a href="/en/upload/c/cd/109_0923.JPG"&gt;109 0923.JPG&lt;/a&gt; . . 885,713 bytes . . &lt;a href="/en/User:Bfalconer" title="User:Bfalconer"&gt;Bfalconer&lt;/a&gt; . . 18:44, 17 November 2005&lt;br /&gt;
        r_images4 = r'(?im)&lt;a href=[^&gt;]+ title="[^:&gt;]+:(?P&lt;filename&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;a href="(?P&lt;url&gt;[^&gt;]+)"&gt;[^&lt;]+&lt;/a&gt;[^&lt;]+&lt;a[^&gt;]+&gt;(?P&lt;uploader&gt;[^&lt;]+)&lt;/a&gt;'
        m = []
        #different mediawiki versions
        if re.search(r_images1, raw):
            m = re.compile(r_images1).finditer(raw)
        elif re.search(r_images2, raw):
            m = re.compile(r_images2).finditer(raw)
        elif re.search(r_images3, raw):
            m = re.compile(r_images3).finditer(raw)
        elif re.search(r_images4, raw):
            m = re.compile(r_images4).finditer(raw)
        
        for i in m:
            url = i.group('url')
            if url[0] == '/' or (not url.startswith('http://') and not url.startswith('https://')): #is it a relative URL?
                if url[0] == '/': #slash is added later
                    url = url[1:]
                domainalone = config['index'].split('://')[1].split('/')[0] #remove from :// (http or https) until the first / after domain
                url = '%s://%s/%s' % (config['index'].split('://')[0], domainalone, url) # concat http(s) + domain + relative url
            url = undoHTMLEntities(text=url)
            #url = urllib.unquote(url) #do not use unquote with url, it break some urls with odd chars
            url = re.sub(' ', '_', url)
            filename = re.sub('_', ' ', i.group('filename'))
            filename = undoHTMLEntities(text=filename)
            filename = urllib.unquote(filename)
            uploader = re.sub('_', ' ', i.group('uploader'))
            uploader = undoHTMLEntities(text=uploader)
            uploader = urllib.unquote(uploader)
            images.append([filename, url, uploader])
            #print filename, url
        
        if re.search(r_next, raw):
            offset = re.findall(r_next, raw)[0]
            retries += 5 # add more retries if we got a page with offset
        else:
            offset = ''
    
    print '    Found %d images' % (len(images))
    images.sort()
    return images

def getImageFilenamesURLAPI(config={}):
    """ Retrieve file list: filename, url, uploader """
    print 'Retrieving image filenames'
    headers = {'User-Agent': getUserAgent()}
    aifrom = '!'
    images = []
    while aifrom:
        sys.stderr.write('.') #progress
        params = {'action': 'query', 'list': 'allimages', 'aiprop': 'url|user', 'aifrom': aifrom, 'format': 'xml', 'ailimit': 500}
        data = urllib.urlencode(params)
        req = urllib2.Request(url=config['api'], data=data, headers=headers)
        try:
            f = urllib2.urlopen(req)
        except:
            try:
                print 'Server is slow... Waiting some seconds and retrying...'
                time.sleep(10)
                f = urllib2.urlopen(req)
            except:
                print 'An error has occurred while retrieving page titles with API'
                print 'Please, resume the dump, --resume'
                sys.exit()
        xml = f.read()
        f.close()
        m = re.findall(r'&lt;allimages aifrom="([^&gt;]+)" /&gt;', xml)
        if m:
            aifrom = undoHTMLEntities(text=m[0]) #&amp;quot; = ", etc
        else:
            aifrom = ''
        m = re.compile(r'(?im)&lt;img name="(?P&lt;filename&gt;[^"]+)"[^&gt;]*user="(?P&lt;uploader&gt;[^"]+)"[^&gt;]* url="(?P&lt;url&gt;[^"]+)"[^&gt;]*/&gt;').finditer(xml) # Retrieves a filename, uploader, url triple from the name, user, url field of the xml line; space before url needed to avoid getting the descriptionurl field instead.
        for i in m:
            url = i.group('url')
            if url[0] == '/' or (not url.startswith('http://') and not url.startswith('https://')): #is it a relative URL?
                if url[0] == '/': #slash is added later
                    url = url[1:]
                domainalone = config['index'].split('://')[1].split('/')[0] #remove from :// (http or https) until the first / after domain
                url = '%s://%s/%s' % (config['index'].split('://')[0], domainalone, url) # concat http(s) + domain + relative url
            url = undoHTMLEntities(text=url)
            #url = urllib.unquote(url) #do not use unquote with url, it break some urls with odd chars
            url = re.sub(' ', '_', url)
            filename = re.sub('_', ' ', i.group('filename'))
            filename = undoHTMLEntities(text=filename)
            filename = urllib.unquote(filename)
            uploader = re.sub('_', ' ', i.group('uploader'))
            uploader = undoHTMLEntities(text=uploader)
            uploader = urllib.unquote(uploader)
            images.append([filename, url, uploader])           
                    
    print '    Found %d images' % (len(images))
    images.sort()
    return images

def undoHTMLEntities(text=''):
    """ Undo some HTML codes """
    text = re.sub('&amp;lt;', '&lt;', text) # i guess only &lt; &gt; &amp; " ' need conversion http://www.w3schools.com/html/html_entities.asp
    text = re.sub('&amp;gt;', '&gt;', text)
    text = re.sub('&amp;amp;', '&amp;', text)
    text = re.sub('&amp;quot;', '"', text)
    text = re.sub('&amp;#039;', '\'', text)
    return text

def generateImageDump(config={}, other={}, images=[], start=''):
    """ Save files and descriptions using a file list """
    #fix use subdirectories md5
    print 'Retrieving images from "%s"' % (start and start or 'start')
    imagepath = '%s/images' % (config['path'])
    if not os.path.isdir(imagepath):
        print 'Creating "%s" directory' % (imagepath)
        os.makedirs(imagepath)
    
    c = 0
    lock = True
    if not start:
        lock = False
    for filename, url, uploader in images:
        if filename == start: #start downloading from start (included)
            lock = False
        if lock:
            continue
        delay(config=config)
        
        #saving file
        #truncate filename if length &gt; 100 (100 + 32 (md5) = 132 &lt; 143 (crash limit). Later .desc is added to filename, so better 100 as max)
        filename2 = filename
        if len(filename2) &gt; other['filenamelimit']:
            # split last . (extension) and then merge
            filename2 = truncateFilename(other=other, filename=filename2)
            print 'Filename is too long, truncating. Now it is:', filename2
        urllib.urlretrieve(url=url, filename='%s/%s' % (imagepath, filename2), data=urllib.urlencode({})) #fix, image request fails on wikipedia (POST neither works?)
        
        #saving description if any
        xmlfiledesc = getXMLFileDesc(config=config, title='Image:%s' % (filename)) # use Image: for backwards compatibility
        f = open('%s/%s.desc' % (imagepath, filename2), 'w')
        if not re.search(r'&lt;/mediawiki&gt;', xmlfiledesc): #&lt;text xml:space="preserve" bytes="36"&gt;Banner featuring SG1, SGA, SGU teams&lt;/text&gt;
            #failure when retrieving desc? then save it as empty .desc
            xmlfiledesc = ''
        f.write(xmlfiledesc)
        f.close()
        c += 1
        if c % 10 == 0:
            print '    Downloaded %d images' % (c)
    print 'Downloaded %d images' % (c)
    
def saveLogs(config={}):
    """ Save Special:Log """
    #get all logs from Special:Log
    """parse
    &lt;select name='type'&gt;
    &lt;option value="block"&gt;Bloqueos de usuarios&lt;/option&gt;
    &lt;option value="rights"&gt;Cambios de perfil de usuario&lt;/option&gt;
    &lt;option value="protect" selected="selected"&gt;Protecciones de pÃ¡ginas&lt;/option&gt;
    &lt;option value="delete"&gt;Registro de borrados&lt;/option&gt;
    &lt;option value="newusers"&gt;Registro de creaciÃ³n de usuarios&lt;/option&gt;
    &lt;option value="merge"&gt;Registro de fusiones&lt;/option&gt;
    &lt;option value="import"&gt;Registro de importaciones&lt;/option&gt;
    &lt;option value="patrol"&gt;Registro de revisiones&lt;/option&gt;
    &lt;option value="move"&gt;Registro de traslados&lt;/option&gt;
    &lt;option value="upload"&gt;Subidas de archivos&lt;/option&gt;
    &lt;option value=""&gt;Todos los registros&lt;/option&gt;
    &lt;/select&gt;
    """
    delay(config=config)

def domain2prefix(config={}):
    """ Convert domain name to a valid prefix filename """
    domain = ''
    if config['api']:
        domain = config['api']
    elif config['index']:
        domain = config['index']
    domain = domain.lower()
    domain = re.sub(r'(https?://|www\.|/index\.php|/api\.php)', '', domain)
    domain = re.sub(r'/', '_', domain)
    domain = re.sub(r'\.', '', domain)
    domain = re.sub(r'[^A-Za-z0-9]', '_', domain)
    return domain

def loadConfig(config={}, configfilename=''):
    """ Load config file """
    try:
        f = open('%s/%s' % (config['path'], configfilename), 'r')
    except:
        print 'There is no config file. we can\'t resume. Start a new dump.'
        sys.exit()
    config = cPickle.load(f)
    f.close()
    return config

def saveConfig(config={}, configfilename=''):
    """ Save config file """
    f = open('%s/%s' % (config['path'], configfilename), 'w')
    cPickle.dump(config, f)
    f.close()
    
def welcome():
    """ Opening message """
    print "#"*73
    print """# Welcome to DumpGenerator 0.1 by WikiTeam (GPL v3)                     #
# More info at: http://code.google.com/p/wikiteam/                      #"""
    print "#"*73
    print ''
    print "#"*73
    print """# Copyright (C) 2011-2013 WikiTeam                                      #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see &lt;http://www.gnu.org/licenses/&gt;. #"""
    print "#"*73
    print ''

def bye():
    """ Closing message """
    print "---&gt; Congratulations! Your dump is complete &lt;---"
    print "If you found any bug, report a new issue here (Google account required): http://code.google.com/p/wikiteam/issues/list"
    print "If this is a public wiki, please, consider publishing this dump. Do it yourself as explained in http://code.google.com/p/wikiteam/wiki/NewTutorial#Publishing_the_dump or contact us at http://code.google.com/p/wikiteam"
    print "Good luck! Bye!"

def usage():
    """  """
    print """Error. You forget mandatory parameters:
    --api or --index: URL to api.php or to index.php, one of them. Examples: --api=http://archiveteam.org/api.php or --index=http://archiveteam.org/index.php
    
And one of these at least:
    --xml: It generates a XML dump. It retrieves full history of all pages (if you want only the current version use --xml --curonly)
           If you want filter by namespace, use the parameter --namespaces=0,1,2,3...
    --images: It generates an image dump

You can resume previous incomplete dumps:
    --resume: It resumes previous incomplete dump. When using --resume, --path is mandatory (path to directory where incomplete dump is).

You can exclude namespaces:
    --exnamespaces: Write the number of the namespaces you want to exclude, split by commas.

You can use authenticaton cookies from a Mozilla cookies.txt file:
    --cookies: Path to a cookies.txt file. Example: --cookies=$HOME/.netscape/cookies.txt

You can be nice with servers using a delay:
    --delay: It adds a delay (in seconds, adding 5 seconds between requests: --delay:5)

Write --help for help."""

def getParameters(params=[]):
    if not params:
        params = sys.argv[1:]
    config = {
        'curonly': False,
        'date': datetime.datetime.now().strftime('%Y%m%d'),
        'api': '',
        'index': '',
        'images': False,
        'logs': False,
        'xml': False,
        'namespaces': ['all'],
        'exnamespaces': [],
        'path': '',
        'cookies': '',
        'delay': 0,
    }
    other = {
        'resume': False,
        'filenamelimit': 100, #do not change
        'force': False,
    }
    #console params
    try:
        opts, args = getopt.getopt(params, "", ["h", "help", "path=", "api=", "index=", "images", "logs", "xml", "curonly", "resume", "cookies=", "delay=", "namespaces=", "exnamespaces=", "force", ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("--path"):
            config["path"] = a
            while len(config["path"])&gt;0:
                if config["path"][-1] == '/': #darÃ¡ problemas con rutas windows?
                    config["path"] = config["path"][:-1]
                else:
                    break
        elif o in ("--api"):
            if not a.startswith('http://') and not a.startswith('https://'):
                print 'api.php must start with http:// or https://'
                sys.exit()
            config['api'] = a
        elif o in ("--index"):
            if not a.startswith('http://') and not a.startswith('https://'):
                print 'index.php must start with http:// or https://'
                sys.exit()
            config["index"] = a
        elif o in ("--images"):
            config["images"] = True
        elif o in ("--logs"):
            config["logs"] = True
        elif o in ("--xml"):
            config["xml"] = True
        elif o in ("--curonly"):
            if not config["xml"]:
                print "If you select --curonly, you must use --xml too"
                sys.exit()
            config["curonly"] = True
        elif o in ("--resume"):
            other["resume"] = True
        elif o in ("--cookies"):
            config["cookies"] = a
        elif o in ("--delay"):
            config["delay"] = int(a)
        elif o in ("--namespaces"):
            if re.search(r'[^\d, \-]', a) and a.lower() != 'all': #fix, why - ?  and... --namespaces= all with a space works?
                print "Invalid namespaces values.\nValid format is integer(s) splitted by commas"
                sys.exit()
            a = re.sub(' ', '', a)
            if a.lower() == 'all':
                config["namespaces"] = ['all']
            else:
                config["namespaces"] = [int(i) for i in a.split(',')]
        elif o in ("--exnamespaces"):
            if re.search(r'[^\d, \-]', a):
                print "Invalid exnamespaces values.\nValid format is integer(s) splitted by commas"
                sys.exit()
            a = re.sub(' ', '', a)
            if a.lower() == 'all':
                print 'You have excluded all namespaces. Error.'
                sys.exit()
            else:
                config["exnamespaces"] = [int(i) for i in a.split(',')]
        elif o in ("--force"):
            other["force"] = True
        else:
            assert False, "unhandled option"

    #missing mandatory params
    #(config['index'] and not re.search('/index\.php', config['index'])) or \ # in EditThis there is no index.php, it is empty editthis.info/mywiki/?title=...
    if (not config['api'] and not config['index']) or \
       (config['api'] and not re.search('/api\.php', config['api'])) or \
       not (config["xml"] or config["images"] or config["logs"]) or \
       (other['resume'] and not config['path']):
        usage()
        sys.exit()
    
    #user chosen --api, --index it is neccesary for special:export, we generate it
    if config['api'] and not config['index']:
        config['index'] = config['api'].split('api.php')[0] + 'index.php'
        #print 'You didn\'t provide a path for index.php, trying to wonder one:', config['index']
    
    if config['cookies']:
        cj = cookielib.MozillaCookieJar()
        cj.load(config['cookies'])
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        urllib2.install_opener(opener)
        print 'Using cookies from %s' % config['cookies']

    if config['api']:
        #check api.php
        if checkAPI(config['api']):
            print 'api.php is OK'
        else:
            print 'Error in api.php, please, provide a correct path to api.php'
            sys.exit()
    
    if config['index']:
        #check index.php
        if checkIndexphp(config['index']):
            print 'index.php is OK'
        else:
            print 'Error in index.php, please, provide a correct path to index.php'
            sys.exit()
    
    #calculating path, if not defined by user with --path=
    if not config['path']:
        config['path'] = './%s-%s-wikidump' % (domain2prefix(config=config), config['date'])
    
    return config, other

def checkAPI(api):
    """ Checking API availability """
    req = urllib2.Request(url=api, headers={'User-Agent': getUserAgent()})
    f = urllib2.urlopen(req)
    raw = f.read()
    f.close()
    print 'Checking api.php...', api
    if re.search(r'action=query', raw):
        return True
    return False

def checkIndexphp(indexphp):
    """ Checking index.php availability """
    req = urllib2.Request(url=indexphp, data=urllib.urlencode({'title': 'Special:Version', }), headers={'User-Agent': getUserAgent()})
    f = urllib2.urlopen(req)
    raw = f.read()
    f.close()
    print 'Checking index.php...', indexphp
    if re.search(r'(This wiki is powered by|&lt;h2 id="mw-version-license"&gt;|meta name="generator" content="MediaWiki)', raw):
        return True
    return False

def removeIP(raw=''):
    """ Remove IP from HTML comments &lt;!-- --&gt; """
    raw = re.sub(r'\d+\.\d+\.\d+\.\d+', '0.0.0.0', raw)
    #http://www.juniper.net/techpubs/software/erx/erx50x/swconfig-routing-vol1/html/ipv6-config5.html
    #weird cases as :: are not included
    raw = re.sub(r'(?i)[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}:[\da-f]{0,4}', '0:0:0:0:0:0:0:0', raw)
    return raw

def checkXMLIntegrity(config={}):
    """ Check XML dump integrity, to detect broken XML chunks """
    return 
    
    #TODO fix, instead relaunch, ask first to user, also it fails on windows without grep : (
    print "Verifying dump..."
    os.chdir(config['path'])
    checktitles = os.system('grep "&lt;title&gt;" *.xml -c &gt; /dev/null')
    checkpageopen = os.system('grep "&lt;page&gt;" *.xml -c &gt; /dev/null')
    checkpageclose = os.system('grep "&lt;/page&gt;" *.xml -c &gt; /dev/null')
    checkrevisionopen = os.system('grep "&lt;revision&gt;" *.xml -c &gt; /dev/null')
    checkrevisionclose = os.system('grep "&lt;/revision&gt;" *.xml -c &gt; /dev/null')
    os.chdir('..')
    if (checktitles == checkpageopen and checktitles == checkpageclose and checkpageopen == checkpageclose):
        xmlisgood = True
    else:
        xmlisgood = False
        print "XML dump is corrupted, regenerating a new dump"
        generateXMLDump(config=config, titles=titles)

def createNewDump(config={}, other={}):
    titles = []
    images = []
    print 'Trying generating a new dump into a new directory...'
    if config['xml']:
        titles += getPageTitles(config=config)
        saveTitles(config=config, titles=titles)
        generateXMLDump(config=config, titles=titles)
        checkXMLIntegrity(config=config)
    if config['images']:
        if config['api']:
            images += getImageFilenamesURLAPI(config=config)
        else:
            images += getImageFilenamesURL(config=config)
        saveImageFilenamesURL(config=config, images=images)
        generateImageDump(config=config, other=other, images=images)
    if config['logs']:
        saveLogs(config=config)

def resumePreviousDump(config={}, other={}):
    titles = []
    images = []
    print 'Resuming previous dump process...'
    if config['xml']:
        #load titles
        lasttitle = ''
        try:
            f = open('%s/%s-%s-titles.txt' % (config['path'], domain2prefix(config=config), config['date']), 'r')
            raw = f.read()
            titles = raw.split('\n')
            lasttitle = titles[-1]
            if not lasttitle: #empty line at EOF ?
                lasttitle = titles[-2]
            f.close()
        except:
            pass #probably file doesnot exists
        if lasttitle == '--END--':
            #titles list is complete
            print 'Title list was completed in the previous session'
        else:
            print 'Title list is incomplete. Reloading...'
            #do not resume, reload, to avoid inconsistences, deleted pages or so
            titles = getPageTitles(config=config)
            saveTitles(config=config, titles=titles)
        #checking xml dump
        xmliscomplete = False
        lastxmltitle = ''
        try:
            f = open('%s/%s-%s-%s.xml' % (config['path'], domain2prefix(config=config), config['date'], config['curonly'] and 'current' or 'history'), 'r')
            for l in f:
                if re.findall('&lt;/mediawiki&gt;', l):
                    #xml dump is complete
                    xmliscomplete = True
                    break
                xmltitles = re.findall(r'&lt;title&gt;([^&lt;]+)&lt;/title&gt;', l) #weird if found more than 1, but maybe
                if xmltitles:
                    lastxmltitle = undoHTMLEntities(text=xmltitles[-1])
            f.close()
        except:
            pass #probably file doesnot exists
        #removing --END-- before getXMLs
        while titles and titles[-1] in ['', '--END--']:
            titles = titles[:-1]
        if xmliscomplete:
            print 'XML dump was completed in the previous session'
        elif lastxmltitle:
            #resuming...
            print 'Resuming XML dump from "%s"' % (lastxmltitle)
            generateXMLDump(config=config, titles=titles, start=lastxmltitle)
        else:
            #corrupt? only has XML header?
            print 'XML is corrupt? Regenerating...'
            generateXMLDump(config=config, titles=titles)
    
    if config['images']:
        #load images
        lastimage = ''
        try:
            f = open('%s/%s-%s-images.txt' % (config['path'], domain2prefix(config=config), config['date']), 'r')
            raw = f.read()
            lines = raw.split('\n')
            for l in lines:
                if re.search(r'\t', l):
                    images.append(l.split('\t'))
            lastimage = lines[-1]
            f.close()
        except:
            pass #probably file doesnot exists
        if lastimage == '--END--':
            print 'Image list was completed in the previous session'
        else:
            print 'Image list is incomplete. Reloading...'
            #do not resume, reload, to avoid inconsistences, deleted images or so
            if config['api']:
                images=getImageFilenamesURLAPI(config=config)
            else:
                images = getImageFilenamesURL(config=config)
            saveImageFilenamesURL(config=config, images=images)
        #checking images directory
        listdir = []
        try:
            listdir = os.listdir('%s/images' % (config['path']))
        except:
            pass #probably directory does not exist
        listdir.sort()
        complete = True
        lastfilename = ''
        lastfilename2 = ''
        c = 0
        for filename, url, uploader in images:
            lastfilename2 = lastfilename
            lastfilename = filename #return always the complete filename, not the truncated
            filename2 = filename
            if len(filename2) &gt; other['filenamelimit']:
                filename2 = truncateFilename(other=other, filename=filename2)
            if filename2 not in listdir:
                complete = False
                break
            c +=1
        print '%d images were found in the directory from a previous session' % (c)
        if complete:
            #image dump is complete
            print 'Image dump was completed in the previous session'
        else:
            generateImageDump(config=config, other=other, images=images, start=lastfilename2) # we resume from previous image, which may be corrupted (or missing .desc)  by the previous session ctrl-c or abort
    
    if config['logs']:
        #fix
        pass

def saveSpecialVersion(config={}):
    #save Special:Version as .html, to preserve extensions details
    if os.path.exists('%s/Special:Version.html' % (config['path'])):
        print 'Special:Version.html exists, do not overwrite'
    else:
        print 'Downloading Special:Version with extensions and other related info'
        req = urllib2.Request(url=config['index'], data=urllib.urlencode({'title': 'Special:Version', }), headers={'User-Agent': getUserAgent()})
        f = urllib2.urlopen(req)
        raw = f.read()
        f.close()
        raw = removeIP(raw=raw)
        f = open('%s/Special:Version.html' % (config['path']), 'w')
        f.write(raw)
        f.close()

def saveIndexPHP(config={}):
    #save index.php as .html, to preserve license details available at the botom of the page
    if os.path.exists('%s/index.html' % (config['path'])):
        print 'index.html exists, do not overwrite'
    else:
        print 'Downloading index.php (Main Page) as index.html'
        req = urllib2.Request(url=config['index'], data=urllib.urlencode({}), headers={'User-Agent': getUserAgent()})
        f = urllib2.urlopen(req)
        raw = f.read()
        f.close()
        raw = removeIP(raw=raw)
        f = open('%s/index.html' % (config['path']), 'w')
        f.write(raw)
        f.close()

def avoidWikimediaProjects(config={}):
    """ Skip Wikimedia projects and redirect to the dumps website """
    #notice about wikipedia dumps
    if re.findall(r'(?i)(wikipedia|wikisource|wiktionary|wikibooks|wikiversity|wikimedia|wikispecies|wikiquote|wikinews|wikidata|wikivoyage)\.org', config['api']+config['index']):
        print 'PLEASE, DO NOT USE THIS SCRIPT TO DOWNLOAD WIKIMEDIA PROJECTS!'
        print 'Download the dumps from http://dumps.wikimedia.org'
        if not other['force']:
            print 'Thanks!'
            sys.exit()

def main(params=[]):
    """ Main function """
    welcome()
    configfilename = 'config.txt'
    config, other = getParameters(params=params)
    avoidWikimediaProjects(config=config)
    print 'Analysing %s' % (config['api'] and config['api'] or config['index'])
    
    #creating path or resuming if desired
    c = 2
    originalpath = config['path'] # to avoid concat blabla-2, blabla-2-3, and so on...
    while not other['resume'] and os.path.isdir(config['path']): #do not enter if resume is requested from begining
        print '\nWarning!: "%s" path exists' % (config['path'])
        reply = ''
        while reply.lower() not in ['yes', 'y', 'no', 'n']:
            reply = raw_input('There is a dump in "%s", probably incomplete.\nIf you choose resume, to avoid conflicts, the parameters you have chosen in the current session will be ignored\nand the parameters available in "%s/%s" will be loaded.\nDo you want to resume ([yes, y], [no, n])? ' % (config['path'], config['path'], configfilename))
        if reply.lower() in ['yes', 'y']:
            if not os.path.isfile('%s/%s' % (config['path'], configfilename)):
                print 'No config file found. I can\'t resume. Aborting.'
                sys.exit()
            print 'You have selected: YES'
            other['resume'] = True
            break
        elif reply.lower() in ['no', 'n']:
            print 'You have selected: NO'
            other['resume'] = False
        config['path'] = '%s-%d' % (originalpath, c)
        print 'Trying to use path "%s"...' % (config['path'])
        c += 1

    if other['resume']:
        print 'Loading config file...'
        config = loadConfig(config=config, configfilename=configfilename)
    else:
        os.mkdir(config['path'])
        saveConfig(config=config, configfilename=configfilename)
    
    if other['resume']:
        resumePreviousDump(config=config, other=other)
    else:
        createNewDump(config=config, other=other)

    saveIndexPHP(config=config)    
    saveSpecialVersion(config=config)
    bye()

if __name__ == "__main__":
    main()
