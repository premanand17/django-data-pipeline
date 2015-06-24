''' Used to fetch publication details from NCBI and generate a JSON. '''

import json
import requests
import xml.etree.ElementTree as ET
import logging
import re
from .exceptions import PublicationDownloadError
import time

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Pubs():

    DUPLICATE_PMIDS = ['22543779']

    @classmethod
    def fetch_details(cls, pmids, filename, disease_code=None, source='auto'):
        ''' Given a list of PMIDs fetch their details from eutils.
        http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi??db=pubmed&retmode=xml&id=<PMIDS>
        Produces a JSON file containing the publications mapping and documents.
        {
          "mapping": {"properties": {...}},
          "docs": [...]
        }
        '''

        # remove known duplicate PMIDs
        pmids = [pmid for pmid in pmids if pmid not in Pubs.DUPLICATE_PMIDS]
        chunk_size = 450
        count = 0
        mapping = {
            "_id": {"type": "integer"},
            "PMID": {"type": "integer"},
            "tags": {"type": "object", "index": "not_analyzed"},
            "journal": {"type": "string"},
            "title": {"type": "string"},
            "date": {"type": "date"},
            "authors": {"type": "object"},
            "abstract": {"type": "string"}
                   }
        mapping_keys = mapping.keys()
        start = time.time()

#         pmids = [25905407, 25905392, 23369186, 24947582, 1476675, 18225448, 10250814, 25743292]
        with open(filename, mode='w', encoding='utf-8') as f:
            f.write('{"mapping": ')
            f.write(json.dumps({"properties": mapping}))
            f.write(',\n"docs":[\n')
            for i in range(0, len(pmids), chunk_size):
                chunk = pmids[i:i+chunk_size]
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' \
                      "?db=pubmed&retmode=xml&id=%s" % \
                      ",".join([str(item) for item in chunk])

                r = requests.get(url, timeout=25)
                tree = ET.fromstring(r.content)
                pubmeds = tree.findall("PubmedArticle") + tree.findall("PubmedBookArticle")
                for pubmed in pubmeds:
                    pub = pubmed.find('MedlineCitation')
                    if pub is None:
                        pub = pubmed.find('BookDocument')
                    if count > 0:
                        f.write(',\n')

                    pub_obj = Pubs._parse_pubmed_record(pub)
                    if disease_code is not None:
                        pub_obj['tags'] = {}
                        pub_obj['tags']['disease'] = [disease_code]
                    if source is not None:
                        if 'tags' not in pub_obj:
                            pub_obj['tags'] = {}
                        pub_obj['tags']['source'] = source

                    keys_not_found = [k for k in mapping_keys if k not in pub_obj]
                    if len(keys_not_found) > 0:
                        logger.warn("PMID: "+pub_obj['PMID']+' not found: '+str(keys_not_found))
                    f.write(json.dumps(pub_obj))
                    count += 1

                time_taken = time.time() - start
                eta = (time_taken / (i+chunk_size)) * (len(pmids) - i+chunk_size)
                logger.debug('Retrieved '+(str(i+chunk_size))+' PMID records of '+str(len(pmids)) +
                             ' :: ETA/s: '+str(int(eta)))

            f.write('\n]}')
        logger.debug("No. publications downloaded "+str(count))
        if count != len(pmids):
            msg = "No. publications does not match the number of requested PMIDs ="+str(len(pmids))
            logger.error(msg)
            raise PublicationDownloadError(msg)

    @classmethod
    def _parse_pubmed_record(cls, pub):
        pmid = pub.find('PMID').text
        pub_obj = {'PMID': pmid, '_id': pmid}
        article = pub.find('Article')
        if article is not None:
            pub_obj['title'] = article.find('ArticleTitle').text
            authors = article.find('AuthorList')
            Pubs.get_authors(pub_obj, authors, pmid)
            Pubs.get_abstract(pub_obj, article)

            journal = article.find('Journal')
            try:
                pub_obj['journal'] = journal.find('ISOAbbreviation').text
            except AttributeError:
                pub_obj['journal'] = journal.find('Title').text

            pub_date = journal.find('JournalIssue').find('PubDate')
            article_date = article.find('ArticleDate')
            if article_date is not None:
                pub_date = article_date
            Pubs.get_date(pub_obj, pub_date)
        elif pub.find('Book') is not None:
            pub_obj['title'] = pub.find('ArticleTitle').text
            authors = pub.find('AuthorList')
            Pubs.get_authors(pub_obj, authors, pmid)
            Pubs.get_abstract(pub_obj, pub)
            pub_date = pub.find('ContributionDate')
            if pub_date is None:
                pub_date = pub.find('Book').find('PubDate')

            Pubs.get_date(pub_obj, pub_date)

        return pub_obj

    @classmethod
    def get_abstract(cls, pub_obj, article):
        ''' Add the abastract to the publication object. '''
        try:
            texts = article.find('Abstract').findall('AbstractText')
            abstract = ''
            for t in texts:
                label = ''
                if hasattr(t, 'Label'):
                    label = t.get('Label') + ': '
                try:
                    abstract = abstract + label + t.text
                except TypeError:
                    continue
            pub_obj['abstract'] = abstract
        except AttributeError:
            return

    @classmethod
    def get_authors(cls, pub_obj, authors, pmid):
        ''' Add the author list to the publication object. '''
        authors_arr = []
        try:
            for author in authors:
                try:
                    author_obj = {'LastName': author.find('LastName').text,
                                  'ForeName': author.find('ForeName').text}
                except AttributeError:
                    if author.find('LastName') is None:
                        continue
                    author_obj = {'LastName': author.find('LastName').text}
                if author.find('Initials') is not None:
                    author_obj.update({'Initials': author.find('Initials').text})
                authors_arr.append(author_obj)
            pub_obj['authors'] = authors_arr
        except TypeError:
            logger.warn('No authors found for PMID:'+pmid)

    # month mappings
    MONTHS = {'jan': '01',
              'feb': '02',
              'mar': '03',
              'apr': '04',
              'may': '05',
              'jun': '06',
              'jul': '07',
              'aug': '08',
              'sep': '09',
              'oct': '10',
              'nov': '11',
              'dec': '12',
              'sum': '06', 'summer': '06',
              'win': '12', 'winter': '12',
              'spr': '03', 'spring': '03',
              'fal': '09', 'fall': '09',
              'aut': '09', 'autumn': '09'}

    @classmethod
    def get_date(cls, pub_obj, pub_date):
        ''' Get the date and save to pub_obj. '''
        if pub_date.find('Month') is not None:
            month = pub_date.find('Month').text
            if month.lower() in Pubs.MONTHS:
                month = Pubs.MONTHS[month.lower()]
            date = pub_date.find('Year').text + '-' + month
            if pub_date.find('Day') is not None:
                day = '%02d' % int(pub_date.find('Day').text)
                date = date + '-' + day
            else:
                date = date + '-01'
            pub_obj['date'] = date
        elif pub_date.find('Year') is not None:
            pub_obj['date'] = pub_date.find('Year').text + '-01-01'
        elif pub_date.find('MedlineDate') is not None:
            date = pub_date.find('MedlineDate').text

            # 1999 May-Jun and 1992 Summer-Fall
            p = re.compile('(\d{4})\s(\w{3,6})\s*-\w+')
            m = p.match(date)
            if m:
                date = m.group(1) + '-' + Pubs.MONTHS[m.group(2).lower()] + '-01'
            else:
                # 2010 May 26-Jun 1
                p = re.compile('(\d{4})\s(\w{3}) (\d{1,2})-')
                m = p.match(date)
                if m:
                    day = '%02d' % int(m.group(3))
                    date = m.group(1) + '-' + Pubs.MONTHS[m.group(2).lower()] + '-' + day
                else:
                    # 1978-1979 and 1981 1st Quart
                    p = re.compile('^(\d{4})\s*(-|1st|2nd|2d|3rd|4th)')
                    m = p.match(date)
                    if m:
                        if m.group(2) == '2nd' or m.group(2) == '2d':
                            date = m.group(1) + '-04-01'
                        elif m.group(2) == '3rd':
                            date = m.group(1) + '-07-01'
                        elif m.group(2) == '4th':
                            date = m.group(1) + '-10-01'
                        else:
                            date = m.group(1) + '-01-01'
                    else:
                        # 2000Jun 8-21
                        p = re.compile('^(\d{4})\s*(\w{3,6})-*')
                        m = p.match(date)
                        if m:
                            date = m.group(1) + '-' + Pubs.MONTHS[m.group(2).lower()] + '-01'

            pub_obj['date'] = date
        else:
            logger.warn("Date not found for PMID:"+pub_obj["PMID"])
