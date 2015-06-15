import json
import requests
import xml.etree.ElementTree as ET
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Pubs():

    @classmethod
    def fetch_details(cls, pmids, filename):
        ''' Given a list of PMIDs fetch their details from eutils.
        http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi??db=pubmed&retmode=xml&id=<PMIDS>
        '''
        chunkSize = 450
        count = 0
        mapping = {
            "PMID": {"type": "integer"},
            "journal": {"type": "string"},
            "title": {"type": "string"},
            "date": {"type": "string", "index": "not_analyzed"},
            "authors": {"type": "object"},
            "abstract": {"type": "string", "index": "no"}
                   }

#         pmids = [25905392, 23369186, 24947582]
        with open(filename, mode='w', encoding='utf-8') as f:
            f.write('{"mapping": ')
            f.write(json.dumps({"properties": mapping}))
            f.write(',\n"docs":[\n')
            for i in range(0, len(pmids), chunkSize):
                chunk = pmids[i:i+chunkSize]
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' \
                      "?db=pubmed&retmode=xml&id=%s" % \
                      ",".join([str(item) for item in chunk])

                r = requests.get(url, timeout=25)
                tree = ET.fromstring(r.content)
                pubmeds = tree.findall("PubmedArticle") + tree.findall("PubmedBookArticle")
                for pubmed in pubmeds:
                    pub = pubmed.find('MedlineCitation')
                    pub_data = pubmed.find('PubmedData')
                    if pub is None:
                        pub = pubmed.find('BookDocument')
                    if count > 0:
                        f.write(',\n')

                    pub_obj = Pubs._parse_pubmed_record(pub)

                    if pub_data is not None:
                        Pubs.get_date_from_pubmeddate(pub_obj, pub_data)
                    f.write(json.dumps(pub_obj))
                    count += 1

            f.write('\n]}')
        logger.debug("No. publications downloaded "+str(count))

    @classmethod
    def _parse_pubmed_record(cls, pub):
        pmid = pub.find('PMID')
        logger.debug(pmid.text)
        pub_obj = {'PMID': pmid.text}
        article = pub.find('Article')
        if article is not None:
            pub_obj['title'] = article.find('ArticleTitle').text
            authors = article.find('AuthorList')
            Pubs.get_authors(pub_obj, authors, pmid)
            Pubs.get_abstract(pub_obj, article, pmid)
            journal = article.find('Journal')
            try:
                pub_obj['journal'] = journal.find('ISOAbbreviation').text
            except AttributeError:
                pub_obj['journal'] = journal.find('Title').text

            pub_date = journal.find('JournalIssue').find('PubDate')
            pub_obj['date'] = Pubs.get_date(pub_date)
        elif pub.find('Book') is not None:
            book = pub.find('Book')
            pub_obj['title'] = pub.find('ArticleTitle').text
            authors = book.find('AuthorList')
            Pubs.get_authors(pub_obj, authors, pmid)
            Pubs.get_abstract(pub_obj, pub, pmid)

        return pub_obj

    @classmethod
    def get_abstract(cls, pub_obj, article, pmid):
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
            logger.warn('No abstract found for PMID:'+pmid.text)

    @classmethod
    def get_authors(cls, pub_obj, authors, pmid):
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
            logger.warn('No authors found for PMID:'+pmid.text)

    @classmethod
    def get_date_from_pubmeddate(cls, pub_obj, pub_date):

        history = pub_date.find('History')
        if history is not None:
            pdates = history.findall('PubMedPubDate')
            for pdate in pdates:
                if hasattr(pdate, 'PubStatus'):
                    if getattr(pdate, 'PubStatus') == 'accepted':
                        pub_obj['date'] = pdate.find('Year') + pdate.find('Month') + pdate.find('Day')
                        return

    @classmethod
    def get_date(cls, pub_date):
        if pub_date.find('Month') is not None:
            date = (pub_date.find('Month').text + ' ' +
                    pub_date.find('Year').text)
            if pub_date.find('Day') is not None:
                date = pub_date.find('Day').text + ' ' + date
            return date
        elif pub_date.find('Year') is not None:
            return pub_date.find('Year').text
        elif pub_date.find('MedlineDate') is not None:
            return pub_date.find('MedlineDate').text
