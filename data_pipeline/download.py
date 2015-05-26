import requests
import ftputil
from urllib.parse import urlparse
import ftplib
import time


class Download:

    def download(self, url, download_dir):
        if url.startswith("ftp://"):
            self.ftp_download(url, download_dir)
        else:
            self.download_file(url, download_dir)

    def ftp_download(self, url, download_dir):
        local_filename = url.split('/')[-1]
        o = urlparse(url)
        ftp_host = ftputil.FTPHost(o.netloc, 'anonymous', '',
                                   session_factory=ftplib.FTP)
        size = ftp_host.path.getsize(o.path)
        ftp_host.download(o.path, local_filename, callback=Monitor(size))

    def download_file(self, url, download_dir):
        local_filename = url.split('/')[-1]
        r = requests.get(url, stream=True)
        monitor = Monitor(r.headers.get('content-length'))
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
                    monitor(chunk)
        return local_filename


class Monitor:

    def __init__(self, size):
        self.size = int(size)
        self.size_progress = 0
        self.previous = 0
        self.start = time.time()

    def __call__(self, chunk):
        self.size_progress += len(chunk)
        percent_progress = int(self.size_progress/self.size * 100)
        if percent_progress != self.previous and percent_progress % 10 == 0:
            time_taken = time.time() - self.start
            eta = (time_taken / self.size_progress) * (self.size - self.size_progress)
            print("\r"+str(self.previous) + "% (eta: "+str(int(eta))+"s)", end="")
            self.previous = percent_progress
