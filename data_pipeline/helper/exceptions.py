

class PublicationDownloadError(Exception):
    ''' Publication download error  '''
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
