from termcolor import colored


class Error(Exception):
    def __init__(self, keyword, message, source=None):
        print(colored(f'ERROR: Wrong \"{keyword}\" argument.', 'red'))
        self.message = message


class InputSDFileError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class FileScreenOutputError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class ClassifierNameError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class InputSMARTSError(Exception):
    def __init__(self, message, source=None):
        super().__init__(message)
        self.source = source


