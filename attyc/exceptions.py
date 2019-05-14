from termcolor import colored


class Error(Exception):
    def __init__(self, keyword, message):
        if keyword == 'smarts':
            print(colored(f'ERROR: Wrong SMARTS file.', 'red'))
        else:
            print(colored(f'ERROR: Wrong \"{keyword}\" argument.', 'red'))
        self.message = message


class InputFileError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class FileScreenOutputError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class ClassifierNameError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class ExternalTypesInputFileError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)


class ClassifierClassError(Error):
    def __init__(self, keyword, message):
        super().__init__(keyword, message)

