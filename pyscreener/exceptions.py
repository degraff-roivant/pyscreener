class BadPDBFileError(Exception):
    pass


class InvalidEnvironmentError(Exception):
    pass


class MissingExecutableError(InvalidEnvironmentError):
    pass


class MissingEnvironmentVariableError(InvalidEnvironmentError):
    pass


class MisconfiguredDirectoryError(InvalidEnvironmentError):
    pass


class NotSimulatedError(Exception):
    pass


class InvalidResultError(Exception):
    pass


class UnsupportedSoftwareError(Exception):
    pass


class ReceptorPreparationError(Exception):
    pass
