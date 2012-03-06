class LibadjointException(Exception):
  pass

class LibadjointErrorException(LibadjointException):
  pass

class LibadjointErrorInvalidInputs(LibadjointErrorException):
  pass

class LibadjointErrorHashFailed(LibadjointErrorException):
  pass

class LibadjointErrorNeedCallback(LibadjointErrorException):
  pass

class LibadjointErrorNeedValue(LibadjointErrorException):
  pass

class LibadjointErrorNotImplemented(LibadjointErrorException):
  pass

class LibadjointErrorDictFailed(LibadjointErrorException):
  pass

class LibadjointErrorToleranceExceeded(LibadjointErrorException):
  pass

class LibadjointErrorMallocFailed(LibadjointErrorException):
  pass

class LibadjointErrorRevolveError(LibadjointErrorException):
  pass

class LibadjointErrorSlepcError(LibadjointErrorException):
  pass

class LibadjointWarnException(LibadjointException):
  pass

class LibadjointWarnAlreadyRecorded(LibadjointWarnException):
  pass

class LibadjointWarnComparisonFailed(LibadjointWarnException):
  pass

class LibadjointWarnUninitialisedValue(LibadjointWarnException):
  pass

class LibadjointWarnNotImplemented(LibadjointWarnException):
  pass

ierr_to_exception = {
    1: LibadjointErrorInvalidInputs,
    2: LibadjointErrorHashFailed,
    3: LibadjointErrorNeedCallback,
    4: LibadjointErrorNeedValue,
    5: LibadjointErrorNotImplemented,
    6: LibadjointErrorDictFailed,
    7: LibadjointErrorToleranceExceeded,
    8: LibadjointErrorMallocFailed,
    9: LibadjointErrorRevolveError,
    10: LibadjointErrorSlepcError
   -1: LibadjointWarnAlreadyRecorded,
   -2: LibadjointWarnComparisonFailed,
   -3: LibadjointWarnUninitialisedValue,
   -4: LibadjointWarnNotImplemented
   }

def get_exception(ierr):
  try:
    return ierr_to_exception[ierr]
  except IndexError:
    return None

