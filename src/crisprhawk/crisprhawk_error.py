""" """


class CrisprHawkError(Exception):
    def __init__(self, value: str):
        # initialize exception object when raised
        self._value = value  # error message or error related info

    def __str__(self):
        return repr(self._value)  # string representation for the exception


class CrisprHawkBedError(CrisprHawkError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkFastaError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkPamError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkBitsetError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkIupacTableError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkGuidesReportError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkVCFError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkEnrichmentError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkHaplotypeError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkGuideError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkVariantMapError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkScoreError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkAzimuthScoreError(CrisprHawkScoreError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkRs3ScoreError(CrisprHawkScoreError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkCfdScoreError(CrisprHawkScoreError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception

class CrisprHawkElevationScoreError(CrisprHawkScoreError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkDeepCpf1ScoreError(CrisprHawkScoreError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkOffTargetsError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkAnnotationError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkConverterError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception


class CrisprHawkPrepareDataError(CrisprHawkError):
    def __init__(self, value):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception
