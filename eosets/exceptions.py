""" EOSets exceptions """


class EOSetsError(Exception):
    """EOSets base error"""

    pass


class IncompatibleProducts(EOSetsError):
    """Incompatible Products for some reason (non overlapping...)"""

    pass
