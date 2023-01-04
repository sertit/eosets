""" EOPairs exceptions """


class EOPairsProducts(Exception):
    """EOPairs error"""

    pass


class IncompatibleProducts(EOPairsProducts):
    """Incompatible Products for some reason (non overlapping...)"""

    pass
