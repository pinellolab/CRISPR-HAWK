"""This module defines the Coordinate class.

The Coordinate class represents a genomic interval. It is used when a BED file is
provided to select genomic regions for haplotype reconstruction.
"""


class Coordinate:
    """Genomic coordinate interval.

    Used when a bed file is provided in the input to select the
    genomic regions to consider when reconstructing haplotypes.

    Attributes:
        contig: Contig name.
        start: Start coordinate, including padding.
        stop: Stop coordinate, including padding.
        padding: Region padding length.
    """

    def __init__(self, contig: str, start: int, stop: int, padding: int) -> None:
        """Initialize a genomic coordinate interval.

        Used when a bed file is provided in the input to select the
        genomic regions to consider when reconstructing haplotypes.

        Attributes:
            contig: Contig name.
            start: Start coordinate, including padding.
            stop: Stop coordinate, including padding.
            padding: Region padding length.
        """
        # initialize a genomic coordinate interval
        # used when a bed file is provided in the input to select the
        # genomic regions to consider when reconstructing haplotypes
        if stop < start:
            raise ValueError("Stop < start coordinate")
        self._contig = contig  # set contig name
        self._start = start  # set start coordinate
        self._startp = max(0, start - padding)  # set start coordinate with padding
        self._stop = stop  # set stop coordinate
        self._stopp = stop + padding  # set stop coordinate with padding
        self._padding = padding  # set region padding length

    def __eq__(self, coord_query: object) -> bool:
        """Check if the current coordinate is equal to another coordinate.

        Checks if the provided coordinate has the same contig, start, and stop
        as the current coordinate.

        Args:
            coord_query: The coordinate to compare to.

        Returns:
            True if the coordinates are identical, False otherwise.
        """
        if not isinstance(coord_query, Coordinate):
            return NotImplemented
        return (
            (self._contig == coord_query.contig)
            and (self._startp == coord_query.start)
            and (self._stopp == coord_query.stop)
        )

    def __hash__(self) -> int:
        """Return a hash value for the Coordinate object.

        The hash is based on the contig, start, and stop values of the coordinate.
        """
        return hash((self.contig, self.start, self.stop))

    def __repr__(self) -> str:
        """Return a string representation of the Coordinate object.

        The representation includes the class name, the coordinate (contig:start-stop),
        and the padding.

        Returns:
            A string representation of the Coordinate object.
        """
        # retrieve original start and stop
        start = self._start + self._padding
        stop = self._stop - self._padding
        return f"<{self.__class__.__name__} object; coordinate={self._contig}:{start}-{stop}; padding={self._padding}>"

    def __str__(self) -> str:
        """Return a string representation of the Coordinate object.

        The representation includes the coordinate (contig:start-stop).

        Returns:
            A string representation of the Coordinate object.
        """
        # retrieve original start and stop
        return f"{self._contig}:{self._start}-{self._stop}"

    def contains(self, query: object) -> bool:
        """Check if the coordinate contains another coordinate.

        Determines whether the provided coordinate is fully contained within the
        current coordinate interval.

        Args:
            query: The coordinate to check for containment.

        Returns:
            True if the query coordinate is contained within this coordinate,
                False otherwise.
        """
        if not isinstance(query, Coordinate):
            return NotImplemented
        return (
            self._contig == query.contig
            and self._start <= query.start
            and self._stop >= query.stop
        )

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._startp

    @property
    def stop(self) -> int:
        return self._stopp

    @property
    def startp(self) -> int:
        return self._start

    @property
    def stopp(self) -> int:
        return self._stop
