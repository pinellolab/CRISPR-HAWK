""" """

from coordinate import Coordinate
from sequence import Sequence

from typing import List, Union


class Region:
    """Represents a genomic region with a sequence.

    Stores a DNA sequence and its associated genomic coordinates (contig, start, stop).

    Attributes:
        _sequence (Sequence): The DNA sequence of the region.
        _coordinates (Coordinate): The genomic coordinates of the region.
    """

    def __init__(self, sequence: Sequence, coord: Coordinate):
        """Initialize a Region object.

        Args:
            sequence: The DNA sequence of the region.
            coord: The genomic coordinates of the region.
        """
        self._sequence = sequence  # store sequence object
        # store contig name, start and stop coordinates of the extracted region
        self._coordinates = coord

    def __len__(self) -> int:
        """Return the length of the sequence.

        Returns the length of the DNA sequence in the region.

        Returns:
            The length of the sequence.
        """
        return len(self._sequence)

    def __eq__(self, region_query: "Region") -> bool:
        """Check if the current region is equal to another region.

        Checks if the provided region has the same sequence and coordinates as the
        current region.

        Args:
            region_query: The region to compare to.

        Returns:
            True if the regions are identical, False otherwise.
        """
        return (self._sequence == region_query.sequence) and (
            self._coordinates == region_query._coordinates
        )
    
    def __hash__(self) -> int:
        return hash((self._sequence.sequence, self._coordinates.contig, self._coordinates.start, self._coordinates.stop))

    def __str__(self) -> str:
        """Return a string representation of the region.

        Returns a FASTA-formatted string representation of the region, including the
        coordinates and the sequence.

        Returns:
            The string representation of the region.
        """
        return f">{str(self._coordinates)}\n{str(self._sequence)}"

    def __repr__(self) -> str:
        """Return a string representation of the region.

        Returns a string representation of the Region object, including the class name and
        the coordinates of the region.

        Returns:
            The string representation of the region.
        """
        return f"<{self.__class__.__name__} object; region={str(self._coordinates)}>"

    def __getitem__(self, idx: Union[int, slice]) -> str:
        """Return the nucleotide at the given index or slice.

        Returns the nucleotide(s) at the specified index or slice in the region's
        sequence.

        Args:
            idx: The index or slice to retrieve.

        Returns:
            The nucleotide or a slice of nucleotides.
        """
        return self._sequence[idx]

    def contain(self, region_query: "Region") -> bool:
        """Check if the current region fully contains another region.

        Checks if the provided region is entirely contained within the current region.
        Both regions must be on the same contig.

        Args:
            region_query: The region to check for containment.

        Returns:
            True if the current region fully contains the query region, False otherwise.

        Raises:
            TypeError: If region_query is not a Region object.
        """
        if not isinstance(region_query, self.__class__):
            raise TypeError(
                f"Full overlap check on input region can only be done on {self.__class__.__name__}"
            )
        # check whether the test region is contained in the current region
        return (
            self.contig == region_query.contig
            and self.start <= region_query.start
            and self.stop >= region_query.stop
        )

    def overlap(self, region_query: "Region") -> bool:
        """Check if the current region overlaps with another region.

        Checks if the provided region overlaps with the current region. Both regions
        must be on the same contig.

        Args:
            region_query: The region to check for overlap.

        Returns:
            True if the current region overlaps with the query region, False otherwise.

        Raises:
            TypeError: If region_query is not a Region object.
        """
        if not isinstance(region_query, self.__class__):
            raise TypeError(
                f"Full overlap check on input region can only be done on {self.__class__.__name__}"
            )
        # check whether the query region overlaps the current region
        return self.contig == region_query.contig and (
            self.start <= region_query.start <= self.stop
            or self.start <= region_query.stop <= self.stop
        )

    @property
    def contig(self) -> str:
        return self._coordinates.contig

    @property
    def start(self) -> int:
        return self._coordinates.start

    @property
    def stop(self) -> int:
        return self._coordinates.stop

    @property
    def sequence(self) -> Sequence:
        return self._sequence

    @property
    def sequence_ref(self) -> Sequence:
        return self._refseq

    @property
    def coordinates(self) -> Coordinate:
        return self._coordinates


class RegionList:
    """Represents a list of genomic regions.

    Stores a list of Region objects, providing methods for accessing,
    manipulating, and iterating over the regions.
    """

    def __init__(self, regions: List[Region]) -> None:
        """Initialize a RegionList object.

        Args:
            regions: A list of Region objects.
        """
        self._regions = regions  # list of region objects

    def __repr__(self) -> str:
        """Return a string representation of the region list.

        Returns a string representation of the RegionList object, including the
        class name and the number of regions in the list.

        Returns:
            The string representation of the region list.
        """
        return f"<{self.__class__.__name__} object; regions={len(self)}>"

    def __str__(self) -> str:
        """Return a string representation of the region list.

        Returns a string representation of the RegionList object, including the
        FASTA-formatted string representation of each region in the list.

        Returns:
            The string representation of the region list.
        """
        return "\n".join(str(region.coordinates) for region in self)

    def __len__(self) -> int:
        """Return the number of regions in the list.

        Returns the number of Region objects stored in the RegionList.

        Returns:
            The number of regions in the list.
        """
        return len(self._regions)

    def __iter__(self) -> "RegionListIterator":
        """Return an iterator over the regions in the list.

        Returns an iterator that allows iterating over the Region objects
        stored in the RegionList.

        Returns:
            An iterator over the regions in the list.
        """
        return RegionListIterator(self)

    def __getitem__(self, idx: Union[int, slice]) -> Region:
        """Return the region at the given index or slice.

        Returns the Region object(s) at the specified index or slice in the RegionList.

        Args:
            idx: The index or slice to retrieve.

        Returns:
            The Region object or a slice of Region objects.

        Raises:
            AttributeError: If the _regions attribute is missing.
            IndexError: If the index is out of range.
        """
        if not hasattr(self, "_regions"):  # always trace this error
            raise AttributeError(
                f"Missing _regions attribute on {self.__class__.__name__}"
            )
        try:
            return self._regions[idx]
        except IndexError as e:
            raise IndexError(f"Index {idx} out of range") from e

    def extend(self, regions: "RegionList") -> None:
        """Extend the region list with another RegionList.

        Extends the current RegionList by adding the Region objects from another
        RegionList to the end.

        Args:
            regions: The RegionList to extend the current list with.

        Raises:
            TypeError: If regions is not a RegionList object.
        """
        if not isinstance(regions, self.__class__):
            raise TypeError(
                f"Cannot extend {self.__class__.__name__} with objects of type {type(regions).__name__}"
            )
        self._regions.extend(regions)  # extend regions list

    def append(self, region: Region) -> None:
        """Append a region to the list.

        Appends a Region object to the end of the RegionList.

        Args:
            region: The Region object to append.

        Raises:
            TypeError: If region is not a Region object.
        """
        if not isinstance(region, Region):
            raise TypeError(
                f"Cannot append to {self.__class__.__name__} objects of type {type(region).__name__}"
            )
        self._regions.append(region)


class RegionListIterator:
    """Iterates over a list of Region objects.

    Provides a way to iterate over the Region objects stored in a RegionList.

    Attributes:
        _regions (RegionList): The RegionList object to iterate over.
        _index (int): The current index of the iterator.
    """

    def __init__(self, regions: RegionList) -> None:
        """Initialize the RegionListIterator.

        Args:
            regions: The RegionList object to iterate over.
        """
        self._regions = regions  # region list object to iterate over
        self._index = 0  # iterator index used over the list

    def __next__(self) -> Region:
        """Return the next Region object in the list.

        Returns the next Region object from the RegionList and advances the iterator.

        Returns:
            The next Region object in the list.

        Raises:
            StopIteration: If there are no more Region objects to iterate over.
        """
        if self._index < len(self._regions):
            result = self._regions[self._index]
            self._index += 1  # go to next position in list
            return result
        raise StopIteration  # stop iteration over regions list
