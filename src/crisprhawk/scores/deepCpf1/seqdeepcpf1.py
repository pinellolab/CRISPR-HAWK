"""
DeepCpf1 sequence-based scoring model and utilities.

This module defines the SeqDeepCpf1 neural network class and associated functions 
for preprocessing guide sequences, loading model weights, and computing DeepCpf1 
activity scores. It is used for predicting CRISPR guide RNA activity using the 
DeepCpf1 model.
"""

from typing import List
from torch import zeros, tensor, no_grad, float32, Tensor
from torch.nn import Linear

import torch.nn as nn

import h5py
import os

NTENCODING = {"A": 0, "C": 1, "G": 2, "T": 3}


class SeqDeepCpf1(nn.Module):
    """
    DeepCpf1 sequence-based neural network model for guide RNA activity prediction.

    This class implements a convolutional neural network architecture for predicting 
    CRISPR guide RNA activity using the DeepCpf1 model. It defines the network layers 
    and the forward pass for processing encoded guide sequences.

    Methods:
        forward(x): Perform a forward pass through the DeepCpf1 model.

    Args:
        x (Tensor): The input tensor representing encoded guide sequences.

    Returns:
        Linear: The output layer of the model producing activity predictions.
    """

    def __init__(self) -> None:
        super().__init__()
        self.conv = nn.Conv1d(4, 80, kernel_size=5)
        self.relu = nn.ReLU()
        self.pool = nn.AvgPool1d(kernel_size=2)
        self.flatten = nn.Flatten()
        self.drop1 = nn.Dropout(0.3)
        self.fc1 = nn.Linear(1200, 80)
        self.drop2 = nn.Dropout(0.3)
        self.fc2 = nn.Linear(80, 40)
        self.drop3 = nn.Dropout(0.3)
        self.fc3 = nn.Linear(40, 40)
        self.drop4 = nn.Dropout(0.3)
        self.output = nn.Linear(40, 1)

    def forward(self, x) -> Linear:
        x = self.conv(x)
        x = self.relu(x)
        x = self.pool(x)
        x = self.flatten(x.transpose(-1, -2))
        x = self.drop1(x)
        x = self.relu(self.fc1(x))
        x = self.drop2(x)
        x = self.relu(self.fc2(x))
        x = self.drop3(x)
        x = self.relu(self.fc3(x))
        x = self.drop4(x)
        return self.output(x)


# TODO: assign NA to guides with Ns
def preprocess(sequences: List[str]) -> Tensor:
    """
    Encode guide RNA sequences into a tensor for DeepCpf1 model input.

    Converts a list of nucleotide sequences into a one-hot encoded tensor suitable 
    for input to the DeepCpf1 neural network. All sequences must be of equal length.

    Args:
        sequences (List[str]): A list of guide RNA sequences.

    Returns:
        Tensor: A tensor of shape (num_sequences, 4, sequence_length) with one-hot 
            encoding.
    """
    seqnum = len(sequences)  # compute number of sequences
    seqlen = len(sequences[0])  # compute sequences length
    assert all(len(s) == seqlen for s in sequences)
    emb_matrix = zeros((seqnum, 4, seqlen), dtype=float32)  # create tensor
    for i, sequence in enumerate(sequences):
        for j, nt in enumerate(sequence):
            emb_matrix[i, NTENCODING[nt.upper()], j] = 1
    return emb_matrix


def load_deepcpf1_weights(model: SeqDeepCpf1) -> None:
    """
    Load pretrained weights into a SeqDeepCpf1 model instance.

    Loads the DeepCpf1 model weights from an HDF5 file and assigns them to the 
    corresponding layers of the provided SeqDeepCpf1 model. This function is 
    required before using the model for prediction.

    Args:
        model (SeqDeepCpf1): The DeepCpf1 model instance to load weights into.

    Returns:
        None
    """
    weightspath = os.path.join(os.path.abspath(os.path.dirname(__file__)), "weights")
    with h5py.File(os.path.join(weightspath, "Seq_deepCpf1_weights.h5"), mode="r") as f:
        w_conv = f["convolution1d_157"]["convolution1d_157_W"][:]  # type: ignore
        b_conv = f["convolution1d_157"]["convolution1d_157_b"][:]  # type: ignore
        model.conv.weight.data.copy_(
            tensor(w_conv).squeeze(1).permute(2, 1, 0).flip(-1)
        )
        model.conv.bias.data.copy_(tensor(b_conv))  # type: ignore
        model.fc1.weight.data.copy_(tensor(f["dense_490"]["dense_490_W"][:]).T)  # type: ignore
        model.fc1.bias.data.copy_(tensor(f["dense_490"]["dense_490_b"][:]))  # type: ignore
        model.fc2.weight.data.copy_(tensor(f["dense_491"]["dense_491_W"][:]).T)  # type: ignore
        model.fc2.bias.data.copy_(tensor(f["dense_491"]["dense_491_b"][:]))  # type: ignore
        model.fc3.weight.data.copy_(tensor(f["dense_492"]["dense_492_W"][:]).T)  # type: ignore
        model.fc3.bias.data.copy_(tensor(f["dense_492"]["dense_492_b"][:]))  # type: ignore
        model.output.weight.data.copy_(tensor(f["dense_493"]["dense_493_W"][:]).T)  # type: ignore
        model.output.bias.data.copy_(tensor(f["dense_493"]["dense_493_b"][:]))  # type: ignore


def compute_deepcpf1(model: SeqDeepCpf1, emb_matrix: Tensor) -> List[float]:
    """
    Compute DeepCpf1 activity scores for encoded guide sequences.

    Runs the provided DeepCpf1 model on the input tensor of encoded guide sequences 
    and returns the predicted activity scores as a list of floats.

    Args:
        model (SeqDeepCpf1): The DeepCpf1 model instance.
        emb_matrix (Tensor): The tensor of one-hot encoded guide sequences.

    Returns:
        List[float]: The predicted activity scores for each guide sequence.
    """
    with no_grad():
        return model(emb_matrix).squeeze().tolist()
