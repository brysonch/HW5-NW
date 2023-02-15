# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nws1s2 = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, seq1_align, seq2_align = nws1s2.align(seq1, seq2)

    assert score == 4.0
    assert seq1_align == "MYQR"
    assert seq2_align == "M-QR"

    ni = -np.inf

    align_check = np.array([[0, ni, ni, ni],
        [ni, 5, -11, -13],
        [ni, -12, 4, -8],
        [ni, -12, -1, 5],
        [ni, -14, -6, 4]]) 

    gapA_check = np.array([[-10, ni, ni, ni],
        [-11, -12,  -6, -7],
        [-12, -13, -14, -7],
        [-13, -14, -15, -12],
        [-14, -15, -16, -17]]) 

    gapB_check = np.array([[-10, -11, -12, -13], 
        [ni, -12, -13, -14], 
        [ni, -6, -14, -15], 
        [ni, -7, -7, -16], 
        [ni, -8, -8, -6]]) 

    assert (nws1s2._align_matrix == align_check).all()
    assert (nws1s2._gapA_matrix == gapA_check).all()
    assert (nws1s2._gapB_matrix == gapB_check).all()


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nws3s4 = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, seq3_align, seq4_align = nws3s4.align(seq3, seq4)

    assert score == 17
    assert seq3_align == "MAVHQLIRRP"
    assert seq4_align == "M---QLIRHP"


