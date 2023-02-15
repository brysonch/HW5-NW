# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    hg_score, h_galign, g_halign = nw.align(hs_seq, gg_seq)
    hm_score, h_malign, m_halign = nw.align(hs_seq, mm_seq)
    hb_score, h_balign, b_halign = nw.align(hs_seq, br_seq)
    ht_score, h_talign, t_halign = nw.align(hs_seq, tt_seq)

    scores = [[hg_score, "Gallus gallus"],
            [hm_score, "Mus musculus"],
            [hb_score, "Balaeniceps rex"],
            [ht_score, "Tursiops truncatus"]]
    sorted_scores = sorted(scores, key=lambda x: (-x[0]))
    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    print("Order of species most similar to least similar aligned with humans: \n")
    for i in sorted_scores:  
        print(i[0], " (score: ", i[1],")\n")
    
    

if __name__ == "__main__":
    main()
