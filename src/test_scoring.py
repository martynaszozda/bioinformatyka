from scoring import sp_score_dna, sp_score_protein

def test_dna():
    w1 = "AGTCGTAG"
    w2 = "A-TCGTCG"
    w3 = "----GTAG"
    w4 = "-GTAG-AG"
    score = sp_score_dna([w1, w2, w3, w4], match=1, mismatch=-1, gap_open=-2, gap_extend=-1)
    print(f"DNA Score: {score} (Expected: -9)")
    assert score == -9, f"Expected -9, got {score}"

def test_protein():
    w1 = "AGTCGTAGNPST"
    w2 = "ASTCGTAG-PST"
    w3 = "-GTPG-AGN-ST"
    w4 = "A-NCGT--NP-T"
    w5 = "-G--ATA--PS-"
    score = sp_score_protein([w1, w2, w3, w4, w5], gap_open=-12, gap_extend=-1)
    print(f"Protein Score: {score} (Expected: -216)")
    assert score == -216, f"Expected -216, got {score}"

if __name__ == '__main__':
    test_dna()
    test_protein()
    print("All tests passed successfully.")
