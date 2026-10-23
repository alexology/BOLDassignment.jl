# Opt-in test against the real BOLD services (BOLDASSIGNMENT_LIVE=1).
# Public COI-5P record CASBI15164-26, Onthophagus taurus.
const LIVE_SEQ = "CCTTATCCTCAAATATTGCCCATGGAGGGGCCTCTGTAGATTTAGCAATTTTTAGACTACATTTGGCAGG" *
    "AATTTCATCAATTTTAGGGGCTGTAAATTTTATTACTACTGTTATTAATATACGAACTACAGGAATAACTTTTGATCGAATAC" *
    "CTTTATTTGCTTGAGCTGTAGCTATTACAGCATTACTCCTTCTTCTTTCCTTACCAGTTTTAGCAGGAGCAATTACAATATTA" *
    "TTAACTGATCGAAATCTTAATACATCTTTTTTTGATCCTATCGGAGGTGGAGATCCAATTTTATATCAACATTTATTT"

live_hits = bold_identify([LIVE_SEQ, "ACGT"^60])
@test nrow(live_hits) > 1
@test maximum(skipmissing(live_hits.similarity)) > 99
@test any(isequal("Onthophagus taurus"), live_hits.species)
tax = bold_assign(DataFrame(ASV = [LIVE_SEQ]); hits = live_hits, min_agreement = 90)
@test tax.genus[1] == "Onthophagus"
records = bold_records(["CASBI15164-26"]; fields = [:processid, :species, :bin_uri])
@test records.species == ["Onthophagus taurus"]
