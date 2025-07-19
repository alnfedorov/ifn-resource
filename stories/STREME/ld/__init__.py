from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
TAGS = RESULTS / "tags.pkl"


class streme:
    comparisons = RESULTS

    target = "target.fa"
    background = "background.fa"
    saveto = "streme"

    markov_order = 2
    seed = 123
    minw = 8
    maxw = 15
    hofract = 0.25
    thresh = 0.05


class single_cell:
    min_expression_tpm = 1

    monocyte = {
        "IFITM3", "APOBEC3A", "SIGLEC1", "CCR1", "SERPING1", "SSB", "TNFSF13B", "SRGN", "CD300E", "FPR3", "JUP",
        "TRIM14", "IL4I1", "CYSTM1", "PML", "ZNFX1", "RAB8A", "CCL8", "GMPR", "KLF6", "LGALS3BP", "AXL", "TOR1B",
        "DHX58", "MARCKS", "IFI27", "LILRB1", "RTCB", "IL1RN", "HESX1", "LYN", "NAPA", "RBX1", "CALM1", "CBR1",
        "LILRB4", "GIMAP6", "ADA2", "TAP1", "CCR5", "RGL1", "MSR1", "OTOF", "MOV10", "TENT5A", "PIK3AP1", "TUBA1C",
        "LILRB2", "FFAR2", "TRIM5", "CXCL10", "CD2AP", "ST3GAL5", "TPM3", "DNAJA1", "ATP10A", "CCL2", "PTP4A1", "SRC",
        "PARP10", "LYSMD2", "C15orf48", "SECTM1", "ATF5", "CNP", "LILRA5", "GIMAP8", "CMTR1", "SAMSN1", "MYD88",
        "PHACTR2", "NAGK", "SAMD4A", "CDKN1A", "CBWD2", "BATF2", "CBWD1", "ATF3", "DOCK8", "C3AR1", "TTC9", "ZC3HAV1",
        "LPAR6", "RRBP1", "MBOAT7", "TAP2", "NAMPT", "MLKL", "S100A11",
    }
    monocyte_and_lymphocytes = {
        "LY6E", "IFI44L", "IFI6", "IFIT1", "ISG15", "RSAD2", "MX1", "IFIT3", "MX2", "IFIT2", "OAS2", "XAF1", "IFI44",
        "OAS1", "EIF2AK2", "IFITM2", "CMPK2", "EPSTI1", "OAS3", "SAMD9", "HERC5", "OASL", "TNFSF10", "PLSCR1", "TYMP",
        "IFIH1", "DDX58", "PARP9", "SAMD9L", "IRF7", "IFI35", "TRIM22", "SPATS2L", "BST2", "RNF213", "DDX60L", "HELZ2",
        "LGALS9", "PARP14", "SP110", "HERC6", "DDX60", "USP18", "DTX3L", "PARP12", "SMCHD1", "GBP1", "LAP3", "ADAR",
        "UBE2L6", "RABGAP1L", "IFIT5", "DRAP1", "IFI16", "ISG20", "PNPT1", "STAT2", "C19orf66", "ZCCHC2", "PHF11",
        "IFITM1", "NMI", "TMEM123", "HLA-B", "SP100", "NT5C3A", "SHISA5", "N4BP1", "HLA-E", "CHMP5", "STAT1", "DYNLT1",
        "HLA-A", "UNC93B1", "HLA-C", "RTP4", "PSME2", "B2M", "MCL1", "PPM1K", "GIMAP4", "HSH2D", "H3F3B", "APOL6",
        "TRIM25",
    }
    lymphocyte = {
        "MT2A", "IRF9", "AC116407.2", "BISPR", "CD164", "PTMA", "GIMAP7", "CD48", "SLFN5", "CCDC85B", "GIMAP2", "NUB1",
        "GIMAP1", "PSMB9", "TRIM56", "USP15", "PSMB8", "SP140", "MYL12A", "TRIM38", "KPNB1", "TMSB10", "SAMHD1",
        "PSME1", "C4orf3", "CHST12", "CD38", "LAG3", "PLAC8", "FCHSD2", "SUB1", "CCSER1", "MNDA", "CACNA1A", "HLA-F",
        "SAT1", "ZBP1", "ALOX5AP", "TMSB4X", "TRANK1", "GPR155", "HLA-DPA1",
    }
