# Loading libraries
library(reticulate)
library(caret)
library(dplyr)

# Reading previously generated dataframe, which contains SMILES notations and extracting them to a separate variable
data_SMILES <- read.csv("smiles_hc50.csv")
smiles_codes <- data_SMILES$SMILES

# Loading a Python interpreter from an Anaconda environment to run 'rdkit' package in R environment
use_python("path to python.exe")

# Importing 'rdkit' modules for chemoinformatics
rdkit <- import("rdkit")
Chem <- rdkit$Chem
Descriptors <- rdkit$Chem$Descriptors
AllChem <- rdkit$Chem$AllChem

# Setting up empty character vectors to store molecular descriptors and MACCS keys fingerprints as a single string
MolarWeights <- character(length(smiles_codes))
logPs <- character(length(smiles_codes))
TPSAs <- character(length(smiles_codes))
Ipcs <- character(length(smiles_codes))
MACCS_keys <- character(length(smiles_codes))
HeavyAtomCounts <- character(length(smiles_codes))
HeavyAtoMolWts <- character(length(smiles_codes))
NHOHCounts <- character(length(smiles_codes))
NOCounts <- character(length(smiles_codes))
NumHAcceptor <- character(length(smiles_codes))
NumHDonor <- character(length(smiles_codes))
NumHeteroatoms <- character(length(smiles_codes))
NumRotatableBond <- character(length(smiles_codes))
NumValenceElectron <- character(length(smiles_codes))
NumRadicalElectron <- character(length(smiles_codes))
NumAmideGroups <- character(length(smiles_codes))
NumAromaticRing <- character(length(smiles_codes))
NumSaturatedRing <- character(length(smiles_codes))
NumAliphaticRing <- character(length(smiles_codes))
RingCounts <- character(length(smiles_codes))
FractionCSP3s <- character(length(smiles_codes))
NumAromaticHeteroCycle <- character(length(smiles_codes))
NumAromaticCarboCycle <- character(length(smiles_codes))
NumSaturatedHeteroCycle <- character(length(smiles_codes))
NumSaturatedCarboCycle <- character(length(smiles_codes))
NumAliphaticHeteroCycle <- character(length(smiles_codes))
NumAliphaticCarboCycle <- character(length(smiles_codes))
FragmentsAlCOO <- character(length(smiles_codes))
FragmentsAlOH <- character(length(smiles_codes))
FragmentsAlOHnotert <- character(length(smiles_codes))
FragmentsAr_N <- character(length(smiles_codes))
FragmentsArCOO <- character(length(smiles_codes))
FragmentsArN <- character(length(smiles_codes))
FragmentsArNH <- character(length(smiles_codes))
FragmentsArOH <- character(length(smiles_codes))
FragmentsCOO <- character(length(smiles_codes))
FragmentsCOO2 <- character(length(smiles_codes))
FragmentsC_O_no_COO <- character(length(smiles_codes))
FragmentsCS <- character(length(smiles_codes))
FragmentsHOCCN <- character(length(smiles_codes))
FragmentsImine <- character(length(smiles_codes))
FragmentsNH0 <- character(length(smiles_codes))
FragmentsNH1 <- character(length(smiles_codes))
FragmentsNH2 <- character(length(smiles_codes))
FragmentsN_O <- character(length(smiles_codes))
FragmentsXCCNR1 <- character(length(smiles_codes))
FragmentsXCCNR2 <- character(length(smiles_codes))
FragmentsSH <- character(length(smiles_codes))
FragmentsAmide <- character(length(smiles_codes))
FragmentsImide <- character(length(smiles_codes))
FragmentsAniline <- character(length(smiles_codes))
FragmentsArylMethyl <- character(length(smiles_codes))
FragmentsAzide <- character(length(smiles_codes))
FragmentsAzo <- character(length(smiles_codes))
FragmentsBarbitur <- character(length(smiles_codes))
FragmentsBenzene <- character(length(smiles_codes))
FragmentsBenzodiazepine <- character(length(smiles_codes))
FragmentsBicyclic <- character(length(smiles_codes))
FragmentsDiazo <- character(length(smiles_codes))
FragmentsDihydropyridine <- character(length(smiles_codes))
FragmentsEpoxideRings <- character(length(smiles_codes))
FragmentsEster <- character(length(smiles_codes))
FragmentsFuran <- character(length(smiles_codes))
FragmentsGuanido <- character(length(smiles_codes))
FragmentsHalogen <- character(length(smiles_codes))
FragmentsHydrazine <- character(length(smiles_codes))
FragmentsHydrazone <- character(length(smiles_codes))
FragmentsImidazole <- character(length(smiles_codes))
FragmentsIsocyan <- character(length(smiles_codes))
FragmentsIsothiocyan <- character(length(smiles_codes))
FragmentsKetone <- character(length(smiles_codes))
FragmentsKetoneTopliss <- character(length(smiles_codes))
FragmentsLactam <- character(length(smiles_codes))
FragmentsLactone <- character(length(smiles_codes))
FragmentsMethoxy <- character(length(smiles_codes))
FragmentsMorpholine <- character(length(smiles_codes))
FragmentsNitrile <- character(length(smiles_codes))
FragmentsNitro <- character(length(smiles_codes))
FragmentsNitroArom <- character(length(smiles_codes))
FragmentsNitroAromNonortho <- character(length(smiles_codes))
FragmentsNitroso <- character(length(smiles_codes))
FragmentsOxazole <- character(length(smiles_codes))
FragmentsOxime <- character(length(smiles_codes))
FragmentsParaHydroxylation <- character(length(smiles_codes))
FragmentsPhenol <- character(length(smiles_codes))
FragmentsPhenolNoOrthoHbonds <- character(length(smiles_codes))
FragmentsPhosAcid <- character(length(smiles_codes))
FragmentsPhosEster <- character(length(smiles_codes))
FragmentsPierdine <- character(length(smiles_codes))
FragmentsPiperzine <- character(length(smiles_codes))
FragmentsPriamide <- character(length(smiles_codes))
FragmentsPrisulfonamd <- character(length(smiles_codes))
FragmentsPyridine <- character(length(smiles_codes))
FragmentsSulfide <- character(length(smiles_codes))
FragmentsSulfonamd <- character(length(smiles_codes))
FragmentsSulfone <- character(length(smiles_codes))
FragmentsTermAcetylene <- character(length(smiles_codes))
FragmentsTetrazole <- character(length(smiles_codes))
FragmentsThiazole <- character(length(smiles_codes))
FragmentsThiocyan <- character(length(smiles_codes))
FragmentsThiophene <- character(length(smiles_codes))
FragmentsUnbrchAlkene <- character(length(smiles_codes))
FragmentsUrea <- character(length(smiles_codes))
FragmentsAldehydes <- character(length(smiles_codes))

# A loop that iterates over SMILES notations present in previously loaded dataframe and converting them to a molecule object. It then calculates specific molecular descriptors and MACCS keys fingerprints for each molecule and at the very end storing them in respective vectors
for (i in 1:length(smiles_codes)) {
  mol <- Chem$MolFromSmiles(smiles_codes[i])
  MolarWeights[i] <- Descriptors$MolWt(mol)
  logPs[i] <- Descriptors$MolLogP(mol)
  TPSAs[i] <- Descriptors$TPSA(mol)
  Ipcs[i] <- Descriptors$Ipc(mol)
  HeavyAtomCounts[i] <- Descriptors$HeavyAtomCount(mol)
  HeavyAtoMolWts[i] <- Descriptors$HeavyAtomMolWt(mol)
  NHOHCounts[i] <- Descriptors$NHOHCount(mol)
  NumHAcceptor[i] <- Descriptors$NumHAcceptors(mol)
  NumHDonor[i] <- Descriptors$NumHDonors(mol)
  NumHeteroatoms[i] <- Descriptors$NumHeteroatoms(mol)
  NumRotatableBond[i] <- Descriptors$NumRotatableBonds(mol)
  NumValenceElectron[i] <- Descriptors$NumValenceElectrons(mol)
  NumRadicalElectron[i] <- Descriptors$NumRadicalElectrons(mol)
  NumAmideGroups[i] <- Descriptors$fr_amidine(mol)
  NumAromaticRing[i] <- Descriptors$NumAromaticRings(mol)
  NumSaturatedRing[i] <- Descriptors$NumSaturatedRings(mol)
  NumAliphaticRing[i] <- Descriptors$NumAliphaticRings(mol)
  RingCounts[i] <- Descriptors$RingCount(mol)
  FractionCSP3s[i] <- Descriptors$FractionCSP3(mol)
  NumAromaticHeteroCycle[i] <- Descriptors$NumAromaticHeterocycles(mol)
  NumAromaticCarboCycle[i] <- Descriptors$NumAromaticCarbocycles(mol)
  NumSaturatedHeteroCycle[i] <- Descriptors$NumSaturatedHeterocycles(mol)
  NumSaturatedCarboCycle[i] <- Descriptors$NumSaturatedCarbocycles(mol)
  NumAliphaticHeteroCycle[i] <- Descriptors$NumAliphaticHeterocycles(mol)
  NumAliphaticCarboCycle[i] <- Descriptors$NumAliphaticCarbocycles(mol)
  FragmentsAlCOO[i] <- Descriptors$fr_Al_COO(mol)
  FragmentsAlOH[i] <- Descriptors$fr_Al_OH(mol)
  FragmentsAlOHnotert[i] <- Descriptors$fr_Al_OH_noTert(mol)
  FragmentsAr_N[i] <- Descriptors$fr_Ar_N(mol)
  FragmentsArCOO[i] <- Descriptors$fr_Ar_COO(mol)
  FragmentsArN[i] <- Descriptors$fr_ArN(mol)
  FragmentsArNH[i] <- Descriptors$fr_Ar_NH(mol)
  FragmentsArOH[i] <- Descriptors$fr_Ar_OH(mol)
  FragmentsCOO[i] <- Descriptors$fr_COO(mol)
  FragmentsCOO2[i] <- Descriptors$fr_COO2(mol)
  FragmentsC_O_no_COO[i] <- Descriptors$fr_C_O_noCOO(mol)
  FragmentsCS[i] <- Descriptors$fr_C_S(mol)
  FragmentsHOCCN[i] <- Descriptors$fr_HOCCN(mol)
  FragmentsImine[i] <- Descriptors$fr_Imine(mol)
  FragmentsNH0[i] <- Descriptors$fr_NH0(mol)
  FragmentsNH1[i] <- Descriptors$fr_NH1(mol)
  FragmentsNH2[i] <- Descriptors$fr_NH2(mol)
  FragmentsN_O[i] <- Descriptors$fr_N_O(mol)
  FragmentsXCCNR1[i] <- Descriptors$fr_Ndealkylation1(mol)
  FragmentsXCCNR2[i] <- Descriptors$fr_Ndealkylation2(mol)
  FragmentsSH[i] <- Descriptors$fr_SH(mol)
  FragmentsAmide[i] <- Descriptors$fr_amide(mol)
  FragmentsImide[i] <- Descriptors$fr_imide(mol)
  FragmentsAniline[i] <- Descriptors$fr_aniline(mol)
  FragmentsArylMethyl[i] <- Descriptors$fr_aryl_methyl(mol)
  FragmentsAzide[i] <- Descriptors$fr_azide(mol)
  FragmentsAzo[i] <- Descriptors$fr_azo(mol)
  FragmentsBarbitur[i] <- Descriptors$fr_barbitur(mol)
  FragmentsBenzene[i] <- Descriptors$fr_benzene(mol)
  FragmentsBenzodiazepine[i] <- Descriptors$fr_benzodiazepine(mol)
  FragmentsBicyclic[i] <- Descriptors$fr_bicyclic(mol)
  FragmentsDiazo[i] <- Descriptors$fr_diazo(mol)
  FragmentsDihydropyridine[i] <- Descriptors$fr_dihydropyridine(mol)
  FragmentsEpoxideRings[i] <- Descriptors$fr_epoxide(mol)
  FragmentsEster[i] <- Descriptors$fr_ester(mol)
  FragmentsFuran[i] <- Descriptors$fr_furan(mol)
  FragmentsGuanido[i] <- Descriptors$fr_guanido(mol)
  FragmentsHalogen[i] <- Descriptors$fr_halogen(mol)
  FragmentsHydrazine[i] <- Descriptors$fr_hdrzine(mol)
  FragmentsHydrazone[i] <- Descriptors$fr_hdrzone(mol)
  FragmentsImidazole[i] <- Descriptors$fr_imidazole(mol)
  FragmentsIsocyan[i] <- Descriptors$fr_isocyan(mol)
  FragmentsIsothiocyan[i] <- Descriptors$fr_isothiocyan(mol)
  FragmentsKetone[i] <- Descriptors$fr_ketone(mol)
  FragmentsKetoneTopliss[i] <- Descriptors$fr_ketone_Topliss(mol)
  FragmentsLactam[i] <- Descriptors$fr_lactam(mol)
  FragmentsLactone[i] <- Descriptors$fr_lactone(mol)
  FragmentsMethoxy[i] <- Descriptors$fr_methoxy(mol)
  FragmentsMorpholine[i] <- Descriptors$fr_morpholine(mol)
  FragmentsNitrile[i] <- Descriptors$fr_nitrile(mol) 
  FragmentsNitro[i] <- Descriptors$fr_nitro(mol)
  FragmentsNitroArom[i] <- Descriptors$fr_nitro_arom(mol)
  FragmentsNitroAromNonortho[i] <- Descriptors$fr_nitro_arom_nonortho(mol)
  FragmentsNitroso[i] <- Descriptors$fr_nitroso(mol)
  FragmentsOxazole[i] <- Descriptors$fr_oxazole(mol)
  FragmentsOxime[i] <- Descriptors$fr_oxime(mol)
  FragmentsParaHydroxylation[i] <- Descriptors$fr_para_hydroxylation(mol)
  FragmentsPhenol[i] <- Descriptors$fr_phenol(mol)
  FragmentsPhenolNoOrthoHbonds[i] <- Descriptors$fr_phenol_noOrthoHbond(mol)
  FragmentsPhosAcid[i] <- Descriptors$fr_phos_acid(mol)
  FragmentsPhosEster[i] <- Descriptors$fr_phos_ester(mol)
  FragmentsPierdine[i] <- Descriptors$fr_piperdine(mol)
  FragmentsPiperzine[i] <- Descriptors$fr_piperzine(mol)
  FragmentsPriamide[i] <- Descriptors$fr_priamide(mol)
  FragmentsPrisulfonamd[i] <- Descriptors$fr_prisulfonamd(mol)
  FragmentsPyridine[i] <- Descriptors$fr_pyridine(mol)
  FragmentsSulfide[i] <- Descriptors$fr_sulfide(mol)
  FragmentsSulfonamd[i] <- Descriptors$fr_sulfonamd(mol)
  FragmentsSulfone[i] <- Descriptors$fr_sulfone(mol)
  FragmentsTermAcetylene[i] <- Descriptors$fr_term_acetylene(mol)
  FragmentsTetrazole[i] <- Descriptors$fr_tetrazole(mol)
  FragmentsThiazole[i] <- Descriptors$fr_thiazole(mol)
  FragmentsThiocyan[i] <- Descriptors$fr_thiocyan(mol)
  FragmentsThiophene[i] <- Descriptors$fr_thiophene(mol)
  FragmentsUnbrchAlkene[i] <- Descriptors$fr_unbrch_alkane(mol)
  FragmentsUrea[i] <- Descriptors$fr_urea(mol)
  FragmentsAldehydes[i] <- Descriptors$fr_aldehyde(mol)
  
  maccs <- AllChem$GetMACCSKeysFingerprint(mol)
  MACCS_keys[i] <- rdkit$DataStructs$BitVectToText(maccs)
}

# Splitting MACCS keys from a single string into separate columns for each bit and storing them in a dataframe (each column named from MACCS1 to MACCS167)
maccs_matrix <- do.call(rbind, strsplit(MACCS_keys, split = "")) #
maccs_df <- as.data.frame(maccs_matrix, stringsAsFactors = FALSE) #
colnames(maccs_df) <- paste0("MACCS", 1:167)

# Generating a dataframe to store SMILES notations, HC50, calculated molecular descriptors and MACCS keys fingerprints
descriptors_dataframe <- data.frame(
  SMILES = smiles_codes,
  HC50.exp = data_SMILES$HC50,
  Molar.Weight = MolarWeights,
  LogP = logPs,
  TPSA = TPSAs,
  Ipc = Ipcs,
  Heavy.Atom.Count = HeavyAtomCounts,
  Heavy.AtomMol.Weight = HeavyAtoMolWts,
  NHOH.Count = NHOHCounts,
  H.Acceptors = NumHAcceptor,
  H.Donors = NumHDonor,
  Heteroatoms = NumHeteroatoms,
  Rotatable.Bonds = NumRotatableBond,
  Valence.Electrons = NumValenceElectron,
  Radical.Electrons = NumRadicalElectron,
  Amide.Groups = NumAmideGroups,
  Aromatic.Rings = NumAromaticRing,
  Aliphatic.Rings = NumAliphaticRing,
  Ring.Count = RingCounts,
  CSP3s = FractionCSP3s,
  Aromatic.Heterocycles = NumAromaticHeteroCycle,
  Aromatic.Carbocycles = NumAromaticCarboCycle,
  Saturated.Heterocycles = NumSaturatedHeteroCycle,
  Saturated.Carbocycles = NumSaturatedCarboCycle,
  Aliphatic.Heterocycles = NumAliphaticHeteroCycle,
  Aliphatic.Carbocycles = NumAliphaticCarboCycle,
  Aliphatic.Carboxylic.Acids = FragmentsAlCOO,
  Aliphatic.Hydroxyl.Groups = FragmentsAlOH,
  Aliphatic.Hydroxyl.Groups.Exc.TertOH = FragmentsAlOHnotert,
  Aromatic.Nitrogens = FragmentsAr_N,
  Aromatic.Carboxylic.Acide = FragmentsArCOO,
  Functional.N.Groups.Att.To.Aromatics = FragmentsArN,
  Aromatic.Amines = FragmentsArNH,
  Aromatic.Hydroxyl.Groups = FragmentsArOH,
  Carboxylic.Acids1 = FragmentsCOO,
  Carboxylic.Acids2 = FragmentsCOO2,
  Carbonyl.O.Exc.COOH = FragmentsC_O_no_COO,
  Thiocarbonyl = FragmentsCS,
  HOCCN.tert = FragmentsHOCCN,
  Imines = FragmentsImine,
  Tertiary.Amines = FragmentsNH0,
  Secondary.Amines = FragmentsNH1,
  Primary.Amines =FragmentsNH2,
  Hydroxylamine.Groups = FragmentsN_O,
  XCCNR1.Groups = FragmentsXCCNR1,
  XCCNR2.Groups = FragmentsXCCNR2,
  Thiol.Groups = FragmentsSH,
  Aldehydes = FragmentsAldehydes,
  Amides = FragmentsAmide,
  Anilines = FragmentsAniline,
  Aryl.Methyl.Sites = FragmentsArylMethyl,
  Azide.Groups = FragmentsAzide,
  Azo.Groups = FragmentsAzo,
  Barbiturate.Groups = FragmentsBarbitur,
  Benzene.Rings = FragmentsBenzene,
  Benzodiazepines.W.O.Additional.Fused.Rings = FragmentsBenzodiazepine,
  Bicyclic = FragmentsBicyclic,
  Diazo.Groups = FragmentsDiazo,
  Dihydropyridines = FragmentsDihydropyridine,
  Epoxide.Groups = FragmentsEpoxideRings,
  Esters = FragmentsEster,
  Furan.Rings = FragmentsFuran,
  Guanidine.Groups = FragmentsGuanido,
  Halogens = FragmentsHalogen,
  Hydrazine.Groups = FragmentsHydrazine,
  Hydrazone.Groups = FragmentsHydrazone,
  Imidazole.Rings = FragmentsImidazole,
  Imide.Groups = FragmentsImide,
  Isocyanates = FragmentsIsocyan,
  Isothiocyanates = FragmentsIsothiocyan,
  Ketones = FragmentsKetone,
  Ketone.Exc = FragmentsKetoneTopliss,
  Beta.Lactams = FragmentsLactam, 
  Cyclic.Esters = FragmentsLactone,
  Methoxy.Groups = FragmentsMethoxy, 
  Morpholine.Rings = FragmentsMorpholine,
  Nitriles = FragmentsNitrile,
  Nitro.Groups = FragmentsNitro,
  Nitro.Benzene.Ring.Substituents = FragmentsNitroArom,
  Non.Ortho.Benzene.Rings.Substituents = FragmentsNitroAromNonortho,
  Nitroso.Groups.Excl.NO2 = FragmentsNitroso,
  Oxazole.Rings = FragmentsOxazole,
  Oxime.Groups = FragmentsOxime,
  Para.Hydroxylation.Sites = FragmentsParaHydroxylation,
  Phenols = FragmentsPhenol, 
  Phenolic.OH.Excl = FragmentsPhenolNoOrthoHbonds,
  Phopshoric.Acid.Groups = FragmentsPhosAcid,
  Phosphoric.Ester.Groups = FragmentsPhosEster,
  Piperzine.Rings = FragmentsPiperzine,
  Primary.Amides = FragmentsPriamide,
  Primary.Sulfonamides = FragmentsPrisulfonamd,
  Pyridine.Rings = FragmentsPyridine,
  Thioether = FragmentsSulfide,
  Sulfonamides = FragmentsSulfonamd,
  Sulfone.Groups = FragmentsSulfone,
  Terminal.Acetylenes = FragmentsTermAcetylene,
  Tetrazole.Rings = FragmentsTetrazole,
  Thiazole.Rings = FragmentsThiazole,
  Thiocyanates = FragmentsThiocyan,
  Thiophene.Rings = FragmentsThiophene,
  Unbranched.Alkanes = FragmentsUnbrchAlkene,
  Urea.Groups = FragmentsUrea
)

#Extracting SMILES notations and experimental values of HC50 into two variables
SMILES <- descriptors_dataframe$SMILES
HC50.exp <- descriptors_dataframe$HC50.exp

#Removing SMILES notations and experimental values of HC50 from a training data set
descriptors_dataframe_s_hc50 <- select(descriptors_dataframe,
                                       -SMILES,
                                       -HC50.exp)

#Converting training data set into 'numeric' type so further calculcations can be performed
descriptors_dataframe_s_hc50 <- as.data.frame(lapply(descriptors_dataframe_s_hc50, 
                                                     as.numeric))

#Preparing model which performs range standardization to a scale between 0 and 1
center_desc_model <- preProcess(descriptors_dataframe_s_hc50,
                                rangeBounds = c(0, 1),
                                method = "range")

#Scaling descriptors to a scale between 0 and 1
descriptors_dataframe_scaled <- predict(center_desc_model,
                                        descriptors_dataframe_s_hc50)

# Combining the previously generated dataframe of scale molecular descriptors and MACCS keys fingerprints into a single dataframe
descriptors_dataframe_scaled_maccs <- cbind(SMILES, 
                                            HC50.exp, 
                                            descriptors_dataframe_scaled, 
                                            maccs_df)

# Saving the previously generated dataframe as a CSV (comma-separated values) file
write.csv(descriptors_dataframe_scaled_maccs, 
          file = "descriptors_dataframe_desc_center_0to1_maccs.csv",
          row.names = TRUE,
          quote = FALSE)