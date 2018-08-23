/* tslint:disable */
/** Generated in 2018-08-22T17:02:01-07:00 */

export enum PdbxLeavingAtomFlag {
  N = "N",
  Y = "Y"
}

export enum PdbxStereoConfig {
  N = "N",
  R = "R",
  S = "S"
}

export enum ExperimentalSupport {
  ASSAY_FOR_OLIGOMERIZATION = "ASSAY_FOR_OLIGOMERIZATION",
  CROSS_LINKING = "CROSS_LINKING",
  EQUILIBRIUM_CENTRIFUGATION = "EQUILIBRIUM_CENTRIFUGATION",
  FLUORESCENCE_RESONANCE_ENERGY_TRANSFER = "FLUORESCENCE_RESONANCE_ENERGY_TRANSFER",
  GEL_FILTRATION = "GEL_FILTRATION",
  HOMOLOGY = "HOMOLOGY",
  IMMUNOPRECIPITATION = "IMMUNOPRECIPITATION",
  ISOTHERMAL_TITRATION_CALORIMETRY = "ISOTHERMAL_TITRATION_CALORIMETRY",
  LIGHT_SCATTERING = "LIGHT_SCATTERING",
  MASS_SPECTROMETRY = "MASS_SPECTROMETRY",
  MICROSCOPY = "MICROSCOPY",
  NATIVE_GEL_ELECTROPHORESIS = "NATIVE_GEL_ELECTROPHORESIS",
  NONE = "NONE",
  SAXS = "SAXS",
  SCANNING_TRANSMISSION_ELECTRON_MICROSCOPY = "SCANNING_TRANSMISSION_ELECTRON_MICROSCOPY",
  SURFACE_PLASMON_RESONANCE = "SURFACE_PLASMON_RESONANCE"
}

export enum SymmetryFeatureType {
  GLOBAL = "GLOBAL",
  LOCAL = "LOCAL",
  PSEUDO = "PSEUDO"
}

export enum UnpublishedFlag {
  N = "N",
  Y = "Y"
}

export enum PdbxMonochromaticOrLaueMl {
  L = "L",
  M = "M"
}

export enum PdbxScatteringType {
  ELECTRON = "ELECTRON",
  NEUTRON = "NEUTRON",
  X_RAY = "X_RAY"
}

export enum RefSpace {
  REAL = "REAL",
  RECIPROCAL = "RECIPROCAL"
}

export enum SymmetryType {
  HELICAL = "HELICAL",
  POINT = "POINT",
  _2_D_CRYSTAL = "_2_D_CRYSTAL",
  _3_D_CRYSTAL = "_3_D_CRYSTAL"
}

export enum AggregationState {
  CELL = "CELL",
  FILAMENT = "FILAMENT",
  HELICAL_ARRAY = "HELICAL_ARRAY",
  PARTICLE = "PARTICLE",
  TISSUE = "TISSUE",
  _2_D_ARRAY = "_2_D_ARRAY",
  _3_D_ARRAY = "_3_D_ARRAY"
}

export enum ReconstructionMethod {
  CRYSTALLOGRAPHY = "CRYSTALLOGRAPHY",
  HELICAL = "HELICAL",
  SINGLE_PARTICLE = "SINGLE_PARTICLE",
  SUBTOMOGRAM_AVERAGING = "SUBTOMOGRAM_AVERAGING",
  TOMOGRAPHY = "TOMOGRAPHY"
}

export enum EmbeddingApplied {
  NO = "NO",
  YES = "YES"
}

export enum StainingApplied {
  NO = "NO",
  YES = "YES"
}

export enum VitrificationApplied {
  NO = "NO",
  YES = "YES"
}

export enum Type {
  NEGATIVE = "NEGATIVE",
  NONE = "NONE",
  POSITIVE = "POSITIVE"
}

export enum SrcMethod {
  MAN = "MAN",
  NAT = "NAT",
  SYN = "SYN"
}

export enum EntityType {
  MACROLIDE = "MACROLIDE",
  NON_POLYMER = "NON_POLYMER",
  POLYMER = "POLYMER",
  WATER = "WATER"
}

export enum EntityPolyType {
  CYCLIC_PSEUDO_PEPTIDE = "CYCLIC_PSEUDO_PEPTIDE",
  OTHER = "OTHER",
  PEPTIDE_NUCLEIC_ACID = "PEPTIDE_NUCLEIC_ACID",
  POLYDEOXYRIBONUCLEOTIDE = "POLYDEOXYRIBONUCLEOTIDE",
  POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID = "POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID",
  POLYPEPTIDE_D = "POLYPEPTIDE_D",
  POLYPEPTIDE_L = "POLYPEPTIDE_L",
  POLYRIBONUCLEOTIDE = "POLYRIBONUCLEOTIDE",
  POLYSACCHARIDE_D = "POLYSACCHARIDE_D",
  POLYSACCHARIDE_L = "POLYSACCHARIDE_L"
}

export enum RcsbHostOrganismSource {
  MMCIF = "MMCIF",
  NCBI = "NCBI",
  UNIPROT = "UNIPROT"
}

export enum RcsbMembraneSource {
  HOMOLOGY = "HOMOLOGY",
  MPSTRUCT = "MPSTRUCT"
}

export enum RcsbOrganismSource {
  MMCIF = "MMCIF",
  NCBI = "NCBI",
  UNIPROT = "UNIPROT"
}

export enum Level {
  _100 = "_100",
  _30 = "_30",
  _40 = "_40",
  _50 = "_50",
  _60 = "_60",
  _70 = "_70",
  _80 = "_80",
  _90 = "_90",
  _95 = "_95"
}

export enum PdbFormatCompatible {
  N = "N",
  Y = "Y"
}

export namespace AssemblySymmetry {
  export type Variables = {
    readonly pdbId: string;
  };

  export type Query = {
    readonly __typename?: "Query";
    readonly assemblies?: ReadonlyArray<Assemblies | null> | null;
  };

  export type Assemblies = {
    readonly __typename?: "CoreAssembly";
    readonly assembly_id: number;
    readonly rcsb_assembly_symmetry?: RcsbAssemblySymmetry | null;
  };

  export type RcsbAssemblySymmetry = {
    readonly __typename?: "RcsbAssemblySymmetry";
    readonly source: string;
    readonly symmetry_features: ReadonlyArray<SymmetryFeatures | null>;
  };

  export type SymmetryFeatures = {
    readonly __typename?: "SymmetryFeature";
    readonly symmetry: Symmetry;
    readonly clusters: ReadonlyArray<Clusters | null>;
    readonly stoichiometry: Stoichiometry;
    readonly symmetry_axes?: ReadonlyArray<SymmetryAxes | null> | null;
    readonly type: SymmetryFeatureType;
  };

  export type Symmetry = {
    readonly __typename?: "QuaternarySymmetry";
    readonly description: string;
    readonly value: string;
  };

  export type Clusters = {
    readonly __typename?: "Cluster";
    readonly members: ReadonlyArray<string | null>;
    readonly avg_rmsd?: number | null;
  };

  export type Stoichiometry = {
    readonly __typename?: "Stoichiometry";
    readonly description: string;
    readonly value: ReadonlyArray<string | null>;
  };

  export type SymmetryAxes = {
    readonly __typename?: "SymmetryAxis";
    readonly start: ReadonlyArray<number | null>;
    readonly end: ReadonlyArray<number | null>;
    readonly order?: number | null;
  };
}
