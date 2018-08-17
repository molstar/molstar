/* tslint:disable */
/** Generated in 2018-08-17T12:05:55-07:00 */

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

export enum Type {
  GLOBAL = "GLOBAL",
  LOCAL = "LOCAL",
  PSEUDO = "PSEUDO"
}

export enum UnpublishedFlag {
  N = "N",
  Y = "Y"
}

export enum PdbxDiffrnProtocol {
  LAUE = "LAUE",
  MAD = "MAD",
  SINGLE_WAVELENGTH = "SINGLE_WAVELENGTH"
}

export enum PdbxMonochromaticOrLaueML {
  L = "L",
  M = "M"
}

export enum PdbxScatteringType {
  ELECTRON = "ELECTRON",
  NEUTRON = "NEUTRON",
  X_RAY = "X_RAY"
}

export enum Source {
  ELECTRON_MICROSCOPE = "ELECTRON_MICROSCOPE",
  FREE_ELECTRON_LASER = "FREE_ELECTRON_LASER",
  LIQUID_ANODE = "LIQUID_ANODE",
  NUCLEAR_REACTOR = "NUCLEAR_REACTOR",
  ROTATING_ANODE = "ROTATING_ANODE",
  SEALED_TUBE = "SEALED_TUBE",
  SPALLATION_SOURCE = "SPALLATION_SOURCE",
  SYNCHROTRON = "SYNCHROTRON"
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

export enum SrcMethod {
  MAN = "MAN",
  NAT = "NAT",
  SYN = "SYN"
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
    readonly pdbId?: string | null;
  };

  export type Query = {
    readonly __typename?: "Query";
    readonly assemblies?: ReadonlyArray<Assemblies | null> | null;
  };

  export type Assemblies = {
    readonly __typename?: "CoreAssembly";
    readonly assembly_id?: number | null;
    readonly rcsb_assembly_symmetry?: RcsbAssemblySymmetry | null;
  };

  export type RcsbAssemblySymmetry = {
    readonly __typename?: "RcsbAssemblySymmetry";
    readonly source?: string | null;
    readonly symmetry_features?: ReadonlyArray<SymmetryFeatures | null> | null;
  };

  export type SymmetryFeatures = {
    readonly __typename?: "SymmetryFeature";
    readonly clusters?: ReadonlyArray<Clusters | null> | null;
    readonly stoichiometry?: Stoichiometry | null;
    readonly symmetry?: Symmetry | null;
    readonly symmetry_axes?: ReadonlyArray<SymmetryAxes | null> | null;
    readonly type?: Type | null;
  };

  export type Clusters = {
    readonly __typename?: "Cluster";
    readonly members?: ReadonlyArray<string | null> | null;
    readonly avg_rmsd?: number | null;
  };

  export type Stoichiometry = {
    readonly __typename?: "Stoichiometry";
    readonly description?: string | null;
    readonly value?: ReadonlyArray<string | null> | null;
  };

  export type Symmetry = {
    readonly __typename?: "Symmetry";
    readonly space_group_name_h_m?: string | null;
  };

  export type SymmetryAxes = {
    readonly __typename?: "SymmetryAxis";
    readonly start?: ReadonlyArray<number | null> | null;
    readonly end?: ReadonlyArray<number | null> | null;
    readonly order?: number | null;
  };
}
