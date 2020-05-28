/* eslint-disable */
export type Maybe<T> = T | null;

// Generated in 2020-05-27T17:35:45-07:00

/** All built-in and custom scalars, mapped to their actual values */
export type Scalars = {
  ID: string;
  String: string;
  Boolean: boolean;
  Int: number;
  Float: number;
  /** Built-in scalar representing an instant in time */
  Date: any;
  /** Unrepresentable type */
  UNREPRESENTABLE: any;
};

export type PdbxDatabasePdbObsSpr = {
  /**
   * The date of replacement.
   * 
   * Examples:
   * 1997-03-30
   */
  readonly date?: Maybe<Scalars['Date']>;
  /** Details related to the replaced or replacing entry. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Identifier for the type of obsolete entry to be added to this entry.
   * 
   * Allowable values:
   * OBSLTE, SPRSDE
   */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * The new PDB identifier for the replaced entry.
   * 
   * Examples:
   * 2ABC
   */
  readonly pdb_id: Scalars['String'];
  /**
   * The PDB identifier for the replaced (OLD) entry/entries.
   * 
   * Examples:
   * 3ABC
   */
  readonly replace_pdb_id: Scalars['String'];
};

export type PdbxStructSpecialSymmetry = {
  /**
   * Part of the identifier for the molecular component.
   * 
   * This data item is a pointer to _atom_site.pdbx_PDB_model_num in the
   * ATOM_SITE category.
   */
  readonly PDB_model_num?: Maybe<Scalars['Int']>;
  /**
   * Part of the identifier for the molecular component.
   * 
   *  This data item is a pointer to _atom_site.auth_seq_id in the
   *  ATOM_SITE category.
   */
  readonly auth_seq_id?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_struct_special_symmetry.id must uniquely identify
   *  each item in the PDBX_STRUCT_SPECIAL_SYMMETRY list.
   * 
   *  This is an integer serial number.
   */
  readonly id: Scalars['Int'];
  /**
   * Part of the identifier for the molecular component.
   * 
   *  This data item is a pointer to _atom_site.label_asym_id in the
   *  ATOM_SITE category.
   */
  readonly label_asym_id?: Maybe<Scalars['String']>;
  /**
   * Part of the identifier for the molecular component.
   * 
   *  This data item is a pointer to _atom_site.label_comp_id in the
   *  ATOM_SITE category.
   */
  readonly label_comp_id?: Maybe<Scalars['String']>;
};

export type PdbxChemCompAudit = {
  /**
   * The action associated with this audit record.
   * 
   * Allowable values:
   * Create component, Initial release, Modify aromatic_flag, Modify atom id, Modify charge, Modify component atom id, Modify component comp_id, Modify coordinates, Modify descriptor, Modify formal charge, Modify formula, Modify identifier, Modify internal type, Modify leaving atom flag, Modify linking type, Modify model coordinates code, Modify name, Modify one letter code, Modify parent residue, Modify processing site, Modify subcomponent list, Modify synonyms, Modify value order, Obsolete component, Other modification
   */
  readonly action_type?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _chem_comp.id in the CHEM_COMP
   *  category.
   */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** The date associated with this audit record. */
  readonly date?: Maybe<Scalars['Date']>;
  /**
   * Additional details decribing this change.
   * 
   * Examples:
   * Added C14 as a leaving atom.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * This data item is an ordinal index for the
   *  PDBX_CHEM_COMP_AUDIT category.
   */
  readonly ordinal: Scalars['Int'];
};

/** Query root */
export type Query = {
  /** Get a PDB polymer entity instance (chain), given the PDB ID and ENTITY INSTANCE ID. Here ENTITY INSTANCE ID identifies structural element in the asymmetric unit, e.g. 'A', 'B', etc. */
  readonly polymer_entity_instance?: Maybe<CorePolymerEntityInstance>;
  /** Get a list of assemblies given the list of ASSEMBLY IDs. Here an ASSEMBLY ID is a compound identifier that includes entry_id and assembly_id separated by '-', e.g. 1XXX-1. */
  readonly assemblies?: Maybe<ReadonlyArray<Maybe<CoreAssembly>>>;
  /** Get a list of PDB non-polymer entities given a list of ENTITY IDs. Here ENTITY ID is a compound identifier that includes entry_id and entity_id separated by '_', e.g. 1XXX_1. */
  readonly nonpolymer_entities?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntity>>>;
  /** Get a PDB non-polymer entity instance (chain), given the PDB ID and ENTITY INSTANCE ID. Here ENTITY INSTANCE ID identifies structural element in the asymmetric unit, e.g. 'A', 'B', etc. */
  readonly nonpolymer_entity_instance?: Maybe<CoreNonpolymerEntityInstance>;
  /** Get a list of PDB polymer entities given a list of ENTITY IDs. Here ENTITY ID is a compound identifier that includes entry_id and entity_id separated by '_', e.g. 1XXX_1. */
  readonly polymer_entities?: Maybe<ReadonlyArray<Maybe<CorePolymerEntity>>>;
  /** Get a PDB polymer entity, given the PDB ID and ENTITY ID. Here ENTITY ID is a '1', '2', '3', etc. */
  readonly polymer_entity?: Maybe<CorePolymerEntity>;
  /** Get a chemical component given the CHEMICAL COMPONENT ID, e.g. 'CFF', 'HEM', 'FE'.For nucleic acid polymer entities, use the one-letter code for the base. */
  readonly chem_comp?: Maybe<CoreChemComp>;
  /** Get PDB entry given the PDB id. */
  readonly entry?: Maybe<CoreEntry>;
  /** Get a list of PDB entries given a list of PDB IDs. */
  readonly entries?: Maybe<ReadonlyArray<Maybe<CoreEntry>>>;
  /** Get literature information from PubMed database given the PubMed identifier. */
  readonly pubmed?: Maybe<CorePubmed>;
  /** Get an assembly given the PDB ID and ASSEMBLY ID. Here ASSEMBLY ID is '1', '2', '3', etc. or 'deposited' for deposited coordinates. */
  readonly assembly?: Maybe<CoreAssembly>;
  /** Get UniProt KB entry given the UniProt primary accession. */
  readonly uniprot?: Maybe<CoreUniprot>;
  /** Get a list of PDB non-polymer entity instances (chains), given the list of ENTITY INSTANCE IDs. Here ENTITY INSTANCE ID identifies structural element in the asymmetric unit, e.g. 'A', 'B', etc. */
  readonly nonpolymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntityInstance>>>;
  /** Get a list of PDB polymer entity instances (chains), given the list of ENTITY INSTANCE IDs. Here ENTITY INSTANCE ID identifies structural element in the asymmetric unit, e.g. 'A', 'B', etc. */
  readonly polymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CorePolymerEntityInstance>>>;
  /** Get a PDB non-polymer entity, given the PDB ID and ENTITY ID. Here ENTITY ID is a '1', '2', '3', etc. */
  readonly nonpolymer_entity?: Maybe<CoreNonpolymerEntity>;
};


/** Query root */
export type QueryPolymer_Entity_InstanceArgs = {
  asym_id: Scalars['String'];
  entry_id: Scalars['String'];
};


/** Query root */
export type QueryAssembliesArgs = {
  assembly_ids: ReadonlyArray<Maybe<Scalars['String']>>;
};


/** Query root */
export type QueryNonpolymer_EntitiesArgs = {
  entity_ids: ReadonlyArray<Scalars['String']>;
};


/** Query root */
export type QueryNonpolymer_Entity_InstanceArgs = {
  asym_id: Scalars['String'];
  entry_id: Scalars['String'];
};


/** Query root */
export type QueryPolymer_EntitiesArgs = {
  entity_ids: ReadonlyArray<Scalars['String']>;
};


/** Query root */
export type QueryPolymer_EntityArgs = {
  entity_id: Scalars['String'];
  entry_id: Scalars['String'];
};


/** Query root */
export type QueryChem_CompArgs = {
  comp_id: Scalars['String'];
};


/** Query root */
export type QueryEntryArgs = {
  entry_id: Scalars['String'];
};


/** Query root */
export type QueryEntriesArgs = {
  entry_ids: ReadonlyArray<Scalars['String']>;
};


/** Query root */
export type QueryPubmedArgs = {
  pubmed_id: Scalars['Int'];
};


/** Query root */
export type QueryAssemblyArgs = {
  assembly_id: Scalars['String'];
  entry_id: Scalars['String'];
};


/** Query root */
export type QueryUniprotArgs = {
  uniprot_id: Scalars['String'];
};


/** Query root */
export type QueryNonpolymer_Entity_InstancesArgs = {
  instance_ids: ReadonlyArray<Maybe<Scalars['String']>>;
};


/** Query root */
export type QueryPolymer_Entity_InstancesArgs = {
  instance_ids: ReadonlyArray<Maybe<Scalars['String']>>;
};


/** Query root */
export type QueryNonpolymer_EntityArgs = {
  entity_id: Scalars['String'];
  entry_id: Scalars['String'];
};

export type EntityPoly = {
  /**
   * A flag to indicate whether the polymer contains at least
   *  one monomer-to-monomer link different from that implied by
   *  _entity_poly.type.
   * 
   * Allowable values:
   * n, no, y, yes
   */
  readonly nstd_linkage?: Maybe<Scalars['String']>;
  /**
   * A flag to indicate whether the polymer contains at least
   *  one monomer that is not considered standard.
   * 
   * Allowable values:
   * n, no, y, yes
   */
  readonly nstd_monomer?: Maybe<Scalars['String']>;
  /**
   * Chemical sequence expressed as string of one-letter
   *  amino acid codes. Modifications and non-standard
   *  amino acids are coded as X.
   * 
   * Examples:
   * HHHH(MSE)AKQRSG or AUCGGAAU, A  for alanine or adenine
   * B  for ambiguous asparagine/aspartic-acid
   * R  for arginine
   * N  for asparagine
   * D  for aspartic-acid
   * C  for cysteine or cystine or cytosine
   * Q  for glutamine
   * E  for glutamic-acid
   * Z  for ambiguous glutamine/glutamic acid
   * G  for glycine or guanine
   * H  for histidine
   * I  for isoleucine
   * L  for leucine
   * K  for lysine
   * M  for methionine
   * F  for phenylalanine
   * P  for proline
   * S  for serine
   * T  for threonine or thymine
   * W  for tryptophan
   * Y  for tyrosine
   * V  for valine
   * U  for uracil
   * O  for water
   * X  for other
   */
  readonly pdbx_seq_one_letter_code?: Maybe<Scalars['String']>;
  /**
   * Cannonical chemical sequence expressed as string of
   *                one-letter amino acid codes. Modifications are coded
   *                as the parent amino acid where possible.
   * 
   * A  for alanine or adenine
   * B  for ambiguous asparagine/aspartic-acid
   * R  for arginine
   * N  for asparagine
   * D  for aspartic-acid
   * C  for cysteine or cystine or cytosine
   * Q  for glutamine
   * E  for glutamic-acid
   * Z  for ambiguous glutamine/glutamic acid
   * G  for glycine or guanine
   * H  for histidine
   * I  for isoleucine
   * L  for leucine
   * K  for lysine
   * M  for methionine
   * F  for phenylalanine
   * P  for proline
   * S  for serine
   * T  for threonine or thymine
   * W  for tryptophan
   * Y  for tyrosine
   * V  for valine
   * U  for uracil
   * 
   * Examples:
   * MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGAAFNVEFD
   */
  readonly pdbx_seq_one_letter_code_can?: Maybe<Scalars['String']>;
  /**
   * The PDB strand/chain id(s) corresponding to this polymer entity.
   * 
   * Examples:
   * A,B, A, B, A,B,C
   */
  readonly pdbx_strand_id?: Maybe<Scalars['String']>;
  /**
   * For Structural Genomics entries, the sequence's target identifier registered at the TargetTrack database.
   * 
   * Examples:
   * JCSG-11211, 356560
   */
  readonly pdbx_target_identifier?: Maybe<Scalars['String']>;
  /**
   * Number of regions in the sample sequence identified as expression tags, linkers, or
   *  cloning artifacts.
   */
  readonly rcsb_artifact_monomer_count?: Maybe<Scalars['Int']>;
  /** Number of monomer conflicts relative to the reference sequence. */
  readonly rcsb_conflict_count?: Maybe<Scalars['Int']>;
  /** Number of monomer deletions relative to the reference sequence. */
  readonly rcsb_deletion_count?: Maybe<Scalars['Int']>;
  /**
   * A coarse-grained polymer entity type.
   * 
   * Allowable values:
   * DNA, NA-hybrid, Other, Protein, RNA
   */
  readonly rcsb_entity_polymer_type?: Maybe<Scalars['String']>;
  /** Number of monomer insertions relative to the reference sequence. */
  readonly rcsb_insertion_count?: Maybe<Scalars['Int']>;
  /** Number of engineered mutations engineered in the sample sequence. */
  readonly rcsb_mutation_count?: Maybe<Scalars['Int']>;
  /** Number of non-standard monomers in the sample sequence. */
  readonly rcsb_non_std_monomer_count?: Maybe<Scalars['Int']>;
  /** Unique list of non-standard monomer chemical component identifiers in the sample sequence. */
  readonly rcsb_non_std_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** For polymer BIRD molecules the BIRD identifier for the entity. */
  readonly rcsb_prd_id?: Maybe<Scalars['String']>;
  /** The monomer length of the sample sequence. */
  readonly rcsb_sample_sequence_length?: Maybe<Scalars['Int']>;
  /**
   * The type of the polymer.
   * 
   * Allowable values:
   * cyclic-pseudo-peptide, other, peptide nucleic acid, polydeoxyribonucleotide, polydeoxyribonucleotide/polyribonucleotide hybrid, polypeptide(D), polypeptide(L), polyribonucleotide
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type CorePolymerEntity = {
  /** Get all unique monomers described in this molecular entity. */
  readonly chem_comp_monomers?: Maybe<ReadonlyArray<Maybe<CoreChemComp>>>;
  /** Get all unique non-standard monomers described in this molecular entity. */
  readonly chem_comp_nstd_monomers?: Maybe<ReadonlyArray<Maybe<CoreChemComp>>>;
  readonly entity_poly?: Maybe<EntityPoly>;
  readonly entity_src_gen?: Maybe<ReadonlyArray<Maybe<EntitySrcGen>>>;
  readonly entity_src_nat?: Maybe<ReadonlyArray<Maybe<EntitySrcNat>>>;
  /** Get PDB entry that contains this molecular entity. */
  readonly entry?: Maybe<CoreEntry>;
  readonly pdbx_entity_src_syn?: Maybe<ReadonlyArray<Maybe<PdbxEntitySrcSyn>>>;
  /** Get all unique Pfam annotations associated with this molecular entity. */
  readonly pfams?: Maybe<ReadonlyArray<Maybe<CorePfam>>>;
  /** Get all unique polymer instances (chains) for this molecular entity. */
  readonly polymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CorePolymerEntityInstance>>>;
  /** Get a BIRD chemical components described in this molecular entity. */
  readonly prd?: Maybe<CoreChemComp>;
  /** Indicates intrinsic flexibility of protein structures determined from structural variations between different depositions and chains in asymmetric units of the same protein in PDB (95% sequence identity). */
  readonly rcsb_cluster_flexibility?: Maybe<RcsbClusterFlexibility>;
  readonly rcsb_cluster_membership?: Maybe<ReadonlyArray<Maybe<RcsbClusterMembership>>>;
  readonly rcsb_entity_host_organism?: Maybe<ReadonlyArray<Maybe<RcsbEntityHostOrganism>>>;
  readonly rcsb_entity_source_organism?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganism>>>;
  readonly rcsb_genomic_lineage?: Maybe<ReadonlyArray<Maybe<RcsbGenomicLineage>>>;
  /**
   * A unique identifier for each object in this entity container formed by
   *  an underscore separated concatenation of entry and entity identifiers.
   */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>;
  /** Members of the membrane protein classification lineage. */
  readonly rcsb_membrane_lineage?: Maybe<ReadonlyArray<Maybe<RcsbMembraneLineage>>>;
  /**
   * Mpstruc keyword denotes original annotation, Homology keyword denotes annotation inferred by homology.
   * 
   * Allowable values:
   * Mpstruc, Homology
   */
  readonly rcsb_membrane_lineage_provenance_code?: Maybe<Scalars['String']>;
  readonly rcsb_polymer_entity?: Maybe<RcsbPolymerEntity>;
  readonly rcsb_polymer_entity_align?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAlign>>>;
  readonly rcsb_polymer_entity_annotation?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAnnotation>>>;
  readonly rcsb_polymer_entity_container_identifiers: RcsbPolymerEntityContainerIdentifiers;
  readonly rcsb_polymer_entity_feature?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeature>>>;
  readonly rcsb_polymer_entity_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeatureSummary>>>;
  readonly rcsb_polymer_entity_keywords?: Maybe<RcsbPolymerEntityKeywords>;
  readonly rcsb_polymer_entity_name_com?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityNameCom>>>;
  readonly rcsb_polymer_entity_name_sys?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityNameSys>>>;
  /** Get all unique UniProt KB annotations associated with this molecular entity. */
  readonly uniprots?: Maybe<ReadonlyArray<Maybe<CoreUniprot>>>;
};

export type PdbxSolnScatterModel = {
  /**
   * A description of the conformer selection criteria
   *  used.
   * 
   * Examples:
   * The modelled scattering curves were assessed by calculation of the 
   *    RG, RSX-1 and RXS-2 values in the same Q ranges 
   *    used in the experimental Guinier fits. models were 
   *    then ranked using a goodness-of-fit R-factor 
   *    defined by analogy with protein crystallography 
   *    and based on the experimental curves in the Q range 
   *    extending to 1.4 nm-1.
   */
  readonly conformer_selection_criteria?: Maybe<Scalars['String']>;
  /**
   * A description of any additional details concerning the experiment.
   * 
   * Examples:
   * Homology models were built for
   *     the 17 SCR domains and energy minimisations were 
   *     performed to improve the connectivity in the fh model.
   *     triantennary complex-type carbohydrate structures
   *     (MAN3GLCNAC6GAL3FUC3NEUNAC1) were added to each of the
   *     N-linked glycosylation sites. a library of linker peptide
   *     conformations was used in domain modelling constrained
   *     by the solution scattering fits. modelling with the
   *     scattering data was also carried out by rotational 
   *     search methods. the x-ray and neutron scattering curve 
   *     I(Q) was calculated assuming a uniform scattering density 
   *     for the spheres using the debye equation as adapted to 
   *     spheres. x-ray curves were calculated from the hydrated 
   *     sphere models without corrections for wavelength spread or
   *     beam divergence, while these corrections were applied for 
   *     the neutron curves but now using unhydrated models.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * A list of the entries used to fit the model
   *  to the scattering data
   * 
   * Examples:
   * PDB CODE 1HFI, 1HCC, 1HFH, 1VCC
   */
  readonly entry_fitting_list?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_soln_scatter_model.id must
   *  uniquely identify the sample in the category PDBX_SOLN_SCATTER_MODEL
   */
  readonly id: Scalars['String'];
  /**
   * A description of the methods used in the modelling
   * 
   * Examples:
   * Constrained scattering fitting of homology models
   */
  readonly method?: Maybe<Scalars['String']>;
  /** The number of model conformers calculated. */
  readonly num_conformers_calculated?: Maybe<Scalars['Int']>;
  /** The number of model conformers submitted in the entry */
  readonly num_conformers_submitted?: Maybe<Scalars['Int']>;
  /** The index of the representative conformer among the submitted conformers for the entry */
  readonly representative_conformer?: Maybe<Scalars['Int']>;
  /** This data item is a pointer to  _pdbx_soln_scatter.id in the  PDBX_SOLN_SCATTER category. */
  readonly scatter_id: Scalars['String'];
  /**
   * A list of the software authors
   * 
   * Examples:
   * MSI
   */
  readonly software_author_list?: Maybe<Scalars['String']>;
  /**
   * A list of the software used in the modeeling
   * 
   * Examples:
   * INSIGHT II, HOMOLOGY, DISCOVERY, BIOPOLYMER, DELPHI
   */
  readonly software_list?: Maybe<Scalars['String']>;
};

export type RcsbUniprotFeature = {
  /** Identifies the version of the feature assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the feature. */
  readonly description?: Maybe<Scalars['String']>;
  /** An identifier for the feature. */
  readonly feature_id?: Maybe<Scalars['String']>;
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeatureFeaturePositions>>>;
  /** A name for the feature. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the feature.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /** Code residue coordinate system for the assigned feature. */
  readonly reference_scheme?: Maybe<Scalars['String']>;
  /**
   * A type or category of the feature.
   * 
   * Allowable values:
   * ACTIVE_SITE, BINDING_SITE, CALCIUM_BINDING_REGION, CHAIN, COMPOSITIONALLY_BIASED_REGION, CROSS_LINK, DNA_BINDING_REGION, DOMAIN, GLYCOSYLATION_SITE, INITIATOR_METHIONINE, LIPID_MOIETY_BINDING_REGION, METAL_ION_BINDING_SITE, MODIFIED_RESIDUE, MUTAGENESIS_SITE, NON_CONSECUTIVE_RESIDUES, NON_TERMINAL_RESIDUE, NUCLEOTIDE_PHOSPHATE_BINDING_REGION, PEPTIDE, PROPEPTIDE, REGION_OF_INTEREST, REPEAT, NON_STANDARD_AMINO_ACID, SEQUENCE_CONFLICT, SEQUENCE_VARIANT, SHORT_SEQUENCE_MOTIF, SIGNAL_PEPTIDE, SITE, SPLICE_VARIANT, TOPOLOGICAL_DOMAIN, TRANSIT_PEPTIDE, TRANSMEMBRANE_REGION, UNSURE_RESIDUE, ZINC_FINGER_REGION, INTRAMEMBRANE_REGION
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type Em3dFittingList = {
  /**
   * The value of _em_3d_fitting_list.3d_fitting_id is a pointer
   *  to  _em_3d_fitting.id in the 3d_fitting category
   */
  readonly _3d_fitting_id: Scalars['String'];
  /** Details about the model used in fitting. */
  readonly details?: Maybe<Scalars['String']>;
  /** This data item is a unique identifier. */
  readonly id: Scalars['String'];
  /**
   * The ID of the biopolymer chain used for fitting, e.g., A.  Please note that
   * only one chain can be specified per instance.  If all chains of a particular
   * structure have been used for fitting, this field can be left blank.
   */
  readonly pdb_chain_id?: Maybe<Scalars['String']>;
  /** The molecular entities represented in this fitting description. */
  readonly pdb_chain_residue_range?: Maybe<Scalars['String']>;
  /**
   * The PDB code for the entry used in fitting.
   * 
   * Examples:
   * PDB entry 1EHZ
   */
  readonly pdb_entry_id?: Maybe<Scalars['String']>;
};

export type EmDiffraction = {
  /**
   * TODO
   * 
   * Examples:
   * 800
   */
  readonly camera_length?: Maybe<Scalars['Float']>;
  /** Primary key */
  readonly id: Scalars['String'];
  /** Foreign key to the EM_IMAGING category */
  readonly imaging_id?: Maybe<Scalars['String']>;
  /**
   * Comma-separated list of tilt angles (in degrees) used in the electron diffraction experiment.
   * 
   * Examples:
   * 20,40,50,55
   */
  readonly tilt_angle_list?: Maybe<Scalars['String']>;
};

export type RcsbChemCompRelated = {
  /**
   * The value of _rcsb_chem_comp_related.comp_id is a reference to
   *  a chemical component definition.
   */
  readonly comp_id: Scalars['String'];
  /**
   * The value of _rcsb_chem_comp_related.ordinal distinguishes
   *  related examples for each chemical component.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The method used to establish the resource correspondence.
   * 
   * Allowable values:
   * assigned by DrugBank resource, assigned by PDB, matching InChIKey in DrugBank, matching InChIKey-prefix in DrugBank, matching by RESID resource
   */
  readonly related_mapping_method?: Maybe<Scalars['String']>;
  /**
   * The resource identifier code for the related chemical reference.
   * 
   * Examples:
   * 124832
   */
  readonly resource_accession_code?: Maybe<Scalars['String']>;
  /**
   * The resource name for the related chemical reference.
   * 
   * Allowable values:
   * CAS, CCDC/CSD, ChEBI, ChEMBL, DrugBank, PubChem, RESID
   */
  readonly resource_name?: Maybe<Scalars['String']>;
};

export type PdbxNmrSpectrometer = {
  /** A text description of the NMR spectrometer. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Select the field strength for protons in MHz.
   * 
   * Examples:
   * 360, 400, 500, 600, 750, 800, 850, 900, 950, 1000
   */
  readonly field_strength?: Maybe<Scalars['Float']>;
  /**
   * The name of the manufacturer of the spectrometer.
   * 
   * Examples:
   * Varian, Bruker, JEOL, GE
   */
  readonly manufacturer?: Maybe<Scalars['String']>;
  /**
   * The model of the NMR spectrometer.
   * 
   * Examples:
   * AVANCE, AVANCE II, AVANCE III, AVANCE III HD, WH, WM, AM, AMX, DMX, DRX, MSL, OMEGA, OMEGA PSG, GX, GSX, A, AL, EC, EX, LA, ECP, VXRS, UNITY, UNITYPLUS, INOVA
   */
  readonly model?: Maybe<Scalars['String']>;
  /**
   * Assign a numerical ID to each instrument.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly spectrometer_id: Scalars['String'];
  /**
   * Select the instrument manufacturer(s) and the model(s) of the NMR(s)
   * used for this work.
   * 
   * Examples:
   * Bruker WH, Bruker WM, Bruker AM, Bruker AMX, Bruker DMX, Bruker DRX, Bruker MSL, Bruker AVANCE, GE Omega, GE Omega PSG, JEOL GX, JEOL GSX, JEOL A, JEOL AL, JEOL EC, JEOL EX, JEOL LA, JEOL ECP, Varian VXRS, Varian UNITY, Varian UNITYplus, Varian INOVA, other
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxAuditRevisionHistory = {
  /**
   * The type of file that the pdbx_audit_revision_history record refers to.
   * 
   * Allowable values:
   * Chemical component, NMR restraints, NMR shifts, Structure factors, Structure model
   */
  readonly data_content_type: Scalars['String'];
  /**
   * The major version number of deposition release.
   * 
   * Examples:
   * 1
   */
  readonly major_revision?: Maybe<Scalars['Int']>;
  /**
   * The minor version number of deposition release.
   * 
   * Examples:
   * 1
   */
  readonly minor_revision?: Maybe<Scalars['Int']>;
  /**
   * A unique identifier for the pdbx_audit_revision_history record.
   * 
   * Examples:
   * 1
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The release date of the revision
   * 
   * Examples:
   * 2017-03-08
   */
  readonly revision_date?: Maybe<Scalars['Date']>;
};

export type DiffrnSource = {
  /** A description of special aspects of the radiation source used. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   */
  readonly diffrn_id: Scalars['String'];
  /** Synchrotron beamline. */
  readonly pdbx_synchrotron_beamline?: Maybe<Scalars['String']>;
  /** Synchrotron site. */
  readonly pdbx_synchrotron_site?: Maybe<Scalars['String']>;
  /** Wavelength of radiation. */
  readonly pdbx_wavelength?: Maybe<Scalars['String']>;
  /**
   * Comma separated list of wavelengths or wavelength range.
   * 
   * Examples:
   * 0.987 or 0.987, 0.988, 1.0 or 0.99-1.5
   */
  readonly pdbx_wavelength_list?: Maybe<Scalars['String']>;
  /**
   * The general class of the radiation source.
   * 
   * Examples:
   * sealed X-ray tube, nuclear reactor, spallation source, electron microscope, rotating-anode X-ray tube, synchrotron
   */
  readonly source?: Maybe<Scalars['String']>;
  /**
   * The make, model or name of the source of radiation.
   * 
   * Examples:
   * NSLS beamline X8C, Rigaku RU200
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxStructAssemblyGen = {
  /**
   * This data item is a pointer to _pdbx_struct_assembly.id in the 
   *  PDBX_STRUCT_ASSEMBLY category.
   */
  readonly assembly_id?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _struct_asym.id in 
   *  the STRUCT_ASYM category.
   * 
   *  This item may be expressed as a comma separated list of identifiers.
   */
  readonly asym_id_list?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Identifies the operation of collection of operations 
   *  from category PDBX_STRUCT_OPER_LIST.  
   * 
   *  Operation expressions may have the forms:
   * 
   *   (1)        the single operation 1
   *   (1,2,5)    the operations 1, 2, 5
   *   (1-4)      the operations 1,2,3 and 4
   *   (1,2)(3,4) the combinations of operations
   *              3 and 4 followed by 1 and 2 (i.e.
   *              the cartesian product of parenthetical
   *              groups applied from right to left)
   * 
   * Examples:
   * (1), (1,2,5), (1-60), (1-60)(61)
   */
  readonly oper_expression?: Maybe<Scalars['String']>;
  /**
   * This data item is an ordinal index for the
   *  PDBX_STRUCT_ASSEMBLY category.
   */
  readonly ordinal: Scalars['Int'];
};

export type PdbxNmrExptl = {
  /**
   * The number to identify the set of sample conditions.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly conditions_id: Scalars['String'];
  /**
   * A numerical ID for each experiment.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly experiment_id: Scalars['String'];
  /**
   * Physical state of the sample either anisotropic or isotropic.
   * 
   * Allowable values:
   * anisotropic, isotropic
   */
  readonly sample_state?: Maybe<Scalars['String']>;
  /**
   * The solution_id from the Experimental Sample to identify the sample
   *  that these conditions refer to. 
   * 
   *  [Remember to save the entries here before returning to the 
   *   Experimental Sample form]
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly solution_id: Scalars['String'];
  /** Pointer to '_pdbx_nmr_spectrometer.spectrometer_id' */
  readonly spectrometer_id?: Maybe<Scalars['Int']>;
  /**
   * The type of NMR experiment.
   * 
   * Examples:
   * 2D NOESY, 3D_15N-separated_NOESY, 3D_13C-separated_NOESY, 4D_13C-separated_NOESY, 4D_13C/15N-separated_NOESY, 3D_15N-separated_ROESY, 3D_13C-separated_ROESY, HNCA-J, HNHA, DQF-COSY, P-COSY, PE-COSY, E-COSY
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbChemCompInfo = {
  /** Chemical component total atom count */
  readonly atom_count?: Maybe<Scalars['Int']>;
  /** Chemical component chiral atom count */
  readonly atom_count_chiral?: Maybe<Scalars['Int']>;
  /** Chemical component heavy atom count */
  readonly atom_count_heavy?: Maybe<Scalars['Int']>;
  /** Chemical component total bond count */
  readonly bond_count?: Maybe<Scalars['Int']>;
  /** Chemical component aromatic bond count */
  readonly bond_count_aromatic?: Maybe<Scalars['Int']>;
  /** The chemical component identifier. */
  readonly comp_id: Scalars['String'];
  /** The initial date the chemical definition was released in the PDB repository. */
  readonly initial_release_date?: Maybe<Scalars['Date']>;
  /**
   * The release status of the chemical definition.
   * 
   * Allowable values:
   * DEL, HOLD, HPUB, OBS, REF_ONLY, REL
   */
  readonly release_status?: Maybe<Scalars['String']>;
  /** The date of last revision of the chemical definition. */
  readonly revision_date?: Maybe<Scalars['Date']>;
};

export type RcsbPolymerEntityInstanceContainerIdentifiers = {
  /** Instance identifier for this container. */
  readonly asym_id: Scalars['String'];
  /** Author instance identifier for this container. */
  readonly auth_asym_id?: Maybe<Scalars['String']>;
  /** Entity identifier for the container. */
  readonly entity_id?: Maybe<Scalars['String']>;
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /**
   * A unique identifier for each object in this entity instance container formed by
   *  an 'dot' (.) separated concatenation of entry and entity instance identifiers.
   */
  readonly rcsb_id?: Maybe<Scalars['String']>;
};

export type EmCtfCorrection = {
  /**
   * Any additional details about CTF correction
   * 
   * Examples:
   * CTF amplitude correction was performed following 3D reconstruction
   */
  readonly details?: Maybe<Scalars['String']>;
  /** Foreign key to the EM_IMAGE_PROCESSING category */
  readonly em_image_processing_id?: Maybe<Scalars['String']>;
  /** Primary key */
  readonly id: Scalars['String'];
  /** Type of CTF correction applied */
  readonly type?: Maybe<Scalars['String']>;
};

export type CoreEntityAlignmentsScores = {
  readonly query_coverage: Scalars['Int'];
  readonly query_length: Scalars['Int'];
  readonly target_coverage: Scalars['Int'];
  readonly target_length: Scalars['Int'];
};

export type RcsbChemCompContainerIdentifiers = {
  /**
   * The Anatomical Therapeutic Chemical (ATC) Classification System identifiers corresponding
   *  to the chemical component.
   */
  readonly atc_codes?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The chemical component identifier. */
  readonly comp_id: Scalars['String'];
  /** The DrugBank identifier corresponding to the chemical component. */
  readonly drugbank_id?: Maybe<Scalars['String']>;
  /** The BIRD definition identifier. */
  readonly prd_id?: Maybe<Scalars['String']>;
  /** A unique identifier for the chemical definition in this container. */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  /** The list of subcomponents contained in this component. */
  readonly subcomponent_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type PdbxAuditSupport = {
  /** The country/region providing the funding support for the entry. */
  readonly country?: Maybe<Scalars['String']>;
  /**
   * The name of the organization providing funding support for the 
   *  entry.
   * 
   * Examples:
   * National Institutes of Health, Welcome Trust, National Institutes of Health/National Institute of General Medical Sciences
   */
  readonly funding_organization?: Maybe<Scalars['String']>;
  /** The grant number associated with this source of support. */
  readonly grant_number?: Maybe<Scalars['String']>;
  /**
   * A unique sequential integer identifier for each source of support for this entry.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly ordinal: Scalars['Int'];
};

export type Entry = {
  /**
   * The value of _entry.id identifies the data block.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   */
  readonly id: Scalars['String'];
};

export type EmSpecimen = {
  /**
   * The concentration (in milligrams per milliliter, mg/ml)
   *  of the complex in the sample.
   * 
   * Examples:
   * 1.35
   */
  readonly concentration?: Maybe<Scalars['Float']>;
  /**
   * A description of any additional details of the specimen preparation.
   * 
   * Examples:
   * This sample was monodisperse., Au was deposited at a 30 degree angle to 15 nm thickness., Colloidal gold particles were deposited by dipping into dilute solution., The specimen was frozen at high pressure using the bal-tec hpm 010 instrument., The embedded sample was sectioned at 100 K to 50 nm final thickness.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * 'YES' indicates that the specimen has been embedded.
   * 
   * Allowable values:
   * NO, YES
   */
  readonly embedding_applied?: Maybe<Scalars['String']>;
  /** Pointer to _em_experiment.id. */
  readonly experiment_id: Scalars['String'];
  /**
   * The item  _em_specimen.id uniquely identifies a specimen along with
   *  its preparation methods.
   */
  readonly id: Scalars['String'];
  /**
   * 'YES' indicates that the specimen has been shadowed.
   * 
   * Allowable values:
   * NO, YES
   */
  readonly shadowing_applied?: Maybe<Scalars['String']>;
  /**
   * 'YES' indicates that the specimen has been stained.
   * 
   * Allowable values:
   * NO, YES
   */
  readonly staining_applied?: Maybe<Scalars['String']>;
  /**
   * 'YES' indicates that the specimen was vitrified by cryopreservation.
   * 
   * Allowable values:
   * NO, YES
   */
  readonly vitrification_applied?: Maybe<Scalars['String']>;
};

export type PdbxChemCompDescriptor = {
  /**
   * This data item is a pointer to _chem_comp.id in the CHEM_COMP
   *  category.
   */
  readonly comp_id: Scalars['String'];
  /**
   * This data item contains the descriptor value for this 
   *  component.
   */
  readonly descriptor?: Maybe<Scalars['String']>;
  /**
   * This data item contains the name of the program
   *  or library used to compute the descriptor.
   * 
   * Examples:
   * OPENEYE, CACTVS, DAYLIGHT, OTHER
   */
  readonly program: Scalars['String'];
  /**
   * This data item contains the version of the program
   *  or library used to compute the descriptor.
   */
  readonly program_version: Scalars['String'];
  /**
   * This data item contains the descriptor type.
   * 
   * Allowable values:
   * InChI, InChIKey, InChI_CHARGE, InChI_FIXEDH, InChI_ISOTOPE, InChI_MAIN, InChI_MAIN_CONNECT, InChI_MAIN_FORMULA, InChI_MAIN_HATOM, InChI_RECONNECT, InChI_STEREO, SMILES, SMILES_CANNONICAL, SMILES_CANONICAL
   */
  readonly type: Scalars['String'];
};

export type RcsbChemCompDescriptor = {
  /**
   * Standard IUPAC International Chemical Identifier (InChI) descriptor for the chemical component.
   * 
   *    InChI, the IUPAC International Chemical Identifier,
   *    by Stephen R Heller, Alan McNaught, Igor Pletnev, Stephen Stein and Dmitrii Tchekhovskoi,
   *    Journal of Cheminformatics, 2015, 7:23;
   */
  readonly InChI?: Maybe<Scalars['String']>;
  /**
   * Standard IUPAC International Chemical Identifier (InChI) descriptor key
   *  for the chemical component
   * 
   *  InChI, the IUPAC International Chemical Identifier,
   *  by Stephen R Heller, Alan McNaught, Igor Pletnev, Stephen Stein and Dmitrii Tchekhovskoi,
   *  Journal of Cheminformatics, 2015, 7:23
   */
  readonly InChIKey?: Maybe<Scalars['String']>;
  /**
   * Simplified molecular-input line-entry system (SMILES) descriptor for the chemical component.
   * 
   *    Weininger D (February 1988). "SMILES, a chemical language and information system. 1.
   *    Introduction to methodology and encoding rules". Journal of Chemical Information and Modeling. 28 (1): 31-6.
   * 
   *    Weininger D, Weininger A, Weininger JL (May 1989).
   *    "SMILES. 2. Algorithm for generation of unique SMILES notation",
   *    Journal of Chemical Information and Modeling. 29 (2): 97-101.
   */
  readonly SMILES?: Maybe<Scalars['String']>;
  /**
   * Simplified molecular-input line-entry system (SMILES) descriptor for the chemical
   *  component including stereochemical features.
   * 
   *  Weininger D (February 1988). "SMILES, a chemical language and information system. 1.
   *  Introduction to methodology and encoding rules".
   *  Journal of Chemical Information and Modeling. 28 (1): 31-6.
   * 
   *  Weininger D, Weininger A, Weininger JL (May 1989).
   *  "SMILES. 2. Algorithm for generation of unique SMILES notation".
   *  Journal of Chemical Information and Modeling. 29 (2): 97-101.
   */
  readonly SMILES_stereo?: Maybe<Scalars['String']>;
  /** The chemical component identifier. */
  readonly comp_id: Scalars['String'];
};

export type PdbxNmrEnsemble = {
  /**
   * The average number of constraint violations on a per residue basis for
   *  the ensemble.
   * 
   * Examples:
   * 0.25
   */
  readonly average_constraint_violations_per_residue?: Maybe<Scalars['Int']>;
  /**
   * The average number of constraints per residue for the ensemble
   * 
   * Examples:
   * 30.2
   */
  readonly average_constraints_per_residue?: Maybe<Scalars['Int']>;
  /**
   * The average distance restraint violation for the ensemble.
   * 
   * Examples:
   * 0.11
   */
  readonly average_distance_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * The average torsion angle constraint violation for the ensemble.
   * 
   * Examples:
   * 2.4
   */
  readonly average_torsion_angle_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * By highlighting the appropriate choice(s), describe how the submitted 
   * conformer (models) were selected.
   * 
   * Examples:
   * structures with the lowest energy, structures with the least restraint violations, structures with acceptable covalent geometry, structures with favorable non-bond energy, target function, back calculated data agree with experimental NOESY spectrum, all calculated structures submitted, The submitted conformer models are the 25 structures with the lowest 
   *     energy., The submitted conformer models are those with the fewest number of 
   *     constraint violations.
   */
  readonly conformer_selection_criteria?: Maybe<Scalars['String']>;
  /**
   * The total number of conformer (models) that were calculated in the final round.
   * 
   * Examples:
   * 40
   */
  readonly conformers_calculated_total_number?: Maybe<Scalars['Int']>;
  /**
   * The number of conformer (models) that are submitted for the ensemble.
   * 
   * Examples:
   * 20
   */
  readonly conformers_submitted_total_number?: Maybe<Scalars['Int']>;
  /**
   * Describe the method used to calculate the distance constraint violation statistics,
   *  i.e. are they calculated over all the distance constraints or calculated for
   *  violations only?
   * 
   * Examples:
   * Statistics were calculated over all of the distance constraints., Statistics were calculated for violations only
   */
  readonly distance_constraint_violation_method?: Maybe<Scalars['String']>;
  /**
   * The maximum distance constraint violation for the ensemble.
   * 
   * Examples:
   * 0.4
   */
  readonly maximum_distance_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * The maximum lower distance constraint violation for the ensemble.
   * 
   * Examples:
   * 0.3
   */
  readonly maximum_lower_distance_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * The maximum torsion angle constraint violation for the ensemble.
   * 
   * Examples:
   * 4
   */
  readonly maximum_torsion_angle_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * The maximum upper distance constraint violation for the ensemble.
   * 
   * Examples:
   * 0.4
   */
  readonly maximum_upper_distance_constraint_violation?: Maybe<Scalars['Float']>;
  /**
   * The number of the conformer identified as most representative.
   * 
   * Examples:
   * 20
   */
  readonly representative_conformer?: Maybe<Scalars['Int']>;
  /**
   * This item describes the method used to calculate the torsion angle constraint violation statistics.
   * i.e. are the entered values based on all torsion angle or calculated for violations only?
   * 
   * Examples:
   * Statistics were calculated over all the torsion angle constraints., Statistics were calculated for torsion angle constraints violations only.
   */
  readonly torsion_angle_constraint_violation_method?: Maybe<Scalars['String']>;
};

export type RefineLsRestr = {
  /**
   * For the given parameter type, the root-mean-square deviation
   *  between the ideal values used as restraints in the least-squares
   *  refinement and the values obtained by refinement. For instance,
   *  bond distances may deviate by 0.018 \%A (r.m.s.) from ideal
   *  values in the current model.
   */
  readonly dev_ideal?: Maybe<Scalars['Float']>;
  /**
   * For the given parameter type, the target root-mean-square
   *  deviation between the ideal values used as restraints in the
   *  least-squares refinement and the values obtained by refinement.
   */
  readonly dev_ideal_target?: Maybe<Scalars['Float']>;
  /**
   * The number of parameters of this type subjected to restraint in
   *  least-squares refinement.
   */
  readonly number?: Maybe<Scalars['Int']>;
  /**
   * This data item uniquely identifies a refinement within an entry.
   *  _refine_ls_restr.pdbx_refine_id can be used to distinguish the results 
   *  of joint refinements.
   */
  readonly pdbx_refine_id: Scalars['String'];
  /**
   * The functional form of the restraint function used in the least-squares
   *  refinement.
   * 
   * Examples:
   * SINUSOIDAL, HARMONIC, SEMIHARMONIC
   */
  readonly pdbx_restraint_function?: Maybe<Scalars['String']>;
  /**
   * The type of the parameter being restrained.
   *  Explicit sets of data values are provided for the programs
   *  PROTIN/PROLSQ (beginning with p_) and RESTRAIN (beginning with
   *  RESTRAIN_). As computer programs change, these data values
   *  are given as examples, not as an enumeration list. Computer
   *  programs that convert a data block to a refinement table will
   *  expect the exact form of the data values given here to be used.
   * 
   * Examples:
   * p_bond_d, p_angle_d, p_planar_d, p_xhbond_d, p_xhangle_d, p_hydrog_d, p_special_d, p_planar, p_chiral, p_singtor_nbd, p_multtor_nbd, p_xyhbond_nbd, p_xhyhbond_nbd, p_special_tor, p_planar_tor, p_staggered_tor, p_orthonormal_tor, p_mcbond_it, p_mcangle_it, p_scbond_it, p_scangle_it, p_xhbond_it, p_xhangle_it, p_special_it, RESTRAIN_Distances < 2.12, RESTRAIN_Distances 2.12 < D < 2.625, RESTRAIN_Distances > 2.625, RESTRAIN_Peptide Planes, RESTRAIN_Ring and other planes, RESTRAIN_rms diffs for Uiso atoms at dist 1.2-1.4, RESTRAIN_rms diffs for Uiso atoms at dist 1.4-1.6, RESTRAIN_rms diffs for Uiso atoms at dist 1.8-2.0, RESTRAIN_rms diffs for Uiso atoms at dist 2.0-2.2, RESTRAIN_rms diffs for Uiso atoms at dist 2.2-2.4, RESTRAIN_rms diffs for Uiso atoms at dist >2.4
   */
  readonly type: Scalars['String'];
  /**
   * The weighting value applied to this type of restraint in
   *  the least-squares refinement.
   */
  readonly weight?: Maybe<Scalars['Float']>;
};

export type PdbxReferenceEntityList = {
  /** The component number of this entity within the molecule. */
  readonly component_id: Scalars['Int'];
  /** Additional details about this entity. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_list.prd_id is a reference
   *  _pdbx_reference_molecule.prd_id in the PDBX_REFERENCE_MOLECULE category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_list.ref_entity_id is a unique identifier 
   *  the a constituent entity within this reference molecule.
   */
  readonly ref_entity_id: Scalars['String'];
  /**
   * Defines the polymer characteristic of the entity.
   * 
   * Allowable values:
   * branched, non-polymer, polymer, polymer-like
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbUniprotAlignments = {
  readonly core_entity_alignments?: Maybe<ReadonlyArray<Maybe<RcsbUniprotAlignmentsCoreEntityAlignments>>>;
};

export type CoreNonpolymerEntity = {
  /** Get PDB entry that contains this molecular entity. */
  readonly entry?: Maybe<CoreEntry>;
  /** Get a non-polymer chemical components described in this molecular entity. */
  readonly nonpolymer_comp?: Maybe<CoreChemComp>;
  /** Get all unique non-polymer instances (chains) for this molecular entity. */
  readonly nonpolymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntityInstance>>>;
  readonly pdbx_entity_nonpoly?: Maybe<PdbxEntityNonpoly>;
  /** Get a BIRD chemical components described in this molecular entity. */
  readonly prd?: Maybe<CoreChemComp>;
  /**
   * A unique identifier for each object in this entity container formed by
   *  an underscore separated concatenation of entry and entity identifiers.
   */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>;
  readonly rcsb_nonpolymer_entity?: Maybe<RcsbNonpolymerEntity>;
  readonly rcsb_nonpolymer_entity_annotation?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityAnnotation>>>;
  readonly rcsb_nonpolymer_entity_container_identifiers?: Maybe<RcsbNonpolymerEntityContainerIdentifiers>;
  readonly rcsb_nonpolymer_entity_feature?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityFeature>>>;
  readonly rcsb_nonpolymer_entity_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityFeatureSummary>>>;
  readonly rcsb_nonpolymer_entity_keywords?: Maybe<RcsbNonpolymerEntityKeywords>;
  readonly rcsb_nonpolymer_entity_name_com?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityNameCom>>>;
};

export type RcsbPolymerEntity = {
  /** A description of special aspects of the entity. */
  readonly details?: Maybe<Scalars['String']>;
  /** Formula mass (KDa) of the entity. */
  readonly formula_weight?: Maybe<Scalars['Float']>;
  /**
   * A description of the polymer entity.
   * 
   * Examples:
   * Green fluorescent protein, DNA (5'-D(*GP*(CH3)CP*GP*(CH3)CP*GP*C)-3'), PROFLAVINE, PROTEIN (DEOXYRIBONUCLEASE I (E.C.3.1.21.1))
   */
  readonly pdbx_description?: Maybe<Scalars['String']>;
  /**
   * Enzyme Commission (EC) number(s)
   * 
   * Examples:
   * 2.7.7.7
   */
  readonly pdbx_ec?: Maybe<Scalars['String']>;
  /**
   * Polymer entity fragment description(s).
   * 
   * Examples:
   * KLENOW FRAGMENT, REPLICASE OPERATOR HAIRPIN, C-TERMINAL DOMAIN
   */
  readonly pdbx_fragment?: Maybe<Scalars['String']>;
  /**
   * Details about any polymer entity mutation(s).
   * 
   * Examples:
   * Y31H, DEL(298-323)
   */
  readonly pdbx_mutation?: Maybe<Scalars['String']>;
  /**
   * The number of molecules of the entity in the entry.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly pdbx_number_of_molecules?: Maybe<Scalars['Int']>;
  readonly rcsb_ec_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbEcLineage>>>;
  readonly rcsb_enzyme_class_combined?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbEnzymeClassCombined>>>;
  readonly rcsb_macromolecular_names_combined?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbMacromolecularNamesCombined>>>;
  /**
   * A code indicating the entity has multiple biological sources.
   * 
   * Allowable values:
   * N, Y
   */
  readonly rcsb_multiple_source_flag?: Maybe<Scalars['String']>;
  /** The number of biological sources for the polymer entity. */
  readonly rcsb_source_part_count?: Maybe<Scalars['Int']>;
  /**
   * The method by which the sample for the polymer entity was produced.
   *  Entities isolated directly from natural sources (tissues, soil
   *  samples etc.) are expected to have further information in the
   *  ENTITY_SRC_NAT category. Entities isolated from genetically
   *  manipulated sources are expected to have further information in
   *  the ENTITY_SRC_GEN category.
   * 
   * Allowable values:
   * man, nat, syn
   */
  readonly src_method?: Maybe<Scalars['String']>;
};

export type RcsbUniprotFeatureFeaturePositions = {
  /** An identifier for the monomer(s) corresponding to the feature assignment. */
  readonly beg_comp_id?: Maybe<Scalars['String']>;
  /** An identifier for the monomer at which this segment of the feature begins. */
  readonly beg_seq_id: Scalars['Int'];
  /** An identifier for the monomer at which this segment of the feature ends. */
  readonly end_seq_id?: Maybe<Scalars['Int']>;
  /** The value for the feature over this monomer segment. */
  readonly value?: Maybe<Scalars['Float']>;
};

export type CoreNonpolymerEntityInstance = {
  /** Get non-polymer entity for this non-polymer entity instance. */
  readonly nonpolymer_entity?: Maybe<CoreNonpolymerEntity>;
  readonly pdbx_struct_special_symmetry?: Maybe<ReadonlyArray<Maybe<PdbxStructSpecialSymmetry>>>;
  /**
   * A unique identifier for each object in this entity instance container formed by
   *  an 'dot' (.) separated concatenation of entry and entity instance identifiers.
   */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>;
  readonly rcsb_nonpolymer_entity_instance_container_identifiers?: Maybe<RcsbNonpolymerEntityInstanceContainerIdentifiers>;
  readonly rcsb_nonpolymer_instance_annotation?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceAnnotation>>>;
  readonly rcsb_nonpolymer_instance_feature?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeature>>>;
  readonly rcsb_nonpolymer_instance_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeatureSummary>>>;
  readonly rcsb_nonpolymer_struct_conn?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerStructConn>>>;
};

export type RcsbUniprotProteinEc = {
  readonly number?: Maybe<Scalars['String']>;
  /** Historical record of the data attribute. */
  readonly provenance_code?: Maybe<Scalars['String']>;
};

export type RcsbChemCompSynonyms = {
  /** The chemical component to which this synonym applies. */
  readonly comp_id: Scalars['String'];
  /** The synonym of this particular chemical component. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * This data item is an ordinal index for the
   *  RCSB_CHEM_COMP_SYNONYMS category.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The provenance of this synonym.
   * 
   * Allowable values:
   * ACDLabs, Author, ChEBI, ChEMBL, DrugBank, GMML, Lexichem, OpenEye OEToolkits, OpenEye/Lexichem, PDB Reference Data, PDB-CARE, PubChem, RESID
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
};

export type PdbxDepositGroup = {
  /** A description of the contents of entries in the collection. */
  readonly group_description?: Maybe<Scalars['String']>;
  /** A unique identifier for a group of entries deposited as a collection. */
  readonly group_id: Scalars['String'];
  /** A title to describe the group of entries deposited in the collection. */
  readonly group_title?: Maybe<Scalars['String']>;
  /**
   * Text to describe a grouping of entries in multiple collections
   * 
   * Allowable values:
   * changed state, ground state, undefined
   */
  readonly group_type?: Maybe<Scalars['String']>;
};

export type Reflns = {
  /**
   * The value of the overall isotropic displacement parameter
   *  estimated from the slope of the Wilson plot.
   */
  readonly B_iso_Wilson_estimate?: Maybe<Scalars['Float']>;
  /**
   * A description of the method by which a subset of reflections was
   *  selected for exclusion from refinement so as to be used in the
   *  calculation of a 'free' R factor.
   * 
   * Examples:
   * The data set was sorted with l varying most
   *                                   rapidly and h varying least rapidly. Every
   *                                   10th reflection in this sorted list was
   *                                   excluded from refinement and included in the
   *                                   calculation of a 'free' R factor.
   */
  readonly R_free_details?: Maybe<Scalars['String']>;
  /**
   * Residual factor Rmerge for all reflections that satisfy the
   *  resolution limits established by _reflns.d_resolution_high
   *  and _reflns.d_resolution_low.
   * 
   *              sum~i~(sum~j~|F~j~ - <F>|)
   *  Rmerge(F) = --------------------------
   *                   sum~i~(sum~j~<F>)
   * 
   *  F~j~ = the amplitude of the jth observation of reflection i
   *  <F>  = the mean of the amplitudes of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_F_all?: Maybe<Scalars['Float']>;
  /**
   * Residual factor Rmerge for reflections that satisfy the
   *  resolution limits established by _reflns.d_resolution_high
   *  and _reflns.d_resolution_low and the observation limit
   *  established by _reflns.observed_criterion.
   * 
   *              sum~i~(sum~j~|F~j~ - <F>|)
   *  Rmerge(F) = --------------------------
   *                   sum~i~(sum~j~<F>)
   * 
   *  F~j~ = the amplitude of the jth observation of reflection i
   *  <F>  = the mean of the amplitudes of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_F_obs?: Maybe<Scalars['Float']>;
  /**
   * The smallest value for the interplanar spacings for
   *  the reflection data. This is called the highest resolution.
   */
  readonly d_resolution_high?: Maybe<Scalars['Float']>;
  /**
   * The largest value for the interplanar spacings for the
   *  reflection data. This is called the lowest resolution.
   */
  readonly d_resolution_low?: Maybe<Scalars['Float']>;
  /**
   * A description of special aspects of the data-reduction
   *  procedures.
   * 
   * Examples:
   * Merging and scaling based on only those
   *                                   reflections with I > sig(I).
   */
  readonly data_reduction_details?: Maybe<Scalars['String']>;
  /**
   * The method used for data reduction.
   * 
   *  Note that this is not the computer program used, which is
   *  described in the SOFTWARE category, but the method
   *  itself.
   * 
   *  This data item should be used to describe significant
   *  methodological options used within the data-reduction programs.
   * 
   * Examples:
   * Profile fitting by method of Kabsch (1987).
   *                                   Scaling used spherical harmonic coefficients.
   */
  readonly data_reduction_method?: Maybe<Scalars['String']>;
  /**
   * A description of reflection data not covered by other data
   *  names. This should include details of the Friedel pairs.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Maximum value of the Miller index h for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_h_max.
   */
  readonly limit_h_max?: Maybe<Scalars['Int']>;
  /**
   * Minimum value of the Miller index h for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_h_min.
   */
  readonly limit_h_min?: Maybe<Scalars['Int']>;
  /**
   * Maximum value of the Miller index k for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_k_max.
   */
  readonly limit_k_max?: Maybe<Scalars['Int']>;
  /**
   * Minimum value of the Miller index k for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_k_min.
   */
  readonly limit_k_min?: Maybe<Scalars['Int']>;
  /**
   * Maximum value of the Miller index l for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_l_max.
   */
  readonly limit_l_max?: Maybe<Scalars['Int']>;
  /**
   * Minimum value of the Miller index l for the reflection data. This
   *  need not have the same value as _diffrn_reflns.limit_l_min.
   */
  readonly limit_l_min?: Maybe<Scalars['Int']>;
  /**
   * The total number of reflections in the REFLN list (not the
   *  DIFFRN_REFLN list). This number may contain Friedel-equivalent
   *  reflections according to the nature of the structure and the
   *  procedures used. The item _reflns.details describes the
   *  reflection data.
   */
  readonly number_all?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections in the REFLN list (not the DIFFRN_REFLN
   *  list) classified as observed (see _reflns.observed_criterion).
   *  This number may contain Friedel-equivalent reflections according
   *  to the nature of the structure and the procedures used.
   */
  readonly number_obs?: Maybe<Scalars['Int']>;
  /**
   * The criterion used to classify a reflection as 'observed'. This
   *  criterion is usually expressed in terms of a sigma(I) or
   *  sigma(F) threshold.
   * 
   * Examples:
   * >2sigma(I)
   */
  readonly observed_criterion?: Maybe<Scalars['String']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as an upper limit for the value of F.
   */
  readonly observed_criterion_F_max?: Maybe<Scalars['Float']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as a lower limit for the value of F.
   */
  readonly observed_criterion_F_min?: Maybe<Scalars['Float']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as an upper limit for the value of I.
   */
  readonly observed_criterion_I_max?: Maybe<Scalars['Float']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as a lower limit for the value of I.
   */
  readonly observed_criterion_I_min?: Maybe<Scalars['Float']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as a multiple of the value of sigma(F).
   */
  readonly observed_criterion_sigma_F?: Maybe<Scalars['Float']>;
  /**
   * The criterion used to classify a reflection as 'observed'
   *  expressed as a multiple of the value of sigma(I).
   */
  readonly observed_criterion_sigma_I?: Maybe<Scalars['Float']>;
  /**
   * The Pearson's correlation coefficient expressed as a decimal value
   *               between the average intensities from randomly selected 
   *               half-datasets.
   * 
   * 	      Ref: Karplus & Diederichs (2012), Science 336, 1030-33
   */
  readonly pdbx_CC_half?: Maybe<Scalars['Float']>;
  /**
   * R split measures the agreement between the sets of intensities created by merging 
   *               odd- and even-numbered images  from the overall data.          
   * 
   * 	      Ref: T. A. White, R. A. Kirian, A. V. Martin, A. Aquila, K. Nass, A. Barty 
   *               and H. N. Chapman (2012), J. Appl. Cryst. 45, 335-341
   */
  readonly pdbx_R_split?: Maybe<Scalars['Float']>;
  /**
   * The R value for merging intensities satisfying the observed
   *  criteria in this data set.
   */
  readonly pdbx_Rmerge_I_obs?: Maybe<Scalars['Float']>;
  /**
   * The precision-indicating merging R factor value Rpim,
   *  for merging all intensities in this data set.
   * 
   *         sum~i~ [1/(N~i~ - 1)]1/2^ sum~j~ | I~j~ - <I~i~> |
   *  Rpim = --------------------------------------------------
   *                       sum~i~ ( sum~j~ I~j~ )
   * 
   *  I~j~   = the intensity of the jth observation of reflection i
   *  <I~i~> = the mean of the intensities of all observations
   *           of reflection i
   *  N~i~   = the redundancy (the number of times reflection i
   *           has been measured).
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection.
   * 
   *  Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
   *       Biol. 4, 269-275.
   *       Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
   *       30, 203-205.
   *       Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
   */
  readonly pdbx_Rpim_I_all?: Maybe<Scalars['Float']>;
  /**
   * The redundancy-independent merging R factor value Rrim,
   *               also denoted Rmeas, for merging all intensities in this
   *               data set.
   * 
   *                      sum~i~ [N~i~/(N~i~ - 1)]1/2^ sum~j~ | I~j~ - <I~i~> |
   *               Rrim = ----------------------------------------------------
   *                                   sum~i~ ( sum~j~ I~j~ )
   * 
   *               I~j~   = the intensity of the jth observation of reflection i
   *               <I~i~> = the mean of the intensities of all observations of
   *                        reflection i
   * 	       N~i~   = the redundancy (the number of times reflection i
   *                        has been measured).
   * 
   *               sum~i~ is taken over all reflections
   *               sum~j~ is taken over all observations of each reflection.
   * 
   *               Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
   *                    Biol. 4, 269-275.
   *                    Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
   *                    30, 203-205.
   *                    Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
   */
  readonly pdbx_Rrim_I_all?: Maybe<Scalars['Float']>;
  /**
   * The R sym value as a decimal number.
   * 
   * Examples:
   * 0.02
   */
  readonly pdbx_Rsym_value?: Maybe<Scalars['Float']>;
  /** Overall  Chi-squared statistic. */
  readonly pdbx_chi_squared?: Maybe<Scalars['Float']>;
  /**
   * An identifier for the diffraction data set for this set of summary statistics.
   * 
   *  Multiple diffraction data sets entered as a comma separated list.
   */
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The ratio of the average intensity to the average uncertainty,
   *  <I>/<sigma(I)>.
   */
  readonly pdbx_netI_over_av_sigmaI?: Maybe<Scalars['Float']>;
  /**
   * The mean of the ratio of the intensities to their
   *  standard uncertainties, <I/sigma(I)>.
   */
  readonly pdbx_netI_over_sigmaI?: Maybe<Scalars['Float']>;
  /**
   * Total number of measured reflections.
   * 
   * Examples:
   * 23000, 140000
   */
  readonly pdbx_number_measured_all?: Maybe<Scalars['Int']>;
  /** An ordinal identifier for this set of reflection statistics. */
  readonly pdbx_ordinal: Scalars['Int'];
  /** Overall redundancy for this data set (%). */
  readonly pdbx_redundancy?: Maybe<Scalars['Float']>;
  /** Number of reflections rejected in scaling operations. */
  readonly pdbx_scaling_rejects?: Maybe<Scalars['Int']>;
  /**
   * The percentage of geometrically possible reflections represented
   *  by reflections that satisfy the resolution limits established
   *  by _reflns.d_resolution_high and _reflns.d_resolution_low and
   *  the observation limit established by
   *  _reflns.observed_criterion.
   */
  readonly percent_possible_obs?: Maybe<Scalars['Float']>;
  /**
   * The value of _reflns.phase_calculation_details describes a
   *  special details about calculation of phases in _refln.phase_calc.
   * 
   * Examples:
   * From model, NCS averaging, Solvent flipping, Solvent flattening, Multiple crystal averaging, Multiple phase modification, Other phase modification
   */
  readonly phase_calculation_details?: Maybe<Scalars['String']>;
};

export type RcsbLatestRevision = {
  /** The major version number of the latest revision. */
  readonly major_revision?: Maybe<Scalars['Int']>;
  /** The minor version number of the latest revision. */
  readonly minor_revision?: Maybe<Scalars['Int']>;
  /** The release date of the latest revision item. */
  readonly revision_date?: Maybe<Scalars['Date']>;
};

export type RcsbPolymerEntityRcsbEcLineage = {
  /** Members of the enzyme classification lineage as parent classification hierarchy depth (1-N). */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the enzyme classification lineage as parent classification codes. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the enzyme classification lineage as parent classification names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbAssemblyInfo = {
  /** Entity identifier for the container. */
  readonly assembly_id?: Maybe<Scalars['String']>;
  /** The assembly atomic coordinate count. */
  readonly atom_count?: Maybe<Scalars['Int']>;
  /** The assembly branched entity atomic coordinate count. */
  readonly branched_atom_count?: Maybe<Scalars['Int']>;
  /** The number of distinct branched entities in the generated assembly. */
  readonly branched_entity_count?: Maybe<Scalars['Int']>;
  /**
   * The number of branched instances in the generated assembly data set.
   *  This is the total count of branched entity instances generated in the assembly coordinate data.
   */
  readonly branched_entity_instance_count?: Maybe<Scalars['Int']>;
  /** The PDB entry accession code. */
  readonly entry_id: Scalars['String'];
  /**
   * The number of modeled polymer monomers in the assembly coordinate data.
   *  This is the total count of monomers with reported coordinate data for all polymer
   *  entity instances in the generated assembly coordinate data.
   */
  readonly modeled_polymer_monomer_count?: Maybe<Scalars['Int']>;
  /**
   * Nucleic acid polymer entity type categories describing the generated assembly.
   * 
   * Allowable values:
   * DNA (only), DNA/RNA (only), NA-hybrid (only), Other, RNA (only)
   */
  readonly na_polymer_entity_types?: Maybe<Scalars['String']>;
  /** The assembly non-polymer entity atomic coordinate count. */
  readonly nonpolymer_atom_count?: Maybe<Scalars['Int']>;
  /** The number of distinct non-polymer entities in the generated assembly exclusive of solvent. */
  readonly nonpolymer_entity_count?: Maybe<Scalars['Int']>;
  /**
   * The number of non-polymer instances in the generated assembly data set exclusive of solvent.
   *  This is the total count of non-polymer entity instances generated in the assembly coordinate data.
   */
  readonly nonpolymer_entity_instance_count?: Maybe<Scalars['Int']>;
  /** The assembly polymer entity atomic coordinate count. */
  readonly polymer_atom_count?: Maybe<Scalars['Int']>;
  /**
   * Categories describing the polymer entity composition for the generated assembly.
   * 
   * Allowable values:
   * DNA, DNA/RNA, NA-hybrid, NA/oligosaccharide, RNA, heteromeric protein, homomeric protein, oligosaccharide, other, other type composition, other type pair, protein/NA, protein/NA/oligosaccharide, protein/oligosaccharide
   */
  readonly polymer_composition?: Maybe<Scalars['String']>;
  /** The number of distinct polymer entities in the generated assembly. */
  readonly polymer_entity_count?: Maybe<Scalars['Int']>;
  /** The number of distinct DNA polymer entities in the generated assembly. */
  readonly polymer_entity_count_DNA?: Maybe<Scalars['Int']>;
  /** The number of distinct RNA polymer entities in the generated assembly. */
  readonly polymer_entity_count_RNA?: Maybe<Scalars['Int']>;
  /** The number of distinct nucleic acid polymer entities (DNA or RNA) in the generated assembly. */
  readonly polymer_entity_count_nucleic_acid?: Maybe<Scalars['Int']>;
  /** The number of distinct hybrid nucleic acid polymer entities in the generated assembly. */
  readonly polymer_entity_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>;
  /** The number of distinct protein polymer entities in the generated assembly. */
  readonly polymer_entity_count_protein?: Maybe<Scalars['Int']>;
  /**
   * The number of polymer instances in the generated assembly data set.
   *  This is the total count of polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count?: Maybe<Scalars['Int']>;
  /**
   * The number of DNA polymer instances in the generated assembly data set.
   *  This is the total count of DNA polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count_DNA?: Maybe<Scalars['Int']>;
  /**
   * The number of RNA polymer instances in the generated assembly data set.
   *  This is the total count of RNA polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count_RNA?: Maybe<Scalars['Int']>;
  /**
   * The number of nucleic acid polymer instances in the generated assembly data set.
   *  This is the total count of nucleic acid polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count_nucleic_acid?: Maybe<Scalars['Int']>;
  /**
   * The number of hybrid nucleic acide polymer instances in the generated assembly data set.
   *  This is the total count of hybrid nucleic acid polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>;
  /**
   * The number of protein polymer instances in the generated assembly data set.
   *  This is the total count of protein polymer entity instances generated in the assembly coordinate data.
   */
  readonly polymer_entity_instance_count_protein?: Maybe<Scalars['Int']>;
  /**
   * The number of polymer monomers in sample entity instances comprising the assembly data set.
   *  This is the total count of monomers for all polymer entity instances
   *  in the generated assembly coordinate data.
   */
  readonly polymer_monomer_count?: Maybe<Scalars['Int']>;
  /**
   * Selected polymer entity type categories describing the generated assembly.
   * 
   * Allowable values:
   * Nucleic acid (only), Other, Protein (only), Protein/NA
   */
  readonly selected_polymer_entity_types?: Maybe<Scalars['String']>;
  /** The assembly solvent atomic coordinate count. */
  readonly solvent_atom_count?: Maybe<Scalars['Int']>;
  /** The number of distinct solvent entities in the generated assembly. */
  readonly solvent_entity_count?: Maybe<Scalars['Int']>;
  /**
   * The number of solvent instances in the generated assembly data set.
   *  This is the total count of solvent entity instances generated in the assembly coordinate data.
   */
  readonly solvent_entity_instance_count?: Maybe<Scalars['Int']>;
  /**
   * The number of unmodeled polymer monomers in the assembly coordinate data. This is
   *  the total count of monomers with unreported coordinate data for all polymer
   *  entity instances in the generated assembly coordinate data.
   */
  readonly unmodeled_polymer_monomer_count?: Maybe<Scalars['Int']>;
};

export type RcsbBindingAffinity = {
  /**
   * Ligand identifier.
   * 
   * Examples:
   * 0WE, SPE, CL
   */
  readonly comp_id: Scalars['String'];
  /** Link to external resource referencing the data. */
  readonly link: Scalars['String'];
  /**
   * The resource name for the related binding affinity reference.
   * 
   * Allowable values:
   * PDBBind, Binding MOAD, BindingDB
   */
  readonly provenance_code: Scalars['String'];
  /**
   * Data point provided by BindingDB. Percent identity between PDB sequence and reference sequence.
   * 
   * Examples:
   * null, null, null
   */
  readonly reference_sequence_identity?: Maybe<Scalars['Int']>;
  /**
   * Binding affinity symbol indicating approximate or precise strength of the binding.
   * 
   * Examples:
   * ~, =, >, <, >=, <=
   */
  readonly symbol?: Maybe<Scalars['String']>;
  /**
   * Binding affinity measurement given in one of the following types:  The concentration constants: IC50: the concentration of ligand that reduces enzyme activity by 50%;  EC50: the concentration of compound that generates a half-maximal response;  The binding constant:  Kd: dissociation constant;  Ka: association constant;  Ki: enzyme inhibition constant;  The thermodynamic parameters:  delta G: Gibbs free energy of binding (for association reaction);  delta H: change in enthalpy associated with a chemical reaction;  delta S: change in entropy associated with a chemical reaction.
   * 
   * Examples:
   * IC50, EC50, Kd, Ka, Ki
   */
  readonly type: Scalars['String'];
  /**
   * Binding affinity unit.  Dissociation constant Kd is normally in molar units (or millimolar , micromolar, nanomolar, etc).  Association constant Ka is normally expressed in inverse molar units (e.g. M-1).
   * 
   * Examples:
   * nM, kJ/mol
   */
  readonly unit: Scalars['String'];
  /** Binding affinity value between a ligand and its target molecule. */
  readonly value: Scalars['Float'];
};

export type EmEntityAssembly = {
  /**
   * Additional details about the component.
   * 
   * Examples:
   * Fab fragment generated by proteolytic cleavage of LA2 IgG antibody.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * macromolecules associated with this component, if defined
   *  as comma separated list of entity ids (integers).
   */
  readonly entity_id_list?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The value of _em_entity_assembly.id identifies
   *  one component of the complex.
   */
  readonly id: Scalars['String'];
  /**
   * Name of this component in the observed assembly.
   * 
   * Examples:
   * Ternary complex of alpha-tubulin with tubulin folding cofactors TBCE and TBCB, 80S Ribosome bound to emetine, messenger RNA, initiation factor 2, GroEL, antibody Fab fragment
   */
  readonly name?: Maybe<Scalars['String']>;
  /** oligomeric details */
  readonly oligomeric_details?: Maybe<Scalars['String']>;
  /**
   * The parent of this assembly.
   *  This data item is an internal category pointer to _em_entity_assembly.id.
   *  By convention, the full assembly (top of hierarchy) is assigned parent id 0 (zero).
   */
  readonly parent_id?: Maybe<Scalars['Int']>;
  /**
   * The assembly type.
   * 
   * Allowable values:
   * MULTIPLE SOURCES, NATURAL, RECOMBINANT
   */
  readonly source?: Maybe<Scalars['String']>;
  /**
   * Alternative name of the component.
   * 
   * Examples:
   * FADV-1
   */
  readonly synonym?: Maybe<Scalars['String']>;
  /**
   * A description of types of components of the
   *  assembly of the biological structure.
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbUniprotProteinSourceOrganism = {
  /** Historical record of the data attribute. */
  readonly provenance_code: Scalars['String'];
  /** The scientific name of the organism in which a protein occurs. */
  readonly scientific_name: Scalars['String'];
  /** NCBI Taxonomy identifier for the organism in which a protein occurs. */
  readonly taxonomy_id?: Maybe<Scalars['Int']>;
};

export type CoreDrugbank = {
  readonly drugbank_container_identifiers?: Maybe<DrugbankContainerIdentifiers>;
  readonly drugbank_info?: Maybe<DrugbankInfo>;
  readonly drugbank_target?: Maybe<ReadonlyArray<Maybe<DrugbankTarget>>>;
};

export type PdbxAuditRevisionGroup = {
  /**
   * The type of file that the pdbx_audit_revision_history record refers to.
   * 
   * Allowable values:
   * Chemical component, NMR restraints, NMR shifts, Structure factors, Structure model
   */
  readonly data_content_type: Scalars['String'];
  /**
   * The collection of categories updated with this revision.
   * 
   * Allowable values:
   * Advisory, Atomic model, Author supporting evidence, Data collection, Data processing, Database references, Derived calculations, Experimental data, Experimental preparation, Initial release, Non-polymer description, Other, Polymer sequence, Refinement description, Source and taxonomy, Structure summary, Version format compliance
   */
  readonly group?: Maybe<Scalars['String']>;
  /**
   * A unique identifier for the pdbx_audit_revision_group record.
   * 
   * Examples:
   * 1
   */
  readonly ordinal: Scalars['Int'];
  /**
   * A pointer to  _pdbx_audit_revision_history.ordinal
   * 
   * Examples:
   * 1
   */
  readonly revision_ordinal: Scalars['Int'];
};

export type RcsbChemCompAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbChemCompAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * ATC, Generating Enzyme, Modification Type, PSI-MOD
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type ExptlCrystalGrow = {
  /**
   * This data item is a pointer to _exptl_crystal.id in the
   *  EXPTL_CRYSTAL category.
   */
  readonly crystal_id: Scalars['String'];
  /**
   * A description of special aspects of the crystal growth.
   * 
   * Examples:
   * Solution 2 was prepared as a well solution and
   *                                   mixed. A droplet containing 2 \ml of solution
   *                                   1 was delivered onto a cover slip; 2 \ml of
   *                                   solution 2 was added to the droplet without
   *                                   mixing., Crystal plates were originally stored at room
   *                                   temperature for 1 week but no nucleation
   *                                   occurred. They were then transferred to 4
   *                                   degrees C, at which temperature well formed
   *                                   single crystals grew in 2 days., The dependence on pH for successful crystal
   *                                   growth is very sharp. At pH 7.4 only showers
   *                                   of tiny crystals grew, at pH 7.5 well formed
   *                                   single crystals grew, at pH 7.6 no
   *                                   crystallization occurred at all.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The method used to grow the crystals.
   * 
   * Examples:
   * MICROBATCH, VAPOR DIFFUSION, HANGING DROP
   */
  readonly method?: Maybe<Scalars['String']>;
  /**
   * The pH at which the crystal was grown. If more than one pH was
   *  employed during the crystallization process, the final pH should
   *  be noted here and the protocol involving multiple pH values
   *  should be described in _exptl_crystal_grow.details.
   * 
   * Examples:
   * 7.4, 7.6, 4.3
   */
  readonly pH?: Maybe<Scalars['Float']>;
  /**
   * Text description of crystal growth procedure.
   * 
   * Examples:
   * PEG 4000, potassium phosphate, magnesium chloride, cacodylate
   */
  readonly pdbx_details?: Maybe<Scalars['String']>;
  /**
   * The range of pH values at which the crystal was grown.   Used when
   *  a point estimate of pH is not appropriate.
   * 
   * Examples:
   * 5.6 - 6.4
   */
  readonly pdbx_pH_range?: Maybe<Scalars['String']>;
  /**
   * The temperature in kelvins at which the crystal was grown.
   *  If more than one temperature was employed during the
   *  crystallization process, the final temperature should be noted
   *  here and the protocol  involving multiple temperatures should be
   *  described in _exptl_crystal_grow.details.
   */
  readonly temp?: Maybe<Scalars['Float']>;
  /**
   * A description of special aspects of temperature control during
   *  crystal growth.
   */
  readonly temp_details?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityAlignAlignedRegions = {
  /** An identifier for the monomer in the entity sequence at which this segment of the alignment begins. */
  readonly entity_beg_seq_id?: Maybe<Scalars['Int']>;
  /** An length of the this segment of the alignment. */
  readonly length?: Maybe<Scalars['Int']>;
  /** An identifier for the monomer in the reference sequence at which this segment of the alignment begins. */
  readonly ref_beg_seq_id?: Maybe<Scalars['Int']>;
};

export type RcsbUniprotAlignmentsCoreEntityAlignments = {
  /** Aligned region */
  readonly aligned_regions?: Maybe<ReadonlyArray<Maybe<CoreEntityAlignmentsAlignedRegions>>>;
  /** core_entity identifiers */
  readonly core_entity_identifiers?: Maybe<CoreEntityAlignmentsCoreEntityIdentifiers>;
  /** Alignment scores */
  readonly scores?: Maybe<CoreEntityAlignmentsScores>;
};

export type RcsbNonpolymerEntityInstanceContainerIdentifiers = {
  /** Instance identifier for this container. */
  readonly asym_id: Scalars['String'];
  /** Author instance identifier for this container. */
  readonly auth_asym_id?: Maybe<Scalars['String']>;
  /** Residue number for non-polymer entity instance. */
  readonly auth_seq_id?: Maybe<Scalars['String']>;
  /** Component identifier for non-polymer entity instance. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** Entity identifier for the container. */
  readonly entity_id?: Maybe<Scalars['String']>;
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /**
   * A unique identifier for each object in this entity instance container formed by
   *  an 'dot' (.) separated concatenation of entry and entity instance identifiers.
   */
  readonly rcsb_id?: Maybe<Scalars['String']>;
};

export type PdbxReferenceMoleculeDetails = {
  /**
   * The value of _pdbx_reference_molecule_details.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_list.family_prd_id' in category PDBX_REFERENCE_MOLECULE_FAMILY.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_molecule_details.ordinal is an ordinal that 
   *  distinguishes each descriptive text for this entity.
   */
  readonly ordinal: Scalars['Int'];
  /** A data source of this information (e.g. PubMed, Merck Index) */
  readonly source?: Maybe<Scalars['String']>;
  /** A identifier within the data source for this information. */
  readonly source_id?: Maybe<Scalars['String']>;
  /** The text of the description of special aspects of the entity. */
  readonly text?: Maybe<Scalars['String']>;
};

export type PdbxReferenceMoleculeSynonyms = {
  /**
   * The value of _pdbx_reference_molecule_synonyms.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_list.family_prd_id in category PDBX_REFERENCE_MOLECULE_FAMILY_LIST.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * A synonym name for the entity.
   * 
   * Examples:
   * thiostrepton
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_molecule_synonyms.ordinal is an ordinal
   * 	       to distinguish synonyms for this entity.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The value of _pdbx_reference_molecule_synonyms.prd_id is a reference
   * 	       _pdbx_reference_molecule.prd_id in the  PDBX_REFERENCE_MOLECULE category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The source of this synonym name for the entity.
   * 
   * Examples:
   * CAS
   */
  readonly source?: Maybe<Scalars['String']>;
};

export type RcsbPolymerInstanceFeatureSummary = {
  /** The feature count. */
  readonly count?: Maybe<Scalars['Int']>;
  /** The fractional feature coverage relative to the full entity sequence. */
  readonly coverage?: Maybe<Scalars['Float']>;
  /** The maximum feature length. */
  readonly maximum_length?: Maybe<Scalars['Int']>;
  /** The maximum feature value. */
  readonly maximum_value?: Maybe<Scalars['Float']>;
  /** The minimum feature length. */
  readonly minimum_length?: Maybe<Scalars['Int']>;
  /** The minimum feature value. */
  readonly minimum_value?: Maybe<Scalars['Float']>;
  /**
   * Type or category of the feature.
   * 
   * Allowable values:
   * ANGLE_OUTLIER, BINDING_SITE, BOND_OUTLIER, CATH, CIS-PEPTIDE, HELIX_P, MOGUL_ANGLE_OUTLIER, MOGUL_BOND_OUTLIER, RAMACHANDRAN_OUTLIER, ROTAMER_OUTLIER, RSRCC_OUTLIER, RSRZ_OUTLIER, SCOP, SHEET, UNASSIGNED_SEC_STRUCT, UNOBSERVED_ATOM_XYZ, UNOBSERVED_RESIDUE_XYZ, ZERO_OCCUPANCY_ATOM_XYZ, ZERO_OCCUPANCY_RESIDUE_XYZ
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type CoreEntry = {
  /** Get all assemblies for this PDB entry. */
  readonly assemblies?: Maybe<ReadonlyArray<Maybe<CoreAssembly>>>;
  readonly audit_author?: Maybe<ReadonlyArray<Maybe<AuditAuthor>>>;
  readonly cell?: Maybe<Cell>;
  readonly citation?: Maybe<ReadonlyArray<Maybe<Citation>>>;
  readonly diffrn?: Maybe<ReadonlyArray<Maybe<Diffrn>>>;
  readonly diffrn_detector?: Maybe<ReadonlyArray<Maybe<DiffrnDetector>>>;
  readonly diffrn_radiation?: Maybe<ReadonlyArray<Maybe<DiffrnRadiation>>>;
  readonly diffrn_source?: Maybe<ReadonlyArray<Maybe<DiffrnSource>>>;
  readonly em_2d_crystal_entity?: Maybe<ReadonlyArray<Maybe<Em2dCrystalEntity>>>;
  readonly em_3d_crystal_entity?: Maybe<ReadonlyArray<Maybe<Em3dCrystalEntity>>>;
  readonly em_3d_fitting?: Maybe<ReadonlyArray<Maybe<Em3dFitting>>>;
  readonly em_3d_fitting_list?: Maybe<ReadonlyArray<Maybe<Em3dFittingList>>>;
  readonly em_3d_reconstruction?: Maybe<ReadonlyArray<Maybe<Em3dReconstruction>>>;
  readonly em_ctf_correction?: Maybe<ReadonlyArray<Maybe<EmCtfCorrection>>>;
  readonly em_diffraction?: Maybe<ReadonlyArray<Maybe<EmDiffraction>>>;
  readonly em_diffraction_shell?: Maybe<ReadonlyArray<Maybe<EmDiffractionShell>>>;
  readonly em_diffraction_stats?: Maybe<ReadonlyArray<Maybe<EmDiffractionStats>>>;
  readonly em_embedding?: Maybe<ReadonlyArray<Maybe<EmEmbedding>>>;
  readonly em_entity_assembly?: Maybe<ReadonlyArray<Maybe<EmEntityAssembly>>>;
  readonly em_experiment?: Maybe<EmExperiment>;
  readonly em_helical_entity?: Maybe<ReadonlyArray<Maybe<EmHelicalEntity>>>;
  readonly em_image_recording?: Maybe<ReadonlyArray<Maybe<EmImageRecording>>>;
  readonly em_imaging?: Maybe<ReadonlyArray<Maybe<EmImaging>>>;
  readonly em_particle_selection?: Maybe<ReadonlyArray<Maybe<EmParticleSelection>>>;
  readonly em_single_particle_entity?: Maybe<ReadonlyArray<Maybe<EmSingleParticleEntity>>>;
  readonly em_software?: Maybe<ReadonlyArray<Maybe<EmSoftware>>>;
  readonly em_specimen?: Maybe<ReadonlyArray<Maybe<EmSpecimen>>>;
  readonly em_staining?: Maybe<ReadonlyArray<Maybe<EmStaining>>>;
  readonly em_vitrification?: Maybe<ReadonlyArray<Maybe<EmVitrification>>>;
  readonly entry?: Maybe<Entry>;
  readonly exptl?: Maybe<ReadonlyArray<Maybe<Exptl>>>;
  readonly exptl_crystal?: Maybe<ReadonlyArray<Maybe<ExptlCrystal>>>;
  readonly exptl_crystal_grow?: Maybe<ReadonlyArray<Maybe<ExptlCrystalGrow>>>;
  /** Get all non-polymer (non-solvent) entities for this PDB entry. */
  readonly nonpolymer_entities?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntity>>>;
  readonly pdbx_SG_project?: Maybe<ReadonlyArray<Maybe<PdbxSgProject>>>;
  readonly pdbx_audit_revision_category?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionCategory>>>;
  readonly pdbx_audit_revision_details?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionDetails>>>;
  readonly pdbx_audit_revision_group?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionGroup>>>;
  readonly pdbx_audit_revision_history?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionHistory>>>;
  readonly pdbx_audit_revision_item?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionItem>>>;
  readonly pdbx_audit_support?: Maybe<ReadonlyArray<Maybe<PdbxAuditSupport>>>;
  readonly pdbx_database_PDB_obs_spr?: Maybe<ReadonlyArray<Maybe<PdbxDatabasePdbObsSpr>>>;
  readonly pdbx_database_related?: Maybe<ReadonlyArray<Maybe<PdbxDatabaseRelated>>>;
  readonly pdbx_database_status?: Maybe<PdbxDatabaseStatus>;
  readonly pdbx_deposit_group?: Maybe<ReadonlyArray<Maybe<PdbxDepositGroup>>>;
  readonly pdbx_molecule_features?: Maybe<ReadonlyArray<Maybe<PdbxMoleculeFeatures>>>;
  readonly pdbx_nmr_details?: Maybe<PdbxNmrDetails>;
  readonly pdbx_nmr_ensemble?: Maybe<PdbxNmrEnsemble>;
  readonly pdbx_nmr_exptl?: Maybe<ReadonlyArray<Maybe<PdbxNmrExptl>>>;
  readonly pdbx_nmr_exptl_sample_conditions?: Maybe<ReadonlyArray<Maybe<PdbxNmrExptlSampleConditions>>>;
  readonly pdbx_nmr_refine?: Maybe<ReadonlyArray<Maybe<PdbxNmrRefine>>>;
  readonly pdbx_nmr_representative?: Maybe<PdbxNmrRepresentative>;
  readonly pdbx_nmr_sample_details?: Maybe<ReadonlyArray<Maybe<PdbxNmrSampleDetails>>>;
  readonly pdbx_nmr_software?: Maybe<ReadonlyArray<Maybe<PdbxNmrSoftware>>>;
  readonly pdbx_nmr_spectrometer?: Maybe<ReadonlyArray<Maybe<PdbxNmrSpectrometer>>>;
  readonly pdbx_serial_crystallography_data_reduction?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographyDataReduction>>>;
  readonly pdbx_serial_crystallography_measurement?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographyMeasurement>>>;
  readonly pdbx_serial_crystallography_sample_delivery?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDelivery>>>;
  readonly pdbx_serial_crystallography_sample_delivery_fixed_target?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDeliveryFixedTarget>>>;
  readonly pdbx_serial_crystallography_sample_delivery_injection?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDeliveryInjection>>>;
  readonly pdbx_soln_scatter?: Maybe<ReadonlyArray<Maybe<PdbxSolnScatter>>>;
  readonly pdbx_soln_scatter_model?: Maybe<ReadonlyArray<Maybe<PdbxSolnScatterModel>>>;
  readonly pdbx_vrpt_summary?: Maybe<PdbxVrptSummary>;
  /** Get all polymer entities for this PDB entry. */
  readonly polymer_entities?: Maybe<ReadonlyArray<Maybe<CorePolymerEntity>>>;
  /** Get literature information from PubMed database. */
  readonly pubmed?: Maybe<CorePubmed>;
  readonly rcsb_accession_info?: Maybe<RcsbAccessionInfo>;
  /** The list of content types associated with this entry. */
  readonly rcsb_associated_holdings?: Maybe<CurrentEntry>;
  readonly rcsb_binding_affinity?: Maybe<ReadonlyArray<Maybe<RcsbBindingAffinity>>>;
  readonly rcsb_entry_container_identifiers: RcsbEntryContainerIdentifiers;
  readonly rcsb_entry_info: RcsbEntryInfo;
  readonly rcsb_external_references?: Maybe<ReadonlyArray<Maybe<RcsbExternalReferences>>>;
  /** A unique identifier for each object in this entry container. */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_primary_citation?: Maybe<RcsbPrimaryCitation>;
  readonly refine?: Maybe<ReadonlyArray<Maybe<Refine>>>;
  readonly refine_analyze?: Maybe<ReadonlyArray<Maybe<RefineAnalyze>>>;
  readonly refine_hist?: Maybe<ReadonlyArray<Maybe<RefineHist>>>;
  readonly refine_ls_restr?: Maybe<ReadonlyArray<Maybe<RefineLsRestr>>>;
  readonly reflns?: Maybe<ReadonlyArray<Maybe<Reflns>>>;
  readonly reflns_shell?: Maybe<ReadonlyArray<Maybe<ReflnsShell>>>;
  readonly software?: Maybe<ReadonlyArray<Maybe<Software>>>;
  readonly struct?: Maybe<Struct>;
  readonly struct_keywords?: Maybe<StructKeywords>;
  readonly symmetry?: Maybe<Symmetry>;
};

export type ReflnsShell = {
  /**
   * Residual factor Rmerge for all reflections that satisfy the
   *  resolution limits established by _reflns_shell.d_res_high and
   *  _reflns_shell.d_res_low.
   * 
   *              sum~i~(sum~j~|F~j~ - <F>|)
   *  Rmerge(F) = --------------------------
   *                   sum~i~(sum~j~<F>)
   * 
   *  F~j~ = the amplitude of the jth observation of reflection i
   *  <F>  = the mean of the amplitudes of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_F_all?: Maybe<Scalars['Float']>;
  /**
   * Residual factor Rmerge for reflections that satisfy the
   *  resolution limits established by _reflns_shell.d_res_high and
   *  _reflns_shell.d_res_low and the observation criterion
   *  established by _reflns.observed_criterion.
   * 
   *              sum~i~(sum~j~|F~j~ - <F>|)
   *  Rmerge(F) = --------------------------
   *                   sum~i~(sum~j~<F>)
   * 
   *  F~j~ = the amplitude of the jth observation of reflection i
   *  <F>  = the mean of the amplitudes of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_F_obs?: Maybe<Scalars['Float']>;
  /**
   * The value of Rmerge(I) for all reflections in a given shell.
   * 
   *              sum~i~(sum~j~|I~j~ - <I>|)
   *  Rmerge(I) = --------------------------
   *                  sum~i~(sum~j~<I>)
   * 
   *  I~j~ = the intensity of the jth observation of reflection i
   *  <I>  = the mean of the intensities of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_I_all?: Maybe<Scalars['Float']>;
  /**
   * The value of Rmerge(I) for reflections classified as 'observed'
   *  (see _reflns.observed_criterion) in a given shell.
   * 
   *              sum~i~(sum~j~|I~j~ - <I>|)
   *  Rmerge(I) = --------------------------
   *                  sum~i~(sum~j~<I>)
   * 
   *  I~j~ = the intensity of the jth observation of reflection i
   *  <I>  = the mean of the intensities of all observations of
   *         reflection i
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection
   */
  readonly Rmerge_I_obs?: Maybe<Scalars['Float']>;
  /**
   * The smallest value in angstroms for the interplanar spacings
   *  for the reflections in this shell. This is called the highest
   *  resolution.
   */
  readonly d_res_high?: Maybe<Scalars['Float']>;
  /**
   * The highest value in angstroms for the interplanar spacings
   *  for the reflections in this shell. This is called the lowest
   *  resolution.
   */
  readonly d_res_low?: Maybe<Scalars['Float']>;
  /**
   * The ratio of the mean of the intensities of all reflections
   *  in this shell to the mean of the standard uncertainties of the
   *  intensities of all reflections in this shell.
   */
  readonly meanI_over_sigI_all?: Maybe<Scalars['Float']>;
  /**
   * The ratio of the mean of the intensities of the reflections
   *  classified as 'observed' (see _reflns.observed_criterion) in
   *  this shell to the mean of the standard uncertainties of the
   *  intensities of the 'observed' reflections in this
   *  shell.
   */
  readonly meanI_over_sigI_obs?: Maybe<Scalars['Float']>;
  /**
   * The ratio of the mean of the intensities of all reflections
   *  in this shell to the mean of the standard uncertainties of the
   *  intensities of all reflections in this shell.
   */
  readonly meanI_over_uI_all?: Maybe<Scalars['Float']>;
  /**
   * The total number of reflections measured for this
   *  shell.
   */
  readonly number_measured_all?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections classified as 'observed'
   *  (see _reflns.observed_criterion) for this
   *  shell.
   */
  readonly number_measured_obs?: Maybe<Scalars['Int']>;
  /**
   * The number of unique reflections it is possible to measure in
   *  this shell.
   */
  readonly number_possible?: Maybe<Scalars['Int']>;
  /**
   * The total number of measured reflections which are symmetry-
   *  unique after merging for this shell.
   */
  readonly number_unique_all?: Maybe<Scalars['Int']>;
  /**
   * The total number of measured reflections classified as 'observed'
   *  (see _reflns.observed_criterion) which are symmetry-unique
   *  after merging for this shell.
   */
  readonly number_unique_obs?: Maybe<Scalars['Int']>;
  /**
   * The Pearson's correlation coefficient expressed as a decimal value
   *               between the average intensities from randomly selected
   *               half-datasets within the resolution shell.
   * 
   * 	      Ref: Karplus & Diederichs (2012), Science 336, 1030-33
   */
  readonly pdbx_CC_half?: Maybe<Scalars['Float']>;
  /**
   * R split measures the agreement between the sets of intensities created by merging 
   *               odd- and even-numbered images from the data within the resolution shell.
   * 
   * 	      Ref: T. A. White, R. A. Kirian, A. V. Martin, A. Aquila, K. Nass, 
   * 	      A. Barty and H. N. Chapman (2012), J. Appl. Cryst. 45, 335-341
   */
  readonly pdbx_R_split?: Maybe<Scalars['Float']>;
  /**
   * The precision-indicating merging R factor value Rpim,
   *  for merging all intensities in a given shell.
   * 
   *         sum~i~ [1/(N~i~ - 1)]1/2^ sum~j~ | I~j~ - <I~i~> |
   *  Rpim = --------------------------------------------------
   *                       sum~i~ ( sum~j~ I~j~ )
   * 
   *  I~j~   = the intensity of the jth observation of reflection i
   *  <I~i~> = the mean of the intensities of all observations of
   *           reflection i
   *  N~i~   = the redundancy (the number of times reflection i
   *           has been measured).
   * 
   *  sum~i~ is taken over all reflections
   *  sum~j~ is taken over all observations of each reflection.
   * 
   *  Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
   *       Biol. 4, 269-275.
   *       Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
   *       30, 203-205.
   *       Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
   */
  readonly pdbx_Rpim_I_all?: Maybe<Scalars['Float']>;
  /**
   * The redundancy-independent merging R factor value Rrim,
   *               also denoted Rmeas, for merging all intensities in a
   *               given shell.
   * 
   *                      sum~i~ [N~i~ /( N~i~ - 1)]1/2^ sum~j~ | I~j~ - <I~i~> |
   *               Rrim = --------------------------------------------------------
   *                                    sum~i~ ( sum~j~ I~j~ )
   * 
   *               I~j~   = the intensity of the jth observation of reflection i
   *               <I~i~> = the mean of the intensities of all observations of
   *                        reflection i
   * 	      N~i~   = the redundancy (the number of times reflection i
   *                        has been measured).
   * 
   *               sum~i~ is taken over all reflections
   *               sum~j~ is taken over all observations of each reflection.
   * 
   *               Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
   *                    Biol. 4, 269-275.
   *                    Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
   *                    30, 203-205.
   *                    Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
   */
  readonly pdbx_Rrim_I_all?: Maybe<Scalars['Float']>;
  /** R sym value in percent. */
  readonly pdbx_Rsym_value?: Maybe<Scalars['Float']>;
  /** Chi-squared statistic for this resolution shell. */
  readonly pdbx_chi_squared?: Maybe<Scalars['Float']>;
  /**
   * An identifier for the diffraction data set corresponding to this resolution shell.
   * 
   *  Multiple diffraction data sets specified as a comma separated list.
   */
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The mean of the ratio of the intensities to their
   *  standard uncertainties of all reflections in the
   *  resolution shell.
   * 
   *  _reflns_shell.pdbx_netI_over_sigmaI_all =  <I/sigma(I)>
   */
  readonly pdbx_netI_over_sigmaI_all?: Maybe<Scalars['Float']>;
  /**
   * The mean of the ratio of the intensities to their
   *  standard uncertainties of observed reflections
   *  (see _reflns.observed_criterion) in the resolution shell.
   * 
   *  _reflns_shell.pdbx_netI_over_sigmaI_obs =  <I/sigma(I)>
   */
  readonly pdbx_netI_over_sigmaI_obs?: Maybe<Scalars['Float']>;
  /** An ordinal identifier for this resolution shell. */
  readonly pdbx_ordinal: Scalars['Int'];
  /** Redundancy for the current shell. */
  readonly pdbx_redundancy?: Maybe<Scalars['Float']>;
  /**
   * The number of rejected reflections in the resolution 
   *  shell.  Reflections may be rejected from scaling
   *  by setting the observation criterion,
   *  _reflns.observed_criterion.
   */
  readonly pdbx_rejects?: Maybe<Scalars['Int']>;
  /**
   * The percentage of geometrically possible reflections represented
   *  by all reflections measured for this shell.
   */
  readonly percent_possible_all?: Maybe<Scalars['Float']>;
  /**
   * The percentage of geometrically possible reflections represented
   *  by reflections classified as 'observed' (see
   *  _reflns.observed_criterion) for this shell.
   */
  readonly percent_possible_obs?: Maybe<Scalars['Float']>;
};

export type RcsbPubmedMeshDescriptorsLineage = {
  /** Hierarchy depth. */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Identifier for MeSH classification term. */
  readonly id?: Maybe<Scalars['String']>;
  /** MeSH classification term. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerInstanceFeatureFeatureValue = {
  /** The chemical component identifier for the instance of the feature value. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** Specific details about the feature. */
  readonly details?: Maybe<Scalars['String']>;
  /** The reference value of the feature. */
  readonly reference?: Maybe<Scalars['Float']>;
  /** The reported value of the feature. */
  readonly reported?: Maybe<Scalars['Float']>;
  /** The estimated uncertainty of the reported feature value. */
  readonly uncertainty_estimate?: Maybe<Scalars['Float']>;
  /**
   * The type of estimated uncertainty for the reported feature value.
   * 
   * Allowable values:
   * Z-Score
   */
  readonly uncertainty_estimate_type?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityNameCom = {
  /**
   * A common name for the polymer entity.
   * 
   * Examples:
   * HIV protease monomer, hemoglobin alpha chain
   */
  readonly name: Scalars['String'];
};

export type EmSingleParticleEntity = {
  /** Unique category label. */
  readonly id: Scalars['String'];
  /** pointer to _em_image_processing.id. */
  readonly image_processing_id: Scalars['String'];
  /** Point symmetry symbol, either Cn, Dn, T, O, or I */
  readonly point_symmetry?: Maybe<Scalars['String']>;
};

export type RcsbEntryInfo = {
  /** The number of assemblies defined for this entry including the deposited assembly. */
  readonly assembly_count?: Maybe<Scalars['Int']>;
  /** The number of distinct branched entities in the structure entry. */
  readonly branched_entity_count?: Maybe<Scalars['Int']>;
  /** The maximum molecular mass (KDa) of a branched entity in the deposited structure entry. */
  readonly branched_molecular_weight_maximum?: Maybe<Scalars['Float']>;
  /** The minimum molecular mass (KDa) of a branched entity in the deposited structure entry. */
  readonly branched_molecular_weight_minimum?: Maybe<Scalars['Float']>;
  /** The number of cis-peptide linkages per deposited structure model. */
  readonly cis_peptide_count?: Maybe<Scalars['Int']>;
  /** The number of heavy atom coordinates records per deposited structure model. */
  readonly deposited_atom_count?: Maybe<Scalars['Int']>;
  /** The number of model structures deposited. */
  readonly deposited_model_count?: Maybe<Scalars['Int']>;
  /**
   * The number of modeled polymer monomers in the deposited coordinate data.
   *  This is the total count of monomers with reported coordinate data for all polymer
   *  entity instances in the deposited coordinate data.
   */
  readonly deposited_modeled_polymer_monomer_count?: Maybe<Scalars['Int']>;
  /**
   * The number of non-polymer instances in the deposited data set.
   *  This is the total count of non-polymer entity instances reported
   *  per deposited structure model.
   */
  readonly deposited_nonpolymer_entity_instance_count?: Maybe<Scalars['Int']>;
  /**
   * The number of polymer instances in the deposited data set.
   *  This is the total count of polymer entity instances reported
   *  per deposited structure model.
   */
  readonly deposited_polymer_entity_instance_count?: Maybe<Scalars['Int']>;
  /**
   * The number of polymer monomers in sample entity instances in the deposited data set.
   *  This is the total count of monomers for all polymer entity instances reported
   *  per deposited structure model.
   */
  readonly deposited_polymer_monomer_count?: Maybe<Scalars['Int']>;
  /**
   * The number of unmodeled polymer monomers in the deposited coordinate data. This is
   *  the total count of monomers with unreported coordinate data for all polymer
   *  entity instances per deposited structure model.
   */
  readonly deposited_unmodeled_polymer_monomer_count?: Maybe<Scalars['Int']>;
  /** The maximum radiation wavelength in angstroms. */
  readonly diffrn_radiation_wavelength_maximum?: Maybe<Scalars['Float']>;
  /** The minimum radiation wavelength in angstroms. */
  readonly diffrn_radiation_wavelength_minimum?: Maybe<Scalars['Float']>;
  /** The number of disulfide bonds per deposited structure model. */
  readonly disulfide_bond_count?: Maybe<Scalars['Int']>;
  /** The number of distinct polymer, non-polymer, branched molecular, and solvent entities per deposited structure model. */
  readonly entity_count?: Maybe<Scalars['Int']>;
  /**
   * The category of experimental method(s) used to determine the structure entry.
   * 
   * Allowable values:
   * EM, Multiple methods, NMR, Neutron, Other, X-ray
   */
  readonly experimental_method?: Maybe<Scalars['String']>;
  /** The number of experimental methods contributing data to the structure determination. */
  readonly experimental_method_count?: Maybe<Scalars['Int']>;
  /** The number of intermolecular covalent bonds. */
  readonly inter_mol_covalent_bond_count?: Maybe<Scalars['Int']>;
  /** The number of intermolecular metalic bonds. */
  readonly inter_mol_metalic_bond_count?: Maybe<Scalars['Int']>;
  /** The molecular mass (KDa) of polymer and non-polymer entities (exclusive of solvent) in the deposited structure entry. */
  readonly molecular_weight?: Maybe<Scalars['Float']>;
  /**
   * Nucleic acid polymer entity type categories describing the entry.
   * 
   * Allowable values:
   * DNA (only), DNA/RNA (only), NA-hybrid (only), Other, RNA (only)
   */
  readonly na_polymer_entity_types?: Maybe<Scalars['String']>;
  /** Bound nonpolymer components in this entry. */
  readonly nonpolymer_bound_components?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The number of distinct non-polymer entities in the structure entry exclusive of solvent. */
  readonly nonpolymer_entity_count?: Maybe<Scalars['Int']>;
  /** The maximum molecular mass (KDa) of a non-polymer entity in the deposited structure entry. */
  readonly nonpolymer_molecular_weight_maximum?: Maybe<Scalars['Float']>;
  /** The minimum molecular mass (KDa) of a non-polymer entity in the deposited structure entry. */
  readonly nonpolymer_molecular_weight_minimum?: Maybe<Scalars['Float']>;
  /**
   * Categories describing the polymer entity composition for the entry.
   * 
   * Allowable values:
   * DNA, DNA/RNA, NA-hybrid, NA/oligosaccharide, RNA, heteromeric protein, homomeric protein, oligosaccharide, other, other type composition, other type pair, protein/NA, protein/NA/oligosaccharide, protein/oligosaccharide
   */
  readonly polymer_composition?: Maybe<Scalars['String']>;
  /** The number of distinct polymer entities in the structure entry. */
  readonly polymer_entity_count?: Maybe<Scalars['Int']>;
  /** The number of distinct DNA polymer entities. */
  readonly polymer_entity_count_DNA?: Maybe<Scalars['Int']>;
  /** The number of distinct RNA polymer entities. */
  readonly polymer_entity_count_RNA?: Maybe<Scalars['Int']>;
  /** The number of distinct nucleic acid polymer entities (DNA or RNA). */
  readonly polymer_entity_count_nucleic_acid?: Maybe<Scalars['Int']>;
  /** The number of distinct hybrid nucleic acid polymer entities. */
  readonly polymer_entity_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>;
  /** The number of distinct protein polymer entities. */
  readonly polymer_entity_count_protein?: Maybe<Scalars['Int']>;
  /** The number of distinct taxonomies represented among the polymer entities in the entry. */
  readonly polymer_entity_taxonomy_count?: Maybe<Scalars['Int']>;
  /** The maximum molecular mass (KDa) of a polymer entity in the deposited structure entry. */
  readonly polymer_molecular_weight_maximum?: Maybe<Scalars['Float']>;
  /** The minimum molecular mass (KDa) of a polymer entity in the deposited structure entry. */
  readonly polymer_molecular_weight_minimum?: Maybe<Scalars['Float']>;
  /** The maximum monomer count of a polymer entity per deposited structure model. */
  readonly polymer_monomer_count_maximum?: Maybe<Scalars['Int']>;
  /** The minimum monomer count of a polymer entity per deposited structure model. */
  readonly polymer_monomer_count_minimum?: Maybe<Scalars['Int']>;
  /** Combined estimates of experimental resolution. */
  readonly resolution_combined?: Maybe<ReadonlyArray<Maybe<Scalars['Float']>>>;
  /**
   * Selected polymer entity type categories describing the entry.
   * 
   * Allowable values:
   * Nucleic acid (only), Other, Protein (only), Protein/NA
   */
  readonly selected_polymer_entity_types?: Maybe<Scalars['String']>;
  /** Combined list of software programs names reported in connection with the production of this entry. */
  readonly software_programs_combined?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The number of distinct solvent entities per deposited structure model. */
  readonly solvent_entity_count?: Maybe<Scalars['Int']>;
};

export type PdbxChemCompFeature = {
  /**
   * The component identifier for this feature.
   * 
   * Examples:
   * ABC, ATP
   */
  readonly comp_id: Scalars['String'];
  /**
   * The information source for the component feature.
   * 
   * Examples:
   * PDB, CHEBI, DRUGBANK, PUBCHEM
   */
  readonly source: Scalars['String'];
  /**
   * The component feature type.
   * 
   * Examples:
   * FUNCTION, ENZYME INHIBITED, STRUCTURE IMAGE URL, CARBOHYDRATE ANOMER, CARBOHYDRATE ISOMER, CARBOHYDRATE RING
   */
  readonly type: Scalars['String'];
  /** The component feature value. */
  readonly value: Scalars['String'];
};

export type RcsbChemCompTarget = {
  /**
   * The value of _rcsb_chem_comp_target.comp_id is a reference to
   *  a chemical component definition.
   */
  readonly comp_id: Scalars['String'];
  /** The type of target interaction. */
  readonly interaction_type?: Maybe<Scalars['String']>;
  /** The target name. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _rcsb_chem_comp_target.ordinal distinguishes
   *  related examples for each chemical component.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * A code indicating the provenance of the target interaction assignment
   * 
   * Allowable values:
   * DrugBank, PDB Primary Data
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * The reference identifier code for the target interaction reference.
   * 
   * Examples:
   * Q9HD40
   */
  readonly reference_database_accession_code?: Maybe<Scalars['String']>;
  /**
   * The reference database name for the target interaction.
   * 
   * Allowable values:
   * UniProt
   */
  readonly reference_database_name?: Maybe<Scalars['String']>;
  /** The actions of the target interaction. */
  readonly target_actions?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type RcsbMembraneLineage = {
  /** Hierarchy depth. */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Identifier for membrane classification term (system generated for internal purposes identifier). */
  readonly id?: Maybe<Scalars['String']>;
  /** Membrane protein classification term. */
  readonly name?: Maybe<Scalars['String']>;
};

export type PdbxNmrRefine = {
  /**
   * Additional details about the NMR refinement.
   * 
   * Examples:
   * Additional comments about the NMR refinement can be placed here, e.g.
   * the structures are based on a total of 3344 restraints, 3167 are NOE-derived
   * distance constraints, 68 dihedral angle restraints,109 distance restraints 
   * from hydrogen bonds.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The method used to determine the structure.
   * 
   * Examples:
   * simulated annealing, distance geometry
   *   simulated annealing
   *   molecular dynamics
   *   matrix relaxation
   *   torsion angle dynamics
   */
  readonly method?: Maybe<Scalars['String']>;
  /** Pointer to _software.ordinal */
  readonly software_ordinal: Scalars['Int'];
};

export type RcsbUniprotProtein = {
  /** Enzyme Commission (EC) number(s). */
  readonly ec?: Maybe<ReadonlyArray<Maybe<RcsbUniprotProteinEc>>>;
  readonly function?: Maybe<RcsbUniprotProteinFunction>;
  /** The name(s) of the gene(s) that code for the protein sequence(s) described in the entry. */
  readonly gene?: Maybe<ReadonlyArray<Maybe<RcsbUniprotProteinGene>>>;
  readonly name?: Maybe<RcsbUniprotProteinName>;
  /** Protein sequence data for canonical protein sequence. */
  readonly sequence?: Maybe<Scalars['String']>;
  /** Taxonomy information on the organism that is the source of the protein sequence. */
  readonly source_organism?: Maybe<RcsbUniprotProteinSourceOrganism>;
};

export type ClustersMembers = {
  /** Internal chain ID used in mmCIF files to uniquely identify structural elements in the asymmetric unit. */
  readonly asym_id: Scalars['String'];
  readonly pdbx_struct_oper_list_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type RcsbPolymerEntityContainerIdentifiersReferenceSequenceIdentifiers = {
  /** Reference database accession code */
  readonly database_accession?: Maybe<Scalars['String']>;
  /** Reference database identifier for the sequence isoform */
  readonly database_isoform?: Maybe<Scalars['String']>;
  /**
   * Reference database name
   * 
   * Allowable values:
   * EMBL, GenBank, NDB, NORINE, PDB, PIR, PRF, RefSeq, UniProt
   */
  readonly database_name?: Maybe<Scalars['String']>;
  /**
   * Source of the reference database assignment
   * 
   * Allowable values:
   * PDB, RCSB, SIFTS
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
};

export type RcsbEntityHostOrganism = {
  /**
   * The beginning polymer sequence position for the polymer section corresponding
   *  to this host organism.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly beg_seq_num?: Maybe<Scalars['Int']>;
  /** The common name of the host organism */
  readonly common_name?: Maybe<Scalars['String']>;
  /**
   * The ending polymer sequence position for the polymer section corresponding
   *  to this host organism.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly end_seq_num?: Maybe<Scalars['Int']>;
  /**
   * Common names associated with this taxonomy code obtained from NCBI Taxonomy Database.
   * 
   *   These names correspond to the taxonomy identifier assigned by the PDB depositor.
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_common_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The parent scientific name in the NCBI taxonomy hierarchy (depth=1) associated with this taxonomy code.
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_parent_scientific_name?: Maybe<Scalars['String']>;
  /**
   * The scientific name associated with this taxonomy code aggregated by the NCBI Taxonomy Database.
   * 
   *   This name corresponds to the taxonomy identifier assigned by the PDB depositor.
   * 
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_scientific_name?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier for the host organism.
   * 
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly ncbi_taxonomy_id?: Maybe<Scalars['Int']>;
  /** An identifier for an entity segment. */
  readonly pdbx_src_id: Scalars['String'];
  /**
   * A code indicating the provenance of the host organism.
   * 
   * Allowable values:
   * PDB Primary Data
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /** The scientific name of the host organism */
  readonly scientific_name?: Maybe<Scalars['String']>;
  readonly taxonomy_lineage?: Maybe<ReadonlyArray<Maybe<RcsbEntityHostOrganismTaxonomyLineage>>>;
};

export type RcsbPolymerEntityAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type PdbxNmrDetails = {
  /**
   * Additional details describing the NMR experiment.
   * 
   * Examples:
   * This structure was determined using standard 2D homonuclear techniques., The structure was determined using triple-resonance NMR spectroscopy.
   */
  readonly text?: Maybe<Scalars['String']>;
};

export type PdbxReferenceEntitySequence = {
  /**
   * A flag to indicate a non-ribosomal entity.
   * 
   * Allowable values:
   * N, Y
   */
  readonly NRP_flag?: Maybe<Scalars['String']>;
  /** The one-letter-code sequence for this entity.  Non-standard monomers are represented as 'X'. */
  readonly one_letter_codes?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_sequence.prd_id is a reference
   * 	       _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_sequence.ref_entity_id is a reference
   *  to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly ref_entity_id: Scalars['String'];
  /**
   * The monomer type for the sequence.
   * 
   * Allowable values:
   * peptide-like, saccharide
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntityFeature = {
  /** Identifies the version of the feature assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** Non-polymer(ligand) chemical component identifier for the entity. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** A description for the feature. */
  readonly description?: Maybe<Scalars['String']>;
  /** An identifier for the feature. */
  readonly feature_id?: Maybe<Scalars['String']>;
  /** A name for the feature. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the feature.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the feature.
   * 
   * Allowable values:
   * SUBJECT_OF_INVESTIGATION
   */
  readonly type?: Maybe<Scalars['String']>;
  /** The feature value. */
  readonly value?: Maybe<Scalars['Float']>;
};

export type PdbxReferenceMoleculeFeatures = {
  /**
   * The value of _pdbx_reference_molecule_features.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_list.family_prd_id in category PDBX_REFERENCE_MOLECULE_FAMILY_LIST.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_molecule_features.ordinal distinguishes 
   * 	       each feature for this entity.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The value of _pdbx_reference_molecule_features.prd_id is a reference
   * 	       _pdbx_reference_molecule.prd_id in the  PDBX_REFERENCE_MOLECULE category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The information source for the component feature.
   * 
   * Examples:
   * PDB, CHEBI, DRUGBANK, PUBCHEM
   */
  readonly source?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_molecule_features.source_ordinal provides
   * 	       the priority order of features from a particular source or database.
   */
  readonly source_ordinal?: Maybe<Scalars['Int']>;
  /**
   * The entity feature type.
   * 
   * Examples:
   * FUNCTION, ENZYME INHIBITED, STRUCTURE IMAGE URL
   */
  readonly type?: Maybe<Scalars['String']>;
  /** The entity feature value. */
  readonly value?: Maybe<Scalars['String']>;
};

export type AuditAuthor = {
  /**
   * The Open Researcher and Contributor ID (ORCID).
   * 
   * Examples:
   * 0000-0002-6681-547X
   */
  readonly identifier_ORCID?: Maybe<Scalars['String']>;
  /**
   * The name of an author of this data block. If there are multiple
   *  authors, _audit_author.name is looped with _audit_author.address.
   *  The family name(s), followed by a comma and including any
   *  dynastic components, precedes the first name(s) or initial(s).
   * 
   * Examples:
   * Jones, T.J., Bleary, Percival R., O'Neil, F.K., Van den Bossche, G., Yang, D.-L., Simonov, Yu.A
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * This data item defines the order of the author's name in the
   *  list of audit authors.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly pdbx_ordinal: Scalars['Int'];
};

export type RcsbPolymerEntityAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * GO, InterPro, Pfam
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RefineAnalyze = {
  /**
   * The estimated coordinate error obtained from the plot of
   *  the R value versus sin(theta)/lambda for the reflections
   *  treated as a test set during refinement.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_coordinate_error_free?: Maybe<Scalars['Float']>;
  /**
   * The estimated coordinate error obtained from the plot of
   *  the R value versus sin(theta)/lambda for reflections classified
   *  as observed.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_coordinate_error_obs?: Maybe<Scalars['Float']>;
  /**
   * The value of the low-resolution cutoff used in constructing the
   *  Luzzati plot for reflections treated as a test set during
   *  refinement.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_d_res_low_free?: Maybe<Scalars['Float']>;
  /**
   * The value of the low-resolution cutoff used in
   *  constructing the Luzzati plot for reflections classified as
   *  observed.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_d_res_low_obs?: Maybe<Scalars['Float']>;
  /**
   * The value of sigma~a~ used in constructing the Luzzati plot for
   *  the reflections treated as a test set during refinement.
   *  Details of the estimation of sigma~a~ can be specified
   *  in _refine_analyze.Luzzati_sigma_a_free_details.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_sigma_a_free?: Maybe<Scalars['Float']>;
  /**
   * The value of sigma~a~ used in constructing the Luzzati plot for
   *  reflections classified as observed. Details of the
   *  estimation of sigma~a~ can be specified in
   *  _refine_analyze.Luzzati_sigma_a_obs_details.
   * 
   *  Ref:  Luzzati, V. (1952). Traitement statistique des erreurs
   *  dans la determination des structures cristallines. Acta
   *  Cryst. 5, 802-810.
   */
  readonly Luzzati_sigma_a_obs?: Maybe<Scalars['Float']>;
  /**
   * The number of discretely disordered residues in the refined
   *  model.
   */
  readonly number_disordered_residues?: Maybe<Scalars['Float']>;
  /**
   * The sum of the occupancies of the hydrogen atoms in the refined
   *  model.
   */
  readonly occupancy_sum_hydrogen?: Maybe<Scalars['Float']>;
  /**
   * The sum of the occupancies of the non-hydrogen atoms in the
   *   refined model.
   */
  readonly occupancy_sum_non_hydrogen?: Maybe<Scalars['Float']>;
  /** record the high resolution for calculating Luzzati statistics. */
  readonly pdbx_Luzzati_d_res_high_obs?: Maybe<Scalars['Float']>;
  /**
   * This data item uniquely identifies a refinement within an entry.
   *  _refine_analyze.pdbx_refine_id can be used to distinguish the results 
   *  of joint refinements.
   */
  readonly pdbx_refine_id: Scalars['String'];
};

export type PdbxEntityNonpoly = {
  /** This data item is a pointer to _chem_comp.id in the CHEM_COMP category. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** This data item is a pointer to _entity.id in the ENTITY category. */
  readonly entity_id: Scalars['String'];
  /** A name for the non-polymer entity */
  readonly name?: Maybe<Scalars['String']>;
  /** For non-polymer BIRD molecules the BIRD identifier for the entity. */
  readonly rcsb_prd_id?: Maybe<Scalars['String']>;
};

export type PdbxDatabaseStatus = {
  /**
   * This code indicates whether the entry belongs to
   *  Structural Genomics Project.
   * 
   * Allowable values:
   * N, Y
   */
  readonly SG_entry?: Maybe<Scalars['String']>;
  /**
   * The site where the file was deposited.
   * 
   * Allowable values:
   * BMRB, BNL, NDB, PDBC, PDBE, PDBJ, RCSB
   */
  readonly deposit_site?: Maybe<Scalars['String']>;
  /**
   * The methods development category in which this
   *  entry has been placed.
   * 
   * Allowable values:
   * CAPRI, CASD-NMR, CASP, D3R, FoldIt, GPCR Dock, RNA-Puzzles
   */
  readonly methods_development_category?: Maybe<Scalars['String']>;
  /**
   * A flag indicating that the entry is compatible with the PDB format.
   * 
   *  A value of 'N' indicates that the no PDB format data file is
   *  corresponding to this entry is available in the PDB archive.
   * 
   * Allowable values:
   * N, Y
   */
  readonly pdb_format_compatible?: Maybe<Scalars['String']>;
  /**
   * The site where the file was deposited.
   * 
   * Allowable values:
   * BNL, NDB, PDBC, PDBE, PDBJ, RCSB
   */
  readonly process_site?: Maybe<Scalars['String']>;
  /**
   * The date of initial deposition.  (The first message for
   *  deposition has been received.)
   * 
   * Examples:
   * 1983-02-21
   */
  readonly recvd_initial_deposition_date?: Maybe<Scalars['Date']>;
  /**
   * Code for status of file.
   * 
   * Allowable values:
   * AUCO, AUTH, BIB, DEL, HOLD, HPUB, OBS, POLC, PROC, REFI, REL, REPL, REV, RMVD, TRSF, UPD, WAIT, WDRN
   */
  readonly status_code?: Maybe<Scalars['String']>;
  /**
   * Code for status of chemical shift data file.
   * 
   * Allowable values:
   * AUTH, HOLD, HPUB, OBS, POLC, PROC, REL, REPL, RMVD, WAIT, WDRN
   */
  readonly status_code_cs?: Maybe<Scalars['String']>;
  /**
   * Code for status of NMR constraints file.
   * 
   * Allowable values:
   * AUTH, HOLD, HPUB, OBS, POLC, PROC, REL, REPL, RMVD, WAIT, WDRN
   */
  readonly status_code_mr?: Maybe<Scalars['String']>;
  /**
   * Code for status of structure factor file.
   * 
   * Allowable values:
   * AUTH, HOLD, HPUB, OBS, POLC, PROC, REL, REPL, RMVD, WAIT, WDRN
   */
  readonly status_code_sf?: Maybe<Scalars['String']>;
};

export type Em2dCrystalEntity = {
  /** Unit-cell angle gamma in degrees. */
  readonly angle_gamma?: Maybe<Scalars['Float']>;
  /** Length used to sample the reciprocal lattice lines in the c-direction. */
  readonly c_sampling_length?: Maybe<Scalars['Float']>;
  /** Unique key for the 2d_crystal_entity category. */
  readonly id: Scalars['String'];
  /** pointer to _em_image_processing.id in the EM_IMAGE_PROCESSING category. */
  readonly image_processing_id: Scalars['String'];
  /** Unit-cell length a in Angstroms. */
  readonly length_a?: Maybe<Scalars['Float']>;
  /** Unit-cell length b in Angstroms. */
  readonly length_b?: Maybe<Scalars['Float']>;
  /** Thickness of 2D crystal */
  readonly length_c?: Maybe<Scalars['Float']>;
  /**
   * There are 17 plane groups classified as oblique, rectangular, square, and hexagonal.
   *  To describe the symmetry of 2D crystals of biological molecules,
   *  plane groups are expanded to equivalent noncentrosymmetric space groups.
   *  The 2D crystal plane corresponds to the 'ab' plane of the space group.
   * 
   *  Enumerated space group descriptions include the plane group number in parentheses,
   *  the H-M plane group symbol, and the plane group class.
   * 
   * Allowable values:
   * C 1 2, C 2 2 2, P 1, P 1 2, P 1 21, P 2, P 2 2 2, P 2 2 21, P 2 21 21, P 3, P 3 1 2, P 3 2 1, P 4, P 4 2 2, P 4 21 2, P 6, P 6 2 2
   */
  readonly space_group_name_H_M?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntityAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** Non-polymer(ligand) chemical component identifier for the entity. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * SUBJECT_OF_INVESTIGATION
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerInstanceFeature = {
  /** Identifies the version of the feature assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** Component identifier for non-polymer entity instance. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** A description for the feature. */
  readonly description?: Maybe<Scalars['String']>;
  /** An identifier for the feature. */
  readonly feature_id?: Maybe<Scalars['String']>;
  readonly feature_value?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeatureFeatureValue>>>;
  /** A name for the feature. */
  readonly name?: Maybe<Scalars['String']>;
  /** Ordinal identifier for this category */
  readonly ordinal: Scalars['Int'];
  /**
   * Code identifying the individual, organization or program that
   *  assigned the feature.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the feature.
   * 
   * Allowable values:
   * HAS_COVALENT_LINKAGE, HAS_METAL_COORDINATION_LINKAGE, MOGUL_ANGLE_OUTLIER, MOGUL_BOND_OUTLIER, RSRCC_OUTLIER, RSRZ_OUTLIER
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerStructConnConnectTarget = {
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.auth_asym_id in the
   *  ATOM_SITE category.
   */
  readonly auth_asym_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.auth_seq_id in the
   *  ATOM_SITE category.
   */
  readonly auth_seq_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_alt_id in the
   *  ATOM_SITE category.
   */
  readonly label_alt_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_asym_id in the
   *  ATOM_SITE category.
   */
  readonly label_asym_id: Scalars['String'];
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_atom_id in the
   *  ATOM_SITE category.
   */
  readonly label_atom_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_comp_id in the
   *  ATOM_SITE category.
   */
  readonly label_comp_id: Scalars['String'];
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.connect_target_label_seq_id in the
   *  ATOM_SITE category.
   */
  readonly label_seq_id?: Maybe<Scalars['Int']>;
  /**
   * Describes the symmetry operation that should be applied to the
   *  atom set specified by _rcsb_nonpolymer_struct_conn.label* to generate the
   *  target of the structure connection.
   * 
   * Examples:
   * 1_555, 7_645
   */
  readonly symmetry?: Maybe<Scalars['String']>;
};

export type Em3dCrystalEntity = {
  /** Unit-cell angle alpha in degrees. */
  readonly angle_alpha?: Maybe<Scalars['Float']>;
  /** Unit-cell angle beta in degrees. */
  readonly angle_beta?: Maybe<Scalars['Float']>;
  /** Unit-cell angle gamma in degrees. */
  readonly angle_gamma?: Maybe<Scalars['Float']>;
  /** Unique key for the em_3d_crystal_entity category. */
  readonly id: Scalars['String'];
  /** pointer to _em_image_processing.id in the EM_IMAGE_PROCESSING category. */
  readonly image_processing_id: Scalars['String'];
  /** Unit-cell length a in Angstroms. */
  readonly length_a?: Maybe<Scalars['Float']>;
  /** Unit-cell length b in Angstroms. */
  readonly length_b?: Maybe<Scalars['Float']>;
  /** Unit-cell length c in Angstroms. */
  readonly length_c?: Maybe<Scalars['Float']>;
  /**
   * Space group name.
   * 
   * Examples:
   * P 1, P 21 21 2, I 4, H 3
   */
  readonly space_group_name?: Maybe<Scalars['String']>;
  /** Space group number. */
  readonly space_group_num?: Maybe<Scalars['Int']>;
};

export type GeneName = {
  /**
   * Allowable values:
   * PRIMARY, SYNONYM, ORDERED_LOCUS, ORF
   */
  readonly type?: Maybe<Scalars['String']>;
  readonly value?: Maybe<Scalars['String']>;
};

export type Citation = {
  /**
   * The International Standard Book Number (ISBN) code assigned to
   *  the book cited; relevant for books or book chapters.
   */
  readonly book_id_ISBN?: Maybe<Scalars['String']>;
  /**
   * The name of the publisher of the citation; relevant
   *  for books or book chapters.
   * 
   * Examples:
   * John Wiley and Sons
   */
  readonly book_publisher?: Maybe<Scalars['String']>;
  /**
   * The location of the publisher of the citation; relevant
   *  for books or book chapters.
   * 
   * Examples:
   * London
   */
  readonly book_publisher_city?: Maybe<Scalars['String']>;
  /**
   * The title of the book in which the citation appeared; relevant
   *  for books or book chapters.
   */
  readonly book_title?: Maybe<Scalars['String']>;
  /**
   * _citation.coordinate_linkage states whether this citation
   *  is concerned with precisely the set of coordinates given in the
   *  data block. If, for instance, the publication described the same
   *  structure, but the coordinates had undergone further refinement
   *  prior to the creation of the data block, the value of this data
   *  item would be 'no'.
   * 
   * Allowable values:
   * n, no, y, yes
   */
  readonly coordinate_linkage?: Maybe<Scalars['String']>;
  /**
   * The country/region of publication; relevant for books
   *  and book chapters.
   */
  readonly country?: Maybe<Scalars['String']>;
  /**
   * The value of _citation.id must uniquely identify a record in the
   *  CITATION list.
   * 
   *  The _citation.id 'primary' should be used to indicate the
   *  citation that the author(s) consider to be the most pertinent to
   *  the contents of the data block.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   * 
   * Examples:
   * primary, 1, 2
   */
  readonly id: Scalars['String'];
  /**
   * Abbreviated name of the cited journal as given in the
   *  Chemical Abstracts Service Source Index.
   * 
   * Examples:
   * J.Mol.Biol., J. Mol. Biol.
   */
  readonly journal_abbrev?: Maybe<Scalars['String']>;
  /**
   * The American Society for Testing and Materials (ASTM) code
   *  assigned to the journal cited (also referred to as the CODEN
   *  designator of the Chemical Abstracts Service); relevant for
   *  journal articles.
   */
  readonly journal_id_ASTM?: Maybe<Scalars['String']>;
  /**
   * The Cambridge Structural Database (CSD) code assigned to the
   *  journal cited; relevant for journal articles. This is also the
   *  system used at the Protein Data Bank (PDB).
   * 
   * Examples:
   * 0070
   */
  readonly journal_id_CSD?: Maybe<Scalars['String']>;
  /**
   * The International Standard Serial Number (ISSN) code assigned to
   *  the journal cited; relevant for journal articles.
   */
  readonly journal_id_ISSN?: Maybe<Scalars['String']>;
  /**
   * Issue number of the journal cited; relevant for journal
   *  articles.
   * 
   * Examples:
   * 2
   */
  readonly journal_issue?: Maybe<Scalars['String']>;
  /**
   * Volume number of the journal cited; relevant for journal
   *  articles.
   * 
   * Examples:
   * 174
   */
  readonly journal_volume?: Maybe<Scalars['String']>;
  /**
   * Language in which the cited article is written.
   * 
   * Examples:
   * German
   */
  readonly language?: Maybe<Scalars['String']>;
  /**
   * The first page of the citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_first?: Maybe<Scalars['String']>;
  /**
   * The last page of the citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_last?: Maybe<Scalars['String']>;
  /**
   * Document Object Identifier used by doi.org to uniquely
   *  specify bibliographic entry.
   * 
   * Examples:
   * 10.2345/S1384107697000225
   */
  readonly pdbx_database_id_DOI?: Maybe<Scalars['String']>;
  /**
   * Ascession number used by PubMed to categorize a specific
   *  bibliographic entry.
   * 
   * Examples:
   * 12627512
   */
  readonly pdbx_database_id_PubMed?: Maybe<Scalars['Int']>;
  /**
   * Names of the authors of the citation; relevant for journal
   *  articles, books and book chapters.  Names are separated by vertical bars.
   * 
   *  The family name(s), followed by a comma and including any
   *  dynastic components, precedes the first name(s) or initial(s).
   */
  readonly rcsb_authors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Flag to indicate a primary citation.
   * 
   * Allowable values:
   * N, Y
   */
  readonly rcsb_is_primary?: Maybe<Scalars['String']>;
  /**
   * Normalized journal abbreviation.
   * 
   * Examples:
   * Nat Struct Mol Biol
   */
  readonly rcsb_journal_abbrev?: Maybe<Scalars['String']>;
  /**
   * The title of the citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * Structure of diferric duck ovotransferrin
   *                                   at 2.35 Angstroms resolution.
   */
  readonly title?: Maybe<Scalars['String']>;
  /**
   * Flag to indicate that this citation will not be published.
   * 
   * Allowable values:
   * N, Y
   */
  readonly unpublished_flag?: Maybe<Scalars['String']>;
  /**
   * The year of the citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * 1984
   */
  readonly year?: Maybe<Scalars['Int']>;
};

export type RcsbEntitySourceOrganismTaxonomyLineage = {
  /** Members of the NCBI Taxonomy lineage as parent taxonomy lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the NCBI Taxonomy lineage as parent taxonomy idcodes. */
  readonly id?: Maybe<Scalars['String']>;
  /** Memebers of the NCBI Taxonomy lineage as parent taxonomy names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbPrimaryCitation = {
  /**
   * The International Standard Book Number (ISBN) code assigned to
   *  the book cited; relevant for books or book chapters.
   */
  readonly book_id_ISBN?: Maybe<Scalars['String']>;
  /**
   * The name of the publisher of the citation; relevant
   *  for books or book chapters.
   * 
   * Examples:
   * John Wiley and Sons
   */
  readonly book_publisher?: Maybe<Scalars['String']>;
  /**
   * The location of the publisher of the citation; relevant
   *  for books or book chapters.
   * 
   * Examples:
   * London
   */
  readonly book_publisher_city?: Maybe<Scalars['String']>;
  /**
   * The title of the book in which the citation appeared; relevant
   *  for books or book chapters.
   */
  readonly book_title?: Maybe<Scalars['String']>;
  /**
   * _rcsb_primary_citation.coordinate_linkage states whether this citation
   *  is concerned with precisely the set of coordinates given in the
   *  data block. If, for instance, the publication described the same
   *  structure, but the coordinates had undergone further refinement
   *  prior to the creation of the data block, the value of this data
   *  item would be 'no'.
   * 
   * Allowable values:
   * n, no, y, yes
   */
  readonly coordinate_linkage?: Maybe<Scalars['String']>;
  /**
   * The country/region of publication; relevant for books
   *  and book chapters.
   */
  readonly country?: Maybe<Scalars['String']>;
  /**
   * The value of _rcsb_primary_citation.id must uniquely identify a record in the
   *  CITATION list.
   * 
   *  The _rcsb_primary_citation.id 'primary' should be used to indicate the
   *  citation that the author(s) consider to be the most pertinent to
   *  the contents of the data block.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   * 
   * Examples:
   * primary
   */
  readonly id: Scalars['String'];
  /**
   * Abbreviated name of the cited journal as given in the
   *  Chemical Abstracts Service Source Index.
   * 
   * Examples:
   * J.Mol.Biol., J. Mol. Biol.
   */
  readonly journal_abbrev?: Maybe<Scalars['String']>;
  /**
   * The American Society for Testing and Materials (ASTM) code
   *  assigned to the journal cited (also referred to as the CODEN
   *  designator of the Chemical Abstracts Service); relevant for
   *  journal articles.
   */
  readonly journal_id_ASTM?: Maybe<Scalars['String']>;
  /**
   * The Cambridge Structural Database (CSD) code assigned to the
   *  journal cited; relevant for journal articles. This is also the
   *  system used at the Protein Data Bank (PDB).
   * 
   * Examples:
   * 0070
   */
  readonly journal_id_CSD?: Maybe<Scalars['String']>;
  /**
   * The International Standard Serial Number (ISSN) code assigned to
   *  the journal cited; relevant for journal articles.
   */
  readonly journal_id_ISSN?: Maybe<Scalars['String']>;
  /**
   * Issue number of the journal cited; relevant for journal
   *  articles.
   * 
   * Examples:
   * 2
   */
  readonly journal_issue?: Maybe<Scalars['String']>;
  /**
   * Volume number of the journal cited; relevant for journal
   *  articles.
   * 
   * Examples:
   * 174
   */
  readonly journal_volume?: Maybe<Scalars['String']>;
  /**
   * Language in which the cited article is written.
   * 
   * Examples:
   * German
   */
  readonly language?: Maybe<Scalars['String']>;
  /**
   * The first page of the citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_first?: Maybe<Scalars['String']>;
  /**
   * The last page of the citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_last?: Maybe<Scalars['String']>;
  /**
   * Document Object Identifier used by doi.org to uniquely
   *  specify bibliographic entry.
   * 
   * Examples:
   * 10.2345/S1384107697000225
   */
  readonly pdbx_database_id_DOI?: Maybe<Scalars['String']>;
  /**
   * Ascession number used by PubMed to categorize a specific
   *  bibliographic entry.
   * 
   * Examples:
   * 12627512
   */
  readonly pdbx_database_id_PubMed?: Maybe<Scalars['Int']>;
  /**
   * Names of the authors of the citation; relevant for journal
   *  articles, books and book chapters.  Names are separated by vertical bars.
   * 
   *  The family name(s), followed by a comma and including any
   *  dynastic components, precedes the first name(s) or initial(s).
   */
  readonly rcsb_authors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Normalized journal abbreviation.
   * 
   * Examples:
   * Nat Struct Mol Biol
   */
  readonly rcsb_journal_abbrev?: Maybe<Scalars['String']>;
  /**
   * The title of the citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * Structure of diferric duck ovotransferrin
   *                                   at 2.35 Angstroms resolution.
   */
  readonly title?: Maybe<Scalars['String']>;
  /**
   * The year of the citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * 1984
   */
  readonly year?: Maybe<Scalars['Int']>;
};

export type PdbxStructOperList = {
  /**
   * This identifier code must uniquely identify a
   *  record in the PDBX_STRUCT_OPER_LIST list.
   */
  readonly id: Scalars['String'];
  /**
   * The [1][1] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_1_1?: Maybe<Scalars['Float']>;
  /**
   * The [1][2] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_1_2?: Maybe<Scalars['Float']>;
  /**
   * The [1][3] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_1_3?: Maybe<Scalars['Float']>;
  /**
   * The [2][1] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_2_1?: Maybe<Scalars['Float']>;
  /**
   * The [2][2] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_2_2?: Maybe<Scalars['Float']>;
  /**
   * The [2][3] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_2_3?: Maybe<Scalars['Float']>;
  /**
   * The [3][1] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_3_1?: Maybe<Scalars['Float']>;
  /**
   * The [3][2] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_3_2?: Maybe<Scalars['Float']>;
  /**
   * The [3][3] element of the 3x3 matrix component of the
   *  transformation operation.
   */
  readonly matrix_3_3?: Maybe<Scalars['Float']>;
  /**
   * A descriptive name for the transformation operation.
   * 
   * Examples:
   * 1_555, two-fold rotation
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The symmetry operation corresponding to the transformation operation.
   * 
   * Examples:
   * x,y,z, x+1/2,y,-z
   */
  readonly symmetry_operation?: Maybe<Scalars['String']>;
  /**
   * A code to indicate the type of operator.
   * 
   * Allowable values:
   * 2D crystal symmetry operation, 3D crystal symmetry operation, build 2D crystal asymmetric unit, build 3D crystal asymmetric unit, build helical asymmetric unit, build point asymmetric unit, crystal symmetry operation, helical symmetry operation, identity operation, point symmetry operation, transform to 2D crystal frame, transform to 3D crystal frame, transform to crystal frame, transform to helical frame, transform to point frame
   */
  readonly type?: Maybe<Scalars['String']>;
  /**
   * The [1] element of the three-element vector component of the
   *  transformation operation.
   */
  readonly vector_1?: Maybe<Scalars['Float']>;
  /**
   * The [2] element of the three-element vector component of the
   *  transformation operation.
   */
  readonly vector_2?: Maybe<Scalars['Float']>;
  /**
   * The [3] element of the three-element vector component of the
   *  transformation operation.
   */
  readonly vector_3?: Maybe<Scalars['Float']>;
};

export type StructKeywords = {
  /**
   * Terms characterizing the macromolecular structure.
   * 
   * Examples:
   * DNA, RNA, T-RNA, DNA/RNA, RIBOZYME, PROTEIN/DNA, PROTEIN/RNA, PEPTIDE NUCLEIC ACID, PEPTIDE NUCLEIC ACID/DNA, DNA-BINDING PROTEIN, RNA-BINDING PROTEIN
   */
  readonly pdbx_keywords?: Maybe<Scalars['String']>;
  /**
   * Keywords describing this structure.
   * 
   * Examples:
   * Inhibitor, Complex, Isomerase..., serine protease, inhibited complex, high-resolution refinement
   */
  readonly text?: Maybe<Scalars['String']>;
};

export type RcsbUniprotProteinFunction = {
  /** General function(s) of a protein. */
  readonly details?: Maybe<Scalars['String']>;
  /** Historical record of the data attribute. */
  readonly provenance_code?: Maybe<Scalars['String']>;
};

export type RcsbPubmedContainerIdentifiers = {
  /** Unique integer value assigned to each PubMed record. */
  readonly pubmed_id?: Maybe<Scalars['Int']>;
};

export type RcsbPolymerInstanceAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /** Ordinal identifier for this category */
  readonly ordinal: Scalars['Int'];
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * CATH, SCOP
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxReferenceEntityPolySeq = {
  /**
   * A flag to indicate that sequence heterogeneity at this monomer position.
   * 
   * Allowable values:
   * N, Y
   */
  readonly hetero: Scalars['String'];
  /** This data item is the chemical component identifier of monomer. */
  readonly mon_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_poly_seq.num must uniquely and sequentially
   *  identify a record in the PDBX_REFERENCE_ENTITY_POLY_SEQ list.
   * 
   *  This value is conforms to author numbering conventions and does not map directly
   *  to the numbering conventions used for _entity_poly_seq.num.
   */
  readonly num: Scalars['Int'];
  /**
   * A flag to indicate that this monomer is observed in the instance example.
   * 
   * Allowable values:
   * N, Y
   */
  readonly observed?: Maybe<Scalars['String']>;
  /** This data item is the chemical component identifier for the parent component corresponding to this monomer. */
  readonly parent_mon_id?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_poly_seq.prd_id is a reference
   * 	       _pdbx_reference_entity_poly.prd_id in the  PDBX_REFERENCE_ENTITY_POLY category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_poly_seq.ref_entity_id is a reference
   *  to _pdbx_reference_entity_poly.ref_entity_id in PDBX_REFERENCE_ENTITY_POLY category.
   */
  readonly ref_entity_id: Scalars['String'];
};

export type RcsbEntityHostOrganismTaxonomyLineage = {
  /** Members of the NCBI Taxonomy lineage as parent taxonomy lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the NCBI Taxonomy lineage as parent taxonomy idcodes. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the NCBI Taxonomy lineage as parent taxonomy names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type CoreEntityAlignmentsAlignedRegions = {
  /** Aligned region length */
  readonly length: Scalars['Int'];
  /** Entity seqeunce start position */
  readonly query_begin: Scalars['Int'];
  /** NCBI sequence start position */
  readonly target_begin: Scalars['Int'];
};

export type RcsbUniprotKeyword = {
  /**
   * A unique keyword identifier.
   * 
   * Examples:
   * KW-0275, KW-0597
   */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * Human-readable keyword term.
   * 
   * Examples:
   * Lipid metabolism, Phosphoprotein, Fatty acid biosynthesis
   */
  readonly value?: Maybe<Scalars['String']>;
};

export type RcsbGenomicLineage = {
  /** Classification hierarchy depth. */
  readonly depth?: Maybe<Scalars['Int']>;
  /**
   * ID that identify taxonomy, chromosome or gene.
   * 
   * Examples:
   * 9606, 568815441, 414325
   */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * A human-readable term name.
   * 
   * Examples:
   * Homo sapiens, 8, defensin beta 103A
   */
  readonly name?: Maybe<Scalars['String']>;
};

export type DiffrnDetector = {
  /** A description of special aspects of the radiation detector. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The general class of the radiation detector.
   * 
   * Examples:
   * photographic film, scintillation counter, CCD plate, BF~3~ counter
   */
  readonly detector?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * The date of data collection.
   * 
   * Examples:
   * 1996-12-25
   */
  readonly pdbx_collection_date?: Maybe<Scalars['Date']>;
  /** The operating frequency of the detector (Hz) used in data collection. */
  readonly pdbx_frequency?: Maybe<Scalars['Float']>;
  /** The make, model or name of the detector device used. */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbSchemaContainerIdentifiers = {
  /** Collection name associated with the data in the container. */
  readonly collection_name: Scalars['String'];
  /** Version string for the schema and collection. */
  readonly collection_schema_version?: Maybe<Scalars['String']>;
  /** Schema name associated with the data in the container. */
  readonly schema_name: Scalars['String'];
};

export type CorePfam = {
  /** Accession number of Pfam entry. */
  readonly rcsb_id: Scalars['String'];
  /**
   * The unique accession code of protein families and domains in the Pfam database.
   * 
   * Examples:
   * PF00621, PF00637, PF00656
   */
  readonly rcsb_pfam_accession: Scalars['String'];
  /** Details of the Pfam clan to which the entity belongs. */
  readonly rcsb_pfam_clan_id?: Maybe<Scalars['String']>;
  /** Textual description of the family. */
  readonly rcsb_pfam_comment?: Maybe<Scalars['String']>;
  readonly rcsb_pfam_container_identifiers: RcsbPfamContainerIdentifiers;
  /**
   * A human-readable name of protein families and domains.
   * 
   * Examples:
   * Lectin like domain, Cell division control protein 24, OB domain 2, Protein of unknown function (DUF722)
   */
  readonly rcsb_pfam_description?: Maybe<Scalars['String']>;
  /**
   * The unique identifier of protein families and domains in the Pfam database.
   * 
   * Examples:
   * RhoGEF, Clathrin, Peptidase_C14
   */
  readonly rcsb_pfam_identifier?: Maybe<Scalars['String']>;
  /**
   * Pfam-A is the manually curated portion of the Pfam database.
   * 
   * Allowable values:
   * Pfam-A
   */
  readonly rcsb_pfam_provenance_code?: Maybe<Scalars['String']>;
  /**
   * Pfam entries are classified into six different categories, depending on the length and nature of the sequence regions included in the entry: family, domain, repeats, motifs, coiled-coil, and disordered.
   * 
   * Allowable values:
   * Family, Domain, Repeat, Motif, Disordered, Coiled-coil
   */
  readonly rcsb_pfam_seed_source?: Maybe<Scalars['String']>;
};

export type EmImageRecording = {
  /**
   * The average exposure time for each image.
   * 
   * Examples:
   * 2.0
   */
  readonly average_exposure_time?: Maybe<Scalars['Float']>;
  /**
   * The electron dose received by the specimen per image (electrons per square angstrom).
   * 
   * Examples:
   * 30.0
   */
  readonly avg_electron_dose_per_image?: Maybe<Scalars['Float']>;
  /**
   * Any additional details about image recording.
   * 
   * Examples:
   * Images were collected in movie-mode at 17 frames per second
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The detector mode used during image recording.
   * 
   * Allowable values:
   * COUNTING, INTEGRATING, OTHER, SUPER-RESOLUTION
   */
  readonly detector_mode?: Maybe<Scalars['String']>;
  /**
   * The detector type used for recording images.
   *  Usually film or CCD camera.
   */
  readonly film_or_detector_model?: Maybe<Scalars['String']>;
  /**
   * The item _em_image_recording.id uniquely identifies
   *  a set of recorded images.
   */
  readonly id: Scalars['String'];
  /** This data item the id of the microscopy settings used in the imaging. */
  readonly imaging_id: Scalars['String'];
  /** The number of diffraction images collected. */
  readonly num_diffraction_images?: Maybe<Scalars['Int']>;
  /** Number of grids in the microscopy session */
  readonly num_grids_imaged?: Maybe<Scalars['Int']>;
  /** The number of micrograph images collected. */
  readonly num_real_images?: Maybe<Scalars['Int']>;
};

export type RcsbBirdCitation = {
  /**
   * The value of _rcsb_bird_citation.id must uniquely identify a record in the
   *  rcsb_bird_citation list.
   * 
   * Examples:
   * 1, 2
   */
  readonly id: Scalars['String'];
  /**
   * Abbreviated name of the cited journal as given in the
   *  Chemical Abstracts Service Source Index.
   * 
   * Examples:
   * J.Mol.Biol., J. Mol. Biol.
   */
  readonly journal_abbrev?: Maybe<Scalars['String']>;
  /**
   * Volume number of the journal cited; relevant for journal
   *  articles.
   * 
   * Examples:
   * 174
   */
  readonly journal_volume?: Maybe<Scalars['String']>;
  /**
   * The first page of the rcsb_bird_citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_first?: Maybe<Scalars['String']>;
  /**
   * The last page of the rcsb_bird_citation; relevant for journal
   *  articles, books and book chapters.
   */
  readonly page_last?: Maybe<Scalars['String']>;
  /**
   * Document Object Identifier used by doi.org to uniquely
   *  specify bibliographic entry.
   * 
   * Examples:
   * 10.2345/S1384107697000225
   */
  readonly pdbx_database_id_DOI?: Maybe<Scalars['String']>;
  /**
   * Ascession number used by PubMed to categorize a specific
   *  bibliographic entry.
   * 
   * Examples:
   * 12627512
   */
  readonly pdbx_database_id_PubMed?: Maybe<Scalars['Int']>;
  /**
   * Names of the authors of the citation; relevant for journal
   *  articles, books and book chapters.  Names are separated by vertical bars.
   * 
   *  The family name(s), followed by a comma and including any
   *  dynastic components, precedes the first name(s) or initial(s).
   */
  readonly rcsb_authors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The title of the rcsb_bird_citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * Structure of diferric duck ovotransferrin
   *                                   at 2.35 Angstroms resolution.
   */
  readonly title?: Maybe<Scalars['String']>;
  /**
   * The year of the rcsb_bird_citation; relevant for journal articles, books
   *  and book chapters.
   * 
   * Examples:
   * 1984
   */
  readonly year?: Maybe<Scalars['Int']>;
};

export type RcsbNonpolymerInstanceAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbUniprotProteinName = {
  /** Historical record of the data attribute. */
  readonly provenance_code: Scalars['String'];
  /** Name that allows to unambiguously identify a protein. */
  readonly value: Scalars['String'];
};

export type EmEmbedding = {
  /**
   * Staining procedure used in the specimen preparation.
   * 
   * Examples:
   * The crystal suspension was injected into the lens of a drop of buffer containing
   *   1 % tannin sitting on a carbon film supported by a molybdenum grid.  An equal volume
   *   of 1% glucose was then added and the solution thoroughly but gently mixed.  The grid
   *   was then blotted, air dried, and frozen in LN2.
   */
  readonly details?: Maybe<Scalars['String']>;
  /** This data item is the primary key of the category. */
  readonly id: Scalars['String'];
  /**
   * The embedding  material.
   * 
   * Examples:
   * tannin and glucose
   */
  readonly material?: Maybe<Scalars['String']>;
  /** Foreign key relationship to the EMD SPECIMEN category */
  readonly specimen_id?: Maybe<Scalars['String']>;
};

export type PdbxReferenceEntitySrcNat = {
  /** The Americal Tissue Culture Collection code for organism from which the entity was isolated. */
  readonly atcc?: Maybe<Scalars['String']>;
  /** The database code for this source information */
  readonly db_code?: Maybe<Scalars['String']>;
  /** The database name for this source information */
  readonly db_name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_src_nat.ordinal distinguishes 
   * 	       source details for this entity.
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The scientific name of the organism from which the entity was isolated.
   * 
   * Examples:
   * Mus musculus
   */
  readonly organism_scientific?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_src_nat.prd_id is a reference
   * 	       _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_src_nat.ref_entity_id is a reference
   *  to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly ref_entity_id: Scalars['String'];
  /** The data source for this information. */
  readonly source?: Maybe<Scalars['String']>;
  /** A identifier within the data source for this information. */
  readonly source_id?: Maybe<Scalars['String']>;
  /** The NCBI TaxId of the organism from which the entity was isolated. */
  readonly taxid?: Maybe<Scalars['String']>;
};

export type EmDiffractionStats = {
  /**
   * Any addition details about the structure factor measurements
   * 
   * Examples:
   * Phases were obtained from micrograph images of the 2D crystals
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Completeness of the structure factor data within the defined space group
   *  at the reported resolution (percent).
   * 
   * Examples:
   * 89.3
   */
  readonly fourier_space_coverage?: Maybe<Scalars['Float']>;
  /**
   * High resolution limit of the structure factor data, in Angstroms
   * 
   * Examples:
   * 7.5
   */
  readonly high_resolution?: Maybe<Scalars['Float']>;
  /** Identifier for this category */
  readonly id: Scalars['String'];
  /** Pointer to _em_image_processing.id */
  readonly image_processing_id?: Maybe<Scalars['String']>;
  /**
   * Total number of diffraction intensities measured (before averaging)
   * 
   * Examples:
   * 1590
   */
  readonly num_intensities_measured?: Maybe<Scalars['Int']>;
  /**
   * Number of structure factors obtained (merged amplitudes + phases)
   * 
   * Examples:
   * 325
   */
  readonly num_structure_factors?: Maybe<Scalars['Int']>;
  /**
   * Overall phase error in degrees
   * 
   * Examples:
   * 17.5
   */
  readonly overall_phase_error?: Maybe<Scalars['Float']>;
  /**
   * Overall phase residual in degrees
   * 
   * Examples:
   * 17.5
   */
  readonly overall_phase_residual?: Maybe<Scalars['Float']>;
  /**
   * Criteria used to reject phases
   * 
   * Examples:
   * Structure factors with phase errors higher than 20 degrees were omitted from refinement
   */
  readonly phase_error_rejection_criteria?: Maybe<Scalars['String']>;
  /**
   * Rmerge value (percent)
   * 
   * Examples:
   * 19.8
   */
  readonly r_merge?: Maybe<Scalars['Float']>;
  /**
   * Rsym value (percent)
   * 
   * Examples:
   * 24.4
   */
  readonly r_sym?: Maybe<Scalars['Float']>;
};

export type PdbxAuditRevisionItem = {
  /**
   * The type of file that the pdbx_audit_revision_history record refers to.
   * 
   * Allowable values:
   * Chemical component, NMR restraints, NMR shifts, Structure factors, Structure model
   */
  readonly data_content_type: Scalars['String'];
  /**
   * A high level explanation the author has provided for submitting a revision.
   * 
   * Examples:
   * _atom_site.type_symbol
   */
  readonly item?: Maybe<Scalars['String']>;
  /**
   * A unique identifier for the pdbx_audit_revision_item record.
   * 
   * Examples:
   * 1
   */
  readonly ordinal: Scalars['Int'];
  /**
   * A pointer to  _pdbx_audit_revision_history.ordinal
   * 
   * Examples:
   * 1
   */
  readonly revision_ordinal: Scalars['Int'];
};

export type PdbxSerialCrystallographySampleDelivery = {
  /**
   * The description of the mechanism by which the specimen in placed in the path
   *  of the source.
   * 
   * Examples:
   * fixed target, electrospin, MESH, CoMESH, gas dynamic virtual nozzle, LCP injector, addressable microarray
   */
  readonly description?: Maybe<Scalars['String']>;
  /**
   * The data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   * 
   * Examples:
   * 1
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * The description of the mechanism by which the specimen in placed in the path
   *  of the source.
   * 
   * Allowable values:
   * fixed target, injection
   */
  readonly method?: Maybe<Scalars['String']>;
};

export type Em3dReconstruction = {
  /**
   * The actual pixel size of projection set of images.
   * 
   * Examples:
   * 2.8, 5.76
   */
  readonly actual_pixel_size?: Maybe<Scalars['Float']>;
  /** The algorithm used project from 2D orientations to 3D map. */
  readonly algorithm?: Maybe<Scalars['String']>;
  /**
   * Any additional details used in the 3d reconstruction.
   * 
   * Examples:
   * a modified version of SPIDER program was used for the reconstruction
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The value of _em_3d_reconstruction.id must
   *  uniquely identify the 3d reconstruction.
   */
  readonly id: Scalars['String'];
  /** Foreign key to the EM_IMAGE_PROCESSING category */
  readonly image_processing_id: Scalars['String'];
  /**
   * The magnification calibration method for the 3d reconstruction.
   * 
   * Examples:
   * TMV images
   */
  readonly magnification_calibration?: Maybe<Scalars['String']>;
  /**
   * The algorithm method used for the 3d-reconstruction.
   * 
   * Examples:
   * cross-common lines, polar Fourier transform (PFT)
   */
  readonly method?: Maybe<Scalars['String']>;
  /**
   * The nominal pixel size of the projection set of images.
   * 
   * Examples:
   * 3.11, 6.78
   */
  readonly nominal_pixel_size?: Maybe<Scalars['Float']>;
  /**
   * This item was correspondence to two type of em dataset
   *  processing_emDataSet_singleParticle.numClassAverages
   *  processing_emDataSet_icosahedral.numClassAverages
   */
  readonly num_class_averages?: Maybe<Scalars['Int']>;
  /** The number of 2D projections or 3D subtomograms used in the 3d reconstruction */
  readonly num_particles?: Maybe<Scalars['Int']>;
  /**
   * type of refinement performed in order to determine map resolution
   * 
   * Allowable values:
   * HALF-MAPS REFINED AGAINST SAME DATA, HALF-MAPS REFINED INDEPENDENTLY, HALF-MAPS REFINED INDEPENDENTLY WITH FREQUENCY RANGE OMITTED, HALF-MAPS REFINED WITH FREQUENCY RANGE OMITTED, OTHER
   */
  readonly refinement_type?: Maybe<Scalars['String']>;
  /**
   * The final resolution (in Angstroms)of the 3D reconstruction.
   * 
   * Examples:
   * 8.9, 10.0
   */
  readonly resolution?: Maybe<Scalars['Float']>;
  /**
   * The  method used to determine the final resolution
   *  of the 3d reconstruction.
   *  The Fourier Shell Correlation criterion as a measure of
   *  resolution is based on the concept of splitting the (2D)
   *  data set into two halves; averaging each and comparing them
   *  using the Fourier Ring Correlation (FRC) technique.
   * 
   * Examples:
   * FSC at 0.5 cut-off
   */
  readonly resolution_method?: Maybe<Scalars['String']>;
  /**
   * The type of symmetry applied to the reconstruction
   * 
   * Allowable values:
   * 2D CRYSTAL, 3D CRYSTAL, HELICAL, POINT
   */
  readonly symmetry_type?: Maybe<Scalars['String']>;
};

export type EmParticleSelection = {
  /**
   * Any additional details used for selecting particles
   * 
   * Examples:
   * negative monitor contrast facilitated particle picking
   */
  readonly details?: Maybe<Scalars['String']>;
  /** Ordinal identifier */
  readonly id: Scalars['String'];
  /**
   * The value of _em_particle_selection.image_processing_id points to
   *  the EM_IMAGE_PROCESSING category.
   */
  readonly image_processing_id: Scalars['String'];
  /**
   * The number of particles selected from the projection set of images.
   * 
   * Examples:
   * 840
   */
  readonly num_particles_selected?: Maybe<Scalars['Int']>;
};

export type CoreChemComp = {
  readonly chem_comp?: Maybe<ChemComp>;
  /** Get DrubBank entry associated with this chemical component. */
  readonly drugbank?: Maybe<CoreDrugbank>;
  readonly pdbx_chem_comp_audit?: Maybe<ReadonlyArray<Maybe<PdbxChemCompAudit>>>;
  readonly pdbx_chem_comp_descriptor?: Maybe<ReadonlyArray<Maybe<PdbxChemCompDescriptor>>>;
  readonly pdbx_chem_comp_feature?: Maybe<ReadonlyArray<Maybe<PdbxChemCompFeature>>>;
  readonly pdbx_chem_comp_identifier?: Maybe<ReadonlyArray<Maybe<PdbxChemCompIdentifier>>>;
  readonly pdbx_family_prd_audit?: Maybe<ReadonlyArray<Maybe<PdbxFamilyPrdAudit>>>;
  readonly pdbx_prd_audit?: Maybe<ReadonlyArray<Maybe<PdbxPrdAudit>>>;
  readonly pdbx_reference_entity_list?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityList>>>;
  readonly pdbx_reference_entity_poly?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPoly>>>;
  readonly pdbx_reference_entity_poly_link?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPolyLink>>>;
  readonly pdbx_reference_entity_poly_seq?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPolySeq>>>;
  readonly pdbx_reference_entity_sequence?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntitySequence>>>;
  readonly pdbx_reference_entity_src_nat?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntitySrcNat>>>;
  readonly pdbx_reference_molecule?: Maybe<PdbxReferenceMolecule>;
  readonly pdbx_reference_molecule_annotation?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeAnnotation>>>;
  readonly pdbx_reference_molecule_details?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeDetails>>>;
  readonly pdbx_reference_molecule_family?: Maybe<PdbxReferenceMoleculeFamily>;
  readonly pdbx_reference_molecule_features?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeFeatures>>>;
  readonly pdbx_reference_molecule_list?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeList>>>;
  readonly pdbx_reference_molecule_related_structures?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeRelatedStructures>>>;
  readonly pdbx_reference_molecule_synonyms?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeSynonyms>>>;
  readonly rcsb_bird_citation?: Maybe<ReadonlyArray<Maybe<RcsbBirdCitation>>>;
  readonly rcsb_chem_comp_annotation?: Maybe<ReadonlyArray<Maybe<RcsbChemCompAnnotation>>>;
  readonly rcsb_chem_comp_container_identifiers?: Maybe<RcsbChemCompContainerIdentifiers>;
  readonly rcsb_chem_comp_descriptor?: Maybe<RcsbChemCompDescriptor>;
  readonly rcsb_chem_comp_info?: Maybe<RcsbChemCompInfo>;
  readonly rcsb_chem_comp_related?: Maybe<ReadonlyArray<Maybe<RcsbChemCompRelated>>>;
  readonly rcsb_chem_comp_synonyms?: Maybe<ReadonlyArray<Maybe<RcsbChemCompSynonyms>>>;
  readonly rcsb_chem_comp_target?: Maybe<ReadonlyArray<Maybe<RcsbChemCompTarget>>>;
  /** A unique identifier for the chemical definition in this container. */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_schema_container_identifiers?: Maybe<ReadonlyArray<Maybe<RcsbSchemaContainerIdentifiers>>>;
};

export type EmExperiment = {
  /**
   * The aggregation/assembly state of the imaged specimen.
   * 
   * Allowable values:
   * 2D ARRAY, 3D ARRAY, CELL, FILAMENT, HELICAL ARRAY, PARTICLE, TISSUE
   */
  readonly aggregation_state?: Maybe<Scalars['String']>;
  /** Foreign key to the EM_ENTITY_ASSEMBLY category */
  readonly entity_assembly_id?: Maybe<Scalars['String']>;
  /** Placeholder ID. */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * The reconstruction method used in the EM experiment.
   * 
   * Allowable values:
   * CRYSTALLOGRAPHY, HELICAL, SINGLE PARTICLE, SUBTOMOGRAM AVERAGING, TOMOGRAPHY
   */
  readonly reconstruction_method?: Maybe<Scalars['String']>;
};

export type PdbxStructAssemblyProp = {
  /** The identifier for the assembly used in category PDBX_STRUCT_ASSEMBLY. */
  readonly assembly_id?: Maybe<Scalars['String']>;
  /** The identifier for the assembly used in category PDBX_STRUCT_ASSEMBLY. */
  readonly biol_id: Scalars['String'];
  /**
   * The property type for the assembly.
   * 
   * Allowable values:
   * ABSA (A^2), MORE, SSA (A^2)
   */
  readonly type: Scalars['String'];
  /** The value of the assembly property. */
  readonly value?: Maybe<Scalars['String']>;
};

export type PdbxNmrSoftware = {
  /**
   * The name of the authors of the software used in this
   *  procedure.
   * 
   * Examples:
   * Brunger, Guentert
   */
  readonly authors?: Maybe<Scalars['String']>;
  /**
   * The purpose of the software.
   * 
   * Examples:
   * collection, processing, data analysis, structure solution, refinement, iterative matrix relaxation
   */
  readonly classification?: Maybe<Scalars['String']>;
  /**
   * The name of the software used for the task.
   * 
   * Examples:
   * ANSIG, AURELIA, AZARA, CHARMM, CoMAND, CORMA, DIANA, DYANA, DSPACE, DISGEO, DGII, DISMAN, DINOSAUR, DISCOVER, FELIX, FT_NMR, GROMOS, IRMA, MARDIGRAS, NMRPipe, SA, UXNMR, VNMR, X-PLOR, XWINNMR
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * An ordinal index for this category
   * 
   * Examples:
   * 1, 2
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The version of the software.
   * 
   * Examples:
   * 940501.3, 2.1
   */
  readonly version?: Maybe<Scalars['String']>;
};

export type RcsbUniprotAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type Diffrn = {
  /**
   * The mean hydrostatic pressure in kilopascals at which the
   *  intensities were measured.
   */
  readonly ambient_pressure?: Maybe<Scalars['Float']>;
  /**
   * The mean temperature in kelvins at which the intensities were
   *  measured.
   */
  readonly ambient_temp?: Maybe<Scalars['Float']>;
  /**
   * A description of special aspects of temperature control during
   *  data collection.
   */
  readonly ambient_temp_details?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _exptl_crystal.id in the
   *  EXPTL_CRYSTAL category.
   */
  readonly crystal_id?: Maybe<Scalars['String']>;
  /**
   * The physical device used to support the crystal during data
   *  collection.
   * 
   * Examples:
   * glass capillary, quartz capillary, fiber, metal loop
   */
  readonly crystal_support?: Maybe<Scalars['String']>;
  /**
   * Special details of the diffraction measurement process. Should
   *  include information about source instability, crystal motion,
   *  degradation and so on.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * This data item uniquely identifies a set of diffraction
   *  data.
   */
  readonly id: Scalars['String'];
  /**
   * Y/N if using serial crystallography experiment in which multiple crystals contribute to each diffraction frame in the experiment.
   * 
   * Examples:
   * Y, N
   */
  readonly pdbx_serial_crystal_experiment?: Maybe<Scalars['String']>;
};

export type PdbxNmrExptlSampleConditions = {
  /**
   * The condition number as defined above.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly conditions_id: Scalars['String'];
  /**
   * General details describing conditions of both the sample and the environment 
   * during measurements.
   * 
   * Examples:
   * The high salinity of the sample may have contributed to overheating of the sample during experiments with long saturation periods like the TOCSY experiments.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The ionic strength at which the NMR data were collected -in lieu of
   *  this enter the concentration and identity of the salt in the sample.
   */
  readonly ionic_strength?: Maybe<Scalars['String']>;
  /**
   * Estimate of the standard error for the value for the sample ionic strength.
   * 
   * Examples:
   * 0.2
   */
  readonly ionic_strength_err?: Maybe<Scalars['Float']>;
  /**
   * Units for the value of the sample condition ionic strength..
   * 
   * Allowable values:
   * M, Not defined, mM
   */
  readonly ionic_strength_units?: Maybe<Scalars['String']>;
  /**
   * A descriptive label that uniquely identifies this set of sample conditions.
   * 
   * Examples:
   * conditions_1
   */
  readonly label?: Maybe<Scalars['String']>;
  /**
   * The pH at which the NMR data were collected.
   * 
   * Examples:
   * 3.1, 7.0
   */
  readonly pH?: Maybe<Scalars['String']>;
  /**
   * Estimate of the standard error for the value for the sample pH.
   * 
   * Examples:
   * 0.05
   */
  readonly pH_err?: Maybe<Scalars['Float']>;
  /**
   * Units for the value of the sample condition pH.
   * 
   * Allowable values:
   * Not defined, pD, pH, pH*
   */
  readonly pH_units?: Maybe<Scalars['String']>;
  /**
   * The pressure at which NMR data were collected.
   * 
   * Examples:
   * 1, ambient, 1atm
   */
  readonly pressure?: Maybe<Scalars['String']>;
  /**
   * Estimate of the standard error for the value for the sample pressure.
   * 
   * Examples:
   * 0.01
   */
  readonly pressure_err?: Maybe<Scalars['Float']>;
  /**
   * The units of pressure at which NMR data were collected.
   * 
   * Examples:
   * Pa, atm, Torr
   */
  readonly pressure_units?: Maybe<Scalars['String']>;
  /**
   * The temperature (in Kelvin) at which NMR data were
   *  collected.
   * 
   * Examples:
   * 298
   */
  readonly temperature?: Maybe<Scalars['String']>;
  /**
   * Estimate of the standard error for the value for the sample temperature.
   * 
   * Examples:
   * 0.2
   */
  readonly temperature_err?: Maybe<Scalars['Float']>;
  /**
   * Units for the value of the sample condition temperature.
   * 
   * Allowable values:
   * C, K, Not defined
   */
  readonly temperature_units?: Maybe<Scalars['String']>;
};

export type Refine = {
  /**
   * The maximum isotropic displacement parameter (B value)
   *  found in the coordinate set.
   */
  readonly B_iso_max?: Maybe<Scalars['Float']>;
  /**
   * The mean isotropic displacement parameter (B value)
   *  for the coordinate set.
   */
  readonly B_iso_mean?: Maybe<Scalars['Float']>;
  /**
   * The minimum isotropic displacement parameter (B value)
   *  found in the coordinate set.
   */
  readonly B_iso_min?: Maybe<Scalars['Float']>;
  /**
   * The [1][1] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_1_1?: Maybe<Scalars['Float']>;
  /**
   * The [1][2] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_1_2?: Maybe<Scalars['Float']>;
  /**
   * The [1][3] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_1_3?: Maybe<Scalars['Float']>;
  /**
   * The [2][2] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_2_2?: Maybe<Scalars['Float']>;
  /**
   * The [2][3] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_2_3?: Maybe<Scalars['Float']>;
  /**
   * The [3][3] element of the matrix that defines the overall
   *  anisotropic displacement model if one was refined for this
   *  structure.
   */
  readonly aniso_B_3_3?: Maybe<Scalars['Float']>;
  /**
   * The correlation coefficient between the observed and
   *              calculated structure factors for reflections included in
   *              the refinement.
   * 
   *              The correlation coefficient is scale-independent and gives
   *              an idea of the quality of the refined model.
   * 
   *                           sum~i~(Fo~i~ Fc~i~ - <Fo><Fc>)
   * R~corr~ = ------------------------------------------------------------
   *           SQRT{sum~i~(Fo~i~)^2^-<Fo>^2^} SQRT{sum~i~(Fc~i~)^2^-<Fc>^2^}
   * 
   *              Fo = observed structure factors
   *              Fc = calculated structure factors
   *              <>   denotes average value
   * 
   *              summation is over reflections included in the refinement
   */
  readonly correlation_coeff_Fo_to_Fc?: Maybe<Scalars['Float']>;
  /**
   * The correlation coefficient between the observed and
   *              calculated structure factors for reflections not included
   *              in the refinement (free reflections).
   * 
   *               The correlation coefficient is scale-independent and gives
   *               an idea of the quality of the refined model.
   * 
   *                           sum~i~(Fo~i~ Fc~i~ - <Fo><Fc>)
   * R~corr~ = ------------------------------------------------------------
   *           SQRT{sum~i~(Fo~i~)^2^-<Fo>^2^} SQRT{sum~i~(Fc~i~)^2^-<Fc>^2^}
   * 
   *               Fo  = observed structure factors
   *               Fc  = calculated structure factors
   *               <>    denotes average value
   * 
   *               summation is over reflections not included
   *               in the refinement (free reflections)
   */
  readonly correlation_coeff_Fo_to_Fc_free?: Maybe<Scalars['Float']>;
  /** Description of special aspects of the refinement process. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Residual factor R for reflections that satisfy the resolution
   *  limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the test
   *  reflections (i.e. were excluded from the refinement) when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   * 
   *      sum|F~obs~ - F~calc~|
   *  R = ---------------------
   *           sum|F~obs~|
   * 
   *  F~obs~  = the observed structure-factor amplitudes
   *  F~calc~ = the calculated structure-factor amplitudes
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_R_factor_R_free?: Maybe<Scalars['Float']>;
  /**
   * The estimated error in _refine.ls_R_factor_R_free.
   *  The method used to estimate the error is described in the
   *  item _refine.ls_R_factor_R_free_error_details.
   */
  readonly ls_R_factor_R_free_error?: Maybe<Scalars['Float']>;
  /**
   * Special aspects of the method used to estimated the error in
   *  _refine.ls_R_factor_R_free.
   */
  readonly ls_R_factor_R_free_error_details?: Maybe<Scalars['String']>;
  /**
   * Residual factor R for reflections that satisfy the resolution
   *  limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the working
   *  reflections (i.e. were included in the refinement)  when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   * 
   *  _refine.ls_R_factor_obs should not be confused with
   *  _refine.ls_R_factor_R_work; the former reports the results of a
   *  refinement in which all observed reflections were used, the
   *  latter a refinement in which a subset of the observed
   *  reflections were excluded from refinement for the calculation
   *  of a 'free' R factor. However, it would be meaningful to quote
   *  both values if a 'free' R factor were calculated for most of
   *  the refinement, but all of the observed reflections were used
   *  in the final rounds of refinement; such a protocol should be
   *  explained in _refine.details.
   * 
   *      sum|F~obs~ - F~calc~|
   *  R = ---------------------
   *           sum|F~obs~|
   * 
   *  F~obs~  = the observed structure-factor amplitudes
   *  F~calc~ = the calculated structure-factor amplitudes
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_R_factor_R_work?: Maybe<Scalars['Float']>;
  /**
   * Residual factor R for all reflections that satisfy the resolution
   *  limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low.
   * 
   *      sum|F~obs~ - F~calc~|
   *  R = ---------------------
   *           sum|F~obs~|
   * 
   *  F~obs~  = the observed structure-factor amplitudes
   *  F~calc~ = the calculated structure-factor amplitudes
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_R_factor_all?: Maybe<Scalars['Float']>;
  /**
   * Residual factor R for reflections that satisfy the resolution
   *  limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low and the observation limit established by
   *  _reflns.observed_criterion.
   * 
   *  _refine.ls_R_factor_obs should not be confused with
   *  _refine.ls_R_factor_R_work; the former reports the results of a
   *  refinement in which all observed reflections were used, the
   *  latter a refinement in which a subset of the observed
   *  reflections were excluded from refinement for the calculation
   *  of a 'free' R factor. However, it would be meaningful to quote
   *  both values if a 'free' R factor were calculated for most of
   *  the refinement, but all of the observed reflections were used
   *  in the final rounds of refinement; such a protocol should be
   *  explained in _refine.details.
   * 
   *      sum|F~obs~ - F~calc~|
   *  R = ---------------------
   *           sum|F~obs~|
   * 
   *  F~obs~  = the observed structure-factor amplitudes
   *  F~calc~ = the calculated structure-factor amplitudes
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_R_factor_obs?: Maybe<Scalars['Float']>;
  /**
   * The smallest value for the interplanar spacings for the
   *  reflection data used in the refinement in angstroms. This is
   *  called the highest resolution.
   */
  readonly ls_d_res_high?: Maybe<Scalars['Float']>;
  /**
   * The largest value for the interplanar spacings for
   *  the reflection data used in the refinement in angstroms.
   *  This is called the lowest resolution.
   */
  readonly ls_d_res_low?: Maybe<Scalars['Float']>;
  /**
   * Type of matrix used to accumulate the least-squares derivatives.
   * 
   * Allowable values:
   * atomblock, diagonal, full, fullcycle, sparse, userblock
   */
  readonly ls_matrix_type?: Maybe<Scalars['String']>;
  /**
   * The number of parameters refined in the least-squares process.
   *  If possible, this number should include some contribution from
   *  the restrained parameters. The restrained parameters are
   *  distinct from the constrained parameters (where one or more
   *  parameters are linearly dependent on the refined value of
   *  another). Least-squares restraints often depend on geometry or
   *  energy considerations and this makes their direct contribution
   *  to this number, and to the goodness-of-fit calculation,
   *  difficult to assess.
   */
  readonly ls_number_parameters?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low
   *  and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the test
   *  reflections (i.e. were excluded from the refinement) when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   */
  readonly ls_number_reflns_R_free?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low
   *  and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the working
   *  reflections (i.e. were included in the refinement) when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   */
  readonly ls_number_reflns_R_work?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low.
   */
  readonly ls_number_reflns_all?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low
   *  and the observation limit established by
   *  _reflns.observed_criterion.
   */
  readonly ls_number_reflns_obs?: Maybe<Scalars['Int']>;
  /**
   * The number of restrained parameters. These are parameters which
   *  are not directly dependent on another refined parameter.
   *  Restrained parameters often involve geometry or energy
   *  dependencies.
   *  See also _atom_site.constraints and _atom_site.refinement_flags.
   *  A general description of refinement constraints may appear in
   *  _refine.details.
   */
  readonly ls_number_restraints?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low
   *  and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the test
   *  reflections (i.e. were excluded from the refinement) when the
   *  refinement included the calculation of a 'free' R factor,
   *  expressed as a percentage of the number of geometrically
   *  observable reflections that satisfy the resolution limits.
   */
  readonly ls_percent_reflns_R_free?: Maybe<Scalars['Float']>;
  /**
   * The number of reflections that satisfy the resolution limits
   *  established by _refine.ls_d_res_high and _refine.ls_d_res_low
   *  and the observation limit established by
   *  _reflns.observed_criterion, expressed as a percentage of the
   *  number of geometrically observable reflections that satisfy
   *  the resolution limits.
   */
  readonly ls_percent_reflns_obs?: Maybe<Scalars['Float']>;
  /**
   * The ratio of the total number of observations of the
   *  reflections that satisfy the resolution limits established by
   *  _refine.ls_d_res_high and _refine.ls_d_res_low to the number
   *  of crystallographically unique reflections that satisfy the
   *  same limits.
   */
  readonly ls_redundancy_reflns_all?: Maybe<Scalars['Float']>;
  /**
   * The ratio of the total number of observations of the
   *  reflections that satisfy the resolution limits established by
   *  _refine.ls_d_res_high and _refine.ls_d_res_low and the
   *  observation limit established by _reflns.observed_criterion to
   *  the number of crystallographically unique reflections that
   *  satisfy the same limits.
   */
  readonly ls_redundancy_reflns_obs?: Maybe<Scalars['Float']>;
  /**
   * Weighted residual factor wR for reflections that satisfy the
   *  resolution limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the test
   *  reflections (i.e. were excluded from the refinement) when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   * 
   *       ( sum|w |Y~obs~ - Y~calc~|^2^| )^1/2^
   *  wR = ( ---------------------------- )
   *       (        sum|w Y~obs~^2^|      )
   * 
   *  Y~obs~  = the observed amplitude specified by
   *            _refine.ls_structure_factor_coef
   *  Y~calc~ = the calculated amplitude specified by
   *            _refine.ls_structure_factor_coef
   *  w       = the least-squares weight
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_wR_factor_R_free?: Maybe<Scalars['Float']>;
  /**
   * Weighted residual factor wR for reflections that satisfy the
   *  resolution limits established by _refine.ls_d_res_high and
   *  _refine.ls_d_res_low and the observation limit established by
   *  _reflns.observed_criterion, and that were used as the working
   *  reflections (i.e. were included in the refinement) when the
   *  refinement included the calculation of a 'free' R factor.
   *  Details of how reflections were assigned to the working and
   *  test sets are given in _reflns.R_free_details.
   * 
   *       ( sum|w |Y~obs~ - Y~calc~|^2^| )^1/2^
   *  wR = ( ---------------------------- )
   *       (        sum|w Y~obs~^2^|      )
   * 
   *  Y~obs~  = the observed amplitude specified by
   *            _refine.ls_structure_factor_coef
   *  Y~calc~ = the calculated amplitude specified by
   *            _refine.ls_structure_factor_coef
   *  w       = the least-squares weight
   * 
   *  sum is taken over the specified reflections
   */
  readonly ls_wR_factor_R_work?: Maybe<Scalars['Float']>;
  /** The maximum value for occupancy found in the coordinate set. */
  readonly occupancy_max?: Maybe<Scalars['Float']>;
  /** The minimum value for occupancy found in the coordinate set. */
  readonly occupancy_min?: Maybe<Scalars['Float']>;
  /**
   * Average figure of merit of phases of reflections not included
   *  in the refinement.
   * 
   *  This value is derived from the likelihood function.
   * 
   *  FOM           = I~1~(X)/I~0~(X)
   * 
   *  I~0~, I~1~     = zero- and first-order modified Bessel functions
   *                  of the first kind
   *  X              = sigma~A~ |E~o~| |E~c~|/SIGMA
   *  E~o~, E~c~     = normalized observed and calculated structure
   *                  factors
   *  sigma~A~       = <cos 2 pi s delta~x~> SQRT(Sigma~P~/Sigma~N~)
   *                  estimated using maximum likelihood
   *  Sigma~P~       = sum~{atoms in model}~ f^2^
   *  Sigma~N~       = sum~{atoms in crystal}~ f^2^
   *  f              = form factor of atoms
   *  delta~x~       = expected error
   *  SIGMA          = (sigma~{E;exp}~)^2^ + epsilon [1-(sigma~A~)^2^]
   *  sigma~{E;exp}~ = uncertainties of normalized observed
   *                  structure factors
   *  epsilon       = multiplicity of the diffracting plane
   * 
   *  Ref: Murshudov, G. N., Vagin, A. A. & Dodson, E. J. (1997).
   *       Acta Cryst. D53, 240-255.
   */
  readonly overall_FOM_free_R_set?: Maybe<Scalars['Float']>;
  /**
   * Average figure of merit of phases of reflections included in
   *  the refinement.
   * 
   *  This value is derived from the likelihood function.
   * 
   *  FOM           = I~1~(X)/I~0~(X)
   * 
   *  I~0~, I~1~     = zero- and first-order modified Bessel functions
   *                  of the first kind
   *  X              = sigma~A~ |E~o~| |E~c~|/SIGMA
   *  E~o~, E~c~     = normalized observed and calculated structure
   *                  factors
   *  sigma~A~       = <cos 2 pi s delta~x~> SQRT(Sigma~P~/Sigma~N~)
   *                  estimated using maximum likelihood
   *  Sigma~P~       = sum~{atoms in model}~ f^2^
   *  Sigma~N~       = sum~{atoms in crystal}~ f^2^
   *  f              = form factor of atoms
   *  delta~x~       = expected error
   *  SIGMA          = (sigma~{E;exp}~)^2^ + epsilon [1-(sigma~A~)^2^]
   *  sigma~{E;exp}~ = uncertainties of normalized observed
   *                  structure factors
   *  epsilon       = multiplicity of the diffracting plane
   * 
   *  Ref: Murshudov, G. N., Vagin, A. A. & Dodson, E. J. (1997).
   *       Acta Cryst. D53, 240-255.
   */
  readonly overall_FOM_work_R_set?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *            of the displacement parameters based on a maximum-likelihood
   *            residual.
   * 
   *            The overall standard uncertainty (sigma~B~)^2^ gives an idea
   *            of the uncertainty in the B values of averagely defined
   *            atoms (atoms with B values equal to the average B value).
   * 
   *                                          N~a~
   * (sigma~B~)^2^ = 8 ----------------------------------------------
   *                   sum~i~ {[1/Sigma - (E~o~)^2^ (1-m^2^)](SUM_AS)s^4^}
   * 
   *            N~a~           = number of atoms
   *            E~o~           = normalized structure factors
   *            m              = figure of merit of phases of reflections
   *                             included in the summation
   *            s              = reciprocal-space vector
   * 
   *            SUM_AS         = (sigma~A~)^2^/Sigma^2^
   *            Sigma          = (sigma~{E;exp}~)^2^ + epsilon [1-(sigma~A~)^2^]
   *            sigma~{E;exp}~  = experimental uncertainties of normalized
   *                             structure factors
   *            sigma~A~        = <cos 2 pi s delta~x~> SQRT(Sigma~P~/Sigma~N~)
   *                             estimated using maximum likelihood
   *            Sigma~P~        = sum~{atoms in model}~ f^2^
   *            Sigma~N~        = sum~{atoms in crystal}~ f^2^
   *            f               = atom form factor 
   *            delta~x~        = expected error
   *            epsilon         = multiplicity of diffracting plane
   * 
   *            summation is over all reflections included in refinement
   * 
   *            Ref: (sigma~A~ estimation) "Refinement of macromolecular
   *                 structures by the maximum-likelihood method",
   *                 Murshudov, G. N., Vagin, A. A. & Dodson, E. J. (1997).
   *                 Acta Cryst. D53, 240-255.
   * 
   *                 (SU B estimation) Murshudov, G. N. & Dodson,
   *                 E. J. (1997). Simplified error estimation a la
   *                 Cruickshank in macromolecular crystallography.
   *                 CCP4 Newsletter on Protein Crystallography, No. 33,
   *                 January 1997, pp. 31-39.
   * 
   *                http://www.ccp4.ac.uk/newsletters/newsletter33/murshudov.html
   */
  readonly overall_SU_B?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *            of the positional parameters based on a maximum likelihood
   *            residual.
   * 
   *            The overall standard uncertainty (sigma~X~)^2^ gives an
   *            idea of the uncertainty in the position of averagely
   *            defined atoms (atoms with B values equal to average B value)
   * 
   *                  3                         N~a~
   * (sigma~X~)^2^  = ---------------------------------------------------------
   *                  8 pi^2^ sum~i~ {[1/Sigma - (E~o~)^2^ (1-m^2^)](SUM_AS)s^2^}
   * 
   *            N~a~           = number of atoms
   *            E~o~           = normalized structure factors
   *            m              = figure of merit of phases of reflections
   *                             included in the summation
   *            s              = reciprocal-space vector
   * 
   *            SUM_AS         = (sigma~A~)^2^/Sigma^2^
   *            Sigma          = (sigma~{E;exp}~)^2^ + epsilon [1-(sigma~A~)^2^]
   *            sigma~{E;exp}~  = experimental uncertainties of normalized
   *                             structure factors
   *            sigma~A~        = <cos 2 pi s delta~x~> SQRT(Sigma~P~/Sigma~N~)
   *                             estimated using maximum likelihood
   *            Sigma~P~        = sum~{atoms in model}~ f^2^
   *            Sigma~N~        = sum~{atoms in crystal}~ f^2^
   *            f               = atom form factor
   *            delta~x~        = expected error
   *            epsilon         = multiplicity of diffracting plane
   * 
   *            summation is over all reflections included in refinement
   * 
   *            Ref: (sigma_A estimation) "Refinement of macromolecular
   *                 structures by the maximum-likelihood method",
   *                 Murshudov, G. N., Vagin, A. A. & Dodson, E. J. (1997).
   *                 Acta Cryst. D53, 240-255.
   * 
   *                 (SU ML estimation) Murshudov, G. N. & Dodson,
   *                 E. J. (1997). Simplified error estimation a la
   *                 Cruickshank in macromolecular crystallography.
   *                 CCP4 Newsletter on Protein Crystallography, No. 33,
   *                 January 1997, pp. 31-39.
   * 
   *                http://www.ccp4.ac.uk/newsletters/newsletter33/murshudov.html
   */
  readonly overall_SU_ML?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *  of the displacement parameters based on the crystallographic
   *  R value, expressed in a formalism known as the dispersion
   *  precision indicator (DPI).
   * 
   *  The overall standard uncertainty (sigma~B~) gives an idea
   *  of the uncertainty in the B values of averagely defined
   *  atoms (atoms with B values equal to the average B value).
   * 
   *                         N~a~
   *  (sigma~B~)^2^ = 0.65 ---------- (R~value~)^2^ (D~min~)^2^ C^-2/3^
   *                       (N~o~-N~p~)
   * 
   * 
   *  N~a~     = number of atoms included in refinement
   *  N~o~     = number of observations
   *  N~p~     = number of parameters refined
   *  R~value~ = conventional crystallographic R value
   *  D~min~   = maximum resolution
   *  C        = completeness of data
   * 
   *  Ref: Cruickshank, D. W. J. (1999). Acta Cryst. D55, 583-601.
   * 
   *       Murshudov, G. N. & Dodson,
   *       E. J. (1997). Simplified error estimation a la
   *       Cruickshank in macromolecular crystallography.
   *       CCP4 Newsletter on Protein Crystallography, No. 33,
   *       January 1997, pp. 31-39.
   * 
   *      http://www.ccp4.ac.uk/newsletters/newsletter33/murshudov.html
   */
  readonly overall_SU_R_Cruickshank_DPI?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *  of the displacement parameters based on the free R value.
   * 
   *  The overall standard uncertainty (sigma~B~) gives an idea
   *  of the uncertainty in the B values of averagely defined
   *  atoms (atoms with B values equal to the average B value).
   * 
   *                         N~a~
   *  (sigma~B~)^2^ = 0.65 ---------- (R~free~)^2^ (D~min~)^2^ C^-2/3^
   *                       (N~o~-N~p~)
   * 
   * 
   *  N~a~     = number of atoms included in refinement
   *  N~o~     = number of observations
   *  N~p~     = number of parameters refined
   *  R~free~  = conventional free crystallographic R value calculated
   *           using reflections not included in refinement
   *  D~min~   = maximum resolution
   *  C        = completeness of data
   * 
   *  Ref: Cruickshank, D. W. J. (1999). Acta Cryst. D55, 583-601.
   * 
   *       Murshudov, G. N. & Dodson,
   *       E. J. (1997). Simplified error estimation a la
   *       Cruickshank in macromolecular crystallography.
   *       CCP4 Newsletter on Protein Crystallography, No. 33,
   *       January 1997, pp. 31-39.
   * 
   *      http://www.ccp4.ac.uk/newsletters/newsletter33/murshudov.html
   */
  readonly overall_SU_R_free?: Maybe<Scalars['Float']>;
  /**
   * Details of the manner in which the cross validation
   *  reflections were selected.
   * 
   * Examples:
   * Random selection
   */
  readonly pdbx_R_Free_selection_details?: Maybe<Scalars['String']>;
  /**
   * A flag for TLS refinements identifying the type of atomic displacement parameters stored
   *  in _atom_site.B_iso_or_equiv.
   * 
   * Allowable values:
   * LIKELY RESIDUAL, UNVERIFIED
   */
  readonly pdbx_TLS_residual_ADP_flag?: Maybe<Scalars['String']>;
  /**
   * Average Fourier Shell Correlation (avgFSC) between model and
   *  observed structure factors for reflections not included in refinement.
   * 
   *  The average FSC is a measure of the agreement between observed
   *  and calculated structure factors.
   * 
   *                   sum(N~i~ FSC~free-i~)
   *  avgFSC~free~   = ---------------------
   *                   sum(N~i~)
   * 
   * 
   *  N~i~          = the number of free reflections in the resolution shell i
   *  FSC~free-i~   = FSC for free reflections in the i-th resolution shell calculated as:
   * 
   *                 (sum(|F~o~| |F~c~| fom cos(phi~c~-phi~o~)))
   *  FSC~free-i~  = -------------------------------------------
   *                 (sum(|F~o~|^2^) (sum(|F~c~|^2^)))^1/2^
   * 
   *  |F~o~|   = amplitude of observed structure factor
   *  |F~c~|   = amplitude of calculated structure factor
   *  phi~o~   = phase of observed structure factor
   *  phi~c~   = phase of calculated structure factor
   *  fom      = figure of merit of the experimental phases.
   * 
   *  Summation of FSC~free-i~ is carried over all free reflections in the resolution shell.
   * 
   *  Summation of avgFSC~free~ is carried over all resolution shells.
   * 
   * 
   *  Ref:  Rosenthal P.B., Henderson R.
   *        "Optimal determination of particle orientation, absolute hand,
   *        and contrast loss in single-particle electron cryomicroscopy.
   *        Journal of Molecular Biology. 2003;333(4):721-745, equation (A6).
   */
  readonly pdbx_average_fsc_free?: Maybe<Scalars['Float']>;
  /**
   * Overall average Fourier Shell Correlation (avgFSC) between model and
   *  observed structure factors for all reflections.
   * 
   *  The average FSC is a measure of the agreement between observed
   *  and calculated structure factors.
   * 
   *             sum(N~i~ FSC~i~)
   *  avgFSC   = ----------------
   *             sum(N~i~)
   * 
   * 
   *  N~i~     = the number of all reflections in the resolution shell i
   *  FSC~i~   = FSC for all reflections in the i-th resolution shell calculated as:
   * 
   *            (sum(|F~o~| |F~c~| fom cos(phi~c~-phi~o~)))
   *  FSC~i~  = -------------------------------------------
   *            (sum(|F~o~|^2^) (sum(|F~c~|^2^)))^1/2^
   * 
   *  |F~o~|   = amplitude of observed structure factor
   *  |F~c~|   = amplitude of calculated structure factor
   *  phi~o~   = phase of observed structure factor
   *  phi~c~   = phase of calculated structure factor
   *  fom      = figure of merit of the experimental phases.
   * 
   *  Summation of FSC~i~ is carried over all reflections in the resolution shell.
   * 
   *  Summation of avgFSC is carried over all resolution shells.
   * 
   * 
   *  Ref:  Rosenthal P.B., Henderson R.
   *        "Optimal determination of particle orientation, absolute hand,
   *        and contrast loss in single-particle electron cryomicroscopy.
   *        Journal of Molecular Biology. 2003;333(4):721-745, equation (A6).
   */
  readonly pdbx_average_fsc_overall?: Maybe<Scalars['Float']>;
  /**
   * Average Fourier Shell Correlation (avgFSC) between model and 
   *  observed structure factors for reflections included in refinement.
   * 
   *  The average FSC is a measure of the agreement between observed 
   *  and calculated structure factors.
   * 
   *                   sum(N~i~ FSC~work-i~) 
   *  avgFSC~work~   = ---------------------
   *                   sum(N~i~)
   * 
   * 
   *  N~i~          = the number of working reflections in the resolution shell i
   *  FSC~work-i~   = FSC for working reflections in the i-th resolution shell calculated as:
   * 
   *                 (sum(|F~o~| |F~c~| fom cos(phi~c~-phi~o~)))
   *  FSC~work-i~  = -------------------------------------------
   *                 (sum(|F~o~|^2^) (sum(|F~c~|^2^)))^1/2^
   * 
   *  |F~o~|   = amplitude of observed structure factor
   *  |F~c~|   = amplitude of calculated structure factor
   *  phi~o~   = phase of observed structure factor
   *  phi~c~   = phase of calculated structure factor
   *  fom      = figure of merit of the experimental phases.
   * 
   *  Summation of FSC~work-i~ is carried over all working reflections in the resolution shell.
   * 
   *  Summation of avgFSC~work~ is carried over all resolution shells.
   * 
   * 
   *  Ref:  Rosenthal P.B., Henderson R.
   *        "Optimal determination of particle orientation, absolute hand,
   *        and contrast loss in single-particle electron cryomicroscopy.
   *        Journal of Molecular Biology. 2003;333(4):721-745, equation (A6).
   */
  readonly pdbx_average_fsc_work?: Maybe<Scalars['Float']>;
  /**
   * Value of F at "high end" of data cutoff.
   * 
   * Examples:
   * 17600
   */
  readonly pdbx_data_cutoff_high_absF?: Maybe<Scalars['Float']>;
  /**
   * Value of RMS |F| used as high data cutoff.
   * 
   * Examples:
   * 205.1
   */
  readonly pdbx_data_cutoff_high_rms_absF?: Maybe<Scalars['Float']>;
  /**
   * Value of F at "low end" of data cutoff.
   * 
   * Examples:
   * 0.30
   */
  readonly pdbx_data_cutoff_low_absF?: Maybe<Scalars['Float']>;
  /**
   * An identifier for the diffraction data set used in this refinement.
   * 
   *  Multiple diffraction data sets specified as a comma separated list.
   */
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Whether the structure was refined with indvidual
   * isotropic, anisotropic or overall temperature factor.
   * 
   * Examples:
   * Isotropic, Overall
   */
  readonly pdbx_isotropic_thermal_model?: Maybe<Scalars['String']>;
  /**
   * Whether the cross validataion method was used through
   * out or only at the end.
   * 
   * Examples:
   * FREE R-VALUE
   */
  readonly pdbx_ls_cross_valid_method?: Maybe<Scalars['String']>;
  /** Data cutoff (SIGMA(F)) */
  readonly pdbx_ls_sigma_F?: Maybe<Scalars['Float']>;
  /** Data cutoff (SIGMA(F^2)) */
  readonly pdbx_ls_sigma_Fsqd?: Maybe<Scalars['Float']>;
  /** Data cutoff (SIGMA(I)) */
  readonly pdbx_ls_sigma_I?: Maybe<Scalars['Float']>;
  /**
   * Method(s) used to determine the structure.
   * 
   * Examples:
   * AB INITIO PHASING, DM, ISAS, ISIR, ISIRAS, MAD, MIR, MIRAS, MR, SIR, SIRAS
   */
  readonly pdbx_method_to_determine_struct?: Maybe<Scalars['String']>;
  /**
   * Overall estimated standard uncertainties of positional
   *  parameters based on R value.
   */
  readonly pdbx_overall_ESU_R?: Maybe<Scalars['Float']>;
  /** Overall estimated standard uncertainties of positional parameters based on R free value. */
  readonly pdbx_overall_ESU_R_Free?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *  of the displacement parameters based on the crystallographic
   *  R value, expressed in a formalism known as the dispersion
   *  precision indicator (DPI).
   * 
   *  Ref: Blow, D (2002) Acta Cryst. D58, 792-797
   */
  readonly pdbx_overall_SU_R_Blow_DPI?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *  of the displacement parameters based on the crystallographic
   *  R-free value, expressed in a formalism known as the dispersion
   *  precision indicator (DPI).
   * 
   *  Ref: Blow, D (2002) Acta Cryst. D58, 792-797
   */
  readonly pdbx_overall_SU_R_free_Blow_DPI?: Maybe<Scalars['Float']>;
  /**
   * The overall standard uncertainty (estimated standard deviation)
   *  of the displacement parameters based on the crystallographic
   *  R-free value, expressed in a formalism known as the dispersion
   *  precision indicator (DPI).
   * 
   *  Ref: Cruickshank, D. W. J. (1999). Acta Cryst. D55, 583-601.
   */
  readonly pdbx_overall_SU_R_free_Cruickshank_DPI?: Maybe<Scalars['Float']>;
  /**
   * The overall phase error for all reflections after refinement using
   *  the current refinement target.
   * 
   * Examples:
   * 0.30
   */
  readonly pdbx_overall_phase_error?: Maybe<Scalars['Float']>;
  /**
   * This data item uniquely identifies a refinement within an entry.
   *  _refine.pdbx_refine_id can be used to distinguish the results of 
   *  joint refinements.
   */
  readonly pdbx_refine_id: Scalars['String'];
  /** For bulk solvent mask calculation, the amount that the ionic radii of atoms, which can be ions, are increased used. */
  readonly pdbx_solvent_ion_probe_radii?: Maybe<Scalars['Float']>;
  /** For bulk solvent mask calculation, amount mask is shrunk after taking away atoms with new radii and a constant value assigned to this new region. */
  readonly pdbx_solvent_shrinkage_radii?: Maybe<Scalars['Float']>;
  /** For bulk solvent mask calculation, the value by which the vdw radii of non-ion atoms (like carbon) are increased and used. */
  readonly pdbx_solvent_vdw_probe_radii?: Maybe<Scalars['Float']>;
  /**
   * Starting model for refinement.  Starting model for
   *  molecular replacement should refer to a previous
   *  structure or experiment.
   * 
   * Examples:
   * 1XYZ, 2XYZ, BDL001
   */
  readonly pdbx_starting_model?: Maybe<Scalars['String']>;
  /**
   * Special case of stereochemistry target values used
   * in SHELXL refinement.
   */
  readonly pdbx_stereochem_target_val_spec_case?: Maybe<Scalars['String']>;
  /** Stereochemistry target values used in refinement. */
  readonly pdbx_stereochemistry_target_values?: Maybe<Scalars['String']>;
  /** Special aspects of the solvent model used during refinement. */
  readonly solvent_model_details?: Maybe<Scalars['String']>;
  /**
   * The value of the BSOL solvent-model parameter describing
   *  the average isotropic displacement parameter of disordered
   *  solvent atoms.
   * 
   *  This is one of the two parameters (the other is
   *  _refine.solvent_model_param_ksol) in Tronrud's method of
   *  modelling the contribution of bulk solvent to the
   *  scattering. The standard scale factor is modified according
   *  to the expression
   * 
   *      k0 exp(-B0 * s^2^)[1-KSOL * exp(-BSOL * s^2^)]
   * 
   *  where k0 and B0 are the scale factors for the protein.
   * 
   *  Ref: Tronrud, D. E. (1997). Methods Enzymol. 277, 243-268.
   */
  readonly solvent_model_param_bsol?: Maybe<Scalars['Float']>;
  /**
   * The value of the KSOL solvent-model parameter describing
   *  the ratio of the electron density in the bulk solvent to the
   *  electron density in the molecular solute.
   * 
   *  This is one of the two parameters (the other is
   *  _refine.solvent_model_param_bsol) in Tronrud's method of
   *  modelling the contribution of bulk solvent to the
   *  scattering. The standard scale factor is modified according
   *  to the expression
   * 
   *      k0 exp(-B0 * s^2^)[1-KSOL * exp(-BSOL * s^2^)]
   * 
   *  where k0 and B0 are the scale factors for the protein.
   * 
   *  Ref: Tronrud, D. E. (1997). Methods Enzymol. 277, 243-268.
   */
  readonly solvent_model_param_ksol?: Maybe<Scalars['Float']>;
};

export type PdbxChemCompIdentifier = {
  /**
   * This data item is a pointer to _chem_comp.id in the CHEM_COMP
   *  category.
   */
  readonly comp_id: Scalars['String'];
  /**
   * This data item contains the identifier value for this 
   *  component.
   */
  readonly identifier?: Maybe<Scalars['String']>;
  /**
   * This data item contains the name of the program
   *  or library used to compute the identifier.
   * 
   * Examples:
   * OPENEYE, DAYLIGHT, ACD, AUTONOM, PUBCHEM_CID, PUBCHEM_SID, OTHER, NONE
   */
  readonly program: Scalars['String'];
  /**
   * This data item contains the version of the program
   *  or library used to compute the identifier.
   */
  readonly program_version: Scalars['String'];
  /**
   * This data item contains the identifier type.
   * 
   * Allowable values:
   * CAS REGISTRY NUMBER, COMMON NAME, CONDENSED IUPAC CARB SYMBOL, CONDENSED IUPAC CARBOHYDRATE SYMBOL, IUPAC CARB SYMBOL, IUPAC CARBOHYDRATE SYMBOL, MDL Identifier, PUBCHEM Identifier, SNFG CARB SYMBOL, SNFG CARBOHYDRATE SYMBOL, SYNONYM, SYSTEMATIC NAME
   */
  readonly type: Scalars['String'];
};

export type RcsbPolymerInstanceFeature = {
  /** Identifies the version of the feature assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the feature. */
  readonly description?: Maybe<Scalars['String']>;
  /** An identifier for the feature. */
  readonly feature_id?: Maybe<Scalars['String']>;
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeatureFeaturePositions>>>;
  /** A name for the feature. */
  readonly name?: Maybe<Scalars['String']>;
  /** Ordinal identifier for this category */
  readonly ordinal: Scalars['Int'];
  /**
   * Code identifying the individual, organization or program that
   *  assigned the feature.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * Code residue coordinate system for the assigned feature.
   * 
   * Allowable values:
   * NCBI, PDB entity, PDB entry, UniProt
   */
  readonly reference_scheme?: Maybe<Scalars['String']>;
  /**
   * A type or category of the feature.
   * 
   * Allowable values:
   * ANGLE_OUTLIER, BINDING_SITE, BOND_OUTLIER, CATH, CIS-PEPTIDE, HELIX_P, MOGUL_ANGLE_OUTLIER, MOGUL_BOND_OUTLIER, RAMACHANDRAN_OUTLIER, ROTAMER_OUTLIER, RSRCC_OUTLIER, RSRZ_OUTLIER, SCOP, SHEET, UNASSIGNED_SEC_STRUCT, UNOBSERVED_ATOM_XYZ, UNOBSERVED_RESIDUE_XYZ, ZERO_OCCUPANCY_ATOM_XYZ, ZERO_OCCUPANCY_RESIDUE_XYZ
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxAuditRevisionCategory = {
  /**
   * The category updated in the pdbx_audit_revision_category record.
   * 
   * Examples:
   * audit_author, citation
   */
  readonly category?: Maybe<Scalars['String']>;
  /**
   * The type of file that the pdbx_audit_revision_history record refers to.
   * 
   * Allowable values:
   * Chemical component, NMR restraints, NMR shifts, Structure factors, Structure model
   */
  readonly data_content_type: Scalars['String'];
  /**
   * A unique identifier for the pdbx_audit_revision_category record.
   * 
   * Examples:
   * 1
   */
  readonly ordinal: Scalars['Int'];
  /**
   * A pointer to  _pdbx_audit_revision_history.ordinal
   * 
   * Examples:
   * 1
   */
  readonly revision_ordinal: Scalars['Int'];
};

export type RcsbPolymerEntityContainerIdentifiers = {
  /** Instance identifiers corresponding to copies of the entity in this container. */
  readonly asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Author instance identifiers corresponding to copies of the entity in this container. */
  readonly auth_asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Unique list of monomer chemical component identifiers in the polymer entity in this container. */
  readonly chem_comp_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Unique list of non-standard monomer chemical component identifiers in the polymer entity in this container. */
  readonly chem_comp_nstd_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Entity identifier for the container. */
  readonly entity_id: Scalars['String'];
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /** The BIRD identifier for the entity in this container. */
  readonly prd_id?: Maybe<Scalars['String']>;
  /**
   * A unique identifier for each object in this entity container formed by
   *  an underscore separated concatenation of entry and entity identifiers.
   */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  readonly reference_sequence_identifiers?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityContainerIdentifiersReferenceSequenceIdentifiers>>>;
  readonly uniprot_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type CoreAssembly = {
  /** Get PDB entry that includes this assembly. */
  readonly entry?: Maybe<CoreEntry>;
  readonly pdbx_struct_assembly?: Maybe<PdbxStructAssembly>;
  readonly pdbx_struct_assembly_auth_evidence?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyAuthEvidence>>>;
  readonly pdbx_struct_assembly_gen?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyGen>>>;
  readonly pdbx_struct_assembly_prop?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyProp>>>;
  readonly pdbx_struct_oper_list?: Maybe<ReadonlyArray<Maybe<PdbxStructOperList>>>;
  readonly rcsb_assembly_container_identifiers: RcsbAssemblyContainerIdentifiers;
  readonly rcsb_assembly_info?: Maybe<RcsbAssemblyInfo>;
  /**
   * A unique identifier for each object in this assembly container formed by
   *  a dash separated concatenation of entry and assembly identifiers.
   */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>;
  readonly rcsb_struct_symmetry?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetry>>>;
  readonly rcsb_struct_symmetry_lineage?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetryLineage>>>;
  /** The title and version of software package used for symmetry calculations. */
  readonly rcsb_struct_symmetry_provenance_code?: Maybe<Scalars['String']>;
};

export type Exptl = {
  /**
   * The total number of crystals used in the  measurement of
   *  intensities.
   */
  readonly crystals_number?: Maybe<Scalars['Int']>;
  /**
   * Any special information about the experimental work prior to the
   *  intensity measurement. See also _exptl_crystal.preparation.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The method used in the experiment.
   * 
   * Allowable values:
   * ELECTRON CRYSTALLOGRAPHY, ELECTRON MICROSCOPY, EPR, FIBER DIFFRACTION, FLUORESCENCE TRANSFER, INFRARED SPECTROSCOPY, NEUTRON DIFFRACTION, POWDER DIFFRACTION, SOLID-STATE NMR, SOLUTION NMR, SOLUTION SCATTERING, THEORETICAL MODEL, X-RAY DIFFRACTION
   */
  readonly method: Scalars['String'];
  /**
   * A description of special aspects of the experimental method.
   * 
   * Examples:
   * 29 structures, minimized average structure
   */
  readonly method_details?: Maybe<Scalars['String']>;
};

export type RcsbPolymerStructConnConnectTarget = {
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.auth_asym_id in the
   *  ATOM_SITE category.
   */
  readonly auth_asym_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.auth_seq_id in the
   *  ATOM_SITE category.
   */
  readonly auth_seq_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_alt_id in the
   *  ATOM_SITE category.
   */
  readonly label_alt_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_asym_id in the
   *  ATOM_SITE category.
   */
  readonly label_asym_id: Scalars['String'];
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_atom_id in the
   *  ATOM_SITE category.
   */
  readonly label_atom_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_comp_id in the
   *  ATOM_SITE category.
   */
  readonly label_comp_id: Scalars['String'];
  /**
   * A component of the identifier for the target of the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.connect_target_label_seq_id in the
   *  ATOM_SITE category.
   */
  readonly label_seq_id?: Maybe<Scalars['Int']>;
  /**
   * Describes the symmetry operation that should be applied to the
   *  atom set specified by _rcsb_polymer_struct_conn.label* to generate the
   *  target of the structure connection.
   * 
   * Examples:
   * 1_555, 7_645
   */
  readonly symmetry?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityAlign = {
  readonly aligned_regions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAlignAlignedRegions>>>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the reference sequence.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /** Reference sequence accession code. */
  readonly reference_database_accession?: Maybe<Scalars['String']>;
  /** Reference sequence isoform identifier. */
  readonly reference_database_isoform?: Maybe<Scalars['String']>;
  /**
   * Reference sequence database name.
   * 
   * Allowable values:
   * EMBL, GenBank, NDB, NORINE, PDB, PIR, PRF, RefSeq, UniProt
   */
  readonly reference_database_name?: Maybe<Scalars['String']>;
};

export type CorePolymerEntityInstance = {
  readonly pdbx_struct_special_symmetry?: Maybe<ReadonlyArray<Maybe<PdbxStructSpecialSymmetry>>>;
  /** Get polymer entity for this polymer entity instance. */
  readonly polymer_entity?: Maybe<CorePolymerEntity>;
  /**
   * A unique identifier for each object in this entity instance container formed by
   *  an 'dot' (.) separated concatenation of entry and entity instance identifiers.
   */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>;
  readonly rcsb_polymer_entity_instance_container_identifiers?: Maybe<RcsbPolymerEntityInstanceContainerIdentifiers>;
  readonly rcsb_polymer_instance_annotation?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceAnnotation>>>;
  readonly rcsb_polymer_instance_feature?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeature>>>;
  readonly rcsb_polymer_instance_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeatureSummary>>>;
  readonly rcsb_polymer_struct_conn?: Maybe<ReadonlyArray<Maybe<RcsbPolymerStructConn>>>;
};

export type RcsbPolymerInstanceAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityRcsbEnzymeClassCombined = {
  /** The enzyme classification hierarchy depth (1-N). */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Combined list of enzyme class assignments. */
  readonly ec?: Maybe<Scalars['String']>;
  /**
   * Combined list of enzyme class associated provenance sources.
   * 
   * Allowable values:
   * PDB Primary Data, UniProt
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityFeatureSummary = {
  /** The feature count. */
  readonly count?: Maybe<Scalars['Int']>;
  /** The fractional feature coverage relative to the full entity sequence. */
  readonly coverage?: Maybe<Scalars['Float']>;
  /** The maximum feature length. */
  readonly maximum_length?: Maybe<Scalars['Int']>;
  /** The maximum feature value. */
  readonly maximum_value?: Maybe<Scalars['Float']>;
  /** The minimum feature length. */
  readonly minimum_length?: Maybe<Scalars['Int']>;
  /** The minimum feature value. */
  readonly minimum_value?: Maybe<Scalars['Float']>;
  /**
   * Type or category of the feature.
   * 
   * Allowable values:
   * artifact, modified_monomer, mutation
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbStructSymmetryClusters = {
  /** Average RMSD between members of a given cluster. */
  readonly avg_rmsd?: Maybe<Scalars['Float']>;
  /** Subunits that belong to the cluster, identified by asym_id and optionally by assembly operator id(s). */
  readonly members: ReadonlyArray<Maybe<ClustersMembers>>;
};

export type RcsbChemCompAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type RcsbPfamContainerIdentifiers = {
  /** Accession number of Pfam entry. */
  readonly pfam_id?: Maybe<Scalars['String']>;
};

export type PdbxSgProject = {
  /**
   * The value identifies the full name of center.
   * 
   * Allowable values:
   * Accelerated Technologies Center for Gene to 3D Structure, Assembly, Dynamics and Evolution of Cell-Cell and Cell-Matrix Adhesions, Atoms-to-Animals: The Immune Function Network, Bacterial targets at IGS-CNRS, France, Berkeley Structural Genomics Center, Center for Eukaryotic Structural Genomics, Center for High-Throughput Structural Biology, Center for Membrane Proteins of Infectious Diseases, Center for Structural Genomics of Infectious Diseases, Center for Structures of Membrane Proteins, Center for the X-ray Structure Determination of Human Transporters, Chaperone-Enabled Studies of Epigenetic Regulation Enzymes, Enzyme Discovery for Natural Product Biosynthesis, GPCR Network, Integrated Center for Structure and Function Innovation, Israel Structural Proteomics Center, Joint Center for Structural Genomics, Marseilles Structural Genomics Program @ AFMB, Medical Structural Genomics of Pathogenic Protozoa, Membrane Protein Structural Biology Consortium, Membrane Protein Structures by Solution NMR, Midwest Center for Macromolecular Research, Midwest Center for Structural Genomics, Mitochondrial Protein Partnership, Montreal-Kingston Bacterial Structural Genomics Initiative, Mycobacterium Tuberculosis Structural Proteomics Project, New York Consortium on Membrane Protein Structure, New York SGX Research Center for Structural Genomics, New York Structural GenomiX Research Consortium, New York Structural Genomics Research Consortium, Northeast Structural Genomics Consortium, Nucleocytoplasmic Transport: a Target for Cellular Control, Ontario Centre for Structural Proteomics, Oxford Protein Production Facility, Paris-Sud Yeast Structural Genomics, Partnership for Nuclear Receptor Signaling Code Biology, Partnership for Stem Cell Biology, Partnership for T-Cell Biology, Program for the Characterization of Secreted Effector Proteins, Protein Structure Factory, RIKEN Structural Genomics/Proteomics Initiative, Scottish Structural Proteomics Facility, Seattle Structural Genomics Center for Infectious Disease, South Africa Structural Targets Annotation Database, Southeast Collaboratory for Structural Genomics, Structural Genomics Consortium, Structural Genomics Consortium for Research on Gene Expression, Structural Genomics of Pathogenic Protozoa Consortium, Structural Proteomics in Europe, Structural Proteomics in Europe 2, Structure 2 Function Project, Structure, Dynamics and Activation Mechanisms of Chemokine Receptors, Structure-Function Analysis of Polymorphic CDI Toxin-Immunity Protein Complexes, Structure-Function Studies of Tight Junction Membrane Proteins, Structures of Mtb Proteins Conferring Susceptibility to Known Mtb Inhibitors, TB Structural Genomics Consortium, Transcontinental EM Initiative for Membrane Protein Structure, Transmembrane Protein Center
   */
  readonly full_name_of_center?: Maybe<Scalars['String']>;
  /**
   * A unique integer identifier for this center
   * 
   * Allowable values:
   * 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
   */
  readonly id: Scalars['Int'];
  /**
   * The value identifies the full name of center.
   * 
   * Allowable values:
   * ATCG3D, BIGS, BSGC, BSGI, CEBS, CELLMAT, CESG, CHSAM, CHTSB, CSGID, CSMP, GPCR, IFN, ISFI, ISPC, JCSG, MCMR, MCSG, MPID, MPP, MPSBC, MPSbyNMR, MSGP, MSGPP, MTBI, NESG, NHRs, NPCXstals, NYCOMPS, NYSGRC, NYSGXRC, NatPro, OCSP, OPPF, PCSEP, PSF, RSGI, S2F, SASTAD, SECSG, SGC, SGCGES, SGPP, SPINE, SPINE-2, SSGCID, SSPF, STEMCELL, TBSGC, TCELL, TEMIMPS, TJMP, TMPC, TransportPDB, UC4CDI, XMTB, YSG
   */
  readonly initial_of_center?: Maybe<Scalars['String']>;
  /**
   * The value identifies the Structural Genomics project.
   * 
   * Allowable values:
   * Enzyme Function Initiative, NIAID, National Institute of Allergy and Infectious Diseases, NPPSFA, National Project on Protein Structural and Functional Analyses, PSI, Protein Structure Initiative, PSI:Biology
   */
  readonly project_name?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityFeatureFeaturePositions = {
  /** An identifier for the monomer(s) corresponding to the feature assignment. */
  readonly beg_comp_id?: Maybe<Scalars['String']>;
  /** An identifier for the monomer at which this segment of the feature begins. */
  readonly beg_seq_id: Scalars['Int'];
  /** An identifier for the monomer at which this segment of the feature ends. */
  readonly end_seq_id?: Maybe<Scalars['Int']>;
  /** The value for the feature over this monomer segment. */
  readonly value?: Maybe<Scalars['Float']>;
};

export type RcsbAssemblyContainerIdentifiers = {
  /** Assembly identifier for the container. */
  readonly assembly_id: Scalars['String'];
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /**
   * A unique identifier for each object in this assembly container formed by
   *  a dash separated concatenation of entry and assembly identifiers.
   */
  readonly rcsb_id?: Maybe<Scalars['String']>;
};

export type RcsbPolymerStructConnConnectPartner = {
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_alt_id in the
   *  ATOM_SITE category.
   */
  readonly label_alt_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_asym_id in the
   *  ATOM_SITE category.
   */
  readonly label_asym_id: Scalars['String'];
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _chem_comp_atom.atom_id in the
   *  CHEM_COMP_ATOM category.
   */
  readonly label_atom_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_comp_id in the
   *  ATOM_SITE category.
   */
  readonly label_comp_id: Scalars['String'];
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_seq_id in the
   *  ATOM_SITE category.
   */
  readonly label_seq_id?: Maybe<Scalars['Int']>;
  /**
   * Describes the symmetry operation that should be applied to the
   *  atom set specified by _rcsb_polymer_struct_conn.connect_partner_label* to generate the
   *  partner in the structure connection.
   * 
   * Examples:
   * 1_555, 7_645
   */
  readonly symmetry?: Maybe<Scalars['String']>;
};

export type RcsbUniprotAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbUniprotAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * disease, phenotype
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxSolnScatter = {
  /**
   * The name of the buffer used for the sample in the solution scattering
   *  experiment.
   * 
   * Examples:
   * acetic acid
   */
  readonly buffer_name?: Maybe<Scalars['String']>;
  /**
   * The concentration range (mg/mL) of the complex in the 
   *  sample used in the solution scattering experiment to
   *  determine the mean radius of structural elongation.
   * 
   * Examples:
   * 0.7 - 14
   */
  readonly concentration_range?: Maybe<Scalars['String']>;
  /**
   * A list of the software used in the data analysis
   * 
   * Examples:
   * SCTPL5 GNOM
   */
  readonly data_analysis_software_list?: Maybe<Scalars['String']>;
  /**
   * A list of the software used in the data reduction
   * 
   * Examples:
   * OTOKO
   */
  readonly data_reduction_software_list?: Maybe<Scalars['String']>;
  /**
   * The particular radiation detector. In general this will be a
   *   manufacturer, description, model number or some combination of
   *   these.
   */
  readonly detector_specific?: Maybe<Scalars['String']>;
  /** The general class of the radiation detector. */
  readonly detector_type?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_soln_scatter.id must
   *  uniquely identify the sample in the category PDBX_SOLN_SCATTER
   */
  readonly id: Scalars['String'];
  /**
   * The maximum mean radius of structural elongation of the sample.
   *  In a given solute-solvent contrast, the radius of gyration
   *  R_G is a measure of structural elongation if the internal
   *  inhomogeneity of scattering densities has no effect. Guiner 
   *  analysis at low Q give the R_G and the forward scattering at 
   *  zero angle I(0).
   * 
   *     lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   *  where 
   *        Q = 4(pi)sin(theta/lamda)
   *        2theta = scattering angle
   *        lamda = wavelength
   * 
   *  The above expression is valid in a QR_G range for extended 
   *  rod-like particles. The relative I(0)/c values ( where
   *   c = sample concentration) for sample measurements in a
   *  constant buffer for a single sample data session, gives the
   *  relative masses of the protein(s) studied when referenced
   *  against a standard.
   * 
   *  see:
   *      O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *       X-ray Scattering, Academic Press, New York.
   *      O.Kratky. (1963). X-ray small angle scattering with
   *       substances of biological interest in diluted solutions.
   *       Prog. Biophys. Chem., 13, 105-173.
   *      G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *       of X-ray and neutron scatter from rigid rods of non-uniform
   *       cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   *  If the structure is elongated, the mean radius of gyration
   *  of the cross-sectional structure R_XS  and the mean cross sectional
   *  intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   *     ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly max_mean_cross_sectional_radii_gyration?: Maybe<Scalars['Float']>;
  /**
   * The estimated standard deviation for the
   * minimum mean radius of structural elongation of the sample.
   * In a given solute-solvent contrast, the radius of gyration
   * R_G is a measure of structural elongation if the internal
   * inhomogeneity of scattering densities has no effect. Guiner 
   * analysis at low Q give the R_G and the forward scattering at 
   * zero angle I(0).
   * 
   *     lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   * where 
   *       Q = 4(pi)sin(theta/lamda)
   *       2theta = scattering angle
   *       lamda = wavelength
   * 
   * The above expression is valid in a QR_G range for extended 
   * rod-like particles. The relative I(0)/c values ( where
   *  c = sample concentration) for sample measurements in a
   * constant buffer for a single sample data session, gives the
   * relative masses of the protein(s) studied when referenced
   * against a standard.
   * 
   * see:  
   *     O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *      X-ray Scattering, Academic Press, New York.
   *     O.Kratky. (1963). X-ray small angle scattering with
   *      substances of biological interest in diluted solutions.
   *      Prog. Biophys. Chem., 13, 105-173.
   *     G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *      of X-ray and neutron scatter from rigid rods of non-uniform
   *      cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   * If the structure is elongated, the mean radius of gyration
   * of the cross-sectional structure R_XS  and the mean cross sectional
   * intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   *    ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly max_mean_cross_sectional_radii_gyration_esd?: Maybe<Scalars['Float']>;
  /**
   * The mean radius of structural elongation of the sample.
   *  In a given solute-solvent contrast, the radius of gyration
   *  R_G is a measure of structural elongation if the internal
   *  inhomogeneity of scattering densities has no effect. Guiner 
   *  analysis at low Q gives the R_G and the forward scattering at 
   *  zero angle I(0).
   * 
   *      lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   *  where 
   *        Q = 4(pi)sin(theta/lamda)
   *        2theta = scattering angle
   *        lamda = wavelength
   * 
   *  The above expression is valid in a QR_G range for extended 
   *  rod-like particles. The relative I(0)/c values ( where
   *   c = sample concentration) for sample measurements in a
   *  constant buffer for a single sample data session, gives the
   *  relative masses of the protein(s) studied when referenced
   *  against a standard.
   * 
   *  see: O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *       X-ray Scattering, Academic Press, New York.
   *       O.Kratky. (1963). X-ray small angle scattering with
   *       substances of biological interest in diluted solutions.
   *       Prog. Biophys. Chem., 13, 105-173.
   * 
   *       G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *       of X-ray and neutron scatter from rigid rods of non-uniform
   *       cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   *  If the structure is elongated, the mean radius of gyration
   *  of the cross-sectional structure R_XS  and the mean cross sectional
   *  intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   * 
   *     ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly mean_guiner_radius?: Maybe<Scalars['Float']>;
  /**
   * The estimated standard deviation for the
   *  mean radius of structural elongation of the sample.
   *  In a given solute-solvent contrast, the radius of gyration
   *  R_G is a measure of structural elongation if the internal
   *  inhomogeneity of scattering densities has no effect. Guiner 
   *  analysis at low Q give the R_G and the forward scattering at 
   *  zero angle I(0).
   * 
   *      lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   *  where 
   *        Q = 4(pi)sin(theta/lamda)
   *        2theta = scattering angle
   *        lamda = wavelength
   * 
   *  The above expression is valid in a QR_G range for extended 
   *  rod-like particles. The relative I(0)/c values ( where
   *   c = sample concentration) for sample measurements in a
   *  constant buffer for a single sample data session, gives the
   *  relative masses of the protein(s) studied when referenced
   *  against a standard.
   * 
   *  see: 
   *      O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *       X-ray Scattering, Academic Press, New York.
   *      O.Kratky. (1963). X-ray small angle scattering with
   *       substances of biological interest in diluted solutions.
   *       Prog. Biophys. Chem., 13, 105-173.
   *      G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *       of X-ray and neutron scatter from rigid rods of non-uniform
   *       cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   *  If the structure is elongated, the mean radius of gyration
   *  of the cross-sectional structure R_XS  and the mean cross sectional
   *  intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   *     ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly mean_guiner_radius_esd?: Maybe<Scalars['Float']>;
  /**
   * The minimum mean radius of structural elongation of the sample.
   * In a given solute-solvent contrast, the radius of gyration
   * R_G is a measure of structural elongation if the internal
   * inhomogeneity of scattering densities has no effect. Guiner 
   * analysis at low Q give the R_G and the forward scattering at 
   * zero angle I(0).
   * 
   *     lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   * where 
   *       Q = 4(pi)sin(theta/lamda)
   *       2theta = scattering angle
   *       lamda = wavelength
   * 
   * The above expression is valid in a QR_G range for extended 
   * rod-like particles. The relative I(0)/c values ( where
   *  c = sample concentration) for sample measurements in a
   * constant buffer for a single sample data session, gives the
   * relative masses of the protein(s) studied when referenced
   * against a standard.
   * 
   * see: 
   *     O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *      X-ray Scattering, Academic Press, New York.
   *     O.Kratky. (1963). X-ray small angle scattering with
   *      substances of biological interest in diluted solutions.
   *      Prog. Biophys. Chem., 13, 105-173.
   *     G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *      of X-ray and neutron scatter from rigid rods of non-uniform
   *      cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   * If the structure is elongated, the mean radius of gyration
   * of the cross-sectional structure R_XS  and the mean cross sectional
   * intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   *    ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly min_mean_cross_sectional_radii_gyration?: Maybe<Scalars['Float']>;
  /**
   * The estimated standard deviation for the
   * minimum mean radius of structural elongation of the sample.
   * In a given solute-solvent contrast, the radius of gyration
   * R_G is a measure of structural elongation if the internal
   * inhomogeneity of scattering densities has no effect. Guiner 
   * analysis at low Q give the R_G and the forward scattering at 
   * zero angle I(0).
   * 
   *    lnl(Q) = lnl(0) - R_G^2Q^2/3
   * 
   * where 
   *       Q = 4(pi)sin(theta/lamda)
   *       2theta = scattering angle
   *       lamda = wavelength
   * 
   * The above expression is valid in a QR_G range for extended 
   * rod-like particles. The relative I(0)/c values ( where
   *  c = sample concentration) for sample measurements in a
   * constant buffer for a single sample data session, gives the
   * relative masses of the protein(s) studied when referenced
   * against a standard.
   * 
   * see: 
   *     O.Glatter & O.Kratky, (1982). Editors of "Small angle
   *      X-ray Scattering, Academic Press, New York.
   *     O.Kratky. (1963). X-ray small angle scattering with
   *      substances of biological interest in diluted solutions.
   *      Prog. Biophys. Chem., 13, 105-173.
   *     G.D.Wignall & F.S.Bates, (1987). The small-angle approximation
   *      of X-ray and neutron scatter from rigid rods of non-uniform
   *      cross section and finite length. J.Appl. Crystallog., 18, 452-460.
   * 
   * If the structure is elongated, the mean radius of gyration
   * of the cross-sectional structure R_XS  and the mean cross sectional
   * intensity at zero angle [I(Q).Q]_Q->0 is obtained from
   * 
   *    ln[I(Q).Q] = ln[l(Q).(Q)]_Q->0 - ((R_XS)^2Q^2)/2
   */
  readonly min_mean_cross_sectional_radii_gyration_esd?: Maybe<Scalars['Float']>;
  /** The number of time frame solution scattering images used. */
  readonly num_time_frames?: Maybe<Scalars['Int']>;
  /**
   * The length (or range) of the protein sample under study.
   * If the solution structure is approximated as an elongated elliptical
   * cyclinder the the length L is determined from,
   * 
   *   L = sqrt [12( (R_G)^2  -  (R_XS)^2 ) ]   
   * 
   * The length should also be given by
   * 
   *   L = pi I(0) / [ I(Q).Q]_Q->0
   */
  readonly protein_length?: Maybe<Scalars['String']>;
  /** The pH value of the buffered sample. */
  readonly sample_pH?: Maybe<Scalars['Float']>;
  /** The beamline name used for the experiment */
  readonly source_beamline?: Maybe<Scalars['String']>;
  /** The instrumentation used on the beamline */
  readonly source_beamline_instrument?: Maybe<Scalars['String']>;
  /**
   * The general class of the radiation source.
   * 
   * Examples:
   * neutron source, synchrotron
   */
  readonly source_class?: Maybe<Scalars['String']>;
  /** The make, model, name or beamline of the source of radiation. */
  readonly source_type?: Maybe<Scalars['String']>;
  /**
   * The temperature in kelvins at which the experiment
   *  was conducted
   */
  readonly temperature?: Maybe<Scalars['Float']>;
  /**
   * The type of solution scattering experiment carried out
   * 
   * Allowable values:
   * modelling, neutron, x-ray
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityRcsbMacromolecularNamesCombined = {
  /** Combined list of macromolecular names. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Combined list of macromolecular names associated provenance code.
   * 
   *  ECO (https://github.com/evidenceontology/evidenceontology)
   */
  readonly provenance_code?: Maybe<Scalars['String']>;
  /**
   * Combined list of macromolecular names associated name source.
   * 
   * Allowable values:
   * PDB Preferred Name, PDB Synonym
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
};

export type EmVitrification = {
  /**
   * The temperature (in degrees Kelvin) of the sample just prior to vitrification.
   * 
   * Examples:
   * 298
   */
  readonly chamber_temperature?: Maybe<Scalars['Float']>;
  /**
   * This is the name of the cryogen.
   * 
   * Allowable values:
   * ETHANE, ETHANE-PROPANE, FREON 12, FREON 22, HELIUM, METHANE, NITROGEN, OTHER, PROPANE
   */
  readonly cryogen_name?: Maybe<Scalars['String']>;
  /**
   * Any additional details relating to vitrification.
   * 
   * Examples:
   * Vitrification carried out in argon atmosphere.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The humidity (%) in the vicinity of the vitrification process.
   * 
   * Examples:
   * 90
   */
  readonly humidity?: Maybe<Scalars['Float']>;
  /**
   * The value of _em_vitrification.id must uniquely identify
   *  the vitrification procedure.
   */
  readonly id: Scalars['String'];
  /**
   * The type of instrument used in the vitrification process.
   * 
   * Allowable values:
   * EMS-002 RAPID IMMERSION FREEZER, FEI VITROBOT MARK I, FEI VITROBOT MARK II, FEI VITROBOT MARK III, FEI VITROBOT MARK IV, GATAN CRYOPLUNGE 3, HOMEMADE PLUNGER, LEICA EM CPC, LEICA EM GP, LEICA KF80, LEICA PLUNGER, REICHERT-JUNG PLUNGER, SPOTITON, ZEISS PLUNGE FREEZER CRYOBOX
   */
  readonly instrument?: Maybe<Scalars['String']>;
  /**
   * The procedure for vitrification.
   * 
   * Examples:
   * plunge freezing
   */
  readonly method?: Maybe<Scalars['String']>;
  /** This data item is a pointer to _em_specimen.id */
  readonly specimen_id: Scalars['String'];
  /**
   * The vitrification temperature (in degrees Kelvin), e.g.,
   *   temperature of the plunge instrument cryogen bath.
   * 
   * Examples:
   * 90
   */
  readonly temp?: Maybe<Scalars['Float']>;
  /**
   * The length of time after an event effecting the sample that
   *  vitrification was induced and a description of the event.
   * 
   * Examples:
   * plunge 30 msec after spraying with effector
   */
  readonly time_resolved_state?: Maybe<Scalars['String']>;
};

export type RefineHist = {
  /**
   * The value of _refine_hist.cycle_id must uniquely identify a
   *  record in the REFINE_HIST list.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   */
  readonly cycle_id: Scalars['String'];
  /**
   * The lowest value for the interplanar spacings for the
   *  reflection data for this cycle of refinement. This is called
   *  the highest resolution.
   */
  readonly d_res_high?: Maybe<Scalars['Float']>;
  /**
   * The highest value for the interplanar spacings for the
   *  reflection data for this cycle of refinement. This is
   *  called the lowest resolution.
   */
  readonly d_res_low?: Maybe<Scalars['Float']>;
  /**
   * The number of solvent atoms that were included in the model at
   *  this cycle of the refinement.
   */
  readonly number_atoms_solvent?: Maybe<Scalars['Int']>;
  /**
   * The total number of atoms that were included in the model at
   *  this cycle of the refinement.
   */
  readonly number_atoms_total?: Maybe<Scalars['Int']>;
  /** Mean isotropic B-value for ligand molecules included in refinement. */
  readonly pdbx_B_iso_mean_ligand?: Maybe<Scalars['Float']>;
  /** Mean isotropic B-value for solvent molecules included in refinement. */
  readonly pdbx_B_iso_mean_solvent?: Maybe<Scalars['Float']>;
  /** Number of ligand atoms included in refinement */
  readonly pdbx_number_atoms_ligand?: Maybe<Scalars['Int']>;
  /** Number of nucleic atoms included in refinement */
  readonly pdbx_number_atoms_nucleic_acid?: Maybe<Scalars['Int']>;
  /** Number of protein atoms included in refinement */
  readonly pdbx_number_atoms_protein?: Maybe<Scalars['Int']>;
  /** Total number of polymer residues included in refinement. */
  readonly pdbx_number_residues_total?: Maybe<Scalars['Int']>;
  /**
   * This data item uniquely identifies a refinement within an entry.
   *  _refine_hist.pdbx_refine_id can be used to distinguish the results 
   *  of joint refinements.
   */
  readonly pdbx_refine_id: Scalars['String'];
};

export type PdbxVrptSummary = {
  /**
   * String for B_factor_type either "PARTIAL" or "FULL".
   * 
   * Allowable values:
   * FULL, PARTIAL
   */
  readonly B_factor_type?: Maybe<Scalars['String']>;
  /**
   * REFMAC scaling parameter as reported in log output line starting 'bulk solvent: scale'.
   * Example:            X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly Babinet_b?: Maybe<Scalars['Float']>;
  /**
   * REFMAC scaling parameter as reported in log output line starting 'bulk solvent: scale'.
   * Example: 
   * X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly Babinet_k?: Maybe<Scalars['Float']>;
  /**
   * The string "yes".
   * 
   * Allowable values:
   * yes
   */
  readonly CA_ONLY?: Maybe<Scalars['String']>;
  /**
   * The overall R-factor from a DCC recalculation of an electron density map.
   * Example:            Currently value is rounded to 2 decimal places.
   * X-ray entry specific, obtained from the DCC program.
   */
  readonly DCC_R?: Maybe<Scalars['Float']>;
  /** Either a decimal number or the string "NotAvailable". */
  readonly DCC_Rfree?: Maybe<Scalars['Float']>;
  /**
   * The pdbx_vrpt_software used by DCC to perform the recaluclation of the electron density maps.
   * &#160;Currently one of "CNS", "REFMAC" or "PHENIX".
   *  Example:            X-ray entry specific, obtained from the DCC program.
   */
  readonly DCC_refinement_program?: Maybe<Scalars['String']>;
  /**
   * The overall R factor from the EDS REFMAC calculation (no free set is used in this).
   * Example:            Currently value is rounded to 2 decimal places.
   * X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly EDS_R?: Maybe<Scalars['Float']>;
  /**
   * The data high resolution diffraction limit, in Angstroms, found in the input structure factor file.
   * Example:            X-ray entry specific, obtained in the eds step.
   */
  readonly EDS_resolution?: Maybe<Scalars['Float']>;
  /**
   * The data low resolution diffraction limit, in Angstroms, found in the input structure factor file.
   * Example:            X-ray entry specific, obtained in the eds step.
   */
  readonly EDS_resolution_low?: Maybe<Scalars['Float']>;
  /**
   * Fo,Fc correlation: The difference between the observed structure factors (Fo) and the 
   * calculated structure factors (Fc) measures the correlation between the PDB_model_num and the i
   * experimental data. 
   * Example:            X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly Fo_Fc_correlation?: Maybe<Scalars['Float']>;
  /**
   * Each reflection has an intensity (I) and an uncertainty in measurement 
   * (sigma(I)), so I/sigma(I) is the signal-to-noise ratio. This
   * ratio decreases at higher resolution. <I/sigma(I)> is the mean of individual I/sigma(I)
   * values. Value for outer resolution shell is given in parentheses. In case
   * structure factor amplitudes are deposited, Xtriage estimates the intensities
   * first and then calculates this metric. When intensities are available in the
   * deposited file, these are converted to amplitudes and then back to intensity
   * estimate before calculating the metric.  
   * Example            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly I_over_sigma?: Maybe<Scalars['String']>;
  /** Either a decimal number or the string "NotAvailable". */
  readonly PDB_R?: Maybe<Scalars['Float']>;
  /** Either a decimal number or the string "NotAvailable". */
  readonly PDB_Rfree?: Maybe<Scalars['Float']>;
  /**
   * Date in yyyy-mm-dd format when structure was deposited to the PDB.
   * Obtained from mmCIF table _database_PDB_rev item _database_PDB_rev.date_original 
   * Reports produced by the validation server or during the initial depositon process should not have this item.
   * If there is a difficulty parsing the item then "unknown" will be given.
   */
  readonly PDB_deposition_date?: Maybe<Scalars['Date']>;
  /** Either a decimal number or the string "NotAvailable". */
  readonly PDB_resolution?: Maybe<Scalars['Float']>;
  /** Either a decimal number or the string "NotAvailable". */
  readonly PDB_resolution_low?: Maybe<Scalars['Float']>;
  /**
   * Date in yyyy-mm-dd format when the structure was last revised by PDB.
   * Obtained from mmCIF table _database_PDB_rev item _database_PDB_rev.date
   * Reports produced by the validation server or during the initial depositon process should not have this item.
   * If there is a difficulty parsing the item then "unknown" will be given.
   */
  readonly PDB_revision_date?: Maybe<Scalars['Date']>;
  /**
   * The last highest number that appears in mmCIF item _database_PDB_rev.num.
   * Data items in the DATABASE_PDB_REV category record details about the history of the data block as archived by the Protein Data Bank (PDB).
   * If the input mmCIF coordinate file lacks the information then a value of -1 is supplied.
   */
  readonly PDB_revision_number?: Maybe<Scalars['Float']>;
  /**
   * The MolProbity conformer-match quality parameter for RNA structures.
   * Low values are worse.
   * Example:           Specific to structures that contain RNA polymers.
   */
  readonly RNA_suiteness?: Maybe<Scalars['Float']>;
  /**
   * Result of absolute likelihood based Wilson scaling, 
   * The anisotropic B value of the data is determined using a likelihood based approach. 
   * The resulting B tensor is reported, the 3 diagonal values are given first, followed
   * by the 3 off diagonal values.
   * A large spread in (especially the diagonal) values indicates anisotropy. 
   * Example: 
   * X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly Wilson_B_aniso?: Maybe<Scalars['String']>;
  /**
   * An estimate of the overall B-value of the structure, calculated from the diffraction data. 
   * Units Angstroms squared.
   * It serves as an indicator of the degree of order in the crystal and the value is usually 
   * not hugely different from the average B-value calculated from the model.
   * Example:            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly Wilson_B_estimate?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_DCC_Rfree?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_RNA_suiteness?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_clashscore?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly absolute_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>;
  /**
   * The number of acentric reflections that Xtriage identifies as outliers on the basis 
   * of Wilson statistics. Note that if pseudo translational symmetry is present, 
   * a large number of 'outliers' will be present.
   * Example:            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly acentric_outliers?: Maybe<Scalars['Int']>;
  /**
   * The overall root mean square of the Z-score for deviations of bond angles in comparison to 
   * "standard geometry" made using the MolProbity dangle program.
   * Standard geometry parameters are taken from Engh and Huber (2001) and Parkinson et al. (1996).
   * This value is for all chains in the structure.
   */
  readonly angles_RMSZ?: Maybe<Scalars['Float']>;
  /**
   * The steps that were attempted by the validation pipeline software. 
   * A step typically involves running a 3rd party validation tool, for instance "mogul"
   * Each step that was successfully completed will result in a pdbx_vrpt_software element in the pdbx_vrpt_sotfware_notused list.
   */
  readonly attempted_validation_steps?: Maybe<Scalars['String']>;
  /**
   * The overall root mean square of the Z-score for deviations of bond lengths in comparison to 
   * "standard geometry" made using the MolProbity dangle program.
   * Standard geometry parameters are taken from Engh and Huber (2001) and Parkinson et al. (1996).
   * This value is for all chains in the structure.
   */
  readonly bonds_RMSZ?: Maybe<Scalars['Float']>;
  /**
   * REFMAC scaling parameter as reported in log output line starting 
   * 'Partial structure    1: scale.'
   * Example: 
   * X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly bulk_solvent_b?: Maybe<Scalars['Float']>;
  /**
   * REFMAC scaling parameter as reported in log output line starting 
   * 'Partial structure    1: scale.'
   * Example:            X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly bulk_solvent_k?: Maybe<Scalars['Float']>;
  /**
   * The version of CCP4 suite pdbx_vrpt_sotfware_notused used in the analysis.
   * Example:            X-ray entry specific, obtained from the eds step.
   */
  readonly ccp4version?: Maybe<Scalars['String']>;
  /**
   * The number of centric reflections that Xtriage identifies as outliers.
   * Example:            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly centric_outliers?: Maybe<Scalars['Float']>;
  /**
   * Overall completeness of the chemical shift assignments for the well-defined 
   * regions of the structure.
   */
  readonly chemical_shift_completeness?: Maybe<Scalars['Float']>;
  /**
   * Overall completeness of the chemical shift assignments for the full 
   * macromolecule or complex as suggested by the molecular description of an entry
   * (whether some portion of it is modelled or not).
   */
  readonly chemical_shift_completeness_full_length?: Maybe<Scalars['Float']>;
  /**
   * The filename for the input chemical shifts file given to the validation pipeline.
   * Not reported for runs at the annotation or release stage.
   */
  readonly chemical_shifts_input_filename?: Maybe<Scalars['String']>;
  /**
   * This score is derived from the number of pairs of atoms in the PDB_model_num that are unusually close to each other. 
   * It is calculated by the MolProbity pdbx_vrpt_software and expressed as the number or such clashes per thousand atoms.
   * For structures determined by NMR the clashscore value here will only consider label_atom_id pairs in the 
   * well-defined (core) residues from ensemble analysis.
   */
  readonly clashscore?: Maybe<Scalars['Float']>;
  /** Only given for structures determined by NMR. The MolProbity pdbx_vrpt_clashes score for all label_atom_id pairs. */
  readonly clashscore_full_length?: Maybe<Scalars['Float']>;
  /**
   * The filename for the input mmCIF coordinate file given to the validation pipeline.
   * Not reported for runs at the annotation or release stage.
   */
  readonly coordinates_input_filename?: Maybe<Scalars['String']>;
  /**
   * Diagnostic message from the wrapper of Cyrange software which identifies the 
   * well-defined cores (domains) of NMR protein structures.
   */
  readonly cyrange_error?: Maybe<Scalars['String']>;
  /** Total number of well-defined cores (domains) identified by Cyrange */
  readonly cyrange_number_of_domains?: Maybe<Scalars['Int']>;
  /** Reference for the Cyrange software. */
  readonly cyrange_version?: Maybe<Scalars['String']>;
  /**
   * The ratio (Bmax &#8209; Bmin) / Bmean where Bmax, Bmin and Bmean are computed from the B-values 
   * associated with the principal axes of the anisotropic thermal ellipsoid. 
   * This ratio is usually less than 0.5; for only 1% of PDB entries it is more than 1.0 (Read et al., 2011).
   * Example:            X-ray entry specific, obtained from the Xtriage program.
   */
  readonly data_anisotropy?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly data_completeness?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-DCC_Rfree". 
   * Note that the "high_resol_relative_percentile_DCC_Rfree" value is numerically smaller than the 
   * corresponding low-* value.
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-RNAsuiteness". 
   * Note that the "high_resol_relative_percentile_RNA_suiteness" value is numerically smaller than the 
   * corresponding low value.
   * Example:            Specific to entries that contain RNA polymers (and have a RNA_suiteness attribute),
   * and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-clashscore". 
   * Note that the "high_resol_relative_percentile_clashscore" value is numerically smaller than the 
   * corresponding low value.
   * Example:            Specific to that have a clashscore attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_clashscore?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-RSRZ-outliers". 
   * Note that the "high_resol_relative_percentile_percent_RSRZ_outliers" value is numerically smaller than the 
   * corresponding low-* value.
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-rama-outliers". 
   * Note that the "high_resol_relative_percentile_percent_ramachandran_outliers" value is numerically smaller than the 
   * corresponding low value.
   * Example:            Specific to structures that have a percent_ramachandran_outliers attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-rota-outliers". 
   * Note that the "high_resol_relative_percentile_percent_rotamer_outliers" value is numerically smaller than the 
   * corresponding low value.
   * Example:            Specific to that have a percent_rotamer_outliers attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly high_resol_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>;
  /**
   * The string "yes".
   * 
   * Allowable values:
   * yes
   */
  readonly ligands_for_buster_report?: Maybe<Scalars['String']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-DCC_Rfree". 
   * Note that the "low_resol_relative_percentile_DCC_Rfree" value is numerically greater than the 
   * corresponding high value.
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-RNAsuiteness". 
   * Note that the "low_resol_relative_percentile_RNA_suiteness" value is numerically greater than the 
   * corresponding high value.
   * Example:            Specific to entries that contain RNA polymers (and have a RNA_suiteness attribute),
   * and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-clashscore". 
   * Note that the "low_resol_relative_percentile_clashscore" value is numerically greater than the 
   * corresponding high value.
   * Example:            Specific to that have a clashscore attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_clashscore?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-RSRZ-outliers". 
   * Note that the "low_resol_relative_percentile_percent_RSRZ_outliers" value is numerically greater than the 
   * corresponding high value.
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-rama-outliers". 
   * Note that the "low_resol_relative_percentile_percent_ramachandran_outliers" value is numerically greater than the 
   * corresponding high value.
   * Example:            Specific to that have a percent_ramachandran_outliers attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>;
  /**
   * The resolution bin limit in Angstroms for PDB depositions used in the comparison when calculating 
   * the attribute "relative-percentile-percent-rota-outliers". 
   * Note that the "low_resol_relative_percentile_percent_rotamer_outliers" value is numerically greater than the 
   * corresponding high value.
   * Example:            Specific to that have a percent_rotamer_outliers attribute and have been determined by X-ray crystallography.
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly low_resol_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>;
  /**
   * For each Cyrange well-defined core ("cyrange_domain") the id of the PDB_model_num which is most 
   * similar to other models as measured by pairwise RMSDs over the domain. 
   * For the whole entry ("Entry"), the medoid PDB_model_num of the largest core is taken as an overall
   * representative of the structure.
   */
  readonly medoid_model?: Maybe<Scalars['Int']>;
  /**
   * A flag indicating if all models in the NMR ensemble contain the exact 
   * same atoms ("True") or if the models differ in this respect ("False").
   */
  readonly nmr_models_consistency_flag?: Maybe<Scalars['String']>;
  /** Diagnostic message from the wrapper of NMRClust software which clusters NMR models. */
  readonly nmrclust_error?: Maybe<Scalars['String']>;
  /** Total number of clusters in the NMR ensemble identified by NMRClust. */
  readonly nmrclust_number_of_clusters?: Maybe<Scalars['Int']>;
  /**
   * Number of models analysed by NMRClust - should in almost all cases be the
   * same as the number of models in the NMR ensemble.
   */
  readonly nmrclust_number_of_models?: Maybe<Scalars['Int']>;
  /** Number of models that do not belong to any cluster as deemed by NMRClust. */
  readonly nmrclust_number_of_outliers?: Maybe<Scalars['Int']>;
  /** Overall representative PDB_model_num of the NMR ensemble as identified by NMRClust. */
  readonly nmrclust_representative_model?: Maybe<Scalars['Int']>;
  /** Reference for the NMRClust software. */
  readonly nmrclust_version?: Maybe<Scalars['String']>;
  /**
   * The string "yes".
   * 
   * Allowable values:
   * yes
   */
  readonly no_ligands_for_buster_report?: Maybe<Scalars['String']>;
  /**
   * The string "yes".
   * 
   * Allowable values:
   * yes
   */
  readonly no_ligands_for_mogul?: Maybe<Scalars['String']>;
  /**
   * Will be set to "true" if no property was found to do percentile analysis on. 
   * Please note that currently due to a bug that this attribute is "true" erronously for NMR structures.
   */
  readonly no_percentile_property?: Maybe<Scalars['String']>;
  /**
   * This is the number of hydrogen atoms added and optimized by the MolProbity reduce pdbx_vrpt_software as part of the
   * all-atom clashscore.
   */
  readonly num_H_reduce?: Maybe<Scalars['Float']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "absolute-percentile-DCC_Rfree".
   * Example: 
   * X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_DCC_Rfree?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "absolute-percentile-RNAsuiteness".
   * Example:            Specific to entries that contain RNA polymers (and have a RNA_suiteness attribute).
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_RNA_suiteness?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "absolute-percentile-clashscore"
   * Example:             Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_clashscore?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "absolute-percentile-percent-RSRZ-outliers".
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "absolute_percentile_percent_ramachandran_outliers" 
   * Example:            Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute
   * "absolute-percentile-percent-rota-outliers"
   * Example:             Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_absolute_percentile_percent_rotamer_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-DCC_Rfree".
   * Example: 
   * X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_DCC_Rfree?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-RNAsuiteness".
   * Example:            Specific to entries that contain RNA polymers (and have a RNA_suiteness attribute).
   * Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_RNA_suiteness?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-clashscore"
   * Example:            Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_clashscore?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-percent-RSRZ-outliers".
   * Example:            X-ray entry specific, produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-percent-rama-outliers"
   * Example:            Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of PDB depositions used in the comparison when calculating the attribute 
   * "relative-percentile-percent-rota-outliers"
   * Example:            Produced by the percentiles step of the validation pipeline software.
   */
  readonly num_PDBids_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Int']>;
  /**
   * The number of bond angless compared to "standard geometry" made using the MolProbity dangle program.
   * Standard geometry parameters are taken from Engh and Huber (2001) and Parkinson et al. (1996).
   * This value is for all chains in the structure.
   */
  readonly num_angles_RMSZ?: Maybe<Scalars['Int']>;
  /**
   * The number of bond lengths compared to "standard geometry" made using the MolProbity dangle program.
   * Standard geometry parameters are taken from Engh and Huber (2001) and Parkinson et al. (1996).
   * This value is for all chains in the structure.
   */
  readonly num_bonds_RMSZ?: Maybe<Scalars['Int']>;
  /**
   * The number of reflections in the free set as defined in the input structure factor file supplied to the validation pipeline. 
   * example:            X-ray entry specific, obtained from the DCC program.
   */
  readonly num_free_reflections?: Maybe<Scalars['Int']>;
  /**
   * The number of Miller Indices reported by the Xtriage program. This should be the same as the
   * number of _refln in the input structure factor file.
   * Example:            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly num_miller_indices?: Maybe<Scalars['Int']>;
  /** Reference for the PANAV software. */
  readonly panav_version?: Maybe<Scalars['String']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_RSRZ_outliers?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_free_reflections?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_ramachandran_outliers?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_ramachandran_outliers_full_length?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_rotamer_outliers?: Maybe<Scalars['Float']>;
  /** A percentage, Normally percent proportion of the total number. Between 0% and 100%. */
  readonly percent_rotamer_outliers_full_length?: Maybe<Scalars['Float']>;
  /**
   * The percentile bins that this structure would contribute to in a recalculation of
   * percentile. The string is a comma separated list.
   * Example: An X-ray entry with a resolution of 1.8 Angstroms, 
   * that would contribute to bins all structures, what ever bin 1.8 Angstrom resolution is in and 
   * the all xray bin.
   * Example #2: An EM entry that would contribute to all and em structures bins.
   * Example #3: A NMR entry hat would contribute to all and nmr structures bins.
   */
  readonly percentilebins?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Comma separated list of the _entity.id's for the molecular entities that are present in the structure
   * and have _entity_poly.type indicating that they are protein, RNA or DNA: that is in the list
   * 'polypeptide(L)', 'polypeptide(D)', 'polyribonucleotide, 'polydeoxyribonucleotide'  or
   * 'polydeoxyribonucleotide/polyribonucleotide hybrid'.
   * Normally the entity.id's are integer numbers but not necessarily so. 
   * Example
   */
  readonly protein_DNA_RNA_entities?: Maybe<Scalars['String']>;
  /** Version and reference of the RCI software */
  readonly rci_version?: Maybe<Scalars['String']>;
  /**
   * The filename for the input mmCIF format reflection file given to the validation pipeline.
   * Not reported for runs at the annotation or release stage.
   */
  readonly reflections_input_filename?: Maybe<Scalars['String']>;
  /**
   * Version of the REFMAC pdbx_vrpt_software used in the EDS step.                               
   * Example:           X-ray entry specific, obtained in the eds step from REFMAC calculation.
   */
  readonly refmac_version?: Maybe<Scalars['String']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_clashscore?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>;
  /**
   * Structures are judged in comparison to previously deposited PDB entries. 
   * The comparison is carried out by calculation of the percentile rank, i.e. the percentage of entries 
   * that are equal or poorer than this structure in terms of a quality indicator.
   * Percentile ranks range from 0 (the worst) to 100 (the best).
   */
  readonly relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>;
  /**
   * The date, time and time-zone that the validation XML file was created. 
   * The string will be formatted like "Feb  7, 2017 -- 12:32 pm GMT".
   */
  readonly report_creation_date?: Maybe<Scalars['String']>;
  /**
   * The data high resolution diffraction limit, in Angstroms, obtained from cif item 
   * _reflns.d_resolution_high.
   * X-ray entry specific.
   */
  readonly resol_high_from_reflectionsfile?: Maybe<Scalars['Float']>;
  /**
   * The data low resolution diffraction limit, in Angstroms, obtained from cif item 
   * _reflns.d_resolution_low.
   * X-ray entry specific.
   */
  readonly resol_low_from_reflectionsfile?: Maybe<Scalars['Float']>;
  /**
   * This is a comma separated list of the residue types whose bond lengths and bond angles have 
   * not been checked against "standard geometry" using the MolProbity dangle program.
   * Standard geometry parameters are taken from Engh and Huber (2001) and Parkinson et al. (1996)
   */
  readonly restypes_notchecked_for_bond_angle_geometry?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Version of the software for chemical shift outlier detection - currently 
   * same as revision number of the validation pipeline.
   */
  readonly shiftchecker_version?: Maybe<Scalars['String']>;
  /**
   * A sentence giving the result of Xtriage&#8217;s analysis on translational NCS.
   * Example: largest off-origin peak in the Patterson function is 8.82% of the height of the origin peak. No significant pseudotranslation is detected."
   * X-ray entry specific, obtained from the Xtriage program.
   */
  readonly trans_NSC?: Maybe<Scalars['String']>;
  /**
   * Padilla and Yeates twinning parameter <|L|>. 
   * Theoretical values is 0.5 in the untwinned case, and 0.375 in the perfectly twinned case.
   * Example:            X-ray entry specific, obtained from the Xtriage program.
   */
  readonly twin_L?: Maybe<Scalars['Float']>;
  /**
   * Padilla and Yeates twinning parameter <|L**2|>. 
   * Theoretical values is 0.333 in the untwinned case, and 0.2 in the perfectly twinned case.
   * Example:            X-ray entry specific, obtained from the Xtriage program.
   */
  readonly twin_L2?: Maybe<Scalars['Float']>;
  /**
   * Estimated twinning fraction for operators as identified by Xtriage. A semicolon separated
   * list of operators with fractions is givens 
   * Example:            X-ray entry specific, obtained from the Xtriage program.
   */
  readonly twin_fraction?: Maybe<Scalars['String']>;
  /**
   * The mmCIF item names of the _refln columns used as input to the Xtriage program.
   * Example            X-ray entry specific, calculated by Phenix Xtriage program.
   */
  readonly xtriage_input_columns?: Maybe<Scalars['String']>;
};

export type RcsbUniprotContainerIdentifiersReferenceSequenceIdentifiers = {
  /** Reference database accession code */
  readonly database_accession?: Maybe<Scalars['String']>;
  /** Reference database identifier for the sequence isoform */
  readonly database_isoform?: Maybe<Scalars['String']>;
  /** Reference database name */
  readonly database_name?: Maybe<Scalars['String']>;
  /** Source of the reference database assignment */
  readonly provenance_source?: Maybe<Scalars['String']>;
};

export type RcsbEntitySourceOrganismRcsbGeneName = {
  /**
   * A code indicating the provenance of the source organism details for the entity
   * 
   * Allowable values:
   * PDB Primary Data, UniProt
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /** Gene name. */
  readonly value?: Maybe<Scalars['String']>;
};

export type RcsbPolymerStructConn = {
  readonly connect_partner?: Maybe<RcsbPolymerStructConnConnectPartner>;
  readonly connect_target?: Maybe<RcsbPolymerStructConnConnectTarget>;
  /**
   * The connection type.
   * 
   * Allowable values:
   * covalent bond, covalent modification of a nucleotide base, covalent modification of a nucleotide phosphate, covalent modification of a nucleotide sugar, covalent residue modification, disulfide bridge, hydrogen bond, ionic interaction, metal coordination, mismatched base pairs
   */
  readonly connect_type?: Maybe<Scalars['String']>;
  /** A description of special details of the connection. */
  readonly description?: Maybe<Scalars['String']>;
  /** Distance value for this contact. */
  readonly dist_value?: Maybe<Scalars['Float']>;
  /**
   * The value of _rcsb_polymer_struct_conn.id is an identifier for connection.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * The value of _rcsb_polymer_struct_conn.id must uniquely identify a record in
   *  the rcsb_polymer_struct_conn list.
   */
  readonly ordinal_id: Scalars['Int'];
  /**
   * The chemical or structural role of the interaction
   * 
   * Allowable values:
   * C-Mannosylation, N-Glycosylation, O-Glycosylation
   */
  readonly role?: Maybe<Scalars['String']>;
  /**
   * The chemical bond order associated with the specified atoms in
   *  this contact.
   * 
   * Allowable values:
   * doub, quad, sing, trip
   */
  readonly value_order?: Maybe<Scalars['String']>;
};

export type EmHelicalEntity = {
  /**
   * The angular rotation per helical subunit in degrees.
   * 
   * Examples:
   * -34.616000
   */
  readonly angular_rotation_per_subunit?: Maybe<Scalars['Float']>;
  /**
   * The axial rise per subunit in the helical assembly.
   * 
   * Examples:
   * 17.400000
   */
  readonly axial_rise_per_subunit?: Maybe<Scalars['Float']>;
  /**
   * Symmetry of the helical axis, either cyclic (Cn) or dihedral (Dn), where n>=1.
   * 
   * Examples:
   * C1, D2, C7
   */
  readonly axial_symmetry?: Maybe<Scalars['String']>;
  /**
   * Any other details regarding the helical assembly
   * 
   * Examples:
   * Dihedral symmetry
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The value of _em_helical_entity.id must uniquely identify
   *  a set of the filament parameters for this assembly component.
   */
  readonly id: Scalars['String'];
  /**
   * The value of _em_helical_entity.reconstruction_id identifies a particular reconstruction.
   * 
   *  This data item is a pointer to _em_image_processing.id.
   */
  readonly image_processing_id: Scalars['String'];
};

export type CoreEntityAlignmentsCoreEntityIdentifiers = {
  readonly entity_id: Scalars['String'];
  readonly entry_id: Scalars['String'];
};

export type PdbxStructAssembly = {
  /**
   * A description of special aspects of the macromolecular assembly.
   * 
   * Examples:
   * The icosahedral virus particle.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_struct_assembly.id must uniquely identify a record in
   *  the PDBX_STRUCT_ASSEMBLY list.
   */
  readonly id: Scalars['String'];
  /**
   * Provides details of the method used to determine or 
   *  compute the assembly.
   */
  readonly method_details?: Maybe<Scalars['String']>;
  /** The number of polymer molecules in the assembly. */
  readonly oligomeric_count?: Maybe<Scalars['Int']>;
  /**
   * Provides the details of the oligomeric state of the assembly.
   * 
   * Examples:
   * monomer, octameric, tetradecameric, eicosameric, 21-meric, 60-meric, 180-meric, helical
   */
  readonly oligomeric_details?: Maybe<Scalars['String']>;
  /**
   * Candidate macromolecular assembly.
   * 
   *  Excludes the following cases classified in pdbx_struct_asembly.details:
   * 
   *  'crystal asymmetric unit', 'crystal asymmetric unit, crystal frame', 'helical asymmetric unit',
   *  'helical asymmetric unit, std helical frame','icosahedral 23 hexamer', 'icosahedral asymmetric unit',
   *  'icosahedral asymmetric unit, std point frame','icosahedral pentamer', 'pentasymmetron capsid unit',
   *  'point asymmetric unit', 'point asymmetric unit, std point frame','trisymmetron capsid unit',
   *   and 'deposited_coordinates'.
   * 
   * Allowable values:
   * N, Y
   */
  readonly rcsb_candidate_assembly?: Maybe<Scalars['String']>;
  /**
   * A filtered description of the macromolecular assembly.
   * 
   * Allowable values:
   * author_and_software_defined_assembly, author_defined_assembly, software_defined_assembly
   */
  readonly rcsb_details?: Maybe<Scalars['String']>;
};

export type RcsbRepositoryHoldingsCurrentEntryContainerIdentifiers = {
  /** The assembly id codes. */
  readonly assembly_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The PDB entry accession code. */
  readonly entry_id: Scalars['String'];
  /** The RCSB entry identifier. */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  /** Identifier for the current data exchange status record. */
  readonly update_id?: Maybe<Scalars['String']>;
};

export type CorePubmed = {
  /** Unique integer value assigned to each PubMed record. */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  /** A concise, accurate and factual mini-version of the paper contents. */
  readonly rcsb_pubmed_abstract_text?: Maybe<Scalars['String']>;
  readonly rcsb_pubmed_affiliation_info?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Unique integer value assigned to each PubMed Central record. */
  readonly rcsb_pubmed_central_id?: Maybe<Scalars['String']>;
  readonly rcsb_pubmed_container_identifiers: RcsbPubmedContainerIdentifiers;
  /** Persistent identifier used to provide a link to an article location on the Internet. */
  readonly rcsb_pubmed_doi?: Maybe<Scalars['String']>;
  readonly rcsb_pubmed_mesh_descriptors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  readonly rcsb_pubmed_mesh_descriptors_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPubmedMeshDescriptorsLineage>>>;
};

export type ChemComp = {
  /**
   * The formula for the chemical component. Formulae are written
   *  according to the following rules:
   * 
   *  (1) Only recognized element symbols may be used.
   * 
   *  (2) Each element symbol is followed by a 'count' number. A count
   *     of '1' may be omitted.
   * 
   *  (3) A space or parenthesis must separate each cluster of
   *     (element symbol + count), but in general parentheses are
   *     not used.
   * 
   *  (4) The order of elements depends on whether carbon is
   *     present or not. If carbon is present, the order should be:
   *     C, then H, then the other elements in alphabetical order
   *     of their symbol. If carbon is not present, the elements
   *     are listed purely in alphabetic order of their symbol. This
   *     is the 'Hill' system used by Chemical Abstracts.
   * 
   * Examples:
   * C18 H19 N7 O8 S
   */
  readonly formula?: Maybe<Scalars['String']>;
  /** Formula mass of the chemical component. */
  readonly formula_weight?: Maybe<Scalars['Float']>;
  /**
   * The value of _chem_comp.id must uniquely identify each item in
   *  the CHEM_COMP list.
   * 
   *  For protein polymer entities, this is the three-letter code for
   *  the amino acid.
   * 
   *  For nucleic acid polymer entities, this is the one-letter code
   *  for the base.
   * 
   * Examples:
   * ALA, VAL, DG, C
   */
  readonly id: Scalars['String'];
  /**
   * The identifier for the parent component of the nonstandard
   *  component. May be be a comma separated list if this component
   *  is derived from multiple components.
   * 
   *  Items in this indirectly point to _chem_comp.id in 
   *  the CHEM_COMP category.
   */
  readonly mon_nstd_parent_comp_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The full name of the component.
   * 
   * Examples:
   * alanine, valine, adenine, cytosine
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * For standard polymer components, the one-letter code for
   *  the component.   For non-standard polymer components, the 
   *  one-letter code for parent component if this exists;
   *  otherwise, the one-letter code should be given as 'X'.
   * 
   *  Components that derived from multiple parents components 
   *  are described by a sequence of one-letter-codes.
   * 
   * Examples:
   * A, B, R, N, D, C, Q, E, Z, G, H, I, L, K, M, F, P, S, T, W, Y, V, U, O, X
   */
  readonly one_letter_code?: Maybe<Scalars['String']>;
  /**
   * A preliminary classification used by PDB to indicate
   *  that the chemistry of this component while described
   *  as clearly as possible is still ambiguous.  Software
   *  tools may not be able to process this component
   *  definition.
   */
  readonly pdbx_ambiguous_flag?: Maybe<Scalars['String']>;
  /**
   * The net integer charge assigned to this component. This is the
   *  formal charge assignment normally found in chemical diagrams.
   */
  readonly pdbx_formal_charge?: Maybe<Scalars['Int']>;
  /** Date component was added to database. */
  readonly pdbx_initial_date?: Maybe<Scalars['Date']>;
  /** Date component was last modified. */
  readonly pdbx_modified_date?: Maybe<Scalars['Date']>;
  /**
   * This data item identifies the deposition site that processed
   *  this chemical component defintion.
   * 
   * Allowable values:
   * EBI, PDBC, PDBE, PDBJ, RCSB
   */
  readonly pdbx_processing_site?: Maybe<Scalars['String']>;
  /**
   * This data item holds the current release status for the component.
   * 
   * Allowable values:
   * DEL, HOLD, HPUB, OBS, REF_ONLY, REL
   */
  readonly pdbx_release_status?: Maybe<Scalars['String']>;
  /**
   * Identifies the _chem_comp.id of the component that
   *  has replaced this component.
   * 
   * Examples:
   * q11, tvx
   */
  readonly pdbx_replaced_by?: Maybe<Scalars['String']>;
  /**
   * Identifies the _chem_comp.id's of the components
   *  which have been replaced by this component.
   *  Multiple id codes should be separated by commas.
   * 
   * Examples:
   * q11, tvx,atv
   */
  readonly pdbx_replaces?: Maybe<Scalars['String']>;
  /**
   * The list of subcomponents contained in this component.
   * 
   * Examples:
   * TSM DPH HIS CHF EMR
   */
  readonly pdbx_subcomponent_list?: Maybe<Scalars['String']>;
  /**
   * For standard polymer components, the common three-letter code for
   *  the component.   Non-standard polymer components and non-polymer
   *  components are also assigned three-letter-codes.
   * 
   *  For ambiguous polymer components three-letter code should 
   *  be given as 'UNK'.  Ambiguous ions are assigned the code 'UNX'.
   *  Ambiguous non-polymer components are assigned the code 'UNL'.
   * 
   * Examples:
   * ALA, ARG, ASN, ASP, ASX, CYS, GLN, GLU, GLY, GLX, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, 1MA, 5MC, OMC, 1MG, 2MG, M2G, 7MG, 0MG, H2U, 5MU, PSU, ACE, FOR, HOH, UNK
   */
  readonly three_letter_code?: Maybe<Scalars['String']>;
  /**
   * For standard polymer components, the type of the monomer.
   *  Note that monomers that will form polymers are of three types:
   *  linking monomers, monomers with some type of N-terminal (or 5')
   *  cap and monomers with some type of C-terminal (or 3') cap.
   * 
   * Allowable values:
   * D-beta-peptide, C-gamma linking, D-gamma-peptide, C-delta linking, D-peptide COOH carboxy terminus, D-peptide NH3 amino terminus, D-peptide linking, D-saccharide, D-saccharide 1,4 and 1,4 linking, D-saccharide 1,4 and 1,6 linking, D-saccharide, alpha linking, D-saccharide, beta linking, DNA OH 3 prime terminus, DNA OH 5 prime terminus, DNA linking, L-DNA linking, L-RNA linking, L-beta-peptide, C-gamma linking, L-gamma-peptide, C-delta linking, L-peptide COOH carboxy terminus, L-peptide NH3 amino terminus, L-peptide linking, L-saccharide, L-saccharide 1,4 and 1,4 linking, L-saccharide 1,4 and 1,6 linking, L-saccharide, alpha linking, L-saccharide, beta linking, RNA OH 3 prime terminus, RNA OH 5 prime terminus, RNA linking, non-polymer, other, peptide linking, peptide-like, saccharide
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbUniprotExternalReference = {
  readonly provenance_source?: Maybe<Scalars['String']>;
  readonly reference_id?: Maybe<Scalars['String']>;
  /**
   * Allowable values:
   * IMPC, GTEX, PHAROS
   */
  readonly reference_name?: Maybe<Scalars['String']>;
};

export type PdbxReferenceEntityPoly = {
  /** The database code for this source information */
  readonly db_code?: Maybe<Scalars['String']>;
  /** The database name for this source information */
  readonly db_name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_entity_poly.prd_id is a reference
   * 	       _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_entity_poly.ref_entity_id is a reference
   *  to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
   */
  readonly ref_entity_id: Scalars['String'];
  /**
   * The type of the polymer.
   * 
   * Allowable values:
   * nucleic-acid-like, peptide-like, polysaccharide-like
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type CurrentEntry = {
  /** The RCSB entry identifier. */
  readonly rcsb_id: Scalars['String'];
  readonly rcsb_repository_holdings_current?: Maybe<RcsbRepositoryHoldingsCurrent>;
  readonly rcsb_repository_holdings_current_entry_container_identifiers?: Maybe<RcsbRepositoryHoldingsCurrentEntryContainerIdentifiers>;
};

export type RcsbStructSymmetry = {
  readonly clusters: ReadonlyArray<Maybe<RcsbStructSymmetryClusters>>;
  /**
   * The granularity at which the symmetry calculation is performed. In 'Global Symmetry' all polymeric subunits in assembly are used. In 'Local Symmetry' only a subset of polymeric subunits is considered. In 'Pseudo Symmetry' the threshold for subunits similarity is relaxed.
   * 
   * Allowable values:
   * Global Symmetry, Pseudo Symmetry, Local Symmetry
   */
  readonly kind: Scalars['String'];
  /** Oligomeric state refers to a composition of polymeric subunits in quaternary structure. Quaternary structure may be composed either exclusively of several copies of identical subunits, in which case they are termed homo-oligomers, or alternatively by at least one copy of different subunits (hetero-oligomers). Quaternary structure composed of a single subunit is denoted as 'Monomer'. */
  readonly oligomeric_state: Scalars['String'];
  /** The orientation of the principal rotation (symmetry) axis. */
  readonly rotation_axes?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetryRotationAxes>>>;
  /** Each type of different subunits is assigned a latter. The number of equivalent subunits is added as a coefficient after each letter (except 1 which is not added explicitly). */
  readonly stoichiometry: ReadonlyArray<Maybe<Scalars['String']>>;
  /** Symmetry symbol refers to point group or helical symmetry of identical polymeric subunits in Schönflies notation. Contains point group symbol (e.g., C2, C5, D2, T, O, I) or H for helical symmetry. */
  readonly symbol: Scalars['String'];
  /**
   * Symmetry type refers to point group or helical symmetry of identical polymeric subunits. Contains point group types (e.g. Cyclic, Dihedral) or Helical for helical symmetry.
   * 
   * Allowable values:
   * Asymmetric, Cyclic, Dihedral, Tetrahedral, Octahedral, Icosahedral, Helical
   */
  readonly type: Scalars['String'];
};

export type RcsbStructSymmetryLineage = {
  /** Hierarchy depth. */
  readonly depth?: Maybe<Scalars['Int']>;
  /** A unique identifier automatically assigned to the symmetry term. */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * A human-readable term describing protein symmetry.
   * 
   * Examples:
   * Asymmetric, Global Symmetry, C1, Hetero 3-mer
   */
  readonly name?: Maybe<Scalars['String']>;
};

export type PdbxReferenceMoleculeAnnotation = {
  /**
   * The value of _pdbx_reference_molecule_annotation.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_list.family_prd_id in category PDBX_REFERENCE_MOLECULE_FAMILY_LIST.
   */
  readonly family_prd_id: Scalars['String'];
  /** This data item distinguishes anotations for this entity. */
  readonly ordinal: Scalars['Int'];
  /**
   * This data item is a pointer to _pdbx_reference_molecule.prd_id in the 
   *  PDB_REFERENCE_MOLECULE category.
   */
  readonly prd_id?: Maybe<Scalars['String']>;
  /**
   * The source of the annoation for this entity.
   * 
   * Examples:
   * depositor provided, from UniProt Entry P200311
   */
  readonly source?: Maybe<Scalars['String']>;
  /**
   * Text describing the annotation for this entity.
   * 
   * Examples:
   * antigen binding, glucose transporter activity
   */
  readonly text?: Maybe<Scalars['String']>;
  /**
   * Type of annotation for this entity.
   * 
   * Examples:
   * Function, Use, Pharmacology, Mechanism_of_Action, Biological_Activity, Inhibitor_Class, Therapeutic_Category, Research_Use, Other_annotation
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxReferenceMoleculeRelatedStructures = {
  /** A link to related reference information in the citation category. */
  readonly citation_id?: Maybe<Scalars['String']>;
  /**
   * The database accession code for the related structure reference.
   * 
   * Examples:
   * 143108
   */
  readonly db_accession?: Maybe<Scalars['String']>;
  /**
   * The database identifier code for the related structure reference.
   * 
   * Examples:
   * QEFHUE
   */
  readonly db_code?: Maybe<Scalars['String']>;
  /**
   * The database name for the related structure reference.
   * 
   * Examples:
   * CCDC
   */
  readonly db_name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_molecule_related_structures.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_list.family_prd_id in category PDBX_REFERENCE_MOLECULE_FAMILY_LIST.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * The formula for the reference entity. Formulae are written
   *  according to the rules:
   * 
   *  1. Only recognised element symbols may be used.
   * 
   *  2. Each element symbol is followed by a 'count' number. A count
   *     of '1' may be omitted.
   * 
   *  3. A space or parenthesis must separate each element symbol and
   *     its count, but in general parentheses are not used.
   * 
   *  4. The order of elements depends on whether or not carbon is
   *     present. If carbon is present, the order should be: C, then
   *     H, then the other elements in alphabetical order of their
   *     symbol. If carbon is not present, the elements are listed
   *     purely in alphabetic order of their symbol. This is the
   *     'Hill' system used by Chemical Abstracts.
   * 
   * Examples:
   * C18 H19 N7 O8 S
   */
  readonly formula?: Maybe<Scalars['String']>;
  /**
   * The chemical name for the structure entry in the related database
   * 
   * Examples:
   * actinomycn
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_molecule_related_structures.ordinal distinguishes
   *  related structural data for each entity.
   */
  readonly ordinal: Scalars['Int'];
};

export type RcsbEntitySourceOrganism = {
  /**
   * The beginning polymer sequence position for the polymer section corresponding
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly beg_seq_num?: Maybe<Scalars['Int']>;
  /** The common name for the source organism assigned by the PDB depositor. */
  readonly common_name?: Maybe<Scalars['String']>;
  /**
   * The ending polymer sequence position for the polymer section corresponding
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly end_seq_num?: Maybe<Scalars['Int']>;
  /**
   * Common names associated with this taxonomy code aggregated by the NCBI Taxonomy Database.
   * 
   *   These name correspond to the taxonomy identifier assigned by the PDB depositor.
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_common_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The parent scientific name in the NCBI taxonomy hierarchy (depth=1) of the source organism assigned by the PDB depositor.
   * 
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_parent_scientific_name?: Maybe<Scalars['String']>;
  /**
   * The scientific name associated with this taxonomy code aggregated by the NCBI Taxonomy Database.
   * 
   *   This name corresponds to the taxonomy identifier assigned by the PDB depositor.
   * 
   * 
   * References:
   * 
   * Sayers EW, Barrett T, Benson DA, Bryant SH, Canese K, Chetvernin V,
   * Church DM, DiCuccio M, Edgar R, Federhen S, Feolo M, Geer LY,
   * Helmberg W, Kapustin Y, Landsman D, Lipman DJ, Madden TL, Maglott DR,
   * Miller V, Mizrachi I, Ostell J, Pruitt KD, Schuler GD, Sequeira E,
   * Sherry ST, Shumway M, Sirotkin K, Souvorov A, Starchenko G,
   * Tatusova TA, Wagner L, Yaschenko E, Ye J (2009). Database resources
   * of the National Center for Biotechnology Information. Nucleic Acids
   * Res. 2009 Jan;37(Database issue):D5-15. Epub 2008 Oct 21.
   * 
   * Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW (2009).
   * GenBank. Nucleic Acids Res. 2009 Jan;37(Database issue):D26-31.
   * Epub 2008 Oct 21.
   */
  readonly ncbi_scientific_name?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier for the gene source organism assigned by the PDB depositor.
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly ncbi_taxonomy_id?: Maybe<Scalars['Int']>;
  /** An identifier for the entity segment. */
  readonly pdbx_src_id: Scalars['String'];
  /**
   * A code indicating the provenance of the source organism details for the entity
   * 
   * Allowable values:
   * PDB Primary Data
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  readonly rcsb_gene_name?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganismRcsbGeneName>>>;
  /** The scientific name of the source organism assigned by the PDB depositor. */
  readonly scientific_name?: Maybe<Scalars['String']>;
  /**
   * The source type for the entity
   * 
   * Allowable values:
   * genetically engineered, natural, synthetic
   */
  readonly source_type?: Maybe<Scalars['String']>;
  readonly taxonomy_lineage?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganismTaxonomyLineage>>>;
};

export type RcsbUniprotContainerIdentifiers = {
  readonly reference_sequence_identifiers?: Maybe<ReadonlyArray<Maybe<RcsbUniprotContainerIdentifiersReferenceSequenceIdentifiers>>>;
  /** Primary accession number of a given UniProtKB entry. */
  readonly uniprot_id?: Maybe<Scalars['String']>;
};

export type DrugbankInfo = {
  /** The DrugBank drug affected organisms. */
  readonly affected_organisms?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The Anatomical Therapeutic Chemical Classification System (ATC) codes. */
  readonly atc_codes?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** DrugBank drug brand names. */
  readonly brand_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * The DrugBank assigned Chemical Abstracts Service identifier.
   * 
   * Examples:
   * 56-65-5
   */
  readonly cas_number?: Maybe<Scalars['String']>;
  /** The DrugBank drug description. */
  readonly description?: Maybe<Scalars['String']>;
  /** The DrugBank drug categories. */
  readonly drug_categories?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The DrugBank drug drug groups. */
  readonly drug_groups?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** The DrugBank accession code */
  readonly drugbank_id: Scalars['String'];
  /**
   * The DrugBank drug indication.
   * 
   * Examples:
   * For nutritional supplementation, also for treating dietary shortage or imbalance
   */
  readonly indication?: Maybe<Scalars['String']>;
  /**
   * The DrugBank drug mechanism of actions.
   * 
   * Examples:
   * ATP is able to store and transport chemical energy within cells.
   */
  readonly mechanism_of_action?: Maybe<Scalars['String']>;
  /** The DrugBank drug name. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The DrugBank drug pharmacology.
   * 
   * Examples:
   * Adenosine triphosphate (ATP) is the nucleotide known in biochemistry as the "molecular currency" of intracellular energy transfer; that is, ATP is able to store and transport chemical energy within cells. ATP also plays an important role in the synthesis of nucleic acids. The total quantity of ATP in the human body is about 0.1 mole. The energy used by human cells requires the hydrolysis of 200 to 300 moles of ATP daily. This means that each ATP molecule is recycled 2000 to 3000 times during a single day. ATP cannot be stored, hence its consumption must closely follow its synthesis.
   */
  readonly pharmacology?: Maybe<Scalars['String']>;
  /** DrugBank drug name synonyms. */
  readonly synonyms?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type RcsbNonpolymerEntityFeatureSummary = {
  /** Non-polymer(ligand) chemical component identifier for the entity. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** The feature count. */
  readonly count?: Maybe<Scalars['Int']>;
  /** The maximum feature length. */
  readonly maximum_length?: Maybe<Scalars['Int']>;
  /** The maximum feature value. */
  readonly maximum_value?: Maybe<Scalars['Float']>;
  /** The minimum feature length. */
  readonly minimum_length?: Maybe<Scalars['Int']>;
  /** The minimum feature value. */
  readonly minimum_value?: Maybe<Scalars['Float']>;
  /**
   * Type or category of the feature.
   * 
   * Allowable values:
   * SUBJECT_OF_INVESTIGATION
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbExternalReferences = {
  /** ID (accession) from external resource linked to this entry. */
  readonly id: Scalars['String'];
  /** Link to this entry in external resource. */
  readonly link: Scalars['String'];
  /**
   * Internal identifier for external resource.
   * 
   * Allowable values:
   * OLDERADO, BMRB, NDB, SB GRID, STORE SYNCHROTRON, PROTEIN DIFFRACTION, EM DATA RESOURCE
   */
  readonly type: Scalars['String'];
};

export type PdbxReferenceMolecule = {
  /**
   * For entities represented as single molecules, the identifier
   *  corresponding to the chemical definition for the molecule.
   * 
   * Examples:
   * 0Z3, CD9
   */
  readonly chem_comp_id?: Maybe<Scalars['String']>;
  /**
   * Broadly defines the function of the entity.
   * 
   * Allowable values:
   * Antagonist, Anthelmintic, Antibiotic, Anticancer, Anticoagulant, Antifungal, Antiinflammatory, Antimicrobial, Antineoplastic, Antiparasitic, Antiretroviral, Antithrombotic, Antitumor, Antiviral, CASPASE inhibitor, Chaperone binding, Enzyme inhibitor, Growth factor, Immunosuppressant, Inhibitor, Lantibiotic, Metabolism, Metal transport, Oxidation-reduction, RNA synthesis Inhibitor, Receptor, Synthetic opioid, Thrombin inhibitor, Toxin, Transition state mimetic, Transport activator, Trypsin inhibitor, Unknown
   */
  readonly class?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Evidence for the assignment of _pdbx_reference_molecule.class */
  readonly class_evidence_code?: Maybe<Scalars['String']>;
  /** Special details about this molecule. */
  readonly compound_details?: Maybe<Scalars['String']>;
  /** Description of this molecule. */
  readonly description?: Maybe<Scalars['String']>;
  /**
   * The formula for the reference entity. Formulae are written
   *  according to the rules:
   * 
   *  1. Only recognised element symbols may be used.
   * 
   *  2. Each element symbol is followed by a 'count' number. A count
   *     of '1' may be omitted.
   * 
   *  3. A space or parenthesis must separate each element symbol and
   *     its count, but in general parentheses are not used.
   * 
   *  4. The order of elements depends on whether or not carbon is
   *     present. If carbon is present, the order should be: C, then
   *     H, then the other elements in alphabetical order of their
   *     symbol. If carbon is not present, the elements are listed
   *     purely in alphabetic order of their symbol. This is the
   *     'Hill' system used by Chemical Abstracts.
   * 
   * Examples:
   * C18 H19 N7 O8 S
   */
  readonly formula?: Maybe<Scalars['String']>;
  /** Formula mass in daltons of the entity. */
  readonly formula_weight?: Maybe<Scalars['Float']>;
  /**
   * A name of the entity.
   * 
   * Examples:
   * thiostrepton
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_reference_molecule.prd_id is the unique identifier 
   *  for the reference molecule in this family.  
   * 
   *  By convention this ID uniquely identifies the reference molecule in 
   *  in the PDB reference dictionary.  
   * 
   *  The ID has the template form PRD_dddddd (e.g. PRD_000001)
   */
  readonly prd_id: Scalars['String'];
  /**
   * Defines the current PDB release status for this molecule definition.
   * 
   * Allowable values:
   * HOLD, OBS, REL, WAIT
   */
  readonly release_status?: Maybe<Scalars['String']>;
  /** Assigns the identifier of the reference molecule that has replaced this molecule. */
  readonly replaced_by?: Maybe<Scalars['String']>;
  /**
   * Assigns the identifier for the reference molecule which have been replaced 
   *  by this reference molecule.
   *  Multiple molecule identifier codes should be separated by commas.
   */
  readonly replaces?: Maybe<Scalars['String']>;
  /**
   * Defines how this entity is represented in PDB data files.
   * 
   * Allowable values:
   * branched, polymer, single molecule
   */
  readonly represent_as?: Maybe<Scalars['String']>;
  /** The PDB accession code for the entry containing a representative example of this molecule. */
  readonly representative_PDB_id_code?: Maybe<Scalars['String']>;
  /**
   * Defines the structural classification of the entity.
   * 
   * Allowable values:
   * Amino acid, Aminoglycoside, Ansamycin, Anthracycline, Anthraquinone, Chalkophore, Chalkophore, Polypeptide, Chromophore, Cyclic depsipeptide, Cyclic lipopeptide, Cyclic peptide, Glycopeptide, Heterocyclic, Imino sugar, Keto acid, Lipoglycopeptide, Lipopeptide, Macrolide, Non-polymer, Nucleoside, Oligopeptide, Oligosaccharide, Peptaibol, Peptide-like, Polycyclic, Polypeptide, Polysaccharide, Quinolone, Siderophore, Thiolactone, Thiopeptide, Tricyclic pentaglycosidic antineoplastic antibiotic, Unknown
   */
  readonly type?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Evidence for the assignment of _pdbx_reference_molecule.type */
  readonly type_evidence_code?: Maybe<Scalars['String']>;
};

export type PdbxReferenceMoleculeFamily = {
  /**
   * The value of _pdbx_reference_entity.family_prd_id must uniquely identify a record in the
   *  PDBX_REFERENCE_MOLECULE_FAMILY list.
   * 
   *  By convention this ID uniquely identifies the reference family in 
   *  in the PDB reference dictionary.  
   * 
   *  The ID has the template form FAM_dddddd (e.g. FAM_000001)
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * The entity family name.
   * 
   * Examples:
   * actinomycin, adriamycin
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Assigns the current PDB release status for this family.
   * 
   * Allowable values:
   * HOLD, OBS, REL, WAIT
   */
  readonly release_status?: Maybe<Scalars['String']>;
  /** Assigns the identifier of the family that has replaced this component. */
  readonly replaced_by?: Maybe<Scalars['String']>;
  /**
   * Assigns the identifier for the family which have been replaced by this family.
   *  Multiple family identifier codes should be separated by commas.
   */
  readonly replaces?: Maybe<Scalars['String']>;
};

export type PdbxEntitySrcSyn = {
  /**
   * A description of special aspects of the source for the
   *  synthetic entity.
   * 
   * Examples:
   * This sequence occurs naturally in humans.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier of the organism from which the sequence of
   *  the synthetic entity was derived.
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly ncbi_taxonomy_id?: Maybe<Scalars['String']>;
  /**
   * The common name of the organism from which the sequence of
   *  the synthetic entity was derived.
   * 
   * Examples:
   * house mouse
   */
  readonly organism_common_name?: Maybe<Scalars['String']>;
  /**
   * The scientific name of the organism from which the sequence of
   *  the synthetic entity was derived.
   * 
   * Examples:
   * synthetic construct, Mus musculus
   */
  readonly organism_scientific?: Maybe<Scalars['String']>;
  /**
   * This data item identifies cases in which an alternative source
   *  modeled.
   * 
   * Allowable values:
   * model, sample
   */
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>;
  /**
   * The beginning polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>;
  /**
   * The ending polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>;
  /** This data item is an ordinal identifier for pdbx_entity_src_syn data records. */
  readonly pdbx_src_id: Scalars['Int'];
};

export type RcsbPolymerEntityFeature = {
  /** Identifies the version of the feature assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** A description for the feature. */
  readonly description?: Maybe<Scalars['String']>;
  /** An identifier for the feature. */
  readonly feature_id?: Maybe<Scalars['String']>;
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeatureFeaturePositions>>>;
  /** A name for the feature. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * Code identifying the individual, organization or program that
   *  assigned the feature.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * Code residue coordinate system for the assigned feature.
   * 
   * Allowable values:
   * NCBI, PDB entity, UniProt
   */
  readonly reference_scheme?: Maybe<Scalars['String']>;
  /**
   * A type or category of the feature.
   * 
   * Allowable values:
   * artifact, modified_monomer, mutation
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type Cell = {
  /**
   * The number of the polymeric chains in a unit cell. In the case
   *  of heteropolymers, Z is the number of occurrences of the most
   *  populous chain.
   * 
   *  This data item is provided for compatibility with the original
   *  Protein Data Bank format, and only for that purpose.
   */
  readonly Z_PDB?: Maybe<Scalars['Int']>;
  /** Unit-cell angle alpha of the reported structure in degrees. */
  readonly angle_alpha?: Maybe<Scalars['Float']>;
  /** Unit-cell angle beta of the reported structure in degrees. */
  readonly angle_beta?: Maybe<Scalars['Float']>;
  /** Unit-cell angle gamma of the reported structure in degrees. */
  readonly angle_gamma?: Maybe<Scalars['Float']>;
  /**
   * The number of the formula units in the unit cell as specified
   *  by _chemical_formula.structural, _chemical_formula.moiety or
   *  _chemical_formula.sum.
   */
  readonly formula_units_Z?: Maybe<Scalars['Int']>;
  /**
   * Unit-cell length a corresponding to the structure reported in
   * angstroms.
   */
  readonly length_a?: Maybe<Scalars['Float']>;
  /**
   * Unit-cell length b corresponding to the structure reported in
   *  angstroms.
   */
  readonly length_b?: Maybe<Scalars['Float']>;
  /**
   * Unit-cell length c corresponding to the structure reported in
   * angstroms.
   */
  readonly length_c?: Maybe<Scalars['Float']>;
  /**
   * To further identify unique axis if necessary.  E.g., P 21 with
   *  an unique C axis will have 'C' in this field.
   */
  readonly pdbx_unique_axis?: Maybe<Scalars['String']>;
  /**
   * Cell volume V in angstroms cubed.
   * 
   *  V = a b c (1 - cos^2^~alpha~ - cos^2^~beta~ - cos^2^~gamma~
   *             + 2 cos~alpha~ cos~beta~ cos~gamma~)^1/2^
   * 
   *  a     = _cell.length_a
   *  b     = _cell.length_b
   *  c     = _cell.length_c
   *  alpha = _cell.angle_alpha
   *  beta  = _cell.angle_beta
   *  gamma = _cell.angle_gamma
   */
  readonly volume?: Maybe<Scalars['Float']>;
};

export type EmDiffractionShell = {
  /** Pointer to EM CRYSTALLOGRAPHY STATS */
  readonly em_diffraction_stats_id?: Maybe<Scalars['String']>;
  /**
   * Completeness of the structure factor data within this resolution shell, in percent
   * 
   * Examples:
   * 93.2
   */
  readonly fourier_space_coverage?: Maybe<Scalars['Float']>;
  /**
   * High resolution limit for this shell (Angstroms)
   * 
   * Examples:
   * 3.0
   */
  readonly high_resolution?: Maybe<Scalars['Float']>;
  /** Unique identifier for the category em_diffraction_shell */
  readonly id: Scalars['String'];
  /**
   * Low resolution limit for this shell (Angstroms)
   * 
   * Examples:
   * 5.5
   */
  readonly low_resolution?: Maybe<Scalars['Float']>;
  /**
   * Multiplicity (average number of measurements) for the structure factors in this resolution shell
   * 
   * Examples:
   * 2.5
   */
  readonly multiplicity?: Maybe<Scalars['Float']>;
  /**
   * Number of measured structure factors in this resolution shell
   * 
   * Examples:
   * 244
   */
  readonly num_structure_factors?: Maybe<Scalars['Int']>;
  /**
   * Phase residual for this resolution shell, in degrees
   * 
   * Examples:
   * 13.5
   */
  readonly phase_residual?: Maybe<Scalars['Float']>;
};

export type PdbxReferenceMoleculeList = {
  /**
   * The value of _pdbx_reference_molecule_list.family_prd_id is a reference to 
   *  _pdbx_reference_molecule_family.family_prd_id' in category PDBX_REFERENCE_MOLECULE_FAMILY.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * The value of _pdbx_reference_molecule_list.prd_id is the unique identifier 
   *  for the reference molecule in this family.  
   * 
   *  By convention this ID uniquely identifies the reference molecule in 
   *  in the PDB reference dictionary.  
   * 
   *  The ID has the template form PRD_dddddd (e.g. PRD_000001)
   */
  readonly prd_id: Scalars['String'];
};

export type ExptlCrystal = {
  /**
   * The colour of the crystal.
   * 
   * Examples:
   * dark green
   */
  readonly colour?: Maybe<Scalars['String']>;
  /**
   * The density of the crystal, expressed as the ratio of the
   *  volume of the asymmetric unit to the molecular mass of a
   *  monomer of the structure, in units of angstroms^3^ per dalton.
   * 
   *  Ref: Matthews, B. W. (1968). J. Mol. Biol. 33, 491-497.
   */
  readonly density_Matthews?: Maybe<Scalars['Float']>;
  /**
   * Density values measured using standard chemical and physical
   *  methods. The units are megagrams per cubic metre (grams per
   *  cubic centimetre).
   */
  readonly density_meas?: Maybe<Scalars['Float']>;
  /**
   * Density value P calculated from the crystal cell and contents,
   *  expressed as per cent solvent.
   * 
   *  P = 1 - (1.23 N MMass) / V
   * 
   *  N     = the number of molecules in the unit cell
   *  MMass = the molecular mass of each molecule (gm/mole)
   *  V     = the volume of the unit cell (A^3^)
   *  1.23  = a conversion factor evaluated as:
   * 
   *          (0.74 cm^3^/g) (10^24^ A^3^/cm^3^)
   *          --------------------------------------
   *               (6.02*10^23^) molecules/mole
   * 
   *          where 0.74 is an assumed value for the partial specific
   *          volume of the molecule
   */
  readonly density_percent_sol?: Maybe<Scalars['Float']>;
  /**
   * A description of the quality and habit of the crystal.
   *  The crystal dimensions should not normally be reported here;
   *  use instead the specific items in the EXPTL_CRYSTAL category
   *  relating to size for the gross dimensions of the crystal and
   *  data items in the EXPTL_CRYSTAL_FACE category to describe the
   *  relationship between individual faces.
   */
  readonly description?: Maybe<Scalars['String']>;
  /**
   * The value of _exptl_crystal.id must uniquely identify a record in
   *  the EXPTL_CRYSTAL list.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   */
  readonly id: Scalars['String'];
  /**
   * The of the distribution of mis-orientation angles specified in degrees
   * of all the unit cells in the crystal. Lower mosaicity indicates better 
   * ordered crystals.
   */
  readonly pdbx_mosaicity?: Maybe<Scalars['Float']>;
  /** The uncertainty in the mosaicity estimate for the crystal. */
  readonly pdbx_mosaicity_esd?: Maybe<Scalars['Float']>;
  /**
   * Details of crystal growth and preparation of the crystal (e.g.
   *  mounting) prior to the intensity measurements.
   * 
   * Examples:
   * mounted in an argon-filled quartz capillary
   */
  readonly preparation?: Maybe<Scalars['String']>;
};

export type EmStaining = {
  /**
   * Staining procedure used in the specimen preparation.
   * 
   * Examples:
   * Negatively stained EM specimens were prepared using a carbon-sandwich technique
   *   and uranyl-formate stain.
   */
  readonly details?: Maybe<Scalars['String']>;
  /** This data item is the primary key of the category. */
  readonly id: Scalars['String'];
  /**
   * The staining  material.
   * 
   * Examples:
   * Uranyl Acetate
   */
  readonly material?: Maybe<Scalars['String']>;
  /** Foreign key relationship to the EMD SPECIMEN category */
  readonly specimen_id?: Maybe<Scalars['String']>;
  /**
   * type of staining
   * 
   * Allowable values:
   * NEGATIVE, NONE, POSITIVE
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type PdbxNmrSampleDetails = {
  /**
   * A complete description of each NMR sample. Include the concentration 
   * and concentration units for each component (include buffers, etc.). For each
   * component describe the isotopic composition, including the % labeling level,
   * if known. 
   * 
   * For example:
   * 1. Uniform (random) labeling with 15N: U-15N
   * 2. Uniform (random) labeling with 13C, 15N at known labeling
   *    levels: U-95% 13C;U-98% 15N
   * 3. Residue selective labeling: U-95% 15N-Thymine
   * 4. Site specific labeling: 95% 13C-Ala18,
   * 5. Natural abundance labeling in an otherwise uniformly labeled
   *    biomolecule is designated by NA: U-13C; NA-K,H
   * 
   * Examples:
   * 2mM Ribonuclease  U-15N,13C; 50mM phosphate buffer NA; 90% H2O, 10% D2O
   */
  readonly contents?: Maybe<Scalars['String']>;
  /**
   * Brief description of the sample providing additional information not captured by other items in the category.
   * 
   * Examples:
   * The added glycerol was used to raise the viscosity of the solution to 1.05 poisson.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * A value that uniquely identifies this sample from the other samples listed 
   * in the entry.
   * 
   * Examples:
   * 15N_sample
   */
  readonly label?: Maybe<Scalars['String']>;
  /**
   * The name (number) of the sample.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly solution_id: Scalars['String'];
  /**
   * The solvent system used for this sample.
   * 
   * Examples:
   * 90% H2O, 10% D2O
   */
  readonly solvent_system?: Maybe<Scalars['String']>;
  /**
   * A descriptive term for the sample that defines the general physical properties 
   * of the sample.
   * 
   * Allowable values:
   * bicelle, emulsion, fiber, fibrous protein, filamentous virus, gel solid, gel solution, liposome, lyophilized powder, membrane, micelle, oriented membrane film, polycrystalline powder, reverse micelle, single crystal, solid, solution
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbPolymerEntityKeywords = {
  /** Keywords describing this polymer entity. */
  readonly text?: Maybe<Scalars['String']>;
};

export type DrugbankContainerIdentifiers = {
  /** The DrugBank accession code */
  readonly drugbank_id: Scalars['String'];
};

export type PdbxDatabaseRelated = {
  /**
   * The identifying content type of the related entry.
   * 
   * Allowable values:
   * associated EM volume, associated NMR restraints, associated SAS data, associated structure factors, complete structure, derivative structure, ensemble, minimized average structure, native structure, other, other EM volume, protein target sequence and/or protocol data, re-refinement, representative structure, split, unspecified
   */
  readonly content_type: Scalars['String'];
  /**
   * The identifying code in the related database.
   * 
   * Examples:
   * 1ABC, BDL001
   */
  readonly db_id: Scalars['String'];
  /**
   * The name of the database containing the related entry.
   * 
   * Examples:
   * PDB  - Protein Databank
   * NDB  - Nucleic Acid Database
   * BMRB - BioMagResBank
   * EMDB - Electron Microscopy Database
   * BMCD - Biological Macromolecule Crystallization Database
   * TargetTrack - Target Registration and Protocol Database
   * SASBDB - Small Angle Scattering Biological Data Bank
   */
  readonly db_name: Scalars['String'];
  /**
   * A description of the related entry.
   * 
   * Examples:
   * 1ABC contains the same protein complexed with Netropsin.
   */
  readonly details?: Maybe<Scalars['String']>;
};

export type RcsbUniprotProteinGene = {
  readonly name?: Maybe<ReadonlyArray<Maybe<GeneName>>>;
};

export type PdbxStructAssemblyAuthEvidence = {
  /** This item references an assembly in pdbx_struct_assembly */
  readonly assembly_id: Scalars['String'];
  /** Provides any additional information regarding the evidence of this assembly */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * Provides the experimental method to determine the state of this assembly
   * 
   * Allowable values:
   * SAXS, assay for oligomerization, cross-linking, equilibrium centrifugation, fluorescence resonance energy transfer, gel filtration, homology, immunoprecipitation, isothermal titration calorimetry, light scattering, mass spectrometry, microscopy, native gel electrophoresis, none, scanning transmission electron microscopy, surface plasmon resonance
   */
  readonly experimental_support?: Maybe<Scalars['String']>;
  /** Identifies a unique record in pdbx_struct_assembly_auth_evidence. */
  readonly id: Scalars['String'];
};

export type RcsbNonpolymerInstanceAnnotation = {
  /** An identifier for the annotation. */
  readonly annotation_id?: Maybe<Scalars['String']>;
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceAnnotationAnnotationLineage>>>;
  /** Identifies the version of the annotation assignment. */
  readonly assignment_version?: Maybe<Scalars['String']>;
  /** Chemical component identifier. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** A description for the annotation. */
  readonly description?: Maybe<Scalars['String']>;
  /** A name for the annotation. */
  readonly name?: Maybe<Scalars['String']>;
  /** Ordinal identifier for this category */
  readonly ordinal: Scalars['Int'];
  /**
   * Code identifying the individual, organization or program that
   *  assigned the annotation.
   */
  readonly provenance_source?: Maybe<Scalars['String']>;
  /**
   * A type or category of the annotation.
   * 
   * Allowable values:
   * HAS_COVALENT_LINKAGE, HAS_METAL_COORDINATION_LINKAGE
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntityNameCom = {
  /** A common name for the nonpolymer entity. */
  readonly name: Scalars['String'];
};

export type Software = {
  /**
   * The classification of the program according to its
   *  major function.
   * 
   * Examples:
   * data collection, data reduction, phasing, model building, refinement, validation, other
   */
  readonly classification?: Maybe<Scalars['String']>;
  /**
   * The recognized contact author of the software. This could be
   *  the original author, someone who has modified the code or
   *  someone who maintains the code.  It should be the person
   *  most commonly associated with the code.
   * 
   * Examples:
   * T. Alwyn Jones, Axel Brunger
   */
  readonly contact_author?: Maybe<Scalars['String']>;
  /**
   * The e-mail address of the person specified in
   *  _software.contact_author.
   * 
   * Examples:
   * bourne@sdsc.edu
   */
  readonly contact_author_email?: Maybe<Scalars['String']>;
  /**
   * The date the software was released.
   * 
   * Examples:
   * 1991-10-01, 1990-04-30
   */
  readonly date?: Maybe<Scalars['String']>;
  /**
   * Description of the software.
   * 
   * Examples:
   * Uses method of restrained least squares
   */
  readonly description?: Maybe<Scalars['String']>;
  /**
   * The major computing language in which the software is
   *  coded.
   * 
   * Allowable values:
   * Ada, Awk, Basic, C, C++, C/C++, Fortran, Fortran 77, Fortran 90, Fortran_77, Java, Java & Fortran, Other, Pascal, Perl, Python, Python/C++, Tcl, assembler, csh, ksh, sh
   */
  readonly language?: Maybe<Scalars['String']>;
  /**
   * The URL for an Internet address at which
   *  details of the software can be found.
   * 
   * Examples:
   * http://rosebud.sdsc.edu/projects/pb/IUCr/software.html, ftp://ftp.sdsc.edu/pub/sdsc/biology/
   */
  readonly location?: Maybe<Scalars['String']>;
  /**
   * The name of the software.
   * 
   * Examples:
   * Merlot, O, Xengen, X-plor
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The name of the operating system under which the software
   *  runs.
   * 
   * Examples:
   * Ultrix, OpenVMS, DOS, Windows 95, Windows NT, Irix, HPUX, DEC Unix
   */
  readonly os?: Maybe<Scalars['String']>;
  /**
   * An ordinal index for this category
   * 
   * Examples:
   * 1, 2
   */
  readonly pdbx_ordinal: Scalars['Int'];
  /**
   * The classification of the software according to the most
   *  common types.
   * 
   * Allowable values:
   * filter, jiffy, library, other, package, program
   */
  readonly type?: Maybe<Scalars['String']>;
  /**
   * The version of the software.
   * 
   * Examples:
   * v1.0, beta, 3.1-2, unknown
   */
  readonly version?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerStructConn = {
  readonly connect_partner?: Maybe<RcsbNonpolymerStructConnConnectPartner>;
  readonly connect_target?: Maybe<RcsbNonpolymerStructConnConnectTarget>;
  /**
   * The connection type.
   * 
   * Allowable values:
   * covalent bond, disulfide bridge, hydrogen bond, ionic interaction, metal coordination, mismatched base pairs
   */
  readonly connect_type?: Maybe<Scalars['String']>;
  /** A description of special details of the connection. */
  readonly description?: Maybe<Scalars['String']>;
  /** Distance value for this contact. */
  readonly dist_value?: Maybe<Scalars['Float']>;
  /**
   * The value of _rcsb_nonpolymer_struct_conn.id is an identifier for connection.
   * 
   *  Note that this item need not be a number; it can be any unique
   *  identifier.
   */
  readonly id?: Maybe<Scalars['String']>;
  /**
   * The value of _rcsb_nonpolymer_struct_conn.id must uniquely identify a record in
   *  the rcsb_nonpolymer_struct_conn list.
   */
  readonly ordinal_id: Scalars['Int'];
  /**
   * The chemical or structural role of the interaction
   * 
   * Allowable values:
   * C-Mannosylation, N-Glycosylation, O-Glycosylation
   */
  readonly role?: Maybe<Scalars['String']>;
  /**
   * The chemical bond order associated with the specified atoms in
   *  this contact.
   * 
   * Allowable values:
   * doub, quad, sing, trip
   */
  readonly value_order?: Maybe<Scalars['String']>;
};

export type PdbxPrdAudit = {
  /**
   * The action associated with this audit record.
   * 
   * Allowable values:
   * Create molecule, Initial release, Modify audit, Modify class, Modify linkage, Modify molecule name, Modify representation, Modify sequence, Modify taxonomy organism, Modify type, Obsolete molecule, Other modification
   */
  readonly action_type: Scalars['String'];
  /**
   * The initials of the annotator creating of modifying the molecule.
   * 
   * Examples:
   * JO, SJ, KB
   */
  readonly annotator?: Maybe<Scalars['String']>;
  /** The date associated with this audit record. */
  readonly date: Scalars['Date'];
  /**
   * Additional details decribing this change.
   * 
   * Examples:
   * Revise molecule sequence.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _pdbx_reference_molecule.prd_id in the 
   * 	       pdbx_reference_molecule category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * An identifier for the wwPDB site creating or modifying the molecule.
   * 
   * Allowable values:
   * BMRB, PDBC, PDBJ, PDBe, RCSB
   */
  readonly processing_site?: Maybe<Scalars['String']>;
};

export type CoreUniprot = {
  /** Primary accession number of a given UniProtKB entry. */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  readonly rcsb_uniprot_accession?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** UniProt pairwise sequence alignments. */
  readonly rcsb_uniprot_alignments?: Maybe<RcsbUniprotAlignments>;
  readonly rcsb_uniprot_annotation?: Maybe<ReadonlyArray<Maybe<RcsbUniprotAnnotation>>>;
  readonly rcsb_uniprot_container_identifiers: RcsbUniprotContainerIdentifiers;
  readonly rcsb_uniprot_entry_name?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  readonly rcsb_uniprot_external_reference?: Maybe<ReadonlyArray<Maybe<RcsbUniprotExternalReference>>>;
  readonly rcsb_uniprot_feature?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeature>>>;
  readonly rcsb_uniprot_keyword?: Maybe<ReadonlyArray<Maybe<RcsbUniprotKeyword>>>;
  readonly rcsb_uniprot_protein?: Maybe<RcsbUniprotProtein>;
};

export type RcsbNonpolymerEntityAnnotationAnnotationLineage = {
  /** Members of the annotation lineage as parent lineage depth (1-N) */
  readonly depth?: Maybe<Scalars['Int']>;
  /** Members of the annotation lineage as parent class identifiers. */
  readonly id?: Maybe<Scalars['String']>;
  /** Members of the annotation lineage as parent class names. */
  readonly name?: Maybe<Scalars['String']>;
};

export type PdbxNmrRepresentative = {
  /**
   * If a member of the ensemble has been selected as a representative
   *  structure, identify it by its model number.
   * 
   * Examples:
   * 15
   */
  readonly conformer_id?: Maybe<Scalars['String']>;
  /**
   * By highlighting the appropriate choice(s), describe the criteria used to 
   * select this structure as a representative structure, or if an average 
   * structure has been calculated describe how this was done.
   * 
   * Examples:
   * The structure closest to the average.
   * The structure with the lowest energy was selected.
   * The structure with the fewest number of violations was selected.
   * A minimized average structure was calculated.
   */
  readonly selection_criteria?: Maybe<Scalars['String']>;
};

export type PdbxSerialCrystallographyDataReduction = {
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of frames collected
   *  in which the crystal was hit.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly crystal_hits?: Maybe<Scalars['Int']>;
  /**
   * The data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   * 
   * Examples:
   * 1
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of frames collected
   *  in which a droplet was hit.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly droplet_hits?: Maybe<Scalars['Int']>;
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of data frames collected
   *  in which the sample was hit.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly frame_hits?: Maybe<Scalars['Int']>;
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of data frames collected
   *  that contained a "hit" but failed to index.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly frames_failed_index?: Maybe<Scalars['Int']>;
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of data frames collected
   *  that were indexed.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly frames_indexed?: Maybe<Scalars['Int']>;
  /**
   * The total number of data frames collected for this
   *  data set.
   * 
   * Examples:
   * 20, 100
   */
  readonly frames_total?: Maybe<Scalars['Int']>;
  /**
   * For experiments in which samples are provided in a
   *  continuous stream, the total number of lattices indexed.
   * 
   * Examples:
   * 1200, 5750
   */
  readonly lattices_indexed?: Maybe<Scalars['Int']>;
  /** For FEL experiments, the number of pulse events in the dataset. */
  readonly xfel_pulse_events?: Maybe<Scalars['Int']>;
  /**
   * For FEL experiments, in which data collection was performed
   * 	       in batches, indicates which subset of the data collected
   *                were used in producing this dataset.
   */
  readonly xfel_run_numbers?: Maybe<Scalars['String']>;
};

export type RcsbClusterFlexibility = {
  /** Average RMSD refer to average pairwise RMSD (Root Mean Square Deviation of C-alpha atoms) between structures in the cluster (95% sequence identity) where a given entity belongs. */
  readonly avg_rmsd?: Maybe<Scalars['Float']>;
  /** Structural flexibility in the cluster (95% sequence identity) where a given entity belongs. */
  readonly label?: Maybe<Scalars['String']>;
  /** Link to the associated PDBFlex database entry. */
  readonly link?: Maybe<Scalars['String']>;
  /** Maximal RMSD refer to maximal pairwise RMSD (Root Mean Square Deviation of C-alpha atoms) between structures in the cluster (95% sequence identity) where a given entity belongs. */
  readonly max_rmsd?: Maybe<Scalars['Float']>;
  /**
   * Allowable values:
   * PDBFlex
   */
  readonly provenance_code?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntityContainerIdentifiers = {
  /** Instance identifiers corresponding to copies of the entity in this container. */
  readonly asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Author instance identifiers corresponding to copies of the entity in this container. */
  readonly auth_asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Entity identifier for the container. */
  readonly entity_id: Scalars['String'];
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /** Non-polymer(ligand) chemical component identifier for the entity in this container. */
  readonly nonpolymer_comp_id?: Maybe<Scalars['String']>;
  /** The BIRD identifier for the entity in this container. */
  readonly prd_id?: Maybe<Scalars['String']>;
  /**
   * A unique identifier for each object in this entity container formed by
   *  an underscore separated concatenation of entry and entity identifiers.
   */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  /**
   * Source of the reference database assignment
   * 
   * Allowable values:
   * PDB, RCSB
   */
  readonly reference_chemical_identifiers_provenance_source?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Reference resource accession code */
  readonly reference_chemical_identifiers_resource_accession?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * Reference resource name
   * 
   * Allowable values:
   * ChEBI, ChEMBL, DrugBank, PubChem
   */
  readonly reference_chemical_identifiers_resource_name?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type RcsbPolymerEntityNameSys = {
  /** The systematic name for the polymer entity. */
  readonly name: Scalars['String'];
  /**
   * The system used to generate the systematic name of the polymer entity.
   * 
   * Examples:
   * Chemical Abstracts conventions
   */
  readonly system?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntityKeywords = {
  /** Keywords describing this non-polymer entity. */
  readonly text?: Maybe<Scalars['String']>;
};

export type PdbxAuditRevisionDetails = {
  /**
   * The type of file that the pdbx_audit_revision_history record refers to.
   * 
   * Allowable values:
   * Chemical component, NMR restraints, NMR shifts, Structure factors, Structure model
   */
  readonly data_content_type: Scalars['String'];
  /** Additional details describing the revision. */
  readonly description?: Maybe<Scalars['String']>;
  /** Further details describing the revision. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * A unique identifier for the pdbx_audit_revision_details record.
   * 
   * Examples:
   * 1
   */
  readonly ordinal: Scalars['Int'];
  /**
   * The provider of the revision.
   * 
   * Allowable values:
   * author, repository
   */
  readonly provider?: Maybe<Scalars['String']>;
  /**
   * A pointer to  _pdbx_audit_revision_history.ordinal
   * 
   * Examples:
   * 1
   */
  readonly revision_ordinal: Scalars['Int'];
  /**
   * A type classification of the revision
   * 
   * Allowable values:
   * Coordinate replacement, Initial release, Obsolete
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type DrugbankTarget = {
  /** The type of target interaction. */
  readonly interaction_type?: Maybe<Scalars['String']>;
  /** The target name. */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _drugbank_target.ordinal distinguishes
   *  related examples for each chemical component.
   */
  readonly ordinal: Scalars['Int'];
  /** The organism common name. */
  readonly organism_common_name?: Maybe<Scalars['String']>;
  /**
   * The reference identifier code for the target interaction reference.
   * 
   * Examples:
   * Q9HD40
   */
  readonly reference_database_accession_code?: Maybe<Scalars['String']>;
  /**
   * The reference database name for the target interaction.
   * 
   * Allowable values:
   * UniProt
   */
  readonly reference_database_name?: Maybe<Scalars['String']>;
  /**
   * Target sequence expressed as string of one-letter amino acid codes.
   * 
   * Examples:
   * MAKQRSG...
   */
  readonly seq_one_letter_code?: Maybe<Scalars['String']>;
  /** The actions of the target interaction. */
  readonly target_actions?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type Struct = {
  /**
   * The item indicates whether the entry is a CASP target, a CASD-NMR target,
   *  or similar target participating in methods development experiments.
   * 
   * Allowable values:
   * N, Y
   */
  readonly pdbx_CASP_flag?: Maybe<Scalars['String']>;
  /**
   * An automatically generated descriptor for an NDB structure or
   *  the unstructured content of the PDB COMPND record.
   * 
   * Examples:
   * 5'-D(*CP*GP*CP*(HYD)AP*AP*AP*TP*TP*TP*GP*CP*G)-3'
   */
  readonly pdbx_descriptor?: Maybe<Scalars['String']>;
  /**
   * Text description of the methodology which produced this
   *  model structure.
   * 
   * Examples:
   * This model was produced from a 10 nanosecond Amber/MD simulation
   * starting from PDB structure ID 1ABC.
   */
  readonly pdbx_model_details?: Maybe<Scalars['String']>;
  /**
   * A description of the type of structure model.
   * 
   * Examples:
   * MINIMIZED AVERAGE
   */
  readonly pdbx_model_type_details?: Maybe<Scalars['String']>;
  /**
   * A title for the data block. The author should attempt to convey
   *  the essence of the structure archived in the CIF in the title,
   *  and to distinguish this structural result from others.
   * 
   * Examples:
   * T4 lysozyme mutant - S32A, 5'-D(*(I)CP*CP*GP*G)-3, T4 lysozyme mutant - S32A, hen egg white lysozyme at -30 degrees C, quail egg white lysozyme at 2 atmospheres
   */
  readonly title?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerStructConnConnectPartner = {
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_alt_id in the
   *  ATOM_SITE category.
   */
  readonly label_alt_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_asym_id in the
   *  ATOM_SITE category.
   */
  readonly label_asym_id: Scalars['String'];
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _chem_comp_atom.atom_id in the
   *  CHEM_COMP_ATOM category.
   */
  readonly label_atom_id?: Maybe<Scalars['String']>;
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_comp_id in the
   *  ATOM_SITE category.
   */
  readonly label_comp_id: Scalars['String'];
  /**
   * A component of the identifier for the partner in the structure
   *  connection.
   * 
   *  This data item is a pointer to _atom_site.label_seq_id in the
   *  ATOM_SITE category.
   */
  readonly label_seq_id?: Maybe<Scalars['Int']>;
  /**
   * Describes the symmetry operation that should be applied to the
   *  atom set specified by _rcsb_nonpolymer_struct_conn.connect_partner_label* to generate the
   *  partner in the structure connection.
   * 
   * Examples:
   * 1_555, 7_645
   */
  readonly symmetry?: Maybe<Scalars['String']>;
};

export type PdbxSerialCrystallographySampleDeliveryFixedTarget = {
  /** The number of crystals per dropplet or pore in fixed target */
  readonly crystals_per_unit?: Maybe<Scalars['Int']>;
  /** For a fixed target sample, a description of sample preparation */
  readonly description?: Maybe<Scalars['String']>;
  /** Any details pertinent to the fixed sample target */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   * 
   * Examples:
   * 1
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * Device used to control movement of the fixed sample
   * 
   * Examples:
   * DMC-4080
   */
  readonly motion_control?: Maybe<Scalars['String']>;
  /**
   * Method to prevent dehydration of sample
   * 
   * Examples:
   * seal, humidifed gas, flash freezing
   */
  readonly sample_dehydration_prevention?: Maybe<Scalars['String']>;
  /**
   * For a fixed target sample, mechanism to hold sample in the beam
   * 
   * Examples:
   * mesh, loop, grid
   */
  readonly sample_holding?: Maybe<Scalars['String']>;
  /** The sample solution content and concentration */
  readonly sample_solvent?: Maybe<Scalars['String']>;
  /**
   * Size of pore in grid supporting sample. Diameter or length in micrometres,
   *  e.g. pore diameter
   */
  readonly sample_unit_size?: Maybe<Scalars['Float']>;
  /**
   * Type of base holding the support
   * 
   * Examples:
   * goniometer
   */
  readonly support_base?: Maybe<Scalars['String']>;
  /** Velocity of sample horizontally relative to a perpendicular beam in millimetres/second */
  readonly velocity_horizontal?: Maybe<Scalars['Float']>;
  /** Velocity of sample vertically relative to a perpendicular beam in millimetres/second */
  readonly velocity_vertical?: Maybe<Scalars['Float']>;
};

export type PdbxReferenceEntityPolyLink = {
  /**
   * The atom identifier/name in the first of the two components making
   *  the linkage.
   */
  readonly atom_id_1?: Maybe<Scalars['String']>;
  /**
   * The atom identifier/name in the second of the two components making
   *  the linkage.
   */
  readonly atom_id_2?: Maybe<Scalars['String']>;
  /**
   * The component identifier in the first of the two components making the 
   *  linkage.
   * 
   *  This data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id 
   *  in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
   */
  readonly comp_id_1?: Maybe<Scalars['String']>;
  /**
   * The component identifier in the second of the two components making the 
   *  linkage.
   * 
   *  This data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id 
   *  in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
   */
  readonly comp_id_2?: Maybe<Scalars['String']>;
  /** The entity component identifier entity containing the linkage. */
  readonly component_id: Scalars['Int'];
  /**
   * For a polymer entity, the sequence number in the first of 
   *  the two components making the linkage. 
   * 
   *  This data item is a pointer to _pdbx_reference_entity_poly_seq.num 
   *  in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
   */
  readonly entity_seq_num_1?: Maybe<Scalars['Int']>;
  /**
   * For a polymer entity, the sequence number in the second of 
   *  the two components making the linkage. 
   * 
   *  This data item is a pointer to _pdbx_reference_entity_poly_seq.num 
   *  in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
   */
  readonly entity_seq_num_2?: Maybe<Scalars['Int']>;
  /**
   * The value of _pdbx_reference_entity_poly_link.link_id uniquely identifies
   *  a linkage within a polymer entity.
   */
  readonly link_id: Scalars['Int'];
  /**
   * The value of _pdbx_reference_entity_poly_link.prd_id is a reference
   *  _pdbx_reference_entity_list.prd_id in the PDBX_REFERENCE_ENTITY_POLY category.
   */
  readonly prd_id: Scalars['String'];
  /**
   * The reference entity id of the polymer entity containing the linkage.
   * 
   *  This data item is a pointer to _pdbx_reference_entity_poly.ref_entity_id 
   *  in the PDBX_REFERENCE_ENTITY_POLY category.
   */
  readonly ref_entity_id: Scalars['String'];
  /**
   * The bond order target for the non-standard linkage.
   * 
   * Allowable values:
   * arom, delo, doub, pi, poly, quad, sing, trip
   */
  readonly value_order?: Maybe<Scalars['String']>;
};

export type RcsbStructSymmetryRotationAxes = {
  /** coordinate */
  readonly end: ReadonlyArray<Maybe<Scalars['Float']>>;
  /** The number of times (order of rotation) that a subunit can be repeated by a rotation operation, being transformed into a new state indistinguishable from its starting state. */
  readonly order?: Maybe<Scalars['Int']>;
  /** coordinate */
  readonly start: ReadonlyArray<Maybe<Scalars['Float']>>;
};

export type Em3dFitting = {
  /**
   * Any additional details regarding fitting of atomic coordinates into
   *  the 3DEM volume, including data and considerations from other
   *  methods used in computation of the model.
   * 
   * Examples:
   * Initial local fitting was done using Chimera and then NMFF was used for flexible fitting.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The value of _em_3d_fitting.id must uniquely identify
   *  a fitting procedure of atomic coordinates
   *  into 3dem reconstructed map volume.
   */
  readonly id: Scalars['String'];
  /**
   * The method used to fit atomic coordinates
   *  into the 3dem reconstructed map.
   */
  readonly method?: Maybe<Scalars['String']>;
  /**
   * The overall B (temperature factor) value for the 3d-em volume.
   * 
   * Examples:
   * 200
   */
  readonly overall_b_value?: Maybe<Scalars['Float']>;
  /**
   * The refinement protocol used.
   * 
   * Allowable values:
   * AB INITIO MODEL, BACKBONE TRACE, FLEXIBLE FIT, OTHER, RIGID BODY FIT
   */
  readonly ref_protocol?: Maybe<Scalars['String']>;
  /**
   * A flag to indicate whether fitting was carried out in real
   *  or reciprocal refinement space.
   * 
   * Allowable values:
   * REAL, RECIPROCAL
   */
  readonly ref_space?: Maybe<Scalars['String']>;
  /**
   * The measure used to assess quality of fit of the atomic coordinates in the
   *  3DEM map volume.
   * 
   * Examples:
   * Cross-correlation coefficient
   */
  readonly target_criteria?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerEntity = {
  /** A description of special aspects of the entity. */
  readonly details?: Maybe<Scalars['String']>;
  /** Formula mass (KDa) of the entity. */
  readonly formula_weight?: Maybe<Scalars['Float']>;
  /** A description of the nonpolymer entity. */
  readonly pdbx_description?: Maybe<Scalars['String']>;
  /**
   * The number of molecules of the nonpolymer entity in the entry.
   * 
   * Examples:
   * 1, 2, 3
   */
  readonly pdbx_number_of_molecules?: Maybe<Scalars['Int']>;
};

export type EntitySrcGen = {
  /**
   * A unique identifier for the expression system. This
   *  should be extracted from a local list of expression
   *  systems.
   */
  readonly expression_system_id?: Maybe<Scalars['String']>;
  /**
   * The common name of the natural organism from which the gene was
   *  obtained.
   * 
   * Examples:
   * man, yeast, bacteria
   */
  readonly gene_src_common_name?: Maybe<Scalars['String']>;
  /**
   * A description of special aspects of the natural organism from
   *  which the gene was obtained.
   */
  readonly gene_src_details?: Maybe<Scalars['String']>;
  /**
   * The genus of the natural organism from which the gene was
   *  obtained.
   * 
   * Examples:
   * Homo, Saccharomyces, Escherichia
   */
  readonly gene_src_genus?: Maybe<Scalars['String']>;
  /**
   * The species of the natural organism from which the gene was
   *  obtained.
   * 
   * Examples:
   * sapiens, cerevisiae, coli
   */
  readonly gene_src_species?: Maybe<Scalars['String']>;
  /**
   * The strain of the natural organism from which the gene was
   *  obtained, if relevant.
   * 
   * Examples:
   * DH5a, BMH 71-18
   */
  readonly gene_src_strain?: Maybe<Scalars['String']>;
  /**
   * The tissue of the natural organism from which the gene was
   *  obtained.
   * 
   * Examples:
   * heart, liver, eye lens
   */
  readonly gene_src_tissue?: Maybe<Scalars['String']>;
  /**
   * The subcellular fraction of the tissue of the natural organism
   *  from which the gene was obtained.
   * 
   * Examples:
   * mitochondria, nucleus, membrane
   */
  readonly gene_src_tissue_fraction?: Maybe<Scalars['String']>;
  /**
   * The common name of the organism that served as host for the
   *  production of the entity.  Where full details of the protein
   *  production are available it would be expected that this item
   *  be derived from _entity_src_gen_express.host_org_common_name
   *  or via _entity_src_gen_express.host_org_tax_id
   * 
   * Examples:
   * yeast, bacteria
   */
  readonly host_org_common_name?: Maybe<Scalars['String']>;
  /**
   * A description of special aspects of the organism that served as
   *  host for the production of the entity. Where full details of
   *  the protein production are available it would be expected that
   *  this item would derived from _entity_src_gen_express.host_org_details
   */
  readonly host_org_details?: Maybe<Scalars['String']>;
  /**
   * The genus of the organism that served as host for the production
   *  of the entity.
   * 
   * Examples:
   * Saccharomyces, Escherichia
   */
  readonly host_org_genus?: Maybe<Scalars['String']>;
  /**
   * The species of the organism that served as host for the
   *  production of the entity.
   * 
   * Examples:
   * cerevisiae, coli
   */
  readonly host_org_species?: Maybe<Scalars['String']>;
  /**
   * This data item identifies cases in which an alternative source
   *  modeled.
   * 
   * Allowable values:
   * model, sample
   */
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>;
  /**
   * The beginning polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>;
  /** Information on the source which is not given elsewhere. */
  readonly pdbx_description?: Maybe<Scalars['String']>;
  /**
   * The ending polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>;
  /**
   * American Type Culture Collection tissue culture number.
   * 
   * Examples:
   * 6051
   */
  readonly pdbx_gene_src_atcc?: Maybe<Scalars['String']>;
  /**
   * Cell type.
   * 
   * Examples:
   * ENDOTHELIAL
   */
  readonly pdbx_gene_src_cell?: Maybe<Scalars['String']>;
  /**
   * The specific line of cells.
   * 
   * Examples:
   * HELA CELLS
   */
  readonly pdbx_gene_src_cell_line?: Maybe<Scalars['String']>;
  /**
   * Identifies the location inside (or outside) the cell.
   * 
   * Examples:
   * CYTOPLASM, NUCLEUS
   */
  readonly pdbx_gene_src_cellular_location?: Maybe<Scalars['String']>;
  /**
   * A domain or fragment of the molecule.
   * 
   * Examples:
   * CYTOPLASM, NUCLEUS
   */
  readonly pdbx_gene_src_fragment?: Maybe<Scalars['String']>;
  /** Identifies the gene. */
  readonly pdbx_gene_src_gene?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier for the gene source organism.
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4 
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly pdbx_gene_src_ncbi_taxonomy_id?: Maybe<Scalars['String']>;
  /**
   * Organized group of tissues that carries on a specialized function.
   * 
   * Examples:
   * KIDNEY, LIVER, PANCREAS
   */
  readonly pdbx_gene_src_organ?: Maybe<Scalars['String']>;
  /**
   * Organized structure within cell.
   * 
   * Examples:
   * MITOCHONDRIA
   */
  readonly pdbx_gene_src_organelle?: Maybe<Scalars['String']>;
  /**
   * Scientific name of the organism.
   * 
   * Examples:
   * Homo sapiens, ESCHERICHIA COLI
   * HOMO SAPIENS
   * SACCHAROMYCES CEREVISIAE
   */
  readonly pdbx_gene_src_scientific_name?: Maybe<Scalars['String']>;
  /**
   * Identifies the variant.
   * 
   * Examples:
   * DELTAH1DELTATRP
   */
  readonly pdbx_gene_src_variant?: Maybe<Scalars['String']>;
  /**
   * Americal Tissue Culture Collection of the expression system. Where
   *  full details of the protein production are available it would
   *  be expected that this item  would be derived from
   *  _entity_src_gen_express.host_org_culture_collection
   */
  readonly pdbx_host_org_atcc?: Maybe<Scalars['String']>;
  /**
   * Cell type from which the gene is derived. Where
   *  entity.target_id is provided this should be derived from
   *  details of the target.
   * 
   * Examples:
   * ENDOTHELIAL
   */
  readonly pdbx_host_org_cell?: Maybe<Scalars['String']>;
  /**
   * A specific line of cells used as the expression system. Where
   *  full details of the protein production are available it would
   *  be expected that this item would be derived from
   *  entity_src_gen_express.host_org_cell_line
   * 
   * Examples:
   * HELA
   */
  readonly pdbx_host_org_cell_line?: Maybe<Scalars['String']>;
  /**
   * Identifies the location inside (or outside) the cell which
   *  expressed the molecule.
   * 
   * Examples:
   * CYTOPLASM, NUCLEUS
   */
  readonly pdbx_host_org_cellular_location?: Maybe<Scalars['String']>;
  /**
   * Culture collection of the expression system. Where
   *  full details of the protein production are available it would
   *  be expected that this item  would be derived somehwere, but
   *  exactly where is not clear.
   */
  readonly pdbx_host_org_culture_collection?: Maybe<Scalars['String']>;
  /**
   * Specific gene which expressed the molecule.
   * 
   * Examples:
   * HIV-1 POL, GLNS7, U1A (2-98, Y31H, Q36R)
   */
  readonly pdbx_host_org_gene?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier for the expression system organism.
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4 
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly pdbx_host_org_ncbi_taxonomy_id?: Maybe<Scalars['String']>;
  /**
   * Specific organ which expressed the molecule.
   * 
   * Examples:
   * KIDNEY
   */
  readonly pdbx_host_org_organ?: Maybe<Scalars['String']>;
  /**
   * Specific organelle which expressed the molecule.
   * 
   * Examples:
   * MITOCHONDRIA
   */
  readonly pdbx_host_org_organelle?: Maybe<Scalars['String']>;
  /**
   * The scientific name of the organism that served as host for the
   *  production of the entity. Where full details of the protein
   *  production are available it would be expected that this item
   *  would be derived from _entity_src_gen_express.host_org_scientific_name
   *  or via _entity_src_gen_express.host_org_tax_id
   * 
   * Examples:
   * ESCHERICHIA COLI, SACCHAROMYCES CEREVISIAE
   */
  readonly pdbx_host_org_scientific_name?: Maybe<Scalars['String']>;
  /**
   * The strain of the organism in which the entity was
   * expressed.
   * 
   * Examples:
   * AR120
   */
  readonly pdbx_host_org_strain?: Maybe<Scalars['String']>;
  /**
   * The specific tissue which expressed the molecule. Where full details
   *  of the protein production are available it would be expected that this
   *  item would be derived from _entity_src_gen_express.host_org_tissue
   * 
   * Examples:
   * heart, liver, eye lens
   */
  readonly pdbx_host_org_tissue?: Maybe<Scalars['String']>;
  /**
   * The fraction of the tissue which expressed the
   * molecule.
   * 
   * Examples:
   * mitochondria, nucleus, membrane
   */
  readonly pdbx_host_org_tissue_fraction?: Maybe<Scalars['String']>;
  /**
   * Variant of the organism used as the expression system. Where
   *  full details of the protein production are available it would
   *  be expected that this item be derived from
   *  entity_src_gen_express.host_org_variant or via
   *  _entity_src_gen_express.host_org_tax_id
   * 
   * Examples:
   * TRP-LAC, LAMBDA DE3
   */
  readonly pdbx_host_org_variant?: Maybe<Scalars['String']>;
  /**
   * Identifies the vector used. Where full details of the protein
   *  production are available it would be expected that this item
   *  would be derived from _entity_src_gen_clone.vector_name.
   * 
   * Examples:
   * PBIT36, PET15B, PUC18
   */
  readonly pdbx_host_org_vector?: Maybe<Scalars['String']>;
  /**
   * Identifies the type of vector used (plasmid, virus, or cosmid).
   *  Where full details of the protein production are available it
   *  would be expected that this item would be derived from
   *  _entity_src_gen_express.vector_type.
   * 
   * Examples:
   * COSMID, PLASMID
   */
  readonly pdbx_host_org_vector_type?: Maybe<Scalars['String']>;
  /**
   * This data item povides additional information about the sequence type.
   * 
   * Allowable values:
   * Biological sequence, C-terminal tag, Linker, N-terminal tag
   */
  readonly pdbx_seq_type?: Maybe<Scalars['String']>;
  /** This data item is an ordinal identifier for entity_src_gen data records. */
  readonly pdbx_src_id: Scalars['Int'];
  /**
   * A description of special aspects of the plasmid that produced the
   *  entity in the host organism. Where full details of the protein
   *  production are available it would be expected that this item
   *  would be derived from _pdbx_construct.details of the construct
   *  pointed to from _entity_src_gen_express.plasmid_id.
   */
  readonly plasmid_details?: Maybe<Scalars['String']>;
  /**
   * The name of the plasmid that produced the entity in the host
   *  organism. Where full details of the protein production are available
   *  it would be expected that this item would be derived from
   *  _pdbx_construct.name of the construct pointed to from
   *  _entity_src_gen_express.plasmid_id.
   * 
   * Examples:
   * pET3C, pT123sab
   */
  readonly plasmid_name?: Maybe<Scalars['String']>;
};

export type Symmetry = {
  /**
   * Space-group number from International Tables for Crystallography
   *  Vol. A (2002).
   */
  readonly Int_Tables_number?: Maybe<Scalars['Int']>;
  /**
   * The cell settings for this space-group symmetry.
   * 
   * Allowable values:
   * cubic, hexagonal, monoclinic, orthorhombic, rhombohedral, tetragonal, triclinic, trigonal
   */
  readonly cell_setting?: Maybe<Scalars['String']>;
  /**
   * Used for PDB space group:
   * 
   *  Example: 'C 1 2 1'  (instead of C 2)
   *           'P 1 2 1'  (instead of P 2)
   *           'P 1 21 1' (instead of P 21)
   *           'P 1 1 21' (instead of P 21 -unique C axis)
   *           'H 3'      (instead of R 3   -hexagonal)
   *           'H 3 2'    (instead of R 3 2 -hexagonal)
   * 
   * Examples:
   * Example: 'C 1 2 1'  (instead of C 2)
   *            'P 1 2 1'  (instead of P 2)
   *            'P 1 21 1' (instead of P 21)
   *            'P 1 1 21' (instead of P 21 -unique C axis)
   *            'H 3'      (instead of R 3   -hexagonal)
   *            'H 3 2'    (instead of R 3 2 -hexagonal)
   */
  readonly pdbx_full_space_group_name_H_M?: Maybe<Scalars['String']>;
  /**
   * Hermann-Mauguin space-group symbol. Note that the
   *  Hermann-Mauguin symbol does not necessarily contain complete
   *  information about the symmetry and the space-group origin. If
   *  used, always supply the FULL symbol from International Tables
   *  for Crystallography Vol. A (2002) and indicate the origin and
   *  the setting if it is not implicit. If there is any doubt that
   *  the equivalent positions can be uniquely deduced from this
   *  symbol, specify the  _symmetry_equiv.pos_as_xyz or
   *  _symmetry.space_group_name_Hall  data items as well. Leave
   *  spaces between symbols referring to
   *  different axes.
   * 
   * Examples:
   * A 1, A 1 2 1, A 2, B 1 1 2, B 2, B 2 21 2, C 2, C 1 2 1, C 21, C 1 21 1, C 2(A 112), C 2 2 2, C 2 2 21, C 4 21 2, F 2 2 2, F 2 3, F 4 2 2, F 4 3 2, F 41 3 2, I 1 2 1, I 1 21 1, I 2, I 2 2 2, I 2 3, I 21, I 21 3, I 21 21 21, I 4, I 4 2 2, I 4 3 2, I 41, I 41/a, I 41 2 2, I 41 3 2, P 1, P 1-, P 2, P 1 2 1, P 1 1 2, P 2 2 2, P 2 3, P 2 2 21, P 2 21 21, P 21, P 1 21 1, P 1 21/c 1, P 1 1 21, P 21(C), P 21 2 21, P 21 3, P 21 21 2, P 21 21 2 A, P 21 21 21, P 3, P 3 1 2, P 3 2 1, P 31, P 31 1 2, P 31 2 1, P 32, P 32 1 2, P 32 2 1, P 4, P 4 2 2, P 4 3 2, P 4 21 2, P 41, P 41 2 2, P 41 3 2, P 41 21 2, P 42, P 42 2 2, P 42 3 2, P 42 21 2, P 43, P 43 2 2, P 43 3 2, P 43 21 2, P 6, P 6 2 2, P 61, P 61 2 2, P 62, P 62 2 2, P 63, P 63 2 2, P 64, P 64 2 2, P 65, P 65 2 2, H 3, R 3, H 3 2, R 3 2
   */
  readonly space_group_name_H_M?: Maybe<Scalars['String']>;
  /**
   * Space-group symbol as described by Hall (1981). This symbol
   *  gives the space-group setting explicitly. Leave spaces between
   *  the separate components of the symbol.
   * 
   *  Ref: Hall, S. R. (1981). Acta Cryst. A37, 517-525; erratum
   *  (1981) A37, 921.
   * 
   * Examples:
   * -P 2ac 2n, -R 3 2", P 61 2 2 (0 0 -1)
   */
  readonly space_group_name_Hall?: Maybe<Scalars['String']>;
};

export type PdbxFamilyPrdAudit = {
  /**
   * The action associated with this audit record.
   * 
   * Allowable values:
   * Add PRD, Create family, Initial release, Modify annotation, Modify citation, Modify family classification, Modify family name, Modify feature, Modify molecule details, Modify related structures, Modify sequence, Modify synonyms, Obsolete family, Obsolete familyt, Other modification, Remove PRD
   */
  readonly action_type: Scalars['String'];
  /**
   * The initials of the annotator creating of modifying the family.
   * 
   * Examples:
   * JO, SJ, KB
   */
  readonly annotator?: Maybe<Scalars['String']>;
  /** The date associated with this audit record. */
  readonly date: Scalars['Date'];
  /**
   * Additional details decribing this change.
   * 
   * Examples:
   * Revise molecule sequence.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _pdbx_reference_molecule_family.family_prd_id in the 
   * 	       pdbx_reference_molecule category.
   */
  readonly family_prd_id: Scalars['String'];
  /**
   * An identifier for the wwPDB site creating or modifying the family.
   * 
   * Examples:
   * RCSB, PDBE, PDBJ, BMRB, PDBC
   */
  readonly processing_site?: Maybe<Scalars['String']>;
};

export type EmSoftware = {
  /**
   * The purpose of the software.
   * 
   * Allowable values:
   * CLASSIFICATION, CRYSTALLOGRAPHY MERGING, CTF CORRECTION, DIFFRACTION INDEXING, FINAL EULER ASSIGNMENT, IMAGE ACQUISITION, INITIAL EULER ASSIGNMENT, LATTICE DISTORTION CORRECTION, LAYERLINE INDEXING, MASKING, MODEL FITTING, MODEL REFINEMENT, MOLECULAR REPLACEMENT, OTHER, PARTICLE SELECTION, RECONSTRUCTION, SERIES ALIGNMENT, SYMMETRY DETERMINATION, VOLUME SELECTION
   */
  readonly category?: Maybe<Scalars['String']>;
  /**
   * Details about the software used.
   * 
   * Examples:
   * EMAN2 e2boxer.py was used to automatically select particle images.
   */
  readonly details?: Maybe<Scalars['String']>;
  /** pointer to _em_3d_fitting.id in the EM_3D_FITTING category. */
  readonly fitting_id?: Maybe<Scalars['String']>;
  /** Unique identifier for each software description. */
  readonly id: Scalars['String'];
  /** pointer to _em_image_processing.id in the EM_IMAGE_PROCESSING category. */
  readonly image_processing_id?: Maybe<Scalars['String']>;
  /** pointer to _em_imaging.id in the EM_IMAGING category. */
  readonly imaging_id?: Maybe<Scalars['String']>;
  /**
   * The name of the software package used, e.g., RELION.  Depositors are strongly
   *   encouraged to provide a value in this field.
   * 
   * Examples:
   * EMAN, Imagic, Spider, Bsoft, UCSF-Chimera
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The version of the software.
   * 
   * Examples:
   * 9.03, 2.1
   */
  readonly version?: Maybe<Scalars['String']>;
};

export type RcsbAccessionInfo = {
  /** The entry deposition date. */
  readonly deposit_date?: Maybe<Scalars['Date']>;
  /**
   * A code indicating the current availibility of experimental data in the repository.
   * 
   * Allowable values:
   * N, Y
   */
  readonly has_released_experimental_data?: Maybe<Scalars['String']>;
  /** The entry initial release date. */
  readonly initial_release_date?: Maybe<Scalars['Date']>;
  /** The latest entry major revision number. */
  readonly major_revision?: Maybe<Scalars['Int']>;
  /** The latest entry minor revision number. */
  readonly minor_revision?: Maybe<Scalars['Int']>;
  /** The latest entry revision date. */
  readonly revision_date?: Maybe<Scalars['Date']>;
  /**
   * The release status for the entry.
   * 
   * Allowable values:
   * AUCO, AUTH, HOLD, HPUB', POLC, PROC, REFI, REL, REPL, WAIT, WDRN
   */
  readonly status_code?: Maybe<Scalars['String']>;
};

export type DiffrnRadiation = {
  /**
   * The collimation or focusing applied to the radiation.
   * 
   * Examples:
   * 0.3 mm double-pinhole, 0.5 mm, focusing mirrors
   */
  readonly collimation?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * The method used to obtain monochromatic radiation. If a mono-
   *  chromator crystal is used, the material and the indices of the
   *  Bragg reflection are specified.
   * 
   * Examples:
   * Zr filter, Ge 220, none, equatorial mounted graphite
   */
  readonly monochromator?: Maybe<Scalars['String']>;
  /**
   * SINGLE WAVELENGTH, LAUE, or MAD.
   * 
   * Examples:
   * SINGLE WAVELENGTH, MONOCHROMATIC, LAUE, MAD, OTHER
   */
  readonly pdbx_diffrn_protocol?: Maybe<Scalars['String']>;
  /**
   * Monochromatic or Laue.
   * 
   * Allowable values:
   * L, M
   */
  readonly pdbx_monochromatic_or_laue_m_l?: Maybe<Scalars['String']>;
  /**
   * The radiation scattering type for this diffraction data set.
   * 
   * Allowable values:
   * electron, neutron, x-ray
   */
  readonly pdbx_scattering_type?: Maybe<Scalars['String']>;
  /** Wavelength of radiation. */
  readonly pdbx_wavelength?: Maybe<Scalars['String']>;
  /** Comma separated list of wavelengths or wavelength range. */
  readonly pdbx_wavelength_list?: Maybe<Scalars['String']>;
  /**
   * The nature of the radiation. This is typically a description
   *  of the X-ray wavelength in Siegbahn notation.
   * 
   * Examples:
   * CuK\a, Cu K\a~1~, Cu K-L~2,3~, white-beam
   */
  readonly type?: Maybe<Scalars['String']>;
  /**
   * This data item is a pointer to _diffrn_radiation_wavelength.id
   *  in the DIFFRN_RADIATION_WAVELENGTH category.
   */
  readonly wavelength_id?: Maybe<Scalars['String']>;
};

export type RcsbNonpolymerInstanceFeatureSummary = {
  /** Component identifier for non-polymer entity instance. */
  readonly comp_id?: Maybe<Scalars['String']>;
  /** The feature count. */
  readonly count?: Maybe<Scalars['Int']>;
  /** The maximum feature length. */
  readonly maximum_length?: Maybe<Scalars['Int']>;
  /** The maximum feature value. */
  readonly maximum_value?: Maybe<Scalars['Float']>;
  /** The minimum feature length. */
  readonly minimum_length?: Maybe<Scalars['Int']>;
  /** The minimum feature value. */
  readonly minimum_value?: Maybe<Scalars['Float']>;
  /**
   * Type or category of the feature.
   * 
   * Allowable values:
   * HAS_COVALENT_LINKAGE, HAS_METAL_COORDINATION_LINKAGE, MOGUL_ANGLE_OUTLIER, MOGUL_BOND_OUTLIER, RSRCC_OUTLIER, RSRZ_OUTLIER
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbClusterMembership = {
  /** Identifier for a cluster at the specified level of sequence identity within the cluster data set. */
  readonly cluster_id?: Maybe<Scalars['Int']>;
  /** Sequence identity expressed as an integer percent value. */
  readonly identity?: Maybe<Scalars['Int']>;
};


export type EntitySrcNat = {
  /**
   * The common name of the organism from which the entity
   *  was isolated.
   * 
   * Examples:
   * man, yeast, bacteria
   */
  readonly common_name?: Maybe<Scalars['String']>;
  /**
   * A description of special aspects of the organism from which the
   *  entity was isolated.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The genus of the organism from which the entity was isolated.
   * 
   * Examples:
   * Homo, Saccharomyces, Escherichia
   */
  readonly genus?: Maybe<Scalars['String']>;
  /**
   * This data item identifies cases in which an alternative source
   *  modeled.
   * 
   * Allowable values:
   * model, sample
   */
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>;
  /**
   * Americal Tissue Culture Collection number.
   * 
   * Examples:
   * 6051
   */
  readonly pdbx_atcc?: Maybe<Scalars['String']>;
  /**
   * The beginning polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>;
  /**
   * A particular cell type.
   * 
   * Examples:
   * BHK-21
   */
  readonly pdbx_cell?: Maybe<Scalars['String']>;
  /**
   * The specific line of cells.
   * 
   * Examples:
   * HELA
   */
  readonly pdbx_cell_line?: Maybe<Scalars['String']>;
  /** Identifies the location inside (or outside) the cell. */
  readonly pdbx_cellular_location?: Maybe<Scalars['String']>;
  /**
   * The ending polymer sequence position for the polymer section corresponding  
   *  to this source.
   * 
   *  A reference to the sequence position in the entity_poly category.
   */
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>;
  /** A domain or fragment of the molecule. */
  readonly pdbx_fragment?: Maybe<Scalars['String']>;
  /**
   * NCBI Taxonomy identifier for the source organism.
   * 
   *  Reference:
   * 
   *  Wheeler DL, Chappey C, Lash AE, Leipe DD, Madden TL, Schuler GD,
   *  Tatusova TA, Rapp BA (2000). Database resources of the National
   *  Center for Biotechnology Information. Nucleic Acids Res 2000 Jan
   *  1;28(1):10-4 
   * 
   *  Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA,
   *  Wheeler DL (2000). GenBank. Nucleic Acids Res 2000 Jan 1;28(1):15-18.
   */
  readonly pdbx_ncbi_taxonomy_id?: Maybe<Scalars['String']>;
  /**
   * Organized group of tissues that carries on a specialized function.
   * 
   * Examples:
   * KIDNEY
   */
  readonly pdbx_organ?: Maybe<Scalars['String']>;
  /**
   * Organized structure within cell.
   * 
   * Examples:
   * MITOCHONDRIA
   */
  readonly pdbx_organelle?: Maybe<Scalars['String']>;
  /**
   * Scientific name of the organism of the natural source.
   * 
   * Examples:
   * Bos taurus, BOS TAURUS, SUS SCROFA, ASPERGILLUS ORYZAE
   */
  readonly pdbx_organism_scientific?: Maybe<Scalars['String']>;
  /**
   * Details about the plasmid.
   * 
   * Examples:
   * PLC28 DERIVATIVE
   */
  readonly pdbx_plasmid_details?: Maybe<Scalars['String']>;
  /**
   * The plasmid containing the gene.
   * 
   * Examples:
   * pB322
   */
  readonly pdbx_plasmid_name?: Maybe<Scalars['String']>;
  /**
   * Identifies the secretion from which the molecule was isolated.
   * 
   * Examples:
   * saliva, urine, venom
   */
  readonly pdbx_secretion?: Maybe<Scalars['String']>;
  /** This data item is an ordinal identifier for entity_src_nat data records. */
  readonly pdbx_src_id: Scalars['Int'];
  /** Identifies the variant. */
  readonly pdbx_variant?: Maybe<Scalars['String']>;
  /**
   * The species of the organism from which the entity was isolated.
   * 
   * Examples:
   * sapiens, cerevisiae, coli
   */
  readonly species?: Maybe<Scalars['String']>;
  /**
   * The strain of the organism from which the entity was isolated.
   * 
   * Examples:
   * DH5a, BMH 71-18
   */
  readonly strain?: Maybe<Scalars['String']>;
  /**
   * The tissue of the organism from which the entity was isolated.
   * 
   * Examples:
   * heart, liver, eye lens
   */
  readonly tissue?: Maybe<Scalars['String']>;
  /**
   * The subcellular fraction of the tissue of the organism from
   *  which the entity was isolated.
   * 
   * Examples:
   * mitochondria, nucleus, membrane
   */
  readonly tissue_fraction?: Maybe<Scalars['String']>;
};

export type EmImaging = {
  /**
   * A value of accelerating voltage (in kV) used for imaging.
   * 
   * Examples:
   * 300
   */
  readonly accelerating_voltage?: Maybe<Scalars['Int']>;
  /**
   * microscope alignment procedure
   * 
   * Allowable values:
   * BASIC, COMA FREE, NONE, OTHER, ZEMLIN TABLEAU
   */
  readonly alignment_procedure?: Maybe<Scalars['String']>;
  /** astigmatism */
  readonly astigmatism?: Maybe<Scalars['String']>;
  /**
   * The open diameter of the c2 condenser lens,
   *  in microns.
   */
  readonly c2_aperture_diameter?: Maybe<Scalars['Float']>;
  /**
   * The maximum defocus value of the objective lens (in nanometers) used
   *  to obtain the recorded images.
   * 
   * Examples:
   * 5000
   */
  readonly calibrated_defocus_max?: Maybe<Scalars['Float']>;
  /**
   * The minimum defocus value of the objective lens (in nanometers) used
   *  to obtain the recorded images.
   * 
   * Examples:
   * 1200
   */
  readonly calibrated_defocus_min?: Maybe<Scalars['Float']>;
  /**
   * The magnification value obtained for a known standard just
   *  prior to, during or just after the imaging experiment.
   * 
   * Examples:
   * 61200
   */
  readonly calibrated_magnification?: Maybe<Scalars['Int']>;
  /**
   * Cryogen type used to maintain the specimen stage temperature during imaging
   *  in the microscope.
   * 
   * Allowable values:
   * HELIUM, NITROGEN
   */
  readonly cryogen?: Maybe<Scalars['String']>;
  /**
   * Date (YYYY-MM-DD) of imaging experiment or the date at which
   *  a series of experiments began.
   * 
   * Examples:
   * 2001-05-08
   */
  readonly date?: Maybe<Scalars['Date']>;
  /**
   * Any additional imaging details.
   * 
   * Examples:
   * Preliminary grid screening was performed manually.
   */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * The camera length (in millimeters). The camera length is the
   *  product of the objective focal length and the combined magnification
   *  of the intermediate and projector lenses when the microscope is
   *  operated in the diffraction mode.
   */
  readonly detector_distance?: Maybe<Scalars['Float']>;
  /** electron beam tilt params */
  readonly electron_beam_tilt_params?: Maybe<Scalars['String']>;
  /** The source of electrons. The electron gun. */
  readonly electron_source?: Maybe<Scalars['String']>;
  /**
   * The value of _em_imaging.id must uniquely identify
   *  each imaging experiment.
   */
  readonly id: Scalars['String'];
  /**
   * The mode of illumination.
   * 
   * Allowable values:
   * FLOOD BEAM, OTHER, SPOT SCAN
   */
  readonly illumination_mode?: Maybe<Scalars['String']>;
  /**
   * The name of the model of microscope.
   * 
   * Allowable values:
   * FEI MORGAGNI, FEI POLARA 300, FEI TALOS ARCTICA, FEI TECNAI 10, FEI TECNAI 12, FEI TECNAI 20, FEI TECNAI ARCTICA, FEI TECNAI F20, FEI TECNAI F30, FEI TECNAI SPHERA, FEI TECNAI SPIRIT, FEI TITAN, FEI TITAN KRIOS, FEI/PHILIPS CM10, FEI/PHILIPS CM12, FEI/PHILIPS CM120T, FEI/PHILIPS CM200FEG, FEI/PHILIPS CM200FEG/SOPHIE, FEI/PHILIPS CM200FEG/ST, FEI/PHILIPS CM200FEG/UT, FEI/PHILIPS CM200T, FEI/PHILIPS CM300FEG/HE, FEI/PHILIPS CM300FEG/ST, FEI/PHILIPS CM300FEG/T, FEI/PHILIPS EM400, FEI/PHILIPS EM420, HITACHI EF2000, HITACHI EF3000, HITACHI H-9500SD, HITACHI H3000 UHVEM, HITACHI H7600, HITACHI HF2000, HITACHI HF3000, JEOL 100B, JEOL 100CX, JEOL 1010, JEOL 1200, JEOL 1200EX, JEOL 1200EXII, JEOL 1230, JEOL 1400, JEOL 2000EX, JEOL 2000EXII, JEOL 2010, JEOL 2010F, JEOL 2010HC, JEOL 2010HT, JEOL 2010UHR, JEOL 2011, JEOL 2100, JEOL 2100F, JEOL 2200FS, JEOL 2200FSC, JEOL 3000SFF, JEOL 3100FEF, JEOL 3100FFC, JEOL 3200FS, JEOL 3200FSC, JEOL 4000, JEOL 4000EX, JEOL CRYO ARM 200, JEOL CRYO ARM 300, JEOL KYOTO-3000SFF, SIEMENS SULEIKA, TFS GLACIOS, TFS KRIOS, TFS TALOS, TFS TALOS F200C, TFS TALOS L120C, ZEISS LEO912, ZEISS LIBRA120PLUS
   */
  readonly microscope_model?: Maybe<Scalars['String']>;
  /**
   * The mode of imaging.
   * 
   * Allowable values:
   * BRIGHT FIELD, DARK FIELD, DIFFRACTION, OTHER
   */
  readonly mode?: Maybe<Scalars['String']>;
  /**
   * The spherical aberration coefficient (Cs) in millimeters,
   *  of the objective lens.
   * 
   * Examples:
   * 2.0
   */
  readonly nominal_cs?: Maybe<Scalars['Float']>;
  /**
   * The maximum defocus value of the objective lens (in nanometers) used
   *  to obtain the recorded images.
   * 
   * Examples:
   * 5000
   */
  readonly nominal_defocus_max?: Maybe<Scalars['Float']>;
  /**
   * The minimum defocus value of the objective lens (in nanometers) used
   *  to obtain the recorded images.
   * 
   * Examples:
   * 1200
   */
  readonly nominal_defocus_min?: Maybe<Scalars['Float']>;
  /**
   * The magnification indicated by the microscope readout.
   * 
   * Examples:
   * 60000
   */
  readonly nominal_magnification?: Maybe<Scalars['Int']>;
  /**
   * The specimen temperature maximum (degrees Kelvin) for the duration
   *  of imaging.
   */
  readonly recording_temperature_maximum?: Maybe<Scalars['Float']>;
  /**
   * The specimen temperature minimum (degrees Kelvin) for the duration
   *  of imaging.
   */
  readonly recording_temperature_minimum?: Maybe<Scalars['Float']>;
  /** residual tilt of the electron beam */
  readonly residual_tilt?: Maybe<Scalars['Float']>;
  /**
   * The name of the model of specimen holder used during imaging.
   * 
   * Allowable values:
   * FEI TITAN KRIOS AUTOGRID HOLDER, FISCHIONE 2550, FISCHIONE INSTRUMENTS DUAL AXIS TOMOGRAPHY HOLDER, GATAN 626 SINGLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER, GATAN 910 MULTI-SPECIMEN SINGLE TILT CRYO TRANSFER HOLDER, GATAN 914 HIGH TILT LIQUID NITROGEN CRYO TRANSFER TOMOGRAPHY HOLDER, GATAN 915 DOUBLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER, GATAN CHDT 3504 DOUBLE TILT HIGH RESOLUTION NITROGEN COOLING HOLDER, GATAN CT3500 SINGLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER, GATAN CT3500TR SINGLE TILT ROTATION LIQUID NITROGEN CRYO TRANSFER HOLDER, GATAN ELSA 698 SINGLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER, GATAN HC 3500 SINGLE TILT HEATING/NITROGEN COOLING HOLDER, GATAN HCHDT 3010 DOUBLE TILT HIGH RESOLUTION HELIUM COOLING HOLDER, GATAN HCHST 3008 SINGLE TILT HIGH RESOLUTION HELIUM COOLING HOLDER, GATAN HELIUM, GATAN LIQUID NITROGEN, GATAN UHRST 3500 SINGLE TILT ULTRA HIGH RESOLUTION NITROGEN COOLING HOLDER, GATAN ULTDT ULTRA LOW TEMPERATURE DOUBLE TILT HELIUM COOLING HOLDER, GATAN ULTST ULTRA LOW TEMPERATURE SINGLE TILT HELIUM COOLING HOLDER, HOME BUILD, JEOL, JEOL 3200FSC CRYOHOLDER, JEOL CRYOSPECPORTER, OTHER, PHILIPS ROTATION HOLDER, SIDE ENTRY, EUCENTRIC
   */
  readonly specimen_holder_model?: Maybe<Scalars['String']>;
  /**
   * The type of specimen holder used during imaging.
   * 
   * Examples:
   * cryo
   */
  readonly specimen_holder_type?: Maybe<Scalars['String']>;
  /** Foreign key to the EM_SPECIMEN category */
  readonly specimen_id?: Maybe<Scalars['String']>;
  /**
   * The mean specimen stage temperature (degrees Kelvin) during imaging
   *  in the microscope.
   * 
   * Examples:
   * 70
   */
  readonly temperature?: Maybe<Scalars['Float']>;
  /**
   * The maximum angle at which the specimen was tilted to obtain
   *  recorded images.
   * 
   * Examples:
   * 70
   */
  readonly tilt_angle_max?: Maybe<Scalars['Float']>;
  /**
   * The minimum angle at which the specimen was tilted to obtain
   *  recorded images.
   * 
   * Examples:
   * -70
   */
  readonly tilt_angle_min?: Maybe<Scalars['Float']>;
};

export type PdbxSerialCrystallographyMeasurement = {
  /**
   * The total number of hours required to measure this data set.
   * 
   * Examples:
   * 120.0
   */
  readonly collection_time_total?: Maybe<Scalars['Float']>;
  /**
   * The collimation or type of focusing optics applied to the radiation.
   * 
   * Examples:
   * Kirkpatrick-Baez mirrors, Beryllium compound refractive lenses, Fresnel zone plates
   */
  readonly collimation?: Maybe<Scalars['String']>;
  /**
   * The data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   * 
   * Examples:
   * 1
   */
  readonly diffrn_id: Scalars['String'];
  /**
   * The focal spot size of the beam
   *  impinging on the sample (micrometres squared).
   */
  readonly focal_spot_size?: Maybe<Scalars['Float']>;
  /** The photons per pulse measured in  (tera photons (10^(12)^)/pulse units). */
  readonly photons_per_pulse?: Maybe<Scalars['Float']>;
  /**
   * The average duration (femtoseconds)
   * 	       of the pulse energy measured at the sample.
   */
  readonly pulse_duration?: Maybe<Scalars['Float']>;
  /** The energy/pulse of the X-ray pulse impacting the sample measured in microjoules. */
  readonly pulse_energy?: Maybe<Scalars['Float']>;
  /** The photon energy of the X-ray pulse measured in KeV. */
  readonly pulse_photon_energy?: Maybe<Scalars['Float']>;
  /** The distance from source to the sample along the optical axis (metres). */
  readonly source_distance?: Maybe<Scalars['Float']>;
  /** The dimension of the source beam measured at the source (micrometres squared). */
  readonly source_size?: Maybe<Scalars['Float']>;
  /** For FEL experiments, the pulse repetition rate measured in cycles per seconds. */
  readonly xfel_pulse_repetition_rate?: Maybe<Scalars['Float']>;
};

export type RcsbRepositoryHoldingsCurrent = {
  /**
   * The list of content types associated with this entry.
   * 
   * Allowable values:
   * 2fo-fc Map, Combined NMR data (NEF), Combined NMR data (NMR-STAR), FASTA sequence, Map Coefficients, NMR chemical shifts, NMR restraints V1, NMR restraints V2, assembly PDB, assembly mmCIF, entry PDB, entry PDB bundle, entry PDBML, entry mmCIF, fo-fc Map, structure factors, validation report
   */
  readonly repository_content_types?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};


export type PdbxMoleculeFeatures = {
  /**
   * Broadly defines the function of the molecule.
   * 
   * Allowable values:
   * Antagonist, Anthelmintic, Antibiotic, Antibiotic, Anthelmintic, Antibiotic, Antimicrobial, Antibiotic, Antineoplastic, Anticancer, Anticoagulant, Anticoagulant, Antithrombotic, Antifungal, Antiinflammatory, Antimicrobial, Antimicrobial, Antiparasitic, Antibiotic, Antimicrobial, Antiretroviral, Antimicrobial, Antitumor, Antineoplastic, Antiparasitic, Antiretroviral, Antithrombotic, Antitumor, Antiviral, CASPASE inhibitor, Chaperone binding, Enzyme inhibitor, Growth factor, Immunosuppressant, Inhibitor, Lantibiotic, Metabolism, Metal transport, Oxidation-reduction, Receptor, Thrombin inhibitor, Thrombin inhibitor, Trypsin inhibitor, Toxin, Transport activator, Trypsin inhibitor, Unknown
   */
  readonly class?: Maybe<Scalars['String']>;
  /** Additional details describing the molecule. */
  readonly details?: Maybe<Scalars['String']>;
  /**
   * A name of the molecule.
   * 
   * Examples:
   * thiostrepton
   */
  readonly name?: Maybe<Scalars['String']>;
  /**
   * The value of _pdbx_molecule_features.prd_id is the PDB accession code for this 
   *  reference molecule.
   */
  readonly prd_id: Scalars['String'];
  /**
   * Defines the structural classification of the molecule.
   * 
   * Allowable values:
   * Amino acid, Aminoglycoside, Ansamycin, Anthracycline, Anthraquinone, Chalkophore, Chalkophore, Polypeptide, Chromophore, Cyclic depsipeptide, Cyclic lipopeptide, Cyclic peptide, Glycopeptide, Heterocyclic, Imino sugar, Keto acid, Lipoglycopeptide, Lipopeptide, Macrolide, Non-polymer, Nucleoside, Oligopeptide, Oligosaccharide, Peptaibol, Peptide-like, Polycyclic, Polypeptide, Polysaccharide, Quinolone, Siderophore, Thiolactone, Thiopeptide, Unknown
   */
  readonly type?: Maybe<Scalars['String']>;
};

export type RcsbEntryContainerIdentifiers = {
  /** List of identifiers for assemblies generated from the entry. */
  readonly assembly_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** List of identifiers for the branched entity constituents for the entry. */
  readonly branched_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /**
   * List of EMDB identifiers for the 3D electron microscopy density maps
   *  used in the production of the structure model.
   */
  readonly emdb_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** List of identifiers or the entity constituents for the entry. */
  readonly entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Entry identifier for the container. */
  readonly entry_id: Scalars['String'];
  /** List of PDB model identifiers for the entry. */
  readonly model_ids?: Maybe<ReadonlyArray<Maybe<Scalars['Int']>>>;
  /** List of identifiers for the non-polymer entity constituents for the entry. */
  readonly non_polymer_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** List of identifiers for the polymer entity constituents for the entry. */
  readonly polymer_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** Unique integer value assigned to each PubMed record. */
  readonly pubmed_id?: Maybe<Scalars['Int']>;
  /** A unique identifier for each object in this entry container. */
  readonly rcsb_id?: Maybe<Scalars['String']>;
  /**
   * List of EMDB identifiers for the 3D electron microscopy density maps
   *  related to the structure model.
   */
  readonly related_emdb_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
  /** List of identifiers for the solvent/water entity constituents for the entry. */
  readonly water_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>;
};

export type PdbxSerialCrystallographySampleDeliveryInjection = {
  /**
   * For continuous sample flow experiments, the carrier buffer used
   *  to move the sample into the beam. Should include protein
   *  concentration.
   * 
   * Examples:
   * LCP, grease, liquid
   */
  readonly carrier_solvent?: Maybe<Scalars['String']>;
  /**
   * For continuous sample flow experiments, the concentration of
   *  crystals in the solution being injected.
   * 
   *  The concentration is measured in million crystals/ml.
   */
  readonly crystal_concentration?: Maybe<Scalars['Float']>;
  /**
   * For continuous sample flow experiments, a description of the injector used
   *  to move the sample into the beam.
   * 
   * Examples:
   * microextrusion injector
   */
  readonly description?: Maybe<Scalars['String']>;
  /**
   * The data item is a pointer to _diffrn.id in the DIFFRN
   *  category.
   * 
   * Examples:
   * 1
   */
  readonly diffrn_id: Scalars['String'];
  /** The size of filter in micrometres in filtering crystals */
  readonly filter_size?: Maybe<Scalars['Float']>;
  /**
   * For continuous sample flow experiments, the flow rate of
   *  solution being injected  measured in ul/min.
   */
  readonly flow_rate?: Maybe<Scalars['Float']>;
  /**
   * For continuous sample flow experiments, the diameter of the
   *  injector in micrometres.
   */
  readonly injector_diameter?: Maybe<Scalars['Float']>;
  /**
   * The type of nozzle to deliver and focus sample jet
   * 
   * Examples:
   * gas, GDVN
   */
  readonly injector_nozzle?: Maybe<Scalars['String']>;
  /**
   * For continuous sample flow experiments, the mean pressure
   *  in kilopascals at which the sample is injected into the beam.
   */
  readonly injector_pressure?: Maybe<Scalars['Float']>;
  /**
   * For continuous sample flow experiments, the temperature in
   *  Kelvins of the speciman injected. This may be different from
   *  the temperature of the sample.
   */
  readonly injector_temperature?: Maybe<Scalars['Float']>;
  /** Diameter in micrometres of jet stream of sample delivery */
  readonly jet_diameter?: Maybe<Scalars['Float']>;
  /**
   * Sample deliver driving force, e.g. Gas, Electronic Potential
   * 
   * Examples:
   * syringe, gas, electronic potential
   */
  readonly power_by?: Maybe<Scalars['String']>;
  /**
   * Details of crystal growth and preparation of the crystals
   * 
   * Examples:
   * Crystals transfered to carrier solvent at room temperature
   */
  readonly preparation?: Maybe<Scalars['String']>;
};

export type RcsbPolymerInstanceFeatureFeaturePositions = {
  /** An identifier for the monomer(s) corresponding to the feature assignment. */
  readonly beg_comp_id?: Maybe<Scalars['String']>;
  /** An identifier for the monomer at which this segment of the feature begins. */
  readonly beg_seq_id: Scalars['Int'];
  /** An identifier for the monomer at which this segment of the feature ends. */
  readonly end_seq_id?: Maybe<Scalars['Int']>;
  /** The value of the feature over the monomer segment. */
  readonly value?: Maybe<Scalars['Float']>;
};

export type AssemblySymmetryQueryVariables = {
  assembly_id: Scalars['String'];
  entry_id: Scalars['String'];
};


export type AssemblySymmetryQuery = { readonly assembly?: Maybe<{ readonly rcsb_struct_symmetry?: Maybe<ReadonlyArray<Maybe<(
      Pick<RcsbStructSymmetry, 'kind' | 'oligomeric_state' | 'stoichiometry' | 'symbol' | 'type'>
      & { readonly clusters: ReadonlyArray<Maybe<(
        Pick<RcsbStructSymmetryClusters, 'avg_rmsd'>
        & { readonly members: ReadonlyArray<Maybe<Pick<ClustersMembers, 'asym_id' | 'pdbx_struct_oper_list_ids'>>> }
      )>>, readonly rotation_axes?: Maybe<ReadonlyArray<Maybe<Pick<RcsbStructSymmetryRotationAxes, 'order' | 'start' | 'end'>>>> }
    )>>> }> };
