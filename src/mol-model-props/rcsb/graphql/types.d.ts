/* eslint-disable */
export type Maybe<T> = T | null;

// Generated in 2020-02-05T13:46:39-08:00

/** All built-in and custom scalars, mapped to their actual values */
export type Scalars = {
  ID: string,
  String: string,
  Boolean: boolean,
  Int: number,
  Float: number,
  Date: any,
  UNREPRESENTABLE: any,
};


export type AuditAuthor = {
  readonly __typename?: 'AuditAuthor',
  readonly identifier_ORCID?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly pdbx_ordinal: Scalars['Int'],
};

export type Cell = {
  readonly __typename?: 'Cell',
  readonly Z_PDB?: Maybe<Scalars['Int']>,
  readonly angle_alpha?: Maybe<Scalars['Float']>,
  readonly angle_beta?: Maybe<Scalars['Float']>,
  readonly angle_gamma?: Maybe<Scalars['Float']>,
  readonly formula_units_Z?: Maybe<Scalars['Int']>,
  readonly length_a?: Maybe<Scalars['Float']>,
  readonly length_b?: Maybe<Scalars['Float']>,
  readonly length_c?: Maybe<Scalars['Float']>,
  readonly pdbx_unique_axis?: Maybe<Scalars['String']>,
  readonly volume?: Maybe<Scalars['Float']>,
};

export type ChemComp = {
  readonly __typename?: 'ChemComp',
  readonly formula?: Maybe<Scalars['String']>,
  readonly formula_weight?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly mon_nstd_parent_comp_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly name?: Maybe<Scalars['String']>,
  readonly one_letter_code?: Maybe<Scalars['String']>,
  readonly pdbx_ambiguous_flag?: Maybe<Scalars['String']>,
  readonly pdbx_formal_charge?: Maybe<Scalars['Int']>,
  readonly pdbx_initial_date?: Maybe<Scalars['Date']>,
  readonly pdbx_modified_date?: Maybe<Scalars['Date']>,
  readonly pdbx_processing_site?: Maybe<Scalars['String']>,
  readonly pdbx_release_status?: Maybe<Scalars['String']>,
  readonly pdbx_replaced_by?: Maybe<Scalars['String']>,
  readonly pdbx_replaces?: Maybe<Scalars['String']>,
  readonly pdbx_subcomponent_list?: Maybe<Scalars['String']>,
  readonly three_letter_code?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type Citation = {
  readonly __typename?: 'Citation',
  readonly book_id_ISBN?: Maybe<Scalars['String']>,
  readonly book_publisher?: Maybe<Scalars['String']>,
  readonly book_publisher_city?: Maybe<Scalars['String']>,
  readonly book_title?: Maybe<Scalars['String']>,
  readonly coordinate_linkage?: Maybe<Scalars['String']>,
  readonly country?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly journal_abbrev?: Maybe<Scalars['String']>,
  readonly journal_id_ASTM?: Maybe<Scalars['String']>,
  readonly journal_id_CSD?: Maybe<Scalars['String']>,
  readonly journal_id_ISSN?: Maybe<Scalars['String']>,
  readonly journal_issue?: Maybe<Scalars['String']>,
  readonly journal_volume?: Maybe<Scalars['String']>,
  readonly language?: Maybe<Scalars['String']>,
  readonly page_first?: Maybe<Scalars['String']>,
  readonly page_last?: Maybe<Scalars['String']>,
  readonly pdbx_database_id_DOI?: Maybe<Scalars['String']>,
  readonly pdbx_database_id_PubMed?: Maybe<Scalars['Int']>,
  readonly rcsb_authors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_is_primary?: Maybe<Scalars['String']>,
  readonly rcsb_journal_abbrev?: Maybe<Scalars['String']>,
  readonly title?: Maybe<Scalars['String']>,
  readonly unpublished_flag?: Maybe<Scalars['String']>,
  readonly year?: Maybe<Scalars['Int']>,
};

export type ClustersMembers = {
  readonly __typename?: 'ClustersMembers',
  readonly asym_id: Scalars['String'],
  readonly pdbx_struct_oper_list_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type CoreAssembly = {
  readonly __typename?: 'CoreAssembly',
  readonly entry?: Maybe<CoreEntry>,
  readonly pdbx_struct_assembly?: Maybe<PdbxStructAssembly>,
  readonly pdbx_struct_assembly_auth_evidence?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyAuthEvidence>>>,
  readonly pdbx_struct_assembly_gen?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyGen>>>,
  readonly pdbx_struct_assembly_prop?: Maybe<ReadonlyArray<Maybe<PdbxStructAssemblyProp>>>,
  readonly pdbx_struct_oper_list?: Maybe<ReadonlyArray<Maybe<PdbxStructOperList>>>,
  readonly rcsb_assembly_container_identifiers: RcsbAssemblyContainerIdentifiers,
  readonly rcsb_assembly_info?: Maybe<RcsbAssemblyInfo>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>,
  readonly rcsb_struct_symmetry?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetry>>>,
  readonly rcsb_struct_symmetry_lineage?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetryLineage>>>,
  readonly rcsb_struct_symmetry_provenance_code?: Maybe<Scalars['String']>,
};

export type CoreChemComp = {
  readonly __typename?: 'CoreChemComp',
  readonly chem_comp?: Maybe<ChemComp>,
  readonly drugbank?: Maybe<CoreDrugbank>,
  readonly pdbx_chem_comp_audit?: Maybe<ReadonlyArray<Maybe<PdbxChemCompAudit>>>,
  readonly pdbx_chem_comp_descriptor?: Maybe<ReadonlyArray<Maybe<PdbxChemCompDescriptor>>>,
  readonly pdbx_chem_comp_feature?: Maybe<ReadonlyArray<Maybe<PdbxChemCompFeature>>>,
  readonly pdbx_chem_comp_identifier?: Maybe<ReadonlyArray<Maybe<PdbxChemCompIdentifier>>>,
  readonly pdbx_family_prd_audit?: Maybe<ReadonlyArray<Maybe<PdbxFamilyPrdAudit>>>,
  readonly pdbx_prd_audit?: Maybe<ReadonlyArray<Maybe<PdbxPrdAudit>>>,
  readonly pdbx_reference_entity_list?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityList>>>,
  readonly pdbx_reference_entity_poly?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPoly>>>,
  readonly pdbx_reference_entity_poly_link?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPolyLink>>>,
  readonly pdbx_reference_entity_poly_seq?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntityPolySeq>>>,
  readonly pdbx_reference_entity_sequence?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntitySequence>>>,
  readonly pdbx_reference_entity_src_nat?: Maybe<ReadonlyArray<Maybe<PdbxReferenceEntitySrcNat>>>,
  readonly pdbx_reference_molecule?: Maybe<PdbxReferenceMolecule>,
  readonly pdbx_reference_molecule_annotation?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeAnnotation>>>,
  readonly pdbx_reference_molecule_details?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeDetails>>>,
  readonly pdbx_reference_molecule_family?: Maybe<PdbxReferenceMoleculeFamily>,
  readonly pdbx_reference_molecule_features?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeFeatures>>>,
  readonly pdbx_reference_molecule_list?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeList>>>,
  readonly pdbx_reference_molecule_related_structures?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeRelatedStructures>>>,
  readonly pdbx_reference_molecule_synonyms?: Maybe<ReadonlyArray<Maybe<PdbxReferenceMoleculeSynonyms>>>,
  readonly rcsb_bird_citation?: Maybe<ReadonlyArray<Maybe<RcsbBirdCitation>>>,
  readonly rcsb_chem_comp_container_identifiers?: Maybe<RcsbChemCompContainerIdentifiers>,
  readonly rcsb_chem_comp_descriptor?: Maybe<RcsbChemCompDescriptor>,
  readonly rcsb_chem_comp_info?: Maybe<RcsbChemCompInfo>,
  readonly rcsb_chem_comp_related?: Maybe<ReadonlyArray<Maybe<RcsbChemCompRelated>>>,
  readonly rcsb_chem_comp_synonyms?: Maybe<ReadonlyArray<Maybe<RcsbChemCompSynonyms>>>,
  readonly rcsb_chem_comp_target?: Maybe<ReadonlyArray<Maybe<RcsbChemCompTarget>>>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_schema_container_identifiers?: Maybe<ReadonlyArray<Maybe<RcsbSchemaContainerIdentifiers>>>,
};

export type CoreDrugbank = {
  readonly __typename?: 'CoreDrugbank',
  readonly drugbank_container_identifiers?: Maybe<DrugbankContainerIdentifiers>,
  readonly drugbank_info?: Maybe<DrugbankInfo>,
  readonly drugbank_target?: Maybe<ReadonlyArray<Maybe<DrugbankTarget>>>,
};

export type CoreEntry = {
  readonly __typename?: 'CoreEntry',
  readonly assemblies?: Maybe<ReadonlyArray<Maybe<CoreAssembly>>>,
  readonly audit_author?: Maybe<ReadonlyArray<Maybe<AuditAuthor>>>,
  readonly cell?: Maybe<Cell>,
  readonly citation?: Maybe<ReadonlyArray<Maybe<Citation>>>,
  readonly diffrn?: Maybe<ReadonlyArray<Maybe<Diffrn>>>,
  readonly diffrn_detector?: Maybe<ReadonlyArray<Maybe<DiffrnDetector>>>,
  readonly diffrn_radiation?: Maybe<ReadonlyArray<Maybe<DiffrnRadiation>>>,
  readonly diffrn_source?: Maybe<ReadonlyArray<Maybe<DiffrnSource>>>,
  readonly em_2d_crystal_entity?: Maybe<ReadonlyArray<Maybe<Em2dCrystalEntity>>>,
  readonly em_3d_crystal_entity?: Maybe<ReadonlyArray<Maybe<Em3dCrystalEntity>>>,
  readonly em_3d_fitting?: Maybe<ReadonlyArray<Maybe<Em3dFitting>>>,
  readonly em_3d_fitting_list?: Maybe<ReadonlyArray<Maybe<Em3dFittingList>>>,
  readonly em_3d_reconstruction?: Maybe<ReadonlyArray<Maybe<Em3dReconstruction>>>,
  readonly em_ctf_correction?: Maybe<ReadonlyArray<Maybe<EmCtfCorrection>>>,
  readonly em_diffraction?: Maybe<ReadonlyArray<Maybe<EmDiffraction>>>,
  readonly em_diffraction_shell?: Maybe<ReadonlyArray<Maybe<EmDiffractionShell>>>,
  readonly em_diffraction_stats?: Maybe<ReadonlyArray<Maybe<EmDiffractionStats>>>,
  readonly em_embedding?: Maybe<ReadonlyArray<Maybe<EmEmbedding>>>,
  readonly em_entity_assembly?: Maybe<ReadonlyArray<Maybe<EmEntityAssembly>>>,
  readonly em_experiment?: Maybe<EmExperiment>,
  readonly em_helical_entity?: Maybe<ReadonlyArray<Maybe<EmHelicalEntity>>>,
  readonly em_image_recording?: Maybe<ReadonlyArray<Maybe<EmImageRecording>>>,
  readonly em_imaging?: Maybe<ReadonlyArray<Maybe<EmImaging>>>,
  readonly em_particle_selection?: Maybe<ReadonlyArray<Maybe<EmParticleSelection>>>,
  readonly em_single_particle_entity?: Maybe<ReadonlyArray<Maybe<EmSingleParticleEntity>>>,
  readonly em_software?: Maybe<ReadonlyArray<Maybe<EmSoftware>>>,
  readonly em_specimen?: Maybe<ReadonlyArray<Maybe<EmSpecimen>>>,
  readonly em_staining?: Maybe<ReadonlyArray<Maybe<EmStaining>>>,
  readonly em_vitrification?: Maybe<ReadonlyArray<Maybe<EmVitrification>>>,
  readonly entry?: Maybe<Entry>,
  readonly exptl?: Maybe<ReadonlyArray<Maybe<Exptl>>>,
  readonly exptl_crystal?: Maybe<ReadonlyArray<Maybe<ExptlCrystal>>>,
  readonly exptl_crystal_grow?: Maybe<ReadonlyArray<Maybe<ExptlCrystalGrow>>>,
  readonly nonpolymer_entities?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntity>>>,
  readonly pdbx_SG_project?: Maybe<ReadonlyArray<Maybe<PdbxSgProject>>>,
  readonly pdbx_audit_revision_category?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionCategory>>>,
  readonly pdbx_audit_revision_details?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionDetails>>>,
  readonly pdbx_audit_revision_group?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionGroup>>>,
  readonly pdbx_audit_revision_history?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionHistory>>>,
  readonly pdbx_audit_revision_item?: Maybe<ReadonlyArray<Maybe<PdbxAuditRevisionItem>>>,
  readonly pdbx_audit_support?: Maybe<ReadonlyArray<Maybe<PdbxAuditSupport>>>,
  readonly pdbx_database_PDB_obs_spr?: Maybe<ReadonlyArray<Maybe<PdbxDatabasePdbObsSpr>>>,
  readonly pdbx_database_related?: Maybe<ReadonlyArray<Maybe<PdbxDatabaseRelated>>>,
  readonly pdbx_database_status?: Maybe<PdbxDatabaseStatus>,
  readonly pdbx_deposit_group?: Maybe<ReadonlyArray<Maybe<PdbxDepositGroup>>>,
  readonly pdbx_molecule_features?: Maybe<ReadonlyArray<Maybe<PdbxMoleculeFeatures>>>,
  readonly pdbx_nmr_details?: Maybe<PdbxNmrDetails>,
  readonly pdbx_nmr_ensemble?: Maybe<PdbxNmrEnsemble>,
  readonly pdbx_nmr_exptl?: Maybe<ReadonlyArray<Maybe<PdbxNmrExptl>>>,
  readonly pdbx_nmr_exptl_sample_conditions?: Maybe<ReadonlyArray<Maybe<PdbxNmrExptlSampleConditions>>>,
  readonly pdbx_nmr_refine?: Maybe<ReadonlyArray<Maybe<PdbxNmrRefine>>>,
  readonly pdbx_nmr_representative?: Maybe<PdbxNmrRepresentative>,
  readonly pdbx_nmr_sample_details?: Maybe<ReadonlyArray<Maybe<PdbxNmrSampleDetails>>>,
  readonly pdbx_nmr_software?: Maybe<ReadonlyArray<Maybe<PdbxNmrSoftware>>>,
  readonly pdbx_nmr_spectrometer?: Maybe<ReadonlyArray<Maybe<PdbxNmrSpectrometer>>>,
  readonly pdbx_serial_crystallography_data_reduction?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographyDataReduction>>>,
  readonly pdbx_serial_crystallography_measurement?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographyMeasurement>>>,
  readonly pdbx_serial_crystallography_sample_delivery?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDelivery>>>,
  readonly pdbx_serial_crystallography_sample_delivery_fixed_target?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDeliveryFixedTarget>>>,
  readonly pdbx_serial_crystallography_sample_delivery_injection?: Maybe<ReadonlyArray<Maybe<PdbxSerialCrystallographySampleDeliveryInjection>>>,
  readonly pdbx_soln_scatter?: Maybe<ReadonlyArray<Maybe<PdbxSolnScatter>>>,
  readonly pdbx_soln_scatter_model?: Maybe<ReadonlyArray<Maybe<PdbxSolnScatterModel>>>,
  readonly pdbx_vrpt_summary?: Maybe<PdbxVrptSummary>,
  readonly polymer_entities?: Maybe<ReadonlyArray<Maybe<CorePolymerEntity>>>,
  readonly pubmed?: Maybe<CorePubmed>,
  readonly rcsb_accession_info?: Maybe<RcsbAccessionInfo>,
  readonly rcsb_associated_holdings?: Maybe<CurrentEntry>,
  readonly rcsb_binding_affinity?: Maybe<ReadonlyArray<Maybe<RcsbBindingAffinity>>>,
  readonly rcsb_entry_container_identifiers: RcsbEntryContainerIdentifiers,
  readonly rcsb_entry_info: RcsbEntryInfo,
  readonly rcsb_external_references?: Maybe<ReadonlyArray<Maybe<RcsbExternalReferences>>>,
  readonly rcsb_id: Scalars['String'],
  readonly refine?: Maybe<ReadonlyArray<Maybe<Refine>>>,
  readonly refine_analyze?: Maybe<ReadonlyArray<Maybe<RefineAnalyze>>>,
  readonly refine_hist?: Maybe<ReadonlyArray<Maybe<RefineHist>>>,
  readonly refine_ls_restr?: Maybe<ReadonlyArray<Maybe<RefineLsRestr>>>,
  readonly reflns?: Maybe<ReadonlyArray<Maybe<Reflns>>>,
  readonly reflns_shell?: Maybe<ReadonlyArray<Maybe<ReflnsShell>>>,
  readonly software?: Maybe<ReadonlyArray<Maybe<Software>>>,
  readonly struct?: Maybe<Struct>,
  readonly struct_keywords?: Maybe<StructKeywords>,
  readonly symmetry?: Maybe<Symmetry>,
};

export type CoreNonpolymerEntity = {
  readonly __typename?: 'CoreNonpolymerEntity',
  readonly nonpolymer_comp?: Maybe<CoreChemComp>,
  readonly nonpolymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntityInstance>>>,
  readonly pdbx_entity_nonpoly?: Maybe<PdbxEntityNonpoly>,
  readonly prd?: Maybe<CoreChemComp>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>,
  readonly rcsb_nonpolymer_entity?: Maybe<RcsbNonpolymerEntity>,
  readonly rcsb_nonpolymer_entity_annotation?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityAnnotation>>>,
  readonly rcsb_nonpolymer_entity_container_identifiers?: Maybe<RcsbNonpolymerEntityContainerIdentifiers>,
  readonly rcsb_nonpolymer_entity_feature?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityFeature>>>,
  readonly rcsb_nonpolymer_entity_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityFeatureSummary>>>,
  readonly rcsb_nonpolymer_entity_keywords?: Maybe<RcsbNonpolymerEntityKeywords>,
  readonly rcsb_nonpolymer_entity_name_com?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityNameCom>>>,
};

export type CoreNonpolymerEntityInstance = {
  readonly __typename?: 'CoreNonpolymerEntityInstance',
  readonly pdbx_struct_special_symmetry?: Maybe<ReadonlyArray<Maybe<PdbxStructSpecialSymmetry>>>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>,
  readonly rcsb_nonpolymer_entity_instance_container_identifiers?: Maybe<RcsbNonpolymerEntityInstanceContainerIdentifiers>,
  readonly rcsb_nonpolymer_instance_annotation?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceAnnotation>>>,
  readonly rcsb_nonpolymer_instance_feature?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeature>>>,
  readonly rcsb_nonpolymer_instance_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeatureSummary>>>,
  readonly rcsb_nonpolymer_struct_conn?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerStructConn>>>,
};

export type CorePfam = {
  readonly __typename?: 'CorePfam',
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_pfam_accession: Scalars['String'],
  readonly rcsb_pfam_clan_id?: Maybe<Scalars['String']>,
  readonly rcsb_pfam_comment?: Maybe<Scalars['String']>,
  readonly rcsb_pfam_container_identifiers: RcsbPfamContainerIdentifiers,
  readonly rcsb_pfam_description?: Maybe<Scalars['String']>,
  readonly rcsb_pfam_identifier?: Maybe<Scalars['String']>,
  readonly rcsb_pfam_provenance_code?: Maybe<Scalars['String']>,
  readonly rcsb_pfam_seed_source?: Maybe<Scalars['String']>,
};

export type CorePolymerEntity = {
  readonly __typename?: 'CorePolymerEntity',
  readonly chem_comp_monomers?: Maybe<ReadonlyArray<Maybe<CoreChemComp>>>,
  readonly chem_comp_nstd_monomers?: Maybe<ReadonlyArray<Maybe<CoreChemComp>>>,
  readonly entity_poly?: Maybe<EntityPoly>,
  readonly entity_src_gen?: Maybe<ReadonlyArray<Maybe<EntitySrcGen>>>,
  readonly entity_src_nat?: Maybe<ReadonlyArray<Maybe<EntitySrcNat>>>,
  readonly pdbx_entity_src_syn?: Maybe<ReadonlyArray<Maybe<PdbxEntitySrcSyn>>>,
  readonly pfams?: Maybe<ReadonlyArray<Maybe<CorePfam>>>,
  readonly polymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CorePolymerEntityInstance>>>,
  readonly prd?: Maybe<CoreChemComp>,
  readonly rcsb_cluster_flexibility?: Maybe<RcsbClusterFlexibility>,
  readonly rcsb_cluster_membership?: Maybe<ReadonlyArray<Maybe<RcsbClusterMembership>>>,
  readonly rcsb_entity_host_organism?: Maybe<ReadonlyArray<Maybe<RcsbEntityHostOrganism>>>,
  readonly rcsb_entity_source_organism?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganism>>>,
  readonly rcsb_genomic_lineage?: Maybe<ReadonlyArray<Maybe<RcsbGenomicLineage>>>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>,
  readonly rcsb_membrane_lineage?: Maybe<ReadonlyArray<Maybe<RcsbMembraneLineage>>>,
  readonly rcsb_membrane_lineage_provenance_code?: Maybe<Scalars['String']>,
  readonly rcsb_polymer_entity?: Maybe<RcsbPolymerEntity>,
  readonly rcsb_polymer_entity_align?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAlign>>>,
  readonly rcsb_polymer_entity_annotation?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAnnotation>>>,
  readonly rcsb_polymer_entity_container_identifiers: RcsbPolymerEntityContainerIdentifiers,
  readonly rcsb_polymer_entity_feature?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeature>>>,
  readonly rcsb_polymer_entity_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeatureSummary>>>,
  readonly rcsb_polymer_entity_keywords?: Maybe<RcsbPolymerEntityKeywords>,
  readonly rcsb_polymer_entity_name_com?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityNameCom>>>,
  readonly rcsb_polymer_entity_name_sys?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityNameSys>>>,
  readonly uniprots?: Maybe<ReadonlyArray<Maybe<CoreUniprot>>>,
};

export type CorePolymerEntityInstance = {
  readonly __typename?: 'CorePolymerEntityInstance',
  readonly pdbx_struct_special_symmetry?: Maybe<ReadonlyArray<Maybe<PdbxStructSpecialSymmetry>>>,
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_latest_revision?: Maybe<RcsbLatestRevision>,
  readonly rcsb_polymer_entity_instance_container_identifiers?: Maybe<RcsbPolymerEntityInstanceContainerIdentifiers>,
  readonly rcsb_polymer_instance_annotation?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceAnnotation>>>,
  readonly rcsb_polymer_instance_feature?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeature>>>,
  readonly rcsb_polymer_instance_feature_summary?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeatureSummary>>>,
  readonly rcsb_polymer_struct_conn?: Maybe<ReadonlyArray<Maybe<RcsbPolymerStructConn>>>,
};

export type CorePubmed = {
  readonly __typename?: 'CorePubmed',
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly rcsb_pubmed_abstract_text?: Maybe<Scalars['String']>,
  readonly rcsb_pubmed_affiliation_info?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_pubmed_central_id?: Maybe<Scalars['String']>,
  readonly rcsb_pubmed_container_identifiers: RcsbPubmedContainerIdentifiers,
  readonly rcsb_pubmed_doi?: Maybe<Scalars['String']>,
  readonly rcsb_pubmed_mesh_descriptors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_pubmed_mesh_descriptors_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPubmedMeshDescriptorsLineage>>>,
};

export type CoreUniprot = {
  readonly __typename?: 'CoreUniprot',
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly rcsb_uniprot_accession?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_uniprot_container_identifiers: RcsbUniprotContainerIdentifiers,
  readonly rcsb_uniprot_entry_name?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_uniprot_feature?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeature>>>,
  readonly rcsb_uniprot_keyword?: Maybe<ReadonlyArray<Maybe<RcsbUniprotKeyword>>>,
  readonly rcsb_uniprot_protein?: Maybe<RcsbUniprotProtein>,
};

export type CurrentEntry = {
  readonly __typename?: 'CurrentEntry',
  readonly rcsb_id: Scalars['String'],
  readonly rcsb_repository_holdings_current?: Maybe<RcsbRepositoryHoldingsCurrent>,
  readonly rcsb_repository_holdings_current_entry_container_identifiers?: Maybe<RcsbRepositoryHoldingsCurrentEntryContainerIdentifiers>,
};


export type Diffrn = {
  readonly __typename?: 'Diffrn',
  readonly ambient_pressure?: Maybe<Scalars['Float']>,
  readonly ambient_temp?: Maybe<Scalars['Float']>,
  readonly ambient_temp_details?: Maybe<Scalars['String']>,
  readonly crystal_id?: Maybe<Scalars['String']>,
  readonly crystal_support?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly pdbx_serial_crystal_experiment?: Maybe<Scalars['String']>,
};

export type DiffrnDetector = {
  readonly __typename?: 'DiffrnDetector',
  readonly details?: Maybe<Scalars['String']>,
  readonly detector?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly pdbx_collection_date?: Maybe<Scalars['Date']>,
  readonly pdbx_frequency?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type DiffrnRadiation = {
  readonly __typename?: 'DiffrnRadiation',
  readonly collimation?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly monochromator?: Maybe<Scalars['String']>,
  readonly pdbx_diffrn_protocol?: Maybe<Scalars['String']>,
  readonly pdbx_monochromatic_or_laue_m_l?: Maybe<Scalars['String']>,
  readonly pdbx_scattering_type?: Maybe<Scalars['String']>,
  readonly pdbx_wavelength?: Maybe<Scalars['String']>,
  readonly pdbx_wavelength_list?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
  readonly wavelength_id?: Maybe<Scalars['String']>,
};

export type DiffrnSource = {
  readonly __typename?: 'DiffrnSource',
  readonly details?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly pdbx_synchrotron_beamline?: Maybe<Scalars['String']>,
  readonly pdbx_synchrotron_site?: Maybe<Scalars['String']>,
  readonly pdbx_wavelength?: Maybe<Scalars['String']>,
  readonly pdbx_wavelength_list?: Maybe<Scalars['String']>,
  readonly source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type DrugbankContainerIdentifiers = {
  readonly __typename?: 'DrugbankContainerIdentifiers',
  readonly drugbank_id: Scalars['String'],
};

export type DrugbankInfo = {
  readonly __typename?: 'DrugbankInfo',
  readonly affected_organisms?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly atc_codes?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly brand_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly cas_number?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly drug_categories?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly drug_groups?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly drugbank_id: Scalars['String'],
  readonly indication?: Maybe<Scalars['String']>,
  readonly mechanism_of_action?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly pharmacology?: Maybe<Scalars['String']>,
  readonly synonyms?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type DrugbankTarget = {
  readonly __typename?: 'DrugbankTarget',
  readonly interaction_type?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly organism_common_name?: Maybe<Scalars['String']>,
  readonly reference_database_accession_code?: Maybe<Scalars['String']>,
  readonly reference_database_name?: Maybe<Scalars['String']>,
  readonly seq_one_letter_code?: Maybe<Scalars['String']>,
  readonly target_actions?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type Em2dCrystalEntity = {
  readonly __typename?: 'Em2dCrystalEntity',
  readonly angle_gamma?: Maybe<Scalars['Float']>,
  readonly c_sampling_length?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
  readonly length_a?: Maybe<Scalars['Float']>,
  readonly length_b?: Maybe<Scalars['Float']>,
  readonly length_c?: Maybe<Scalars['Float']>,
  readonly space_group_name_H_M?: Maybe<Scalars['String']>,
};

export type Em3dCrystalEntity = {
  readonly __typename?: 'Em3dCrystalEntity',
  readonly angle_alpha?: Maybe<Scalars['Float']>,
  readonly angle_beta?: Maybe<Scalars['Float']>,
  readonly angle_gamma?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
  readonly length_a?: Maybe<Scalars['Float']>,
  readonly length_b?: Maybe<Scalars['Float']>,
  readonly length_c?: Maybe<Scalars['Float']>,
  readonly space_group_name?: Maybe<Scalars['String']>,
  readonly space_group_num?: Maybe<Scalars['Int']>,
};

export type Em3dFitting = {
  readonly __typename?: 'Em3dFitting',
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly method?: Maybe<Scalars['String']>,
  readonly overall_b_value?: Maybe<Scalars['Float']>,
  readonly ref_protocol?: Maybe<Scalars['String']>,
  readonly ref_space?: Maybe<Scalars['String']>,
  readonly target_criteria?: Maybe<Scalars['String']>,
};

export type Em3dFittingList = {
  readonly __typename?: 'Em3dFittingList',
  readonly _3d_fitting_id: Scalars['String'],
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly pdb_chain_id?: Maybe<Scalars['String']>,
  readonly pdb_chain_residue_range?: Maybe<Scalars['String']>,
  readonly pdb_entry_id?: Maybe<Scalars['String']>,
};

export type Em3dReconstruction = {
  readonly __typename?: 'Em3dReconstruction',
  readonly actual_pixel_size?: Maybe<Scalars['Float']>,
  readonly algorithm?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
  readonly magnification_calibration?: Maybe<Scalars['String']>,
  readonly method?: Maybe<Scalars['String']>,
  readonly nominal_pixel_size?: Maybe<Scalars['Float']>,
  readonly num_class_averages?: Maybe<Scalars['Int']>,
  readonly num_particles?: Maybe<Scalars['Int']>,
  readonly refinement_type?: Maybe<Scalars['String']>,
  readonly resolution?: Maybe<Scalars['Float']>,
  readonly resolution_method?: Maybe<Scalars['String']>,
  readonly symmetry_type?: Maybe<Scalars['String']>,
};

export type EmCtfCorrection = {
  readonly __typename?: 'EmCtfCorrection',
  readonly details?: Maybe<Scalars['String']>,
  readonly em_image_processing_id?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type EmDiffraction = {
  readonly __typename?: 'EmDiffraction',
  readonly camera_length?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly imaging_id?: Maybe<Scalars['String']>,
  readonly tilt_angle_list?: Maybe<Scalars['String']>,
};

export type EmDiffractionShell = {
  readonly __typename?: 'EmDiffractionShell',
  readonly em_diffraction_stats_id?: Maybe<Scalars['String']>,
  readonly fourier_space_coverage?: Maybe<Scalars['Float']>,
  readonly high_resolution?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly low_resolution?: Maybe<Scalars['Float']>,
  readonly multiplicity?: Maybe<Scalars['Float']>,
  readonly num_structure_factors?: Maybe<Scalars['Int']>,
  readonly phase_residual?: Maybe<Scalars['Float']>,
};

export type EmDiffractionStats = {
  readonly __typename?: 'EmDiffractionStats',
  readonly details?: Maybe<Scalars['String']>,
  readonly fourier_space_coverage?: Maybe<Scalars['Float']>,
  readonly high_resolution?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly image_processing_id?: Maybe<Scalars['String']>,
  readonly num_intensities_measured?: Maybe<Scalars['Int']>,
  readonly num_structure_factors?: Maybe<Scalars['Int']>,
  readonly overall_phase_error?: Maybe<Scalars['Float']>,
  readonly overall_phase_residual?: Maybe<Scalars['Float']>,
  readonly phase_error_rejection_criteria?: Maybe<Scalars['String']>,
  readonly r_merge?: Maybe<Scalars['Float']>,
  readonly r_sym?: Maybe<Scalars['Float']>,
};

export type EmEmbedding = {
  readonly __typename?: 'EmEmbedding',
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly material?: Maybe<Scalars['String']>,
  readonly specimen_id?: Maybe<Scalars['String']>,
};

export type EmEntityAssembly = {
  readonly __typename?: 'EmEntityAssembly',
  readonly details?: Maybe<Scalars['String']>,
  readonly entity_id_list?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly id: Scalars['String'],
  readonly name?: Maybe<Scalars['String']>,
  readonly oligomeric_details?: Maybe<Scalars['String']>,
  readonly parent_id?: Maybe<Scalars['Int']>,
  readonly source?: Maybe<Scalars['String']>,
  readonly synonym?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type EmExperiment = {
  readonly __typename?: 'EmExperiment',
  readonly aggregation_state?: Maybe<Scalars['String']>,
  readonly entity_assembly_id?: Maybe<Scalars['String']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly reconstruction_method?: Maybe<Scalars['String']>,
};

export type EmHelicalEntity = {
  readonly __typename?: 'EmHelicalEntity',
  readonly angular_rotation_per_subunit?: Maybe<Scalars['Float']>,
  readonly axial_rise_per_subunit?: Maybe<Scalars['Float']>,
  readonly axial_symmetry?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
};

export type EmImageRecording = {
  readonly __typename?: 'EmImageRecording',
  readonly average_exposure_time?: Maybe<Scalars['Float']>,
  readonly avg_electron_dose_per_image?: Maybe<Scalars['Float']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly detector_mode?: Maybe<Scalars['String']>,
  readonly film_or_detector_model?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly imaging_id: Scalars['String'],
  readonly num_diffraction_images?: Maybe<Scalars['Int']>,
  readonly num_grids_imaged?: Maybe<Scalars['Int']>,
  readonly num_real_images?: Maybe<Scalars['Int']>,
};

export type EmImaging = {
  readonly __typename?: 'EmImaging',
  readonly accelerating_voltage?: Maybe<Scalars['Int']>,
  readonly alignment_procedure?: Maybe<Scalars['String']>,
  readonly astigmatism?: Maybe<Scalars['String']>,
  readonly c2_aperture_diameter?: Maybe<Scalars['Float']>,
  readonly calibrated_defocus_max?: Maybe<Scalars['Float']>,
  readonly calibrated_defocus_min?: Maybe<Scalars['Float']>,
  readonly calibrated_magnification?: Maybe<Scalars['Int']>,
  readonly cryogen?: Maybe<Scalars['String']>,
  readonly date?: Maybe<Scalars['Date']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly detector_distance?: Maybe<Scalars['Float']>,
  readonly electron_beam_tilt_params?: Maybe<Scalars['String']>,
  readonly electron_source?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly illumination_mode?: Maybe<Scalars['String']>,
  readonly microscope_model?: Maybe<Scalars['String']>,
  readonly mode?: Maybe<Scalars['String']>,
  readonly nominal_cs?: Maybe<Scalars['Float']>,
  readonly nominal_defocus_max?: Maybe<Scalars['Float']>,
  readonly nominal_defocus_min?: Maybe<Scalars['Float']>,
  readonly nominal_magnification?: Maybe<Scalars['Int']>,
  readonly recording_temperature_maximum?: Maybe<Scalars['Float']>,
  readonly recording_temperature_minimum?: Maybe<Scalars['Float']>,
  readonly residual_tilt?: Maybe<Scalars['Float']>,
  readonly specimen_holder_model?: Maybe<Scalars['String']>,
  readonly specimen_holder_type?: Maybe<Scalars['String']>,
  readonly specimen_id?: Maybe<Scalars['String']>,
  readonly temperature?: Maybe<Scalars['Float']>,
  readonly tilt_angle_max?: Maybe<Scalars['Float']>,
  readonly tilt_angle_min?: Maybe<Scalars['Float']>,
};

export type EmParticleSelection = {
  readonly __typename?: 'EmParticleSelection',
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
  readonly num_particles_selected?: Maybe<Scalars['Int']>,
};

export type EmSingleParticleEntity = {
  readonly __typename?: 'EmSingleParticleEntity',
  readonly id: Scalars['String'],
  readonly image_processing_id: Scalars['String'],
  readonly point_symmetry?: Maybe<Scalars['String']>,
};

export type EmSoftware = {
  readonly __typename?: 'EmSoftware',
  readonly category?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly fitting_id?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly image_processing_id?: Maybe<Scalars['String']>,
  readonly imaging_id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly version?: Maybe<Scalars['String']>,
};

export type EmSpecimen = {
  readonly __typename?: 'EmSpecimen',
  readonly concentration?: Maybe<Scalars['Float']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly embedding_applied?: Maybe<Scalars['String']>,
  readonly experiment_id: Scalars['String'],
  readonly id: Scalars['String'],
  readonly shadowing_applied?: Maybe<Scalars['String']>,
  readonly staining_applied?: Maybe<Scalars['String']>,
  readonly vitrification_applied?: Maybe<Scalars['String']>,
};

export type EmStaining = {
  readonly __typename?: 'EmStaining',
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly material?: Maybe<Scalars['String']>,
  readonly specimen_id?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type EmVitrification = {
  readonly __typename?: 'EmVitrification',
  readonly chamber_temperature?: Maybe<Scalars['Float']>,
  readonly cryogen_name?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly humidity?: Maybe<Scalars['Float']>,
  readonly id: Scalars['String'],
  readonly instrument?: Maybe<Scalars['String']>,
  readonly method?: Maybe<Scalars['String']>,
  readonly specimen_id: Scalars['String'],
  readonly temp?: Maybe<Scalars['Float']>,
  readonly time_resolved_state?: Maybe<Scalars['String']>,
};

export type EntityPoly = {
  readonly __typename?: 'EntityPoly',
  readonly nstd_linkage?: Maybe<Scalars['String']>,
  readonly nstd_monomer?: Maybe<Scalars['String']>,
  readonly pdbx_seq_one_letter_code?: Maybe<Scalars['String']>,
  readonly pdbx_seq_one_letter_code_can?: Maybe<Scalars['String']>,
  readonly pdbx_strand_id?: Maybe<Scalars['String']>,
  readonly pdbx_target_identifier?: Maybe<Scalars['String']>,
  readonly rcsb_artifact_monomer_count?: Maybe<Scalars['Int']>,
  readonly rcsb_conflict_count?: Maybe<Scalars['Int']>,
  readonly rcsb_deletion_count?: Maybe<Scalars['Int']>,
  readonly rcsb_entity_polymer_type?: Maybe<Scalars['String']>,
  readonly rcsb_insertion_count?: Maybe<Scalars['Int']>,
  readonly rcsb_mutation_count?: Maybe<Scalars['Int']>,
  readonly rcsb_non_std_monomer_count?: Maybe<Scalars['Int']>,
  readonly rcsb_non_std_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly rcsb_prd_id?: Maybe<Scalars['String']>,
  readonly rcsb_sample_sequence_length?: Maybe<Scalars['Int']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type EntitySrcGen = {
  readonly __typename?: 'EntitySrcGen',
  readonly expression_system_id?: Maybe<Scalars['String']>,
  readonly gene_src_common_name?: Maybe<Scalars['String']>,
  readonly gene_src_details?: Maybe<Scalars['String']>,
  readonly gene_src_genus?: Maybe<Scalars['String']>,
  readonly gene_src_species?: Maybe<Scalars['String']>,
  readonly gene_src_strain?: Maybe<Scalars['String']>,
  readonly gene_src_tissue?: Maybe<Scalars['String']>,
  readonly gene_src_tissue_fraction?: Maybe<Scalars['String']>,
  readonly host_org_common_name?: Maybe<Scalars['String']>,
  readonly host_org_details?: Maybe<Scalars['String']>,
  readonly host_org_genus?: Maybe<Scalars['String']>,
  readonly host_org_species?: Maybe<Scalars['String']>,
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>,
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_description?: Maybe<Scalars['String']>,
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_gene_src_atcc?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_cell?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_cell_line?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_cellular_location?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_fragment?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_gene?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_ncbi_taxonomy_id?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_organ?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_organelle?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_scientific_name?: Maybe<Scalars['String']>,
  readonly pdbx_gene_src_variant?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_atcc?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_cell?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_cell_line?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_cellular_location?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_culture_collection?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_gene?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_ncbi_taxonomy_id?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_organ?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_organelle?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_scientific_name?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_strain?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_tissue?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_tissue_fraction?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_variant?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_vector?: Maybe<Scalars['String']>,
  readonly pdbx_host_org_vector_type?: Maybe<Scalars['String']>,
  readonly pdbx_seq_type?: Maybe<Scalars['String']>,
  readonly pdbx_src_id: Scalars['Int'],
  readonly plasmid_details?: Maybe<Scalars['String']>,
  readonly plasmid_name?: Maybe<Scalars['String']>,
};

export type EntitySrcNat = {
  readonly __typename?: 'EntitySrcNat',
  readonly common_name?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly genus?: Maybe<Scalars['String']>,
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>,
  readonly pdbx_atcc?: Maybe<Scalars['String']>,
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_cell?: Maybe<Scalars['String']>,
  readonly pdbx_cell_line?: Maybe<Scalars['String']>,
  readonly pdbx_cellular_location?: Maybe<Scalars['String']>,
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_fragment?: Maybe<Scalars['String']>,
  readonly pdbx_ncbi_taxonomy_id?: Maybe<Scalars['String']>,
  readonly pdbx_organ?: Maybe<Scalars['String']>,
  readonly pdbx_organelle?: Maybe<Scalars['String']>,
  readonly pdbx_organism_scientific?: Maybe<Scalars['String']>,
  readonly pdbx_plasmid_details?: Maybe<Scalars['String']>,
  readonly pdbx_plasmid_name?: Maybe<Scalars['String']>,
  readonly pdbx_secretion?: Maybe<Scalars['String']>,
  readonly pdbx_src_id: Scalars['Int'],
  readonly pdbx_variant?: Maybe<Scalars['String']>,
  readonly species?: Maybe<Scalars['String']>,
  readonly strain?: Maybe<Scalars['String']>,
  readonly tissue?: Maybe<Scalars['String']>,
  readonly tissue_fraction?: Maybe<Scalars['String']>,
};

export type Entry = {
  readonly __typename?: 'Entry',
  readonly id: Scalars['String'],
};

export type Exptl = {
  readonly __typename?: 'Exptl',
  readonly crystals_number?: Maybe<Scalars['Int']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly method: Scalars['String'],
  readonly method_details?: Maybe<Scalars['String']>,
};

export type ExptlCrystal = {
  readonly __typename?: 'ExptlCrystal',
  readonly colour?: Maybe<Scalars['String']>,
  readonly density_Matthews?: Maybe<Scalars['Float']>,
  readonly density_meas?: Maybe<Scalars['Float']>,
  readonly density_percent_sol?: Maybe<Scalars['Float']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly pdbx_mosaicity?: Maybe<Scalars['Float']>,
  readonly pdbx_mosaicity_esd?: Maybe<Scalars['Float']>,
  readonly preparation?: Maybe<Scalars['String']>,
};

export type ExptlCrystalGrow = {
  readonly __typename?: 'ExptlCrystalGrow',
  readonly crystal_id: Scalars['String'],
  readonly details?: Maybe<Scalars['String']>,
  readonly method?: Maybe<Scalars['String']>,
  readonly pH?: Maybe<Scalars['Float']>,
  readonly pdbx_details?: Maybe<Scalars['String']>,
  readonly pdbx_pH_range?: Maybe<Scalars['String']>,
  readonly temp?: Maybe<Scalars['Float']>,
  readonly temp_details?: Maybe<Scalars['String']>,
};

export type GeneName = {
  readonly __typename?: 'GeneName',
  readonly type?: Maybe<Scalars['String']>,
  readonly value?: Maybe<Scalars['String']>,
};

export type PdbxAuditRevisionCategory = {
  readonly __typename?: 'PdbxAuditRevisionCategory',
  readonly category?: Maybe<Scalars['String']>,
  readonly data_content_type: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly revision_ordinal: Scalars['Int'],
};

export type PdbxAuditRevisionDetails = {
  readonly __typename?: 'PdbxAuditRevisionDetails',
  readonly data_content_type: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly provider?: Maybe<Scalars['String']>,
  readonly revision_ordinal: Scalars['Int'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxAuditRevisionGroup = {
  readonly __typename?: 'PdbxAuditRevisionGroup',
  readonly data_content_type: Scalars['String'],
  readonly group?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly revision_ordinal: Scalars['Int'],
};

export type PdbxAuditRevisionHistory = {
  readonly __typename?: 'PdbxAuditRevisionHistory',
  readonly data_content_type: Scalars['String'],
  readonly major_revision?: Maybe<Scalars['Int']>,
  readonly minor_revision?: Maybe<Scalars['Int']>,
  readonly ordinal: Scalars['Int'],
  readonly revision_date?: Maybe<Scalars['Date']>,
};

export type PdbxAuditRevisionItem = {
  readonly __typename?: 'PdbxAuditRevisionItem',
  readonly data_content_type: Scalars['String'],
  readonly item?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly revision_ordinal: Scalars['Int'],
};

export type PdbxAuditSupport = {
  readonly __typename?: 'PdbxAuditSupport',
  readonly country?: Maybe<Scalars['String']>,
  readonly funding_organization?: Maybe<Scalars['String']>,
  readonly grant_number?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
};

export type PdbxChemCompAudit = {
  readonly __typename?: 'PdbxChemCompAudit',
  readonly action_type?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly date?: Maybe<Scalars['Date']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
};

export type PdbxChemCompDescriptor = {
  readonly __typename?: 'PdbxChemCompDescriptor',
  readonly comp_id: Scalars['String'],
  readonly descriptor?: Maybe<Scalars['String']>,
  readonly program: Scalars['String'],
  readonly program_version: Scalars['String'],
  readonly type: Scalars['String'],
};

export type PdbxChemCompFeature = {
  readonly __typename?: 'PdbxChemCompFeature',
  readonly comp_id: Scalars['String'],
  readonly source: Scalars['String'],
  readonly type: Scalars['String'],
  readonly value: Scalars['String'],
};

export type PdbxChemCompIdentifier = {
  readonly __typename?: 'PdbxChemCompIdentifier',
  readonly comp_id: Scalars['String'],
  readonly identifier?: Maybe<Scalars['String']>,
  readonly program: Scalars['String'],
  readonly program_version: Scalars['String'],
  readonly type: Scalars['String'],
};

export type PdbxDatabasePdbObsSpr = {
  readonly __typename?: 'PdbxDatabasePDBObsSpr',
  readonly date?: Maybe<Scalars['Date']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly pdb_id: Scalars['String'],
  readonly replace_pdb_id: Scalars['String'],
};

export type PdbxDatabaseRelated = {
  readonly __typename?: 'PdbxDatabaseRelated',
  readonly content_type: Scalars['String'],
  readonly db_id: Scalars['String'],
  readonly db_name: Scalars['String'],
  readonly details?: Maybe<Scalars['String']>,
};

export type PdbxDatabaseStatus = {
  readonly __typename?: 'PdbxDatabaseStatus',
  readonly SG_entry?: Maybe<Scalars['String']>,
  readonly deposit_site?: Maybe<Scalars['String']>,
  readonly methods_development_category?: Maybe<Scalars['String']>,
  readonly pdb_format_compatible?: Maybe<Scalars['String']>,
  readonly process_site?: Maybe<Scalars['String']>,
  readonly recvd_initial_deposition_date?: Maybe<Scalars['Date']>,
  readonly status_code?: Maybe<Scalars['String']>,
  readonly status_code_cs?: Maybe<Scalars['String']>,
  readonly status_code_mr?: Maybe<Scalars['String']>,
  readonly status_code_sf?: Maybe<Scalars['String']>,
};

export type PdbxDepositGroup = {
  readonly __typename?: 'PdbxDepositGroup',
  readonly group_description?: Maybe<Scalars['String']>,
  readonly group_id: Scalars['String'],
  readonly group_title?: Maybe<Scalars['String']>,
  readonly group_type?: Maybe<Scalars['String']>,
};

export type PdbxEntityNonpoly = {
  readonly __typename?: 'PdbxEntityNonpoly',
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly entity_id: Scalars['String'],
  readonly name?: Maybe<Scalars['String']>,
  readonly rcsb_prd_id?: Maybe<Scalars['String']>,
};

export type PdbxEntitySrcSyn = {
  readonly __typename?: 'PdbxEntitySrcSyn',
  readonly details?: Maybe<Scalars['String']>,
  readonly ncbi_taxonomy_id?: Maybe<Scalars['String']>,
  readonly organism_common_name?: Maybe<Scalars['String']>,
  readonly organism_scientific?: Maybe<Scalars['String']>,
  readonly pdbx_alt_source_flag?: Maybe<Scalars['String']>,
  readonly pdbx_beg_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_end_seq_num?: Maybe<Scalars['Int']>,
  readonly pdbx_src_id: Scalars['Int'],
};

export type PdbxFamilyPrdAudit = {
  readonly __typename?: 'PdbxFamilyPrdAudit',
  readonly action_type: Scalars['String'],
  readonly annotator?: Maybe<Scalars['String']>,
  readonly date: Scalars['Date'],
  readonly details?: Maybe<Scalars['String']>,
  readonly family_prd_id: Scalars['String'],
  readonly processing_site?: Maybe<Scalars['String']>,
};

export type PdbxMoleculeFeatures = {
  readonly __typename?: 'PdbxMoleculeFeatures',
  readonly class?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxNmrDetails = {
  readonly __typename?: 'PdbxNmrDetails',
  readonly text?: Maybe<Scalars['String']>,
};

export type PdbxNmrEnsemble = {
  readonly __typename?: 'PdbxNmrEnsemble',
  readonly average_constraint_violations_per_residue?: Maybe<Scalars['Int']>,
  readonly average_constraints_per_residue?: Maybe<Scalars['Int']>,
  readonly average_distance_constraint_violation?: Maybe<Scalars['Float']>,
  readonly average_torsion_angle_constraint_violation?: Maybe<Scalars['Float']>,
  readonly conformer_selection_criteria?: Maybe<Scalars['String']>,
  readonly conformers_calculated_total_number?: Maybe<Scalars['Int']>,
  readonly conformers_submitted_total_number?: Maybe<Scalars['Int']>,
  readonly distance_constraint_violation_method?: Maybe<Scalars['String']>,
  readonly maximum_distance_constraint_violation?: Maybe<Scalars['Float']>,
  readonly maximum_lower_distance_constraint_violation?: Maybe<Scalars['Float']>,
  readonly maximum_torsion_angle_constraint_violation?: Maybe<Scalars['Float']>,
  readonly maximum_upper_distance_constraint_violation?: Maybe<Scalars['Float']>,
  readonly representative_conformer?: Maybe<Scalars['Int']>,
  readonly torsion_angle_constraint_violation_method?: Maybe<Scalars['String']>,
};

export type PdbxNmrExptl = {
  readonly __typename?: 'PdbxNmrExptl',
  readonly conditions_id: Scalars['String'],
  readonly experiment_id: Scalars['String'],
  readonly sample_state?: Maybe<Scalars['String']>,
  readonly solution_id: Scalars['String'],
  readonly spectrometer_id?: Maybe<Scalars['Int']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxNmrExptlSampleConditions = {
  readonly __typename?: 'PdbxNmrExptlSampleConditions',
  readonly conditions_id: Scalars['String'],
  readonly details?: Maybe<Scalars['String']>,
  readonly ionic_strength?: Maybe<Scalars['String']>,
  readonly ionic_strength_err?: Maybe<Scalars['Float']>,
  readonly ionic_strength_units?: Maybe<Scalars['String']>,
  readonly label?: Maybe<Scalars['String']>,
  readonly pH?: Maybe<Scalars['String']>,
  readonly pH_err?: Maybe<Scalars['Float']>,
  readonly pH_units?: Maybe<Scalars['String']>,
  readonly pressure?: Maybe<Scalars['String']>,
  readonly pressure_err?: Maybe<Scalars['Float']>,
  readonly pressure_units?: Maybe<Scalars['String']>,
  readonly temperature?: Maybe<Scalars['String']>,
  readonly temperature_err?: Maybe<Scalars['Float']>,
  readonly temperature_units?: Maybe<Scalars['String']>,
};

export type PdbxNmrRefine = {
  readonly __typename?: 'PdbxNmrRefine',
  readonly details?: Maybe<Scalars['String']>,
  readonly method?: Maybe<Scalars['String']>,
  readonly software_ordinal: Scalars['Int'],
};

export type PdbxNmrRepresentative = {
  readonly __typename?: 'PdbxNmrRepresentative',
  readonly conformer_id?: Maybe<Scalars['String']>,
  readonly selection_criteria?: Maybe<Scalars['String']>,
};

export type PdbxNmrSampleDetails = {
  readonly __typename?: 'PdbxNmrSampleDetails',
  readonly contents?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly label?: Maybe<Scalars['String']>,
  readonly solution_id: Scalars['String'],
  readonly solvent_system?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxNmrSoftware = {
  readonly __typename?: 'PdbxNmrSoftware',
  readonly authors?: Maybe<Scalars['String']>,
  readonly classification?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly version?: Maybe<Scalars['String']>,
};

export type PdbxNmrSpectrometer = {
  readonly __typename?: 'PdbxNmrSpectrometer',
  readonly details?: Maybe<Scalars['String']>,
  readonly field_strength?: Maybe<Scalars['Float']>,
  readonly manufacturer?: Maybe<Scalars['String']>,
  readonly model?: Maybe<Scalars['String']>,
  readonly spectrometer_id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxPrdAudit = {
  readonly __typename?: 'PdbxPrdAudit',
  readonly action_type: Scalars['String'],
  readonly annotator?: Maybe<Scalars['String']>,
  readonly date: Scalars['Date'],
  readonly details?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly processing_site?: Maybe<Scalars['String']>,
};

export type PdbxReferenceEntityList = {
  readonly __typename?: 'PdbxReferenceEntityList',
  readonly component_id: Scalars['Int'],
  readonly details?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxReferenceEntityPoly = {
  readonly __typename?: 'PdbxReferenceEntityPoly',
  readonly db_code?: Maybe<Scalars['String']>,
  readonly db_name?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxReferenceEntityPolyLink = {
  readonly __typename?: 'PdbxReferenceEntityPolyLink',
  readonly atom_id_1?: Maybe<Scalars['String']>,
  readonly atom_id_2?: Maybe<Scalars['String']>,
  readonly comp_id_1?: Maybe<Scalars['String']>,
  readonly comp_id_2?: Maybe<Scalars['String']>,
  readonly component_id: Scalars['Int'],
  readonly entity_seq_num_1?: Maybe<Scalars['Int']>,
  readonly entity_seq_num_2?: Maybe<Scalars['Int']>,
  readonly link_id: Scalars['Int'],
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
  readonly value_order?: Maybe<Scalars['String']>,
};

export type PdbxReferenceEntityPolySeq = {
  readonly __typename?: 'PdbxReferenceEntityPolySeq',
  readonly hetero: Scalars['String'],
  readonly mon_id: Scalars['String'],
  readonly num: Scalars['Int'],
  readonly observed?: Maybe<Scalars['String']>,
  readonly parent_mon_id?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
};

export type PdbxReferenceEntitySequence = {
  readonly __typename?: 'PdbxReferenceEntitySequence',
  readonly NRP_flag?: Maybe<Scalars['String']>,
  readonly one_letter_codes?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxReferenceEntitySrcNat = {
  readonly __typename?: 'PdbxReferenceEntitySrcNat',
  readonly atcc?: Maybe<Scalars['String']>,
  readonly db_code?: Maybe<Scalars['String']>,
  readonly db_name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly organism_scientific?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly ref_entity_id: Scalars['String'],
  readonly source?: Maybe<Scalars['String']>,
  readonly source_id?: Maybe<Scalars['String']>,
  readonly taxid?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMolecule = {
  readonly __typename?: 'PdbxReferenceMolecule',
  readonly chem_comp_id?: Maybe<Scalars['String']>,
  readonly class?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly class_evidence_code?: Maybe<Scalars['String']>,
  readonly compound_details?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly formula?: Maybe<Scalars['String']>,
  readonly formula_weight?: Maybe<Scalars['Float']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly prd_id: Scalars['String'],
  readonly release_status?: Maybe<Scalars['String']>,
  readonly replaced_by?: Maybe<Scalars['String']>,
  readonly replaces?: Maybe<Scalars['String']>,
  readonly represent_as?: Maybe<Scalars['String']>,
  readonly representative_PDB_id_code?: Maybe<Scalars['String']>,
  readonly type?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly type_evidence_code?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMoleculeAnnotation = {
  readonly __typename?: 'PdbxReferenceMoleculeAnnotation',
  readonly family_prd_id: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly prd_id?: Maybe<Scalars['String']>,
  readonly source?: Maybe<Scalars['String']>,
  readonly text?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMoleculeDetails = {
  readonly __typename?: 'PdbxReferenceMoleculeDetails',
  readonly family_prd_id: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly source?: Maybe<Scalars['String']>,
  readonly source_id?: Maybe<Scalars['String']>,
  readonly text?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMoleculeFamily = {
  readonly __typename?: 'PdbxReferenceMoleculeFamily',
  readonly family_prd_id: Scalars['String'],
  readonly name?: Maybe<Scalars['String']>,
  readonly release_status?: Maybe<Scalars['String']>,
  readonly replaced_by?: Maybe<Scalars['String']>,
  readonly replaces?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMoleculeFeatures = {
  readonly __typename?: 'PdbxReferenceMoleculeFeatures',
  readonly family_prd_id: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly prd_id: Scalars['String'],
  readonly source?: Maybe<Scalars['String']>,
  readonly source_ordinal?: Maybe<Scalars['Int']>,
  readonly type?: Maybe<Scalars['String']>,
  readonly value?: Maybe<Scalars['String']>,
};

export type PdbxReferenceMoleculeList = {
  readonly __typename?: 'PdbxReferenceMoleculeList',
  readonly family_prd_id: Scalars['String'],
  readonly prd_id: Scalars['String'],
};

export type PdbxReferenceMoleculeRelatedStructures = {
  readonly __typename?: 'PdbxReferenceMoleculeRelatedStructures',
  readonly citation_id?: Maybe<Scalars['String']>,
  readonly db_accession?: Maybe<Scalars['String']>,
  readonly db_code?: Maybe<Scalars['String']>,
  readonly db_name?: Maybe<Scalars['String']>,
  readonly family_prd_id: Scalars['String'],
  readonly formula?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
};

export type PdbxReferenceMoleculeSynonyms = {
  readonly __typename?: 'PdbxReferenceMoleculeSynonyms',
  readonly family_prd_id: Scalars['String'],
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly prd_id: Scalars['String'],
  readonly source?: Maybe<Scalars['String']>,
};

export type PdbxSerialCrystallographyDataReduction = {
  readonly __typename?: 'PdbxSerialCrystallographyDataReduction',
  readonly crystal_hits?: Maybe<Scalars['Int']>,
  readonly diffrn_id: Scalars['String'],
  readonly droplet_hits?: Maybe<Scalars['Int']>,
  readonly frame_hits?: Maybe<Scalars['Int']>,
  readonly frames_failed_index?: Maybe<Scalars['Int']>,
  readonly frames_indexed?: Maybe<Scalars['Int']>,
  readonly frames_total?: Maybe<Scalars['Int']>,
  readonly lattices_indexed?: Maybe<Scalars['Int']>,
  readonly xfel_pulse_events?: Maybe<Scalars['Int']>,
  readonly xfel_run_numbers?: Maybe<Scalars['String']>,
};

export type PdbxSerialCrystallographyMeasurement = {
  readonly __typename?: 'PdbxSerialCrystallographyMeasurement',
  readonly collection_time_total?: Maybe<Scalars['Float']>,
  readonly collimation?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly focal_spot_size?: Maybe<Scalars['Float']>,
  readonly photons_per_pulse?: Maybe<Scalars['Float']>,
  readonly pulse_duration?: Maybe<Scalars['Float']>,
  readonly pulse_energy?: Maybe<Scalars['Float']>,
  readonly pulse_photon_energy?: Maybe<Scalars['Float']>,
  readonly source_distance?: Maybe<Scalars['Float']>,
  readonly source_size?: Maybe<Scalars['Float']>,
  readonly xfel_pulse_repetition_rate?: Maybe<Scalars['Float']>,
};

export type PdbxSerialCrystallographySampleDelivery = {
  readonly __typename?: 'PdbxSerialCrystallographySampleDelivery',
  readonly description?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly method?: Maybe<Scalars['String']>,
};

export type PdbxSerialCrystallographySampleDeliveryFixedTarget = {
  readonly __typename?: 'PdbxSerialCrystallographySampleDeliveryFixedTarget',
  readonly crystals_per_unit?: Maybe<Scalars['Int']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly motion_control?: Maybe<Scalars['String']>,
  readonly sample_dehydration_prevention?: Maybe<Scalars['String']>,
  readonly sample_holding?: Maybe<Scalars['String']>,
  readonly sample_solvent?: Maybe<Scalars['String']>,
  readonly sample_unit_size?: Maybe<Scalars['Float']>,
  readonly support_base?: Maybe<Scalars['String']>,
  readonly velocity_horizontal?: Maybe<Scalars['Float']>,
  readonly velocity_vertical?: Maybe<Scalars['Float']>,
};

export type PdbxSerialCrystallographySampleDeliveryInjection = {
  readonly __typename?: 'PdbxSerialCrystallographySampleDeliveryInjection',
  readonly carrier_solvent?: Maybe<Scalars['String']>,
  readonly crystal_concentration?: Maybe<Scalars['Float']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly diffrn_id: Scalars['String'],
  readonly filter_size?: Maybe<Scalars['Float']>,
  readonly flow_rate?: Maybe<Scalars['Float']>,
  readonly injector_diameter?: Maybe<Scalars['Float']>,
  readonly injector_nozzle?: Maybe<Scalars['String']>,
  readonly injector_pressure?: Maybe<Scalars['Float']>,
  readonly injector_temperature?: Maybe<Scalars['Float']>,
  readonly jet_diameter?: Maybe<Scalars['Float']>,
  readonly power_by?: Maybe<Scalars['String']>,
  readonly preparation?: Maybe<Scalars['String']>,
};

export type PdbxSgProject = {
  readonly __typename?: 'PdbxSGProject',
  readonly full_name_of_center?: Maybe<Scalars['String']>,
  readonly id: Scalars['Int'],
  readonly initial_of_center?: Maybe<Scalars['String']>,
  readonly project_name?: Maybe<Scalars['String']>,
};

export type PdbxSolnScatter = {
  readonly __typename?: 'PdbxSolnScatter',
  readonly buffer_name?: Maybe<Scalars['String']>,
  readonly concentration_range?: Maybe<Scalars['String']>,
  readonly data_analysis_software_list?: Maybe<Scalars['String']>,
  readonly data_reduction_software_list?: Maybe<Scalars['String']>,
  readonly detector_specific?: Maybe<Scalars['String']>,
  readonly detector_type?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly max_mean_cross_sectional_radii_gyration?: Maybe<Scalars['Float']>,
  readonly max_mean_cross_sectional_radii_gyration_esd?: Maybe<Scalars['Float']>,
  readonly mean_guiner_radius?: Maybe<Scalars['Float']>,
  readonly mean_guiner_radius_esd?: Maybe<Scalars['Float']>,
  readonly min_mean_cross_sectional_radii_gyration?: Maybe<Scalars['Float']>,
  readonly min_mean_cross_sectional_radii_gyration_esd?: Maybe<Scalars['Float']>,
  readonly num_time_frames?: Maybe<Scalars['Int']>,
  readonly protein_length?: Maybe<Scalars['String']>,
  readonly sample_pH?: Maybe<Scalars['Float']>,
  readonly source_beamline?: Maybe<Scalars['String']>,
  readonly source_beamline_instrument?: Maybe<Scalars['String']>,
  readonly source_class?: Maybe<Scalars['String']>,
  readonly source_type?: Maybe<Scalars['String']>,
  readonly temperature?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type PdbxSolnScatterModel = {
  readonly __typename?: 'PdbxSolnScatterModel',
  readonly conformer_selection_criteria?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly entry_fitting_list?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly method?: Maybe<Scalars['String']>,
  readonly num_conformers_calculated?: Maybe<Scalars['Int']>,
  readonly num_conformers_submitted?: Maybe<Scalars['Int']>,
  readonly representative_conformer?: Maybe<Scalars['Int']>,
  readonly scatter_id: Scalars['String'],
  readonly software_author_list?: Maybe<Scalars['String']>,
  readonly software_list?: Maybe<Scalars['String']>,
};

export type PdbxStructAssembly = {
  readonly __typename?: 'PdbxStructAssembly',
  readonly details?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
  readonly method_details?: Maybe<Scalars['String']>,
  readonly oligomeric_count?: Maybe<Scalars['Int']>,
  readonly oligomeric_details?: Maybe<Scalars['String']>,
  readonly rcsb_candidate_assembly?: Maybe<Scalars['String']>,
  readonly rcsb_details?: Maybe<Scalars['String']>,
};

export type PdbxStructAssemblyAuthEvidence = {
  readonly __typename?: 'PdbxStructAssemblyAuthEvidence',
  readonly assembly_id: Scalars['String'],
  readonly details?: Maybe<Scalars['String']>,
  readonly experimental_support?: Maybe<Scalars['String']>,
  readonly id: Scalars['String'],
};

export type PdbxStructAssemblyGen = {
  readonly __typename?: 'PdbxStructAssemblyGen',
  readonly assembly_id?: Maybe<Scalars['String']>,
  readonly asym_id_list?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly oper_expression?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
};

export type PdbxStructAssemblyProp = {
  readonly __typename?: 'PdbxStructAssemblyProp',
  readonly assembly_id?: Maybe<Scalars['String']>,
  readonly biol_id: Scalars['String'],
  readonly type: Scalars['String'],
  readonly value?: Maybe<Scalars['String']>,
};

export type PdbxStructOperList = {
  readonly __typename?: 'PdbxStructOperList',
  readonly id: Scalars['String'],
  readonly matrix_1_1?: Maybe<Scalars['Float']>,
  readonly matrix_1_2?: Maybe<Scalars['Float']>,
  readonly matrix_1_3?: Maybe<Scalars['Float']>,
  readonly matrix_2_1?: Maybe<Scalars['Float']>,
  readonly matrix_2_2?: Maybe<Scalars['Float']>,
  readonly matrix_2_3?: Maybe<Scalars['Float']>,
  readonly matrix_3_1?: Maybe<Scalars['Float']>,
  readonly matrix_3_2?: Maybe<Scalars['Float']>,
  readonly matrix_3_3?: Maybe<Scalars['Float']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly symmetry_operation?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
  readonly vector_1?: Maybe<Scalars['Float']>,
  readonly vector_2?: Maybe<Scalars['Float']>,
  readonly vector_3?: Maybe<Scalars['Float']>,
};

export type PdbxStructSpecialSymmetry = {
  readonly __typename?: 'PdbxStructSpecialSymmetry',
  readonly PDB_model_num?: Maybe<Scalars['Int']>,
  readonly auth_seq_id?: Maybe<Scalars['String']>,
  readonly id: Scalars['Int'],
  readonly label_asym_id?: Maybe<Scalars['String']>,
  readonly label_comp_id?: Maybe<Scalars['String']>,
};

export type PdbxVrptSummary = {
  readonly __typename?: 'PdbxVrptSummary',
  readonly B_factor_type?: Maybe<Scalars['String']>,
  readonly Babinet_b?: Maybe<Scalars['Float']>,
  readonly Babinet_k?: Maybe<Scalars['Float']>,
  readonly CA_ONLY?: Maybe<Scalars['String']>,
  readonly DCC_R?: Maybe<Scalars['Float']>,
  readonly DCC_Rfree?: Maybe<Scalars['Float']>,
  readonly DCC_refinement_program?: Maybe<Scalars['String']>,
  readonly EDS_R?: Maybe<Scalars['Float']>,
  readonly EDS_resolution?: Maybe<Scalars['Float']>,
  readonly EDS_resolution_low?: Maybe<Scalars['Float']>,
  readonly Fo_Fc_correlation?: Maybe<Scalars['Float']>,
  readonly I_over_sigma?: Maybe<Scalars['String']>,
  readonly PDB_R?: Maybe<Scalars['Float']>,
  readonly PDB_Rfree?: Maybe<Scalars['Float']>,
  readonly PDB_deposition_date?: Maybe<Scalars['Date']>,
  readonly PDB_resolution?: Maybe<Scalars['Float']>,
  readonly PDB_resolution_low?: Maybe<Scalars['Float']>,
  readonly PDB_revision_date?: Maybe<Scalars['Date']>,
  readonly PDB_revision_number?: Maybe<Scalars['Float']>,
  readonly RNA_suiteness?: Maybe<Scalars['Float']>,
  readonly Wilson_B_aniso?: Maybe<Scalars['String']>,
  readonly Wilson_B_estimate?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_DCC_Rfree?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_RNA_suiteness?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_clashscore?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>,
  readonly absolute_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>,
  readonly acentric_outliers?: Maybe<Scalars['Int']>,
  readonly angles_RMSZ?: Maybe<Scalars['Float']>,
  readonly attempted_validation_steps?: Maybe<Scalars['String']>,
  readonly bonds_RMSZ?: Maybe<Scalars['Float']>,
  readonly bulk_solvent_b?: Maybe<Scalars['Float']>,
  readonly bulk_solvent_k?: Maybe<Scalars['Float']>,
  readonly ccp4version?: Maybe<Scalars['String']>,
  readonly centric_outliers?: Maybe<Scalars['Float']>,
  readonly chemical_shift_completeness?: Maybe<Scalars['Float']>,
  readonly chemical_shift_completeness_full_length?: Maybe<Scalars['Float']>,
  readonly chemical_shifts_input_filename?: Maybe<Scalars['String']>,
  readonly clashscore?: Maybe<Scalars['Float']>,
  readonly clashscore_full_length?: Maybe<Scalars['Float']>,
  readonly coordinates_input_filename?: Maybe<Scalars['String']>,
  readonly cyrange_error?: Maybe<Scalars['String']>,
  readonly cyrange_number_of_domains?: Maybe<Scalars['Int']>,
  readonly cyrange_version?: Maybe<Scalars['String']>,
  readonly data_anisotropy?: Maybe<Scalars['Float']>,
  readonly data_completeness?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_clashscore?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>,
  readonly high_resol_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>,
  readonly ligands_for_buster_report?: Maybe<Scalars['String']>,
  readonly low_resol_relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>,
  readonly low_resol_relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>,
  readonly low_resol_relative_percentile_clashscore?: Maybe<Scalars['Float']>,
  readonly low_resol_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>,
  readonly low_resol_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>,
  readonly low_resol_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>,
  readonly medoid_model?: Maybe<Scalars['Int']>,
  readonly nmr_models_consistency_flag?: Maybe<Scalars['String']>,
  readonly nmrclust_error?: Maybe<Scalars['String']>,
  readonly nmrclust_number_of_clusters?: Maybe<Scalars['Int']>,
  readonly nmrclust_number_of_models?: Maybe<Scalars['Int']>,
  readonly nmrclust_number_of_outliers?: Maybe<Scalars['Int']>,
  readonly nmrclust_representative_model?: Maybe<Scalars['Int']>,
  readonly nmrclust_version?: Maybe<Scalars['String']>,
  readonly no_ligands_for_buster_report?: Maybe<Scalars['String']>,
  readonly no_ligands_for_mogul?: Maybe<Scalars['String']>,
  readonly no_percentile_property?: Maybe<Scalars['String']>,
  readonly num_H_reduce?: Maybe<Scalars['Float']>,
  readonly num_PDBids_absolute_percentile_DCC_Rfree?: Maybe<Scalars['Int']>,
  readonly num_PDBids_absolute_percentile_RNA_suiteness?: Maybe<Scalars['Int']>,
  readonly num_PDBids_absolute_percentile_clashscore?: Maybe<Scalars['Int']>,
  readonly num_PDBids_absolute_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Int']>,
  readonly num_PDBids_absolute_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Int']>,
  readonly num_PDBids_absolute_percentile_percent_rotamer_outliers?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_DCC_Rfree?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_RNA_suiteness?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_clashscore?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Int']>,
  readonly num_PDBids_relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Int']>,
  readonly num_angles_RMSZ?: Maybe<Scalars['Int']>,
  readonly num_bonds_RMSZ?: Maybe<Scalars['Int']>,
  readonly num_free_reflections?: Maybe<Scalars['Int']>,
  readonly num_miller_indices?: Maybe<Scalars['Int']>,
  readonly panav_version?: Maybe<Scalars['String']>,
  readonly percent_RSRZ_outliers?: Maybe<Scalars['Float']>,
  readonly percent_free_reflections?: Maybe<Scalars['Float']>,
  readonly percent_ramachandran_outliers?: Maybe<Scalars['Float']>,
  readonly percent_ramachandran_outliers_full_length?: Maybe<Scalars['Float']>,
  readonly percent_rotamer_outliers?: Maybe<Scalars['Float']>,
  readonly percent_rotamer_outliers_full_length?: Maybe<Scalars['Float']>,
  readonly percentilebins?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly protein_DNA_RNA_entities?: Maybe<Scalars['String']>,
  readonly rci_version?: Maybe<Scalars['String']>,
  readonly reflections_input_filename?: Maybe<Scalars['String']>,
  readonly refmac_version?: Maybe<Scalars['String']>,
  readonly relative_percentile_DCC_Rfree?: Maybe<Scalars['Float']>,
  readonly relative_percentile_RNA_suiteness?: Maybe<Scalars['Float']>,
  readonly relative_percentile_clashscore?: Maybe<Scalars['Float']>,
  readonly relative_percentile_percent_RSRZ_outliers?: Maybe<Scalars['Float']>,
  readonly relative_percentile_percent_ramachandran_outliers?: Maybe<Scalars['Float']>,
  readonly relative_percentile_percent_rotamer_outliers?: Maybe<Scalars['Float']>,
  readonly report_creation_date?: Maybe<Scalars['String']>,
  readonly resol_high_from_reflectionsfile?: Maybe<Scalars['Float']>,
  readonly resol_low_from_reflectionsfile?: Maybe<Scalars['Float']>,
  readonly restypes_notchecked_for_bond_angle_geometry?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly shiftchecker_version?: Maybe<Scalars['String']>,
  readonly trans_NSC?: Maybe<Scalars['String']>,
  readonly twin_L?: Maybe<Scalars['Float']>,
  readonly twin_L2?: Maybe<Scalars['Float']>,
  readonly twin_fraction?: Maybe<Scalars['String']>,
  readonly xtriage_input_columns?: Maybe<Scalars['String']>,
};

export type Query = {
  readonly __typename?: 'Query',
  readonly polymer_entity_instance?: Maybe<CorePolymerEntityInstance>,
  readonly assemblies?: Maybe<ReadonlyArray<Maybe<CoreAssembly>>>,
  readonly nonpolymer_entities?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntity>>>,
  readonly nonpolymer_entity_instance?: Maybe<CoreNonpolymerEntityInstance>,
  readonly polymer_entities?: Maybe<ReadonlyArray<Maybe<CorePolymerEntity>>>,
  readonly polymer_entity?: Maybe<CorePolymerEntity>,
  readonly chem_comp?: Maybe<CoreChemComp>,
  readonly entry?: Maybe<CoreEntry>,
  readonly entries?: Maybe<ReadonlyArray<Maybe<CoreEntry>>>,
  readonly pubmed?: Maybe<CorePubmed>,
  readonly assembly?: Maybe<CoreAssembly>,
  readonly uniprot?: Maybe<CoreUniprot>,
  readonly nonpolymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CoreNonpolymerEntityInstance>>>,
  readonly polymer_entity_instances?: Maybe<ReadonlyArray<Maybe<CorePolymerEntityInstance>>>,
  readonly nonpolymer_entity?: Maybe<CoreNonpolymerEntity>,
};


export type QueryPolymer_Entity_InstanceArgs = {
  asym_id: Scalars['String'],
  entry_id: Scalars['String']
};


export type QueryAssembliesArgs = {
  assembly_ids: ReadonlyArray<Maybe<Scalars['String']>>
};


export type QueryNonpolymer_EntitiesArgs = {
  entity_ids: ReadonlyArray<Scalars['String']>
};


export type QueryNonpolymer_Entity_InstanceArgs = {
  asym_id: Scalars['String'],
  entry_id: Scalars['String']
};


export type QueryPolymer_EntitiesArgs = {
  entity_ids: ReadonlyArray<Scalars['String']>
};


export type QueryPolymer_EntityArgs = {
  entity_id: Scalars['String'],
  entry_id: Scalars['String']
};


export type QueryChem_CompArgs = {
  comp_id: Scalars['String']
};


export type QueryEntryArgs = {
  entry_id: Scalars['String']
};


export type QueryEntriesArgs = {
  entry_ids: ReadonlyArray<Scalars['String']>
};


export type QueryPubmedArgs = {
  pubmed_id: Scalars['Int']
};


export type QueryAssemblyArgs = {
  assembly_id: Scalars['String'],
  entry_id: Scalars['String']
};


export type QueryUniprotArgs = {
  uniprot_id: Scalars['String']
};


export type QueryNonpolymer_Entity_InstancesArgs = {
  instance_ids: ReadonlyArray<Maybe<Scalars['String']>>
};


export type QueryPolymer_Entity_InstancesArgs = {
  instance_ids: ReadonlyArray<Maybe<Scalars['String']>>
};


export type QueryNonpolymer_EntityArgs = {
  entity_id: Scalars['String'],
  entry_id: Scalars['String']
};

export type RcsbAccessionInfo = {
  readonly __typename?: 'RcsbAccessionInfo',
  readonly deposit_date?: Maybe<Scalars['Date']>,
  readonly initial_release_date?: Maybe<Scalars['Date']>,
  readonly major_revision?: Maybe<Scalars['Int']>,
  readonly minor_revision?: Maybe<Scalars['Int']>,
  readonly revision_date?: Maybe<Scalars['Date']>,
  readonly status_code?: Maybe<Scalars['String']>,
};

export type RcsbAssemblyContainerIdentifiers = {
  readonly __typename?: 'RcsbAssemblyContainerIdentifiers',
  readonly assembly_id: Scalars['String'],
  readonly entry_id: Scalars['String'],
  readonly rcsb_id?: Maybe<Scalars['String']>,
};

export type RcsbAssemblyInfo = {
  readonly __typename?: 'RcsbAssemblyInfo',
  readonly assembly_id?: Maybe<Scalars['String']>,
  readonly atom_count?: Maybe<Scalars['Int']>,
  readonly branched_atom_count?: Maybe<Scalars['Int']>,
  readonly branched_entity_count?: Maybe<Scalars['Int']>,
  readonly branched_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly entry_id: Scalars['String'],
  readonly modeled_polymer_monomer_count?: Maybe<Scalars['Int']>,
  readonly na_polymer_entity_types?: Maybe<Scalars['String']>,
  readonly nonpolymer_atom_count?: Maybe<Scalars['Int']>,
  readonly nonpolymer_entity_count?: Maybe<Scalars['Int']>,
  readonly nonpolymer_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly polymer_atom_count?: Maybe<Scalars['Int']>,
  readonly polymer_composition?: Maybe<Scalars['String']>,
  readonly polymer_entity_count?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_DNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_RNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_nucleic_acid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_protein?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count_DNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count_RNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count_nucleic_acid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_instance_count_protein?: Maybe<Scalars['Int']>,
  readonly polymer_monomer_count?: Maybe<Scalars['Int']>,
  readonly selected_polymer_entity_types?: Maybe<Scalars['String']>,
  readonly solvent_atom_count?: Maybe<Scalars['Int']>,
  readonly solvent_entity_count?: Maybe<Scalars['Int']>,
  readonly solvent_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly unmodeled_polymer_monomer_count?: Maybe<Scalars['Int']>,
};

export type RcsbBindingAffinity = {
  readonly __typename?: 'RcsbBindingAffinity',
  readonly comp_id: Scalars['String'],
  readonly display_order: Scalars['Int'],
  readonly display_value: Scalars['String'],
  readonly link: Scalars['String'],
  readonly provenance_code: Scalars['String'],
  readonly reference_sequence_identity?: Maybe<Scalars['Int']>,
  readonly symbol?: Maybe<Scalars['String']>,
  readonly type: Scalars['String'],
  readonly unit: Scalars['String'],
  readonly value: Scalars['Float'],
};

export type RcsbBirdCitation = {
  readonly __typename?: 'RcsbBirdCitation',
  readonly id: Scalars['String'],
  readonly journal_abbrev?: Maybe<Scalars['String']>,
  readonly journal_volume?: Maybe<Scalars['String']>,
  readonly page_first?: Maybe<Scalars['String']>,
  readonly page_last?: Maybe<Scalars['String']>,
  readonly pdbx_database_id_DOI?: Maybe<Scalars['String']>,
  readonly pdbx_database_id_PubMed?: Maybe<Scalars['Int']>,
  readonly rcsb_authors?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly title?: Maybe<Scalars['String']>,
  readonly year?: Maybe<Scalars['Int']>,
};

export type RcsbChemCompContainerIdentifiers = {
  readonly __typename?: 'RcsbChemCompContainerIdentifiers',
  readonly atc_codes?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly comp_id: Scalars['String'],
  readonly drugbank_id?: Maybe<Scalars['String']>,
  readonly prd_id?: Maybe<Scalars['String']>,
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly subcomponent_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbChemCompDescriptor = {
  readonly __typename?: 'RcsbChemCompDescriptor',
  readonly InChI?: Maybe<Scalars['String']>,
  readonly InChIKey?: Maybe<Scalars['String']>,
  readonly SMILES?: Maybe<Scalars['String']>,
  readonly SMILES_stereo?: Maybe<Scalars['String']>,
  readonly comp_id: Scalars['String'],
};

export type RcsbChemCompInfo = {
  readonly __typename?: 'RcsbChemCompInfo',
  readonly atom_count?: Maybe<Scalars['Int']>,
  readonly atom_count_chiral?: Maybe<Scalars['Int']>,
  readonly atom_count_heavy?: Maybe<Scalars['Int']>,
  readonly bond_count?: Maybe<Scalars['Int']>,
  readonly bond_count_aromatic?: Maybe<Scalars['Int']>,
  readonly comp_id: Scalars['String'],
  readonly initial_release_date?: Maybe<Scalars['Date']>,
  readonly release_status?: Maybe<Scalars['String']>,
  readonly revision_date?: Maybe<Scalars['Date']>,
};

export type RcsbChemCompRelated = {
  readonly __typename?: 'RcsbChemCompRelated',
  readonly comp_id: Scalars['String'],
  readonly ordinal: Scalars['Int'],
  readonly related_mapping_method?: Maybe<Scalars['String']>,
  readonly resource_accession_code?: Maybe<Scalars['String']>,
  readonly resource_lineage?: Maybe<ReadonlyArray<Maybe<RcsbChemCompRelatedResourceLineage>>>,
  readonly resource_name?: Maybe<Scalars['String']>,
};

export type RcsbChemCompRelatedResourceLineage = {
  readonly __typename?: 'RcsbChemCompRelatedResourceLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbChemCompSynonyms = {
  readonly __typename?: 'RcsbChemCompSynonyms',
  readonly comp_id: Scalars['String'],
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
};

export type RcsbChemCompTarget = {
  readonly __typename?: 'RcsbChemCompTarget',
  readonly comp_id: Scalars['String'],
  readonly interaction_type?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly reference_database_accession_code?: Maybe<Scalars['String']>,
  readonly reference_database_name?: Maybe<Scalars['String']>,
  readonly target_actions?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbClusterFlexibility = {
  readonly __typename?: 'RcsbClusterFlexibility',
  readonly avg_rmsd?: Maybe<Scalars['Float']>,
  readonly label?: Maybe<Scalars['String']>,
  readonly link?: Maybe<Scalars['String']>,
  readonly max_rmsd?: Maybe<Scalars['Float']>,
  readonly provenance_code?: Maybe<Scalars['String']>,
};

export type RcsbClusterMembership = {
  readonly __typename?: 'RcsbClusterMembership',
  readonly cluster_id: Scalars['Int'],
  readonly identity: Scalars['Int'],
};

export type RcsbEntityHostOrganism = {
  readonly __typename?: 'RcsbEntityHostOrganism',
  readonly beg_seq_num?: Maybe<Scalars['Int']>,
  readonly common_name?: Maybe<Scalars['String']>,
  readonly end_seq_num?: Maybe<Scalars['Int']>,
  readonly ncbi_common_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly ncbi_parent_scientific_name?: Maybe<Scalars['String']>,
  readonly ncbi_scientific_name?: Maybe<Scalars['String']>,
  readonly ncbi_taxonomy_id?: Maybe<Scalars['Int']>,
  readonly pdbx_src_id: Scalars['String'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly scientific_name?: Maybe<Scalars['String']>,
  readonly taxonomy_lineage?: Maybe<ReadonlyArray<Maybe<RcsbEntityHostOrganismTaxonomyLineage>>>,
};

export type RcsbEntityHostOrganismTaxonomyLineage = {
  readonly __typename?: 'RcsbEntityHostOrganismTaxonomyLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbEntitySourceOrganism = {
  readonly __typename?: 'RcsbEntitySourceOrganism',
  readonly beg_seq_num?: Maybe<Scalars['Int']>,
  readonly common_name?: Maybe<Scalars['String']>,
  readonly end_seq_num?: Maybe<Scalars['Int']>,
  readonly ncbi_common_names?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly ncbi_parent_scientific_name?: Maybe<Scalars['String']>,
  readonly ncbi_scientific_name?: Maybe<Scalars['String']>,
  readonly ncbi_taxonomy_id?: Maybe<Scalars['Int']>,
  readonly pdbx_src_id: Scalars['String'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly rcsb_gene_name?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganismRcsbGeneName>>>,
  readonly scientific_name?: Maybe<Scalars['String']>,
  readonly source_type?: Maybe<Scalars['String']>,
  readonly taxonomy_lineage?: Maybe<ReadonlyArray<Maybe<RcsbEntitySourceOrganismTaxonomyLineage>>>,
};

export type RcsbEntitySourceOrganismRcsbGeneName = {
  readonly __typename?: 'RcsbEntitySourceOrganismRcsbGeneName',
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly value?: Maybe<Scalars['String']>,
};

export type RcsbEntitySourceOrganismTaxonomyLineage = {
  readonly __typename?: 'RcsbEntitySourceOrganismTaxonomyLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbEntryContainerIdentifiers = {
  readonly __typename?: 'RcsbEntryContainerIdentifiers',
  readonly assembly_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly branched_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly emdb_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly entry_id: Scalars['String'],
  readonly non_polymer_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly polymer_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly pubmed_id?: Maybe<Scalars['Int']>,
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly related_emdb_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly water_entity_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbEntryInfo = {
  readonly __typename?: 'RcsbEntryInfo',
  readonly assembly_count?: Maybe<Scalars['Int']>,
  readonly branched_entity_count?: Maybe<Scalars['Int']>,
  readonly branched_molecular_weight_maximum?: Maybe<Scalars['Float']>,
  readonly branched_molecular_weight_minimum?: Maybe<Scalars['Float']>,
  readonly cis_peptide_count?: Maybe<Scalars['Int']>,
  readonly deposited_atom_count?: Maybe<Scalars['Int']>,
  readonly deposited_model_count?: Maybe<Scalars['Int']>,
  readonly deposited_modeled_polymer_monomer_count?: Maybe<Scalars['Int']>,
  readonly deposited_nonpolymer_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly deposited_polymer_entity_instance_count?: Maybe<Scalars['Int']>,
  readonly deposited_polymer_monomer_count?: Maybe<Scalars['Int']>,
  readonly deposited_unmodeled_polymer_monomer_count?: Maybe<Scalars['Int']>,
  readonly disulfide_bond_count?: Maybe<Scalars['Int']>,
  readonly entity_count?: Maybe<Scalars['Int']>,
  readonly experimental_method?: Maybe<Scalars['String']>,
  readonly experimental_method_count?: Maybe<Scalars['Int']>,
  readonly inter_mol_covalent_bond_count?: Maybe<Scalars['Int']>,
  readonly inter_mol_metalic_bond_count?: Maybe<Scalars['Int']>,
  readonly molecular_weight?: Maybe<Scalars['Float']>,
  readonly na_polymer_entity_types?: Maybe<Scalars['String']>,
  readonly nonpolymer_bound_components?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly nonpolymer_entity_count?: Maybe<Scalars['Int']>,
  readonly nonpolymer_molecular_weight_maximum?: Maybe<Scalars['Float']>,
  readonly nonpolymer_molecular_weight_minimum?: Maybe<Scalars['Float']>,
  readonly polymer_composition?: Maybe<Scalars['String']>,
  readonly polymer_entity_count?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_DNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_RNA?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_nucleic_acid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_nucleic_acid_hybrid?: Maybe<Scalars['Int']>,
  readonly polymer_entity_count_protein?: Maybe<Scalars['Int']>,
  readonly polymer_entity_taxonomy_count?: Maybe<Scalars['Int']>,
  readonly polymer_molecular_weight_maximum?: Maybe<Scalars['Float']>,
  readonly polymer_molecular_weight_minimum?: Maybe<Scalars['Float']>,
  readonly polymer_monomer_count_maximum?: Maybe<Scalars['Int']>,
  readonly polymer_monomer_count_minimum?: Maybe<Scalars['Int']>,
  readonly resolution_combined?: Maybe<ReadonlyArray<Maybe<Scalars['Float']>>>,
  readonly selected_polymer_entity_types?: Maybe<Scalars['String']>,
  readonly software_programs_combined?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly solvent_entity_count?: Maybe<Scalars['Int']>,
};

export type RcsbExternalReferences = {
  readonly __typename?: 'RcsbExternalReferences',
  readonly id: Scalars['String'],
  readonly link: Scalars['String'],
  readonly type: Scalars['String'],
};

export type RcsbGenomicLineage = {
  readonly __typename?: 'RcsbGenomicLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbLatestRevision = {
  readonly __typename?: 'RcsbLatestRevision',
  readonly major_revision?: Maybe<Scalars['Int']>,
  readonly minor_revision?: Maybe<Scalars['Int']>,
  readonly revision_date?: Maybe<Scalars['Date']>,
};

export type RcsbMembraneLineage = {
  readonly __typename?: 'RcsbMembraneLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntity = {
  readonly __typename?: 'RcsbNonpolymerEntity',
  readonly details?: Maybe<Scalars['String']>,
  readonly formula_weight?: Maybe<Scalars['Float']>,
  readonly pdbx_description?: Maybe<Scalars['String']>,
  readonly pdbx_number_of_molecules?: Maybe<Scalars['Int']>,
};

export type RcsbNonpolymerEntityAnnotation = {
  readonly __typename?: 'RcsbNonpolymerEntityAnnotation',
  readonly annotation_id?: Maybe<Scalars['String']>,
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerEntityAnnotationAnnotationLineage>>>,
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntityAnnotationAnnotationLineage = {
  readonly __typename?: 'RcsbNonpolymerEntityAnnotationAnnotationLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntityContainerIdentifiers = {
  readonly __typename?: 'RcsbNonpolymerEntityContainerIdentifiers',
  readonly asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly auth_asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly entity_id: Scalars['String'],
  readonly entry_id: Scalars['String'],
  readonly nonpolymer_comp_id?: Maybe<Scalars['String']>,
  readonly prd_id?: Maybe<Scalars['String']>,
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly reference_chemical_identifiers_provenance_source?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly reference_chemical_identifiers_resource_accession?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly reference_chemical_identifiers_resource_name?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbNonpolymerEntityFeature = {
  readonly __typename?: 'RcsbNonpolymerEntityFeature',
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly feature_id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
  readonly value?: Maybe<Scalars['Float']>,
};

export type RcsbNonpolymerEntityFeatureSummary = {
  readonly __typename?: 'RcsbNonpolymerEntityFeatureSummary',
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly count?: Maybe<Scalars['Int']>,
  readonly maximum_length?: Maybe<Scalars['Int']>,
  readonly maximum_value?: Maybe<Scalars['Float']>,
  readonly minimum_length?: Maybe<Scalars['Int']>,
  readonly minimum_value?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntityInstanceContainerIdentifiers = {
  readonly __typename?: 'RcsbNonpolymerEntityInstanceContainerIdentifiers',
  readonly asym_id: Scalars['String'],
  readonly auth_asym_id?: Maybe<Scalars['String']>,
  readonly auth_seq_id?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly entity_id?: Maybe<Scalars['String']>,
  readonly entry_id: Scalars['String'],
  readonly rcsb_id?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntityKeywords = {
  readonly __typename?: 'RcsbNonpolymerEntityKeywords',
  readonly text?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerEntityNameCom = {
  readonly __typename?: 'RcsbNonpolymerEntityNameCom',
  readonly name: Scalars['String'],
};

export type RcsbNonpolymerInstanceAnnotation = {
  readonly __typename?: 'RcsbNonpolymerInstanceAnnotation',
  readonly annotation_id?: Maybe<Scalars['String']>,
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceAnnotationAnnotationLineage>>>,
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerInstanceAnnotationAnnotationLineage = {
  readonly __typename?: 'RcsbNonpolymerInstanceAnnotationAnnotationLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerInstanceFeature = {
  readonly __typename?: 'RcsbNonpolymerInstanceFeature',
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly auth_seq_id?: Maybe<Scalars['String']>,
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly feature_id?: Maybe<Scalars['String']>,
  readonly feature_value?: Maybe<ReadonlyArray<Maybe<RcsbNonpolymerInstanceFeatureFeatureValue>>>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerInstanceFeatureFeatureValue = {
  readonly __typename?: 'RcsbNonpolymerInstanceFeatureFeatureValue',
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly reference?: Maybe<Scalars['Float']>,
  readonly reported?: Maybe<Scalars['Float']>,
  readonly uncertainty_estimate?: Maybe<Scalars['Float']>,
  readonly uncertainty_estimate_type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerInstanceFeatureSummary = {
  readonly __typename?: 'RcsbNonpolymerInstanceFeatureSummary',
  readonly count?: Maybe<Scalars['Int']>,
  readonly maximum_length?: Maybe<Scalars['Int']>,
  readonly maximum_value?: Maybe<Scalars['Float']>,
  readonly minimum_length?: Maybe<Scalars['Int']>,
  readonly minimum_value?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerStructConn = {
  readonly __typename?: 'RcsbNonpolymerStructConn',
  readonly connect_partner?: Maybe<RcsbNonpolymerStructConnConnectPartner>,
  readonly connect_target?: Maybe<RcsbNonpolymerStructConnConnectTarget>,
  readonly connect_type?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly dist_value?: Maybe<Scalars['Float']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly ordinal_id: Scalars['Int'],
  readonly value_order?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerStructConnConnectPartner = {
  readonly __typename?: 'RcsbNonpolymerStructConnConnectPartner',
  readonly label_alt_id?: Maybe<Scalars['String']>,
  readonly label_asym_id: Scalars['String'],
  readonly label_atom_id?: Maybe<Scalars['String']>,
  readonly label_comp_id: Scalars['String'],
  readonly label_seq_id?: Maybe<Scalars['Int']>,
  readonly symmetry?: Maybe<Scalars['String']>,
};

export type RcsbNonpolymerStructConnConnectTarget = {
  readonly __typename?: 'RcsbNonpolymerStructConnConnectTarget',
  readonly auth_asym_id?: Maybe<Scalars['String']>,
  readonly auth_seq_id?: Maybe<Scalars['String']>,
  readonly label_alt_id?: Maybe<Scalars['String']>,
  readonly label_asym_id: Scalars['String'],
  readonly label_atom_id?: Maybe<Scalars['String']>,
  readonly label_comp_id: Scalars['String'],
  readonly label_seq_id?: Maybe<Scalars['Int']>,
  readonly symmetry?: Maybe<Scalars['String']>,
};

export type RcsbPfamContainerIdentifiers = {
  readonly __typename?: 'RcsbPfamContainerIdentifiers',
  readonly pfam_id?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntity = {
  readonly __typename?: 'RcsbPolymerEntity',
  readonly details?: Maybe<Scalars['String']>,
  readonly formula_weight?: Maybe<Scalars['Float']>,
  readonly pdbx_description?: Maybe<Scalars['String']>,
  readonly pdbx_ec?: Maybe<Scalars['String']>,
  readonly pdbx_fragment?: Maybe<Scalars['String']>,
  readonly pdbx_mutation?: Maybe<Scalars['String']>,
  readonly pdbx_number_of_molecules?: Maybe<Scalars['Int']>,
  readonly rcsb_ec_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbEcLineage>>>,
  readonly rcsb_enzyme_class_combined?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbEnzymeClassCombined>>>,
  readonly rcsb_macromolecular_names_combined?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityRcsbMacromolecularNamesCombined>>>,
  readonly rcsb_multiple_source_flag?: Maybe<Scalars['String']>,
  readonly rcsb_source_part_count?: Maybe<Scalars['Int']>,
  readonly src_method?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityAlign = {
  readonly __typename?: 'RcsbPolymerEntityAlign',
  readonly aligned_regions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAlignAlignedRegions>>>,
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly reference_database_accession?: Maybe<Scalars['String']>,
  readonly reference_database_isoform?: Maybe<Scalars['String']>,
  readonly reference_database_name?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityAlignAlignedRegions = {
  readonly __typename?: 'RcsbPolymerEntityAlignAlignedRegions',
  readonly entity_beg_seq_id?: Maybe<Scalars['Int']>,
  readonly length?: Maybe<Scalars['Int']>,
  readonly ref_beg_seq_id?: Maybe<Scalars['Int']>,
};

export type RcsbPolymerEntityAnnotation = {
  readonly __typename?: 'RcsbPolymerEntityAnnotation',
  readonly annotation_id?: Maybe<Scalars['String']>,
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityAnnotationAnnotationLineage>>>,
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityAnnotationAnnotationLineage = {
  readonly __typename?: 'RcsbPolymerEntityAnnotationAnnotationLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityContainerIdentifiers = {
  readonly __typename?: 'RcsbPolymerEntityContainerIdentifiers',
  readonly asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly auth_asym_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly chem_comp_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly chem_comp_nstd_monomers?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly entity_id: Scalars['String'],
  readonly entry_id: Scalars['String'],
  readonly prd_id?: Maybe<Scalars['String']>,
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly reference_sequence_identifiers?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityContainerIdentifiersReferenceSequenceIdentifiers>>>,
  readonly uniprot_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbPolymerEntityContainerIdentifiersReferenceSequenceIdentifiers = {
  readonly __typename?: 'RcsbPolymerEntityContainerIdentifiersReferenceSequenceIdentifiers',
  readonly database_accession?: Maybe<Scalars['String']>,
  readonly database_isoform?: Maybe<Scalars['String']>,
  readonly database_name?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityFeature = {
  readonly __typename?: 'RcsbPolymerEntityFeature',
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly feature_id?: Maybe<Scalars['String']>,
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerEntityFeatureFeaturePositions>>>,
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly reference_scheme?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityFeatureFeaturePositions = {
  readonly __typename?: 'RcsbPolymerEntityFeatureFeaturePositions',
  readonly beg_comp_id?: Maybe<Scalars['String']>,
  readonly beg_seq_id: Scalars['Int'],
  readonly end_seq_id?: Maybe<Scalars['Int']>,
  readonly value?: Maybe<Scalars['Float']>,
};

export type RcsbPolymerEntityFeatureSummary = {
  readonly __typename?: 'RcsbPolymerEntityFeatureSummary',
  readonly count?: Maybe<Scalars['Int']>,
  readonly coverage?: Maybe<Scalars['Float']>,
  readonly maximum_length?: Maybe<Scalars['Int']>,
  readonly maximum_value?: Maybe<Scalars['Float']>,
  readonly minimum_length?: Maybe<Scalars['Int']>,
  readonly minimum_value?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityInstanceContainerIdentifiers = {
  readonly __typename?: 'RcsbPolymerEntityInstanceContainerIdentifiers',
  readonly asym_id: Scalars['String'],
  readonly auth_asym_id?: Maybe<Scalars['String']>,
  readonly entity_id?: Maybe<Scalars['String']>,
  readonly entry_id: Scalars['String'],
  readonly rcsb_id?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityKeywords = {
  readonly __typename?: 'RcsbPolymerEntityKeywords',
  readonly text?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityNameCom = {
  readonly __typename?: 'RcsbPolymerEntityNameCom',
  readonly name: Scalars['String'],
};

export type RcsbPolymerEntityNameSys = {
  readonly __typename?: 'RcsbPolymerEntityNameSys',
  readonly name: Scalars['String'],
  readonly system?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityRcsbEcLineage = {
  readonly __typename?: 'RcsbPolymerEntityRcsbEcLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityRcsbEnzymeClassCombined = {
  readonly __typename?: 'RcsbPolymerEntityRcsbEnzymeClassCombined',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly ec?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
};

export type RcsbPolymerEntityRcsbMacromolecularNamesCombined = {
  readonly __typename?: 'RcsbPolymerEntityRcsbMacromolecularNamesCombined',
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_code?: Maybe<Scalars['String']>,
  readonly provenance_source?: Maybe<Scalars['String']>,
};

export type RcsbPolymerInstanceAnnotation = {
  readonly __typename?: 'RcsbPolymerInstanceAnnotation',
  readonly annotation_id?: Maybe<Scalars['String']>,
  readonly annotation_lineage?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceAnnotationAnnotationLineage>>>,
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerInstanceAnnotationAnnotationLineage = {
  readonly __typename?: 'RcsbPolymerInstanceAnnotationAnnotationLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbPolymerInstanceFeature = {
  readonly __typename?: 'RcsbPolymerInstanceFeature',
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly feature_id?: Maybe<Scalars['String']>,
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbPolymerInstanceFeatureFeaturePositions>>>,
  readonly name?: Maybe<Scalars['String']>,
  readonly ordinal: Scalars['Int'],
  readonly provenance_source?: Maybe<Scalars['String']>,
  readonly reference_scheme?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerInstanceFeatureFeaturePositions = {
  readonly __typename?: 'RcsbPolymerInstanceFeatureFeaturePositions',
  readonly beg_comp_id?: Maybe<Scalars['String']>,
  readonly beg_seq_id: Scalars['Int'],
  readonly end_seq_id?: Maybe<Scalars['Int']>,
  readonly value?: Maybe<Scalars['Float']>,
};

export type RcsbPolymerInstanceFeatureSummary = {
  readonly __typename?: 'RcsbPolymerInstanceFeatureSummary',
  readonly count?: Maybe<Scalars['Int']>,
  readonly coverage?: Maybe<Scalars['Float']>,
  readonly maximum_length?: Maybe<Scalars['Int']>,
  readonly maximum_value?: Maybe<Scalars['Float']>,
  readonly minimum_length?: Maybe<Scalars['Int']>,
  readonly minimum_value?: Maybe<Scalars['Float']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbPolymerStructConn = {
  readonly __typename?: 'RcsbPolymerStructConn',
  readonly connect_partner?: Maybe<RcsbPolymerStructConnConnectPartner>,
  readonly connect_target?: Maybe<RcsbPolymerStructConnConnectTarget>,
  readonly connect_type?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly dist_value?: Maybe<Scalars['Float']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly ordinal_id: Scalars['Int'],
  readonly value_order?: Maybe<Scalars['String']>,
};

export type RcsbPolymerStructConnConnectPartner = {
  readonly __typename?: 'RcsbPolymerStructConnConnectPartner',
  readonly label_alt_id?: Maybe<Scalars['String']>,
  readonly label_asym_id: Scalars['String'],
  readonly label_atom_id?: Maybe<Scalars['String']>,
  readonly label_comp_id: Scalars['String'],
  readonly label_seq_id?: Maybe<Scalars['Int']>,
  readonly symmetry?: Maybe<Scalars['String']>,
};

export type RcsbPolymerStructConnConnectTarget = {
  readonly __typename?: 'RcsbPolymerStructConnConnectTarget',
  readonly auth_asym_id?: Maybe<Scalars['String']>,
  readonly auth_seq_id?: Maybe<Scalars['String']>,
  readonly label_alt_id?: Maybe<Scalars['String']>,
  readonly label_asym_id: Scalars['String'],
  readonly label_atom_id?: Maybe<Scalars['String']>,
  readonly label_comp_id: Scalars['String'],
  readonly label_seq_id?: Maybe<Scalars['Int']>,
  readonly symmetry?: Maybe<Scalars['String']>,
};

export type RcsbPubmedContainerIdentifiers = {
  readonly __typename?: 'RcsbPubmedContainerIdentifiers',
  readonly pubmed_id?: Maybe<Scalars['Int']>,
};

export type RcsbPubmedMeshDescriptorsLineage = {
  readonly __typename?: 'RcsbPubmedMeshDescriptorsLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbRepositoryHoldingsCurrent = {
  readonly __typename?: 'RcsbRepositoryHoldingsCurrent',
  readonly repository_content_types?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
};

export type RcsbRepositoryHoldingsCurrentEntryContainerIdentifiers = {
  readonly __typename?: 'RcsbRepositoryHoldingsCurrentEntryContainerIdentifiers',
  readonly assembly_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly entry_id: Scalars['String'],
  readonly rcsb_id?: Maybe<Scalars['String']>,
  readonly update_id?: Maybe<Scalars['String']>,
};

export type RcsbSchemaContainerIdentifiers = {
  readonly __typename?: 'RcsbSchemaContainerIdentifiers',
  readonly collection_name: Scalars['String'],
  readonly collection_schema_version?: Maybe<Scalars['String']>,
  readonly schema_name: Scalars['String'],
};

export type RcsbStructSymmetry = {
  readonly __typename?: 'RcsbStructSymmetry',
  readonly clusters: ReadonlyArray<Maybe<RcsbStructSymmetryClusters>>,
  readonly kind: Scalars['String'],
  readonly oligomeric_state: Scalars['String'],
  readonly rotation_axes?: Maybe<ReadonlyArray<Maybe<RcsbStructSymmetryRotationAxes>>>,
  readonly stoichiometry: ReadonlyArray<Maybe<Scalars['String']>>,
  readonly symbol: Scalars['String'],
  readonly type: Scalars['String'],
};

export type RcsbStructSymmetryClusters = {
  readonly __typename?: 'RcsbStructSymmetryClusters',
  readonly avg_rmsd?: Maybe<Scalars['Float']>,
  readonly members: ReadonlyArray<Maybe<ClustersMembers>>,
};

export type RcsbStructSymmetryLineage = {
  readonly __typename?: 'RcsbStructSymmetryLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbStructSymmetryRotationAxes = {
  readonly __typename?: 'RcsbStructSymmetryRotationAxes',
  readonly end: ReadonlyArray<Maybe<Scalars['Float']>>,
  readonly order?: Maybe<Scalars['Int']>,
  readonly start: ReadonlyArray<Maybe<Scalars['Float']>>,
};

export type RcsbUniprotContainerIdentifiers = {
  readonly __typename?: 'RcsbUniprotContainerIdentifiers',
  readonly ensembl_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly go_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly pfam_ids?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly uniprot_id?: Maybe<Scalars['String']>,
};

export type RcsbUniprotFeature = {
  readonly __typename?: 'RcsbUniprotFeature',
  readonly assignment_version?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly feature_class_lineage?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeatureFeatureClassLineage>>>,
  readonly feature_id?: Maybe<Scalars['String']>,
  readonly feature_positions?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeatureFeaturePositions>>>,
  readonly feature_ranges?: Maybe<ReadonlyArray<Maybe<RcsbUniprotFeatureFeatureRanges>>>,
  readonly name?: Maybe<Scalars['String']>,
  readonly provenance_code?: Maybe<Scalars['String']>,
  readonly reference_scheme?: Maybe<Scalars['String']>,
  readonly type?: Maybe<Scalars['String']>,
};

export type RcsbUniprotFeatureFeatureClassLineage = {
  readonly __typename?: 'RcsbUniprotFeatureFeatureClassLineage',
  readonly depth?: Maybe<Scalars['Int']>,
  readonly id?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
};

export type RcsbUniprotFeatureFeaturePositions = {
  readonly __typename?: 'RcsbUniprotFeatureFeaturePositions',
  readonly comp_id?: Maybe<Scalars['String']>,
  readonly seq_id?: Maybe<Scalars['Int']>,
  readonly value?: Maybe<Scalars['Float']>,
};

export type RcsbUniprotFeatureFeatureRanges = {
  readonly __typename?: 'RcsbUniprotFeatureFeatureRanges',
  readonly beg_seq_id?: Maybe<Scalars['Int']>,
  readonly end_seq_id?: Maybe<Scalars['Int']>,
  readonly value?: Maybe<Scalars['Float']>,
};

export type RcsbUniprotKeyword = {
  readonly __typename?: 'RcsbUniprotKeyword',
  readonly id?: Maybe<Scalars['String']>,
  readonly value?: Maybe<Scalars['String']>,
};

export type RcsbUniprotProtein = {
  readonly __typename?: 'RcsbUniprotProtein',
  readonly ec?: Maybe<ReadonlyArray<Maybe<RcsbUniprotProteinEc>>>,
  readonly function?: Maybe<RcsbUniprotProteinFunction>,
  readonly gene?: Maybe<ReadonlyArray<Maybe<RcsbUniprotProteinGene>>>,
  readonly name?: Maybe<RcsbUniprotProteinName>,
  readonly sequence?: Maybe<Scalars['String']>,
  readonly source_organism?: Maybe<RcsbUniprotProteinSourceOrganism>,
};

export type RcsbUniprotProteinEc = {
  readonly __typename?: 'RcsbUniprotProteinEc',
  readonly number?: Maybe<Scalars['String']>,
  readonly provenance_code?: Maybe<Scalars['String']>,
};

export type RcsbUniprotProteinFunction = {
  readonly __typename?: 'RcsbUniprotProteinFunction',
  readonly details?: Maybe<Scalars['String']>,
  readonly provenance_code?: Maybe<Scalars['String']>,
};

export type RcsbUniprotProteinGene = {
  readonly __typename?: 'RcsbUniprotProteinGene',
  readonly name?: Maybe<ReadonlyArray<Maybe<GeneName>>>,
};

export type RcsbUniprotProteinName = {
  readonly __typename?: 'RcsbUniprotProteinName',
  readonly provenance_code: Scalars['String'],
  readonly value: Scalars['String'],
};

export type RcsbUniprotProteinSourceOrganism = {
  readonly __typename?: 'RcsbUniprotProteinSourceOrganism',
  readonly provenance_code: Scalars['String'],
  readonly scientific_name: Scalars['String'],
  readonly taxonomy_id?: Maybe<Scalars['Int']>,
};

export type Refine = {
  readonly __typename?: 'Refine',
  readonly B_iso_max?: Maybe<Scalars['Float']>,
  readonly B_iso_mean?: Maybe<Scalars['Float']>,
  readonly B_iso_min?: Maybe<Scalars['Float']>,
  readonly aniso_B_1_1?: Maybe<Scalars['Float']>,
  readonly aniso_B_1_2?: Maybe<Scalars['Float']>,
  readonly aniso_B_1_3?: Maybe<Scalars['Float']>,
  readonly aniso_B_2_2?: Maybe<Scalars['Float']>,
  readonly aniso_B_2_3?: Maybe<Scalars['Float']>,
  readonly aniso_B_3_3?: Maybe<Scalars['Float']>,
  readonly correlation_coeff_Fo_to_Fc?: Maybe<Scalars['Float']>,
  readonly correlation_coeff_Fo_to_Fc_free?: Maybe<Scalars['Float']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly ls_R_factor_R_free?: Maybe<Scalars['Float']>,
  readonly ls_R_factor_R_free_error?: Maybe<Scalars['Float']>,
  readonly ls_R_factor_R_free_error_details?: Maybe<Scalars['String']>,
  readonly ls_R_factor_R_work?: Maybe<Scalars['Float']>,
  readonly ls_R_factor_all?: Maybe<Scalars['Float']>,
  readonly ls_R_factor_obs?: Maybe<Scalars['Float']>,
  readonly ls_d_res_high?: Maybe<Scalars['Float']>,
  readonly ls_d_res_low?: Maybe<Scalars['Float']>,
  readonly ls_matrix_type?: Maybe<Scalars['String']>,
  readonly ls_number_parameters?: Maybe<Scalars['Int']>,
  readonly ls_number_reflns_R_free?: Maybe<Scalars['Int']>,
  readonly ls_number_reflns_R_work?: Maybe<Scalars['Int']>,
  readonly ls_number_reflns_all?: Maybe<Scalars['Int']>,
  readonly ls_number_reflns_obs?: Maybe<Scalars['Int']>,
  readonly ls_number_restraints?: Maybe<Scalars['Int']>,
  readonly ls_percent_reflns_R_free?: Maybe<Scalars['Float']>,
  readonly ls_percent_reflns_obs?: Maybe<Scalars['Float']>,
  readonly ls_redundancy_reflns_all?: Maybe<Scalars['Float']>,
  readonly ls_redundancy_reflns_obs?: Maybe<Scalars['Float']>,
  readonly ls_wR_factor_R_free?: Maybe<Scalars['Float']>,
  readonly ls_wR_factor_R_work?: Maybe<Scalars['Float']>,
  readonly occupancy_max?: Maybe<Scalars['Float']>,
  readonly occupancy_min?: Maybe<Scalars['Float']>,
  readonly overall_FOM_free_R_set?: Maybe<Scalars['Float']>,
  readonly overall_FOM_work_R_set?: Maybe<Scalars['Float']>,
  readonly overall_SU_B?: Maybe<Scalars['Float']>,
  readonly overall_SU_ML?: Maybe<Scalars['Float']>,
  readonly overall_SU_R_Cruickshank_DPI?: Maybe<Scalars['Float']>,
  readonly overall_SU_R_free?: Maybe<Scalars['Float']>,
  readonly pdbx_R_Free_selection_details?: Maybe<Scalars['String']>,
  readonly pdbx_TLS_residual_ADP_flag?: Maybe<Scalars['String']>,
  readonly pdbx_average_fsc_free?: Maybe<Scalars['Float']>,
  readonly pdbx_average_fsc_overall?: Maybe<Scalars['Float']>,
  readonly pdbx_average_fsc_work?: Maybe<Scalars['Float']>,
  readonly pdbx_data_cutoff_high_absF?: Maybe<Scalars['Float']>,
  readonly pdbx_data_cutoff_high_rms_absF?: Maybe<Scalars['Float']>,
  readonly pdbx_data_cutoff_low_absF?: Maybe<Scalars['Float']>,
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly pdbx_isotropic_thermal_model?: Maybe<Scalars['String']>,
  readonly pdbx_ls_cross_valid_method?: Maybe<Scalars['String']>,
  readonly pdbx_ls_sigma_F?: Maybe<Scalars['Float']>,
  readonly pdbx_ls_sigma_Fsqd?: Maybe<Scalars['Float']>,
  readonly pdbx_ls_sigma_I?: Maybe<Scalars['Float']>,
  readonly pdbx_method_to_determine_struct?: Maybe<Scalars['String']>,
  readonly pdbx_overall_ESU_R?: Maybe<Scalars['Float']>,
  readonly pdbx_overall_ESU_R_Free?: Maybe<Scalars['Float']>,
  readonly pdbx_overall_SU_R_Blow_DPI?: Maybe<Scalars['Float']>,
  readonly pdbx_overall_SU_R_free_Blow_DPI?: Maybe<Scalars['Float']>,
  readonly pdbx_overall_SU_R_free_Cruickshank_DPI?: Maybe<Scalars['Float']>,
  readonly pdbx_overall_phase_error?: Maybe<Scalars['Float']>,
  readonly pdbx_refine_id: Scalars['String'],
  readonly pdbx_solvent_ion_probe_radii?: Maybe<Scalars['Float']>,
  readonly pdbx_solvent_shrinkage_radii?: Maybe<Scalars['Float']>,
  readonly pdbx_solvent_vdw_probe_radii?: Maybe<Scalars['Float']>,
  readonly pdbx_starting_model?: Maybe<Scalars['String']>,
  readonly pdbx_stereochem_target_val_spec_case?: Maybe<Scalars['String']>,
  readonly pdbx_stereochemistry_target_values?: Maybe<Scalars['String']>,
  readonly solvent_model_details?: Maybe<Scalars['String']>,
  readonly solvent_model_param_bsol?: Maybe<Scalars['Float']>,
  readonly solvent_model_param_ksol?: Maybe<Scalars['Float']>,
};

export type RefineAnalyze = {
  readonly __typename?: 'RefineAnalyze',
  readonly Luzzati_coordinate_error_free?: Maybe<Scalars['Float']>,
  readonly Luzzati_coordinate_error_obs?: Maybe<Scalars['Float']>,
  readonly Luzzati_d_res_low_free?: Maybe<Scalars['Float']>,
  readonly Luzzati_d_res_low_obs?: Maybe<Scalars['Float']>,
  readonly Luzzati_sigma_a_free?: Maybe<Scalars['Float']>,
  readonly Luzzati_sigma_a_obs?: Maybe<Scalars['Float']>,
  readonly number_disordered_residues?: Maybe<Scalars['Float']>,
  readonly occupancy_sum_hydrogen?: Maybe<Scalars['Float']>,
  readonly occupancy_sum_non_hydrogen?: Maybe<Scalars['Float']>,
  readonly pdbx_Luzzati_d_res_high_obs?: Maybe<Scalars['Float']>,
  readonly pdbx_refine_id: Scalars['String'],
};

export type RefineHist = {
  readonly __typename?: 'RefineHist',
  readonly cycle_id: Scalars['String'],
  readonly d_res_high?: Maybe<Scalars['Float']>,
  readonly d_res_low?: Maybe<Scalars['Float']>,
  readonly number_atoms_solvent?: Maybe<Scalars['Int']>,
  readonly number_atoms_total?: Maybe<Scalars['Int']>,
  readonly pdbx_B_iso_mean_ligand?: Maybe<Scalars['Float']>,
  readonly pdbx_B_iso_mean_solvent?: Maybe<Scalars['Float']>,
  readonly pdbx_number_atoms_ligand?: Maybe<Scalars['Int']>,
  readonly pdbx_number_atoms_nucleic_acid?: Maybe<Scalars['Int']>,
  readonly pdbx_number_atoms_protein?: Maybe<Scalars['Int']>,
  readonly pdbx_number_residues_total?: Maybe<Scalars['Int']>,
  readonly pdbx_refine_id: Scalars['String'],
};

export type RefineLsRestr = {
  readonly __typename?: 'RefineLsRestr',
  readonly dev_ideal?: Maybe<Scalars['Float']>,
  readonly dev_ideal_target?: Maybe<Scalars['Float']>,
  readonly number?: Maybe<Scalars['Int']>,
  readonly pdbx_refine_id: Scalars['String'],
  readonly pdbx_restraint_function?: Maybe<Scalars['String']>,
  readonly type: Scalars['String'],
  readonly weight?: Maybe<Scalars['Float']>,
};

export type Reflns = {
  readonly __typename?: 'Reflns',
  readonly B_iso_Wilson_estimate?: Maybe<Scalars['Float']>,
  readonly R_free_details?: Maybe<Scalars['String']>,
  readonly Rmerge_F_all?: Maybe<Scalars['Float']>,
  readonly Rmerge_F_obs?: Maybe<Scalars['Float']>,
  readonly d_resolution_high?: Maybe<Scalars['Float']>,
  readonly d_resolution_low?: Maybe<Scalars['Float']>,
  readonly data_reduction_details?: Maybe<Scalars['String']>,
  readonly data_reduction_method?: Maybe<Scalars['String']>,
  readonly details?: Maybe<Scalars['String']>,
  readonly limit_h_max?: Maybe<Scalars['Int']>,
  readonly limit_h_min?: Maybe<Scalars['Int']>,
  readonly limit_k_max?: Maybe<Scalars['Int']>,
  readonly limit_k_min?: Maybe<Scalars['Int']>,
  readonly limit_l_max?: Maybe<Scalars['Int']>,
  readonly limit_l_min?: Maybe<Scalars['Int']>,
  readonly number_all?: Maybe<Scalars['Int']>,
  readonly number_obs?: Maybe<Scalars['Int']>,
  readonly observed_criterion?: Maybe<Scalars['String']>,
  readonly observed_criterion_F_max?: Maybe<Scalars['Float']>,
  readonly observed_criterion_F_min?: Maybe<Scalars['Float']>,
  readonly observed_criterion_I_max?: Maybe<Scalars['Float']>,
  readonly observed_criterion_I_min?: Maybe<Scalars['Float']>,
  readonly observed_criterion_sigma_F?: Maybe<Scalars['Float']>,
  readonly observed_criterion_sigma_I?: Maybe<Scalars['Float']>,
  readonly pdbx_CC_half?: Maybe<Scalars['Float']>,
  readonly pdbx_R_split?: Maybe<Scalars['Float']>,
  readonly pdbx_Rmerge_I_obs?: Maybe<Scalars['Float']>,
  readonly pdbx_Rpim_I_all?: Maybe<Scalars['Float']>,
  readonly pdbx_Rrim_I_all?: Maybe<Scalars['Float']>,
  readonly pdbx_Rsym_value?: Maybe<Scalars['Float']>,
  readonly pdbx_chi_squared?: Maybe<Scalars['Float']>,
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly pdbx_netI_over_av_sigmaI?: Maybe<Scalars['Float']>,
  readonly pdbx_netI_over_sigmaI?: Maybe<Scalars['Float']>,
  readonly pdbx_number_measured_all?: Maybe<Scalars['Int']>,
  readonly pdbx_ordinal: Scalars['Int'],
  readonly pdbx_redundancy?: Maybe<Scalars['Float']>,
  readonly pdbx_scaling_rejects?: Maybe<Scalars['Int']>,
  readonly percent_possible_obs?: Maybe<Scalars['Float']>,
  readonly phase_calculation_details?: Maybe<Scalars['String']>,
};

export type ReflnsShell = {
  readonly __typename?: 'ReflnsShell',
  readonly Rmerge_F_all?: Maybe<Scalars['Float']>,
  readonly Rmerge_F_obs?: Maybe<Scalars['Float']>,
  readonly Rmerge_I_all?: Maybe<Scalars['Float']>,
  readonly Rmerge_I_obs?: Maybe<Scalars['Float']>,
  readonly d_res_high?: Maybe<Scalars['Float']>,
  readonly d_res_low?: Maybe<Scalars['Float']>,
  readonly meanI_over_sigI_all?: Maybe<Scalars['Float']>,
  readonly meanI_over_sigI_obs?: Maybe<Scalars['Float']>,
  readonly meanI_over_uI_all?: Maybe<Scalars['Float']>,
  readonly number_measured_all?: Maybe<Scalars['Int']>,
  readonly number_measured_obs?: Maybe<Scalars['Int']>,
  readonly number_possible?: Maybe<Scalars['Int']>,
  readonly number_unique_all?: Maybe<Scalars['Int']>,
  readonly number_unique_obs?: Maybe<Scalars['Int']>,
  readonly pdbx_CC_half?: Maybe<Scalars['Float']>,
  readonly pdbx_R_split?: Maybe<Scalars['Float']>,
  readonly pdbx_Rpim_I_all?: Maybe<Scalars['Float']>,
  readonly pdbx_Rrim_I_all?: Maybe<Scalars['Float']>,
  readonly pdbx_Rsym_value?: Maybe<Scalars['Float']>,
  readonly pdbx_chi_squared?: Maybe<Scalars['Float']>,
  readonly pdbx_diffrn_id?: Maybe<ReadonlyArray<Maybe<Scalars['String']>>>,
  readonly pdbx_netI_over_sigmaI_all?: Maybe<Scalars['Float']>,
  readonly pdbx_netI_over_sigmaI_obs?: Maybe<Scalars['Float']>,
  readonly pdbx_ordinal: Scalars['Int'],
  readonly pdbx_redundancy?: Maybe<Scalars['Float']>,
  readonly pdbx_rejects?: Maybe<Scalars['Int']>,
  readonly percent_possible_all?: Maybe<Scalars['Float']>,
  readonly percent_possible_obs?: Maybe<Scalars['Float']>,
};

export type Software = {
  readonly __typename?: 'Software',
  readonly classification?: Maybe<Scalars['String']>,
  readonly contact_author?: Maybe<Scalars['String']>,
  readonly contact_author_email?: Maybe<Scalars['String']>,
  readonly date?: Maybe<Scalars['String']>,
  readonly description?: Maybe<Scalars['String']>,
  readonly language?: Maybe<Scalars['String']>,
  readonly location?: Maybe<Scalars['String']>,
  readonly name?: Maybe<Scalars['String']>,
  readonly os?: Maybe<Scalars['String']>,
  readonly pdbx_ordinal: Scalars['Int'],
  readonly type?: Maybe<Scalars['String']>,
  readonly version?: Maybe<Scalars['String']>,
};

export type Struct = {
  readonly __typename?: 'Struct',
  readonly pdbx_CASP_flag?: Maybe<Scalars['String']>,
  readonly pdbx_descriptor?: Maybe<Scalars['String']>,
  readonly pdbx_model_details?: Maybe<Scalars['String']>,
  readonly pdbx_model_type_details?: Maybe<Scalars['String']>,
  readonly title?: Maybe<Scalars['String']>,
};

export type StructKeywords = {
  readonly __typename?: 'StructKeywords',
  readonly pdbx_keywords?: Maybe<Scalars['String']>,
  readonly text?: Maybe<Scalars['String']>,
};

export type Symmetry = {
  readonly __typename?: 'Symmetry',
  readonly Int_Tables_number?: Maybe<Scalars['Int']>,
  readonly cell_setting?: Maybe<Scalars['String']>,
  readonly pdbx_full_space_group_name_H_M?: Maybe<Scalars['String']>,
  readonly space_group_name_H_M?: Maybe<Scalars['String']>,
  readonly space_group_name_Hall?: Maybe<Scalars['String']>,
};


export type AssemblySymmetryQueryVariables = {
  assembly_id: Scalars['String'],
  entry_id: Scalars['String']
};


export type AssemblySymmetryQuery = (
  { readonly __typename?: 'Query' }
  & { readonly assembly: Maybe<(
    { readonly __typename?: 'CoreAssembly' }
    & { readonly rcsb_struct_symmetry: Maybe<ReadonlyArray<Maybe<(
      { readonly __typename?: 'RcsbStructSymmetry' }
      & Pick<RcsbStructSymmetry, 'kind' | 'oligomeric_state' | 'stoichiometry' | 'symbol' | 'type'>
      & { readonly clusters: ReadonlyArray<Maybe<(
        { readonly __typename?: 'RcsbStructSymmetryClusters' }
        & Pick<RcsbStructSymmetryClusters, 'avg_rmsd'>
        & { readonly members: ReadonlyArray<Maybe<(
          { readonly __typename?: 'ClustersMembers' }
          & Pick<ClustersMembers, 'asym_id' | 'pdbx_struct_oper_list_ids'>
        )>> }
      )>>, readonly rotation_axes: Maybe<ReadonlyArray<Maybe<(
        { readonly __typename?: 'RcsbStructSymmetryRotationAxes' }
        & Pick<RcsbStructSymmetryRotationAxes, 'order' | 'start' | 'end'>
      )>>> }
    )>>> }
  )> }
);
