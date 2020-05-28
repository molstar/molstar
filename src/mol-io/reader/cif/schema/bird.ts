/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'BIRD' schema file. Dictionary versions: mmCIF 5.326, IHM 1.09, CARB draft.
 *
 * @author molstar/ciftools package
 */

import { Database, Column } from '../../../../mol-data/db';

import Schema = Column.Schema;

const str = Schema.str;
const float = Schema.float;
const Aliased = Schema.Aliased;
const int = Schema.int;

export const BIRD_Schema = {
    /**
     * Data items in the PDBX_REFERENCE_MOLECULE category record
     * reference information about small polymer molecules.
     */
    pdbx_reference_molecule: {
        /**
         * The value of _pdbx_reference_molecule.prd_id is the unique identifier
         * for the reference molecule in this family.
         *
         * By convention this ID uniquely identifies the reference molecule in
         * in the PDB reference dictionary.
         *
         * The ID has the template form PRD_dddddd (e.g. PRD_000001)
         */
        prd_id: str,
        /**
         * Formula mass in daltons of the entity.
         */
        formula_weight: float,
        /**
         * The formula for the reference entity. Formulae are written
         * according to the rules:
         *
         * 1. Only recognised element symbols may be used.
         *
         * 2. Each element symbol is followed by a 'count' number. A count
         * of '1' may be omitted.
         *
         * 3. A space or parenthesis must separate each element symbol and
         * its count, but in general parentheses are not used.
         *
         * 4. The order of elements depends on whether or not carbon is
         * present. If carbon is present, the order should be: C, then
         * H, then the other elements in alphabetical order of their
         * symbol. If carbon is not present, the elements are listed
         * purely in alphabetic order of their symbol. This is the
         * 'Hill' system used by Chemical Abstracts.
         */
        formula: str,
        /**
         * Defines the structural classification of the entity.
         */
        type: Aliased<'Amino acid' | 'Aminoglycoside' | 'Anthracycline' | 'Anthraquinone' | 'Ansamycin' | 'Chalkophore' | 'Chromophore' | 'Glycopeptide' | 'Cyclic depsipeptide' | 'Cyclic lipopeptide' | 'Cyclic peptide' | 'Heterocyclic' | 'Imino sugar' | 'Keto acid' | 'Lipoglycopeptide' | 'Lipopeptide' | 'Macrolide' | 'Non-polymer' | 'Nucleoside' | 'Oligopeptide' | 'Oligosaccharide' | 'Peptaibol' | 'Peptide-like' | 'Polycyclic' | 'Polypeptide' | 'Polysaccharide' | 'Quinolone' | 'Thiolactone' | 'Thiopeptide' | 'Siderophore' | 'Unknown' | 'Chalkophore, Polypeptide'>(str),
        /**
         * Evidence for the assignment of _pdbx_reference_molecule.type
         */
        type_evidence_code: str,
        /**
         * Broadly defines the function of the entity.
         */
        class: Aliased<'Antagonist' | 'Antibiotic' | 'Anticancer' | 'Anticoagulant' | 'Antifungal' | 'Antigen' | 'Antiinflammatory' | 'Antimicrobial' | 'Antineoplastic' | 'Antiparasitic' | 'Antiretroviral' | 'Anthelmintic' | 'Antithrombotic' | 'Antitumor' | 'Antiviral' | 'CASPASE inhibitor' | 'Chaperone binding' | 'Enzyme inhibitor' | 'Drug delivery' | 'Glycan component' | 'Growth factor' | 'Immunosuppressant' | 'Inducer' | 'Inhibitor' | 'Lantibiotic' | 'Metabolism' | 'Metal transport' | 'Nutrient' | 'Oxidation-reduction' | 'Protein binding' | 'Receptor' | 'Substrate analog' | 'Thrombin inhibitor' | 'Trypsin inhibitor' | 'Toxin' | 'Unknown' | 'Water retention' | 'Anticoagulant, Antithrombotic' | 'Antibiotic, Antimicrobial' | 'Antibiotic, Anthelmintic' | 'Antibiotic, Antineoplastic' | 'Antimicrobial, Antiretroviral' | 'Antimicrobial, Antitumor' | 'Antimicrobial, Antiparasitic, Antibiotic' | 'Thrombin inhibitor, Trypsin inhibitor'>(str),
        /**
         * Evidence for the assignment of _pdbx_reference_molecule.class
         */
        class_evidence_code: str,
        /**
         * A name of the entity.
         */
        name: str,
        /**
         * Defines how this entity is represented in PDB data files.
         */
        represent_as: Aliased<'polymer' | 'single molecule' | 'branched'>(str),
        /**
         * For entities represented as single molecules, the identifier
         * corresponding to the chemical definition for the molecule.
         */
        chem_comp_id: str,
        /**
         * Special details about this molecule.
         */
        compound_details: str,
        /**
         * Description of this molecule.
         */
        description: str,
        /**
         * The PDB accession code for the entry containing a representative example of this molecule.
         */
        representative_PDB_id_code: str,
        /**
         * Defines the current PDB release status for this molecule definition.
         */
        release_status: Aliased<'REL' | 'HOLD' | 'OBS' | 'WAIT'>(str),
        /**
         * Assigns the identifier for the reference molecule which have been replaced
         * by this reference molecule.
         * Multiple molecule identifier codes should be separated by commas.
         */
        replaces: str,
        /**
         * Assigns the identifier of the reference molecule that has replaced this molecule.
         */
        replaced_by: str,
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_LIST category record
     * the list of entities within each reference molecule.
     */
    pdbx_reference_entity_list: {
        /**
         * The value of _pdbx_reference_entity_list.prd_id is a reference
         * _pdbx_reference_molecule.prd_id in the PDBX_REFERENCE_MOLECULE category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_list.ref_entity_id is a unique identifier
         * the a constituent entity within this reference molecule.
         */
        ref_entity_id: str,
        /**
         * Defines the polymer characteristic of the entity.
         */
        type: Aliased<'polymer' | 'polymer-like' | 'non-polymer' | 'branched'>(str),
        /**
         * Additional details about this entity.
         */
        details: str,
        /**
         * The component number of this entity within the molecule.
         */
        component_id: int,
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_NONPOLY category record
     * the list of entities within each reference molecule.
     */
    pdbx_reference_entity_nonpoly: {
        /**
         * The value of _pdbx_reference_entity_nonpoly.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the PDBX_REFERENCE_ENTITY_LIST category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_nonpoly.ref_entity_id is a reference
         * to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id: str,
        /**
         * A name of the non-polymer entity.
         */
        name: str,
        /**
         * For non-polymer entities, the identifier corresponding
         * to the chemical definition for the molecule.
         */
        chem_comp_id: str,
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_LINK category give details about
     * the linkages between entities within reference molecules.
     */
    pdbx_reference_entity_link: {
        /**
         * The value of _pdbx_reference_entity_link.link_id uniquely identifies
         * linkages between entities with a molecule.
         */
        link_id: int,
        /**
         * The value of _pdbx_reference_entity_link.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the PDBX_REFERENCE_ENTITY_LIST category.
         */
        prd_id: str,
        /**
         * A description of special aspects of a linkage between
         * chemical components in the structure.
         */
        details: str,
        /**
         * The reference entity id of the first of the two entities joined by the
         * linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_list.ref_entity_id
         * in the PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id_1: str,
        /**
         * The reference entity id of the second of the two entities joined by the
         * linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_list.ref_entity_id
         * in the PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id_2: str,
        /**
         * For a polymer entity, the sequence number in the first of
         * the two entities containing the linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.num
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        entity_seq_num_1: int,
        /**
         * For a polymer entity, the sequence number in the second of
         * the two entities containing the linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.num
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        entity_seq_num_2: int,
        /**
         * The component identifier in the first of the two entities containing the linkage.
         *
         * For polymer entities, this data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         *
         * For non-polymer entities, this data item is a pointer to
         * _pdbx_reference_entity_nonpoly.chem_comp_id in the
         * PDBX_REFERENCE_ENTITY_NONPOLY category.
         */
        comp_id_1: str,
        /**
         * The component identifier in the second of the two entities containing the linkage.
         *
         * For polymer entities, this data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         *
         * For non-polymer entities, this data item is a pointer to
         * _pdbx_reference_entity_nonpoly.chem_comp_id in the
         * PDBX_REFERENCE_ENTITY_NONPOLY category.
         */
        comp_id_2: str,
        /**
         * The atom identifier/name in the first of the two entities containing the linkage.
         */
        atom_id_1: str,
        /**
         * The atom identifier/name in the second of the two entities containing the linkage.
         */
        atom_id_2: str,
        /**
         * The bond order target for the chemical linkage.
         */
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
        /**
         * The entity component identifier for the first of two entities containing the linkage.
         */
        component_1: int,
        /**
         * The entity component identifier for the second of two entities containing the linkage.
         */
        component_2: int,
        /**
         * A code indicating the entity types involved in the linkage.
         */
        link_class: Aliased<'PP' | 'PN' | 'NP' | 'NN'>(str),
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_POLY_LINK category give details about
     * polymer linkages including both standard and non-standard linkages between
     * polymer componnents.
     */
    pdbx_reference_entity_poly_link: {
        /**
         * The value of _pdbx_reference_entity_poly_link.link_id uniquely identifies
         * a linkage within a polymer entity.
         */
        link_id: int,
        /**
         * The value of _pdbx_reference_entity_poly_link.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the PDBX_REFERENCE_ENTITY_POLY category.
         */
        prd_id: str,
        /**
         * The reference entity id of the polymer entity containing the linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly.ref_entity_id
         * in the PDBX_REFERENCE_ENTITY_POLY category.
         */
        ref_entity_id: str,
        /**
         * The entity component identifier entity containing the linkage.
         */
        component_id: int,
        /**
         * For a polymer entity, the sequence number in the first of
         * the two components making the linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.num
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        entity_seq_num_1: int,
        /**
         * For a polymer entity, the sequence number in the second of
         * the two components making the linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.num
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        entity_seq_num_2: int,
        /**
         * The component identifier in the first of the two components making the
         * linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        comp_id_1: str,
        /**
         * The component identifier in the second of the two components making the
         * linkage.
         *
         * This data item is a pointer to _pdbx_reference_entity_poly_seq.mon_id
         * in the PDBX_REFERENCE_ENTITY_POLY_SEQ category.
         */
        comp_id_2: str,
        /**
         * The atom identifier/name in the first of the two components making
         * the linkage.
         */
        atom_id_1: str,
        /**
         * The atom identifier/name in the second of the two components making
         * the linkage.
         */
        atom_id_2: str,
        /**
         * The bond order target for the non-standard linkage.
         */
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_POLY category record details about
     * the polymer, such as the type of the polymer, the number of
     * monomers and whether it has nonstandard features.
     */
    pdbx_reference_entity_poly: {
        /**
         * The value of _pdbx_reference_entity_poly.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_poly.ref_entity_id is a reference
         * to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id: str,
        /**
         * The type of the polymer.
         */
        type: Aliased<'peptide-like' | 'nucleic-acid-like' | 'polysaccharide-like' | 'oligosaccharide'>(str),
        /**
         * The database code for this source information
         */
        db_code: str,
        /**
         * The database name for this source information
         */
        db_name: str,
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_POLY_SEQ category specify the sequence
     * of monomers in a polymer.
     */
    pdbx_reference_entity_poly_seq: {
        /**
         * The value of _pdbx_reference_entity_poly_seq.prd_id is a reference
         * _pdbx_reference_entity_poly.prd_id in the  PDBX_REFERENCE_ENTITY_POLY category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_poly_seq.ref_entity_id is a reference
         * to _pdbx_reference_entity_poly.ref_entity_id in PDBX_REFERENCE_ENTITY_POLY category.
         */
        ref_entity_id: str,
        /**
         * This data item is the chemical component identifier of monomer.
         */
        mon_id: str,
        /**
         * This data item is the chemical component identifier for the parent component corresponding to this monomer.
         */
        parent_mon_id: str,
        /**
         * The value of _pdbx_reference_entity_poly_seq.num must uniquely and sequentially
         * identify a record in the PDBX_REFERENCE_ENTITY_POLY_SEQ list.
         *
         * This value is conforms to author numbering conventions and does not map directly
         * to the numbering conventions used for _entity_poly_seq.num.
         */
        num: int,
        /**
         * A flag to indicate that this monomer is observed in the instance example.
         */
        observed: Aliased<'Y' | 'N'>(str),
        /**
         * A flag to indicate that sequence heterogeneity at this monomer position.
         */
        hetero: Aliased<'Y' | 'N'>(str),
    },
    /**
     * Additional features associated with the reference entity.
     */
    pdbx_reference_entity_sequence: {
        /**
         * The value of _pdbx_reference_entity_sequence.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_sequence.ref_entity_id is a reference
         * to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id: str,
        /**
         * The monomer type for the sequence.
         */
        type: Aliased<'peptide-like' | 'saccharide'>(str),
        /**
         * A flag to indicate a non-ribosomal entity.
         */
        NRP_flag: Aliased<'Y' | 'N'>(str),
        /**
         * The one-letter-code sequence for this entity.  Non-standard monomers are represented as 'X'.
         */
        one_letter_codes: str,
    },
    /**
     * Data items in the PDBX_REFERENCE_ENTITY_SRC_NAT category record
     * details of the source from which the entity was obtained.
     */
    pdbx_reference_entity_src_nat: {
        /**
         * The value of _pdbx_reference_entity_src_nat.prd_id is a reference
         * _pdbx_reference_entity_list.prd_id in the  PDBX_REFERENCE_ENTITY_LIST category.
         */
        prd_id: str,
        /**
         * The value of _pdbx_reference_entity_src_nat.ref_entity_id is a reference
         * to _pdbx_reference_entity_list.ref_entity_id in PDBX_REFERENCE_ENTITY_LIST category.
         */
        ref_entity_id: str,
        /**
         * The value of _pdbx_reference_entity_src_nat.ordinal distinguishes
         * source details for this entity.
         */
        ordinal: int,
        /**
         * The scientific name of the organism from which the entity was isolated.
         */
        organism_scientific: str,
        /**
         * The NCBI TaxId of the organism from which the entity was isolated.
         */
        taxid: str,
        /**
         * The database code for this source information
         */
        db_code: str,
        /**
         * The database name for this source information
         */
        db_name: str,
    },
    /**
     * Data items in the PDBX_PRD_AUDIT category records
     * the status and tracking information for this molecule.
     */
    pdbx_prd_audit: {
        /**
         * This data item is a pointer to _pdbx_reference_molecule.prd_id in the
         * pdbx_reference_molecule category.
         */
        prd_id: str,
        /**
         * The date associated with this audit record.
         */
        date: str,
        /**
         * An identifier for the wwPDB site creating or modifying the molecule.
         */
        processing_site: Aliased<'RCSB' | 'PDBe' | 'PDBJ' | 'BMRB' | 'PDBC'>(str),
        /**
         * The action associated with this audit record.
         */
        action_type: Aliased<'Initial release' | 'Create molecule' | 'Modify type' | 'Modify class' | 'Modify molecule name' | 'Modify representation' | 'Modify sequence' | 'Modify linkage' | 'Modify taxonomy organism' | 'Modify audit' | 'Other modification' | 'Obsolete molecule'>(str),
    },
};

export type BIRD_Schema = typeof BIRD_Schema;
export interface BIRD_Database extends Database<BIRD_Schema> {};