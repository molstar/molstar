/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'BIRD' schema file
 *
 * @author mol-star package (src/apps/schema-generator/generate)
 */

import { Database, Column } from 'mol-data/db'

import Schema = Column.Schema

const str = Schema.str;
const int = Schema.int;
const float = Schema.float;
// const coord = Schema.coord;

const Aliased = Schema.Aliased;
// const Matrix = Schema.Matrix;
// const Vector = Schema.Vector;
// const List = Schema.List;

export const BIRD_Schema = {
    pdbx_reference_molecule: {
        prd_id: str,
        formula_weight: float,
        formula: str,
        type: Aliased<'Amino acid' | 'Aminoglycoside' | 'Anthracycline' | 'Anthraquinone' | 'Ansamycin' | 'Chalkophore' | 'Chromophore' | 'Glycopeptide' | 'Cyclic depsipeptide' | 'Cyclic lipopeptide' | 'Cyclic peptide' | 'Heterocyclic' | 'Imino sugar' | 'Keto acid' | 'Lipoglycopeptide' | 'Lipopeptide' | 'Macrolide' | 'Non-polymer' | 'Nucleoside' | 'Oligopeptide' | 'Oligosaccharide' | 'Peptaibol' | 'Peptide-like' | 'Polycyclic' | 'Polypeptide' | 'Polysaccharide' | 'Quinolone' | 'Thiolactone' | 'Thiopeptide' | 'Siderophore' | 'Unknown' | 'Chalkophore, Polypeptide'>(str),
        type_evidence_code: str,
        class: Aliased<'Antagonist' | 'Antibiotic' | 'Anticancer' | 'Anticoagulant' | 'Antifungal' | 'Antiinflammatory' | 'Antimicrobial' | 'Antineoplastic' | 'Antiparasitic' | 'Antiretroviral' | 'Anthelmintic' | 'Antithrombotic' | 'Antitumor' | 'Antiviral' | 'CASPASE inhibitor' | 'Chaperone binding' | 'Enzyme inhibitor' | 'Growth factor' | 'Immunosuppressant' | 'Inhibitor' | 'Lantibiotic' | 'Metabolism' | 'Metal transport' | 'Oxidation-reduction' | 'Receptor' | 'Thrombin inhibitor' | 'Trypsin inhibitor' | 'Toxin' | 'Unknown' | 'Anticoagulant, Antithrombotic' | 'Antibiotic, Antimicrobial' |'Antibiotic, Anthelmintic' | 'Antibiotic, Antineoplastic' | 'Antimicrobial, Antiretroviral' | 'Antimicrobial, Antitumor' | 'Antimicrobial, Antiparasitic, Antibiotic' | 'Thrombin inhibitor, Trypsin inhibitor'>(str),
        class_evidence_code: str,
        name: str,
        represent_as: Aliased<'polymer' | 'single molecule'>(str),
        chem_comp_id: str,
        description: str,
        representative_PDB_id_code: str,
        release_status: Aliased<'REL' | 'HOLD' | 'OBS' | 'WAIT'>(str),
        replaces: str,
        replaced_by: str,
    },
    pdbx_reference_entity_list: {
        prd_id: str,
        ref_entity_id: str,
        type: str,
        details: str,
        component_id: int,
    },
    pdbx_reference_entity_nonpoly: {
        prd_id: str,
        ref_entity_id: str,
        name: str,
        chem_comp_id: str,
    },
    pdbx_reference_entity_link: {
        link_id: int,
        prd_id: str,
        details: str,
        ref_entity_id_1: str,
        ref_entity_id_2: str,
        entity_seq_num_1: int,
        entity_seq_num_2: int,
        comp_id_1: str,
        comp_id_2: str,
        atom_id_1: str,
        atom_id_2: str,
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
        component_1: int,
        component_2: int,
        link_class: Aliased<'PP' | 'PN' | 'NP' | 'NN'>(str),
    },
    pdbx_reference_entity_poly_link: {
        link_id: int,
        prd_id: str,
        ref_entity_id: str,
        component_id: int,
        entity_seq_num_1: int,
        entity_seq_num_2: int,
        comp_id_1: str,
        comp_id_2: str,
        atom_id_1: str,
        atom_id_2: str,
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
    },
    pdbx_reference_entity_poly: {
        prd_id: str,
        ref_entity_id: str,
        type: Aliased<'peptide-like' | 'nucleic-acid-like' | 'polysaccharide-like'>(str),
        db_code: str,
        db_name: str,
    },
    pdbx_reference_entity_poly_seq: {
        prd_id: str,
        ref_entity_id: str,
        mon_id: str,
        parent_mon_id: str,
        num: int,
        observed: Aliased<'Y' | 'N'>(str),
        hetero: Aliased<'Y' | 'N'>(str),
    },
    pdbx_reference_entity_sequence: {
        prd_id: str,
        ref_entity_id: str,
        type: str,
        NRP_flag: Aliased<'Y' | 'N'>(str),
        one_letter_codes: str,
    },
    pdbx_reference_entity_src_nat: {
        prd_id: str,
        ref_entity_id: str,
        ordinal: int,
        organism_scientific: str,
        taxid: str,
        db_code: str,
        db_name: str,
    },
    pdbx_prd_audit: {
        prd_id: str,
        date: str,
        processing_site: str,
        action_type: Aliased<'Initial release' | 'Create molecule' | 'Modify type' | 'Modify class' | 'Modify molecule name' | 'Modify representation' | 'Modify sequence' | 'Modify linkage' | 'Modify taxonomy organism' | 'Modify audit' | 'Other modification' | 'Obsolete molecule'>(str),
    },
}

export type BIRD_Schema = typeof BIRD_Schema;
export type BIRD_Database = Database<BIRD_Schema>