/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'CCD' schema file. Dictionary versions: mmCIF 5.326, IHM 1.09, CARB draft.
 *
 * @author molstar/ciftools package
 */

import { Database, Column } from '../../../../mol-data/db';

import Schema = Column.Schema;

const str = Schema.str;
const float = Schema.float;
const List = Schema.List;
const Aliased = Schema.Aliased;
const int = Schema.int;
const coord = Schema.coord;

export const CCD_Schema = {
    /**
     * Data items in the CHEM_COMP category give details about each
     * of the chemical components from which the relevant chemical
     * structures can be constructed, such as name, mass or charge.
     *
     * The related categories CHEM_COMP_ATOM, CHEM_COMP_BOND,
     * CHEM_COMP_ANGLE etc. describe the detailed geometry of these
     * chemical components.
     */
    chem_comp: {
        /**
         * The formula for the chemical component. Formulae are written
         * according to the following rules:
         *
         * (1) Only recognized element symbols may be used.
         *
         * (2) Each element symbol is followed by a 'count' number. A count
         * of '1' may be omitted.
         *
         * (3) A space or parenthesis must separate each cluster of
         * (element symbol + count), but in general parentheses are
         * not used.
         *
         * (4) The order of elements depends on whether carbon is
         * present or not. If carbon is present, the order should be:
         * C, then H, then the other elements in alphabetical order
         * of their symbol. If carbon is not present, the elements
         * are listed purely in alphabetic order of their symbol. This
         * is the 'Hill' system used by Chemical Abstracts.
         */
        formula: str,
        /**
         * Formula mass in daltons of the chemical component.
         */
        formula_weight: float,
        /**
         * The value of _chem_comp.id must uniquely identify each item in
         * the CHEM_COMP list.
         *
         * For protein polymer entities, this is the three-letter code for
         * the amino acid.
         *
         * For nucleic acid polymer entities, this is the one-letter code
         * for the base.
         */
        id: str,
        /**
         * The identifier for the parent component of the nonstandard
         * component. May be be a comma separated list if this component
         * is derived from multiple components.
         *
         * Items in this indirectly point to _chem_comp.id in
         * the CHEM_COMP category.
         */
        mon_nstd_parent_comp_id: List(',', x => x),
        /**
         * The full name of the component.
         */
        name: str,
        /**
         * For standard polymer components, the one-letter code for
         * the component.   For non-standard polymer components, the
         * one-letter code for parent component if this exists;
         * otherwise, the one-letter code should be given as 'X'.
         *
         * Components that derived from multiple parents components
         * are described by a sequence of one-letter-codes.
         */
        one_letter_code: str,
        /**
         * For standard polymer components, the common three-letter code for
         * the component.   Non-standard polymer components and non-polymer
         * components are also assigned three-letter-codes.
         *
         * For ambiguous polymer components three-letter code should
         * be given as 'UNK'.  Ambiguous ions are assigned the code 'UNX'.
         * Ambiguous non-polymer components are assigned the code 'UNL'.
         */
        three_letter_code: str,
        /**
         * For standard polymer components, the type of the monomer.
         * Note that monomers that will form polymers are of three types:
         * linking monomers, monomers with some type of N-terminal (or 5')
         * cap and monomers with some type of C-terminal (or 3') cap.
         */
        type: Aliased<'D-peptide linking' | 'L-peptide linking' | 'D-peptide NH3 amino terminus' | 'L-peptide NH3 amino terminus' | 'D-peptide COOH carboxy terminus' | 'L-peptide COOH carboxy terminus' | 'DNA linking' | 'RNA linking' | 'L-RNA linking' | 'L-DNA linking' | 'DNA OH 5 prime terminus' | 'RNA OH 5 prime terminus' | 'DNA OH 3 prime terminus' | 'RNA OH 3 prime terminus' | 'D-saccharide 1,4 and 1,4 linking' | 'L-saccharide 1,4 and 1,4 linking' | 'D-saccharide 1,4 and 1,6 linking' | 'L-saccharide 1,4 and 1,6 linking' | 'D-saccharide, beta linking' | 'D-saccharide, alpha linking' | 'L-saccharide, beta linking' | 'L-saccharide, alpha linking' | 'L-saccharide' | 'D-saccharide' | 'saccharide' | 'non-polymer' | 'peptide linking' | 'peptide-like' | 'L-gamma-peptide, C-delta linking' | 'D-gamma-peptide, C-delta linking' | 'L-beta-peptide, C-gamma linking' | 'D-beta-peptide, C-gamma linking' | 'other'>(str),
        /**
         * Synonym list for the component.
         */
        pdbx_synonyms: List(';', x => x),
        /**
         * A preliminary classification used by PDB.
         */
        pdbx_type: str,
        /**
         * A preliminary classification used by PDB to indicate
         * that the chemistry of this component while described
         * as clearly as possible is still ambiguous.  Software
         * tools may not be able to process this component
         * definition.
         */
        pdbx_ambiguous_flag: str,
        /**
         * Identifies the _chem_comp.id of the component that
         * has replaced this component.
         */
        pdbx_replaced_by: str,
        /**
         * Identifies the _chem_comp.id's of the components
         * which have been replaced by this component.
         * Multiple id codes should be separated by commas.
         */
        pdbx_replaces: str,
        /**
         * The net integer charge assigned to this component. This is the
         * formal charge assignment normally found in chemical diagrams.
         */
        pdbx_formal_charge: int,
        /**
         * This data item provides additional details about the model coordinates
         * in the component definition.
         */
        pdbx_model_coordinates_details: str,
        /**
         * This data item identifies the PDB database code from which the heavy
         * atom model coordinates were obtained.
         */
        pdbx_model_coordinates_db_code: str,
        /**
         * This data item identifies the source of the ideal coordinates in the
         * component definition.
         */
        pdbx_ideal_coordinates_details: str,
        /**
         * This data item identifies if ideal coordinates are missing in this definition.
         */
        pdbx_ideal_coordinates_missing_flag: Aliased<'Y' | 'N'>(str),
        /**
         * This data item identifies if model coordinates are missing in this definition.
         */
        pdbx_model_coordinates_missing_flag: Aliased<'Y' | 'N'>(str),
        /**
         * Date component was added to database.
         */
        pdbx_initial_date: str,
        /**
         * Date component was last modified.
         */
        pdbx_modified_date: str,
        /**
         * This data item holds the current release status for the component.
         */
        pdbx_release_status: Aliased<'REL' | 'HOLD' | 'HPUB' | 'OBS' | 'DEL' | 'REF_ONLY'>(str),
        /**
         * This data item identifies the deposition site that processed
         * this chemical component defintion.
         */
        pdbx_processing_site: Aliased<'PDBE' | 'EBI' | 'PDBJ' | 'PDBC' | 'RCSB'>(str),
    },
    /**
     * Data items in the CHEM_COMP_ATOM category record details about
     * the atoms in a chemical component. Specifying the atomic
     * coordinates for the components in this category is an
     * alternative to specifying the structure of the component
     * via bonds, angles, planes etc. in the appropriate
     * CHEM_COMP subcategories.
     */
    chem_comp_atom: {
        /**
         * An alternative identifier for the atom. This data item would be
         * used in cases where alternative nomenclatures exist for labelling
         * atoms in a group.
         */
        alt_atom_id: str,
        /**
         * The value of _chem_comp_atom.atom_id must uniquely identify
         * each atom in each monomer in the CHEM_COMP_ATOM list.
         *
         * The atom identifiers need not be unique over all atoms in the
         * data block; they need only be unique for each atom in a
         * component.
         *
         * Note that this item need not be a number; it can be any unique
         * identifier.
         */
        atom_id: str,
        /**
         * The net integer charge assigned to this atom. This is the
         * formal charge assignment normally found in chemical diagrams.
         */
        charge: int,
        /**
         * The x component of the coordinates for this atom in this
         * component specified as orthogonal angstroms. The choice of
         * reference axis frame for the coordinates is arbitrary.
         *
         * The set of coordinates input for the entity here is intended to
         * correspond to the atomic model used to generate restraints for
         * structure refinement, not to atom sites in the ATOM_SITE
         * list.
         */
        model_Cartn_x: coord,
        /**
         * The y component of the coordinates for this atom in this
         * component specified as orthogonal angstroms. The choice of
         * reference axis frame for the coordinates is arbitrary.
         *
         * The set of coordinates input for the entity here is intended to
         * correspond to the atomic model used to generate restraints for
         * structure refinement, not to atom sites in the ATOM_SITE
         * list.
         */
        model_Cartn_y: coord,
        /**
         * The z component of the coordinates for this atom in this
         * component specified as orthogonal angstroms. The choice of
         * reference axis frame for the coordinates is arbitrary.
         *
         * The set of coordinates input for the entity here is intended to
         * correspond to the atomic model used to generate restraints for
         * structure refinement, not to atom sites in the ATOM_SITE
         * list.
         */
        model_Cartn_z: coord,
        /**
         * This data item is a pointer to _chem_comp.id in the CHEM_COMP
         * category.
         */
        comp_id: str,
        /**
         * The code used to identify the atom species representing
         * this atom type. Normally this code is the element
         * symbol.
         */
        type_symbol: str,
        /**
         * Atom name alignment offset in PDB atom field.
         */
        pdbx_align: int,
        /**
         * Ordinal index for the component atom list.
         */
        pdbx_ordinal: int,
        /**
         * An alternative x component of the coordinates for this atom in this
         * component specified as orthogonal angstroms.
         */
        pdbx_model_Cartn_x_ideal: coord,
        /**
         * An alternative y component of the coordinates for this atom in this
         * component specified as orthogonal angstroms.
         */
        pdbx_model_Cartn_y_ideal: coord,
        /**
         * An alternative z component of the coordinates for this atom in this
         * component specified as orthogonal angstroms.
         */
        pdbx_model_Cartn_z_ideal: coord,
        /**
         * The chiral configuration of the atom that is a chiral center.
         */
        pdbx_stereo_config: Aliased<'R' | 'S' | 'N'>(str),
        /**
         * A flag indicating an aromatic atom.
         */
        pdbx_aromatic_flag: Aliased<'Y' | 'N'>(str),
        /**
         * A flag indicating a leaving atom.
         */
        pdbx_leaving_atom_flag: Aliased<'Y' | 'N'>(str),
    },
    /**
     * Data items in the CHEM_COMP_BOND category record details about
     * the bonds between atoms in a chemical component. Target values
     * may be specified as bond orders, as a distance between the two
     * atoms, or both.
     */
    chem_comp_bond: {
        /**
         * The ID of the first of the two atoms that define the bond.
         *
         * This data item is a pointer to _chem_comp_atom.atom_id in the
         * CHEM_COMP_ATOM category.
         */
        atom_id_1: str,
        /**
         * The ID of the second of the two atoms that define the bond.
         *
         * This data item is a pointer to _chem_comp_atom.atom_id in the
         * CHEM_COMP_ATOM category.
         */
        atom_id_2: str,
        /**
         * This data item is a pointer to _chem_comp.id in the CHEM_COMP
         * category.
         */
        comp_id: str,
        /**
         * The value that should be taken as the target for the chemical
         * bond associated with the specified atoms, expressed as a bond
         * order.
         */
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
        /**
         * Ordinal index for the component bond list.
         */
        pdbx_ordinal: int,
        /**
         * Stereochemical configuration across a double bond.
         */
        pdbx_stereo_config: Aliased<'E' | 'Z' | 'N'>(str),
        /**
         * A flag indicating an aromatic bond.
         */
        pdbx_aromatic_flag: Aliased<'Y' | 'N'>(str),
    },
    /**
     * Data items in the CHEM_COMP_DESCRIPTOR category provide
     * string descriptors of component chemical structure.
     */
    pdbx_chem_comp_descriptor: {
        /**
         * This data item is a pointer to _chem_comp.id in the CHEM_COMP
         * category.
         */
        comp_id: str,
        /**
         * This data item contains the descriptor value for this
         * component.
         */
        descriptor: str,
        /**
         * This data item contains the descriptor type.
         */
        type: Aliased<'SMILES_CANNONICAL' | 'SMILES_CANONICAL' | 'SMILES' | 'SMILES' | 'InChI' | 'InChI_MAIN' | 'InChI_MAIN_FORMULA' | 'InChI_MAIN_CONNECT' | 'InChI_MAIN_HATOM' | 'InChI_CHARGE' | 'InChI_STEREO' | 'InChI_ISOTOPE' | 'InChI_FIXEDH' | 'InChI_RECONNECT' | 'InChIKey'>(str),
        /**
         * This data item contains the name of the program
         * or library used to compute the descriptor.
         */
        program: str,
        /**
         * This data item contains the version of the program
         * or library used to compute the descriptor.
         */
        program_version: str,
    },
    /**
     * Data items in the CHEM_COMP_IDENTIFIER category provide
     * identifiers for chemical components.
     */
    pdbx_chem_comp_identifier: {
        /**
         * This data item is a pointer to _chem_comp.id in the CHEM_COMP
         * category.
         */
        comp_id: str,
        /**
         * This data item contains the identifier value for this
         * component.
         */
        identifier: str,
        /**
         * This data item contains the identifier type.
         */
        type: Aliased<'COMMON NAME' | 'SYSTEMATIC NAME' | 'CAS REGISTRY NUMBER' | 'PUBCHEM Identifier' | 'MDL Identifier' | 'SYNONYM' | 'CONDENSED IUPAC CARB SYMBOL' | 'IUPAC CARB SYMBOL' | 'SNFG CARB SYMBOL' | 'CONDENSED IUPAC CARBOHYDRATE SYMBOL' | 'IUPAC CARBOHYDRATE SYMBOL' | 'SNFG CARBOHYDRATE SYMBOL'>(str),
        /**
         * This data item contains the name of the program
         * or library used to compute the identifier.
         */
        program: str,
        /**
         * This data item contains the version of the program
         * or library used to compute the identifier.
         */
        program_version: str,
    },
};

export type CCD_Schema = typeof CCD_Schema;
export interface CCD_Database extends Database<CCD_Schema> {};