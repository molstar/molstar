/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'CifCore' schema file. Dictionary versions: CifCore 3.3.0.
 *
 * @author molstar/ciftools package
 */

import { Database, Column } from '../../../../mol-data/db';

import Schema = Column.Schema;

const float = Schema.float;
const int = Schema.int;
const str = Schema.str;
const Matrix = Schema.Matrix;

export const CifCore_Schema = {
    /**
     * The category of data items used to describe the parameters of
     * the crystal unit cell.
     */
    cell: {
        /**
         * The angle between the bounding cell axes.
         */
        angle_alpha: float,
        /**
         * The angle between the bounding cell axes.
         */
        angle_beta: float,
        /**
         * The angle between the bounding cell axes.
         */
        angle_gamma: float,
        /**
         * The number of the formula units in the unit cell as specified
         * by _chemical_formula.structural, _chemical_formula.moiety or
         * _chemical_formula.sum.
         */
        formula_units_z: int,
        /**
         * The length of each cell axis.
         */
        length_a: float,
        /**
         * The length of each cell axis.
         */
        length_b: float,
        /**
         * The length of each cell axis.
         */
        length_c: float,
        /**
         * Volume of the crystal unit cell.
         */
        volume: float,
    },
    /**
     * The category of data items which describe the composition and
     * chemical properties of the compound under study. The formula data
     * items must be consistent with the density, unit-cell and Z values.
     */
    chemical: {
        /**
         * The temperature at which a crystalline solid changes to a liquid.
         */
        melting_point: float,
        /**
         * Trivial name by which the compound is commonly known.
         */
        name_common: str,
        /**
         * IUPAC or Chemical Abstracts full name of compound.
         */
        name_systematic: str,
    },
    /**
     * The category of data items which specify the composition and chemical
     * properties of the compound. The formula data items must agree
     * with those that specify the density, unit-cell and Z values.
     *
     * The following rules apply to the construction of the data items
     * _chemical_formula.analytical, *.structural and *.sum. For the
     * data item *.moiety the formula construction is broken up into
     * residues or moieties, i.e. groups of atoms that form a molecular
     * unit or molecular ion. The rules given below apply within each
     * moiety, but different requirements apply to the way that moieties
     * are connected (see _chemical_formula.moiety).
     *
     * 1. Only recognized element symbols may be used.
     *
     * 2. Each element symbol is followed by a 'count' number. A count of
     * '1' may be omitted.
     *
     * 3. A space or parenthesis must separate each cluster of (element
     * symbol + count).
     *
     * 4. Where a group of elements is enclosed in parentheses, the
     * multiplier for the group must follow the closing parentheses.
     * That is, all element and group multipliers are assumed to be
     * printed as subscripted numbers. [An exception to this rule
     * exists for *.moiety formulae where pre- and post-multipliers
     * are permitted for molecular units].
     *
     * 5. Unless the elements are ordered in a manner that corresponds to
     * their chemical structure, as in _chemical_formula.structural,
     * the order of the elements within any group or moiety
     * depends on whether or not carbon is present. If carbon is
     * present, the order should be: C, then H, then the other
     * elements in alphabetical order of their symbol. If carbon is
     * not present, the elements are listed purely in alphabetical order
     * of their symbol. This is the 'Hill' system used by Chemical
     * Abstracts. This ordering is used in _chemical_formula.moiety
     * and _chemical_formula.sum.
     *
     * _chemical_formula.IUPAC      '[Mo (C O)4 (C18 H33 P)2]'
     * _chemical_formula.moiety     'C40 H66 Mo O4 P2'
     * _chemical_formula.structural '((C O)4 (P (C6 H11)3)2)Mo'
     * _chemical_formula.sum         'C40 H66 Mo O4 P2'
     * _chemical_formula.weight      768.81
     */
    chemical_formula: {
        /**
         * Formula with each discrete bonded residue or ion shown as a
         * separate moiety. See above CHEMICAL_FORMULA for rules
         * for writing chemical formulae. In addition to the general
         * formulae requirements, the following rules apply:
         * 1. Moieties are separated by commas ','.
         * 2. The order of elements within a moiety follows general rule
         * 5 in CHEMICAL_FORMULA.
         * 3. Parentheses are not used within moieties but may surround
         * a moiety. Parentheses may not be nested.
         * 4. Charges should be placed at the end of the moiety. The
         * charge '+' or '-' may be preceded by a numerical multiplier
         * and should be separated from the last (element symbol +
         * count) by a space. Pre- or post-multipliers may be used for
         * individual moieties.
         */
        moiety: str,
        /**
         * Chemical formulae in which all discrete bonded residues and ions are
         * summed over the constituent elements, following the ordering given
         * in rule 5 of the category description. Parentheses normally not used.
         */
        sum: str,
        /**
         * Mass corresponding to the formulae _chemical_formula.structural,
         * *_IUPAC, *_moiety or *_sum and, together with the Z value and cell
         * parameters yield the density given as _exptl_crystal.density_diffrn.
         */
        weight: float,
    },
    /**
     * The category of data items used to specify space group
     * information about the crystal used in the diffraction measurements.
     *
     * Space-group types are identified by their number as listed in
     * International Tables for Crystallography Volume A, or by their
     * Schoenflies symbol. Specific settings of the space groups can
     * be identified by their Hall symbol, by specifying their
     * symmetry operations or generators, or by giving the
     * transformation that relates the specific setting to the
     * reference setting based on International Tables Volume A and
     * stored in this dictionary.
     *
     * The commonly used Hermann-Mauguin symbol determines the
     * space-group type uniquely, but several different Hermann-Mauguin
     * symbols may refer to the same space-group type. A
     * Hermann-Mauguin symbol contains information on the choice of
     * the basis, but not on the choice of origin.
     *
     * Ref: International Tables for Crystallography (2002). Volume A,
     * Space-group symmetry, edited by Th. Hahn, 5th ed.
     * Dordrecht: Kluwer Academic Publishers.
     */
    space_group: {
        /**
         * The name of the system of geometric crystal classes of space
         * groups (crystal system) to which the space group belongs.
         * Note that rhombohedral space groups belong to the
         * trigonal system.
         */
        crystal_system: str,
        /**
         * The number as assigned in International Tables for Crystallography
         * Vol. A, specifying the proper affine class (i.e. the orientation
         * preserving affine class) of space groups (crystallographic space
         * group type) to which the space group belongs. This number defines
         * the space group type but not the coordinate system expressed.
         */
        it_number: int,
        /**
         * The full international Hermann-Mauguin space-group symbol as
         * defined in Section 2.2.3 and given as the second item of the
         * second line of each of the space-group tables of Part 7 of
         * International Tables for Crystallography Volume A (2002).
         *
         * Each component of the space-group name is separated by a
         * space or an underscore character. The use of a space is
         * strongly recommended.  The underscore is only retained
         * because it was used in old CIFs. It should not be used in
         * new CIFs.
         *
         * Subscripts should appear without special symbols. Bars should
         * be given as negative signs before the numbers to which they
         * apply. The commonly used Hermann-Mauguin symbol determines the
         * space-group type uniquely, but a given space-group type may
         * be described by more than one Hermann-Mauguin symbol. The
         * space-group type is best described using
         * _space_group.IT_number or _space_group.name_Schoenflies. The
         * full international Hermann-Mauguin symbol contains information
         * about the choice of basis for monoclinic and orthorhombic
         * space groups, but does not give information about the choice
         * of origin. To define the setting uniquely use
         * _space_group.name_Hall, or list the symmetry operations
         * or generators.
         *
         * Ref: International Tables for Crystallography (2002). Volume A,
         * Space-group symmetry, edited by Th. Hahn, 5th ed.
         * Dordrecht: Kluwer Academic Publishers.
         */
        'name_h-m_full': str,
    },
    /**
     * The category of data items used to describe symmetry equivalent sites
     * in the crystal unit cell.
     */
    space_group_symop: {
        /**
         * A parsable string giving one of the symmetry operations of the
         * space group in algebraic form.  If W is a matrix representation
         * of the rotational part of the symmetry operation defined by the
         * positions and signs of x, y and z, and w is a column of
         * translations defined by fractions, an equivalent position
         * X' is generated from a given position X by the equation
         *
         * X' = WX + w.
         *
         * (Note: X is used to represent bold_italics_x in International
         * Tables for Crystallography Vol. A, Part 5.)
         *
         * When a list of symmetry operations is given, it must contain
         * a complete set of coordinate representatives which generates
         * all the operations of the space group by the addition of
         * all primitive translations of the space group. Such
         * representatives are to be found as the coordinates of
         * the general-equivalent position in International Tables for
         * Crystallography Vol. A (2002), to which it is necessary to
         * add any centring translations shown above the
         * general-equivalent position.
         *
         * That is to say, it is necessary to list explicitly all the
         * symmetry operations required to generate all the atoms in
         * the unit cell defined by the setting used.
         */
        operation_xyz: str,
    },
    /**
     * The category of data items used to specify the geometry bonds in the
     * structural model as derived from the atomic sites.
     */
    geom_bond: {
        /**
         * This label is a unique identifier for a particular site in the
         * asymmetric unit of the crystal unit cell.
         */
        atom_site_label_1: str,
        /**
         * This label is a unique identifier for a particular site in the
         * asymmetric unit of the crystal unit cell.
         */
        atom_site_label_2: str,
        /**
         * Intramolecular bond distance between the sites identified
         * by _geom_bond.id
         */
        distance: float,
        /**
         * This code signals whether the angle is referred to in a
         * publication or should be placed in a table of significant angles.
         */
        publ_flag: str,
        /**
         * Data item specifying the symmetry operation codes applied to the atom
         * sites involved in a specific geometric configuration or other correlated
         * behaviour.
         *
         * The symmetry code of each atom site as the symmetry-equivalent position
         * number 'n' and the cell translation number 'pqr'. These numbers are
         * combined to form the code 'n pqr' or n_pqr.
         *
         * The character string n_pqr is composed as follows:
         *
         * n refers to the symmetry operation that is applied to the
         * coordinates stored in _atom_site.fract_xyz. It must match
         * a number given in _space_group_symop.id (or one of its
         * aliases, such as _symmetry_equiv_pos_site_id).
         *
         * p, q and r refer to the translations that are subsequently
         * applied to the symmetry transformed coordinates to generate
         * the related atom position. These translations (x,y,z) are related
         * to (p,q,r) by the relations
         * p = 5 + x
         * q = 5 + y
         * r = 5 + z
         */
        site_symmetry_1: str,
        /**
         * Data item specifying the symmetry operation codes applied to the atom
         * sites involved in a specific geometric configuration or other correlated
         * behaviour.
         *
         * The symmetry code of each atom site as the symmetry-equivalent position
         * number 'n' and the cell translation number 'pqr'. These numbers are
         * combined to form the code 'n pqr' or n_pqr.
         *
         * The character string n_pqr is composed as follows:
         *
         * n refers to the symmetry operation that is applied to the
         * coordinates stored in _atom_site.fract_xyz. It must match
         * a number given in _space_group_symop.id (or one of its
         * aliases, such as _symmetry_equiv_pos_site_id).
         *
         * p, q and r refer to the translations that are subsequently
         * applied to the symmetry transformed coordinates to generate
         * the related atom position. These translations (x,y,z) are related
         * to (p,q,r) by the relations
         * p = 5 + x
         * q = 5 + y
         * r = 5 + z
         */
        site_symmetry_2: str,
        /**
         * Bond valence calculated from the bond distance.
         */
        valence: float,
    },
    /**
     * The category of data items used to record details about the
     * creation and subsequent updating of the data block.
     */
    audit: {
        /**
         * The digital object identifier (DOI) registered to identify
         * the data set publication represented by the current
         * data block. This can be used as a unique identifier for
         * the data block so long as the code used is a valid DOI
         * (i.e. begins with a valid publisher prefix assigned by a
         * Registration Agency and a suffix guaranteed to be unique
         * by the publisher) and has had its metadata deposited
         * with a DOI Registration Agency.
         *
         * A DOI is a unique character string identifying any
         * object of intellectual property. It provides a
         * persistent identifier for an object on a digital network
         * and permits the association of related current data in a
         * structured extensible way. A DOI is an implementation
         * of the Internet concepts of Uniform Resource Name and
         * Universal Resource Locator managed according to the
         * specifications of the International DOI Foundation
         * (see https://www.doi.org/).
         */
        block_doi: str,
    },
    /**
     * The category of data items recording database deposition. These data items
     * are assigned by database managers and should only appear in a CIF if they
     * originate from that source.
     */
    database_code: {
        /**
         * Code assigned by the Crystallography Open Database (COD).
         */
        cod: str,
        /**
         * Code assigned by the Cambridge Structural Database.
         */
        csd: str,
        /**
         * Deposition numbers assigned by the Cambridge Crystallographic
         * Data Centre (CCDC) to files containing structural information
         * archived by the CCDC.
         */
        depnum_ccdc_archive: str,
        /**
         * Deposition numbers assigned by the Fachinformationszentrum
         * Karlsruhe (FIZ) to files containing structural information
         * archived by the Cambridge Crystallographic Data Centre (CCDC).
         */
        depnum_ccdc_fiz: str,
        /**
         * Code assigned by the Inorganic Crystal Structure Database.
         */
        icsd: str,
        /**
         * Code assigned in the Metals Data File.
         */
        mdf: str,
        /**
         * Code assigned by the NBS (NIST) Crystal Data Database.
         */
        nbs: str,
    },
    /**
     * The category of data items used to describe atom site information
     * used in crystallographic structure studies.
     */
    atom_site: {
        /**
         * Code for type of atomic displacement parameters used for the site.
         */
        adp_type: str,
        /**
         * A standard code to signal if the site coordinates have been
         * determined from the intensities or calculated from the geometry
         * of surrounding sites, or have been assigned dummy coordinates.
         */
        calc_flag: str,
        /**
         * A code which identifies a cluster of atoms that show long range disorder
         * but are locally ordered. Within each such cluster of atoms,
         * _atom_site.disorder_group is used to identify the sites that are
         * simultaneously occupied. This field is only needed if there is more than
         * one cluster of disordered atoms showing independent local order.
         */
        disorder_assembly: str,
        /**
         * A code that identifies a group of disordered atom sites that are locally
         * simultaneously occupied. Atoms that are positionally disordered over two or
         * more sites (e.g. the H atoms of a methyl group that exists in two
         * orientations) should be assigned to two or more groups. Similarly, atoms
         * that describe a specific alternative composition of a compositionally
         * disordered site should be assigned to a distinct disorder group (e.g. a site
         * that is partially occupied by Mg and Mn atoms should be described by
         * assigning the Mg atom to one group and the Mn atom to another group). Sites
         * belonging to the same group are simultaneously occupied, but those belonging
         * to different groups are not. A minus prefix (e.g. "-1") is used to indicate
         * sites disordered about a special position.
         */
        disorder_group: str,
        /**
         * Atom site coordinates as fractions of the cell length values.
         */
        fract_x: float,
        /**
         * Atom site coordinates as fractions of the cell length values.
         */
        fract_y: float,
        /**
         * Atom site coordinates as fractions of the cell length values.
         */
        fract_z: float,
        /**
         * This label is a unique identifier for a particular site in the
         * asymmetric unit of the crystal unit cell. It is made up of
         * components, _atom_site.label_component_0 to *_6, which may be
         * specified as separate data items. Component 0 usually matches one
         * of the specified _atom_type.symbol codes. This is not mandatory
         * if an _atom_site.type_symbol item is included in the atom site
         * list. The _atom_site.type_symbol always takes precedence over
         * an _atom_site.label in the identification of the atom type. The
         * label components 1 to 6 are optional, and normally only
         * components 0 and 1 are used. Note that components 0 and 1 are
         * concatenated, while all other components, if specified, are
         * separated by an underline character. Underline separators are
         * only used if higher-order components exist. If an intermediate
         * component is not used it may be omitted provided the underline
         * separators are inserted. For example the label 'C233__ggg' is
         * acceptable and represents the components C, 233, '', and ggg.
         * Each label may have a different number of components.
         */
        label: str,
        /**
         * The fraction of the atom type present at this site.
         * The sum of the occupancies of all the atom types at this site
         * may not significantly exceed 1.0 unless it is a dummy site. The
         * value must lie in the 99.97% Gaussian confidence interval
         * -3u =< x =< 1 + 3u. The _enumeration.range of 0.0:1.0 is thus
         * correctly interpreted as meaning (0.0 - 3u) =< x =< (1.0 + 3u).
         */
        occupancy: float,
        /**
         * A concatenated series of single-letter codes which indicate the
         * refinement restraints or constraints applied to this site. This
         * item should not be used. It has been replaced by
         * _atom_site.refinement_flags_posn, _ADP and _occupancy. It is
         * retained in this dictionary only to provide compatibility with
         * legacy CIFs.
         */
        refinement_flags: str,
        /**
         * The number of different sites that are generated by the
         * application of the space-group symmetry to the coordinates
         * given for this site. It is equal to the multiplicity given
         * for this Wyckoff site in International Tables for Cryst.
         * Vol. A (2002). It is equal to the multiplicity of the general
         * position divided by the order of the site symmetry given in
         * _atom_site.site_symmetry_order.
         *
         * The _atom_site_symmetry_multiplicity form of this data name is
         * deprecated because of historical inconsistencies in practice among
         * structure refinement software packages and should not be used.
         */
        site_symmetry_multiplicity: int,
        /**
         * A code to identify the atom specie(s) occupying this site.
         * This code must match a corresponding _atom_type.symbol. The
         * specification of this code is optional if component_0 of the
         * _atom_site.label is used for this purpose. See _atom_type.symbol.
         */
        type_symbol: str,
        /**
         * Isotropic atomic displacement parameter, or equivalent isotropic
         * atomic displacement parameter, U(equiv), in angstroms squared,
         * calculated from anisotropic atomic displacement parameters.
         *
         * U(equiv) = (1/3) sum~i~[sum~j~(U^ij^ a*~i~ a*~j~ a~i~.a~j~)]
         *
         * a  = the real-space cell vectors
         * a* = the reciprocal-space cell lengths
         * Ref: Fischer, R. X. & Tillmanns, E. (1988). Acta Cryst. C44, 775-776.
         */
        u_iso_or_equiv: float,
    },
    /**
     * The category of data items used to describe the anisotropic atomic
     * displacement parameters of the atomic sites in a crystal structure.
     */
    atom_site_aniso: {
        /**
         * Anisotropic atomic displacement parameters are usually looped in
         * a separate list. If this is the case, this code must match the
         * _atom_site.label of the associated atom in the atom coordinate
         * list and conform with the same rules described in _atom_site.label.
         */
        label: str,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_11: float,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u: Matrix(3, 3),
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_12: float,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_13: float,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_22: float,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_23: float,
        /**
         * These are the standard anisotropic atomic displacement components in
         * angstroms squared which appear in the structure-factor term
         *
         * T = exp{ -2π^2^ sum~i~ [ sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }
         *
         * h = the Miller indices
         * a* = the reciprocal-space cell lengths
         *
         * The unique elements of the real symmetric matrix are entered by row.
         *
         * The IUCr Commission on Nomenclature recommends the use of U for reporting
         * atomic displacement parameters.
         *
         * Note that U^ij^ = β^ij^/(2 π^2^ a*~i~ a*~j~) = B^ij^/(8 π^2^).
         *
         * Ref: Trueblood, K. N. et al. (1996). Acta Cryst. A52, 770-781.
         */
        u_33: float,
    },
    /**
     * The category of data items used to describe atomic type information
     * used in crystallographic structure studies.
     */
    atom_type: {
        /**
         * A description of the atom(s) designated by this atom type. In
         * most cases this will be the element name and oxidation state of
         * a single atom species. For disordered or nonstoichiometric
         * structures it will describe a combination of atom species.
         */
        description: str,
        /**
         * The identity of the atom specie(s) representing this atom type.
         * Normally this code is the element symbol followed by the charge
         * if there is one. The symbol may be composed of any character except
         * an underline or a blank, with the proviso that digits designate an
         * oxidation state and must be followed by a + or - character.
         */
        symbol: str,
    },
    /**
     * The category of data items used to describe atomic scattering
     * information used in crystallographic structure studies.
     */
    atom_type_scat: {
        /**
         * The imaginary component of the anomalous dispersion scattering factors
         * for this atom type and radiation by _diffrn_radiation_wavelength.value
         */
        dispersion_imag: float,
        /**
         * The real component of the anomalous dispersion scattering factors
         * for this atom type and radiation by _diffrn_radiation_wavelength.value
         */
        dispersion_real: float,
        /**
         * Reference to source of scattering factors used for this atom type.
         */
        source: str,
    },
};

export const CifCore_Aliases = {
    'cell.formula_units_z': [
        'cell_formula_units_Z',
    ],
    'space_group.it_number': [
        'space_group_IT_number',
        'symmetry_Int_Tables_number',
    ],
    'space_group.name_h-m_full': [
        'symmetry_space_group_name_H-M',
    ],
    'space_group_symop.operation_xyz': [
        'symmetry_equiv_pos_as_xyz',
    ],
    'geom_bond.atom_site_label_1': [
        'geom_bond_atom_site_id_1',
    ],
    'geom_bond.atom_site_label_2': [
        'geom_bond_atom_site_id_2',
    ],
    'geom_bond.distance': [
        'geom_bond_dist',
    ],
    'audit.block_doi': [
        'audit_block_DOI',
    ],
    'database_code.cod': [
        'database_code_COD',
    ],
    'database_code.csd': [
        'database_code_CSD',
    ],
    'database_code.depnum_ccdc_archive': [
        'database_code_depnum_CCDC_archive',
    ],
    'database_code.depnum_ccdc_fiz': [
        'database_code_depnum_CCDC_fiz',
    ],
    'database_code.icsd': [
        'database_code_ICSD',
    ],
    'database_code.mdf': [
        'database_code_MDF',
    ],
    'database_code.nbs': [
        'database_code_NBS',
    ],
    'atom_site.adp_type': [
        'atom_site_ADP_type',
        'atom_site_thermal_displace_type',
    ],
    'atom_site.label': [
        'atom_site_id',
    ],
    'atom_site.site_symmetry_multiplicity': [
        'atom_site_symmetry_multiplicity',
    ],
    'atom_site.u_iso_or_equiv': [
        'atom_site_U_iso_or_equiv',
    ],
    'atom_site_aniso.label': [
        'atom_site_anisotrop_id',
    ],
    'atom_site_aniso.u_11': [
        'atom_site_aniso_U_11',
        'atom_site_anisotrop_U_11',
    ],
    'atom_site_aniso.u_12': [
        'atom_site_aniso_U_12',
        'atom_site_anisotrop_U_12',
    ],
    'atom_site_aniso.u_13': [
        'atom_site_aniso_U_13',
        'atom_site_anisotrop_U_13',
    ],
    'atom_site_aniso.u_22': [
        'atom_site_aniso_U_22',
        'atom_site_anisotrop_U_22',
    ],
    'atom_site_aniso.u_23': [
        'atom_site_aniso_U_23',
        'atom_site_anisotrop_U_23',
    ],
    'atom_site_aniso.u_33': [
        'atom_site_aniso_U_33',
        'atom_site_anisotrop_U_33',
    ],
};

export type CifCore_Schema = typeof CifCore_Schema;
export interface CifCore_Database extends Database<CifCore_Schema> {};