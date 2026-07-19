/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'SF' schema file. Dictionary versions: mmCIF 5.416, IHM 1.28, MA 1.4.9.
 *
 * @author molstar/ciftools package
 */

import { Database, Column } from '../../../../mol-data/db';

import Schema = Column.Schema;

const float = Schema.float;
const str = Schema.str;
const int = Schema.int;
const lstr = Schema.lstr;
const Aliased = Schema.Aliased;

export const SF_Schema = {
    /**
     * Data items in the CELL category record details about the
     * crystallographic cell parameters.
     */
    cell: {
        /**
         * Unit-cell angle alpha of the reported structure in degrees.
         */
        angle_alpha: float,
        /**
         * Unit-cell angle beta of the reported structure in degrees.
         */
        angle_beta: float,
        /**
         * Unit-cell angle gamma of the reported structure in degrees.
         */
        angle_gamma: float,
        /**
         * This data item is a pointer to _entry.id in the ENTRY category.
         */
        entry_id: str,
        /**
         * Unit-cell length a corresponding to the structure reported in
         * angstroms.
         */
        length_a: float,
        /**
         * Unit-cell length b corresponding to the structure reported in
         * angstroms.
         */
        length_b: float,
        /**
         * Unit-cell length c corresponding to the structure reported in
         * angstroms.
         */
        length_c: float,
        /**
         * The number of the polymeric chains in a unit cell. In the case
         * of heteropolymers, Z is the number of occurrences of the most
         * populous chain.
         *
         * This data item is provided for compatibility with the original
         * Protein Data Bank format, and only for that purpose.
         */
        Z_PDB: int,
        /**
         * To further identify unique axis if necessary.  E.g., P 21 with
         * an unique C axis will have 'C' in this field.
         */
        pdbx_unique_axis: str,
    },
    /**
     * Data items in the DIFFRN category record details about the
     * diffraction data and their measurement.
     */
    diffrn: {
        /**
         * The mean temperature in kelvins at which the intensities were
         * measured.
         */
        ambient_temp: float,
        /**
         * This data item is a pointer to _exptl_crystal.id in the
         * EXPTL_CRYSTAL category.
         */
        crystal_id: str,
        /**
         * Remarks about how the crystal was treated prior to intensity
         * measurement. Particularly relevant when intensities were
         * measured at low temperature.
         */
        crystal_treatment: str,
        /**
         * Special details of the diffraction measurement process. Should
         * include information about source instability, crystal motion,
         * degradation and so on.
         */
        details: str,
        /**
         * This data item uniquely identifies a set of diffraction
         * data.
         */
        id: str,
    },
    /**
     * Data items in the DIFFRN_RADIATION_WAVELENGTH category
     * describe the wavelength of the radiation used to measure the
     * diffraction intensities. Items may be looped to identify
     * and assign weights to distinct components of a
     * polychromatic beam.
     */
    diffrn_radiation_wavelength: {
        /**
         * The code identifying each value of
         * _diffrn_radiation_wavelength.wavelength.
         * Items in the DIFFRN_RADIATION_WAVELENGTH category are looped
         * when multiple wavelengths are used.
         *
         * This code is used to link with the DIFFRN_REFLN category.
         * The _diffrn_refln.wavelength_id codes must match one of
         * the codes defined in this category.
         */
        id: str,
        /**
         * The radiation wavelength in angstroms.
         */
        wavelength: float,
    },
    /**
     * Data items in the DIFFRN_REFLNS category record details about
     * the set of intensities measured in the diffraction experiment.
     *
     * The DIFFRN_REFLN data items refer to individual intensity
     * measurements and must be included in looped lists.
     *
     * The DIFFRN_REFLNS data items specify the parameters that apply
     * to all intensity measurements in a diffraction data set.
     */
    diffrn_reflns: {
        /**
         * This data item is a pointer to _diffrn.id in the DIFFRN
         * category.
         */
        diffrn_id: str,
        /**
         * The maximum value of the Miller index h for the
         * reflection data specified by _diffrn_refln.index_h.
         */
        limit_h_max: int,
        /**
         * The minimum value of the Miller index h for the
         * reflection data specified by _diffrn_refln.index_h.
         */
        limit_h_min: int,
        /**
         * The maximum value of the Miller index k for the
         * reflection data specified by _diffrn_refln.index_k.
         */
        limit_k_max: int,
        /**
         * The minimum value of the Miller index k for the
         * reflection data specified by _diffrn_refln.index_k.
         */
        limit_k_min: int,
        /**
         * The maximum value of the Miller index l for the
         * reflection data specified by _diffrn_refln.index_l.
         */
        limit_l_max: int,
        /**
         * The minimum value of the Miller index l for the
         * reflection data specified by _diffrn_refln.index_l.
         */
        limit_l_min: int,
        /**
         * The total number of measured intensities, excluding reflections
         * that are classified as systematically absent.
         */
        number: int,
        /**
         * The lowest resolution for the interplanar spacings in the
         * reflection data set. This is the largest d value.
         */
        pdbx_d_res_low: float,
        /**
         * The highest resolution for the interplanar spacings in the
         * reflection data set. This is the smallest d value.
         */
        pdbx_d_res_high: float,
        /**
         * The  number of reflections satisfying the observation criterion
         * as in _diffrn_reflns.pdbx_observed_criterion
         */
        pdbx_number_obs: int,
    },
    /**
     * There is only one item in the ENTRY category, _entry.id. This
     * data item gives a name to this entry and is indirectly a key to
     * the categories (such as CELL, GEOM, EXPTL) that describe
     * information pertinent to the entire data block.
     */
    entry: {
        /**
         * The value of _entry.id identifies the data block.
         *
         * Note that this item need not be a number; it can be any unique
         * identifier.
         */
        id: str,
    },
    /**
     * Data items in the EXPTL_CRYSTAL category record the results of
     * experimental measurements on the crystal or crystals used,
     * such as shape, size or density.
     */
    exptl_crystal: {
        /**
         * The value of _exptl_crystal.id must uniquely identify a record in
         * the EXPTL_CRYSTAL list.
         *
         * Note that this item need not be a number; it can be any unique
         * identifier.
         */
        id: str,
    },
    /**
     * Data items in the REFLN category record details about the
     * reflection data used to determine the ATOM_SITE data items.
     *
     * The REFLN data items refer to individual reflections and must
     * be included in looped lists.
     *
     * The REFLNS data items specify the parameters that apply to all
     * reflections. The REFLNS data items are not looped.
     */
    refln: {
        /**
         * This data item is a pointer to _exptl_crystal.id in the
         * EXPTL_CRYSTAL category.
         */
        crystal_id: str,
        /**
         * The measured value of the structure factor in arbitrary units.
         */
        F_meas_au: float,
        /**
         * The standard uncertainty (estimated standard deviation) of
         * _refln.F_meas_au in arbitrary units.
         */
        F_meas_sigma_au: float,
        /**
         * The figure of merit m for this reflection.
         *
         * int P~alpha~ exp(i*alpha) dalpha
         * m = --------------------------------
         * int P~alpha~ dalpha
         *
         * P~a~ = the probability that the phase angle a is correct
         *
         * int is taken over the range alpha = 0 to 2 pi.
         */
        fom: float,
        /**
         * Miller index h of the reflection. The values of the Miller
         * indices in the REFLN category must correspond to the cell
         * defined by cell lengths and cell angles in the CELL category.
         */
        index_h: int,
        /**
         * Miller index k of the reflection. The values of the Miller
         * indices in the REFLN category must correspond to the cell
         * defined by cell lengths and cell angles in the CELL category.
         */
        index_k: int,
        /**
         * Miller index l of the reflection. The values of the Miller
         * indices in the REFLN category must correspond to the cell
         * defined by cell lengths and cell angles in the CELL category.
         */
        index_l: int,
        /**
         * The measured value of the intensity.
         */
        intensity_meas: float,
        /**
         * The standard uncertainty (derived from measurement) of the
         * intensity in the same units as _refln.intensity_meas.
         */
        intensity_sigma: float,
        /**
         * Classification of a reflection so as to indicate its status with
         * respect to inclusion in the refinement and the calculation of
         * R factors.
         */
        status: Aliased<'o' | '<' | '-' | 'x' | 'h' | 'l' | 'f'>(lstr),
        /**
         * This data item is a pointer to _reflns_scale.group_code in the
         * REFLNS_SCALE category.
         */
        scale_group_code: str,
        /**
         * This data item is a pointer to _diffrn_radiation.wavelength_id in
         * the DIFFRN_RADIATION category.
         */
        wavelength_id: str,
        /**
         * The weighted structure factor amplitude for the 2mFo-DFc map.
         */
        pdbx_FWT: float,
        /**
         * The weighted phase for the 2mFo-DFc map.
         */
        pdbx_PHWT: float,
        /**
         * The weighted structure factor amplitude for the mFo-DFc map.
         */
        pdbx_DELFWT: float,
        /**
         * The weighted phase for the mFo-DFc map.
         */
        pdbx_DELPHWT: float,
        /**
         * The R-free flag originally assigned to the reflection.  The convention used for
         * labeling the work and test sets differs depending on choice of data processing
         * software and refinement program.
         */
        pdbx_r_free_flag: int,
    },
    /**
     * Data items in the REFLNS_SCALE category record details about
     * the structure-factor scales. They are referenced from within
     * the REFLN list through _refln.scale_group_code.
     */
    reflns_scale: {
        /**
         * The code identifying a scale _reflns_scale.meas_F,
         * _reflns_scale.meas_F_squared or _reflns_scale.meas_intensity.
         * These are linked to the REFLN list by the
         * _refln.scale_group_code. These codes
         * need not correspond to those in the DIFFRN_SCALE list.
         */
        group_code: str,
    },
    /**
     * Data items in the SYMMETRY category record details about the
     * space-group symmetry.
     */
    symmetry: {
        /**
         * This data item is a pointer to _entry.id in the ENTRY category.
         */
        entry_id: str,
        /**
         * The cell settings for this space-group symmetry.
         */
        cell_setting: Aliased<'triclinic' | 'monoclinic' | 'orthorhombic' | 'tetragonal' | 'rhombohedral' | 'trigonal' | 'hexagonal' | 'cubic'>(lstr),
        /**
         * Space-group number from International Tables for Crystallography
         * Vol. A (2002).
         */
        Int_Tables_number: int,
        /**
         * Space-group symbol as described by Hall (1981). This symbol
         * gives the space-group setting explicitly. Leave spaces between
         * the separate components of the symbol.
         *
         * Ref: Hall, S. R. (1981). Acta Cryst. A37, 517-525; erratum
         * (1981) A37, 921.
         */
        space_group_name_Hall: str,
        /**
         * Hermann-Mauguin space-group symbol. Note that the
         * Hermann-Mauguin symbol does not necessarily contain complete
         * information about the symmetry and the space-group origin. If
         * used, always supply the FULL symbol from International Tables
         * for Crystallography Vol. A (2002) and indicate the origin and
         * the setting if it is not implicit. If there is any doubt that
         * the equivalent positions can be uniquely deduced from this
         * symbol, specify the  _symmetry_equiv.pos_as_xyz or
         * _symmetry.space_group_name_Hall  data items as well. Leave
         * spaces between symbols referring to
         * different axes.
         */
        'space_group_name_H-M': str,
    },
};

export type SF_Schema = typeof SF_Schema;
export interface SF_Database extends Database<SF_Schema> {};