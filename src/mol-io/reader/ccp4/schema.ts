/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface Ccp4Header {
    /** number of columns (fastest changing) */
    NC: number
    /** number of rows */
    NR: number
    /** number of sections (slowest changing) */
    NS: number
    /**
     * 0 image : signed 8-bit bytes range -128 to 127
     * 1 image : 16-bit halfwords
     * 2 image : 32-bit reals
     * 3 transform : complex 16-bit integers
     * 4 transform : complex 32-bit reals
     * 6 image : unsigned 16-bit range 0 to 65535
     * 16 image: unsigned char * 3 (for rgb data, non-standard)
     *
     * Note: Mode 2 is the normal mode used in the CCP4 programs.
     * This parser only supports modes 0 and 2
     */
    MODE: number
    /** first column */
    NCSTART: number
    /** first row */
    NRSTART: number
    /** first section */
    NSSTART: number
    /** intervals along x */
    NX: number
    /** intervals along y */
    NY: number
    /** intervals along z */
    NZ: number
    /** x axis cell length (Angstroms in CCP4) */
    xLength: number
    /** y axis cell length (Angstroms in CCP4) */
    yLength: number
    /** z axis cell length (Angstroms in CCP4) */
    zLength: number
    /** alpha cell angle (Degrees) */
    alpha: number
    /** beta cell angle (Degrees) */
    beta: number
    /** gamma cell angle (Degrees) */
    gamma: number
    /** axis corresponds to columns (1,2,3 for X,Y,Z) */
    MAPC: number
    /** axis corresponds to rows (1,2,3 for X,Y,Z) */
    MAPR: number
    /** axis corresponds to sections (1,2,3 for X,Y,Z) */
    MAPS: number
    /** minimum density value */
    AMIN: number
    /** maximum density value */
    AMAX: number
    /** mean/average density value */
    AMEAN: number
    /** space group number */
    ISPG: number
    /** number of bytes used for storing symmetry operators */
    NSYMBT: number
    /** flag for skew transformation, =0 none, =1 if foll */
    LSKFLG: number
    /**
     * Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0
     * May be used in CCP4 but not in MRC
     */
    SKWMAT: number[]
    /**
     * Skew translation t if LSKFLG != 0
     * Skew transformation is from standard orthogonal
     * coordinate frame (as used for atoms) to orthogonal
     * map frame, as Xo(map) = S * (Xo(atoms) - t)
     *
     * May be used in CCP4 but not in MRC
     */
    SKWTRN: number[]
    /** see https://github.com/uglymol/uglymol/blob/master/tools/mapmode2to0#L69 */
    userFlag1: number,
    userFlag2: number,
    /** x axis origin transformation (not used in CCP4) */
    originX: number
    /** y axis origin transformation (not used in CCP4) */
    originY: number
    /** z axis origin transformation (not used in CCP4) */
    originZ: number
    /** Character string 'MAP ' to identify file type */
    MAP: string
    /**
     * Machine stamp indicating machine type which wrote file,
     * 17 and 17 for big-endian or 68 and 65 for little-endian.
     */
    MACHST: number[]
    /** rms deviation of map from mean density */
    ARMS: number
}

/**
 * MRC
 * http://ami.scripps.edu/software/mrctools/mrc_specification.php
 * http://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/image-processing-software/#image
 * http://bio3d.colorado.edu/imod/doc/mrc_format.txt
 *
 * CCP4 (MAP)
 * http://www.ccp4.ac.uk/html/maplib.html
 *
 * MRC format does not use the skew transformation header records (words 25-37)
 * CCP4 format does not use the ORIGIN header records (words 50-52)
 */
export interface Ccp4File {
    name: string
    header: Ccp4Header
    values: Float32Array | Int16Array | Int8Array | Uint16Array
}