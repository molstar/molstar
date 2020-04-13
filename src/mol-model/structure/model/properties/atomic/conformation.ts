/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../../../mol-data/db';
import UUID from '../../../../../mol-util/uuid';

export interface AtomicConformation {
    id: UUID,

    // ID is part of conformation because mmCIF is a leaky abstraction
    // that assigns different atom ids to corresponding atoms in different models
    // ... go figure.

    /**
     * Uniquely identifies an atom, i.e. even accross models in a mmCIF file.
     */
    atomId: Column<number>,

    /**
     * The fraction of the atom type present at this site. The sum of the occupancies of all the atom types
     * at this site may not significantly exceed 1.0 unless it is a dummy site.
     */
    occupancy: Column<number>,
    /**
     * Isotropic atomic displacement parameter, or equivalent isotropic atomic displacement parameter,
     * B~eq~, calculated from the anisotropic displacement parameters.
     */
    B_iso_or_equiv: Column<number>

    // Coordinates. Generally, not to be accessed directly because the coordinate might be
    // transformed by an operator. Use Unit.getPosition instead.

    /**
     * Are xyz coordinates defined?
     */
    xyzDefined: boolean,
    /**
     * The x coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    x: ArrayLike<number>,
    /**
     * The y coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    y: ArrayLike<number>,
    /**
     * The z coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    z: ArrayLike<number>
}