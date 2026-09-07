/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { isSffPrimaryDescriptor } from '../parser';

describe('SFF identification', () => {
    it('accepts the three descriptors defined by EMDB-SFF', () => {
        expect(isSffPrimaryDescriptor('mesh_list')).toBe(true);
        expect(isSffPrimaryDescriptor('three_d_volume')).toBe(true);
        expect(isSffPrimaryDescriptor('shape_primitive_list')).toBe(true);
    });

    it('rejects an HDF5 that carries no primary_descriptor', () => {
        // any other HDF5-based format - NeXus, EMAN2, Imaris, BigDataViewer
        expect(isSffPrimaryDescriptor(undefined)).toBe(false);
    });

    it('rejects an unknown primary_descriptor', () => {
        expect(isSffPrimaryDescriptor('NXentry')).toBe(false);
        expect(isSffPrimaryDescriptor('')).toBe(false);
    });
});
