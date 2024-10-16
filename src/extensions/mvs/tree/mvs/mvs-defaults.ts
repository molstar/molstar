/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultsForTree } from '../generic/tree-schema';
import { MVSTreeSchema } from './mvs-tree';


/** Default values for params in `MVSTree` */
export const MVSDefaults = {
    root: {},
    download: {
    },
    parse: {
    },
    structure: {
        block_header: null,
        block_index: 0,
        model_index: 0,
        assembly_id: null,
        radius: 5,
        ijk_min: [-1, -1, -1],
        ijk_max: [1, 1, 1],
    },
    component: {
        selector: 'all' as const,
    },
    component_from_uri: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'component',
        field_values: null,
    },
    component_from_source: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'component',
        field_values: null,
    },
    representation: {
    },
    color: {
        selector: 'all' as const,
    },
    color_from_uri: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'color',
    },
    color_from_source: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'color',
    },
    transparency: {
    },
    label: {
    },
    label_from_uri: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'label',
    },
    label_from_source: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'label',
    },
    tooltip: {
    },
    tooltip_from_uri: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'tooltip',
    },
    tooltip_from_source: {
        block_header: null,
        block_index: 0,
        category_name: null,
        field_name: 'tooltip',
    },
    focus: {
        direction: [0, 0, -1],
        up: [0, 1, 0],
    },
    transform: {
        rotation: [1, 0, 0, 0, 1, 0, 0, 0, 1], // 3x3 identitity matrix
        translation: [0, 0, 0],
    },
    canvas: {
    },
    camera: {
        up: [0, 1, 0],
    },
    primitives: {
        color: null,
        label_color: null,
        tooltip: null,
        transparency: null,
        label_transparency: null,
        instances: null,
    },
    primitives_from_uri: {
        references: null,
    },
    primitive: { },
} satisfies DefaultsForTree<typeof MVSTreeSchema>;

/** Color to be used e.g. for representations without 'color' node */
export const DefaultColor = 'white';
