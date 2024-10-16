/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { deepClone, pickObjectKeys } from '../../../../mol-util/object';
import { MVSData } from '../../mvs-data';
import { CustomProps } from '../generic/tree-schema';
import { MVSDefaults } from './mvs-defaults';
import { MVSKind, MVSNode, MVSSubtree, MVSNodeParams, MVSTreeSchema } from './mvs-tree';


/** Create a new MolViewSpec builder containing only a root node. Example of MVS builder usage:
 *
 * ```
 * const builder = createMVSBuilder();
 * builder.canvas({ background_color: 'white' });
 * const struct = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og2_updated.cif' }).parse({ format: 'mmcif' }).modelStructure();
 * struct.component().representation().color({ color: '#3050F8' });
 * console.log(JSON.stringify(builder.getState()));
 * ```
 */
export function createMVSBuilder(params: CustomAndRef = {}) {
    return new Root(params);
}


/** Base class for MVS builder pointing to anything */
class _Base<TKind extends MVSKind> {
    protected constructor(
        protected readonly _root: Root,
        protected readonly _node: MVSSubtree<TKind>,
    ) { }
    /** Create a new node, append as child to current _node, and return the new node */
    protected addChild<TChildKind extends MVSKind>(kind: TChildKind, params_: MVSNodeParams<TChildKind> & CustomAndRef) {
        const { params, custom, ref } = splitParams<MVSNodeParams<TChildKind>>(params_);
        const allowedParamNames = Object.keys(MVSTreeSchema.nodes[kind].params) as (keyof MVSNodeParams<TChildKind>)[];
        const node = {
            kind,
            params: pickObjectKeys(params, allowedParamNames) as unknown,
            custom,
            ref,
        } as MVSSubtree<TChildKind>;
        this._node.children ??= [];
        this._node.children.push(node);
        return node;
    }
}


/** MVS builder pointing to the 'root' node */
export class Root extends _Base<'root'> {
    constructor(params_: CustomAndRef) {
        const { custom, ref } = params_;
        const node: MVSNode<'root'> = { kind: 'root', custom, ref };
        super(undefined as any, node);
        (this._root as Root) = this;
    }
    /** Return the current state of the builder as object in MVS format. */
    getState(metadata?: Partial<Pick<MVSData['metadata'], 'title' | 'description' | 'description_format'>>): MVSData {
        return {
            root: deepClone(this._node),
            metadata: {
                ...metadata,
                version: `${MVSData.SupportedVersion}`,
                timestamp: utcNowISO(),
            },
        };
    }
    // omitting `saveState`, filesystem operations are responsibility of the caller code (platform-dependent)

    /** Add a 'camera' node and return builder pointing to the root. 'camera' node instructs to set the camera position and orientation. */
    camera(params: MVSNodeParams<'camera'> & CustomAndRef): Root {
        this.addChild('camera', params);
        return this;
    }
    /** Add a 'canvas' node and return builder pointing to the root. 'canvas' node sets canvas properties. */
    canvas(params: MVSNodeParams<'canvas'> & CustomAndRef): Root {
        this.addChild('canvas', params);
        return this;
    }
    /** Add a 'download' node and return builder pointing to it. 'download' node instructs to retrieve a data resource. */
    download(params: MVSNodeParams<'download'> & CustomAndRef): Download {
        return new Download(this._root, this.addChild('download', params));
    }
}


/** MVS builder pointing to a 'download' node */
export class Download extends _Base<'download'> {
    /** Add a 'parse' node and return builder pointing to it. 'parse' node instructs to parse a data resource. */
    parse(params: MVSNodeParams<'parse'> & CustomAndRef) {
        return new Parse(this._root, this.addChild('parse', params));
    }
}


/** Subsets of 'structure' node params which will be passed to individual builder functions. */
const StructureParamsSubsets = {
    model: ['block_header', 'block_index', 'model_index'],
    assembly: ['block_header', 'block_index', 'model_index', 'assembly_id'],
    symmetry: ['block_header', 'block_index', 'model_index', 'ijk_min', 'ijk_max'],
    symmetry_mates: ['block_header', 'block_index', 'model_index', 'radius'],
} satisfies { [kind in MVSNodeParams<'structure'>['type']]: (keyof MVSNodeParams<'structure'>)[] };


/** MVS builder pointing to a 'parse' node */
export class Parse extends _Base<'parse'> {
    /** Add a 'structure' node representing a "model structure", i.e. includes all coordinates from the original model without applying any transformations.
     * Return builder pointing to the new node. */
    modelStructure(params: Pick<MVSNodeParams<'structure'>, typeof StructureParamsSubsets['model'][number]> & CustomAndRef = {}): Structure {
        return new Structure(this._root, this.addChild('structure', {
            type: 'model',
            ...pickObjectKeys(params, [...StructureParamsSubsets.model]),
            custom: params.custom,
            ref: params.ref,
        }));
    }
    /** Add a 'structure' node representing an "assembly structure", i.e. may apply filters and symmetry operators to the original model coordinates.
     * Return builder pointing to the new node. */
    assemblyStructure(params: Pick<MVSNodeParams<'structure'>, typeof StructureParamsSubsets['assembly'][number]> & CustomAndRef = {}): Structure {
        return new Structure(this._root, this.addChild('structure', {
            type: 'assembly',
            ...pickObjectKeys(params, StructureParamsSubsets.assembly),
            custom: params.custom,
            ref: params.ref,
        }));
    }
    /** Add a 'structure' node representing a "symmetry structure", i.e. applies symmetry operators to build crystal unit cells within given Miller indices.
     * Return builder pointing to the new node. */
    symmetryStructure(params: Pick<MVSNodeParams<'structure'>, typeof StructureParamsSubsets['symmetry'][number]> & CustomAndRef = {}): Structure {
        return new Structure(this._root, this.addChild('structure', {
            type: 'symmetry',
            ...pickObjectKeys(params, StructureParamsSubsets.symmetry),
            custom: params.custom,
            ref: params.ref,
        }));
    }
    /** Add a 'structure' node representing a "symmetry mates structure", i.e. applies symmetry operators to build asymmetric units within a radius from the original model.
     * Return builder pointing to the new node. */
    symmetryMatesStructure(params: Pick<MVSNodeParams<'structure'>, typeof StructureParamsSubsets['symmetry_mates'][number]> & CustomAndRef = {}): Structure {
        return new Structure(this._root, this.addChild('structure', {
            type: 'symmetry_mates',
            ...pickObjectKeys(params, StructureParamsSubsets.symmetry_mates),
            custom: params.custom,
            ref: params.ref,
        }));
    }
}


/** MVS builder pointing to a 'structure' node */
export class Structure extends _Base<'structure'> {
    /** Add a 'component' node and return builder pointing to it. 'component' node instructs to create a component (i.e. a subset of the parent structure). */
    component(params: Partial<MVSNodeParams<'component'>> & CustomAndRef = {}): Component {
        const fullParams = { ...params, selector: params.selector ?? MVSDefaults.component.selector };
        return new Component(this._root, this.addChild('component', fullParams));
    }
    /** Add a 'component_from_uri' node and return builder pointing to it. 'component_from_uri' node instructs to create a component defined by an external annotation resource. */
    componentFromUri(params: MVSNodeParams<'component_from_uri'> & CustomAndRef): Component {
        return new Component(this._root, this.addChild('component_from_uri', params));
    }
    /** Add a 'component_from_source' node and return builder pointing to it. 'component_from_source' node instructs to create a component defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
    componentFromSource(params: MVSNodeParams<'component_from_source'> & CustomAndRef): Component {
        return new Component(this._root, this.addChild('component_from_source', params));
    }
    /** Add a 'label_from_uri' node and return builder pointing back to the structure node. 'label_from_uri' node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an external annotation resource. */
    labelFromUri(params: MVSNodeParams<'label_from_uri'> & CustomAndRef): Structure {
        this.addChild('label_from_uri', params);
        return this;
    }
    /** Add a 'label_from_source' node and return builder pointing back to the structure node. 'label_from_source' node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
    labelFromSource(params: MVSNodeParams<'label_from_source'> & CustomAndRef): Structure {
        this.addChild('label_from_source', params);
        return this;
    }
    /** Add a 'tooltip_from_uri' node and return builder pointing back to the structure node. 'tooltip_from_uri' node instructs to add tooltips to parts of a structure. The tooltips are defined by an external annotation resource. */
    tooltipFromUri(params: MVSNodeParams<'tooltip_from_uri'> & CustomAndRef): Structure {
        this.addChild('tooltip_from_uri', params);
        return this;
    }
    /** Add a 'tooltip_from_source' node and return builder pointing back to the structure node. 'tooltip_from_source' node instructs to add tooltips to parts of a structure. The tooltips are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
    tooltipFromSource(params: MVSNodeParams<'tooltip_from_source' & CustomAndRef>): Structure {
        this.addChild('tooltip_from_source', params);
        return this;
    }
    /** Add a 'transform' node and return builder pointing back to the structure node. 'transform' node instructs to rotate and/or translate structure coordinates. */
    transform(params: MVSNodeParams<'transform'> & CustomAndRef = {}): Structure {
        if (params.rotation && params.rotation.length !== 9) {
            throw new Error('ValueError: `rotation` parameter must be an array of 9 numbers');
        }
        this.addChild('transform', params);
        return this;
    }
}


/** MVS builder pointing to a 'component' or 'component_from_uri' or 'component_from_source' node */
export class Component extends _Base<'component' | 'component_from_uri' | 'component_from_source'> {
    /** Add a 'representation' node and return builder pointing to it. 'representation' node instructs to create a visual representation of a component. */
    representation(params: Partial<MVSNodeParams<'representation'>> & CustomAndRef = {}): Representation {
        const fullParams: MVSNodeParams<'representation'> = { ...params, type: params.type ?? 'cartoon' };
        return new Representation(this._root, this.addChild('representation', fullParams));
    }
    /** Add a 'label' node and return builder pointing back to the component node. 'label' node instructs to add a label (textual visual representation) to a component. */
    label(params: MVSNodeParams<'label'> & CustomAndRef): Component {
        this.addChild('label', params);
        return this;
    }
    /** Add a 'tooltip' node and return builder pointing back to the component node. 'tooltip' node instructs to add a text which is not a part of the visualization but should be presented to the users when they interact with the component (typically, the tooltip will be shown somewhere on the screen when the user hovers over a visual representation of the component). */
    tooltip(params: MVSNodeParams<'tooltip'> & CustomAndRef): Component {
        this.addChild('tooltip', params);
        return this;
    }
    /** Add a 'focus' node and return builder pointing back to the component node. 'focus' node instructs to set the camera focus to a component (zoom in). */
    focus(params: MVSNodeParams<'focus'> & CustomAndRef = {}): Component {
        this.addChild('focus', params);
        return this;
    }
}


/** MVS builder pointing to a 'representation' node */
export class Representation extends _Base<'representation'> {
    /** Add a 'color' node and return builder pointing back to the representation node. 'color' node instructs to apply color to a visual representation. */
    color(params: MVSNodeParams<'color'> & CustomAndRef): Representation {
        this.addChild('color', params);
        return this;
    }
    /** Add a 'color_from_uri' node and return builder pointing back to the representation node. 'color_from_uri' node instructs to apply colors to a visual representation. The colors are defined by an external annotation resource. */
    colorFromUri(params: MVSNodeParams<'color_from_uri'> & CustomAndRef): Representation {
        this.addChild('color_from_uri', params);
        return this;
    }
    /** Add a 'color_from_source' node and return builder pointing back to the representation node. 'color_from_source' node instructs to apply colors to a visual representation. The colors are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
    colorFromSource(params: MVSNodeParams<'color_from_source'> & CustomAndRef): Representation {
        this.addChild('color_from_source', params);
        return this;
    }
}


/** Demonstration of usage of MVS builder */
export function builderDemo() {
    const builder = createMVSBuilder();
    builder.canvas({ background_color: 'white' });
    const struct = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og2_updated.cif' }).parse({ format: 'mmcif' }).modelStructure();
    struct.component().representation().color({ color: 'white' });
    struct.component({ selector: 'ligand' }).representation({ type: 'ball_and_stick', custom: { repr_quality: 'high' }, ref: 'Ligand' })
        .color({ color: '#555555' })
        .color({ selector: { type_symbol: 'N' }, color: '#3050F8' })
        .color({ selector: { type_symbol: 'O' }, color: '#FF0D0D' })
        .color({ selector: { type_symbol: 'S' }, color: '#FFFF30' })
        .color({ selector: { type_symbol: 'FE' }, color: '#E06633' });
    builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og5_updated.cif' }).parse({ format: 'mmcif' }).assemblyStructure({ assembly_id: '1' }).component().representation().color({ color: 'cyan' });
    builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og5_updated.cif' }).parse({ format: 'mmcif' }).assemblyStructure({ assembly_id: '2' }).component().representation().color({ color: 'blue' });
    const cif = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1wrf_updated.cif' }).parse({ format: 'mmcif' });

    cif.modelStructure({ model_index: 0 }).component().representation().color({ color: '#CC0000' });
    cif.modelStructure({ model_index: 1 }).component().representation().color({ color: '#EE7700' });
    cif.modelStructure({ model_index: 2 }).component().representation().color({ color: '#FFFF00' });

    cif.modelStructure({ model_index: 0 }).transform({ translation: [30, 0, 0] }).component().representation().color({ color: '#ff88bb' });
    cif.modelStructure({ model_index: 0 as any }).transform({ translation: [60, 0, 0], rotation: [0, 1, 0, -1, 0, 0, 0, 0, 1] }).component().representation().color({ color: '#aa0077' });

    return builder.getState();
}

/** Return the current universal time, in ISO format, e.g. '2023-11-24T10:45:49.873Z' */
function utcNowISO(): string {
    return new Date().toISOString();
}

export interface CustomAndRef {
    custom?: CustomProps,
    ref?: string,
};

function splitParams<TParams extends {}>(params_custom_ref: TParams & CustomAndRef): { params: TParams, custom?: CustomProps, ref?: string } {
    const { custom, ref, ...params } = params_custom_ref;
    return { params: params as TParams, custom, ref };
}
