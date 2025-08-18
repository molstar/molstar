/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, dict, float, int, list, literal, nullable, OptionalField, RequiredField, str, tuple, union } from '../generic/field-schema';
import { SimpleParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema, TreeSchemaWithAllRequired } from '../generic/tree-schema';
import { MVSClipParams, MVSRepresentationParams, MVSVolumeRepresentationParams } from './mvs-tree-representations';
import { MVSPrimitiveParams } from './mvs-tree-primitives';
import { ColorT, ComponentExpressionT, ComponentSelectorT, Matrix, Palette, ParseFormatT, SchemaFormatT, SchemaT, StrList, StructureTypeT, Vector3 } from './param-types';


const _DataFromUriParams = {
    /** URL of the annotation resource. */
    uri: RequiredField(str, 'URL of the annotation resource.'),
    /** Format of the annotation resource. */
    format: RequiredField(SchemaFormatT, 'Format of the annotation resource.'),
    /** Annotation schema defines what fields in the annotation will be taken into account. */
    schema: RequiredField(SchemaT, 'Annotation schema defines what fields in the annotation will be taken into account.'),
    /** Header of the CIF block to read annotation from (only applies when `format` is `"cif"` or `"bcif"`). If `null`, block is selected based on `block_index`. */
    block_header: OptionalField(nullable(str), null, 'Header of the CIF block to read annotation from (only applies when `format` is `"cif"` or `"bcif"`). If `null`, block is selected based on `block_index`.'),
    /** 0-based index of the CIF block to read annotation from (only applies when `format` is `"cif"` or `"bcif"` and `block_header` is `null`). */
    block_index: OptionalField(int, 0, '0-based index of the CIF block to read annotation from (only applies when `format` is `"cif"` or `"bcif"` and `block_header` is `null`).'),
    /** Name of the CIF category to read annotation from (only applies when `format` is `"cif"` or `"bcif"`). If `null`, the first category in the block is used. */
    category_name: OptionalField(nullable(str), null, 'Name of the CIF category to read annotation from (only applies when `format` is `"cif"` or `"bcif"`). If `null`, the first category in the block is used.'),
    /** Name of the column in CIF or field name (key) in JSON that contains the dependent variable (color/label/tooltip/component_id...). The default value is 'color'/'label'/'tooltip'/'component' depending on the node type */
    field_name: RequiredField(str, 'Name of the column in CIF or field name (key) in JSON that contains the dependent variable (color/label/tooltip/component_id...).'),
    /** Optional remapping of annotation field names `{ standardName1: actualName1, ... }`. Use `{ "label_asym_id": "X" }` to load actual field "X" as "label_asym_id". Use `{ "label_asym_id": null }` to ignore actual field "label_asym_id". Fields not mentioned here are mapped implicitely (i.e. actual name = standard name). */
    field_remapping: OptionalField(dict(str, nullable(str)), {}, 'Optional remapping of annotation field names `{ standardName1: actualName1, ... }`. Use `{ "label_asym_id": "X" }` to load actual field "X" as "label_asym_id". Use `{ "label_asym_id": null }` to ignore actual field "label_asym_id". Fields not mentioned here are mapped implicitely (i.e. actual name = standard name).'),
};

const _DataFromSourceParams = {
    /** Annotation schema defines what fields in the annotation will be taken into account. */
    schema: RequiredField(SchemaT, 'Annotation schema defines what fields in the annotation will be taken into account.'),
    /** Header of the CIF block to read annotation from. If `null`, block is selected based on `block_index`. */
    block_header: OptionalField(nullable(str), null, 'Header of the CIF block to read annotation from. If `null`, block is selected based on `block_index`.'),
    /** 0-based index of the CIF block to read annotation from (only applies when `block_header` is `null`). */
    block_index: OptionalField(int, 0, '0-based index of the CIF block to read annotation from (only applies when `block_header` is `null`).'),
    /** Name of the CIF category to read annotation from. If `null`, the first category in the block is used. */
    category_name: OptionalField(nullable(str), null, 'Name of the CIF category to read annotation from. If `null`, the first category in the block is used.'),
    /** Name of the column in CIF or field name (key) in JSON that contains the dependent variable (color/label/tooltip/component_id...). The default value is 'color'/'label'/'tooltip'/'component' depending on the node type */
    field_name: RequiredField(str, 'Name of the column in CIF or field name (key) in JSON that contains the dependent variable (color/label/tooltip/component_id...).'),
    /** Optional remapping of annotation field names `{ standardName1: actualName1, ... }`. Use `{ "label_asym_id": "X" }` to load actual field "X" as "label_asym_id". Use `{ "label_asym_id": null }` to ignore actual field "label_asym_id". Fields not mentioned here are mapped implicitely (i.e. actual name = standard name). */
    field_remapping: OptionalField(dict(str, nullable(str)), {}, 'Optional remapping of annotation field names `{ standardName1: actualName1, ... }`. Use `{ "label_asym_id": "X" }` to load actual field "X" as "label_asym_id". Use `{ "label_asym_id": null }` to ignore actual field "label_asym_id". Fields not mentioned here are mapped implicitely (i.e. actual name = standard name).'),
};

/** Color to be used e.g. for representations without 'color' node */
export const DefaultColor = 'white';

const LabelAttachments = literal('bottom-left', 'bottom-center', 'bottom-right', 'middle-left', 'middle-center', 'middle-right', 'top-left', 'top-center', 'top-right');

const TransformParams = SimpleParamsSchema({
    /** Rotation matrix (3x3 matrix flattened in column major format (j*3+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. The default value is the identity matrix (corresponds to no rotation). */
    rotation: OptionalField(Matrix, [1, 0, 0, 0, 1, 0, 0, 0, 1], 'Rotation matrix (3x3 matrix flattened in column major format (j*3+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. The default value is the identity matrix (corresponds to no rotation).'),
    /** Translation vector, applied to the structure coordinates after rotation. The default value is the zero vector (corresponds to no translation). */
    translation: OptionalField(Vector3, [0, 0, 0], 'Translation vector, applied to the structure coordinates after rotation. The default value is the zero vector (corresponds to no translation).'),
    /** Point to rotate the object around. Can be either a 3D vector or dynamically computed object centroid. */
    rotation_center: OptionalField(nullable(union(Vector3, literal('centroid'))), null, 'Point to rotate the object around. Can be either a 3D vector or dynamically computed object centroid.'),
    /** Transform matrix (4x4 matrix flattened in column major format (j*4+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. Takes precedence over `rotation` and `translation`. */
    matrix: OptionalField(nullable(Matrix), null, 'Transform matrix (4x4 matrix flattened in column major format (j*4+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. Takes precedence over `rotation` and `translation`.'),
});

/** Schema for `MVSTree` (MolViewSpec tree) */
export const MVSTreeSchema = TreeSchema({
    rootKind: 'root',
    nodes: {
        /** Auxiliary node kind that only appears as the tree root. */
        root: {
            description: 'Auxiliary node kind that only appears as the tree root.',
            parent: [],
            params: SimpleParamsSchema({
            }),
        },
        /** This node instructs to retrieve a data resource. */
        download: {
            description: 'This node instructs to retrieve a data resource.',
            parent: ['root'],
            params: SimpleParamsSchema({
                /** URL of the data resource. */
                url: RequiredField(str, 'URL of the data resource.'),
            }),
        },
        /** This node instructs to parse a data resource. */
        parse: {
            description: 'This node instructs to parse a data resource.',
            parent: ['download'],
            params: SimpleParamsSchema({
                /** Format of the input data resource. */
                format: RequiredField(ParseFormatT, 'Format of the input data resource.'),
            }),
        },
        /** This node instructs to retrieve molecular coordinates from a parsed data resource. */
        coordinates: {
            description: 'This node instructs to retrieve molecular coordinates from a parsed data resource.',
            parent: ['parse'],
            params: SimpleParamsSchema({}),
        },
        /** This node instructs to create a structure from a parsed data resource. "Structure" refers to an internal representation of molecular coordinates without any visual representation. */
        structure: {
            description: 'This node instructs to create a structure from a parsed data resource. "Structure" refers to an internal representation of molecular coordinates without any visual representation.',
            parent: ['parse'],
            params: SimpleParamsSchema({
                /** Type of structure to be created (`"model"` for original model coordinates, `"assembly"` for assembly structure, `"symmetry"` for a set of crystal unit cells based on Miller indices, `"symmetry_mates"` for a set of asymmetric units within a radius from the original model). */
                type: RequiredField(StructureTypeT, 'Type of structure to be created (`"model"` for original model coordinates, `"assembly"` for assembly structure, `"symmetry"` for a set of crystal unit cells based on Miller indices, `"symmetry_mates"` for a set of asymmetric units within a radius from the original model).'),
                /** Header of the CIF block to read coordinates from (only applies when the input data are from CIF or BinaryCIF). If `null`, block is selected based on `block_index`. */
                block_header: OptionalField(nullable(str), null, 'Header of the CIF block to read coordinates from (only applies when the input data are from CIF or BinaryCIF). If `null`, block is selected based on `block_index`.'),
                /** 0-based index of the CIF block to read coordinates from (only applies when the input data are from CIF or BinaryCIF and `block_header` is `null`). */
                block_index: OptionalField(int, 0, '0-based index of the CIF block to read coordinates from (only applies when the input data are from CIF or BinaryCIF and `block_header` is `null`).'),
                /** 0-based index of model in case the input data contain multiple models. */
                model_index: OptionalField(int, 0, '0-based index of model in case the input data contain multiple models.'),
                /** Assembly identifier (only applies when `kind` is `"assembly"`). If `null`, the first assembly is selected. */
                assembly_id: OptionalField(nullable(str), null, 'Assembly identifier (only applies when `kind` is `"assembly"`). If `null`, the first assembly is selected.'),
                /** Distance (in Angstroms) from the original model in which asymmetric units should be included (only applies when `kind` is `"symmetry_mates"`). */
                radius: OptionalField(float, 5, 'Distance (in Angstroms) from the original model in which asymmetric units should be included (only applies when `kind` is `"symmetry_mates"`).'),
                /** Miller indices of the bottom-left unit cell to be included (only applies when `kind` is `"symmetry"`). */
                ijk_min: OptionalField(tuple([int, int, int]), [-1, -1, -1], 'Miller indices of the bottom-left unit cell to be included (only applies when `kind` is `"symmetry"`).'),
                /** Miller indices of the top-right unit cell to be included (only applies when `kind` is `"symmetry"`). */
                ijk_max: OptionalField(tuple([int, int, int]), [1, 1, 1], 'Miller indices of the top-right unit cell to be included (only applies when `kind` is `"symmetry"`).'),
                /** Reference to a specific set of coordinates. */
                coordinates_ref: OptionalField(nullable(str), null, 'Reference to a specific set of coordinates.')
            }),
        },
        /** This node instructs to rotate and/or translate structure coordinates. */
        transform: {
            description: 'This node instructs to rotate and/or translate coordinates OR provide a transformation matrix.',
            parent: ['structure', 'component', 'volume'],
            params: TransformParams,
        },
        /** This node allows instantiation using the provided transformation parameters. */
        instance: {
            description: 'This node allows instantiation using the provided transformation parameters.',
            parent: ['structure', 'component', 'volume'],
            params: TransformParams,
        },
        /** This node instructs to create a component (i.e. a subset of the parent structure). */
        component: {
            description: 'This node instructs to create a component (i.e. a subset of the parent structure).',
            parent: ['structure'],
            params: SimpleParamsSchema({
                /** Defines what part of the parent structure should be included in this component. */
                selector: RequiredField(union(ComponentSelectorT, ComponentExpressionT, list(ComponentExpressionT)), 'Defines what part of the parent structure should be included in this component.'),
            }),
        },
        /** This node instructs to create a component defined by an external annotation resource. */
        component_from_uri: {
            description: 'This node instructs to create a component defined by an external annotation resource.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromUriParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the component identifier. */
                field_name: OptionalField(str, 'component', 'Name of the column in CIF or field name (key) in JSON that contains the component identifier.'),
                /** List of component identifiers (i.e. values in the field given by `field_name`) which should be included in this component. If `null`, component identifiers are ignored (all annotation rows are included), and `field_name` field can be dropped from the annotation. */
                field_values: OptionalField(nullable(list(str)), null, 'List of component identifiers (i.e. values in the field given by `field_name`) which should be included in this component. If `null`, component identifiers are ignored (all annotation rows are included), and `field_name` field can be dropped from the annotation.'),
            }),
        },
        /** This node instructs to create a component defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
        component_from_source: {
            description: 'This node instructs to create a component defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromSourceParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the component identifier. */
                field_name: OptionalField(str, 'component', 'Name of the column in CIF or field name (key) in JSON that contains the component identifier.'),
                /** List of component identifiers (i.e. values in the field given by `field_name`) which should be included in this component. If `null`, component identifiers are ignored (all annotation rows are included), and `field_name` field can be dropped from the annotation. */
                field_values: OptionalField(nullable(list(str)), null, 'List of component identifiers (i.e. values in the field given by `field_name`) which should be included in this component. If `null`, component identifiers are ignored (all annotation rows are included), and `field_name` field can be dropped from the annotation.'),
            }),
        },
        /** This node instructs to create a visual representation of a component. */
        representation: {
            description: 'This node instructs to create a visual representation of a component.',
            parent: ['component', 'component_from_uri', 'component_from_source'],
            params: MVSRepresentationParams,
        },
        /** This node instructs to create a volume from a parsed data resource. "Volume" refers to an internal representation of volumetric data without any visual representation. */
        volume: {
            description: 'This node instructs to create a volume from a parsed data resource. "Volume" refers to an internal representation of volumetric data without any visual representation.',
            parent: ['parse'],
            params: SimpleParamsSchema({
                channel_id: OptionalField(nullable(str), null, 'Channel identifier (only applies when the input data contain multiple channels).'),
            }),
        },
        /** This node instructs to create a visual representation of a volume. */
        volume_representation: {
            description: 'This node instructs to create a visual representation of a volume.',
            parent: ['volume'],
            params: MVSVolumeRepresentationParams,
        },
        /** This node instructs to apply color to a visual representation. */
        color: {
            description: 'This node instructs to apply color to a visual representation.',
            parent: ['representation', 'volume_representation'],
            params: SimpleParamsSchema({
                /** Color to apply to the representation. Can be either an X11 color name (e.g. `"red"`) or a hexadecimal code (e.g. `"#FF0011"`). */
                color: OptionalField(ColorT, DefaultColor, 'Color to apply to the representation. Can be either an X11 color name (e.g. `"red"`) or a hexadecimal code (e.g. `"#FF0011"`).'),
                /** Defines to what part of the representation this color should be applied. */
                selector: OptionalField(union(ComponentSelectorT, ComponentExpressionT, list(ComponentExpressionT)), 'all', 'Defines to what part of the representation this color should be applied.'),
            }),
        },
        /** This node instructs to apply colors to a visual representation. The colors are defined by an external annotation resource. */
        color_from_uri: {
            description: 'This node instructs to apply colors to a visual representation. The colors are defined by an external annotation resource.',
            parent: ['representation'],
            params: SimpleParamsSchema({
                ..._DataFromUriParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the color. */
                field_name: OptionalField(str, 'color', 'Name of the column in CIF or field name (key) in JSON that contains the color.'),
                /** Customize mapping of annotation values to colors. */
                palette: OptionalField(nullable(Palette), null, 'Customize mapping of annotation values to colors.'),
            }),
        },
        /** This node instructs to apply colors to a visual representation. The colors are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
        color_from_source: {
            description: 'This node instructs to apply colors to a visual representation. The colors are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file.',
            parent: ['representation'],
            params: SimpleParamsSchema({
                ..._DataFromSourceParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the color. */
                field_name: OptionalField(str, 'color', 'Name of the column in CIF or field name (key) in JSON that contains the color.'),
                /** Customize mapping of annotation values to colors. */
                palette: OptionalField(nullable(Palette), null, 'Customize mapping of annotation values to colors.'),
            }),
        },
        /** This node instructs to apply clipping to a visual representation. */
        clip: {
            description: 'This node instructs to apply clipping to a visual representation.',
            parent: ['representation', 'volume_representation'],
            params: MVSClipParams,
        },
        /** This node instructs to apply opacity/transparency to a visual representation. */
        opacity: {
            description: 'This node instructs to apply opacity/transparency to a visual representation.',
            parent: ['representation', 'volume_representation'],
            params: SimpleParamsSchema({
                /** Opacity of a representation. 0.0: fully transparent, 1.0: fully opaque. */
                opacity: RequiredField(float, 'Opacity of a representation. 0.0: fully transparent, 1.0: fully opaque.'),
            }),
        },
        /** This node instructs to add a label (textual visual representation) to a component. */
        label: {
            description: 'This node instructs to add a label (textual visual representation) to a component.',
            parent: ['component', 'component_from_uri', 'component_from_source'],
            params: SimpleParamsSchema({
                /** Content of the shown label. */
                text: RequiredField(str, 'Content of the shown label.'),
            }),
        },
        /** This node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an external annotation resource. */
        label_from_uri: {
            description: 'This node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an external annotation resource.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromUriParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the label text. */
                field_name: OptionalField(str, 'label', 'Name of the column in CIF or field name (key) in JSON that contains the label text.'),
            }),
        },
        /** This node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
        label_from_source: {
            description: 'This node instructs to add labels (textual visual representations) to parts of a structure. The labels are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromSourceParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the label text. */
                field_name: OptionalField(str, 'label', 'Name of the column in CIF or field name (key) in JSON that contains the label text.'),
            }),
        },
        /** This node instructs to add a tooltip to a component. "Tooltip" is a text which is not a part of the visualization but should be presented to the users when they interact with the component (typically, the tooltip will be shown somewhere on the screen when the user hovers over a visual representation of the component). */
        tooltip: {
            description: 'This node instructs to add a tooltip to a component. "Tooltip" is a text which is not a part of the visualization but should be presented to the users when they interact with the component (typically, the tooltip will be shown somewhere on the screen when the user hovers over a visual representation of the component).',
            parent: ['component', 'component_from_uri', 'component_from_source'],
            params: SimpleParamsSchema({
                /** Content of the shown tooltip. */
                text: RequiredField(str, 'Content of the shown tooltip.'),
            }),
        },
        /** This node instructs to add tooltips to parts of a structure. The tooltips are defined by an external annotation resource. */
        tooltip_from_uri: {
            description: 'This node instructs to add tooltips to parts of a structure. The tooltips are defined by an external annotation resource.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromUriParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the tooltip text. */
                field_name: OptionalField(str, 'tooltip', 'Name of the column in CIF or field name (key) in JSON that contains the tooltip text.'),
            }),
        },
        /** This node instructs to add tooltips to parts of a structure. The tooltips are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file. */
        tooltip_from_source: {
            description: 'This node instructs to add tooltips to parts of a structure. The tooltips are defined by an annotation resource included in the same file this structure was loaded from. Only applicable if the structure was loaded from an mmCIF or BinaryCIF file.',
            parent: ['structure'],
            params: SimpleParamsSchema({
                ..._DataFromSourceParams,
                /** Name of the column in CIF or field name (key) in JSON that contains the tooltip text. */
                field_name: OptionalField(str, 'tooltip', 'Name of the column in CIF or field name (key) in JSON that contains the tooltip text.'),
            }),
        },
        /** This node instructs to set the camera focus to a component (zoom in). */
        focus: {
            description: 'This node instructs to set the camera focus to a component (zoom in).',
            parent: ['root', 'component', 'component_from_uri', 'component_from_source', 'primitives', 'primitives_from_uri', 'volume', 'volume_representation'],
            params: SimpleParamsSchema({
                /** Vector describing the direction of the view (camera position -> focused target). */
                direction: OptionalField(Vector3, [0, 0, -1], 'Vector describing the direction of the view (camera position -> focused target).'),
                /** Vector which will be aligned with the screen Y axis. */
                up: OptionalField(Vector3, [0, 1, 0], 'Vector which will be aligned with the screen Y axis.'),
                /** Radius of the focused sphere (overrides `radius_factor` and `radius_extra`. */
                radius: OptionalField(nullable(float), null, 'Radius of the focused sphere (overrides `radius_factor` and `radius_extra`).'),
                /** Radius of the focused sphere relative to the radius of parent component (default: 1). Focused radius = component_radius * radius_factor + radius_extent. */
                radius_factor: OptionalField(float, 1, 'Radius of the focused sphere relative to the radius of parent component (default: 1). Focused radius = component_radius * radius_factor + radius_extent.'),
                /** Addition to the radius of the focused sphere, if computed from the radius of parent component (default: 0). Focused radius = component_radius * radius_factor + radius_extent. */
                radius_extent: OptionalField(float, 0, 'Addition to the radius of the focused sphere, if computed from the radius of parent component (default: 0). Focused radius = component_radius * radius_factor + radius_extent.'),
            }),
        },
        /** This node instructs to set the camera position and orientation. */
        camera: {
            description: 'This node instructs to set the camera position and orientation.',
            parent: ['root'],
            params: SimpleParamsSchema({
                /** Coordinates of the point in space at which the camera is pointing. */
                target: RequiredField(Vector3, 'Coordinates of the point in space at which the camera is pointing.'),
                /** Coordinates of the camera. */
                position: RequiredField(Vector3, 'Coordinates of the camera.'),
                /** Vector which will be aligned with the screen Y axis. */
                up: OptionalField(Vector3, [0, 1, 0], 'Vector which will be aligned with the screen Y axis.'),
            }),
        },
        /** This node sets canvas properties. */
        canvas: {
            description: 'This node sets canvas properties.',
            parent: ['root'],
            params: SimpleParamsSchema({
                /** Color of the canvas background. Can be either an X11 color name (e.g. `"red"`) or a hexadecimal code (e.g. `"#FF0011"`). */
                background_color: OptionalField(ColorT, 'white', 'Color of the canvas background. Can be either an X11 color name (e.g. `"red"`) or a hexadecimal code (e.g. `"#FF0011"`). Defaults to white.'),
            }),
        },
        primitives: {
            description: 'This node groups a list of geometrical primitives',
            parent: ['structure', 'root'],
            params: SimpleParamsSchema({
                /** Default color for primitives in this group. */
                color: OptionalField(ColorT, 'white', 'Default color for primitives in this group.'),
                /** Default label color for primitives in this group. */
                label_color: OptionalField(ColorT, 'white', 'Default label color for primitives in this group.'),
                /** Default tooltip for primitives in this group. */
                tooltip: OptionalField(nullable(str), null, 'Default tooltip for primitives in this group.'),
                /** Opacity of primitive geometry in this group. */
                opacity: OptionalField(float, 1, 'Opacity of primitive geometry in this group.'),
                /** Opacity of primitive labels in this group. */
                label_opacity: OptionalField(float, 1, 'Opacity of primitive labels in this group.'),
                /** Whether to show a tether line between the label and the target. Defaults to false. */
                label_show_tether: OptionalField(bool, false, 'Whether to show a tether line between the label and the target. Defaults to false.'),
                /** Length of the tether line between the label and the target. Defaults to 1 (Angstrom). */
                label_tether_length: OptionalField(float, 1, 'Length of the tether line between the label and the target. Defaults to 1 (Angstrom).'),
                /** How to attach the label to the target. Defaults to "middle-center". */
                label_attachment: OptionalField(LabelAttachments, 'middle-center', 'How to attach the label to the target. Defaults to "middle-center".'),
                /** Background color of the label. Defaults to none/transparent. */
                label_background_color: OptionalField(nullable(ColorT), null, 'Background color of the label. Defaults to none/transparent.'),
                /** Load snapshot with the provided key when interacting with this primitives group. */
                snapshot_key: OptionalField(nullable(str), null, 'Load snapshot with the provided key when interacting with this primitives group.'),
                /** Instances of this primitive group defined as 4x4 column major (j * 4 + i indexing) transformation matrices. */
                instances: OptionalField(nullable(list(Matrix)), null, 'Instances of this primitive group defined as 4x4 column major (j * 4 + i indexing) transformation matrices.'),
            }),
        },
        primitives_from_uri: {
            description: 'This node loads a list of primitives from URI',
            parent: ['structure', 'root'],
            params: SimpleParamsSchema({
                /** Location of the resource. */
                uri: RequiredField(str, 'Location of the resource.'),
                /** Format of the data. */
                format: RequiredField(literal('mvs-node-json'), 'Format of the data.'),
                /** List of nodes the data are referencing. */
                references: OptionalField(StrList, [], 'List of nodes the data are referencing.'),
            }),
        },
        primitive: {
            description: 'This node represents a geometrical primitive',
            parent: ['primitives'],
            params: MVSPrimitiveParams,
        },
    }
});

/** Node kind in a `MVSTree` */
export type MVSKind = keyof typeof MVSTreeSchema.nodes

/** Node in a `MVSTree` */
export type MVSNode<TKind extends MVSKind = MVSKind> = NodeFor<typeof MVSTreeSchema, TKind>

/** Params for a specific node kind in a `MVSTree` */
export type MVSNodeParams<TKind extends MVSKind> = ParamsOfKind<MVSTree, TKind>

/** MolViewSpec tree */
export type MVSTree = TreeFor<typeof MVSTreeSchema>

/** Any subtree in a `MVSTree` (e.g. its root doesn't need to be 'root') */
export type MVSSubtree<TKind extends MVSKind = MVSKind> = SubtreeOfKind<MVSTree, TKind>

/** Schema for `MVSTree` (MolViewSpec tree with all params provided) */
export const FullMVSTreeSchema = TreeSchemaWithAllRequired(MVSTreeSchema);

/** MolViewSpec tree with all params provided */
export type FullMVSTree = TreeFor<typeof FullMVSTreeSchema>