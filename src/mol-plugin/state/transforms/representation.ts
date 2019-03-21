/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
import { PluginContext } from 'mol-plugin/context';
import { RepresentationProvider } from 'mol-repr/representation';
import { BuiltInStructureRepresentationsName } from 'mol-repr/structure/registry';
import { StructureParams } from 'mol-repr/structure/representation';
import { BuiltInVolumeRepresentationsName } from 'mol-repr/volume/registry';
import { VolumeParams } from 'mol-repr/volume/representation';
import { StateTransformer } from 'mol-state';
import { Task } from 'mol-task';
import { BuiltInColorThemeName, ColorTheme } from 'mol-theme/color';
import { BuiltInSizeThemeName, SizeTheme } from 'mol-theme/size';
import { createTheme, ThemeRegistryContext } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { Text } from 'mol-geo/geometry/text/text';
import { ColorNames } from 'mol-util/color/tables';
import { getLabelRepresentation } from 'mol-plugin/util/structure-labels';
import { ShapeRepresentation } from 'mol-repr/shape/representation';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';
import { unwindStructureAssembly, explodeStructure, getStructureOverpaint } from '../animation/helpers';
import { Color } from 'mol-util/color';
import { Overpaint } from 'mol-theme/overpaint';
import { BaseGeometry } from 'mol-geo/geometry/base';

export { StructureRepresentation3D }
export { StructureRepresentation3DHelpers }
export { StructureLabels3D}
export { ExplodeStructureRepresentation3D }
export { UnwindStructureAssemblyRepresentation3D }
export { OverpaintStructureRepresentation3D as ColorStructureRepresentation3D }
export { VolumeRepresentation3D }

namespace StructureRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structure: Structure, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);

        const themeDataCtx = { structure };
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const structureDefaultParams = PD.getDefaultValues(type.getParams(ctx.structureRepresentation.themeCtx, structure))
        return ({
            type: { name, params: structureParams ? { ...structureDefaultParams, ...structureParams } : structureDefaultParams },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function createParams<R extends RepresentationProvider<Structure, any, any>, C extends ColorTheme.Provider<any>, S extends SizeTheme.Provider<any>>(
            ctx: PluginContext, structure: Structure, params: {
            repr?: R | [R, (r: R, ctx: ThemeRegistryContext, s: Structure) => RepresentationProvider.ParamValues<R>],
            color?: C | [C, (c: C, ctx: ThemeRegistryContext) => ColorTheme.ParamValues<C>],
            size?: S | [S, (c: S, ctx: ThemeRegistryContext) => SizeTheme.ParamValues<S>]
        }): StateTransformer.Params<StructureRepresentation3D> {

        const themeCtx = ctx.structureRepresentation.themeCtx

        const repr = params.repr
            ? params.repr instanceof Array ? params.repr[0] : params.repr
            : ctx.structureRepresentation.registry.default.provider;
        const reprParams = params.repr instanceof Array
            ? params.repr[1](repr as R, themeCtx, structure)
            : PD.getDefaultValues(repr.getParams(themeCtx, structure));

        const color = params.color
            ? params.color instanceof Array ? params.color[0] : params.color
            : themeCtx.colorThemeRegistry.get(repr.defaultColorTheme);
        const colorParams = params.color instanceof Array
            ? params.color[1](color as C, themeCtx)
            : PD.getDefaultValues(color.getParams(themeCtx));

        const size = params.size
            ? params.size instanceof Array ? params.size[0] : params.size
            : themeCtx.sizeThemeRegistry.get(repr.defaultSizeTheme);
        const sizeParams = params.size instanceof Array
            ? params.size[1](size as S, themeCtx)
            : PD.getDefaultValues(size.getParams(themeCtx));

        return ({
            type: { name: ctx.structureRepresentation.registry.getName(repr), params: reprParams },
            colorTheme: { name: themeCtx.colorThemeRegistry.getName(color), params: colorParams },
            sizeTheme: { name: themeCtx.sizeThemeRegistry.getName(size), params: sizeParams }
        })
    }

    export function getDefaultParamsWithTheme(ctx: PluginContext, reprName: BuiltInStructureRepresentationsName, colorName: BuiltInColorThemeName | undefined, structure: Structure, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(reprName);

        const themeDataCtx = { structure };
        const color = colorName || type.defaultColorTheme;
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(color).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const structureDefaultParams = PD.getDefaultValues(type.getParams(ctx.structureRepresentation.themeCtx, structure))
        return ({
            type: { name: reprName, params: structureParams ? { ...structureDefaultParams, ...structureParams } : structureDefaultParams },
            colorTheme: { name: color, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).defaultValues;
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).defaultValues
        return ({
            type: { name, params: structureParams ? { ...type.defaultValues, ...structureParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: colorParams },
            sizeTheme: { name: type.defaultSizeTheme, params: sizeParams }
        })
    }
}

type StructureRepresentation3D = typeof StructureRepresentation3D
const StructureRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'structure-representation-3d',
    display: '3D Representation',
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure.Representation3D,
    params: (a, ctx: PluginContext) => {
        const { registry, themeCtx } = ctx.structureRepresentation
        const type = registry.get(registry.default.name);

        if (!a) {
            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, Structure.Empty))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ structure: Structure.Empty }))
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ structure: Structure.Empty }))
                )
            }
        }

        const dataCtx = { structure: a.data }
        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.types,
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme,
                themeCtx.sizeThemeRegistry.types,
                name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx))
            )
        })
    }
})({
    canAutoUpdate({ a, oldParams, newParams }) {
        // TODO: other criteria as well?
        return a.data.elementCount < 10000 || oldParams.type.name === newParams.type.name;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const provider = plugin.structureRepresentation.registry.get(params.type.name)
            const props = params.type.params || {}
            const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.structureRepresentation.themeCtx }, provider.getParams)
            repr.setTheme(createTheme(plugin.structureRepresentation.themeCtx, { structure: a.data }, params))
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Molecule.Structure.Representation3D({ repr, source: a } , { label: provider.label });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;
            const props = { ...b.data.repr.props, ...newParams.type.params }
            b.data.repr.setTheme(createTheme(plugin.structureRepresentation.themeCtx, { structure: a.data }, newParams));
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});


type StructureLabels3D = typeof StructureLabels3D
const StructureLabels3D = PluginStateTransform.BuiltIn({
    name: 'structure-labels-3d',
    display: '3D Labels',
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure.Representation3D,
    params: {
        // TODO: other targets
        target: PD.MappedStatic('residues', {
            'elements': PD.Group({ }),
            'residues': PD.Group({ }),
            'static-text': PD.Group({
                value: PD.Text(''),
                size: PD.Optional(PD.Numeric(1, { min: 1, max: 1000, step: 0.1 })),
                // TODO: this changes the position while rotated etc... fix
                position: PD.Optional(Text.Params.attachment)
            }, { isFlat: true })
        }),
        options: PD.Group({
            ...Text.Params,

            background: PD.Boolean(true),
            backgroundMargin: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
            backgroundColor: PD.Color(ColorNames.snow),
            backgroundOpacity: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }),
        })
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return (oldParams.target.name === 'static-text' && newParams.target.name === 'static-text' && oldParams.target.params.value === newParams.target.params.value)
            || newParams.target.name === oldParams.target.name;
    },
    apply({ a, params }) {
        return Task.create('Structure Labels', async ctx => {
            const repr = await getLabelRepresentation(ctx, a.data, params);
            return new SO.Molecule.Structure.Representation3D({ repr, source: a }, { label: `Labels`, description: params.target.name });
        });
    },
    update({ a, b, newParams }) {
        return Task.create('Structure Labels', async ctx => {
            await getLabelRepresentation(ctx, a.data, newParams, b.data.repr as ShapeRepresentation<any, any, any>);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});

type UnwindStructureAssemblyRepresentation3D = typeof UnwindStructureAssemblyRepresentation3D
const UnwindStructureAssemblyRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'unwind-structure-assembly-representation-3d',
    display: 'Unwind Assembly 3D Representation',
    from: SO.Molecule.Structure.Representation3D,
    to: SO.Molecule.Structure.Representation3DState,
    params: { t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }) }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }) {
        const structure = a.data.source.data;
        const unitTransforms = new StructureUnitTransforms(structure);
        unwindStructureAssembly(structure, unitTransforms, params.t);
        return new SO.Molecule.Structure.Representation3DState({
            state: { unitTransforms },
            initialState: { unitTransforms: new StructureUnitTransforms(structure) },
            info: structure,
            source: a
        }, { label: `Unwind T = ${params.t.toFixed(2)}` });
    },
    update({ a, b, newParams, oldParams }) {
        const structure = b.data.info as Structure;
        if (a.data.source.data !== structure) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.t === newParams.t) return StateTransformer.UpdateResult.Unchanged;
        const unitTransforms = b.data.state.unitTransforms!;
        unwindStructureAssembly(structure, unitTransforms, newParams.t);
        b.label = `Unwind T = ${newParams.t.toFixed(2)}`;
        b.data.source = a;
        return StateTransformer.UpdateResult.Updated;
    }
});


type ExplodeStructureRepresentation3D = typeof ExplodeStructureRepresentation3D
const ExplodeStructureRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'explode-structure-representation-3d',
    display: 'Explode 3D Representation',
    from: SO.Molecule.Structure.Representation3D,
    to: SO.Molecule.Structure.Representation3DState,
    params: { t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }) }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params, spine }) {
        const rootStructure = spine.getRootOfType(SO.Molecule.Structure)!.data;
        const structure = a.data.source.data;
        const unitTransforms = new StructureUnitTransforms(rootStructure);
        explodeStructure(structure, unitTransforms, params.t);
        return new SO.Molecule.Structure.Representation3DState({
            state: { unitTransforms },
            initialState: { unitTransforms: new StructureUnitTransforms(rootStructure) },
            info: rootStructure,
            source: a
        }, { label: `Explode T = ${params.t.toFixed(2)}` });
    },
    update({ a, b, newParams, oldParams, spine }) {
        const rootStructure = spine.getRootOfType(SO.Molecule.Structure)!.data;
        if (b.data.info !== rootStructure) return StateTransformer.UpdateResult.Recreate;
        if (oldParams.t === newParams.t) return StateTransformer.UpdateResult.Unchanged;
        const unitTransforms = b.data.state.unitTransforms!;
        explodeStructure(rootStructure, unitTransforms, newParams.t);
        b.label = `Explode T = ${newParams.t.toFixed(2)}`;
        b.data.source = a;
        return StateTransformer.UpdateResult.Updated;
    }
});

type OverpaintStructureRepresentation3D = typeof OverpaintStructureRepresentation3D
const OverpaintStructureRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'overpaint-structure-representation-3d',
    display: 'Overpaint 3D Representation',
    from: SO.Molecule.Structure.Representation3D,
    to: SO.Molecule.Structure.Representation3DState,
    params: {
        layers: PD.ObjectList({
            script: PD.ScriptExpression({ language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname LYS))' }),
            color: PD.Color(ColorNames.blueviolet)
        }, e => `${Color.toRgbString(e.color)}`, {
            defaultValue: [
                {
                    script: {
                        language: 'mol-script',
                        expression: '(sel.atom.atom-groups :residue-test (= atom.resname LYS))'
                    },
                    color: ColorNames.blueviolet
                },
                {
                    script: {
                        language: 'mol-script',
                        expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))'
                    },
                    color: ColorNames.chartreuse
                }
            ]
        }),
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity' }),
    }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }) {
        const structure = a.data.source.data
        const overpaint = getStructureOverpaint(structure, params.layers, params.alpha)

        return new SO.Molecule.Structure.Representation3DState({
            state: { overpaint },
            initialState: { overpaint: Overpaint.Empty },
            info: structure,
            source: a
        }, { label: `Overpaint (${overpaint.layers.length} Layers)` })
    },
    update({ a, b, newParams, oldParams }) {
        const structure = b.data.info as Structure
        if (a.data.source.data !== structure) return StateTransformer.UpdateResult.Recreate
        const oldOverpaint = b.data.state.overpaint!
        const newOverpaint = getStructureOverpaint(structure, newParams.layers, newParams.alpha)
        if (oldParams.alpha === newParams.alpha && Overpaint.areEqual(oldOverpaint, newOverpaint)) return StateTransformer.UpdateResult.Unchanged

        b.data.state.overpaint = newOverpaint
        b.data.source = a
        b.label = `Overpaint (${newOverpaint.layers.length} Layers)`
        return StateTransformer.UpdateResult.Updated
    }
});

//

export namespace VolumeRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInVolumeRepresentationsName, volume: VolumeData, volumeParams?: Partial<PD.Values<VolumeParams>>): StateTransformer.Params<VolumeRepresentation3D> {
        const type = ctx.volumeRepresentation.registry.get(name);

        const themeDataCtx = { volume };
        const colorParams = ctx.volumeRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.volumeRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const volumeDefaultParams = PD.getDefaultValues(type.getParams(ctx.volumeRepresentation.themeCtx, volume))
        return ({
            type: { name, params: volumeParams ? { ...volumeDefaultParams, ...volumeParams } : volumeDefaultParams },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInVolumeRepresentationsName, volumeParams?: Partial<PD.Values<PD.Params>>, colorName?: BuiltInColorThemeName, colorParams?: Partial<ColorTheme.Props>, sizeName?: BuiltInSizeThemeName, sizeParams?: Partial<SizeTheme.Props>): StateTransformer.Params<VolumeRepresentation3D> {
        const type = ctx.volumeRepresentation.registry.get(name);
        const colorType = ctx.volumeRepresentation.themeCtx.colorThemeRegistry.get(colorName || type.defaultColorTheme);
        const sizeType = ctx.volumeRepresentation.themeCtx.sizeThemeRegistry.get(sizeName || type.defaultSizeTheme);
        return ({
            type: { name, params: volumeParams ? { ...type.defaultValues, ...volumeParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: colorParams ? { ...colorType.defaultValues, ...colorParams } : colorType.defaultValues },
            sizeTheme: { name: type.defaultSizeTheme, params: sizeParams ? { ...sizeType.defaultValues, ...sizeParams } : sizeType.defaultValues }
        })
    }

    export function getDescription(props: any) {
        return props.isoValue && VolumeIsoValue.toString(props.isoValue)
    }
}
type VolumeRepresentation3D = typeof VolumeRepresentation3D
const VolumeRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'volume-representation-3d',
    display: '3D Representation',
    from: SO.Volume.Data,
    to: SO.Volume.Representation3D,
    params: (a, ctx: PluginContext) => {
        const { registry, themeCtx } = ctx.volumeRepresentation
        const type = registry.get(registry.default.name);

        if (!a) {
            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, VolumeData.One ))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ volume: VolumeData.One }))
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ volume: VolumeData.One }))
                )
            }
        }

        const dataCtx = { volume: a.data }
        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.types,
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme,
                themeCtx.sizeThemeRegistry.types,
                name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx))
            )
        })
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        // TODO: allow for small molecules
        return oldParams.type.name === newParams.type.name;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Volume Representation', async ctx => {
            const provider = plugin.volumeRepresentation.registry.get(params.type.name)
            const props = params.type.params || {}
            const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.volumeRepresentation.themeCtx }, provider.getParams)
            repr.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: a.data }, params))
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Volume.Representation3D({ repr, source: a }, { label: provider.label, description: VolumeRepresentation3DHelpers.getDescription(props) });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Volume Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;
            const props = { ...b.data.repr.props, ...newParams.type.params }
            b.data.repr.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: a.data }, newParams))
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            b.description = VolumeRepresentation3DHelpers.getDescription(props)
            return StateTransformer.UpdateResult.Updated;
        });
    }
});

//

export { ShapeRepresentation3D }
type ShapeRepresentation3D = typeof ShapeRepresentation3D
const ShapeRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'shape-representation-3d',
    display: '3D Representation',
    from: SO.Shape.Provider,
    to: SO.Shape.Representation3D,
    params: (a, ctx: PluginContext) => {
        return a ? a.data.params : BaseGeometry.Params
    }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Shape Representation', async ctx => {
            const props = { ...PD.getDefaultValues(a.data.params), params }
            const repr = ShapeRepresentation(a.data.getShape, a.data.geometryUtils)
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data.data).runInContext(ctx);
            return new SO.Shape.Representation3D({ repr, source: a }, { label: a.data.label });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Shape Representation', async ctx => {
            const props = { ...b.data.repr.props, ...newParams }
            await b.data.repr.createOrUpdate(props, a.data.data).runInContext(ctx);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});