/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject as PSO, PluginStateTransform } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { Clip } from '../../../mol-util/clip';
import { escapeRegExp, stringToWords } from '../../../mol-util/string';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { ParamMapping } from '../../../mol-util/param-mapping';
import { EntityNode } from '../ui/entities';
import { DistinctColorsProps, distinctColors } from '../../../mol-util/color/distinct';
import { Sphere3D } from '../../../mol-math/geometry';
import { Hcl } from '../../../mol-util/color/spaces/hcl';
import { StateObjectCell, StateObjectRef, StateSelection } from '../../../mol-state';
import { ShapeRepresentation3D, StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { SpacefillRepresentationProvider } from '../../../mol-repr/structure/representation/spacefill';
import { assertUnreachable } from '../../../mol-util/type-helpers';
import { MesoscaleExplorerState } from '../app';
import { saturate } from '../../../mol-math/interpolate';
import { Material } from '../../../mol-util/material';

function getHueRange(hue: number, variability: number) {
    let min = hue - variability;
    const minOverflow = (min < 0 ? -min : 0);
    let max = hue + variability;
    if (max > 360) min -= max - 360;
    max += minOverflow;
    return [Math.max(0, min), Math.min(360, max)] as [number, number];
}

function getGrayscaleColors(count: number, luminance: number, variability: number) {
    const out: Color[] = [];
    for (let i = 0; i < count; ++ i) {
        const l = saturate(luminance / 100);
        const v = saturate(variability / 180) * Math.random();
        const s = Math.random() > 0.5 ? 1 : -1;
        const d = Math.abs(l + s * v) % 1;
        out[i] = Color.fromNormalizedRgb(d, d, d);
    }
    return out;
}

export function getDistinctGroupColors(count: number, color: Color, variability: number, shift: number, props?: Partial<DistinctColorsProps>) {
    const hcl = Hcl.fromColor(Hcl(), color);
    if (isNaN(hcl[0])) {
        return getGrayscaleColors(count, hcl[2], variability);
    }

    if (count === 1) {
        hcl[1] = 65;
        hcl[2] = 55;
        return [Hcl.toColor(hcl)];
    }

    const colors = distinctColors(count, {
        hue: getHueRange(hcl[0], variability),
        chroma: [30, 100],
        luminance: [50, 100],
        clusteringStepCount: 0,
        minSampleCount: 1000,
        sampleCountFactor: 100,
        sort: 'none',
        ...props,
    });

    if (shift !== 0) {
        const offset = Math.floor(shift / 100 * count);
        return [...colors.slice(offset), ...colors.slice(0, offset)];
    } else {
        return colors;
    }
}

const Colors = [0x377eb8, 0xe41a1c, 0x4daf4a, 0x984ea3, 0xff7f00, 0xffff33, 0xa65628, 0xf781bf] as Color[];

export function getDistinctBaseColors(count: number, shift: number, props?: Partial<DistinctColorsProps>): Color[] {
    let colors: Color[];
    if (count <= Colors.length) {
        colors = Colors.slice(0, count).map(e => Array.isArray(e) ? e[0] : e);
    } else {
        colors = distinctColors(count, {
            hue: [1, 360],
            chroma: [25, 100],
            luminance: [30, 100],
            clusteringStepCount: 0,
            minSampleCount: 1000,
            sampleCountFactor: 100,
            sort: 'none',
            ...props,
        });
    }

    if (shift !== 0) {
        const offset = Math.floor(shift / 100 * count);
        return [...colors.slice(offset), ...colors.slice(0, offset)];
    } else {
        return colors;
    }
}

export const ColorParams = {
    type: PD.Select('generate', PD.arrayToOptions(['generate', 'uniform', 'custom'])),
    illustrative: PD.Boolean(false, { description: 'Illustrative style', hideIf: p => p.type === 'custom' }),
    value: PD.Color(Color(0xFFFFFF), { hideIf: p => p.type === 'custom' }),
    variability: PD.Numeric(20, { min: 1, max: 180, step: 1 }, { hideIf: p => p.type !== 'generate' }),
    shift: PD.Numeric(0, { min: 0, max: 100, step: 1 }, { hideIf: p => !p.type.includes('generate') }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }, { hideIf: p => p.type === 'custom' }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
    emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
};
export type ColorProps = PD.Values<typeof ColorParams>

export const ColorValueParam = PD.Color(Color(0xFFFFFF));

export const RootParams = {
    type: PD.Select('custom', PD.arrayToOptions(['group-generate', 'group-uniform', 'generate', 'uniform', 'custom'])),
    illustrative: PD.Boolean(false, { description: 'Illustrative style', hideIf: p => p.type === 'custom' }),
    value: PD.Color(Color(0xFFFFFF), { hideIf: p => p.type !== 'uniform' }),
    variability: PD.Numeric(20, { min: 1, max: 180, step: 1 }, { hideIf: p => p.type !== 'group-generate' }),
    shift: PD.Numeric(0, { min: 0, max: 100, step: 1 }, { hideIf: p => !p.type.includes('generate') }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }, { hideIf: p => p.type === 'custom' }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
    emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
};

export const LightnessParams = {
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
export const DimLightness = 6;

export const IllustrativeParams = {
    illustrative: PD.Boolean(false, { description: 'Illustrative style' }),
};

export const OpacityParams = {
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
};

export const EmissiveParams = {
    emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
};

export const celShaded = {
    celShaded: PD.Boolean(false, { description: 'Cel Shading light for stylized rendering of representations' })
};

export type celShadedProps = PD.Values<typeof celShaded>;


export const PatternParams = {
    frequency: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    amplitude: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
};

export const StyleParams = {
    ignoreLight: PD.Boolean(false, { description: 'Ignore light for stylized rendering of representations' }),
    materialStyle: Material.getParam(),
    celShaded: PD.Boolean(false, { description: 'Cel Shading light for stylized rendering of representations' }),
};

export const LodParams = {
    lodLevels: Spheres.Params.lodLevels,
    cellSize: Spheres.Params.cellSize,
    batchSize: Spheres.Params.batchSize,
    approximate: Spheres.Params.approximate,
};

export const SimpleClipParams = {
    type: PD.Select('none', PD.objectToOptions(Clip.Type, t => stringToWords(t))),
    invert: PD.Boolean(false),
    position: PD.Group({
        x: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(0, { min: -100, max: 100, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => g.type === 'none', isExpanded: true }),
    rotation: PD.Group({
        axis: PD.Vec3(Vec3.create(1, 0, 0)),
        angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => g.type === 'none', isExpanded: true }),
    scale: PD.Group({
        x: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { immediateUpdate: true }),
    }, { hideIf: g => ['none', 'plane'].includes(g.type), isExpanded: true }),
};
export type SimpleClipParams = typeof SimpleClipParams
export type SimpleClipProps = PD.Values<SimpleClipParams>

export function getClipObjects(values: SimpleClipProps, boundingSphere: Sphere3D): Clip.Props['objects'] {
    const { center, radius } = boundingSphere;

    const position = Vec3.clone(center);
    Vec3.add(position, position, Vec3.create(
        radius * values.position.x / 100,
        radius * values.position.y / 100,
        radius * values.position.z / 100
    ));

    const scale = Vec3.create(values.scale.x, values.scale.y, values.scale.z);
    Vec3.scale(scale, scale, 2 * radius / 100);

    return [{
        type: values.type,
        invert: values.invert,
        position,
        scale,
        rotation: values.rotation,
        transform: Mat4.identity(),
    }];
}

export function createClipMapping(node: EntityNode) {
    return ParamMapping({
        params: SimpleClipParams,
        target: (ctx: PluginContext) => {
            return node.clipValue;
        }
    })({
        values(props, ctx) {
            if (!props || props.objects.length === 0) {
                return {
                    type: 'none',
                    invert: false,
                    position: { x: 0, y: 0, z: 0 },
                    rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
                    scale: { x: 100, y: 100, z: 100 },
                };
            }

            const { center, radius } = node.plugin.canvas3d!.boundingSphere;
            const { invert, position, scale, rotation, type } = props.objects[0];

            const p = Vec3.clone(position);
            Vec3.sub(p, p, center);
            Vec3.scale(p, p, 100 / radius);
            Vec3.round(p, p);

            const s = Vec3.clone(scale);
            Vec3.scale(s, s, 100 / radius / 2);
            Vec3.round(s, s);

            return {
                type,
                invert,
                position: { x: p[0], y: p[1], z: p[2] },
                rotation,
                scale: { x: s[0], y: s[1], z: s[2] },
            };
        },
        update: (s, props) => {
            if (!props) return;

            const clipObjects = getClipObjects(s, node.plugin.canvas3d!.boundingSphere);
            props.objects = clipObjects;
        },
        apply: async (props, ctx) => {
            if (props) node.updateClip(props);
        }
    });
}

export const MesoscaleGroupParams = {
    root: PD.Value<boolean>(false, { isHidden: true }),
    index: PD.Value<number>(-1, { isHidden: true }),
    tag: PD.Value<string>('', { isHidden: true }),
    label: PD.Value<string>('', { isHidden: true }),
    description: PD.Value<string>('', { isHidden: true }),
    hidden: PD.Boolean(false),
    color: PD.Group(RootParams),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    lod: PD.Group(LodParams),
    clip: PD.Group(SimpleClipParams),
};
export type MesoscaleGroupProps = PD.Values<typeof MesoscaleGroupParams>;

export class MesoscaleGroupObject extends PSO.Create({ name: 'Mesoscale Group', typeClass: 'Object' }) { }

export const MesoscaleGroup = PluginStateTransform.BuiltIn({
    name: 'mesoscale-group',
    display: { name: 'Mesoscale Group' },
    from: [PSO.Root, MesoscaleGroupObject],
    to: MesoscaleGroupObject,
    params: MesoscaleGroupParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Apply Mesoscale Group', async () => {
            return new MesoscaleGroupObject({}, { label: params.label, description: params.description });
        });
    },
});

export function getMesoscaleGroupParams(graphicsMode: GraphicsMode): MesoscaleGroupProps {
    const groupParams = PD.getDefaultValues(MesoscaleGroupParams);
    if (graphicsMode === 'custom') return groupParams;

    return {
        ...groupParams,
        lod: {
            ...groupParams.lod,
            ...getGraphicsModeProps(graphicsMode),
        }
    };
}

//

export type LodLevels = typeof SpacefillRepresentationProvider.defaultValues['lodLevels']

export function getLodLevels(graphicsMode: Exclude<GraphicsMode, 'custom'>): LodLevels {
    switch (graphicsMode) {
        case 'performance':
            return [
                { minDistance: 1, maxDistance: 300, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 300, maxDistance: 2000, overlap: 0, stride: 40, scaleBias: 3 },
                { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 150, scaleBias: 3 },
                { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 300, scaleBias: 2.5 },
            ];
        case 'balanced':
            return [
                { minDistance: 1, maxDistance: 500, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 500, maxDistance: 2000, overlap: 0, stride: 15, scaleBias: 3 },
                { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 70, scaleBias: 2.7 },
                { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.5 },
            ];
        case 'quality':
            return [
                { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 1000, maxDistance: 4000, overlap: 0, stride: 10, scaleBias: 3 },
                { minDistance: 4000, maxDistance: 10000, overlap: 0, stride: 50, scaleBias: 2.7 },
                { minDistance: 10000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.3 },
            ];
        case 'ultra':
            return [
                { minDistance: 1, maxDistance: 5000, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 5000, maxDistance: 10000, overlap: 0, stride: 10, scaleBias: 3 },
                { minDistance: 10000, maxDistance: 30000, overlap: 0, stride: 50, scaleBias: 2.5 },
                { minDistance: 30000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2 },
            ];
        default:
            assertUnreachable(graphicsMode);
    }
}

export type GraphicsMode = 'ultra' | 'quality' | 'balanced' | 'performance' | 'custom';

export function getGraphicsModeProps(graphicsMode: Exclude<GraphicsMode, 'custom'>) {
    return {
        lodLevels: getLodLevels(graphicsMode),
        approximate: graphicsMode !== 'quality' && graphicsMode !== 'ultra',
        alphaThickness: graphicsMode === 'performance' ? 15 : 12,
    };
}

export function setGraphicsCanvas3DProps(ctx: PluginContext, graphics: GraphicsMode) {
    const pixelScale = graphics === 'balanced' ? 0.75
        : graphics === 'performance' ? 0.5 : 1;

    ctx.canvas3dContext?.setProps({ pixelScale });

    ctx.canvas3d?.setProps({
        postprocessing: {
            sharpening: pixelScale < 1 ? {
                name: 'on',
                params: { sharpness: 0.5, denoise: true }
            } : { name: 'off', params: {} }
        }
    });
}

//

export const MesoscaleStateParams = {
    filter: PD.Value<string>('', { isHidden: true }),
    graphics: PD.Select('quality', PD.arrayToOptions(['ultra', 'quality', 'balanced', 'performance', 'custom'] as GraphicsMode[])),
    description: PD.Value<string>('', { isHidden: true }),
    focusInfo: PD.Value<string>('', { isHidden: true }),
    link: PD.Value<string>('', { isHidden: true }),
    textSizeDescription: PD.Numeric(14, { min: 1, max: 100, step: 1 }, { isHidden: true }),
    index: PD.Value<number>(-1, { isHidden: true })
};

export class MesoscaleStateObject extends PSO.Create<MesoscaleState>({ name: 'Mesoscale State', typeClass: 'Object' }) { }

const MesoscaleStateTransform = PluginStateTransform.BuiltIn({
    name: 'mesoscale-state',
    display: { name: 'Mesoscale State' },
    from: PSO.Root,
    to: MesoscaleStateObject,
    params: MesoscaleStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Apply Mesoscale State', async () => {
            return new MesoscaleStateObject(params);
        });
    },
});

export { MesoscaleState };
type MesoscaleState = PD.Values<typeof MesoscaleStateParams>;
const MesoscaleState = {
    async init(ctx: PluginContext) {
        const cell = ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
        if (cell) throw new Error('MesoscaleState already initialized');

        const customState = ctx.customState as MesoscaleExplorerState;
        const state = await ctx.state.data.build().toRoot().apply(MesoscaleStateTransform, {
            filter: '',
            graphics: customState.graphicsMode,
        }).commit();
        customState.stateRef = state.ref;
    },
    get(ctx: PluginContext): MesoscaleState {
        const ref = this.ref(ctx);
        return ctx.state.data.tryGetCellData<MesoscaleStateObject>(ref);
    },
    async set(ctx: PluginContext, props: Partial<MesoscaleState>) {
        const ref = this.ref(ctx);
        await ctx.state.data.build().to(ref).update(MesoscaleStateTransform, old => Object.assign(old, props)).commit();
    },
    ref(ctx: PluginContext): string {
        const ref = (ctx.customState as MesoscaleExplorerState).stateRef;
        if (!ref) throw new Error('MesoscaleState not initialized');
        return ref;
    },
    has(ctx: PluginContext): boolean {
        const ref = (ctx.customState as MesoscaleExplorerState).stateRef || '';
        return ctx.state.data.cells.has(ref) ? true : false;
    },
};

//

export function getRoots(plugin: PluginContext): StateSelection.CellSeq<StateObjectCell<MesoscaleGroupObject>> {
    const s = plugin.customState as MesoscaleExplorerState;
    if (!s.stateCache.roots) {
        s.stateCache.roots = plugin.state.data.select(StateSelection.Generators.rootsOfType(MesoscaleGroupObject));
    }
    return s.stateCache.roots;
}

export function getGroups(plugin: PluginContext, tag?: string): StateSelection.CellSeq<StateObjectCell<MesoscaleGroupObject>> {
    const s = plugin.customState as MesoscaleExplorerState;
    const k = `groups-${tag || ''}`;
    if (!s.stateCache[k]) {
        const selector = tag !== undefined
            ? StateSelection.Generators.ofTransformer(MesoscaleGroup).withTag(tag)
            : StateSelection.Generators.ofTransformer(MesoscaleGroup);
        s.stateCache[k] = plugin.state.data.select(selector);
    }
    return s.stateCache[k];
}

function _getAllGroups(plugin: PluginContext, tag: string | undefined, list: StateObjectCell[]) {
    const groups = getGroups(plugin, tag);
    list.push(...groups);
    for (const g of groups) {
        _getAllGroups(plugin, g.params?.values.tag, list);
    }
    return list;
}

export function getAllGroups(plugin: PluginContext, tag?: string) {
    return _getAllGroups(plugin, tag, []);
}

export function getAllLeafGroups(plugin: PluginContext, tag: string) {
    const allGroups = getAllGroups(plugin, tag);
    allGroups.sort((a, b) => a.params?.values.index - b.params?.values.index);
    return allGroups.filter(g => {
        return getEntities(plugin, g.params?.values.tag).length > 0;
    });
}

type EntityCells = StateSelection.CellSeq<StateObjectCell<PSO.Molecule.Structure.Representation3D | PSO.Shape.Representation3D>>

export function getEntities(plugin: PluginContext, tag?: string): EntityCells {
    const s = plugin.customState as MesoscaleExplorerState;
    const k = `entities-${tag || ''}`;
    if (!s.stateCache[k]) {
        const structureSelector = tag !== undefined
            ? StateSelection.Generators.ofTransformer(StructureRepresentation3D).withTag(tag)
            : StateSelection.Generators.ofTransformer(StructureRepresentation3D);
        const shapeSelector = tag !== undefined
            ? StateSelection.Generators.ofTransformer(ShapeRepresentation3D).withTag(tag)
            : StateSelection.Generators.ofTransformer(ShapeRepresentation3D);
        s.stateCache[k] = [
            ...plugin.state.data.select(structureSelector).filter(c => c.obj!.data.sourceData.elementCount > 0),
            ...plugin.state.data.select(shapeSelector),
        ];
    }
    return s.stateCache[k];
}

function getFilterMatcher(filter: string) {
    return filter.startsWith('"') && filter.endsWith('"')
        ? new RegExp(`^${escapeRegExp(filter.substring(1, filter.length - 1))}$`, 'g')
        : new RegExp(escapeRegExp(filter), 'gi');
}

export function getFilteredEntities(plugin: PluginContext, tag: string, filter: string) {
    if (!filter) return getEntities(plugin, tag);
    const matcher = getFilterMatcher(filter);
    return getEntities(plugin, tag).filter(c => getEntityLabel(plugin, c).match(matcher) !== null);
}

function _getAllEntities(plugin: PluginContext, tag: string | undefined, list: EntityCells) {
    list.push(...getEntities(plugin, tag));
    for (const g of getGroups(plugin, tag)) {
        _getAllEntities(plugin, g.params?.values.tag, list);
    }
    return list;
}

export function getAllEntities(plugin: PluginContext, tag?: string) {
    return _getAllEntities(plugin, tag, []);
}

export function getAllFilteredEntities(plugin: PluginContext, tag: string, filter: string) {
    if (!filter) return getAllEntities(plugin, tag);
    const matcher = getFilterMatcher(filter);
    return getAllEntities(plugin, tag).filter(c => getEntityLabel(plugin, c).match(matcher) !== null);
}

export function getEveryEntity(plugin: PluginContext, filter?: string, tag?: string) {
    if (filter) {
        const matcher = getFilterMatcher(filter);
        return getAllEntities(plugin, tag).filter(c => getEntityLabel(plugin, c).match(matcher) !== null);
    } else {
        return getAllEntities(plugin, tag);
    }
}

export function getEntityLabel(plugin: PluginContext, cell: StateObjectCell) {
    return StateObjectRef.resolve(plugin.state.data, cell.transform.parent)?.obj?.label || 'Entity';
}

export function getCellDescription(cell: StateObjectCell) {
    // markdown style for description
    return '**' + cell?.obj?.label + '**\n\n' + cell?.obj?.description;
}

export function getEntityDescription(plugin: PluginContext, cell: StateObjectCell) {
    const s = StateObjectRef.resolve(plugin.state.data, cell.transform.parent);
    const d = getCellDescription(s!);
    return d;
}

export async function updateStyle(plugin: PluginContext, options: { ignoreLight: boolean, material: Material, celShaded: boolean, illustrative: boolean }) {
    const update = plugin.state.data.build();
    const { ignoreLight, material, celShaded, illustrative } = options;

    const entities = getAllEntities(plugin);

    for (let j = 0; j < entities.length; ++j) {
        update.to(entities[j]).update(old => {
            if (old.type) {
                const value = old.colorTheme.name === 'illustrative'
                    ? old.colorTheme.params.style.params.value
                    : old.colorTheme.params.value;
                const lightness = old.colorTheme.name === 'illustrative'
                    ? old.colorTheme.params.style.params.lightness
                    : old.colorTheme.params.lightness;
                if (illustrative) {
                    old.colorTheme = { name: 'illustrative', params: { style: { name: 'uniform', params: { value, lightness } } } };
                } else {
                    old.colorTheme = { name: 'uniform', params: { value, lightness } };
                }
                old.type.params.ignoreLight = ignoreLight;
                old.type.params.material = material;
                old.type.params.celShaded = celShaded;
            }
        });
    }

    await update.commit();
};

export async function updateColors(plugin: PluginContext, values: PD.Values, tag: string, filter: string) {
    const update = plugin.state.data.build();
    const { type, illustrative, value, shift, lightness, alpha, emissive } = values;
    if (type === 'group-generate' || type === 'group-uniform') {
        const leafGroups = getAllLeafGroups(plugin, tag);
        const rootLeafGroups = getRoots(plugin).filter(g => g.params?.values.tag === tag && getEntities(plugin, g.params?.values.tag).length > 0);
        const groups = [...leafGroups, ...rootLeafGroups];
        const baseColors = getDistinctBaseColors(groups.length, shift);

        for (let i = 0; i < groups.length; ++i) {
            const g = groups[i];
            const entities = getFilteredEntities(plugin, g.params?.values.tag, filter);
            let groupColors: Color[] = [];

            if (type === 'group-generate') {
                const c = g.params?.values.color;
                groupColors = getDistinctGroupColors(entities.length, baseColors[i], c.variability, c.shift);
            }

            for (let j = 0; j < entities.length; ++j) {
                const c = type === 'group-generate' ? groupColors[j] : baseColors[i];
                update.to(entities[j]).update(old => {
                    if (old.type) {
                        if (illustrative) {
                            old.colorTheme = { name: 'illustrative', params: { style: { name: 'uniform', params: { value: c, lightness: lightness } } } };
                        } else {
                            old.colorTheme = { name: 'uniform', params: { value: c, lightness: lightness } };
                        }
                        old.type.params.alpha = alpha;
                        old.type.params.xrayShaded = alpha < 1 ? 'inverted' : false;
                        old.type.params.emissive = emissive;
                    } else if (old.coloring) {
                        old.coloring.params.color = c;
                        old.coloring.params.lightness = lightness;
                        old.alpha = alpha;
                        old.xrayShaded = alpha < 1 ? true : false;
                        old.emissive = emissive;
                    }
                });
            }

            update.to(g.transform.ref).update(old => {
                old.color.type = type === 'group-generate' ? 'generate' : 'uniform';
                old.color.illustrative = illustrative;
                old.color.value = baseColors[i];
                old.color.lightness = lightness;
                old.color.alpha = alpha;
                old.color.emissive = emissive;
            });
        }
    } else if (type === 'generate' || type === 'uniform') {
        const entities = getAllFilteredEntities(plugin, tag, filter);
        let groupColors: Color[] = [];

        if (type === 'generate') {
            groupColors = getDistinctBaseColors(entities.length, shift);
        }

        for (let j = 0; j < entities.length; ++j) {
            const c = type === 'generate' ? groupColors[j] : value;
            update.to(entities[j]).update(old => {
                if (old.type) {
                    if (illustrative) {
                        old.colorTheme = { name: 'illustrative', params: { style: { name: 'uniform', params: { value: c, lightness: lightness } } } };
                    } else {
                        old.colorTheme = { name: 'uniform', params: { value: c, lightness: lightness } };
                    }
                    old.type.params.alpha = alpha;
                    old.type.params.xrayShaded = alpha < 1 ? 'inverted' : false;
                    old.type.params.emissive = emissive;
                } else if (old.coloring) {
                    old.coloring.params.color = c;
                    old.coloring.params.lightness = lightness;
                    old.alpha = alpha;
                    old.xrayShaded = alpha < 1 ? true : false;
                    old.emissive = emissive;
                }
            });
        }

        const others = getAllLeafGroups(plugin, tag);
        for (const o of others) {
            update.to(o).update(old => {
                old.color.type = type === 'generate' ? 'custom' : 'uniform';
                old.color.illustrative = illustrative;
                old.color.value = value;
                old.color.lightness = lightness;
                old.color.alpha = alpha;
                old.color.emissive = emissive;
            });
        }
    }

    await update.commit();
};

export function expandAllGroups(plugin: PluginContext) {
    for (const g of getAllGroups(plugin)) {
        if (g.state.isCollapsed) {
            plugin.state.data.updateCellState(g.transform.ref, { isCollapsed: false });
        }
    }
};

