/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ParamMapping } from '../../../mol-util/param-mapping';
import { EntityNode } from '../ui/entities';
import { DistinctColorsProps, distinctColors } from '../../../mol-util/color/distinct';
import { Sphere3D } from '../../../mol-math/geometry';
import { Hcl } from '../../../mol-util/color/spaces/hcl';
import { StateObjectCell, StateObjectRef, StateSelection } from '../../../mol-state';
import { StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { SpacefillRepresentationProvider } from '../../../mol-repr/structure/representation/spacefill';
import { assertUnreachable } from '../../../mol-util/type-helpers';
import { MesoscaleExplorerState } from '../app';

export function getDistinctGroupColors(count: number, color: Color, variablity: number, props?: Partial<DistinctColorsProps>) {
    const hcl = Hcl.fromColor(Hcl(), color);
    const hue = color === 0
        ? [1, 360] as [number, number]
        : [Math.max(1, hcl[0] - variablity), Math.min(360, hcl[0] + variablity)] as [number, number];
    return distinctColors(count, {
        hue,
        chroma: [30, 80],
        luminance: [15, 85],
        clusteringStepCount: 50,
        minSampleCount: 800,
        ...props,
    });
}

export function getDistinctBaseColors(count: number, props?: Partial<DistinctColorsProps>) {
    return distinctColors(count, {
        hue: [1, 360],
        chroma: [40, 70],
        luminance: [15, 85],
        clusteringStepCount: 50,
        minSampleCount: 800,
        ...props,
    });
}

export const ColorParams = {
    type: PD.Select('generate', PD.arrayToOptions(['generate', 'uniform', 'custom'])),
    value: PD.Color(Color(0xFFFFFF), { hideIf: p => p.type === 'custom' }),
    variablity: PD.Numeric(35, { min: 1, max: 360, step: 1 }, { hideIf: p => p.type !== 'generate' }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }, { hideIf: p => p.type === 'custom' }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
};
export type ColorProps = PD.Values<typeof ColorParams>

export const ColorValueParam = PD.Color(Color(0xFFFFFF));

export const RootParams = {
    type: PD.Select('custom', PD.arrayToOptions(['group-generate', 'group-uniform', 'generate', 'uniform', 'custom'])),
    value: PD.Color(Color(0xFFFFFF), { hideIf: p => p.type !== 'uniform' }),
    variablity: PD.Numeric(35, { min: 1, max: 360, step: 1 }, { hideIf: p => p.type !== 'group-generate' }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }, { hideIf: p => p.type === 'custom' }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { hideIf: p => p.type === 'custom' }),
};

export const LightnessParams = {
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
export const DimLightness = 6;

export const OpacityParams = {
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
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
        rotation: values.rotation
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
    hidden: PD.Boolean(false),
    color: PD.Group(RootParams),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
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
            return new MesoscaleGroupObject({}, { label: params.label });
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
                { minDistance: 300, maxDistance: 2000, overlap: 40, stride: 40, scaleBias: 3 },
                { minDistance: 2000, maxDistance: 6000, overlap: 200, stride: 150, scaleBias: 2.5 },
                { minDistance: 6000, maxDistance: 10000000, overlap: 600, stride: 300, scaleBias: 2 },
            ];
        case 'balanced':
            return [
                { minDistance: 1, maxDistance: 500, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 500, maxDistance: 2000, overlap: 50, stride: 15, scaleBias: 3 },
                { minDistance: 2000, maxDistance: 6000, overlap: 200, stride: 70, scaleBias: 2.5 },
                { minDistance: 6000, maxDistance: 10000000, overlap: 600, stride: 200, scaleBias: 2 },
            ];
        case 'quality':
            return [
                { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1, scaleBias: 1 },
                { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10, scaleBias: 3 },
                { minDistance: 4000, maxDistance: 10000, overlap: 500, stride: 50, scaleBias: 2.5 },
                { minDistance: 10000, maxDistance: 10000000, overlap: 1000, stride: 200, scaleBias: 2 },
            ];
        default:
            assertUnreachable(graphicsMode);
    }
}

export type GraphicsMode = 'quality' | 'balanced' | 'performance' | 'custom';

export function getGraphicsModeProps(graphicsMode: Exclude<GraphicsMode, 'custom'>) {
    return {
        lodLevels: getLodLevels(graphicsMode),
        approximate: graphicsMode !== 'quality',
        alphaThickness: graphicsMode === 'performance' ? 15 : 12,
    };
}

//

export const MesoscaleStateParams = {
    filter: PD.Value<string>('', { isHidden: true }),
    graphics: PD.Select('quality', PD.arrayToOptions(['quality', 'balanced', 'performance', 'custom'] as GraphicsMode[])),
};

class MesoscaleStateObject extends PSO.Create<MesoscaleState>({ name: 'Mesoscale State', typeClass: 'Object' }) { }

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

function getMesoscaleStateCell(ctx: PluginContext) {
    const cell = ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
    if (!cell) throw new Error('MesoscaleState not initialized');

    return cell;
}

export { MesoscaleState };
type MesoscaleState = PD.Values<typeof MesoscaleStateParams>;
const MesoscaleState = {
    async init(ctx: PluginContext) {
        const cell = ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
        if (cell) throw new Error('MesoscaleState already initialized');

        const customState = ctx.customState as MesoscaleExplorerState;
        await ctx.state.data.build().toRoot().apply(MesoscaleStateTransform, {
            filter: '',
            graphics: customState.graphicsMode,
        }).commit();
    },
    get(ctx: PluginContext): MesoscaleState {
        const cell = getMesoscaleStateCell(ctx);
        return cell.obj!.data;
    },
    async set(ctx: PluginContext, props: Partial<MesoscaleState>) {
        const cell = getMesoscaleStateCell(ctx);
        await ctx.state.data.build().to(cell).update(MesoscaleStateTransform, old => Object.assign(old, props)).commit();
    },
    ref(ctx: PluginContext): string {
        const cell = getMesoscaleStateCell(ctx);
        return cell.transform.ref;
    },
    has(ctx: PluginContext): boolean {
        return !!ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
    },
};

//

export function getRoots(plugin: PluginContext) {
    return plugin.state.data.select(StateSelection.Generators.rootsOfType(MesoscaleGroupObject));
}

export function getGroups(plugin: PluginContext, tag?: string) {
    const selector = tag !== undefined
        ? StateSelection.Generators.ofTransformer(MesoscaleGroup).withTag(tag)
        : StateSelection.Generators.ofTransformer(MesoscaleGroup);
    return plugin.state.data.select(selector);
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
        return plugin.state.data.select(StateSelection.Generators.ofTransformer(StructureRepresentation3D).withTag(g.params?.values.tag)).length > 0;
    });
}

export function getEntities(plugin: PluginContext, tag?: string) {
    const selector = tag !== undefined
        ? StateSelection.Generators.ofTransformer(StructureRepresentation3D).withTag(tag)
        : StateSelection.Generators.ofTransformer(StructureRepresentation3D);
    return plugin.state.data.select(selector).filter(c => c.obj!.data.sourceData.elementCount > 0);
}

export function getFilteredEntities(plugin: PluginContext, tag: string, filter: string) {
    const reFilter = new RegExp(escapeRegExp(filter), 'gi');
    return getEntities(plugin, tag).filter(c => getEntityLabel(plugin, c).match(reFilter) !== null);
}

function _getAllEntities(plugin: PluginContext, tag: string | undefined, list: StateObjectCell[]) {
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
    const reFilter = new RegExp(escapeRegExp(filter), 'gi');
    return getAllEntities(plugin, tag).filter(c => getEntityLabel(plugin, c).match(reFilter) !== null);
}

export function getEntityLabel(plugin: PluginContext, cell: StateObjectCell) {
    return StateObjectRef.resolve(plugin.state.data, cell.transform.parent)?.obj?.label || 'Entity';
}

//

export async function updateColors(plugin: PluginContext, values: PD.Values, tag: string, filter: string) {
    const update = plugin.state.data.build();
    const { type, value, lightness, alpha } = values;

    if (type === 'group-generate' || type === 'group-uniform') {
        const groups = getAllLeafGroups(plugin, tag);
        const baseColors = getDistinctBaseColors(groups.length);

        for (let i = 0; i < groups.length; ++i) {
            const g = groups[i];
            const entities = getFilteredEntities(plugin, g.params?.values.tag, filter);
            let groupColors: Color[] = [];

            if (type === 'group-generate') {
                groupColors = getDistinctGroupColors(entities.length, baseColors[i], g.params?.values.color.variablity);
            }

            for (let j = 0; j < entities.length; ++j) {
                const c = type === 'group-generate' ? groupColors[j] : baseColors[i];
                update.to(entities[j]).update(old => {
                    old.colorTheme.params.value = c;
                    old.colorTheme.params.lightness = lightness;
                    old.type.params.alpha = alpha;
                    old.type.params.xrayShaded = alpha < 1 ? 'inverted' : false;
                });
            }

            update.to(g.transform.ref).update(old => {
                old.color.type = type === 'group-generate' ? 'generate' : 'uniform';
                old.color.value = baseColors[i];
                old.color.lightness = lightness;
                old.color.alpha = alpha;
            });
        }
    } else if (type === 'generate' || type === 'uniform') {
        const entities = getAllFilteredEntities(plugin, tag, filter);
        let groupColors: Color[] = [];

        if (type === 'generate') {
            groupColors = getDistinctBaseColors(entities.length);
        }

        for (let j = 0; j < entities.length; ++j) {
            const c = type === 'generate' ? groupColors[j] : value;
            update.to(entities[j]).update(old => {
                old.colorTheme.params.value = c;
                old.colorTheme.params.lightness = lightness;
                old.type.params.alpha = alpha;
                old.type.params.xrayShaded = alpha < 1 ? 'inverted' : false;
            });
        }

        const others = getAllLeafGroups(plugin, tag);
        for (const o of others) {
            update.to(o).update(old => {
                old.color.type = type === 'generate' ? 'custom' : 'uniform';
                old.color.value = value;
                old.color.lightness = lightness;
                old.color.alpha = alpha;
            });
        }
    }

    await update.commit();
};
