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
import { stringToWords } from '../../../mol-util/string';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ParamMapping } from '../../../mol-util/param-mapping';
import { EntityNode } from '../ui/entities';
import { DistinctColorsProps, distinctColors } from '../../../mol-util/color/distinct';
import { Sphere3D } from '../../../mol-math/geometry';
import { Hcl } from '../../../mol-util/color/spaces/hcl';

export function getDistinctGroupColors(count: number, color: Color, props?: Partial<DistinctColorsProps>) {
    const hcl = Hcl.fromColor(Hcl(), color);
    const hue = color === 0
        ? [1, 360] as [number, number]
        : [Math.max(1, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number];
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
    value: PD.Color(Color(0xFFFFFF)),
    type: PD.Select('generate', PD.arrayToOptions(['generate', 'uniform', 'custom'])),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
};
export type ColorProps = PD.Values<typeof ColorParams>

export const ColorValueParam = PD.Color(Color(0xFFFFFF));

export const RootParams = {
    type: PD.Select('generate', PD.arrayToOptions(['generate', 'uniform', 'custom'])),
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

export function getClipProps(values: SimpleClipProps, boundingSphere: Sphere3D): Clip.Props {
    const { center, radius } = boundingSphere;

    const position = Vec3.clone(center);
    Vec3.add(position, position, Vec3.create(
        radius * values.position.x / 100,
        radius * values.position.y / 100,
        radius * values.position.z / 100
    ));

    const scale = Vec3.create(values.scale.x, values.scale.y, values.scale.z);
    Vec3.scale(scale, scale, 2 * radius / 100);

    return {
        variant: 'instance',
        objects: [{
            type: values.type,
            invert: values.invert,
            position,
            scale,
            rotation: values.rotation
        }],
    };
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

            const clip = getClipProps(s, node.plugin.canvas3d!.boundingSphere);
            props.variant = clip.variant;
            props.objects = clip.objects;
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
    color: PD.Group(ColorParams),
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

//

const MesoscaleStateParams = {
    filter: PD.Value<string>('', { isHidden: true }),
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

        await ctx.state.data.build().toRoot().apply(MesoscaleStateTransform, {}).commit();
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
