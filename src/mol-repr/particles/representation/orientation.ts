/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';
import { Task } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../representation';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../mol-geo/geometry/spheres/spheres-builder';
import { MarkerAction, MarkerActions } from '../../../mol-util/marker-action';
import { ThemeRegistryContext, Theme } from '../../../mol-theme/theme';
import { Particle, ParticleList, getParticleTransforms } from '../../../mol-model/particles/particle-list';
import { GraphicsRenderObject, createRenderObject, getNextMaterialId } from '../../../mol-gl/render-object';
import { createTransform } from '../../../mol-geo/geometry/transform-data';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { UniformColorTheme } from '../../../mol-theme/color/uniform';
import { UniformSizeTheme } from '../../../mol-theme/size/uniform';
import { ColorTheme } from '../../../mol-theme/color';
import { SizeTheme } from '../../../mol-theme/size';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Loci as ModelLoci, EmptyLoci, isEveryLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { Visual } from '../../visual';
import { ParticleRepresentation, ParticleRepresentationProvider } from '../representation';

const AxisColorByGroup = [ColorNames.red, ColorNames.green, ColorNames.blue] as const;

const AxisLengthOptions = { min: 0.1, max: 1000, step: 0.1 } as const;
const PointSizeOptions = { min: 0.1, max: 100, step: 0.1 } as const;

export const BaseOrientationParticlesParams = {
    ...Spheres.Params,
    ...Lines.Params,
    pointSize: PD.Numeric(1, PointSizeOptions, { description: 'Radius used for the particle position marker.' }),
    axisLength: PD.Numeric(10, AxisLengthOptions, { description: 'Length of the particle orientation axes.' }),
    positionColor: PD.Color(ColorNames.white),
    xColor: PD.Color(ColorNames.red),
    yColor: PD.Color(ColorNames.green),
    zColor: PD.Color(ColorNames.blue),
};

const OrientationParticlesVisualKinds = ['position', 'orientation'] as const;
type OrientationParticlesVisualKind = typeof OrientationParticlesVisualKinds[number]

export const OrientationParticlesParams = {
    ...BaseOrientationParticlesParams,
    visuals: PD.MultiSelect(['position', 'orientation'] as OrientationParticlesVisualKind[], PD.arrayToOptions(OrientationParticlesVisualKinds as unknown as OrientationParticlesVisualKind[])),
};

export type OrientationParticlesParams = typeof OrientationParticlesParams
export type OrientationParticlesProps = PD.Values<OrientationParticlesParams>

export function getOrientationParticlesParams(ctx: ThemeRegistryContext, data: ParticleList) {
    const hasRotations = !!data.rotations;
    const visualKinds: OrientationParticlesVisualKind[] = hasRotations
        ? ['position', 'orientation']
        : ['position'];
    const defaultVisuals: OrientationParticlesVisualKind[] = hasRotations
        ? ['position', 'orientation']
        : ['position'];
    return {
        ...BaseOrientationParticlesParams,
        visuals: PD.MultiSelect(defaultVisuals, PD.arrayToOptions(visualKinds)),
    };
}

function getAxisColor(props: OrientationParticlesProps, groupId: number) {
    switch (groupId) {
        case 0: return props.xColor;
        case 1: return props.yColor;
        case 2: return props.zColor;
        default: return AxisColorByGroup[groupId % AxisColorByGroup.length];
    }
}

function createPositionGeometry(spheres?: Spheres) {
    const builder = SpheresBuilder.create(1, 1, spheres);
    builder.add(0, 0, 0, 0);
    return builder.getSpheres();
}

function createOrientationGeometry(axisLength: number, lines?: Lines) {
    const builder = LinesBuilder.create(3, 3, lines);
    builder.add(0, 0, 0, axisLength, 0, 0, 0);
    builder.add(0, 0, 0, 0, axisLength, 0, 1);
    builder.add(0, 0, 0, 0, 0, axisLength, 2);
    return builder.getLines();
}

function particleTransformArray(data: ParticleList): Float32Array {
    const transforms = getParticleTransforms(data);
    const arr = new Float32Array(transforms.length * 16);
    for (let i = 0; i < transforms.length; ++i) arr.set(transforms[i], i * 16);
    return arr;
}

type CarriedLocation = Particle.Location & { groupIndex: number }

function createCarriedLocation(data: ParticleList): CarriedLocation {
    return Object.assign(Particle.Location(data, 0), { groupIndex: 0 });
}

function createParticleLocationIterator(groupCount: number, instanceCount: number, location: CarriedLocation) {
    return LocationIterator(groupCount, instanceCount, 1, (groupIndex, instanceIndex) => {
        location.index = instanceIndex;
        location.groupIndex = groupIndex;
        return location;
    });
}

function axisColorTheme(props: OrientationParticlesProps): ColorTheme<{}> {
    const c0 = getAxisColor(props, 0);
    const c1 = getAxisColor(props, 1);
    const c2 = getAxisColor(props, 2);
    const theme: ColorTheme<{}> = {
        factory: (() => theme) as any,
        granularity: 'group',
        color: (location) => {
            const g = (location as any).groupIndex as number | undefined;
            if (g === 0) return c0;
            if (g === 1) return c1;
            if (g === 2) return c2;
            return Color(0xCCCCCC);
        },
        props: {},
    };
    return theme;
}

function uniformColorTheme(color: Color): ColorTheme<any> {
    return UniformColorTheme({} as any, { value: color, lightness: 0, saturation: 0 });
}

function uniformSize(value: number) {
    return UniformSizeTheme({} as any, { value });
}

export type OrientationParticlesRepresentation = ParticleRepresentation<OrientationParticlesParams>

interface VisualEntry {
    kind: OrientationParticlesVisualKind
    renderObject: GraphicsRenderObject
    geometryVersion: number
}

export function OrientationParticlesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, OrientationParticlesParams>): OrientationParticlesRepresentation {
    let version = 0;
    const updated = new Subject<number>();
    const geometryState = new Representation.GeometryState();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = Representation.createState();
    let _theme = Theme.createEmpty();
    _state.markerActions = MarkerActions.Highlighting;

    let _data: ParticleList | undefined = undefined;
    let _params: OrientationParticlesParams = OrientationParticlesParams;
    let _props: OrientationParticlesProps = PD.getDefaultValues(OrientationParticlesParams);
    const entries = new Map<OrientationParticlesVisualKind, VisualEntry>();
    let dirty = true;

    function disposeEntries() {
        for (const e of entries.values()) {
            e.renderObject.state.disposed = true;
        }
        entries.clear();
    }

    function buildPositionEntry(data: ParticleList, props: OrientationParticlesProps): VisualEntry | undefined {
        const spheres = createPositionGeometry();
        const transformArr = particleTransformArray(data);
        const instanceCount = transformArr.length / 16 | 0;
        if (instanceCount === 0) return undefined;
        const transform = createTransform(transformArr, instanceCount, spheres.boundingSphere, props.cellSize, props.batchSize);
        const carry = createCarriedLocation(data);
        const locationIt = createParticleLocationIterator(1, instanceCount, carry);
        const colorFromTheme = _theme.color.factory !== ColorTheme.EmptyFactory;
        const sizeFromTheme = _theme.size.factory !== SizeTheme.EmptyFactory;
        const theme: Theme = {
            color: colorFromTheme ? _theme.color : uniformColorTheme(props.positionColor),
            size: sizeFromTheme ? _theme.size : uniformSize(props.pointSize),
        };
        const values = Spheres.Utils.createValues(spheres, transform, locationIt, theme, props);
        const state = Spheres.Utils.createRenderableState(props);
        const renderObject = createRenderObject('spheres', values, state, materialId);
        return { kind: 'position', renderObject, geometryVersion: 0 };
    }

    function buildOrientationEntry(data: ParticleList, props: OrientationParticlesProps): VisualEntry | undefined {
        const lines = createOrientationGeometry(props.axisLength);
        const transformArr = particleTransformArray(data);
        const instanceCount = transformArr.length / 16 | 0;
        if (instanceCount === 0) return undefined;
        const transform = createTransform(transformArr, instanceCount, lines.boundingSphere, props.cellSize, props.batchSize);
        const carry = createCarriedLocation(data);
        const locationIt = createParticleLocationIterator(3, instanceCount, carry);
        const theme: Theme = {
            color: axisColorTheme(props),
            size: uniformSize(1),
        };
        const values = Lines.Utils.createValues(lines, transform, locationIt, theme, props);
        const state = Lines.Utils.createRenderableState(props);
        const renderObject = createRenderObject('lines', values, state, materialId);
        return { kind: 'orientation', renderObject, geometryVersion: 0 };
    }

    function rebuild() {
        disposeEntries();
        renderObjects.length = 0;
        if (!_data) return;
        const enabled = new Set(_props.visuals as OrientationParticlesVisualKind[]);
        if (enabled.has('position')) {
            const e = buildPositionEntry(_data, _props);
            if (e) entries.set('position', e);
        }
        if (enabled.has('orientation') && _data.rotations) {
            const e = buildOrientationEntry(_data, _props);
            if (e) entries.set('orientation', e);
        }
        for (const e of entries.values()) {
            renderObjects.push(e.renderObject);
            geometryState.add(e.renderObject.id, e.geometryVersion);
        }
        geometryState.snapshot();
        applyCurrentState();
    }

    function applyCurrentState() {
        for (const e of entries.values()) {
            Visual.setVisibility(e.renderObject, _state.visible);
            Visual.setAlphaFactor(e.renderObject, _state.alphaFactor);
            Visual.setPickable(e.renderObject, _state.pickable);
            Visual.setColorOnly(e.renderObject, _state.colorOnly);
            Visual.setThemeStrength(e.renderObject, _state.themeStrength);
            Visual.setTransform(e.renderObject, _state.transform);
        }
    }

    function createOrUpdate(props: Partial<OrientationParticlesProps> = {}, data?: ParticleList) {
        if (data && data !== _data) {
            _params = getParams(ctx, data);
            _data = data;
            dirty = true;
        }
        const next = { ..._props, ...props } as OrientationParticlesProps;
        // Any prop change triggers a rebuild; the geometry is trivial (1 sphere + 3 lines).
        if (JSON.stringify(next) !== JSON.stringify(_props)) dirty = true;
        _props = next;

        return Task.create('Creating or updating Orientation Particles representation', async _runtime => {
            if (dirty && _data) {
                rebuild();
                dirty = false;
            }
            updated.next(version++);
        });
    }

    function lociApplyFor(blockSize: number): Visual.LociApply {
        return (l, apply) => {
            if (Particle.isLoci(l)) {
                let c = false;
                OrderedSet.forEach(l.indices, p => {
                    const start = p * blockSize;
                    if (apply(Interval.ofBounds(start, start + blockSize))) c = true;
                });
                return c;
            }
            return false;
        };
    }

    function mark(loci: ModelLoci, action: MarkerAction) {
        if (entries.size === 0) return false;
        if (!MarkerActions.is(_state.markerActions, action)) return false;

        let remapped: ModelLoci = loci;
        if (Particle.isLoci(loci) && _data) {
            remapped = Particle.remapLoci(loci, _data);
            if (Particle.isLociEmpty(remapped as Particle.Loci)) return false;
        } else if (!isEveryLoci(loci)) {
            return false;
        }

        let changed = false;
        for (const e of entries.values()) {
            const values: any = e.renderObject.values;
            const instanceGranularity: boolean = values.instanceGranularity.ref.value;
            const groupCount: number = values.uGroupCount.ref.value;
            const blockSize = instanceGranularity ? 1 : groupCount;
            if (Visual.mark(e.renderObject, remapped, action, lociApplyFor(blockSize))) changed = true;
        }
        return changed;
    }

    function setState(state: Partial<Representation.State>) {
        Representation.updateState(_state, state);
        applyCurrentState();
    }

    function setTheme(t: Theme) {
        _theme = t;
        dirty = true;
    }

    function destroy() {
        disposeEntries();
        renderObjects.length = 0;
        _data = undefined;
    }

    return {
        label: 'Orientation Particles',
        updated,
        get groupCount() {
            let n = 0;
            for (const e of entries.values()) {
                const v: any = e.renderObject.values;
                n += v.uGroupCount.ref.value * v.instanceCount.ref.value;
            }
            return n;
        },
        get renderObjects() { return renderObjects; },
        get props() { return _props; },
        get params() { return _params; },
        get state() { return _state; },
        get theme() { return _theme; },
        get geometryVersion() { return geometryState.version; },
        createOrUpdate,
        setState,
        setTheme,
        getLoci: (pickingId: PickingId): ModelLoci => {
            if (!_data) return EmptyLoci;
            for (const e of entries.values()) {
                if (e.renderObject.id !== pickingId.objectId) continue;
                const idx = pickingId.instanceId;
                if (idx < 0 || idx >= _data.count) return EmptyLoci;
                return Particle.Loci(_data, OrderedSet.ofSingleton(idx));
            }
            return EmptyLoci;
        },
        getAllLoci: (): ModelLoci[] => {
            if (!_data) return [];
            const { count } = _data;
            if (count === 0) return [];
            return [Particle.Loci(_data, OrderedSet.ofBounds(0, count))];
        },
        eachLocation: (cb) => {
            if (!_data) return;
            const { count } = _data;
            const location = Particle.Location(_data, 0);
            for (let i = 0; i < count; ++i) {
                location.index = i;
                cb(location, false);
            }
        },
        mark,
        destroy,
    };
}

export const OrientationParticlesRepresentationProvider = ParticleRepresentationProvider({
    name: 'orientation',
    label: 'Orientation',
    description: 'Displays particles as position markers and orientation axes.',
    factory: OrientationParticlesRepresentation,
    getParams: getOrientationParticlesParams,
    defaultValues: PD.getDefaultValues(OrientationParticlesParams),
    defaultColorTheme: { name: 'particle-index' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (data: ParticleList) => data.count > 0,
});
