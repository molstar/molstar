/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Representation, RepresentationContext, RepresentationParamsGetter, RepresentationProvider } from '../../../representation';
import { ParticleList, Particle, ParticleTarget } from '../../../../mol-model/particles/particle-list';
import { Theme } from '../../../../mol-theme/theme';
import { GraphicsRenderObject, getNextMaterialId } from '../../../../mol-gl/render-object';
import { Subject } from 'rxjs';
import { Task } from '../../../../mol-task';
import { Loci as ModelLoci, EmptyLoci, isEveryLoci, EveryLoci } from '../../../../mol-model/loci';
import { MarkerAction, MarkerActions } from '../../../../mol-util/marker-action';
import { LocationCallback, getQualityProps } from '../../../util';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { OrderedSet } from '../../../../mol-data/int';
import { ParticleTargetRepresentationParams, ParticleTargetRepresentationProps, createTargetVisual, makeFilteredParticleList, TargetVisual } from './common';

export { ParticleTargetRepresentationParams };
export type { ParticleTargetRepresentationProps };

export interface ParticleTargetRepresentation extends Representation<ParticleList, ParticleTargetRepresentationParams> { }

export type ParticleTargetRepresentationProvider<Id extends string = string> = RepresentationProvider<ParticleList, ParticleTargetRepresentationParams, Representation.State, Id>

export function ParticleTargetRepresentation(
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<ParticleList, ParticleTargetRepresentationParams>
): ParticleTargetRepresentation {
    const { webgl } = ctx;
    const updated = new Subject<number>();
    const geometryState = new Representation.GeometryState();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = Representation.createState();
    let _theme = Theme.createEmpty();

    let version = 0;
    let _particles: ParticleList | undefined;
    let _params: ParticleTargetRepresentationParams;
    let _props: ParticleTargetRepresentationProps;
    // Map from targetId → TargetVisual; re-used across updates.
    const targetVisuals = new Map<number, TargetVisual>();

    function createOrUpdate(props: Partial<ParticleTargetRepresentationProps> = {}, input?: ParticleList) {
        if (input) {
            _particles = input;
            _params = getParams(ctx, _particles);
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        const qualityProps = getQualityProps(Object.assign({}, _props, props));
        Object.assign(_props, props, qualityProps);

        return Task.create('ParticleTargetRepresentation', async runtime => {
            if (!_particles) return;

            const particles = _particles;
            const targets = Particle.getParticleTargets(particles) ?? new Map<number, ParticleTarget>();
            const { targets: targetIds } = particles;

            // Collect unique target IDs present in this particle list that have a target object.
            const presentTargets = new Set<number>();
            for (let i = 0; i < particles.count; i++) {
                const tid = targetIds[i];
                if (targets.has(tid)) presentTargets.add(tid);
            }

            renderObjects.length = 0;

            for (const targetId of presentTargets) {
                const target = targets.get(targetId)!;
                const { list: filtered, indices: filteredIndices } = makeFilteredParticleList(particles, targetId);
                if (filtered.count === 0) continue;

                let visual = targetVisuals.get(targetId);
                if (!visual) {
                    visual = createTargetVisual(targetId, getNextMaterialId(), webgl);
                    targetVisuals.set(targetId, visual);
                }

                await visual.createOrUpdate({ webgl, runtime }, _theme, _props, particles, OrderedSet.ofSortedArray(filteredIndices), target, filtered);

                if (visual.renderObject) {
                    renderObjects.push(visual.renderObject);
                    geometryState.add(visual.renderObject.id, visual.geometryVersion);
                }
            }

            // Remove visuals for targets no longer present.
            for (const [tid, visual] of targetVisuals) {
                if (!presentTargets.has(tid)) {
                    visual.destroy();
                    targetVisuals.delete(tid);
                }
            }

            geometryState.snapshot();
            updated.next(version++);
        });
    }

    function setState(state: Partial<Representation.State>) {
        Representation.updateState(_state, state);
        for (const visual of targetVisuals.values()) {
            if (!visual.renderObject) continue;
            if (state.visible !== undefined) visual.renderObject.state.visible = state.visible;
            if (state.alphaFactor !== undefined) visual.renderObject.state.alphaFactor = state.alphaFactor;
            if (state.pickable !== undefined) visual.renderObject.state.pickable = state.pickable;
        }
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function getLoci(pickingId: PickingId): ModelLoci {
        if (!_particles) return EmptyLoci;
        const { objectId, instanceId } = pickingId;
        for (const visual of targetVisuals.values()) {
            if (visual.renderObject && visual.renderObject.id === objectId) {
                // instanceId is local to this target visual; map to global particle index
                const indices = visual.particleIndices;
                if (!indices || instanceId < 0 || instanceId >= OrderedSet.size(indices)) return EmptyLoci;
                return Particle.Loci(_particles, OrderedSet.ofSingleton(OrderedSet.getAt(indices, instanceId)));
            }
        }
        return EmptyLoci;
    }

    function getAllLoci(): ModelLoci[] {
        if (!_particles) return [];
        return [Particle.Loci(_particles, OrderedSet.ofRange(0, _particles.count - 1))];
    }

    function eachLocation(cb: LocationCallback) {
        if (!_particles) return;
        const particles = _particles;
        const loc = Particle.Location(particles, 0);
        for (let i = 0; i < particles.count; i++) {
            loc.index = i;
            cb(loc, false);
        }
    }

    function mark(loci: ModelLoci, action: MarkerAction): boolean {
        if (!_particles) return false;
        if (!MarkerActions.is(_state.markerActions, action)) return false;
        if (!isEveryLoci(loci) && (!Particle.isLoci(loci) || loci.particles !== _particles)) return false;
        if (Particle.isLoci(loci) && Particle.lociSize(loci) === _particles.count) {
            // Change to `EveryLoci` to allow for downstream optimizations
            loci = EveryLoci;
        }
        let changed = false;
        for (const visual of targetVisuals.values()) {
            changed = visual.mark(loci, action) || changed;
        }
        return changed;
    }

    function destroy() {
        for (const visual of targetVisuals.values()) {
            visual.destroy();
        }
        targetVisuals.clear();
        renderObjects.length = 0;
    }

    return {
        label: 'Particle Target',
        get updated() { return updated; },
        get renderObjects(): ReadonlyArray<GraphicsRenderObject> { return renderObjects; },
        get props(): Readonly<ParticleTargetRepresentationProps> { return _props; },
        get params(): Readonly<ParticleTargetRepresentationParams> { return _params; },
        get state(): Readonly<Representation.State> { return _state; },
        get theme(): Readonly<Theme> { return _theme; },
        get groupCount() {
            return renderObjects.reduce((sum, ro) => sum + ((ro.values as any).uGroupCount?.ref?.value || 0), 0);
        },
        get geometryVersion() { return geometryState.version; },
        createOrUpdate,
        setState,
        setTheme,
        getLoci,
        getAllLoci,
        eachLocation,
        mark,
        destroy,
    };
}

export const ParticleTargetRepresentationProvider: ParticleTargetRepresentationProvider = {
    name: 'target',
    label: 'Target',
    description: 'Displays each particle as an instanced reference structure or shape.',
    factory: ParticleTargetRepresentation,
    getParams: (_ctx: RepresentationContext, _particles: ParticleList) => ParticleTargetRepresentationParams,
    defaultValues: PD.getDefaultValues(ParticleTargetRepresentationParams),
    defaultColorTheme: { name: 'particle-entity' },
    defaultSizeTheme: { name: 'uniform', props: { value: 1.6 } },
    isApplicable: (data: any) => (data as ParticleList).count > 0,
};
