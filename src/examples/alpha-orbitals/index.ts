import { BehaviorSubject } from 'rxjs';
import { debounceTime, skip } from 'rxjs/operators';
/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Basis, computeIsocontourValues } from '../../extensions/alpha-orbitals/cubes';
import { SphericalBasisOrder } from '../../extensions/alpha-orbitals/orbitals';
import { createPluginAsync, DefaultPluginSpec } from '../../mol-plugin';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition } from '../../mol-util/param-definition';
import { mountControls } from './controls';
import { DemoMoleculeSDF, DemoOrbitals } from './example-data';
import './index.html';
import { CreateOrbitalVolume, StaticBasisAndOrbitals } from './transforms';
require('mol-plugin-ui/skin/light.scss');

interface DemoInput {
    moleculeSdf: string,
    basis: Basis,
    order: SphericalBasisOrder,
    orbitals: {
        energy: number,
        alpha: number[]
    }[]
}

interface Params {
    orbitalIndex: number,
    isoValue: number,
    staticIsovalues: boolean
}

export class AlphaOrbitalsExample {
    plugin: PluginContext;

    async init(target: string | HTMLElement) {
        this.plugin = await createPluginAsync(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: { left: 'none', right: 'none', top: 'none', bottom: 'none' },
            },
            config: [
                [PluginConfig.Viewport.ShowExpand, false],
                [PluginConfig.Viewport.ShowControls, false],
                [PluginConfig.Viewport.ShowSelectionMode, false],
                [PluginConfig.Viewport.ShowAnimation, false],
            ]
        });

        this.plugin.managers.interactivity.setProps({ granularity: 'element' });

        this.load({
            moleculeSdf: DemoMoleculeSDF,
            ...DemoOrbitals
        });

        mountControls(this, document.getElementById('controls')!);
    }

    readonly params = new BehaviorSubject<ParamDefinition.For<Params>>({ } as any);
    readonly state = new BehaviorSubject<Params>({ orbitalIndex: 32, isoValue: 1, staticIsovalues: true });

    private volume?: StateObjectSelector<PluginStateObject.Volume.Data, typeof CreateOrbitalVolume>;
    private positive?: StateObjectSelector<PluginStateObject.Volume.Representation3D, typeof StateTransforms.Representation.VolumeRepresentation3D>;
    private negative?: StateObjectSelector<PluginStateObject.Volume.Representation3D, typeof StateTransforms.Representation.VolumeRepresentation3D>;
    private isovalues: { negative?: number, positive?: number } = { negative: void 0, positive: void 0 }
    private currentParams: Params = { ...this.state.value };

    private async setIndex() {
        if (!this.volume?.isOk) return;
        const state = this.state.value;
        await this.plugin.build().to(this.volume).update(CreateOrbitalVolume, () => ({ index: state.orbitalIndex })).commit();
        this.currentParams.orbitalIndex = this.state.value.orbitalIndex;
        this.currentParams.staticIsovalues = this.state.value.staticIsovalues;
        this.currentParams.isoValue = this.state.value.isoValue;

        if (!state.staticIsovalues) {
            this.isovalues = computeIsocontourValues(this.volume.data!.grid.cells.data as any, 0.85);
            await this.setIsovalue();
        }
    }

    private setIsovalue() {
        const { positive, negative } = this.isovalues;
        this.currentParams.isoValue = this.state.value.isoValue;
        this.currentParams.staticIsovalues = this.state.value.staticIsovalues;
        const update = this.plugin.build();
        update.to(this.positive!).update(this.volumeParams(positive, ColorNames.blue));
        update.to(this.negative!).update(this.volumeParams(negative, ColorNames.red));
        return update.commit();
    }

    private volumeParams(value: number | undefined, color: Color) {
        return createVolumeRepresentationParams(this.plugin, this.volume!.data!, {
            // type: 'isosurface',
            // typeParams: { isoValue: { kind: 'absolute', absoluteValue: (value ?? 1000) * this.state.value.isoValue } },
            // color: 'uniform',
            // colorParams: { value: ColorNames.blue }
            type: 'direct-volume',
            typeParams: {
                renderMode: {
                    name: 'isosurface',
                    params: { isoValue: { kind: 'absolute', absoluteValue: (value ?? 1000) * this.state.value.isoValue }, singleLayer: true }
                }
            },
            color: 'uniform',
            colorParams: { value: color }
        });
    }

    async load(input: DemoInput) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.rawData({ data: input.moleculeSdf }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mol');
        const model = await this.plugin.builders.structure.createModel(trajectory);
        const structure = await this.plugin.builders.structure.createStructure(model);

        const all = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'all');
        if (all) await this.plugin.builders.structure.representation.addRepresentation(all, { type: 'ball-and-stick', color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } } });

        const state = this.state.value;

        this.volume = await this.plugin.build().toRoot()
            .apply(StaticBasisAndOrbitals, { basis: input.basis, order: input.order, orbitals: input.orbitals })
            .apply(CreateOrbitalVolume, { index: state.orbitalIndex })
            .commit();

        if (!this.volume.isOk) {
            this.volume = void 0;
            return;
        }

        // TODO: the isovalues are being computed twice. Need to add more flexible support to Volume object
        //       for controlling them
        this.isovalues = computeIsocontourValues(this.volume.data!.grid.cells.data as any, 0.85);
        const { positive, negative } = this.isovalues;

        const repr = this.plugin.build().to(this.volume);

        this.positive = repr.apply(StateTransforms.Representation.VolumeRepresentation3D, this.volumeParams(positive, ColorNames.blue)).selector;
        this.negative = repr.apply(StateTransforms.Representation.VolumeRepresentation3D, this.volumeParams(negative, ColorNames.red)).selector;

        await repr.commit();

        this.params.next({
            orbitalIndex: ParamDefinition.Numeric(this.currentParams.orbitalIndex, { min: 0, max: input.orbitals.length - 1 }, { immediateUpdate: true, isEssential: true }),
            isoValue: ParamDefinition.Numeric(this.currentParams.isoValue, { min: 0.5, max: 3, step: 0.1 }, { immediateUpdate: true, isEssential: false }),
            staticIsovalues: ParamDefinition.Boolean(this.currentParams.staticIsovalues)
        });

        this.state.pipe(skip(1), debounceTime(1000 / 24)).subscribe(async params => {
            if (params.orbitalIndex !== this.currentParams.orbitalIndex || params.staticIsovalues !== this.currentParams.staticIsovalues) {
                this.setIndex();
            } else if (params.isoValue !== this.currentParams.isoValue) {
                this.setIsovalue();
            }
        });
    }
}

(window as any).AlphaOrbitalsExample = new AlphaOrbitalsExample();