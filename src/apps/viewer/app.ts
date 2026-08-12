/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { MVSData } from '../../extensions/mvs';
import { loadMVS, MVSLoadOptions } from '../../extensions/mvs/load';
import { applyStructureInteractivity, StructureInteractivityOptions } from '../../extensions/plugin/interactivity';
import * as loaders from '../../extensions/plugin/loaders';
import { Volume } from '../../mol-model/volume';
import { PluginComponent } from '../../mol-plugin-state/component';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { BuildInVolumeFormat } from '../../mol-plugin-state/formats/volume';
import { createPluginUI } from '../../mol-plugin-ui';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { PluginState } from '../../mol-plugin/state';
import { Color } from '../../mol-util/color';
import { decodeColor } from '../../mol-util/color/utils';
import '../../mol-util/polyfill';
import { DefaultViewerOptions, ViewerOptions } from './options';
import { createViewerSpec } from './plugin-spec';
import { ViewerAutoPreset } from './presets';

export { PLUGIN_VERSION as version } from '../../mol-plugin/version';
export { consoleStats, isDebugMode, isProductionMode, isTimingMode, setDebugMode, setProductionMode, setTimingMode } from '../../mol-util/debug';

// re-export for backwards compatibility, but these should ideally be imported from the plugin extension directly
// TODO: consider removing these in v6.0
export type { LoadStructureOptions, LoadTrajectoryParams, VolumeIsovalueInfo } from '../../extensions/plugin/loaders';

export class Viewer {
    private _events = new PluginComponent();
    public readonly plugin: PluginUIContext;

    constructor(plugin: PluginUIContext) {
        this.plugin = plugin;
    }

    static async create(elementOrId: string | HTMLElement, options: Partial<ViewerOptions> = {}) {
        const spec = createViewerSpec(options);

        const element = typeof elementOrId === 'string'
            ? document.getElementById(elementOrId)
            : elementOrId;
        if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);
        const plugin = await createPluginUI({
            target: element,
            spec,
            render: renderReact18,
            onBeforeUIRender: plugin => {
                // the preset needs to be added before the UI renders otherwise
                // "Download Structure" wont be able to pick it up
                plugin.builders.structure.representation.registerPreset(ViewerAutoPreset);
            }
        });

        plugin.canvas3d?.setProps({ illumination: { enabled: options.illumination ?? DefaultViewerOptions.illumination } });
        if (options.viewportBackgroundColor ?? DefaultViewerOptions.viewportBackgroundColor) {
            const backgroundColor = decodeColor(options.viewportBackgroundColor ?? DefaultViewerOptions.viewportBackgroundColor);
            if (typeof backgroundColor === 'number') {
                plugin.canvas3d?.setProps({ renderer: { backgroundColor } });
            }
        }
        setTimeout(async () => {
            const mvs = builderDemo2();
            console.log(MVSData.toMVSJ(mvs))
            await loadMVS(plugin, mvs);
        }, 1000)
        return new Viewer(plugin);
    }

    /**
     * Allows subscribing to rxjs observables in the context of the viewer.
     * All subscriptions will be disposed of when the viewer is destroyed.
     */
    subscribe = this._events.subscribe.bind(this._events);

    setRemoteSnapshot(id: string) {
        return loaders.setRemoteSnapshot(this.plugin, id);
    }

    loadSnapshotFromUrl(url: string, type: PluginState.SnapshotType) {
        return loaders.loadSnapshotFromUrl(this.plugin, url, type);
    }

    loadStructureFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: loaders.LoadStructureOptions & { label?: string }) {
        return loaders.loadStructureFromUrl(this.plugin, url, format, isBinary, options);
    }

    loadAllModelsOrAssemblyFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: loaders.LoadStructureOptions) {
        return loaders.loadAllModelsOrAssemblyFromUrl(this.plugin, url, format, isBinary, options);
    }

    loadStructureFromData(data: string | number[], format: BuiltInTrajectoryFormat, options?: { dataLabel?: string }) {
        return loaders.loadStructureFromData(this.plugin, data, format, options);
    }

    loadPdb(pdb: string, options?: loaders.LoadStructureOptions) {
        return loaders.loadPdb(this.plugin, pdb, options);
    }

    /**
     * @deprecated Scheduled for removal in v5. Use {@link loadPdbIhm | loadPdbIhm(pdbIhm: string)} instead.
     */
    loadPdbDev(pdbDev: string) {
        return this.loadPdbIhm(pdbDev);
    }

    loadPdbIhm(pdbIhm: string) {
        return loaders.loadPdbIhm(this.plugin, pdbIhm);
    }

    loadEmdb(emdb: string, options?: { detail?: number }) {
        return loaders.loadEmdb(this.plugin, emdb, options);
    }

    loadAlphaFoldDb(afdb: string) {
        return loaders.loadAlphaFoldDb(this.plugin, afdb);
    }

    loadModelArchive(id: string) {
        return loaders.loadModelArchive(this.plugin, id);
    }

    /**
     * @example Load X-ray density from volume server
        viewer.loadVolumeFromUrl({
            url: 'https://www.ebi.ac.uk/pdbe/densities/x-ray/1tqn/cell?detail=3',
            format: 'dscif',
            isBinary: true
        }, [{
            type: 'relative',
            value: 1.5,
            color: 0x3362B2
        }, {
            type: 'relative',
            value: 3,
            color: 0x33BB33,
            volumeIndex: 1
        }, {
            type: 'relative',
            value: -3,
            color: 0xBB3333,
            volumeIndex: 1
        }], {
            entryId: ['2FO-FC', 'FO-FC'],
            isLazy: true
        });
     * *********************
     * @example Load EM density from volume server
        viewer.loadVolumeFromUrl({
            url: 'https://maps.rcsb.org/em/emd-30210/cell?detail=6',
            format: 'dscif',
            isBinary: true
        }, [{
            type: 'relative',
            value: 1,
            color: 0x3377aa
        }], {
            entryId: 'EMD-30210',
            isLazy: true
        });
     */
    loadVolumeFromUrl({ url, format, isBinary }: { url: string, format: BuildInVolumeFormat, isBinary: boolean }, isovalues: loaders.VolumeIsovalueInfo[], options?: { entryId?: string | string[], isLazy?: boolean }) {
        return loaders.loadVolumeFromUrl(this.plugin, { url, format, isBinary }, isovalues, options);
    }

    loadFullResolutionEMDBMap(emdbId: string, options: { isoValue: Volume.IsoValue, color?: Color }) {
        return loaders.loadFullResolutionEMDBMap(this.plugin, emdbId, options);
    }

    /**
     * @example
     *  viewer.loadTrajectory({
     *      model: { kind: 'model-url', url: 'villin.gro', format: 'gro' },
     *      coordinates: { kind: 'coordinates-url', url: 'villin.xtc', format: 'xtc', isBinary: true },
     *      preset: 'all-models' // or 'default'
     *  });
     */
    loadTrajectory(params: loaders.LoadTrajectoryParams) {
        return loaders.loadTrajectory(this.plugin, params);
    }

    loadMvsFromUrl(url: string, format: 'mvsj' | 'mvsx', options?: MVSLoadOptions) {
        return loaders.loadMVSFromUrl(this.plugin, url, format, options);
    }

    /** Load MolViewSpec from `data`.
     * If `format` is 'mvsj', `data` must be a string or a Uint8Array containing a UTF8-encoded string.
     * If `format` is 'mvsx', `data` must be a Uint8Array or a string containing base64-encoded binary data prefixed with 'base64,'. */
    loadMvsData(data: string | Uint8Array<ArrayBuffer>, format: 'mvsj' | 'mvsx', options?: MVSLoadOptions) {
        return loaders.loadMvsData(this.plugin, data, format, options);
    }

    loadFiles(files: File[]) {
        return loaders.loadFiles(this.plugin, files);
    }

    loadUrl(url: string, format: string, isBinary = false) {
        return loaders.loadUrl(this.plugin, url, format, isBinary);
    }

    handleResize() {
        this.plugin.layout.events.updated.next(void 0);
    }

    /**
     * Triggers structure element selection or highlighting based on the provided
     * MolScript expression or StructureElement schema. Focus action will only apply to the
     * first structure that matches the criteria.
     *
     * If neither `expression` nor `elements` are provided, all selections/highlights
     * will be cleared based on the specified `action`.
     */
    structureInteractivity(options: StructureInteractivityOptions) {
        return applyStructureInteractivity(this.plugin, options);
    }

    dispose() {
        this._events.dispose();
        this.plugin.dispose();
    }
}

/** Demonstration of usage of MVS builder */
export function builderDemo2() {
    const builder = MVSData.createBuilder();
    builder.canvas({ background_color: 'white' });
    const struct = builder
        .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/5hwa.bcif' })
        .parse({ format: 'bcif' })
        .modelStructure();
    const cartoon = struct.component().representation({ type: 'cartoon', size_factor: 0.5 });
    const carb = struct.component({ selector: 'branched' }).representation({ type: 'carbohydrate' });

    carb.color({ color: 'green/lightgreen', ref: 'color1' });
    // carb.color({ color: 'red/orange', selector: { auth_seq_id: 1 }, ref: 'color2' });

    const anim = builder.animation({ autoplay: true });
    anim
        .interpolate({
            target_ref: 'color1',
            kind: 'color',
            property: 'color',
            start_ms: 1000,
            duration_ms: 1000,
            start: 'green/lightgreen',
            end: 'yellow',
        }).interpolate({
            target_ref: 'color1',
            kind: 'color',
            property: 'color',
            start_ms: 3000,
            duration_ms: 1000,
            end: '#ff0000/#ffffff',
        })
        .interpolate({
            target_ref: 'color1',
            kind: 'color',
            property: 'color',
            start_ms: 5000,
            duration_ms: 4000,
            palette: { kind: 'continuous', colors: ['red/white', 'orange/white', 'yellow/white', 'yellowgreen/white', 'green/black', 'cyan/white', 'blue/black', 'blue/white'] },
        })
        .interpolate({
            target_ref: 'color1',
            kind: 'color',
            property: 'color',
            start_ms: 9000,
            duration_ms: 10_000,
            alternate_direction: true,
            frequency: 10,
            easing: 'sin-in-out',
            palette: { kind: 'continuous', colors: ['blue/white', 'white/blue'] },
        });

    // const annotationUri = `data:text/plain,
    //     data_annot
    //     loop_
    //     _annot.auth_seq_id
    //     _annot.color
    //     _annot.elem
    //     1 green/cyan  S
    //     2 blue/yellow N
    //     3 blue/white  N
    //     4 blue  O
    //     `;

    // // Direct
    // carb.colorFromUri({
    //     uri: annotationUri,
    //     format: 'cif',
    //     schema: 'all_atomic',
    //     category_name: 'annot',
    //     field_name: 'color',
    // });

    // // Categorical with dict
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'categorical',
    //         colors: { 1: 'red', 2: 'orange/cyan' },
    //         missing_color: 'green/lightgreen',
    //     },
    // });

    // // Categorical with named dict
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'label_comp_id',
    //     palette: {
    //         kind: 'categorical',
    //         colors: 'CarbohydrateSymbol',
    //         missing_color: 'lightgreen/green',
    //     },
    // });

    // // Categorical with list
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'categorical',
    //         colors: ['red', 'pink/cyan'],
    //         missing_color: 'green/lightgreen',
    //     },
    // });

    // // Categorical with named list
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'categorical',
    //         colors: 'Blues',
    //         missing_color: 'green/lightgreen',
    //     },
    // });

    // // Discrete
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'discrete',
    //         colors: ['red', 'pink/cyan'],
    //         value_domain: [1, 4],
    //         // colors: 'Blues',
    //         // reverse: true,
    //     },
    // });

    // // Discrete with checkpoints
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'discrete',
    //         colors: [['pink/cyan', null, null], [null, 2, 3], ['orange', 3, 3]],
    //         mode: 'absolute',
    //     },
    // });

    // // Continuous
    // carb.colorFromSource({
    //     schema: 'all_atomic',
    //     category_name: 'atom_site',
    //     field_name: 'auth_seq_id',
    //     palette: {
    //         kind: 'continuous',
    //         // colors: ['red/darkblue', 'yellow/lightgreen',  'white'],
    //         // value_domain: [1, 4],
    //         // colors: 'Blues',
    //         colors: [['red/darkblue', 2], ['yellow/lightgreen', 3]], mode: 'absolute',
    //         underflow_color: 'orange/#ffffdd',
    //         overflow_color: 'cyan/#ddffff',
    //         // reverse: true,
    //     },
    // });

    // carb.color({ color: 'blue/green' });
    const snapshot = builder.getSnapshot({ linger_duration_ms: 3000 });
    return MVSData.createMultistate([snapshot, snapshot]);
}

// TODO: test animations
