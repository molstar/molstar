/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { Representation } from 'mol-repr/representation';
import { StructureRepresentationRegistry } from 'mol-repr/structure/registry';
import { State, Transform, Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { LogEntry } from 'mol-util/log-entry';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { merge } from 'rxjs';
import { BuiltInPluginBehaviors } from './behavior';
import { PluginCommand, PluginCommands } from './command';
import { PluginSpec } from './spec';
import { PluginState } from './state';
import { TaskManager } from './util/task-manager';
import { Color } from 'mol-util/color';
import { LociLabelEntry, LociLabelManager } from './util/loci-label-manager';
import { ajaxGet } from 'mol-util/data-source';
import { CustomPropertyRegistry } from './util/custom-prop-registry';
import { VolumeRepresentationRegistry } from 'mol-repr/volume/registry';
import { PLUGIN_VERSION, PLUGIN_VERSION_DATE } from './version';
import { PluginLayout } from './layout';
import { List } from 'immutable';
import { StateTransformParameters } from './ui/state/common';
import { DataFormatRegistry } from './state/actions/basic';

export class PluginContext {
    private disposed = false;
    private ev = RxEventHelper.create();
    private tasks = new TaskManager();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommand.Manager();

    readonly events = {
        state: {
            cell: {
                stateUpdated: merge(this.state.dataState.events.cell.stateUpdated, this.state.behaviorState.events.cell.stateUpdated),
                created: merge(this.state.dataState.events.cell.created, this.state.behaviorState.events.cell.created),
                removed: merge(this.state.dataState.events.cell.removed, this.state.behaviorState.events.cell.removed),
            },
            object: {
                created: merge(this.state.dataState.events.object.created, this.state.behaviorState.events.object.created),
                removed: merge(this.state.dataState.events.object.removed, this.state.behaviorState.events.object.removed),
                updated: merge(this.state.dataState.events.object.updated, this.state.behaviorState.events.object.updated)
            },
            cameraSnapshots: this.state.cameraSnapshots.events,
            snapshots: this.state.snapshots.events,
        },
        log: this.ev<LogEntry>(),
        task: this.tasks.events,
        labels: {
            highlight: this.ev<{ entries: ReadonlyArray<LociLabelEntry> }>()
        },
        canvad3d: {
            settingsUpdated: this.ev()
        }
    };

    readonly behaviors = {
        canvas: {
            highlightLoci: this.ev.behavior<{ loci: Loci, repr?: Representation.Any }>({ loci: EmptyLoci }),
            selectLoci: this.ev.behavior<{ loci: Loci, repr?: Representation.Any }>({ loci: EmptyLoci }),
        },
        command: this.commands.behaviour
    };

    readonly canvas3d: Canvas3D;
    readonly layout: PluginLayout = new PluginLayout(this);

    readonly lociLabels: LociLabelManager;

    readonly structureRepresentation = {
        registry: new StructureRepresentationRegistry(),
        themeCtx: { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() } as ThemeRegistryContext
    }

    readonly volumeRepresentation = {
        registry: new VolumeRepresentationRegistry(),
        themeCtx: { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() } as ThemeRegistryContext
    }

    readonly dataFormat = {
        registry: new DataFormatRegistry()
    }

    readonly customModelProperties = new CustomPropertyRegistry();
    readonly customParamEditors = new Map<string, StateTransformParameters.Class>();

    initViewer(canvas: HTMLCanvasElement, container: HTMLDivElement) {
        try {
            this.layout.setRoot(container);
            if (this.spec.initialLayout) this.layout.updateState(this.spec.initialLayout);
            (this.canvas3d as Canvas3D) = Canvas3D.create(canvas, container);
            PluginCommands.Canvas3D.SetSettings.dispatch(this, { settings: { backgroundColor: Color(0xFCFBF9) } });
            this.canvas3d.animate();
            return true;
        } catch (e) {
            this.log.error('' + e);
            console.error(e);
            return false;
        }
    }

    readonly log = {
        entries: List<LogEntry>(),
        entry: (e: LogEntry) => this.events.log.next(e),
        error: (msg: string) => this.events.log.next(LogEntry.error(msg)),
        message: (msg: string) => this.events.log.next(LogEntry.message(msg)),
        info: (msg: string) => this.events.log.next(LogEntry.info(msg)),
        warn: (msg: string) => this.events.log.next(LogEntry.warning(msg)),
    };

    /**
     * This should be used in all transform related request so that it could be "spoofed" to allow
     * "static" access to resources.
     */
    fetch(url: string, type: 'string' | 'binary' = 'string', body?: string): Task<string | Uint8Array> {
        return ajaxGet({ url, type, body });
        // const req = await fetch(url, { referrerPolicy: 'origin-when-cross-origin' });
        // return type === 'string' ? await req.text() : new Uint8Array(await req.arrayBuffer());
    }

    runTask<T>(task: Task<T>) {
        return this.tasks.run(task);
    }

    dispose() {
        if (this.disposed) return;
        this.commands.dispose();
        this.canvas3d.dispose();
        this.ev.dispose();
        this.state.dispose();
        this.tasks.dispose();
        this.disposed = true;
    }

    applyTransform(state: State, a: Transform.Ref, transformer: Transformer, params: any) {
        const tree = state.tree.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    updateTransform(state: State, a: Transform.Ref, params: any) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    private initBuiltInBehavior() {
        BuiltInPluginBehaviors.State.registerDefault(this);
        BuiltInPluginBehaviors.Representation.registerDefault(this);
        BuiltInPluginBehaviors.Camera.registerDefault(this);
        BuiltInPluginBehaviors.Misc.registerDefault(this);

        merge(this.state.dataState.events.log, this.state.behaviorState.events.log).subscribe(e => this.events.log.next(e));
    }

    private async initBehaviors() {
        const tree = this.state.behaviorState.tree.build();

        for (const b of this.spec.behaviors) {
            tree.toRoot().apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
        }

        await this.runTask(this.state.behaviorState.updateTree(tree, true));
    }

    private initDataActions() {
        for (const a of this.spec.actions) {
            this.state.dataState.actions.add(a.action);
        }
    }

    private initCustomParamEditors() {
        if (!this.spec.customParamEditors) return;

        for (const [t, e] of this.spec.customParamEditors) {
            this.customParamEditors.set(t.id, e);
        }
    }

    constructor(public spec: PluginSpec) {
        this.events.log.subscribe(e => this.log.entries = this.log.entries.push(e));

        this.initBuiltInBehavior();

        this.initBehaviors();
        this.initDataActions();
        this.initCustomParamEditors();

        this.lociLabels = new LociLabelManager(this);

        this.log.message(`Mol* Plugin ${PLUGIN_VERSION} [${PLUGIN_VERSION_DATE.toLocaleString()}]`);
    }
}