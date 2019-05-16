/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { List } from 'immutable';
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { CustomPropertyRegistry } from 'mol-model-props/common/custom-property-registry';
import { StructureRepresentationRegistry } from 'mol-repr/structure/registry';
import { VolumeRepresentationRegistry } from 'mol-repr/volume/registry';
import { State, StateTransform, StateTransformer } from 'mol-state';
import { Task } from 'mol-task';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Color } from 'mol-util/color';
import { ajaxGet } from 'mol-util/data-source';
import { LogEntry } from 'mol-util/log-entry';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { merge } from 'rxjs';
import { BuiltInPluginBehaviors } from './behavior';
import { PluginBehavior } from './behavior/behavior';
import { PluginCommand, PluginCommands } from './command';
import { PluginLayout } from './layout';
import { PluginSpec } from './spec';
import { PluginState } from './state';
import { DataFormatRegistry } from './state/actions/data-format';
import { StateTransformParameters } from './ui/state/common';
import { LociLabelEntry, LociLabelManager } from './util/loci-label-manager';
import { TaskManager } from './util/task-manager';
import { PLUGIN_VERSION, PLUGIN_VERSION_DATE } from './version';
import { StructureElementSelectionManager } from './util/structure-element-selection';
import { SubstructureParentHelper } from './util/substructure-parent-helper';
import { Representation } from 'mol-repr/representation';
import { ModifiersKeys } from 'mol-util/input/input-observer';
import { isProductionMode, isDebugMode } from 'mol-util/debug';
import { Model, Structure } from 'mol-model/structure';

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
        canvas3d: {
            settingsUpdated: this.ev()
        }
    };

    readonly behaviors = {
        state: {
            isAnimating: this.ev.behavior<boolean>(false),
            isUpdating: this.ev.behavior<boolean>(false)
        },
        canvas3d: {
            highlight: this.ev.behavior<Canvas3D.HighlightEvent>({ current: Representation.Loci.Empty, prev: Representation.Loci.Empty }),
            click: this.ev.behavior<Canvas3D.ClickEvent>({ current: Representation.Loci.Empty, modifiers: ModifiersKeys.None, buttons: 0 })
        },
        labels: {
            highlight: this.ev.behavior<{ entries: ReadonlyArray<LociLabelEntry> }>({ entries: [] })
        }
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

    readonly customModelProperties = new CustomPropertyRegistry<Model>();
    readonly customStructureProperties = new CustomPropertyRegistry<Structure>();
    readonly customParamEditors = new Map<string, StateTransformParameters.Class>();

    readonly helpers = {
        structureSelection: new StructureElementSelectionManager(this),
        substructureParent: new SubstructureParentHelper(this)
    };

    initViewer(canvas: HTMLCanvasElement, container: HTMLDivElement) {
        try {
            this.layout.setRoot(container);
            if (this.spec.layout && this.spec.layout.initial) this.layout.setProps(this.spec.layout.initial);
            (this.canvas3d as Canvas3D) = Canvas3D.fromCanvas(canvas);
            const renderer = this.canvas3d.props.renderer;
            PluginCommands.Canvas3D.SetSettings.dispatch(this, { settings: { renderer: { ...renderer, backgroundColor: Color(0xFCFBF9) } } });
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
    readonly fetch = ajaxGet

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
        this.layout.dispose();
        this.disposed = true;
    }

    applyTransform(state: State, a: StateTransform.Ref, transformer: StateTransformer, params: any) {
        const tree = state.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    updateTransform(state: State, a: StateTransform.Ref, params: any) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    private initBehaviorEvents() {
        merge(this.state.dataState.events.isUpdating, this.state.behaviorState.events.isUpdating).subscribe(u => {
            this.behaviors.state.isUpdating.next(u);
        });
    }

    private initBuiltInBehavior() {
        BuiltInPluginBehaviors.State.registerDefault(this);
        BuiltInPluginBehaviors.Representation.registerDefault(this);
        BuiltInPluginBehaviors.Camera.registerDefault(this);
        BuiltInPluginBehaviors.Misc.registerDefault(this);

        merge(this.state.dataState.events.log, this.state.behaviorState.events.log).subscribe(e => this.events.log.next(e));
    }

    private async initBehaviors() {
        const tree = this.state.behaviorState.build();

        for (const cat of Object.keys(PluginBehavior.Categories)) {
            tree.toRoot().apply(PluginBehavior.CreateCategory, { label: (PluginBehavior.Categories as any)[cat] }, { ref: cat, state: { isLocked: true } });
        }

        for (const b of this.spec.behaviors) {
            tree.to(PluginBehavior.getCategoryId(b.transformer)).apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
        }

        await this.runTask(this.state.behaviorState.updateTree(tree, { doNotUpdateCurrent: true, doNotLogTiming: true }));
    }

    private initDataActions() {
        for (const a of this.spec.actions) {
            this.state.dataState.actions.add(a.action);
        }
    }

    private initAnimations() {
        if (!this.spec.animations) return;
        for (const anim of this.spec.animations) {
            this.state.animation.register(anim);
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

        this.initBehaviorEvents();
        this.initBuiltInBehavior();

        this.initBehaviors();
        this.initDataActions();
        this.initAnimations();
        this.initCustomParamEditors();

        this.lociLabels = new LociLabelManager(this);

        this.log.message(`Mol* Plugin ${PLUGIN_VERSION} [${PLUGIN_VERSION_DATE.toLocaleString()}]`);
        if (!isProductionMode) this.log.message(`Development mode enabled`);
        if (isDebugMode) this.log.message(`Debug mode enabled`);
    }
}