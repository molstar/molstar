/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { setAutoFreeze } from 'immer';
import { List } from 'immutable';
import { merge } from 'rxjs';
import { Canvas3D, DefaultCanvas3DParams } from '../mol-canvas3d/canvas3d';
import { CustomProperty } from '../mol-model-props/common/custom-property';
import { Model, Structure } from '../mol-model/structure';
import { DataBuilder } from '../mol-plugin-state/builder/data';
import { StructureBuilder } from '../mol-plugin-state/builder/structure';
import { DataFormatRegistry } from '../mol-plugin-state/formats/registry';
import { StructureSelectionQueryRegistry } from '../mol-plugin-state/helpers/structure-selection-query';
import { CameraManager } from '../mol-plugin-state/manager/camera';
import { InteractivityManager } from '../mol-plugin-state/manager/interactivity';
import { LociLabel, LociLabelManager } from '../mol-plugin-state/manager/loci-label';
import { StructureComponentManager } from '../mol-plugin-state/manager/structure/component';
import { StructureFocusManager } from '../mol-plugin-state/manager/structure/focus';
import { StructureHierarchyManager } from '../mol-plugin-state/manager/structure/hierarchy';
import { StructureHierarchyRef } from '../mol-plugin-state/manager/structure/hierarchy-state';
import { StructureMeasurementManager } from '../mol-plugin-state/manager/structure/measurement';
import { StructureSelectionManager } from '../mol-plugin-state/manager/structure/selection';
import { PluginUIComponent } from '../mol-plugin-ui/base';
import { StateTransformParameters } from '../mol-plugin-ui/state/common';
import { Representation } from '../mol-repr/representation';
import { StructureRepresentationRegistry } from '../mol-repr/structure/registry';
import { VolumeRepresentationRegistry } from '../mol-repr/volume/registry';
import { StateTransform } from '../mol-state';
import { Progress, Task, RuntimeContext } from '../mol-task';
import { ColorTheme } from '../mol-theme/color';
import { SizeTheme } from '../mol-theme/size';
import { ThemeRegistryContext } from '../mol-theme/theme';
import { Color } from '../mol-util/color';
import { ajaxGet } from '../mol-util/data-source';
import { isDebugMode, isProductionMode } from '../mol-util/debug';
import { ModifiersKeys } from '../mol-util/input/input-observer';
import { LogEntry } from '../mol-util/log-entry';
import { RxEventHelper } from '../mol-util/rx-event-helper';
import { BuiltInPluginBehaviors } from './behavior';
import { PluginBehavior } from './behavior/behavior';
import { PluginCommandManager } from './command';
import { PluginCommands } from './commands';
import { PluginConfig, PluginConfigManager } from './config';
import { LeftPanelTabName, PluginLayout } from './layout';
import { PluginSpec } from './spec';
import { PluginState } from './state';
import { SubstructureParentHelper } from './util/substructure-parent-helper';
import { TaskManager } from './util/task-manager';
import { PluginToastManager } from './util/toast';
import { ViewportScreenshotHelper } from './util/viewport-screenshot';
import { PLUGIN_VERSION, PLUGIN_VERSION_DATE } from './version';
import { AssetManager } from '../mol-util/assets';
import { PluginStateSnapshotManager } from '../mol-plugin-state/manager/snapshots';
import { PluginAnimationManager } from '../mol-plugin-state/manager/animation';
import { objectForEach } from '../mol-util/object';
import { VolumeHierarchyManager } from '../mol-plugin-state/manager/volume/hierarchy';
import { filter, take } from 'rxjs/operators';

export class PluginContext {
    runTask = <T>(task: Task<T>) => this.tasks.run(task);

    private disposed = false;
    private ev = RxEventHelper.create();
    private tasks = new TaskManager();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommandManager();

    readonly events = {
        log: this.ev<LogEntry>(),
        task: this.tasks.events,
        canvas3d: {
            settingsUpdated: this.ev(),
        }
    } as const

    readonly config = new PluginConfigManager(this.spec.config);

    private canvas3dInit = this.ev.behavior<boolean>(false);
    readonly behaviors = {
        state: {
            isAnimating: this.ev.behavior<boolean>(false),
            isUpdating: this.ev.behavior<boolean>(false),
            // TODO: should there be separate "updated" event?
            //   Often, this is used to indicate that the state has updated
            //   and it might not be the best way to react to state updates.
            isBusy: this.ev.behavior<boolean>(false)
        },
        interaction: {
            hover: this.ev.behavior<InteractivityManager.HoverEvent>({ current: Representation.Loci.Empty, modifiers: ModifiersKeys.None, buttons: 0, button: 0 }),
            click: this.ev.behavior<InteractivityManager.ClickEvent>({ current: Representation.Loci.Empty, modifiers: ModifiersKeys.None, buttons: 0, button: 0 }),
            selectionMode: this.ev.behavior<boolean>(false)
        },
        labels: {
            highlight: this.ev.behavior<{ labels: ReadonlyArray<LociLabel> }>({ labels: [] })
        },
        layout: {
            leftPanelTabName: this.ev.behavior<LeftPanelTabName>('root')
        },
        canvas3d: {
            initialized: this.canvas3dInit.pipe(filter(v => !!v), take(1))
        }
    } as const;

    readonly canvas3d: Canvas3D | undefined;
    readonly layout = new PluginLayout(this);

    readonly representation = {
        structure: {
            registry: new StructureRepresentationRegistry(),
            themes: { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() } as ThemeRegistryContext,
        },
        volume: {
            registry: new VolumeRepresentationRegistry(),
            themes: { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() } as ThemeRegistryContext
        }
    } as const;

    readonly query = {
        structure: {
            registry: new StructureSelectionQueryRegistry()
        }
    } as const;

    readonly dataFormats = new DataFormatRegistry();

    readonly builders = {
        data: new DataBuilder(this),
        structure: void 0 as any as StructureBuilder
    };

    build() {
        return this.state.data.build();
    }

    readonly helpers = {
        substructureParent: new SubstructureParentHelper(this),
        viewportScreenshot: void 0 as ViewportScreenshotHelper | undefined
    } as const;

    readonly managers = {
        structure: {
            hierarchy: new StructureHierarchyManager(this),
            component: new StructureComponentManager(this),
            measurement: new StructureMeasurementManager(this),
            selection: new StructureSelectionManager(this),
            focus: new StructureFocusManager(this),
        },
        volume: {
            hierarchy: new VolumeHierarchyManager(this)
        },
        interactivity: void 0 as any as InteractivityManager,
        camera: new CameraManager(this),
        animation: new PluginAnimationManager(this),
        snapshot: new PluginStateSnapshotManager(this),
        lociLabels: void 0 as any as LociLabelManager,
        toast: new PluginToastManager(this),
        asset: new AssetManager()
    } as const

    readonly customModelProperties = new CustomProperty.Registry<Model>();
    readonly customStructureProperties = new CustomProperty.Registry<Structure>();
    readonly customParamEditors = new Map<string, StateTransformParameters.Class>();

    readonly customStructureControls = new Map<string, { new(): PluginUIComponent<any, any, any> }>();
    readonly genericRepresentationControls = new Map<string, (selection: StructureHierarchyManager['selection']) => [StructureHierarchyRef[], string]>();

    /**
     * Used to store application specific custom state which is then available
     * to State Actions and similar constructs via the PluginContext.
     */
    readonly customState: unknown = Object.create(null);

    initViewer(canvas: HTMLCanvasElement, container: HTMLDivElement) {
        try {
            this.layout.setRoot(container);
            if (this.spec.layout && this.spec.layout.initial) this.layout.setProps(this.spec.layout.initial);

            (this.canvas3d as Canvas3D) = Canvas3D.fromCanvas(canvas);
            this.canvas3dInit.next(true);
            const renderer = this.canvas3d!.props.renderer;
            PluginCommands.Canvas3D.SetSettings(this, { settings: { renderer: { ...renderer, backgroundColor: Color(0xFCFBF9) } } });
            this.canvas3d!.animate();
            (this.helpers.viewportScreenshot as ViewportScreenshotHelper) = new ViewportScreenshotHelper(this);
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

    /** return true is animating or updating */
    get isBusy() {
        return this.behaviors.state.isAnimating.value || this.behaviors.state.isUpdating.value;
    }

    get selectionMode() {
        return this.behaviors.interaction.selectionMode.value;
    }

    set selectionMode(mode: boolean) {
        this.behaviors.interaction.selectionMode.next(mode);
    }

    dataTransaction(f: (ctx: RuntimeContext) => Promise<void> | void, options?: { canUndo?: string | boolean }) {
        return this.runTask(this.state.data.transaction(f, options));
    }

    requestTaskAbort(progress: Progress, reason?: string) {
        this.tasks.requestAbort(progress, reason);
    }

    clear(resetViewportSettings = false) {
        if (resetViewportSettings) this.canvas3d?.setProps(DefaultCanvas3DParams);
        return PluginCommands.State.RemoveObject(this, { state: this.state.data, ref: StateTransform.RootRef });
    }

    dispose() {
        if (this.disposed) return;
        this.commands.dispose();
        this.canvas3d?.dispose();
        this.ev.dispose();
        this.state.dispose();
        this.tasks.dispose();
        this.layout.dispose();
        this.helpers.substructureParent.dispose();

        objectForEach(this.managers, m => (m as any)?.dispose?.());
        objectForEach(this.managers.structure, m => (m as any)?.dispose?.());

        this.disposed = true;
    }

    private initBehaviorEvents() {
        merge(this.state.data.behaviors.isUpdating, this.state.behaviors.behaviors.isUpdating).subscribe(u => {
            if (this.behaviors.state.isUpdating.value !== u) this.behaviors.state.isUpdating.next(u);
        });

        const timeoutMs = this.config.get(PluginConfig.General.IsBusyTimeoutMs) || 750;
        const isBusy = this.behaviors.state.isBusy;

        let timeout: any = void 0;
        const setBusy = () => {
            if (!isBusy.value) isBusy.next(true);
        };
        const reset = () => {
            if (timeout !== void 0) clearTimeout(timeout);
            timeout = void 0;
        };

        merge(this.behaviors.state.isUpdating, this.behaviors.state.isAnimating).subscribe(v => {
            const isUpdating = this.behaviors.state.isUpdating.value;
            const isAnimating = this.behaviors.state.isAnimating.value;

            if (isUpdating || isAnimating) {
                if (!isBusy.value) {
                    reset();
                    timeout = setTimeout(setBusy, timeoutMs);
                }
            } else {
                reset();
                isBusy.next(false);
            }
        });

        this.behaviors.interaction.selectionMode.subscribe(v => {
            if (!v) {
                this.managers.interactivity?.lociSelects.deselectAll();
            }
        });
    }

    private initBuiltInBehavior() {
        BuiltInPluginBehaviors.State.registerDefault(this);
        BuiltInPluginBehaviors.Representation.registerDefault(this);
        BuiltInPluginBehaviors.Camera.registerDefault(this);
        BuiltInPluginBehaviors.Misc.registerDefault(this);

        merge(this.state.data.events.log, this.state.behaviors.events.log).subscribe(e => this.events.log.next(e));
    }

    private async initBehaviors() {
        let tree = this.state.behaviors.build();

        for (const cat of Object.keys(PluginBehavior.Categories)) {
            tree.toRoot().apply(PluginBehavior.CreateCategory, { label: (PluginBehavior.Categories as any)[cat] }, { ref: cat, state: { isLocked: true } });
        }

        // Init custom properties 1st
        for (const b of this.spec.behaviors) {
            const cat = PluginBehavior.getCategoryId(b.transformer);
            if (cat !== 'custom-props') continue;

            tree.to(PluginBehavior.getCategoryId(b.transformer)).apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
        }
        await this.runTask(this.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true, doNotLogTiming: true }));

        tree = this.state.behaviors.build();
        for (const b of this.spec.behaviors) {
            const cat = PluginBehavior.getCategoryId(b.transformer);
            if (cat === 'custom-props') continue;

            tree.to(PluginBehavior.getCategoryId(b.transformer)).apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
        }
        await this.runTask(this.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true, doNotLogTiming: true }));
    }

    private initDataActions() {
        for (const a of this.spec.actions) {
            this.state.data.actions.add(a.action);
        }
    }

    private initAnimations() {
        if (!this.spec.animations) return;
        for (const anim of this.spec.animations) {
            this.managers.animation.register(anim);
        }
    }

    private initCustomParamEditors() {
        if (!this.spec.customParamEditors) return;

        for (const [t, e] of this.spec.customParamEditors) {
            this.customParamEditors.set(t.id, e);
        }
    }

    constructor(public spec: PluginSpec) {
        // the reason for this is that sometimes, transform params get modified inline (i.e. palette.valueLabel)
        // and freezing the params object causes "read-only exception"
        // TODO: is this the best place to do it?
        setAutoFreeze(false);

        this.events.log.subscribe(e => this.log.entries = this.log.entries.push(e));

        this.initBehaviorEvents();
        this.initBuiltInBehavior();

        this.initBehaviors();
        this.initDataActions();
        this.initAnimations();
        this.initCustomParamEditors();

        (this.managers.interactivity as InteractivityManager) = new InteractivityManager(this);
        (this.managers.lociLabels as LociLabelManager) = new LociLabelManager(this);

        (this.builders.structure as StructureBuilder) = new StructureBuilder(this);

        this.log.message(`Mol* Plugin ${PLUGIN_VERSION} [${PLUGIN_VERSION_DATE.toLocaleString()}]`);
        if (!isProductionMode) this.log.message(`Development mode enabled`);
        if (isDebugMode) this.log.message(`Debug mode enabled`);
    }
}