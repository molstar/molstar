/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { produce } from '@molstar/core/util/produce';
import { List } from 'immutable';
import { merge, Subscription } from 'rxjs';
import { debounceTime, filter, take, throttleTime } from 'rxjs/operators';
import { Canvas3D, Canvas3DContext, DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { resizeCanvas } from '@molstar/graphics/canvas3d/util';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { Model, Structure } from '@molstar/model/model/structure';
import { DataBuilder } from '@molstar/plugin/state/builder/data';
import { StructureBuilder } from '@molstar/plugin/state/builder/structure';
import { DataFormatRegistry } from '@molstar/plugin/state/formats/registry';
import { StructureSelectionQueryRegistry } from '@molstar/plugin/state/queries/structure/registry';
import { PluginAnimationManager } from '@molstar/plugin/state/manager/animation';
import { CameraManager } from '@molstar/plugin/state/manager/camera';
import { InteractivityManager } from '@molstar/plugin/state/manager/interactivity';
import { type LociLabel, LociLabelManager } from '@molstar/plugin/state/manager/loci-label';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { StructureComponentManager } from '@molstar/plugin/state/manager/structure/component';
import { StructureFocusManager } from '@molstar/plugin/state/manager/structure/focus';
import { StructureHierarchyManager } from '@molstar/plugin/state/manager/structure/hierarchy';
import type { StructureHierarchyRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { StructureMeasurementManager } from '@molstar/plugin/state/manager/structure/measurement';
import { StructureSelectionManager } from '@molstar/plugin/state/manager/structure/selection';
import { ParticleHierarchyManager } from '@molstar/plugin/state/manager/particles/hierarchy';
import { VolumeHierarchyManager } from '@molstar/plugin/state/manager/volume/hierarchy';
import { MarkdownExtensionManager } from '@molstar/plugin/state/manager/markdown-extensions';
import { type LeftPanelTabName, PluginLayout } from '@molstar/plugin/layout';
import { Representation } from '@molstar/graphics/repr/representation';
import { ParticleRepresentationRegistry } from '@molstar/graphics/repr/particles/registry';
import { StructureRepresentationRegistry } from '@molstar/graphics/repr/structure/registry';
import { VolumeRepresentationRegistry } from '@molstar/graphics/repr/volume/registry';
import { StateTransform } from '@molstar/core/state';
import { RuntimeContext, Scheduler, Task } from '@molstar/core/task';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { AssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { ajaxGet } from '@molstar/core/util/data-source';
import { isDebugMode, isProductionMode } from '@molstar/core/util/debug';
import { EmptyKeyInput, type KeyInput, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { LogEntry } from '@molstar/core/util/log-entry';
import { objectForEach } from '@molstar/core/util/object';
import { RxEventHelper } from '@molstar/core/util/rx-event-helper';
import { PluginAnimationLoop } from '@molstar/plugin/animation-loop';
import { BuiltInPluginBehaviors } from '@molstar/plugin/behavior/built-in';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginCommandManager } from '@molstar/plugin/command';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginConfig, PluginConfigManager } from '@molstar/plugin/config';
import type { PluginRegistryEntry, PluginSpec } from '@molstar/plugin/spec';
import { registerEntries } from '@molstar/plugin/registry-entry';
import { PluginState } from '@molstar/plugin/state';
import { SubstructureParentHelper } from '@molstar/plugin/util/substructure-parent-helper';
import { TaskManager } from '@molstar/plugin/util/task-manager';
import { PluginToastManager } from '@molstar/plugin/util/toast';
import { ViewportScreenshotHelper } from '@molstar/plugin/util/viewport-screenshot';
import { PLUGIN_VERSION } from '@molstar/plugin/version';
import { setSaccharideCompIdMapType } from '@molstar/model/model/structure/structure/carbohydrates/constants';
import { DragAndDropManager } from '@molstar/plugin/state/manager/drag-and-drop';
import { ErrorContext } from '@molstar/core/util/error-context';
import { PluginContainer } from '@molstar/plugin/container';
import type { Volume } from '@molstar/model/model/volume';

export type PluginInitializedState = { kind: 'no' } | { kind: 'yes' } | { kind: 'error'; error: any };

/** `PluginSpec` keys removed in 6.0; the constructor rejects them with a pointer to registry entries. */
const REMOVED_SPEC_KEYS = ['actions', 'animations', 'customFormats'] as const;

export class PluginContext {
  runTask = <T>(task: Task<T>, params?: { useOverlay?: boolean }) => this.managers.task.run(task, params);
  resolveTask = <T>(object: Task<T> | T | undefined) => {
    if (!object) return void 0;
    if (Task.is(object)) return this.runTask(object);
    return object;
  };

  protected subs: Subscription[] = [];
  private initCanvas3dPromiseCallbacks: [res: () => void, rej: (err: any) => void] = [() => {}, () => {}];
  private _isInitialized = false;
  private initializedPromiseCallbacks: [res: () => void, rej: (err: any) => void] = [() => {}, () => {}];

  private disposed = false;
  private container: PluginContainer | undefined = void 0;
  private ev = RxEventHelper.create();

  readonly config = new PluginConfigManager(this.spec.config); // needed to init state
  readonly state = new PluginState(this);
  readonly commands = new PluginCommandManager();

  private canvas3dInit = this.ev.behavior<boolean>(false);
  readonly behaviors = {
    state: {
      isAnimating: this.ev.behavior<boolean>(false),
      isUpdating: this.ev.behavior<boolean>(false),
      // TODO: should there be separate "updated" event?
      //   Often, this is used to indicate that the state has updated
      //   and it might not be the best way to react to state updates.
      isBusy: this.ev.behavior<boolean>(false),
    },
    interaction: {
      hover: this.ev.behavior<InteractivityManager.HoverEvent>({
        current: Representation.Loci.Empty,
        modifiers: ModifiersKeys.None,
        buttons: 0,
        button: 0,
      }),
      click: this.ev.behavior<InteractivityManager.ClickEvent>({
        current: Representation.Loci.Empty,
        modifiers: ModifiersKeys.None,
        buttons: 0,
        button: 0,
      }),
      drag: this.ev.behavior<InteractivityManager.DragEvent>({
        current: Representation.Loci.Empty,
        modifiers: ModifiersKeys.None,
        buttons: 0,
        button: 0,
        pageStart: Vec2(),
        pageEnd: Vec2(),
      }),
      key: this.ev.behavior<KeyInput>(EmptyKeyInput),
      keyReleased: this.ev.behavior<KeyInput>(EmptyKeyInput),
      selectionMode: this.ev.behavior<boolean>(false),
    },
    labels: {
      highlight: this.ev.behavior<{ labels: ReadonlyArray<LociLabel> }>({ labels: [] }),
    },
    layout: {
      leftPanelTabName: this.ev.behavior<LeftPanelTabName>('root'),
    },
    canvas3d: {
      // TODO: remove in 4.0?
      initialized: this.canvas3dInit.pipe(
        filter((v) => !!v),
        take(1),
      ),
    },
  } as const;

  readonly canvas3dInitialized = new Promise<void>((res, rej) => {
    this.initCanvas3dPromiseCallbacks = [res, rej];
  });

  readonly initialized = new Promise<void>((res, rej) => {
    this.initializedPromiseCallbacks = [res, rej];
  });

  get isInitialized() {
    return this._isInitialized;
  }

  readonly canvas3dContext: Canvas3DContext | undefined;
  readonly canvas3d: Canvas3D | undefined;
  readonly layout = new PluginLayout(this);
  readonly animationLoop = new PluginAnimationLoop(this);

  readonly representation = {
    structure: {
      registry: new StructureRepresentationRegistry(),
      themes: {
        colorThemeRegistry: ColorTheme.createRegistry(),
        sizeThemeRegistry: SizeTheme.createRegistry(),
      } as ThemeRegistryContext,
    },
    volume: {
      registry: new VolumeRepresentationRegistry(),
      themes: {
        colorThemeRegistry: ColorTheme.createRegistry(),
        sizeThemeRegistry: SizeTheme.createRegistry(),
      } as ThemeRegistryContext,
    },
    particles: {
      registry: new ParticleRepresentationRegistry(),
      themes: {
        colorThemeRegistry: ColorTheme.createRegistry(),
        sizeThemeRegistry: SizeTheme.createRegistry(),
      } as ThemeRegistryContext,
    },
  } as const;

  readonly query = {
    structure: {
      registry: new StructureSelectionQueryRegistry(),
    },
  } as const;

  readonly dataFormats = new DataFormatRegistry();

  readonly builders = {
    data: new DataBuilder(this),
    structure: void 0 as any as StructureBuilder,
  };

  build() {
    return this.state.data.build();
  }

  readonly helpers = {
    substructureParent: new SubstructureParentHelper(this),
    viewportScreenshot: void 0 as ViewportScreenshotHelper | undefined,
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
      hierarchy: new VolumeHierarchyManager(this),
    },
    particles: {
      hierarchy: new ParticleHierarchyManager(this),
    },
    interactivity: void 0 as any as InteractivityManager,
    camera: new CameraManager(this),
    animation: new PluginAnimationManager(this),
    snapshot: new PluginStateSnapshotManager(this),
    lociLabels: void 0 as any as LociLabelManager,
    toast: new PluginToastManager(this),
    asset: new AssetManager(),
    task: new TaskManager(),
    markdownExtensions: new MarkdownExtensionManager(this),
    dragAndDrop: new DragAndDropManager(this),
  } as const;

  readonly events = {
    log: this.ev<LogEntry>(),
    task: this.managers.task.events,
    canvas3d: {
      settingsUpdated: this.ev(),
    },
  } as const;

  readonly customModelProperties = new CustomProperty.Registry<Model>();
  readonly customStructureProperties = new CustomProperty.Registry<Structure>();
  readonly customVolumeProperties = new CustomProperty.Registry<Volume>();

  readonly customStructureControls = new Map<
    string,
    { new (): any /* constructible react components with <action.customControl /> */ }
  >();
  readonly customImportControls = new Map<
    string,
    { new (): any /* constructible react components with <action.customControl /> */ }
  >();
  readonly genericRepresentationControls = new Map<
    string,
    (selection: StructureHierarchyManager['selection']) => [StructureHierarchyRef[], string]
  >();

  /**
   * A helper for collecting and notifying errors
   * in async contexts such as custom properties.
   *
   * Individual extensions are responsible for using this
   * context and displaying the errors in appropriate ways.
   */
  readonly errorContext = new ErrorContext();

  /**
   * Used to store application specific custom state which is then available
   * to State Actions and similar constructs via the PluginContext.
   */
  readonly customState: unknown = Object.create(null);

  async initViewerAsync(canvas: HTMLCanvasElement, container: HTMLDivElement, canvas3dContext?: Canvas3DContext) {
    return this._initViewer(canvas, container, canvas3dContext);
  }

  async initContainerAsync(options?: { canvas3dContext?: Canvas3DContext; checkeredCanvasBackground?: boolean }) {
    return this._initContainer(options);
  }

  async mountAsync(
    target: HTMLElement,
    initOptions?: { canvas3dContext?: Canvas3DContext; checkeredCanvasBackground?: boolean },
  ) {
    return this._mount(target, initOptions);
  }

  private _initContainer(options?: { canvas3dContext?: Canvas3DContext; checkeredCanvasBackground?: boolean }) {
    if (this.container) return true;
    const container = new PluginContainer({
      checkeredCanvasBackground: options?.checkeredCanvasBackground,
      canvas: options?.canvas3dContext?.canvas,
    });
    if (!this._initViewer(container.canvas, container.parent, options?.canvas3dContext)) {
      return false;
    }
    if (options?.checkeredCanvasBackground) {
      this.canvas3d?.setProps({ checkeredTransparentBackground: true });
    }
    this.container = container;
    return true;
  }

  /**
   * Mount the plugin into the target element (assumes the target has "relative"-like positioninig).
   * If initContainer wasn't called separately before, initOptions will be passed to it.
   */
  private _mount(
    target: HTMLElement,
    initOptions?: { canvas3dContext?: Canvas3DContext; checkeredCanvasBackground?: boolean },
  ) {
    if (this.disposed) throw new Error('Cannot mount a disposed context');

    if (!this._initContainer(initOptions)) return false;
    this.container?.mount(target);
    this.handleResize();
    return true;
  }

  unmount() {
    this.container?.unmount();
  }

  private _initViewer(canvas: HTMLCanvasElement, container: HTMLDivElement, canvas3dContext?: Canvas3DContext) {
    try {
      this.layout.setRoot(container);
      if (this.spec.layout && this.spec.layout.initial) this.layout.setProps(this.spec.layout.initial);

      if (!canvas3dContext) {
        canvas3dContext = Canvas3DContext.fromCanvas(
          canvas,
          this.managers.asset,
          {
            antialias: !(this.config.get(PluginConfig.General.DisableAntialiasing) ?? false),
            preserveDrawingBuffer: !(this.config.get(PluginConfig.General.DisablePreserveDrawingBuffer) ?? false),
            preferWebGl1: this.config.get(PluginConfig.General.PreferWebGl1) || false,
            failIfMajorPerformanceCaveat: !(this.config.get(PluginConfig.General.AllowMajorPerformanceCaveat) ?? false),
            powerPreference: this.config.get(PluginConfig.General.PowerPreference) || 'high-performance',
            handleResize: this.handleResize,
          },
          {
            pixelScale: this.config.get(PluginConfig.General.PixelScale) || 1,
            pickScale: this.config.get(PluginConfig.General.PickScale) || 0.25,
            transparency: this.config.get(PluginConfig.General.Transparency) || 'wboit',
            resolutionMode: this.config.get(PluginConfig.General.ResolutionMode) || 'auto',
          },
        );
      }
      (this.canvas3dContext as Canvas3DContext) = canvas3dContext;
      (this.canvas3d as Canvas3D) = Canvas3D.create(this.canvas3dContext!);
      this.canvas3dInit.next(true);
      let props = this.spec.canvas3d;

      const backgroundColor = Color(0xfcfbf9);
      if (!props) {
        this.canvas3d?.setProps({ renderer: { backgroundColor } });
      } else {
        if (props.renderer?.backgroundColor === void 0) {
          props = produce(props, (p) => {
            if (p.renderer) p.renderer.backgroundColor = backgroundColor;
            else p.renderer = { backgroundColor };
          });
        }
        this.canvas3d?.setProps(props);
      }
      this.animationLoop.start();
      (this.helpers.viewportScreenshot as ViewportScreenshotHelper) = new ViewportScreenshotHelper(this);

      this.subs.push(this.canvas3d!.interaction.click.subscribe((e) => this.behaviors.interaction.click.next(e)));
      this.subs.push(this.canvas3d!.interaction.drag.subscribe((e) => this.behaviors.interaction.drag.next(e)));
      this.subs.push(this.canvas3d!.interaction.hover.subscribe((e) => this.behaviors.interaction.hover.next(e)));
      this.subs.push(
        this.canvas3d!.input.resize.pipe(
          debounceTime(50),
          throttleTime(100, undefined, { leading: false, trailing: true }),
        ).subscribe(() => this.handleResize()),
      );
      this.subs.push(this.canvas3d!.input.keyDown.subscribe((e) => this.behaviors.interaction.key.next(e)));
      this.subs.push(this.canvas3d!.input.keyUp.subscribe((e) => this.behaviors.interaction.keyReleased.next(e)));
      this.subs.push(
        this.canvas3d!.xr.isPresenting.subscribe((e) => this.log.info(`WebXR ${e ? 'enabled' : 'disabled'}`)),
      );
      this.subs.push(this.canvas3d!.xr.requestFailed.subscribe((e) => this.log.error(`WebXR request failed: ${e}`)));
      this.subs.push(this.layout.events.updated.subscribe(() => requestAnimationFrame(() => this.handleResize())));

      this.handleResize();

      Scheduler.setImmediate(() => this.initCanvas3dPromiseCallbacks[0]());
      return true;
    } catch (e) {
      this.log.error('' + e);
      console.error(e);
      Scheduler.setImmediate(() => this.initCanvas3dPromiseCallbacks[1](e));
      return false;
    }
  }

  handleResize = () => {
    const canvas = this.canvas3dContext?.canvas;
    const container = this.layout.root;
    if (container && canvas) {
      resizeCanvas(canvas, container, this.canvas3dContext.pixelScale);
      this.canvas3dContext.syncPixelScale();
      this.canvas3d?.requestResize();
    }
  };

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
  readonly fetch = ajaxGet;

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

  dataTransaction(
    f: (ctx: RuntimeContext) => Promise<void> | void,
    options?: { canUndo?: string | boolean; rethrowErrors?: boolean },
  ) {
    return this.runTask(this.state.data.transaction(f, options));
  }

  clear(resetViewportSettings = false) {
    if (resetViewportSettings) this.canvas3d?.setProps(DefaultCanvas3DParams);
    return PluginCommands.State.RemoveObject(this, { state: this.state.data, ref: StateTransform.RootRef });
  }

  dispose(options?: { doNotForceWebGLContextLoss?: boolean; doNotDisposeCanvas3DContext?: boolean }) {
    if (this.disposed) return;

    for (const s of this.subs) {
      s.unsubscribe();
    }
    this.subs = [];

    this.layout.dispose();
    this.managers.markdownExtensions.audio.dispose();
    this.animationLoop.stop();
    this.commands.dispose();
    this.canvas3d?.dispose();
    if (!options?.doNotDisposeCanvas3DContext) {
      this.canvas3dContext?.dispose(options);
    }
    this.ev.dispose();
    this.state.dispose();
    this.helpers.substructureParent.dispose();

    objectForEach(this.managers, (m) => (m as any)?.dispose?.());
    objectForEach(this.managers.structure, (m) => (m as any)?.dispose?.());
    objectForEach(this.managers.volume, (m) => (m as any)?.dispose?.());

    this.unmount();
    this.container = undefined;
    (this.customState as any) = {};

    this.disposed = true;
  }

  private initBehaviorEvents() {
    this.subs.push(
      merge(this.state.data.behaviors.isUpdating, this.state.behaviors.behaviors.isUpdating).subscribe((u) => {
        if (this.behaviors.state.isUpdating.value !== u) this.behaviors.state.isUpdating.next(u);
      }),
    );

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

    this.subs.push(
      merge(this.behaviors.state.isUpdating, this.behaviors.state.isAnimating).subscribe((v) => {
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
      }),
    );

    this.subs.push(
      this.behaviors.interaction.selectionMode.subscribe((v) => {
        if (!v) {
          this.managers.interactivity?.lociSelects.deselectAll();
        }
      }),
    );
  }

  private initBuiltInBehavior() {
    BuiltInPluginBehaviors.State.registerDefault(this);
    BuiltInPluginBehaviors.Representation.registerDefault(this);
    BuiltInPluginBehaviors.Camera.registerDefault(this);
    BuiltInPluginBehaviors.Misc.registerDefault(this);

    this.subs.push(
      merge(this.state.data.events.log, this.state.behaviors.events.log).subscribe((e) => this.events.log.next(e)),
    );
  }

  private async initBehaviors() {
    let tree = this.state.behaviors.build();

    for (const cat of Object.keys(PluginBehavior.Categories)) {
      tree
        .toRoot()
        .apply(
          PluginBehavior.CreateCategory,
          { label: (PluginBehavior.Categories as any)[cat] },
          { ref: cat, state: { isLocked: true } },
        );
    }

    // Init custom properties 1st
    for (const b of this.spec.behaviors) {
      const cat = PluginBehavior.getCategoryId(b.transformer);
      if (cat !== 'custom-props') continue;

      tree
        .to(PluginBehavior.getCategoryId(b.transformer))
        .apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
    }
    await this.runTask(this.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true, doNotLogTiming: true }));

    tree = this.state.behaviors.build();
    for (const b of this.spec.behaviors) {
      const cat = PluginBehavior.getCategoryId(b.transformer);
      if (cat === 'custom-props') continue;

      tree
        .to(PluginBehavior.getCategoryId(b.transformer))
        .apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
    }
    await this.runTask(this.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true, doNotLogTiming: true }));
  }

  /**
   * Registers the providers of one or more entries, in the fixed order of the registry contract, and returns an
   * idempotent function that removes exactly what this call registered.
   *
   * Checks every provider against existing registrations and against each other first; on a conflict it throws one
   * error listing all of them and changes nothing. Must be called after `init()` has created the managers, which is
   * the case for `spec.registry`, behaviors, and everything that runs later.
   */
  register(entry: PluginRegistryEntry | readonly PluginRegistryEntry[]): () => void {
    if (!this.managers.interactivity || !this.managers.lociLabels || !this.builders.structure) {
      throw new Error('PluginContext.register called before init()');
    }
    return registerEntries(this, entry);
  }

  /**
   * Transitional (removed with the constructor preloads in plugin composition step 3): the format registry still
   * preloads every built-in provider, so a registry entry that lists a different provider under a built-in name, such
   * as the Viewer's `customFormats` override, would conflict with the preloaded one. Drop the preloaded provider so the
   * listed one replaces it, as `customFormats` did in 5.x.
   */
  private dropOverriddenPreloadedFormats(registry: readonly PluginRegistryEntry[]) {
    for (const entry of registry) {
      for (const provider of entry.formats ?? []) {
        if (this.dataFormats.has(provider.name) && this.dataFormats.get(provider.name) !== provider) {
          this.dataFormats.remove(provider.name);
        }
      }
    }
  }

  async init() {
    try {
      this.subs.push(this.events.log.subscribe((e) => (this.log.entries = this.log.entries.push(e))));

      this.initBehaviorEvents();
      this.initBuiltInBehavior();

      (this.managers.interactivity as InteractivityManager) = new InteractivityManager(this);
      (this.managers.lociLabels as LociLabelManager) = new LociLabelManager(this);
      (this.builders.structure as StructureBuilder) = new StructureBuilder(this);

      if (this.spec.registry?.length) {
        this.dropOverriddenPreloadedFormats(this.spec.registry);
        this.register(this.spec.registry);
      }

      await this.initBehaviors();

      this.log.message(`Mol* Plugin ${PLUGIN_VERSION}`);
      if (!isProductionMode) this.log.message(`Development mode enabled`);
      if (isDebugMode) this.log.message(`Debug mode enabled`);

      this._isInitialized = true;
      this.initializedPromiseCallbacks[0]();
    } catch (err) {
      this.initializedPromiseCallbacks[1](err);
      throw err;
    }
  }

  constructor(public spec: PluginSpec) {
    for (const key of REMOVED_SPEC_KEYS) {
      if ((spec as Record<string, unknown>)[key] !== undefined) {
        throw new Error(`PluginSpec.${key} was removed in 6.0; use registry entries (see the migration guide)`);
      }
    }
    setSaccharideCompIdMapType(this.config.get(PluginConfig.Structure.SaccharideCompIdMapType) ?? 'default');
  }
}
