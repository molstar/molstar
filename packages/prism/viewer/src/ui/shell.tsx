/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  BandFadePx,
  CanvasStackFadeMs,
  type CanvasBackgroundProps,
  CommandMenu,
  commandItems,
  type CommandMenuAssistant,
  canvasInk,
  domView,
  FullscreenModel,
  For,
  formatChord,
  HoverLabel,
  isMac,
  LogView,
  observableState,
  render,
  ShellModel,
  ShellView,
  Show,
  ThemeModel,
  PlaybackDeck,
  TrajectoryRow,
  type UITheme,
  UIThemes,
  ViewHost,
} from '@molstar/prism-ui';
import { ownsKey } from '@molstar/prism-ui/dom/events';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import type { HoverLine } from '@molstar/prism-ui/model/playback';
import { type SceneNode, sequenceNodeId } from '../model/scene.js';
import type { SelectionAction } from '../model/selection.js';
import { pickHoverLabel } from '../sources/hover-lines.js';
// The reset and the theme tokens: every prism-ui root needs them (MspElement imports them; the shell is not one).
import '@molstar/prism-ui/styles/base';
import { adoptStyles, defineStyles, releaseStyles } from '@molstar/prism-ui/styles/registry';
import {
  combineLatest,
  debounceTime,
  distinctUntilChanged,
  map,
  merge,
  type Observable,
  of,
  shareReplay,
  skip,
  switchMap,
  timer,
} from 'rxjs';
import { extensionUIPort } from './extensions.js';
import { historyTab, keepSessionRecord } from './history.js';
import { exportCard } from './export-card.js';
import { isLazyExtension } from '../spec.js';
import { settingsTab } from './settings.js';
import { keepLayout, layoutCommands, layoutStore, restoreLayout } from './layout.js';
import { ModeSwitch, modeMenuRows, modeSwitchProps, PickLevel, SelectionBar, selectionModel } from './interaction.js';
import { createStyleForms } from '../sources/command-forms.js';
import { styleCommandForms } from '../views/command-forms/style-forms.js';
import type { ViewerCore } from '../viewer.js';
import { createTrajectoryControls } from './trajectory.js';
import { downloadBlob } from './download.js';
import { decodeColor } from '@molstar/core/util/color/utils';
import { FrameGlideMs, frameBelowTop } from '../sources/framing.js';
import { watchToasts } from '../sources/toasts.js';
import { createXrSource } from '../sources/xr.js';

import { createRoot, getScope, runInScope } from '@molstar/prism-ui/dom/scope';
import { frame, on } from '@molstar/prism-ui/dom/scoped';
import { derive } from '@molstar/prism-ui/reactive/state';
// The viewer on a page (viewer-spec.md §6): the shell over the viewer's managers, rendered into the target's
// shadow root (README decision 12), with the plugin's canvas mounted under it. Keys go to this viewer only (decision
// 29); styles are adopted into its shadow root (library.md §3).

defineStyles(
  'prism-viewer.host',
  // No `contain`: it would make the host the containing block of fixed popups (tooltips, menus) placed in viewport
  // coordinates.
  '@layer msp.reset { :host { display: block; position: relative; overflow: hidden; } }',
);

/** Files that open as a session rather than as data. */
const SessionFile = /\.(molx|molj)$/i;

/** The assistant's tab (`@molstar/prism-assistant`'s), stood in for by a placeholder until it loads. */
const AssistantTabId = 'assistant.tab.assistant';

export interface MountedShell {
  readonly shell: ShellModel;
  readonly root: HTMLElement;
  dispose(): void;
}

export async function mountShell(core: ViewerCore, target: HTMLElement): Promise<MountedShell> {
  const { managers, sources, spec } = core;
  const person = core.api.as({ kind: 'person' });

  const shadow = target.shadowRoot ?? target.attachShadow({ mode: 'open' });
  adoptStyles(shadow);
  const root = document.createElement('div');
  root.className = 'msp-root';
  root.style.cssText = 'position: absolute; inset: 0';
  shadow.appendChild(root);
  // The shell's own attachments (its keys, the first frame: rendering.md §2.4) in a root scope of their own, ended with
  // the shell.
  const attached = createRoot(() => getScope()!);
  // The shell's teardown, defined first: on a failure partway through, a step whose part was never made is passed
  // over (each one guarded), so what was made goes; `dispose` runs it whole.
  let teardown: (() => void) | undefined;
  try {
    const steps: (() => void)[] = [
      () => attached.dispose(),
      () => {
        for (const u of unregister.reverse()) u();
      },
      () => disposeRender(),
      () => shell.unmount(),
      () => shell.retire(),
      () => inkSub.unsubscribe(),
      () => describeCanvas.unsubscribe(),
      () => ink.dispose(),
      () => stopTheme(),
      () => keepTheme.unsubscribe(),
      () => askedFullscreen.unsubscribe(),
      () => browserFullscreen.unsubscribe(),
      () => fullscreen.dispose(),
      () => stopToasts(),
      () => xr.dispose(),
      () => root.remove(),
      () => releaseStyles(shadow),
    ];
    teardown = () => {
      for (const step of steps) {
        try {
          step();
        } catch {
          // Not made, or already gone.
        }
      }
    };

    const layout = spec.ui?.layout ?? {};
    const shell = new ShellModel({
      ...layout,
      viewState: managers.viewState,
      navigation: managers.navigation,
    });
    if (layout.compactRail) shell.setCompactRail(true);
    if (layout.mode) shell.setMode(layout.mode);
    // The layout the person left (journal.md §6, `ui/layout.ts`), before the tabs register so each lands where it was;
    // the host's parts of it (mode, compact rail, buttons) stay the spec's.
    const layoutPreferences = layoutStore(spec);
    if (layoutPreferences) restoreLayout(managers, layoutPreferences);
    // The theme the person chose last, kept in the preferences when the spec gives them (journal.md §6).
    const preferences = spec.persist?.preferences;
    const kept = preferences?.get('theme');
    const theme = new ThemeModel(UIThemes.some((t) => t.id === kept) ? (kept as UITheme) : (spec.ui?.theme ?? 'light'));
    const stopTheme = theme.applyTo(root);
    // Kept when the person changes it, not the spec's or the default one the viewer starts with.
    const keepTheme = theme.theme.pipe(skip(1)).subscribe((t) => preferences?.set('theme', t));
    const fullscreen = new FullscreenModel(() => target);
    // `view.fullscreen` asks through `navigation.fullscreen`; the browser's fullscreen follows it, and what the browser
    // did (Escape, a refusal) is written back.
    const askedFullscreen = managers.navigation.fullscreen.subscribe((want) => {
      if (want === fullscreen.active.value) return;
      void fullscreen.toggle().then(() => managers.navigation.setFullscreen(fullscreen.active.value));
    });
    const browserFullscreen = fullscreen.active.subscribe((on) => managers.navigation.setFullscreen(on));
    // The plugin's toasts in the status column (parity.md's "toasts" row).
    const stopToasts = watchToasts(core.plugin, managers.notifications);
    const xr = createXrSource(core.plugin);
    // The command menu's own forms (commands.md §5.1): Show as… and Color by theme… as pickers.
    const commandForms = styleCommandForms(createStyleForms(core.plugin));

    // The canvas: the plugin is mounted here once the shell is in the page and has a size.
    // Its text alternative (accessibility.md): what is loaded and what is selected, in words, for a screen reader; kept
    // in step with the scene, coalesced (one update per settle), and never announced (a description, not a live
    // region).
    const canvasText = (
      <div
        id="msp-canvas-text"
        style="position: absolute; width: 1px; height: 1px; overflow: hidden; clip-path: inset(50%)"
      />
    ) as HTMLElement;
    const canvasHost = (
      <div
        style="position: absolute; inset: 0"
        data-msp-canvas=""
        role="img"
        aria-label="3D view"
        aria-describedby="msp-canvas-text"
      />
    ) as HTMLElement;
    canvasHost.appendChild(canvasText);
    const describeCanvas = combineLatest([sources.structures, sources.interaction.selection, sources.interaction.focus])
      .pipe(debounceTime(250))
      .subscribe(([structures, selected, focused]) => {
        const loaded = structures.length
          ? `${structures.map((st) => st.label).join(', ')} (${structures.reduce((n, st) => n + st.atoms, 0)} atoms).`
          : 'Nothing is loaded.';
        canvasText.textContent = [
          loaded,
          selected ? `Selected: ${selected}.` : '',
          focused ? `Focused: ${focused}.` : '',
        ]
          .filter(Boolean)
          .join(' ');
      });

    // The ink follows the plugin's canvas props (the background); white until the canvas exists.
    const inkProps = observableState<CanvasBackgroundProps>({ renderer: { backgroundColor: 0xffffff } });
    const inkSub = sources.canvas.subscribe((c) => {
      if (c) inkProps.set(c.props as unknown as CanvasBackgroundProps);
    });
    const ink = canvasInk(inkProps.asReadonly);

    const openFiles = (files: readonly File[]) => {
      const sessions = files.filter((f) => SessionFile.test(f.name));
      const data = files.filter((f) => !SessionFile.test(f.name));
      if (sessions.length > 1) {
        managers.notifications.post(`Opened ${sessions[0].name} only: a session replaces what is loaded.`, 'warning');
      }
      if (sessions.length) void person.session.open({ file: sessions[0] });
      if (data.length) void person.load.files({ files: data });
    };

    const commandMenuOpen = observableState(false);
    const commandMenuGroup = observableState<string | undefined>(undefined);
    // Opened for one command from elsewhere: Esc closes it (no list to go back to).
    let commandMenuDirect = false;
    const openCommandMenu = (group?: string, direct = false) => {
      commandMenuGroup.set(group);
      commandMenuDirect = direct;
      commandMenuOpen.set(true);
    };

    const unregister: (() => void)[] = [];
    // The core's own parts answer to overrides too (viewer-spec.md §4, `extensions/overrides.ts`): `core.tab.home`,
    // `core.toolbar.reset`, … removed with null, or their fields replaced.
    const own = <T extends { id: string }>(kind: string, parts: readonly T[]): T[] => core.overrides.core(kind, parts);
    unregister.push(
      shell.registerTabs(
        own('tab', [
          // History and the log in one panel; `options.ui.showHistory` (viewer-spec.md §2.2) leaves History out: the
          // log alone.
          core.options.ui.showHistory
            ? historyTab(core, person, managers.notifications)
            : {
                id: 'log',
                title: 'Log',
                icon: 'scroll-text',
                order: 2.5,
                place: 'foot',
                render: domView(() => <LogView status={managers.notifications} />),
              },
          settingsTab(sources, managers, person, core.storage),
        ]),
      ),
    );
    // The built-in features' panels (viewer-spec.md §3.1): Home, Structure, Selection, Help, with the core's ids; each
    // only when the spec lists its feature.
    for (const e of core.spec.extensions ?? []) {
      if (isLazyExtension(e) || !e.builtIn || !e.shell) continue;
      unregister.push(
        e.shell({ shell, core, person, own, openCommandMenu: (group) => openCommandMenu(group), openFiles }),
      );
    }

    const hasData = derive(sources.empty, (e) => !e);
    // The canvas stack's own: its buttons stay while it fades out as the canvas empties (`CanvasStackFadeMs`).
    const stackData = hasData.pipe(
      switchMap((data) => (data ? of(true) : timer(CanvasStackFadeMs).pipe(map(() => false)))),
      distinctUntilChanged(),
      shareReplay({ bufferSize: 1, refCount: true }),
    );
    const sequence = core.views.sequence;
    const seqLines = observableState<readonly string[]>([]);
    const hoverLines = observableState<readonly HoverLine[]>([]);
    // The label is Mol*'s whenever there is one: a sequence hover the link maps into the structure lights it in 3D, and
    // Mol*'s label providers label that (sequence.md §10). Seq*'s own label is for what does not map (a gap, an
    // annotation track); while a mapped hover waits for Mol*'s label, nothing new is shown (the label holds a frame).
    // Seq*'s lines come a frame after the mapped flag (sequence-control.tsx hands them on per frame): lines seen while
    // mapped stay Mol*'s after the flag drops, so leaving a mapped residue never shows Seq*'s label of it.
    let mappedSeq: readonly string[] | undefined;
    const mergeHover = () => {
      const mol = sources.pointer.hoverLabel.value;
      const seq = seqLines.value;
      const mapped = sequence?.source.hoverMapped.value ?? false;
      if (mapped) mappedSeq = seq;
      hoverLines.set(pickHoverLabel(mol, seq, mapped || seq === mappedSeq));
    };
    const seqHoverSub = seqLines.asReadonly.subscribe(mergeHover);
    const pointerHoverSub = sources.pointer.hoverLabel.subscribe(mergeHover);
    const mappedSub = sequence ? sequence.source.hoverMapped.subscribe(mergeHover) : undefined;
    unregister.push(() => {
      seqHoverSub.unsubscribe();
      pointerHoverSub.unsubscribe();
      if (mappedSub) mappedSub.unsubscribe();
    });
    // `options.ui.showSequence` leaves the band out; the sequence (its API, its state) stays. Its views come with its
    // runtime (its own chunk), once the core has loaded it (a structure shown).
    if (sequence && core.options.ui.showSequence) {
      const inkToSeq = ink.subscribe((i) => sequence.setInk(i));
      unregister.push(() => inkToSeq.unsubscribe());
      let dropUI: (() => void) | undefined;
      let gone = false;
      const uiSub = sequence.viewers.subscribe((viewers) => {
        if (!viewers || dropUI) return;
        void sequence.load().then(({ runtime }) => {
          if (gone || dropUI) return;
          dropUI = runtime.sequenceShellUI({
            viewers,
            source: sequence.source,
            shell,
            person,
            registerCommands: (defs) => core.registerCommands(defs),
            own,
            openCommand: (id) => void managers.commands.open(id),
            onHover: (lines) => seqLines.set(lines),
            marks: sources.marks,
            registerGroups: (groups) => core.managers.commands.registerGroups(groups),
            openCommandMenu: (group) => openCommandMenu(group),
            ...('scene' in core.views && sequence.tree
              ? {
                  reveal: (structure: string) => {
                    shell.show('structure');
                    core.views.scene.reveal?.(sequenceNodeId(structure));
                  },
                }
              : {}),
          });
        });
      });
      unregister.push(() => {
        gone = true;
        uiSub.unsubscribe();
        dropUI?.();
      });
    }
    // The core's home for features' export sections: the Export card (ui/export-card.tsx).
    const exporting = exportCard(shell);
    unregister.push(shell.registerOverlays(own('overlay', [exporting.overlay])));
    unregister.push(shell.registerToolbar(own('toolbar', [exporting.toolbar])));
    // The trajectory bar over the canvas while the structure in view has several models (`trajectory.*`).
    // Only with the `trajectory` feature (viewer-spec.md §3.1).
    const trajectory = person.area('trajectory') ? createTrajectoryControls(core.plugin, person) : undefined;
    if (trajectory) unregister.push(trajectory.dispose);
    const modeProps = modeSwitchProps(sources.interaction, person);
    // The mode choices as rows of the collapsed toolbar's menu (ui shell.md §2); the same buttons as `ModeSwitch`.
    const modeRows = modeMenuRows(modeProps);
    // The canvas's own controls that `options.ui.buttons` can hide (not toolbar items): while there is data, and not off.
    const shown = (id: string, data: Observable<boolean> = hasData) =>
      combineLatest([data, managers.viewState.layout]).pipe(
        map(([data, layout]) => data && layout.buttons[id] !== false),
        distinctUntilChanged(),
      );
    // The followed trajectory's row on screen in the Structure panel (not under a card in the panel's place): the bar
    // steps back to a compact chip then (the row has the controls). Any other way (another tab, the panel closed, the
    // structure collapsed or searched away, the canvas expanded) the full bar is there, playing or not.
    const rowOpen = trajectory
      ? combineLatest([trajectory.rowShown, shell.overlay]).pipe(
          map(([row, overlay]) => row && !overlay),
          distinctUntilChanged(),
        )
      : undefined;
    // The chip stays while it has the focus (its Pause pressed from the keyboard does not drop the focus to the page).
    const chipFocused = observableState(false);
    // The selection bar's menu: the viewer's rows, the extensions' (`ui.selectionActions`, through the port below)
    const extraSelectionActions = observableState<readonly SelectionAction[]>([]);
    // Show in structure tree: the Structure panel with its tree filtered by the selection, the first component that
    // holds it picked (viewer-spec.md §3.1, its search's non-text part)
    const showInTree = () => {
      if (!('scene' in core.views)) return;
      const scene = core.views.scene;
      scene.selectionFilter?.set(true);
      shell.show('structure');
      const first = (nodes: readonly SceneNode[]): SceneNode | undefined => {
        for (const n of nodes) {
          if (n.kind === 'component' && scene.selectionFilter?.holds(n)) return n;
          const below = n.kind === 'structure' || n.kind === 'components' ? first(n.children.value) : undefined;
          if (below) return below;
        }
        return undefined;
      };
      const hit = first(scene.roots.value);
      if (hit) scene.reveal?.(hit.id);
    };
    const selection = selectionModel(sources.interaction, person, {
      openPanel: () => shell.show('selection'),
      ...('scene' in core.views ? { showInTree } : {}),
      extra: extraSelectionActions.asReadonly,
    });
    unregister.push(
      shell.registerToolbar(
        own('toolbar', [
          {
            id: 'theme',
            priority: 100,
            icon: 'sun',
            title: 'Theme',
            face: derive(theme.theme, (t) => {
              const current = UIThemes.find((u) => u.id === t)!;
              return { icon: current.icon, title: `Theme: ${current.title}` };
            }),
            // A choice of the themes beside the rail (five is too many to cycle through); under More, a click cycles.
            menu: UIThemes.map((t) => ({
              id: t.id,
              label: t.title,
              icon: t.icon,
              choice: true,
              pressed: derive(theme.theme, (x) => x === t.id),
              run: () => theme.set(t.id),
            })),
            run: () => theme.next(),
          },
          {
            id: 'palette',
            priority: -1,
            icon: 'command',
            title: 'Command menu',
            caption: formatChord({ key: 'k', mod: true }),
            keepsCard: true,
            run: () => openCommandMenu(),
          },
          {
            id: 'reset',
            place: 'canvas',
            group: 'canvas',
            visible: stackData,
            icon: 'focus',
            title: 'Reset camera',
            run: () => void person.camera.reset(),
            // v5's Reset menu (parity.md §4): the principal axes, or x, y and z. Entries by command id (commands.md §7):
            // whether they can run and their keys are the commands'; a command an override removed takes its entry.
            menu: commandItems(managers.commands, [
              {
                command: 'camera.orient',
                id: 'orient-axes',
                label: 'Orient axes',
                hint: 'Principal axes along the screen',
              },
              {
                command: 'camera.reset-axes',
                id: 'reset-axes',
                label: 'Reset axes',
                hint: 'x, y and z along the screen',
              },
            ]),
          },
          // Expand (the canvas fills the viewer, `view.expand`), or fullscreen where expanding would only fill an iframe
          // (FullscreenModel.mode); the other beside it in the menu.
          {
            id: 'expand',
            place: 'canvas',
            group: 'canvas',
            visible: stackData,
            icon: fullscreen.mode === 'fullscreen' ? 'fullscreen' : 'maximize-2',
            title: fullscreen.mode === 'fullscreen' ? 'Fullscreen' : 'Expand',
            pressed: derive([managers.navigation.expanded, fullscreen.active], (e, f) => e || f),
            run: () =>
              void (fullscreen.mode === 'fullscreen'
                ? person.view.fullscreen({ on: !managers.navigation.fullscreen.value })
                : person.view.expand({ on: !managers.navigation.expanded.value })),
            menu: [
              {
                id: 'expand-viewport',
                label: 'Expand',
                icon: 'maximize-2',
                hint: 'Hide the panels; Escape brings them back',
                pressed: managers.navigation.expanded,
                run: () => void person.view.expand({ on: !managers.navigation.expanded.value }),
              },
              {
                id: 'fullscreen',
                label: 'Fullscreen',
                icon: 'fullscreen',
                hint: 'The viewer fills the screen',
                pressed: fullscreen.active,
                ...(fullscreen.unavailableReason ? { disabled: fullscreen.unavailableReason } : {}),
                run: () => void person.view.fullscreen({ on: !managers.navigation.fullscreen.value }),
              },
            ],
          },
          // Rendering (parity.md §4's toolbar Rendering row): illumination, pressed while on; XR only where WebXR is.
          {
            id: 'illumination',
            place: 'canvas',
            group: 'rendering',
            visible: combineLatest([stackData, sources.canvas]).pipe(
              map(([data, c]) => data && c !== undefined && person.describe('canvas.lighting') !== undefined),
              distinctUntilChanged(),
            ),
            icon: 'lightbulb',
            title: 'Illumination',
            pressed: derive(sources.canvas, (c) => !!c?.illumination),
            // The lighting's illumination switch, set to the other of what it is (api.md §1).
            run: () =>
              void person.canvas.lighting({
                illumination: { name: sources.canvas.value?.illumination ? 'off' : 'on' },
              }),
          },
          {
            id: 'xr',
            place: 'canvas',
            group: 'rendering',
            // `options.ui.xr` (viewer-spec.md §2.2): auto, only where WebXR is; always, with a note where it is not.
            // Either way only while `canvas.xr` is there to run.
            visible: combineLatest([stackData, xr.state]).pipe(
              map(([data, x]) => {
                const mode = core.options.ui.xr;
                if (!data || mode === 'never' || person.describe('canvas.xr') === undefined) return false;
                return mode === 'always' || x.supported;
              }),
              distinctUntilChanged(),
            ),
            icon: 'glasses',
            title: 'XR',
            pressed: derive(xr.state, (x) => x.presenting),
            run: () => {
              if (!xr.state.value.supported) {
                managers.notifications.post(
                  'XR is not available: this browser or device does not support WebXR.',
                  'warning',
                );
                return;
              }
              void person.canvas.xr({ on: !xr.state.value.presenting });
            },
          },
        ]),
      ),
    );
    unregister.push(
      core.registerCommands([
        {
          id: 'record.keep',
          title: "Keep this session's record",
          alias: ['keep record', 'keep history', 'save record'],
          category: 'Session',
          icon: 'save',
          when: () => (core.record.kept ? 'The record is already kept in this browser' : true),
          run: async () => {
            // At the bound (journal.md §5) History asks before forgetting the oldest: the command does not forget one.
            const oldest = await core.keptRecords.overBound();
            if (oldest) {
              shell.show('history');
              return {
                message: `${core.keptRecords.max} records are kept: keep this one from History, which asks before forgetting “${oldest.title}”.`,
              };
            }
            return { message: (await keepSessionRecord(core)) ? "This session's record is kept." : '' };
          },
        },
        {
          id: 'view.theme',
          title: 'Next theme',
          alias: ['theme', 'dark', 'light'],
          category: 'View',
          icon: 'sun-moon',
          run: () => {
            theme.next();
            return { message: '' };
          },
        },
        ...layoutCommands(shell),
      ]),
    );

    // Extensions' views: the eager ones' now, a lazy one's placeholders until it loads (viewer-spec.md §3).
    const port = extensionUIPort(
      shell,
      managers.commands,
      {
        api: person,
        sources,
        managers,
        ...(core.spec.persist?.preferences ? { preferences: core.spec.persist.preferences } : {}),
      },
      extraSelectionActions,
      'scene' in core.views ? core.views.scene.nodeActions : undefined,
    );
    unregister.push(core.attachExtensionUI(port));
    // An extension that did not load, or whose views did not: said, not silent (its placeholder offers Retry).
    const extensionFailed = core.extensions.events.failed.subscribe((f) =>
      managers.notifications.post(`${f.id} could not load: ${f.message}`, 'error'),
    );
    unregister.push(() => extensionFailed.unsubscribe());

    // The assistant, when the spec has one (viewer-spec.md §2.4): loaded on first use, with the agent's code (both
    // chunks of their own, with the AI SDK), then its tab, the command menu's Ask row, and its turn card under the
    // status line. Until then its tab is a placeholder whose opening loads it, and the Ask row loads it and then asks.
    const assistantMenu = observableState<CommandMenuAssistant | undefined>(undefined);
    const agentView = observableState<(() => Child) | undefined>(undefined);
    let closed = false;
    unregister.push(() => {
      closed = true;
    });
    const agentSpec = core.spec.agent;
    if (agentSpec && agentSpec.assistant) {
      const load = agentSpec.assistant;
      let stub: (() => void) | undefined;
      const removeStub = () => {
        stub?.();
        stub = undefined;
      };
      let loading: Promise<CommandMenuAssistant | undefined> | undefined;
      const loadAssistant = () =>
        (loading ??= (async (): Promise<CommandMenuAssistant | undefined> => {
          // The assistant's module, with the agent's code (its own chunk; the viewer does not import it).
          let loaded: Awaited<ReturnType<typeof load>>;
          try {
            loaded = await load();
          } catch (e) {
            // A chunk that does not load (a deploy since the page loaded, journal.md §7): the next use tries again.
            loading = undefined;
            console.error('[prism-viewer] the assistant did not load:', e);
            managers.notifications.post(
              'The assistant could not load. A newer version of the viewer is available: reload.',
              'error',
            );
            return undefined;
          }
          const { createAssistant, createViewerAgent } = loaded;
          if (closed) return undefined;
          // The agents it made, disposed here if it fails to start (their subscriptions to the viewer would outlive it).
          const agents: { dispose(): void }[] = [];
          try {
            const assistant = await createAssistant({
              viewer: core,
              createAgent: (models, options) => {
                const agent = createViewerAgent({
                  viewer: core,
                  models,
                  ...(options?.rules ? { rules: options.rules } : {}),
                });
                agents.push(agent);
                return agent;
              },
            });
            if (closed) {
              assistant.dispose();
              return undefined;
            }
            // Its real tab takes the placeholder's place (same id).
            removeStub();
            unregister.push(port.register('assistant', assistant.ui));
            if (assistant.commands?.length) unregister.push(core.registerCommands(assistant.commands));
            unregister.push(() => assistant.dispose());
            assistantMenu.set(assistant.commandMenu);
            agentView.set(assistant.agentView);
            return assistant.commandMenu;
          } catch (e) {
            for (const a of agents.splice(0)) a.dispose();
            managers.notifications.post(
              `The assistant could not start: ${e instanceof Error ? e.message : e}`,
              'error',
            );
            // Not kept: the placeholder's Retry, or the next ask, tries again.
            loading = undefined;
            return undefined;
          }
        })());
      // The placeholder tab, where the assistant's will be: opening it loads the assistant, whose tab then shows.
      stub = port.registerStub(
        'assistant',
        { parts: [{ id: AssistantTabId, kind: 'tab', title: 'Assistant', icon: 'bot', order: 2.65, place: 'foot' }] },
        () => loadAssistant().then((a) => !!a),
      );
      unregister.push(removeStub);
      // The Ask row before it loads: whether a model is set up is not known yet, so it offers to ask; asking loads it,
      // then asks (or opens the setup when no model is).
      assistantMenu.set({
        ready: () => true,
        propose: () => undefined,
        ask: (text) =>
          void loadAssistant().then((a) => {
            if (!a) return;
            if (a.ready()) a.ask(text);
            else a.setup();
          }),
        setup: () => void loadAssistant().then((a) => a?.setup()),
      });
    }

    // A file the agent made for the person (an export, a saved session: a result's `{ blob, filename }`) is saved here, the
    // page's job; the person's own commands and views save theirs (agent.md §4: such a call asked first).
    const saving = core.api.events.finished.subscribe((r) => {
      // The outer run only (a stub's replay of the real one is nested): one file saved.
      if (r.origin.kind !== 'agent' || r.nested || !r.result?.ok) return;
      const file = r.result.value as { blob?: unknown; filename?: unknown } | undefined;
      if (file?.blob instanceof Blob && typeof file.filename === 'string') downloadBlob(file.blob, file.filename);
    });
    unregister.push(() => saving.unsubscribe());

    // A command opened from elsewhere (a component's action, a shortcut of a command with a form, an extension) shows in
    // the menu: an open command with no menu around it would be invisible.
    const showActive = managers.commands.active.subscribe((invocation) => {
      if (invocation && !commandMenuOpen.value) openCommandMenu(undefined, true);
    });
    unregister.push(() => showActive.unsubscribe());

    // ⌘K / ⇧⌘P open the command menu; commands' own shortcuts run them; only this viewer's keys.
    const onKey = (e: KeyboardEvent) => {
      if (!ownsKey(root, e)) return;
      const mod = isMac() ? e.metaKey : e.ctrlKey;
      if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
        e.preventDefault();
        openCommandMenu();
      } else if (!commandMenuOpen.value && managers.commands.handleKey(e, isMac())) e.preventDefault();
    };
    runInScope(attached, () => on(document, 'keydown', onKey));

    // It starts clean: no panel open (the rail and the canvas tools are there), unless the spec names a tab, or
    // `homeWhenEmpty` opens Home while nothing is loaded (where to load from); not in the embedded layout, whose rail
    // starts hidden and whose page drives what loads (shell.md §3).
    const homeAtStart = core.options.ui.homeWhenEmpty && sources.empty.value && layout.mode !== 'embedded';
    const start = layout.start ?? (homeAtStart ? 'home' : undefined);
    // A tab the person hid does not open at start (`show` would bring it back).
    const startHidden = !!start && managers.viewState.layout.value.hiddenTabs.includes(start);
    if (start && !startHidden && !core.options.ui.collapsed) shell.show(start);
    else managers.navigation.setDock('right', { collapsed: true });
    shell.mount();
    // From here on, what changes is the person's: kept a moment after each change. Stopped first on dispose, before
    // the tabs unregister (that is not a change to keep).
    if (layoutPreferences) steps.unshift(keepLayout(managers, layoutPreferences));

    const disposeRender = render(
      () => (
        <ShellView
          shell={shell}
          status={managers.notifications}
          canvas={canvasHost}
          overlays={
            <Show when={hasData}>
              {() => (
                // The viewport, not the canvas: the sequence band floats over it and feeds the same label
                // (sequence.md §10).
                <HoverLabel
                  model={{ lines: hoverLines.asReadonly }}
                  area={() => canvasHost.closest<HTMLElement>('[data-msp-viewport]') ?? canvasHost.parentElement}
                  // Over the canvas or the sequence, not a panel over them (a tree row has its own card).
                  over={(t) => !!t.closest?.('[data-msp-canvas], [data-msp-hover-source]')}
                  placement={core.options.ui.hoverLabel}
                />
              )}
            </Show>
          }
          playback={
            // One card: the extensions' rows (the snapshots', their transition's), then the trajectory's: a compact
            // chip while frames play and its row is on screen (`rowOpen`), else the full bar, playing or not.
            <PlaybackDeck>
              <For each={shell.sections('playback')}>{(s) => <ViewHost renderer={s.render} ctx={{ shell }} />}</For>
              {trajectory ? (
                <>
                  <Show
                    when={shown('trajectory', combineLatest([hasData, rowOpen!]).pipe(map(([d, open]) => d && !open)))}
                  >
                    {() => <TrajectoryRow model={trajectory.model} />}
                  </Show>
                  <Show
                    when={shown(
                      'trajectory',
                      combineLatest([hasData, rowOpen!, trajectory.model.isPlaying, chipFocused.asReadonly]).pipe(
                        map(([d, open, playing, focused]) => d && open && (playing || focused)),
                      ),
                    )}
                  >
                    {() => (
                      <div
                        style={{ display: 'contents' }}
                        onFocusIn={() => chipFocused.set(true)}
                        onFocusOut={(e) => {
                          if (!(e.currentTarget as HTMLElement).contains(e.relatedTarget as Node | null))
                            chipFocused.set(false);
                        }}
                      >
                        <TrajectoryRow model={trajectory.model} compact />
                      </div>
                    )}
                  </Show>
                </>
              ) : null}
            </PlaybackDeck>
          }
          agent={<Show when={agentView}>{(view) => view()}</Show>}
          ink={ink}
          empty={sources.empty}
          splash={core.options.ui.splash !== 'off'}
          {...(core.options.ui.splash !== 'off' ? { splashField: core.options.ui.splash } : {})}
          busy={core.state.working}
          context={() => <Show when={hasData}>{() => <SelectionBar selection={selection} />}</Show>}
          canvasLead={<Show when={shown('selection-mode', stackData)}>{() => <ModeSwitch {...modeProps} />}</Show>}
          canvasLeadMenu={shown('selection-mode', stackData).pipe(map((on) => (on ? modeRows : [])))}
          canvasAside={<Show when={stackData}>{() => <PickLevel {...modeProps} />}</Show>}
          onDropFiles={openFiles}
          dropFormats={sources.formats}
          onOpenLog={() => (core.options.ui.showHistory ? shell.show('history', 'log') : shell.show('log'))}
          commandMenu={
            <Show when={commandMenuOpen}>
              {() => {
                const menu = (assistant: CommandMenuAssistant | undefined) => (
                  <CommandMenu
                    manager={managers.commands}
                    forms={commandForms}
                    {...(assistant ? { assistant } : {})}
                    group={commandMenuGroup.value}
                    direct={commandMenuDirect}
                    prefixes={core.extensions.prefixes.value}
                    onClose={() => {
                      commandMenuOpen.set(false);
                      commandMenuGroup.set(undefined);
                    }}
                  />
                );
                // Keyed on the assistant: one that loads while the menu is open gets its Ask row there.
                return (
                  <Show when={assistantMenu} fallback={() => menu(undefined)}>
                    {(assistant) => menu(assistant)}
                  </Show>
                );
              }}
            </Show>
          }
        />
      ),
      root,
    );

    // The canvas once the shell has laid out (the plugin sizes its canvas from its container).
    // The next frame, or a moment later where there are no frames (a background tab): creation never waits on visibility.
    await new Promise<void>((resolve) => {
      runInScope(attached, () => frame(() => resolve()));
      setTimeout(resolve, 50);
    });
    const mounted = await core.mountCanvas(canvasHost);
    if (!mounted) managers.notifications.post('The canvas could not start: WebGL is not available.', 'error');
    // The canvas takes the theme's color (`--msp-canvas-bg`) as the theme is chosen, if its background is a plain color
    // (no image, gradient or skybox style); at the start too, unless the options gave a background.
    const themeCanvas = () => {
      const c3 = core.plugin.canvas3d;
      if (!c3 || c3.props.postprocessing.background.variant.name !== 'off') return;
      const color = decodeColor(getComputedStyle(root).getPropertyValue('--msp-canvas-bg').trim());
      if (color === undefined || c3.props.renderer.backgroundColor === color) return;
      core.plugin.setCanvas3DProps({ renderer: { backgroundColor: color } });
    };
    if (!core.options.canvas.background) themeCanvas();
    const followTheme = theme.theme.pipe(skip(1)).subscribe(themeCanvas);
    steps.push(() => followTheme.unsubscribe());
    // Resets frame the scene below the sequence band (a frameless top tab; its faded bottom lets the scene through), not
    // under it (sequence.md §8, README decision 28), Mol*'s first framing of a scene too. The band comes a few frames
    // after the scene (Seq* makes its document first): the first framing counts it in while it is coming, so its arrival
    // moves nothing. It is coming while it has not been shown yet (Seq* shows it by itself only the first time: closed,
    // or the sequence left out of the UI, it does not come), and the structure it follows has a polymer to show.
    const canvas3d = core.plugin.canvas3d;
    if (canvas3d) {
      const top = shell.dock('top');
      let bandShown = false;
      const bandComing = () => {
        if (bandShown || !sequence || !core.options.ui.showSequence) return false;
        const host = sequence.source.structure.value;
        return !!host && host.structure.polymerResidueCount > 0;
      };
      const bandPx = () => {
        const active = top.active.value;
        if (active?.frame === 'none' || (!active && bandComing()))
          return Math.max(0, top.state.value.size - BandFadePx);
        return 0;
      };
      const framing = frameBelowTop(canvas3d, canvasHost, bandPx);
      let last = bandPx();
      // a scene framed before this hold (a load the options started) with the band already there
      framing.reframe();
      // The band shown or closed otherwise (by the person, its menu): the scene glides to its new frame while the
      // camera is where the last framing left it; a resize (a drag of its edge) moves nothing.
      const follow = merge(top.state, top.active).subscribe(() => {
        if (top.active.value) bandShown = true;
        const now = bandPx();
        if (now > 0 === last > 0) {
          last = now;
          return;
        }
        last = now;
        framing.reframe(FrameGlideMs);
      });
      steps.push(
        () => follow.unsubscribe(),
        () => framing.undo(),
      );
    }

    return {
      shell,
      root,
      dispose() {
        teardown?.();
      },
    };
  } catch (e) {
    teardown?.();
    root.remove();
    releaseStyles(shadow);
    throw e;
  }
}
