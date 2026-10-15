/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type CanvasBackgroundProps,
  CommandMenu,
  type CommandMenuAssistant,
  canvasInk,
  domView,
  FullscreenModel,
  formatChord,
  HelpPanel,
  HoverLabel,
  isMac,
  mapState,
  LogView,
  MolstarWordmark,
  observableState,
  render,
  SectionList,
  ShellModel,
  ShellView,
  Show,
  ThemeModel,
  TrajectoryBar,
  type UITheme,
  UIThemes,
} from '@molstar/prism-ui';
import { ownsKey } from '@molstar/prism-ui/dom/events';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
// The reset and the theme tokens: every prism-ui root needs them (MspElement imports them; the shell is not one).
import '@molstar/prism-ui/styles/base';
import { adoptStyles, defineStyles, releaseStyles } from '@molstar/prism-ui/styles/registry';
import { combineLatest, debounceTime, distinctUntilChanged, map } from 'rxjs';
import { createAnnotationOffer, type HomeModel } from '../model/home.js';
import { HomeTab } from '../views/home/home-tab.js';
import { extensionUIPort } from './extensions.js';
import { historyTab, keepSessionRecord } from './history.js';
import { exportCard, selectionSections, selectionTab } from './panels.js';
import { settingsTab } from './settings.js';
import { LegendCard } from '../sequence/legend-card.js';
import { LegendOverlay, sequenceCommands, sequenceMenu } from '../sequence/commands.js';
import { SequenceControl } from '../sequence/sequence-control.js';
import { ModeSwitch, modeSwitchProps, PickLevel, SelectionBar, selectionModel } from './interaction.js';
import { StructureTab } from '../views/structure/structure-tab.js';
import { HelpProse } from '../views/help/help-prose.js';
import type { ViewerCore } from '../viewer.js';
import { CanvasBindings } from './bindings.js';
import { createTrajectoryControls } from './trajectory.js';
import { downloadBlob } from './download.js';
import { watchToasts } from '../sources/toasts.js';
import { createXrSource } from '../sources/xr.js';

// The viewer on a page (viewer-spec.md §6, plan 1.4): the shell over the viewer's managers, rendered into the target's
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

  const layout = spec.ui?.layout ?? {};
  const shell = new ShellModel({
    ...layout,
    viewState: managers.viewState,
    navigation: managers.navigation,
  });
  if (layout.compactRail) shell.setCompactRail(true);
  if (layout.mode) shell.setMode(layout.mode);
  // The theme the person chose last, kept in the preferences when the spec gives them (journal.md §6).
  const preferences = spec.persist?.preferences;
  const kept = preferences?.get('theme');
  const theme = new ThemeModel(UIThemes.some((t) => t.id === kept) ? (kept as UITheme) : (spec.ui?.theme ?? 'light'));
  const stopTheme = theme.applyTo(root);
  const keepTheme = theme.theme.subscribe((t) => preferences?.set('theme', t));
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

  // The canvas: the plugin is mounted here once the shell is in the page and has a size.
  // Its text alternative (README §3.7): what is loaded and what is selected, in words, for a screen reader; kept in
  // step with the scene, coalesced (one update per settle), and never announced (a description, not a live region).
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
      canvasText.textContent = [loaded, selected ? `Selected: ${selected}.` : '', focused ? `Focused: ${focused}.` : '']
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

  const home: HomeModel = {
    actions: core.views.homeActions,
    // After a load, what a structure carries beyond the built-ins, offered in one click (parity.md §6a).
    annotationOffer: createAnnotationOffer(person, sources.structures),
    formats: sources.formats,
    openFiles,
    session: {
      isEmpty: sources.empty,
      summary: mapState(sources.structures, summaryOf, Object.is),
      clear: () => void person.session.clear(),
    },
  };

  const commandMenuOpen = observableState(false);
  const commandMenuGroup = observableState<string | undefined>(undefined);
  const openCommandMenu = (group?: string) => {
    commandMenuGroup.set(group);
    commandMenuOpen.set(true);
  };

  const unregister: (() => void)[] = [];
  // The core's own parts answer to overrides too (viewer-spec.md §4, `extensions/overrides.ts`): `core.tab.home`,
  // `core.toolbar.reset`, … removed with null, or their fields replaced.
  const own = <T extends { id: string }>(kind: string, parts: readonly T[]): T[] => core.overrides.core(kind, parts);
  unregister.push(
    shell.registerTabs(
      own('tab', [
        {
          id: 'home',
          compactBar: true,
          title: 'Home',
          heading: () => <MolstarWordmark live />,
          icon: 'prism',
          order: 0,
          place: 'head',
          render: domView(({ ctx }) => (
            <HomeTab
              home={home}
              onOpenCommandMenu={() => openCommandMenu()}
              // Features' import sections (the loaders' database and URL forms) after the registry's actions.
              more={() => <SectionList shell={ctx.shell} slot="home.import" />}
            />
          )),
        },
        {
          id: 'structure',
          compactBar: true,
          title: 'Structure',
          icon: 'structure',
          order: 1,
          render: domView(({ ctx }) => (
            <StructureTab
              shell={ctx.shell}
              scene={core.views.scene}
              presets={core.views.presets}
              options={core.views.componentOptions}
              quickStyles={core.views.quickStyles}
              onAction={(action) => {
                if (action === 'add') managers.commands.open('component.add');
              }}
            />
          )),
        },
        {
          id: 'log',
          title: 'Log',
          icon: 'scroll-text',
          order: 2.5,
          place: 'foot',
          render: domView(() => <LogView status={managers.notifications} />),
        },
        // `options.ui.showHistory` (viewer-spec.md §2.2) leaves History out.
        ...(core.options.ui.showHistory ? [historyTab(core, person)] : []),
        settingsTab(sources, managers, person, core.storage),
        {
          id: 'help',
          title: 'Help',
          icon: 'circle-help',
          order: 2.8,
          place: 'foot',
          render: domView(() => (
            <HelpPanel commands={managers.commands} bindings={CanvasBindings} intro={() => <HelpProse />} />
          )),
        },
      ]),
    ),
  );

  const hasData = sources.empty.pipe(map((e) => !e));
  const sequence = core.views.sequence;
  const seqLines = observableState<readonly string[]>([]);
  const hoverLines = observableState<readonly string[]>([]);
  const mergeHover = () => hoverLines.set(seqLines.value.length ? seqLines.value : sources.pointer.hoverLines.value);
  const seqHoverSub = seqLines.asReadonly.subscribe(mergeHover);
  const pointerHoverSub = sources.pointer.hoverLines.subscribe(mergeHover);
  unregister.push(() => {
    seqHoverSub.unsubscribe();
    pointerHoverSub.unsubscribe();
  });
  // `options.ui.showSequence` leaves the band out; the sequence (its API, its state) stays.
  if (sequence && core.options.ui.showSequence) {
    const viewer = sequence.viewers.get();
    const inkToSeq = ink.subscribe((i) => sequence.setInk(i));
    unregister.push(() => inkToSeq.unsubscribe());
    unregister.push(core.registerCommands(sequenceCommands(viewer, shell, person)));
    unregister.push(
      shell.registerOverlays([
        ...own('overlay', [
          { id: LegendOverlay, title: 'Sequence colors', render: domView(() => <LegendCard viewer={viewer} />) },
        ]),
      ]),
    );
    // The band, there while a structure is shown (prism-sequence.md §2.2); the first one shows it.
    let dropTab: (() => void) | undefined;
    let shownBefore = false;
    const bandSub = viewer.structure.subscribe((shown) => {
      if (!shown) {
        dropTab?.();
        dropTab = undefined;
        return;
      }
      if (dropTab) return;
      const tabs = own<Parameters<typeof shell.registerTabs>[0][number]>('tab', [
        {
          id: 'sequence',
          title: 'Sequence',
          icon: 'sequence',
          dock: 'top',
          frame: 'none',
          render: domView(({ ctx }) => (
            <SequenceControl
              viewer={viewer}
              view3d={(s) => sequence.source.view3d(s, () => void person.camera.reset())}
              shell={ctx.shell}
              dock={ctx.dock ?? 'top'}
              onHover={(lines) => seqLines.set(lines)}
            />
          )),
        },
      ]);
      if (!tabs.length) return;
      dropTab = shell.registerTabs(tabs);
      if (!shownBefore) shell.show('sequence');
      shownBefore = true;
    });
    unregister.push(() => bandSub.unsubscribe());
    unregister.push(() => dropTab?.());
    unregister.push(
      shell.registerToolbar(
        own('toolbar', [
          {
            id: 'sequence',
            place: 'canvas',
            group: 'view',
            icon: 'sequence',
            title: 'Sequence',
            visible: viewer.structure.pipe(map((s) => s !== undefined)),
            menu: sequenceMenu(viewer, shell, (id) => managers.commands.open(id)),
          },
        ]),
      ),
    );
  }
  // The core's homes for features' sections: the Selection tab and the Export card (ui/panels.tsx).
  unregister.push(shell.registerTabs(own('tab', [selectionTab(sources.interaction.selection)])));
  unregister.push(shell.registerSections(own('section', selectionSections(person, sources.interaction.selection))));
  const exporting = exportCard(shell);
  unregister.push(shell.registerOverlays(own('overlay', [exporting.overlay])));
  unregister.push(shell.registerToolbar(own('toolbar', [exporting.toolbar])));
  // The trajectory bar over the canvas while the structure in view has several models (`trajectory.*`).
  const trajectory = createTrajectoryControls(core.plugin, person);
  unregister.push(trajectory.dispose);
  const modeProps = modeSwitchProps(sources.interaction, person);
  // The canvas's own controls that `options.ui.buttons` can hide (not toolbar items): while there is data, and not off.
  const shown = (id: string) =>
    combineLatest([hasData, managers.viewState.layout]).pipe(
      map(([data, layout]) => data && layout.buttons[id] !== false),
      distinctUntilChanged(),
    );
  const selection = selectionModel(sources.interaction, person, () => shell.show('selection'));
  unregister.push(
    shell.registerToolbar(
      own('toolbar', [
        {
          id: 'theme',
          priority: 100,
          icon: 'sun',
          title: 'Theme',
          face: theme.theme.pipe(
            map((t) => {
              const current = UIThemes.find((u) => u.id === t)!;
              return { icon: current.icon, title: `Theme: ${current.title}` };
            }),
          ),
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
          visible: hasData,
          icon: 'focus',
          title: 'Reset camera',
          run: () => void person.camera.reset(),
          // v5's Reset menu (parity.md §4): the principal axes, or x, y and z.
          menu: [
            {
              id: 'orient-axes',
              label: 'Orient axes',
              icon: 'axis-3d',
              hint: 'Principal axes along the screen',
              run: () => void person.camera.orient(),
            },
            {
              id: 'reset-axes',
              label: 'Reset axes',
              icon: 'rotate-ccw',
              hint: 'x, y and z along the screen',
              run: () => void person.camera.resetAxes(),
            },
          ],
        },
        // Expand (the canvas fills the viewer, `view.expand`), or fullscreen where expanding would only fill an iframe
        // (FullscreenModel.mode); the other beside it in the menu.
        {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: fullscreen.mode === 'fullscreen' ? 'fullscreen' : 'maximize-2',
          title: fullscreen.mode === 'fullscreen' ? 'Fullscreen' : 'Expand',
          pressed: combineLatest([managers.navigation.expanded, fullscreen.active]).pipe(map(([e, f]) => e || f)),
          run: () =>
            void (fullscreen.mode === 'fullscreen'
              ? person.view.fullscreen({ on: 'toggle' })
              : person.view.expand({ on: 'toggle' })),
          menu: [
            {
              id: 'expand-viewport',
              label: 'Expand',
              icon: 'maximize-2',
              hint: 'Hide the panels; Escape brings them back',
              pressed: managers.navigation.expanded,
              run: () => void person.view.expand({ on: 'toggle' }),
            },
            {
              id: 'fullscreen',
              label: 'Fullscreen',
              icon: 'fullscreen',
              hint: 'The viewer fills the screen',
              pressed: fullscreen.active,
              ...(fullscreen.unavailableReason ? { disabled: fullscreen.unavailableReason } : {}),
              run: () => void person.view.fullscreen({ on: 'toggle' }),
            },
          ],
        },
        // Rendering (parity.md §4's toolbar Rendering row): illumination, pressed while on; XR only where WebXR is.
        {
          id: 'illumination',
          place: 'canvas',
          group: 'rendering',
          visible: combineLatest([hasData, sources.canvas]).pipe(
            map(([data, c]) => data && c !== undefined && person.canvas.illumination !== undefined),
            distinctUntilChanged(),
          ),
          icon: 'lightbulb',
          title: 'Illumination',
          pressed: sources.canvas.pipe(
            map((c) => !!c?.illumination),
            distinctUntilChanged(),
          ),
          run: () => void person.canvas.illumination?.({ on: !sources.canvas.value?.illumination }),
        },
        {
          id: 'xr',
          place: 'canvas',
          group: 'rendering',
          visible: combineLatest([hasData, xr.state]).pipe(
            map(([data, x]) => data && x.supported && person.canvas.xr !== undefined),
            distinctUntilChanged(),
          ),
          icon: 'glasses',
          title: 'XR',
          pressed: xr.state.pipe(
            map((x) => x.presenting),
            distinctUntilChanged(),
          ),
          run: () => void person.canvas.xr?.({ on: 'toggle' }),
        },
      ]),
    ),
  );
  unregister.push(
    core.registerCommands([
      {
        id: 'help.show',
        title: 'Keyboard and mouse',
        alias: ['help', 'shortcuts', 'keys', 'bindings', 'controls'],
        category: 'Help',
        icon: 'circle-help',
        shortcut: { key: '?' },
        run: () => {
          shell.show('help');
          return { message: '' };
        },
      },
      {
        id: 'record.keep',
        title: "Keep this session's record",
        alias: ['keep record', 'keep history', 'save record'],
        category: 'Session',
        icon: 'save',
        when: () => (core.record.kept ? 'The record is already kept in this browser' : true),
        run: async () => ({ message: (await keepSessionRecord(core)) ? "This session's record is kept." : '' }),
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
    ]),
  );

  // Extensions' views: the eager ones' now, a lazy one's placeholders until it loads (viewer-spec.md §3).
  const port = extensionUIPort(shell, managers.commands, { api: person, sources, managers });
  unregister.push(core.attachExtensionUI(port));

  // The assistant, when the spec has one (viewer-spec.md §2.4): loaded on first use, with the agent's code (both chunks
  // of their own, with the AI SDK: G8 perf #5), then its tab, the command menu's Ask row, and its turn card under the
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
        let chunks: [Awaited<ReturnType<typeof load>>, typeof import('../agent/session.js')];
        try {
          chunks = await Promise.all([load(), import('../agent/session.js')]);
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
        const [{ createAssistant }, { createViewerAgent }] = chunks;
        if (closed) return undefined;
        // The agents it made, disposed here if it fails to start (their subscriptions to the viewer would outlive it).
        const agents: { dispose(): void }[] = [];
        try {
          const assistant = await createAssistant({
            viewer: core,
            createAgent: (models) => {
              const agent = createViewerAgent({ viewer: core, models });
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
          unregister.push(() => assistant.dispose());
          assistantMenu.set(assistant.commandMenu);
          agentView.set(assistant.agentView);
          return assistant.commandMenu;
        } catch (e) {
          for (const a of agents.splice(0)) a.dispose();
          managers.notifications.post(`The assistant could not start: ${e instanceof Error ? e.message : e}`, 'error');
          return undefined;
        }
      })());
    // The placeholder tab, where the assistant's will be: opening it loads the assistant, whose tab then shows.
    stub = port.registerStub(
      'assistant',
      { parts: [{ id: AssistantTabId, kind: 'tab', title: 'Assistant', icon: 'bot', order: 2.6 }] },
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

  // ⌘K / ⇧⌘P open the command menu; commands' own shortcuts run them; only this viewer's keys.
  const onKey = (e: KeyboardEvent) => {
    if (!ownsKey(root, e)) return;
    const mod = isMac() ? e.metaKey : e.ctrlKey;
    if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
      e.preventDefault();
      openCommandMenu();
    } else if (!commandMenuOpen.value && managers.commands.handleKey(e, isMac())) e.preventDefault();
  };
  document.addEventListener('keydown', onKey);

  // It starts clean: no panel open (the rail and the canvas tools are there), unless the spec names a tab.
  const start = layout.start;
  if (start && !core.options.ui.collapsed) shell.show(start);
  else managers.navigation.setDock('right', { collapsed: true });
  shell.mount();

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
              // (prism-sequence.md §6).
              <HoverLabel
                model={{ lines: hoverLines.asReadonly }}
                area={() => canvasHost.closest<HTMLElement>('[data-msp-viewport]') ?? canvasHost.parentElement}
              />
            )}
          </Show>
        }
        playback={<Show when={shown('trajectory')}>{() => <TrajectoryBar model={trajectory.model} />}</Show>}
        agent={<Show when={agentView}>{(view) => view()}</Show>}
        ink={ink}
        empty={sources.empty}
        context={() => <Show when={hasData}>{() => <SelectionBar selection={selection} />}</Show>}
        canvasLead={<Show when={shown('selection-mode')}>{() => <ModeSwitch {...modeProps} />}</Show>}
        canvasAside={<Show when={hasData}>{() => <PickLevel {...modeProps} />}</Show>}
        onDropFiles={openFiles}
        dropFormats={home.formats}
        onOpenLog={() => shell.show('log')}
        commandMenu={
          <Show when={commandMenuOpen}>
            {() => {
              const menu = (assistant: CommandMenuAssistant | undefined) => (
                <CommandMenu
                  manager={managers.commands}
                  {...(assistant ? { assistant } : {})}
                  group={commandMenuGroup.value}
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
    requestAnimationFrame(() => resolve());
    setTimeout(resolve, 50);
  });
  const mounted = await core.mountCanvas(canvasHost);
  if (!mounted) managers.notifications.post('The canvas could not start: WebGL is not available.', 'error');

  return {
    shell,
    root,
    dispose() {
      document.removeEventListener('keydown', onKey);
      for (const u of unregister.reverse()) u();
      disposeRender();
      shell.unmount();
      shell.retire();
      inkSub.unsubscribe();
      describeCanvas.unsubscribe();
      ink.dispose();
      stopTheme();
      keepTheme.unsubscribe();
      askedFullscreen.unsubscribe();
      browserFullscreen.unsubscribe();
      fullscreen.dispose();
      stopToasts();
      xr.dispose();
      root.remove();
      releaseStyles(shadow);
    },
  };
}

function summaryOf(list: readonly { label: string; atoms: number }[]): string {
  if (!list.length) return 'Nothing loaded';
  if (list.length === 1) return `${list[0].label} · ${list[0].atoms.toLocaleString('en')} atoms`;
  return `${list.length} structures · ${list.reduce((n, s) => n + s.atoms, 0).toLocaleString('en')} atoms`;
}
