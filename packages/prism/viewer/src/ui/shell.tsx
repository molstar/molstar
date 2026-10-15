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
  UIThemes,
} from '@molstar/prism-ui';
import { ownsKey } from '@molstar/prism-ui/dom/events';
// The reset and the theme tokens: every prism-ui root needs them (MspElement imports them; the shell is not one).
import '@molstar/prism-ui/styles/base';
import { adoptStyles, defineStyles, releaseStyles } from '@molstar/prism-ui/styles/registry';
import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import type { HomeModel } from '../model/home.js';
import { HomeTab } from '../views/home/home-tab.js';
import { extensionUIPort } from './extensions.js';
import { historyTab } from './history.js';
import { exportCard, selectionSections, selectionTab } from './panels.js';
import { settingsTab } from './settings.js';
import { LegendCard } from '../sequence/legend-card.js';
import { LegendOverlay, sequenceCommands, sequenceMenu } from '../sequence/commands.js';
import { SequenceControl } from '../sequence/sequence-control.js';
import { ModeSwitch, modeSwitchProps, PickLevel, SelectionBar, selectionModel } from './interaction.js';
import { StructureTab } from '../views/structure/structure-tab.js';
import type { ViewerCore } from '../viewer.js';
import { CanvasBindings } from './bindings.js';
import { createTrajectoryControls } from './trajectory.js';

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
  const theme = new ThemeModel(spec.ui?.theme ?? 'light');
  const stopTheme = theme.applyTo(root);
  const fullscreen = new FullscreenModel(() => target);

  // The canvas: the plugin is mounted here once the shell is in the page and has a size.
  const canvasHost = (<div style="position: absolute; inset: 0" data-msp-canvas="" />) as HTMLElement;

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
        ...(core.options.ui.showHistory ? [historyTab(managers, person)] : []),
        settingsTab(sources, managers, person),
        {
          id: 'help',
          title: 'Help',
          icon: 'circle-help',
          order: 2.8,
          place: 'foot',
          render: domView(() => <HelpPanel commands={managers.commands} bindings={CanvasBindings} />),
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
        {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'fullscreen',
          title: 'Fullscreen',
          pressed: fullscreen.active,
          run: () => fullscreen.toggle(),
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
  unregister.push(
    core.attachExtensionUI(extensionUIPort(shell, managers.commands, { api: person, sources, managers })),
  );

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
            {() => <HoverLabel model={{ lines: hoverLines.asReadonly }} area={() => canvasHost.parentElement} />}
          </Show>
        }
        playback={<Show when={shown('trajectory')}>{() => <TrajectoryBar model={trajectory.model} />}</Show>}
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
            {() => (
              <CommandMenu
                manager={managers.commands}
                group={commandMenuGroup.value}
                onClose={() => {
                  commandMenuOpen.set(false);
                  commandMenuGroup.set(undefined);
                }}
              />
            )}
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
      ink.dispose();
      stopTheme();
      fullscreen.dispose();
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
