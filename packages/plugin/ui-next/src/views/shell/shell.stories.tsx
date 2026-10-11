/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { combineLatest, map, of } from 'rxjs';
import { story } from '../../_stories/story-root.js';
import { formatChord, isMac } from '../../commands/keys.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { Show } from '../../dom/control-flow.js';
import { onCleanup } from '../../dom/scope.js';
import { domView } from '../../dom/view.js';
import { parseColor, themeForColor } from '../../model/canvas-theme.js';
import { observableState } from '../../reactive/observable-state.js';
import type { HoverLabelModel } from '../../model/playback.js';
import { FullscreenModel } from '../../shell/fullscreen.js';
import { ThemeModel, type UITheme, UIThemes } from '../../shell/theme.js';
import { ShellModel } from '../../shell/shell.js';
import { StatusModel } from '../../shell/status.js';
import { HoverLabel } from '../canvas/hover-label.js';
import {
  type CanvasMode,
  type Granularity,
  ModeSwitch,
  type ModeSwitchProps,
  PickLevel,
} from '../canvas/mode-switch.js';
import { SnapshotCard, TrajectoryBar } from '../canvas/playback.js';
import { SelectionBar } from '../canvas/selection-bar.js';
import { CommandMenu } from '../command-menu/command-menu.js';
import { type SequenceSelector, SequenceStrip } from '../sequence/sequence-strip.js';
import { CanvasPlaceholder } from './canvas-placeholder.js';
import { LogView } from './log-view.js';
import { MockAgentModel, mockAssistant } from '../../_stories/mocks/agent.js';
import { AgentTurnCard } from '../agent/agent-turn.js';
import { AssistantTab } from '../agent/assistant-tab.js';
import { HistoryTab } from '../agent/history-tab.js';
import { MolstarWordmark } from '../../components/icon/wordmark.js';
import { HomeTab } from '../home/home-tab.js';
import { ImageExport } from '../export/image-export.js';
import { RemoteStatesPanel } from '../state/remote-states-panel.js';
import { SnapshotsPanel } from '../state/snapshots-panel.js';
import { StateTreeTab } from '../state-tree/state-tree-tab.js';
import { SettingsTab } from '../settings/settings-tab.js';
import { mockSettings } from '../../_stories/mocks/settings.js';
import { mockMeasurements, registerMockMeasurements } from '../../_stories/mocks/measurements.js';
import { mockSuperposition, registerMockSuperposition } from '../../_stories/mocks/superposition.js';
import { mockStateTree } from '../../_stories/mocks/state-tree.js';
import { mockImageExport } from '../../_stories/mocks/image-export.js';
import { mockSnapshotsModel } from '../../_stories/mocks/snapshots.js';
import { SectionList } from './section-list.js';
import { ShellView } from './shell-view.js';
import { StructureTab } from '../structure/structure-tab.js';
import { mockHome } from '../../_stories/mocks/home.js';
import { mockSnapshots, mockTrajectory } from '../../_stories/mocks/playback.js';
import { registerMockAnimation } from '../../_stories/mocks/animation.js';
import { mockSelection } from '../../_stories/mocks/selection.js';
import { mockCommands, mockPrefixes } from '../../_stories/mocks/commands.js';
import { registerMockCanvasSettings } from '../../_stories/mocks/canvas-settings.js';
import {
  type MockScenePreset,
  mockComponentOptions,
  mockPresets,
  mockQuickStyles,
  mockScene,
} from '../../_stories/mocks/scene.js';

const meta: Meta = { title: 'Shell/Shell', parameters: { layout: 'fullscreen' } };
export default meta;

const Sequence =
  'MALYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGFCMFDMECHKKYGKVWGFYDGQQPVLAITDPDMIKTVLVKECYSVFTNRRPFGPVGFMKSAISIAEDEEWKRLRSLLSPTFTSGKLKEMVPIIAQYGDVLVRNLRREAETGKPVTLKDVFGAYSMDVITSTSFGVNIDSLNNPQDPFVENTKKLLRFDFLDPFFLSITVFPFLIPILEVLNICVFPREVTNFLRKSVKRMKESRLEDTQKHRVDFLQLMIDSQNSKETESHKALSDLELVAQSIIFIFAGYETTSSVLSFIMYELATHPDVQQKLQEEIDAVLPNKAPPTYDTVLQMEYLDMVVNETLRLFPIAMRLERVCKKDVEINGMFIPKGVVVMIPSYALHRDPKYWTEPEKFLPERFSKKNKDNIDPYIYTPFGSGPRNCIGMRFALMNMKLALIRVLQNFSFKPCKETQIPLKLSLGGLLQPEKPVVLKVESRDGTVSGA';
const ThreeLetter: Record<string, string> = {
  A: 'ALA',
  R: 'ARG',
  N: 'ASN',
  D: 'ASP',
  C: 'CYS',
  Q: 'GLN',
  E: 'GLU',
  G: 'GLY',
  H: 'HIS',
  I: 'ILE',
  L: 'LEU',
  K: 'LYS',
  M: 'MET',
  F: 'PHE',
  P: 'PRO',
  S: 'SER',
  T: 'THR',
  W: 'TRP',
  Y: 'TYR',
  V: 'VAL',
}; // prettier-ignore

/**
 * The Export panel: one rail button, its sections contributed to the `export` slot. Built-in sections (image, session)
 * and extension sections (animation, geometry, models) look alike; each is also a command that opens the panel at it.
 */
function mockExport(shell: ShellModel, status: StatusModel) {
  const imageExport = mockImageExport((t, k) => status.post(t, k));
  const rendering = observableState<number | undefined>(undefined);
  // Each section: choose in a form, then one action in the footer (Mol*'s action footer).
  const form = 'display: grid; gap: 8px';
  const sessionForm = new ParamsFormModel({
    kind: PD.Select('session', [
      ['session', 'Session (.molx) · with data'],
      ['state', 'State (.molj) · description only'],
    ]),
  });
  const mp4Form = new ParamsFormModel({
    animation: PD.Select('spin', [
      ['spin', 'Camera spin'],
      ['rock', 'Camera rock'],
      ['snapshots', 'State snapshots'],
    ]),
    resolution: PD.Select('full-hd', [
      ['hd', 'HD (1280 × 720)'],
      ['full-hd', 'Full HD (1920 × 1080)'],
      ['ultra-hd', 'Ultra HD (3840 × 2160)'],
    ]),
    transparent: PD.Boolean(false),
  });
  const geometryForm = new ParamsFormModel({
    format: PD.Select('glb', [
      ['glb', 'GLB (glTF binary)'],
      ['stl', 'STL'],
      ['obj', 'OBJ'],
      ['usdz', 'USDZ (AR)'],
    ]),
  });
  const modelsForm = new ParamsFormModel({
    format: PD.Select('cif', [
      ['cif', 'mmCIF'],
      ['bcif', 'BinaryCIF'],
    ]),
  });
  shell.registerOverlays([
    {
      id: 'export',
      title: 'Export',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" focus={ctx.focus} />),
    },
  ]);
  shell.registerSections([
    {
      id: 'image',
      slot: 'export',
      title: 'Image',
      icon: 'camera',
      order: 0,
      render: domView(() => <ImageExport model={imageExport} />),
    },
    {
      id: 'session',
      slot: 'export',
      title: 'Session & state',
      icon: 'save',
      order: 1,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={sessionForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() =>
              status.post(
                (sessionForm.values.value as any).kind === 'session' ? 'Saved state.molx' : 'Saved state.molj',
                'success',
              )
            }
          />
        </div>
      )),
    },
    // Contributed by the mp4-export extension.
    {
      id: 'animation',
      slot: 'export',
      // Named for what it makes, not the format: more formats may follow.
      title: 'Animation',
      icon: 'film',
      keywords: ['mp4', 'video', 'movie'],
      order: 10,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={mp4Form} />
          <ActionFooter
            icon="film"
            label={rendering.asReadonly.pipe(
              map((r) => (r === undefined ? 'Render' : `Rendering ${Math.round(r * 100)}%`)),
            )}
            busy={rendering.asReadonly.pipe(map((r) => r !== undefined))}
            onReset={() => mp4Form.reset()}
            onApply={() => {
              // Rendering holds the plugin: an overlay task, with Cancel.
              let p = 0;
              rendering.set(0);
              const done = (message: string) => {
                clearInterval(t);
                rendering.set(undefined);
                status.setTasks(status.tasks.value.filter((x) => x.id !== 'mp4'));
                status.post(message, 'success');
              };
              const task = () => ({
                id: 'mp4',
                label: 'Rendering animation.mp4',
                progress: p,
                overlay: true,
                abort: () => done('Rendering cancelled'),
              });
              const put = () => status.setTasks([...status.tasks.value.filter((x) => x.id !== 'mp4'), task()]);
              put();
              const t = setInterval(() => {
                p = Math.min(1, p + 0.04);
                rendering.set(p);
                put();
                if (p >= 1) done('Rendered animation.mp4 (6 s)');
              }, 120);
            }}
          />
        </div>
      )),
    },
    // Contributed by the geo-export extension.
    {
      id: 'geometry',
      slot: 'export',
      title: 'Geometry',
      icon: 'box',
      keywords: ['glb', 'stl', 'obj', 'usdz', 'mesh', '3d print', 'ar'],
      order: 11,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={geometryForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() => status.post(`Exported scene.${(geometryForm.values.value as any).format}`, 'success')}
          />
        </div>
      )),
    },
    // Contributed by the model-export extension.
    {
      id: 'models',
      slot: 'export',
      title: 'Models',
      icon: 'file',
      keywords: ['mmcif', 'bcif', 'cif', 'coordinates'],
      order: 12,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={modelsForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() => status.post(`Saved 1tqn.${(modelsForm.values.value as any).format}`, 'success')}
          />
        </div>
      )),
    },
  ]);
}

// Mol*'s granularities (interactivity), labelled and drawn as picking levels.
const Granularities: readonly Granularity[] = [
  { value: 'element', label: 'Atom', icon: 'pick-atom' },
  { value: 'residue', label: 'Residue', icon: 'pick-residue' },
  { value: 'chain', label: 'Chain', icon: 'pick-chain' },
  { value: 'entity', label: 'Entity', icon: 'pick-entity' },
  { value: 'model', label: 'Model', icon: 'pick-model' },
  { value: 'operator', label: 'Operator', icon: 'copy' },
  { value: 'structure', label: 'Structure', icon: 'structure' },
  { value: 'elementInstances', label: 'Atom', icon: 'pick-atom', group: 'All instances' },
  { value: 'residueInstances', label: 'Residue', icon: 'pick-residue', group: 'All instances' },
  { value: 'chainInstances', label: 'Chain', icon: 'pick-chain', group: 'All instances' },
];

/** The sequence view's choices for 1TQN (Mol*'s sequence controls): structure, mode, then entity and chain. */
function sequenceSelectors(post: (text: string) => void): SequenceSelector[] {
  const selector = (
    id: string,
    title: string,
    initial: string,
    options: [string, string][],
    extra: Partial<SequenceSelector> = {},
  ): SequenceSelector => {
    const value = observableState(initial);
    return {
      id,
      title,
      options,
      value: value.asReadonly,
      set: (v) => {
        value.set(v);
        post(`Sequence: ${title.toLowerCase()} ${options.find((o) => o[0] === v)?.[1] ?? v}`);
      },
      ...extra,
    };
  };
  const mode = selector('mode', 'Mode', 'single', [
    ['single', 'Chain'],
    ['polymers', 'Polymers'],
    ['all', 'Everything'],
  ]);
  const single = mode.value.pipe(map((m) => m === 'single'));
  return [
    selector('structure', 'Structure', '1tqn', [['1tqn', '1TQN | Crystal structure of human microsomal P450 3A4']], {
      display: () => '1TQN',
      strong: true,
    }),
    mode,
    selector(
      'entity',
      'Entity',
      '1',
      [
        ['1', '1: Cytochrome P450 3A4'],
        ['2', '2: Protoporphyrin IX containing Fe'],
        ['3', '3: Sulfate ion'],
        ['4', '4: Water'],
      ],
      {
        display: (v) =>
          ['Cytochrome P450 3A4', 'Protoporphyrin IX containing Fe', 'Sulfate ion', 'Water'][Number(v) - 1],
        visible: single,
      },
    ),
    selector('chain', 'Chain', 'A', [['A', 'A']], { display: (v) => `Chain ${v}`, visible: single }),
  ];
}

/** A mock of everything around the canvas, with the shell's real model. */
function mockShell(options: {
  loaded: boolean;
  preset?: MockScenePreset;
  pinned?: boolean;
  remoteStates?: boolean | number;
  ai?: { configured?: boolean };
}) {
  const shell = new ShellModel({ docks: { right: { pinned: !!options.pinned } } });
  const status = new StatusModel();
  // The agent module, only when the app includes it (design/agent.md §6): the Assistant tab, History, the Ask row, the
  // turn card. A mock agent with scripted turns stands in for any-harness.
  const agent = options.ai
    ? new MockAgentModel((text, kind) => status.post(text, kind), {
        configured: options.ai.configured,
        history: options.loaded,
      })
    : undefined;
  const canvasMode = observableState<CanvasMode>('interact');
  const commandMenuOpen = observableState(false);
  // The group the command menu opens in (Preset, Styles); none for ⌘K.
  const commandMenuGroup = observableState<string | undefined>(undefined);
  const granularity = observableState('residue');
  const illumination = observableState(false);
  // The plugin root goes fullscreen (in the stories, the story page).
  const fullscreen = new FullscreenModel(() => document.documentElement);
  // The UI theme, starting from Storybook's (its toolbar sets the same attribute on the page).
  const theme = new ThemeModel((document.documentElement.dataset.mspTheme as UITheme | undefined) ?? 'light');
  theme.applyTo(document.documentElement);
  const hoverLines = observableState<readonly string[]>([]);
  // The scene; hovering a node highlights it (the hover label stands in for the 3D highlight here).
  const scene = mockScene({
    preset: options.preset,
    onHighlight: (label) => hoverLines.set(label ? [label] : []),
    onAction: (text) => status.post(text),
    onMeasure: () => shell.show('selection', 'measure'),
  });
  // Superposing, with two or more structures: the Selection tab's Superpose section.
  registerMockSuperposition(
    shell,
    mockSuperposition(scene, (t) => status.post(t)),
  );
  // Adding measurements: the Selection tab's Measure section.
  registerMockMeasurements(
    shell,
    mockMeasurements(scene, (t) => status.post(t)),
  );
  const hover: HoverLabelModel = { lines: hoverLines.asReadonly };
  const home = mockHome((text, kind) => status.post(text, kind), { remoteStates: options.remoteStates });
  const snapshotsModel = mockSnapshotsModel((text, kind) => status.post(text, kind));
  const stateTree = mockStateTree((text) => status.post(text));
  const settings = mockSettings();
  const snapshotStepper = mockSnapshots((t) => status.post(t));
  const commands = mockCommands(shell, () => scene.roots.value.filter((n) => n.kind === 'structure').length);

  shell.registerTabs([
    {
      id: 'home',
      compactBar: true,
      title: 'Home',
      // The Mol* wordmark over Home (the rail's tooltip stays "Home").
      heading: () => <MolstarWordmark live />,
      icon: 'molstar',
      order: 0,
      place: 'head',
      render: domView(() => <HomeTab home={home} onOpenCommandMenu={() => commandMenuOpen.set(true)} />),
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
          scene={scene}
          presets={mockPresets((t) => status.post(t))}
          options={mockComponentOptions()}
          quickStyles={mockQuickStyles((t) => status.post(t))}
          openCommands={(group) => {
            commandMenuGroup.set(group);
            commandMenuOpen.set(true);
          }}
          onAction={(a) =>
            status.post(
              `${a === 'preset' ? 'Apply preset' : a === 'styles' ? 'Styles' : 'Add component'}… opens in the command menu`,
            )
          }
        />
      )),
    },
    {
      id: 'snapshots',
      title: 'Snapshots',
      icon: 'images',
      order: 1.5,
      render: domView(() => <SnapshotsPanel model={snapshotsModel} />),
    },
    {
      id: 'state-tree',
      title: 'State tree',
      icon: 'list-tree',
      order: 2,
      render: domView(() => <StateTreeTab model={stateTree} />),
    },
    {
      // With the agent module: History (the commands, the assistant's turns, and the messages under Messages).
      id: 'log',
      title: agent ? 'History' : 'Log',
      icon: 'scroll-text',
      order: 2.5,
      // Not a working panel: with Settings at the rail's foot.
      place: 'foot',
      render: domView(() => (agent ? <HistoryTab agent={agent} status={status} /> : <LogView status={status} />)),
    },
    ...(agent
      ? [
          {
            id: 'assistant',
            // The rail's tooltip says what it is for; the panel is headed Assistant.
            title: 'Set up assistant',
            heading: 'Assistant',
            icon: 'sparkles' as const,
            order: 2.7,
            place: 'foot' as const,
            render: domView(() => <AssistantTab agent={agent} />),
          },
        ]
      : []),
    {
      id: 'settings',
      title: 'Settings',
      icon: 'settings',
      order: 3,
      place: 'foot',
      render: domView(() => <SettingsTab model={settings} />),
    },
  ]);
  if (options.loaded) {
    shell.registerTabs([
      {
        id: 'sequence',
        title: 'Sequence',
        icon: 'sequence',
        dock: 'top',
        render: domView(({ ctx }) => (
          <SequenceStrip
            chain={{ structure: '1TQN', name: 'Chain A', info: 'Cytochrome P450 3A4', sequence: Sequence, start: 1 }}
            selectors={sequenceSelectors((t) => status.post(t))}
            onCollapse={() => ctx.shell.toggle(ctx.dock ?? 'top')}
            onHover={(n) =>
              hoverLines.set(n === undefined ? [] : [`${ThreeLetter[Sequence[n - 1]] ?? '?'} ${n}`, 'Chain A · 1TQN'])
            }
            onToggle={(n, on) => status.post(`${on ? 'Selected' : 'Deselected'} residue ${n}`)}
          />
        )),
      },
    ]);
    shell.show('structure');
  }

  // Remote states (a state server): a card next to Export, only when the spec turns them on.
  if (home.remoteStates) {
    const remote = home.remoteStates;
    shell.registerOverlays([
      {
        id: 'remote-states',
        title: 'Remote states',
        render: domView(() => <RemoteStatesPanel model={remote} />),
      },
    ]);
  }

  // The canvas stack's camera and rendering tools act on what is shown: nothing loaded, nothing to show them for.
  const hasData = of(options.loaded);

  shell.registerToolbar([
    // The UI theme, by Settings at the rail's foot: a click goes to the next one; the icon is the current one.
    {
      id: 'theme',
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
      icon: 'command',
      title: 'Command menu',
      caption: formatChord({ key: 'k', mod: true }),
      run: () => commandMenuOpen.set(true),
    },
    {
      id: 'reset',
      place: 'canvas',
      group: 'canvas',
      visible: hasData,
      icon: 'focus',
      title: 'Reset camera',
      run: () => status.post('Camera reset'),
      menu: [
        {
          id: 'reset-zoom',
          label: 'Reset zoom',
          icon: 'zoom-out',
          hint: 'Fit everything, keep the orientation',
          run: () => status.post('Zoom reset'),
        },
        {
          id: 'orient-axes',
          label: 'Orient axes',
          icon: 'axis-3d',
          hint: 'Principal axes along the screen',
          run: () => status.post('Oriented along the principal axes'),
        },
        {
          id: 'reset-axes',
          label: 'Reset axes',
          icon: 'rotate-ccw',
          hint: 'Back to the default orientation',
          run: () => status.post('Axes reset'),
        },
      ],
    },
    // The scene as a whole: how it plays, how the canvas renders it.
    { id: 'animation', group: 'scene', icon: 'clapperboard', title: 'Animation', overlay: 'animation' },
    {
      id: 'canvas-settings',
      group: 'scene',
      icon: 'monitor-cog',
      title: 'Canvas settings',
      overlay: 'canvas-settings',
    },
    { id: 'export', group: 'export', icon: 'download', title: 'Export…', overlay: 'export' },
    ...(home.remoteStates
      ? [
          {
            id: 'remote-states',
            group: 'export',
            icon: 'cloud-upload' as const,
            title: 'Remote states',
            overlay: 'remote-states',
          },
        ]
      : []),
    // Rendering modes under one button: a click opens the menu (no action of its own); shown pressed while one is on.
    {
      id: 'rendering',
      place: 'canvas',
      group: 'canvas',
      visible: hasData,
      icon: 'lightbulb',
      title: 'Illumination & XR',
      pressed: illumination.asReadonly,
      menu: [
        {
          id: 'illumination',
          label: 'Illumination',
          icon: 'lightbulb',
          hint: 'Path-traced global illumination',
          pressed: illumination.asReadonly,
          run: () => illumination.set(!illumination.value),
        },
        {
          id: 'xr',
          label: 'Enter XR',
          icon: 'glasses',
          hint: 'View on a headset',
          run: () => status.post('XR is not available in this browser.', 'warning'),
        },
      ],
    },
    // Expand by context (FullscreenModel.mode): on a page, Expand fills it and Fullscreen is beside it; in an iframe
    // that allows it, the button is Fullscreen (expanding would only fill the frame); in one that does not, Expand fills
    // the frame and Fullscreen shows why it is not possible.
    fullscreen.mode === 'fullscreen'
      ? {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'fullscreen' as const,
          title: 'Fullscreen',
          pressed: fullscreen.active,
          run: () => fullscreen.toggle(),
        }
      : {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'maximize-2' as const,
          title: 'Expand',
          run: () => status.post('Expanded'),
          menu: [
            {
              id: 'fullscreen',
              label: 'Fullscreen',
              icon: 'fullscreen' as const,
              hint: 'The whole screen, not just the page',
              disabled: fullscreen.unavailableReason,
              pressed: fullscreen.active,
              run: () => fullscreen.toggle(),
            },
          ],
        },
    // Panels over the canvas: brought back from the stack when hidden (hidden, they take no room at all).
    {
      id: 'sequence',
      place: 'canvas',
      group: 'view',
      icon: 'sequence',
      title: 'Sequence',
      // Only while there is a sequence and its card is hidden (shown, it carries its own collapse button).
      visible: combineLatest([shell.dock('top').tabs, shell.dock('top').active]).pipe(
        map(([tabs, active]) => tabs.length > 0 && !active),
      ),
      run: () => shell.toggle('top'),
    },
  ]);

  // One Export button; built-in and extension exports are sections of its panel.
  mockExport(shell, status);
  registerMockCanvasSettings(shell);
  registerMockAnimation(shell, (text) => status.post(text));
  const selection = mockSelection(shell, (text) => status.post(text));

  if (options.loaded) {
    status.post('Loaded 1TQN (3,624 atoms) in 0.4 s', 'success');
    // Two tasks at once: they stack.
    const drop = (id: string) => () => status.setTasks(status.tasks.value.filter((t) => t.id !== id));
    status.setTasks([
      { id: 'density', label: 'Downloading density', progress: 0.7, abort: drop('density') },
      { id: 'surface', label: 'Computing surface', progress: 0.4, abort: drop('surface') },
    ]);
  }

  // The canvas mode lives on the canvas, not in the rail: what a click does is visible before the first click.
  const modeProps: ModeSwitchProps = {
    mode: canvasMode.asReadonly,
    setMode: (m) => {
      if (m === canvasMode.value) return;
      canvasMode.set(m);
      status.post(m === 'select' ? 'Select: click to select, Shift+click to extend' : 'Interact: click to focus');
    },
    granularity: granularity.asReadonly,
    granularities: Granularities,
    setGranularity: (g) => granularity.set(g),
  };
  const modeSwitch = <ModeSwitch {...modeProps} />;
  const pickLevel = <PickLevel {...modeProps} />;

  return {
    shell,
    status,
    hover,
    hoverLines,
    modeSwitch,
    pickLevel,
    selection,
    commandMenuOpen,
    commandMenuGroup,
    commands,
    prefixes: mockPrefixes({ manager: commands, shell, scene, snapshots: snapshotStepper }),
    snapshots: snapshotStepper,
    trajectory: mockTrajectory(),
    home,
    agent,
  };
}

function ShellStory(props: {
  loaded: boolean;
  preset?: MockScenePreset;
  background?: string;
  openExportAt?: string;
  /** The right panel pinned beside the canvas. */
  pinned?: boolean;
  /** The spec turns remote states (a state server) on; a number fills the server with that many. */
  remoteStates?: boolean | number;
  /** Opens with the command menu open. */
  commandMenu?: boolean;
  /** What the loaded data has to play: a session's snapshots, a trajectory's models (plain 1TQN has neither). */
  playback?: { snapshots?: boolean; trajectory?: boolean };
  /** The app includes the agent module (set up with a model, or not yet). */
  ai?: { configured?: boolean };
}) {
  const {
    shell,
    status,
    hover,
    hoverLines,
    modeSwitch,
    pickLevel,
    selection,
    commandMenuOpen,
    commandMenuGroup,
    commands,
    prefixes,
    snapshots,
    trajectory,
    home,
    agent,
  } = mockShell(props);
  // What commands did goes to the status line; ⌘K / ⇧⌘P open the command menu; commands' own shortcuts run them.
  const outcomes = commands.outcomes.subscribe((o) => {
    if (o.kind === 'done' && o.result.message) status.post(o.result.message);
    // The person's commands go to History, one row each.
    if (o.kind === 'done') agent?.record(o.result.message ?? o.command.title, o.command.id);
    else if (o.kind === 'failed' || o.kind === 'blocked') status.post(o.message, 'warning');
  });
  const onKey = (e: KeyboardEvent) => {
    const mod = isMac() ? e.metaKey : e.ctrlKey;
    if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
      e.preventDefault();
      commandMenuOpen.set(true);
    } else if (!commandMenuOpen.value && commands.handleKey(e, isMac())) e.preventDefault();
  };
  window.addEventListener('keydown', onKey);
  onCleanup(() => {
    outcomes.unsubscribe();
    window.removeEventListener('keydown', onKey);
  });
  shell.mount();
  status.mount();
  if (props.openExportAt) shell.openOverlay('export', undefined, props.openExportAt);
  if (props.commandMenu) commandMenuOpen.set(true);
  onCleanup(() => {
    shell.unmount();
    status.unmount();
  });
  return (
    <ShellView
      shell={shell}
      status={status}
      canvas={
        <CanvasPlaceholder
          empty={!props.loaded}
          background={props.background}
          onHover={(i) =>
            hoverLines.set(
              i === undefined ? [] : [`${ThreeLetter[Sequence[i * 4]] ?? 'GLY'} ${i * 4 + 1}`, 'Chain A · 1TQN'],
            )
          }
        />
      }
      overlays={
        props.loaded ? (
          <>
            <HoverLabel model={hover} area={() => document.querySelector('main')} />
          </>
        ) : null
      }
      playback={
        props.loaded && (props.playback?.snapshots || props.playback?.trajectory) ? (
          <>
            {props.playback.snapshots ? <SnapshotCard model={snapshots} /> : null}
            {props.playback.trajectory ? <TrajectoryBar model={trajectory} /> : null}
          </>
        ) : null
      }
      canvasColor={props.background ?? '#ffffff'}
      canvasTheme={themeForColor(parseColor(props.background ?? '#ffffff') ?? 0xffffff)}
      context={props.loaded ? () => <SelectionBar selection={selection} /> : undefined}
      commandMenu={
        <Show when={commandMenuOpen}>
          {() => (
            <CommandMenu
              manager={commands}
              prefixes={prefixes}
              group={commandMenuGroup.value}
              assistant={agent ? mockAssistant(agent, () => shell.show('assistant')) : undefined}
              onClose={() => {
                commandMenuOpen.set(false);
                commandMenuGroup.set(undefined);
              }}
            />
          )}
        </Show>
      }
      agent={agent ? <AgentTurnCard agent={agent} /> : undefined}
      canvasLead={props.loaded ? modeSwitch : null}
      canvasAside={props.loaded ? pickLevel : null}
      onOpenLog={() => shell.show('log')}
      empty={of(!props.loaded)}
      onDropFiles={(files) => home.openFiles(files)}
      dropFormats={home.formats}
    />
  );
}

/** The command menu open (⌘K / Ctrl+K anywhere, or the rail's ⌘): everything by name. */
export const CommandMenuOpen: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} commandMenu={true} />}</div>),
};

/** Nothing loaded: Home is the open tab. */
export const Empty: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={false} />}</div>),
};

/** A structure loaded: the Structure tab, a status message, a running task and the selection chip. */
export const Loaded: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} />}</div>),
};

/** Several structures: four CYP3A4 structures superposed on 1TQN, collapsed but the first. */
/**
 * A session with snapshots and a multi-model trajectory: playback appears at the top left only because the data has
 * something to play (plain 1TQN has neither).
 */
export const SessionWithPlayback: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} playback={{ snapshots: true, trajectory: true }} />
      </div>
    )),
};

export const ManyStructures: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} preset="many" />}</div>),
};

/** An EMDB entry: the map is a top-level volume next to the fitted model (EMD-2984 and 5A1A). */
export const EmdbMap: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} preset="em" />}</div>),
};

/**
 * The chrome over the canvas follows the canvas background automatically (light glass on a light background, dark on
 * a dark one), whatever the UI theme. Pick a background color in the controls.
 */
export const CanvasBackground: StoryObj<{ background: string }> = {
  args: { background: '#14181f' },
  argTypes: { background: { control: 'color' } },
  render: (args) =>
    story(() => <div style="height: 100vh">{<ShellStory loaded={true} background={args.background} />}</div>),
};

/** A command opened the Export panel at a section ("Export animation…" in the command menu). */
export const ExportAtSection: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} openExportAt="animation" />
      </div>
    )),
};

/**
 * The panel pinned (the pin in its header): it sits beside the canvas, which makes room, instead of floating over it;
 * for workflows that keep a panel open the whole time.
 */
export const Pinned: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} pinned={true} />}</div>),
};

/** Remote states turned on in the spec: a Remote states card next to Export, and the server's states on Home. */
export const RemoteStates: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} remoteStates={true} />}</div>),
};

/** A busy state server (50 states): the lists scroll on their own, on Home and in the Remote states card. */
export const ManyRemoteStates: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={false} remoteStates={50} />}</div>),
};

/** A narrow host: the panel becomes a bottom sheet with its tabs along the bottom. */
export const Compact: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh; display: grid; place-items: center; background: var(--msp-color-border)">
        <div style="width: 390px; height: 760px; border-radius: 16px; overflow: hidden; box-shadow: var(--msp-shadow-overlay)">
          <ShellStory loaded={true} />
        </div>
      </div>
    )),
};

/**
 * Expand by context, in real frames of the Loaded story (hover the canvas's Expand button, at the top right):
 * - an iframe that allows fullscreen (`allow="fullscreen"`): the button is Fullscreen, since expanding would only fill
 *   the frame; Escape leaves it;
 * - an iframe that does not: Expand fills the frame, and Fullscreen beside it says why it is not possible;
 * - a page (open it below): Expand fills the page, Fullscreen beside it.
 */
export const ExpandByContext: StoryObj = {
  parameters: { layout: 'fullscreen' },
  render: () =>
    story(() => {
      const src = 'iframe.html?id=shell-shell--loaded&viewMode=story';
      const caption = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);';
      const frame =
        'width: 100%; height: 100%; border: 0; border-radius: var(--msp-radius-lg); box-shadow: inset 0 0 0 1px var(--msp-glass-edge);';
      return (
        <div style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-3); height: 100vh; padding: var(--msp-space-4); box-sizing: border-box; background: var(--msp-color-bg);">
          <div style={`${caption} display: flex; gap: var(--msp-space-4);`}>
            <span>Hover Expand at the top right of each canvas.</span>
            <a href={src} target="_blank" rel="noopener" style="color: var(--msp-color-accent)">
              Open as a page (Expand, Fullscreen beside it)
            </a>
          </div>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: var(--msp-space-4); min-height: 0;">
            <figure style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-2); margin: 0; min-height: 0;">
              <figcaption style={caption}>In an iframe that allows fullscreen: the button is Fullscreen</figcaption>
              <iframe title="Allows fullscreen" src={src} allow="fullscreen" style={frame} />
            </figure>
            <figure style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-2); margin: 0; min-height: 0;">
              <figcaption style={caption}>In an iframe that does not: Expand, Fullscreen not possible</figcaption>
              {/* Same-origin frames may go fullscreen by default ('self'): forbidden explicitly, as a page on another origin is. */}
              <iframe title="Does not allow fullscreen" src={src} allow="fullscreen 'none'" style={frame} />
            </figure>
          </div>
        </div>
      );
    }),
};

/**
 * The agent module included and set up (a mock agent stands in for any-harness): ask from the command menu (⌘K, type
 * a request, Ask the assistant; “color it red” opens the filled sentence, “show the active site”, “load 1w0e and
 * superpose it”, “try something broken” run turns), watch the turn card, Allow or Deny a download, Stop, Undo turn;
 * History at the rail's foot; the Assistant tab with the model and keys.
 */
export const Assistant: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} ai={{ configured: true }} />}</div>),
};

/** The agent module included but no model yet: the Assistant tab asks for a key; the Ask row opens it. */
export const AssistantSetup: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} ai={{}} />}</div>),
};
