/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';
import { MmcifRuntime } from '@molstar/seqx-core/runtimes/mmcif';
import * as addressing from '@molstar/seqx-examples/services/addressing';
import * as composed from '@molstar/seqx-examples/services/composed';
import { serviceUrl } from '@molstar/seqx-examples/services/endpoints';
import * as predictions from '@molstar/seqx-examples/services/predictions';
import { parseSifts, type SiftsMapping } from '@molstar/seqx-examples/services/sifts';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { type FakeLink, linkFake } from '../features/structure-link/compose.js';
import { colorRunner } from '../features/structure-link/compose.js';
import type { ChainTable } from './fake-structure/chain-table.js';
import { ColorByEntry, colorByBehavior } from './host-adapter.js';
import { type PortClock, type ServiceMode, servicesPort, type ServiceRead } from './services-port.js';

// An application composing services, DOM-free (the stories and their headless twins call this): the viewer, the services'
// assets port in the mode the application is in, SIFTS read by the application into maps, and the fake structure viewer
// linked once the structure and SIFTS are both here. A scenario is one view of it (the composition, the predictions, the
// addressing); the mode is a switch the application owns: it starts again, so the spec is the same and only the answers
// differ (recorded files, or the services through the asset policy's hosts).

/** One view of the application: what it builds from SIFTS's answer, and the structure it links. */
export interface Scenario {
  readonly id: 'composed' | 'predictions' | 'addressing';
  readonly build: (options: { readonly sifts?: readonly SiftsMapping[] }) => Node;
  /** The PDB entry whose SIFTS mappings the spec waits for. */
  readonly entry: string;
  /** The fake's chain table of the entry (relative to the examples). */
  readonly tableFile: string;
  /** The host's header entry "Color the structure by this track" on every track. */
  readonly colorBy?: boolean;
}

export const Scenarios: Readonly<Record<Scenario['id'], Scenario>> = {
  composed: {
    id: 'composed',
    build: composed.build,
    entry: composed.Entry,
    tableFile: 'fixtures/p53/2OCJ-chains.json',
  },
  predictions: {
    id: 'predictions',
    build: predictions.build,
    entry: composed.Entry,
    tableFile: 'fixtures/p53/2OCJ-chains.json',
    colorBy: true,
  },
  addressing: {
    id: 'addressing',
    build: addressing.build,
    entry: composed.Entry,
    tableFile: 'fixtures/p53/2OCJ-chains.json',
  },
};

/** Where the application is: the structure and the data are on their way, SIFTS is, it is whole, or a service failed. */
export type Phase = 'loading' | 'sifts' | 'ready' | 'failed';

export interface ServicesAppOptions {
  readonly scenario: Scenario;
  /** The chain table of the scenario's entry. */
  readonly table: ChainTable;
  readonly fixture: (file: string) => Promise<Uint8Array>;
  readonly mode?: ServiceMode;
  /** A service's latency (recorded answers wait this long, on the clock). */
  readonly delayMs?: number;
  readonly clock?: PortClock;
  /** Where the canvas mounts: none in a headless twin. */
  readonly target?: HTMLElement;
  /** The theme the structure is colored for (the host's). */
  readonly scheme?: () => 'light' | 'dark';
  /** The live reader (a twin's refuses the network); the policy's by default. */
  readonly live?: Parameters<typeof servicesPort>[0]['live'];
  /** Extra behaviors of the viewer (a host's choices). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
}

const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));

export class ServicesApp {
  readonly mode: State<ServiceMode>;
  readonly phase = new State<Phase>('loading');
  /** What the services were asked, newest last (the log the story prints). */
  readonly reads = new State<readonly ServiceRead[]>([]);
  /** What went wrong or what is awaited, a sentence. */
  readonly message = new State<string>('');
  readonly sifts = new State<readonly SiftsMapping[]>([]);
  /** Moves whenever the viewer or the link is another one: what a view of the application binds again. */
  readonly objects = new State(0);
  viewer: SeqViewer | undefined;
  link: FakeLink | undefined;
  private _port: ReturnType<typeof servicesPort> | undefined;
  private _run = 0;
  private _gone = false;

  constructor(private readonly options: ServicesAppOptions) {
    this.mode = new State<ServiceMode>(options.mode ?? 'fixtures');
  }

  /** Starts (or starts again, in the mode now): the viewer without SIFTS, then with it, the fake linked. */
  async start(): Promise<void> {
    const run = ++this._run;
    const o = this.options;
    this._teardown();
    this.reads.set([]);
    this.sifts.set([]);
    this.phase.set('loading');
    this.message.set('');
    const port = (this._port = servicesPort({
      mode: () => this.mode.value,
      fixture: o.fixture,
      ...(o.delayMs !== undefined ? { delayMs: o.delayMs } : {}),
      ...(o.clock ? { clock: o.clock } : {}),
      ...(o.live ? { live: o.live } : {}),
      onRead: (r) => this.reads.set([...this.reads.value, r]),
    }));
    const base = DefaultSeqViewerSpec();
    const colorBy = !!o.scenario.colorBy;
    const viewer: SeqViewer = (this.viewer = createSeqViewer(
      {
        ...base,
        ui: o.target ? { ...base.ui, toolbar: mountToolbar } : undefined,
        registry: [
          ...(base.registry ?? []),
          { kind: 'format-runtime', runtime: MmcifRuntime },
          ...(colorBy ? [ColorByEntry] : []),
        ],
        behaviors: [
          ...base.behaviors,
          ...(o.behaviors ?? []),
          ...(colorBy
            ? [
                colorByBehavior((layer) =>
                  colorRunner(
                    viewer,
                    () => this.link,
                    () => o.scheme?.() ?? 'light',
                  )(layer),
                ),
              ]
            : []),
        ],
        ports: { ...base.ports, assets: () => port },
      },
      o.target,
    ));
    this.objects.set(this.objects.value + 1);
    if (!o.target) viewer.resize(1200);
    const stale = () => this._gone || run !== this._run;
    try {
      await viewer.load(o.scenario.build({}));
      await viewer.settle();
      if (stale()) return;
      this.phase.set('sifts');
      this.message.set(`Waiting for SIFTS's mappings of ${o.scenario.entry}`);
      const text = await port.reader.readText(serviceUrl.sifts(o.scenario.entry), new AbortController().signal);
      if (stale()) return;
      const sifts = parseSifts(text);
      this.sifts.set(sifts);
      await viewer.load(o.scenario.build({ sifts }));
      await viewer.settle();
      if (stale()) return;
      try {
        this.link = linkFake(viewer, o.table, composed.ServiceAddresses.complex);
        this.objects.set(this.objects.value + 1);
        this.phase.set('ready');
        this.message.set('');
      } catch (e) {
        // the structure did not arrive (its file failed): the data is drawn, the structure is not linked
        this.phase.set('failed');
        this.message.set(`The structure is not linked: ${messageOf(e)}`);
      }
    } catch (e) {
      if (stale()) return;
      this.phase.set('failed');
      this.message.set(messageOf(e));
    }
  }

  /** Switches the mode and starts again. */
  async setMode(mode: ServiceMode): Promise<void> {
    if (mode === this.mode.value && this.phase.value !== 'failed') return;
    this.mode.set(mode);
    await this.start();
  }

  private _teardown() {
    this.link?.dispose();
    this.link = undefined;
    this.viewer?.retire();
    this.viewer = undefined;
    this._port?.dispose?.();
    this._port = undefined;
  }

  dispose() {
    this._gone = true;
    this._run++;
    this._teardown();
  }
}
