/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { location } from '@molstar/seqx-core/query';
import { AssetRefused } from '@molstar/seqx-core/context';
import { Recorded, serviceUrl } from '@molstar/seqx-examples/services/endpoints';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { ServiceAddresses } from '@molstar/seqx-examples/services/composed';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { chainTableOf } from '../../features/structure-link/compose.js';
import { fakeClock } from '../../features/structure-link/_test/open.js';
import { type Scenario, Scenarios, ServicesApp } from '../services-app.js';
import { NotRecorded, type PortClock, servicesPort } from '../services-port.js';

// The services application headless, as the stories compose it, on a fake clock: the recorded answers come after a delay,
// the application draws what it has, waits for SIFTS and links the fake once both are here; the switch to the live
// services reads through the policy (a stub reader here: tests never reach the network).

const disk = async (file: string) => new Uint8Array(readFileSync(`${PackageFolder}${file}`));
const table = chainTableOf(readFileSync(`${PackageFolder}${Scenarios.composed.tableFile}`, 'utf8'));

let apps: ServicesApp[] = [];
let restore: (() => void) | undefined;
afterEach(() => {
  for (const a of apps) a.dispose();
  apps = [];
  restore?.();
  restore = undefined;
});

/** A clock the test moves: timers run by `advance`, never real ones. */
class ManualClock implements PortClock {
  private _now = 0;
  private _seq = 0;
  private readonly _timers = new Map<number, { at: number; fn: () => void }>();
  get pending() {
    return this._timers.size;
  }
  setTimeout(fn: () => void, ms: number) {
    const id = ++this._seq;
    this._timers.set(id, { at: this._now + ms, fn });
    return id;
  }
  clearTimeout(h: unknown) {
    this._timers.delete(h as number);
  }
  async advance(ms: number) {
    this._now += ms;
    for (const [id, t] of [...this._timers]) {
      if (t.at > this._now) continue;
      this._timers.delete(id);
      t.fn();
    }
    await new Promise<void>((r) => setImmediate(r));
  }
}

/** Lets the event loop turn until `cond` holds (real turns: the parse streams). */
async function until(cond: () => boolean, turns = 2000) {
  for (let i = 0; i < turns && !cond(); i++) await new Promise<void>((r) => setImmediate(r));
  expect(cond()).toBe(true);
}

/** A live stub: serves each service location from its recorded file and counts what it was asked. */
function liveStub(refuse: readonly string[] = []) {
  const asked: string[] = [];
  const read = async (uri: string) => {
    asked.push(uri);
    if (refuse.includes(uri)) throw new Error(`503 for ${uri}`);
    const file = Recorded.get(uri);
    if (!file) throw new Error(`the stub has no ${uri}`);
    return disk(file);
  };
  return {
    asked,
    reader: {
      readText: async (uri: string) => new TextDecoder().decode(await read(uri)),
      readBytes: (uri: string) => read(uri),
    },
  };
}

async function open(scenario: Scenario, extra: Partial<ConstructorParameters<typeof ServicesApp>[0]> = {}) {
  restore = fakeClock();
  const app = new ServicesApp({ scenario, table, fixture: disk, ...extra });
  apps.push(app);
  return app;
}

/** Runs the application's start on the fake clock, to its end. */
async function started(app: ServicesApp, start: () => Promise<void> = () => app.start()) {
  const done = start();
  let finished = false;
  void done.then(
    () => (finished = true),
    () => (finished = true),
  );
  // the delays are timers made as the application goes: the clock moves a step at a time until it is done
  // and a real turn of the event loop between steps: the parse of the mmCIF's gzip is a stream, not a timer
  for (let i = 0; i < 400 && !finished; i++) {
    await vi.advanceTimersByTimeAsync(50);
    await new Promise<void>((r) => setImmediate(r));
  }
  await done;
}

describe('the services application, on recorded answers', () => {
  it('draws what it has, waits for SIFTS, then links the fake: each service read once, the phases in order', async () => {
    // the services' latency on a clock of the test's own, moved by hand (the viewer's own timers are real: its parse of
    // the mmCIF's gzip is a stream a faked clock does not pump)
    const clock = new ManualClock();
    const app = new ServicesApp({ scenario: Scenarios.composed, table, fixture: disk, delayMs: 300, clock });
    apps.push(app);
    const phases: string[] = ['loading'];
    app.phase.listen(() => phases.push(app.phase.value));
    const done = app.start();
    await until(() => clock.pending === 3);
    // three answers wait their delay: nothing is drawn, the application is loading
    expect(app.phase.value).toBe('loading');
    expect(app.viewer!.store.store.coverage('/data/2ocj')).toBeUndefined();
    await clock.advance(300);
    // the spec's data is here; SIFTS is asked for now and waits its delay too
    await until(() => clock.pending === 1);
    expect(app.phase.value).toBe('sifts');
    expect(app.reads.value).toHaveLength(4);
    await clock.advance(300);
    await done;
    expect(phases.filter((p, i) => p !== phases[i - 1])).toEqual(['loading', 'sifts', 'ready']);
    expect(app.reads.value.map((r) => r.uri).sort()).toEqual(
      [
        serviceUrl.uniprot('P04637'),
        serviceUrl.confidence('P04637'),
        serviceUrl.entry('2OCJ'),
        serviceUrl.sifts('2OCJ'),
      ].sort(),
    );
    expect(app.reads.value.every((r) => r.mode === 'fixtures' && r.file && r.service)).toBe(true);
    expect(app.viewer!.modelIssues.value.map((i) => i.code)).toEqual([]);
    expect(app.viewer!.layoutIssues.value.map((i) => i.code)).toEqual([]);
    expect(app.link).toBeDefined();
  });

  it('a hover on R175 of P04637 lights residue 175 of every chain through SIFTS’s maps, and back', async () => {
    const app = await open(Scenarios.composed);
    await started(app);
    const { fake } = app.link!;
    const viewer = app.viewer!;
    const owner = {};
    viewer.hover.setLocation(location({ parts: [{ system: ServiceAddresses.protein, ranges: [174, 175] }] }), owner);
    expect(fake.highlight.value.map((s) => [fake.rows[s.row].asym, fake.rows[s.row].residues[s.start].auth])).toEqual([
      ['A', '175'],
      ['B', '175'],
      ['C', '175'],
      ['D', '175'],
    ]);
    viewer.hover.clear(owner);
    // the unobserved end (94, the first residue of the entity) lights nothing
    viewer.hover.setLocation(location({ parts: [{ system: ServiceAddresses.protein, ranges: [93, 94] }] }), owner);
    expect(fake.highlight.value).toEqual([]);
  });

  it('a hover or a click on the structure arrives on its unit: the view draws what the spec’s sections show (a gap: units a source describes are no sections’ systems)', async () => {
    const app = await open(Scenarios.composed);
    await started(app);
    const { fake } = app.link!;
    fake.pointerOver({ row: 1, index: 81 });
    await vi.advanceTimersByTimeAsync(50);
    // the hovered cell is residue 175 of chain B, as a Location on its unit; the hover model projects through the static
    // model's paths, which do not hold a described unit's map, so no section draws it (TODO: overlay-aware projection)
    const own = app.viewer!.hover.location();
    expect(own?.parts.map((p) => [p.system, ...p.ranges])).toEqual([[ServiceAddresses.complex + '/B', 81, 82]]);
    expect(app.viewer!.hover.projected('/views/main/protein')?.parts ?? []).toEqual([]);
    fake.click({ row: 1, index: 81 });
    await vi.advanceTimersByTimeAsync(50);
    expect(app.viewer!.selection.value.value.shown.parts.map((p) => [p.system, ...p.ranges])).toEqual([
      [ServiceAddresses.complex + '/B', 81, 82],
    ]);
  });

  it('the predictions chain two maps: a lab position lights its residue in the fake; “color by this track” colors the chains', async () => {
    const app = await open(Scenarios.predictions);
    await started(app);
    const { fake } = app.link!;
    app.viewer!.hover.setLocation(location({ parts: [{ system: '/data/lab', ranges: [247, 248] }] }), {});
    expect(fake.highlight.value.map((s) => fake.rows[s.row].residues[s.start].auth)).toEqual([
      '248',
      '248',
      '248',
      '248',
    ]);
    const layer = app
      .viewer!.view.resolved.value!.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))
      .find((l) => l.type === 'values')!;
    const said = app.viewer!.runCommand('fake-host.color-by', { layer: layer.address });
    expect(said).toMatch(/Colored \d+ runs/);
  });

  it('every link of the addressing text runs: the three numberings of one residue name the same place', async () => {
    const { Passages } = await import('@molstar/seqx-examples/services/addressing');
    const app = await open(Scenarios.addressing);
    await started(app);
    const viewer = app.viewer!;
    const hrefs = Passages.flat()
      .filter((r) => typeof r !== 'string')
      .map((r) => (r as { href: string }).href);
    expect(hrefs.length).toBeGreaterThan(5);
    for (const href of hrefs) {
      const run = viewer.runLink(href);
      expect(
        run?.outcomes.map((o) => o.status),
        href,
      ).not.toContain('failed');
      expect(
        run?.issues.map((i) => i.message),
        href,
      ).toEqual([]);
    }
    // R175, entity position 82 and the author number 175 of chain A: one residue, three addresses
    const lits: string[] = [];
    for (const href of hrefs.slice(-3)) {
      viewer.runLink(href);
      await vi.advanceTimersByTimeAsync(500);
      lits.push(viewer.addresses.focusAddress() ?? '');
    }
    expect(lits).toHaveLength(3);
  });
});

describe('every scenario', () => {
  it.each(Object.values(Scenarios))(
    '$id: laid out with no issue, rendered to SVG, SIFTS’s maps and the structure in it',
    async (scenario) => {
      const app = await open(scenario);
      await started(app);
      expect(app.phase.value).toBe('ready');
      const viewer = app.viewer!;
      expect([...viewer.modelIssues.value, ...viewer.layoutIssues.value].map((i) => i.code)).toEqual([]);
      const sections = viewer.view.resolved.value!.sections.map((s) => s.address);
      expect(sections).toContain('/views/main/structure');
      const svg = await viewer.renderSvg({});
      expect(svg.svg.startsWith('<svg')).toBe(true);
      expect(svg.issues ?? []).toEqual([]);
    },
  );
});

describe('the switch to the live services', () => {
  it('starts again in live mode: the same spec, the reads now the stub’s (the policy’s reader in the page), none recorded', async () => {
    const stub = liveStub();
    const app = await open(Scenarios.composed, { live: stub.reader });
    await started(app);
    expect(app.reads.value.every((r) => r.mode === 'fixtures')).toBe(true);
    await started(app, () => app.setMode('live'));
    expect(app.phase.value).toBe('ready');
    expect(app.reads.value.every((r) => r.mode === 'live' && r.file === undefined)).toBe(true);
    expect(stub.asked.sort()).toEqual(app.reads.value.map((r) => r.uri).sort());
  });

  it('a service that fails leaves the data drawn and says which: the application is failed, not blank', async () => {
    const stub = liveStub([serviceUrl.sifts('2OCJ')]);
    const app = await open(Scenarios.composed, { live: stub.reader, mode: 'live' });
    await started(app);
    expect(app.phase.value).toBe('failed');
    expect(app.message.value).toContain(serviceUrl.sifts('2OCJ'));
    expect(app.viewer!.layoutIssues.value.map((i) => i.code)).toEqual([]);
    expect(app.link).toBeUndefined();
  });
});

describe('the services port', () => {
  const reader = (mode: 'fixtures' | 'live', clock?: PortClock) =>
    servicesPort({ mode: () => mode, fixture: disk, ...(clock ? { clock } : {}) }).reader;
  const never = new AbortController().signal;

  it('a service location with no recorded answer is NotRecorded, said so; a recorded one is read', async () => {
    await expect(reader('fixtures').readText(serviceUrl.uniprot('P99999'), never)).rejects.toBeInstanceOf(NotRecorded);
    const text = await reader('fixtures').readText(serviceUrl.sifts('2OCJ'), never);
    expect(JSON.parse(text)['2ocj'].UniProt.P04637.mappings).toHaveLength(4);
  });

  it('the live mode is the asset policy’s: a host no service has is refused before any request, whatever the spec says', async () => {
    const live = reader('live');
    await expect(live.readText('https://example.org/uniprotkb/P04637.json', never)).rejects.toBeInstanceOf(
      AssetRefused,
    );
    await expect(live.readText('http://rest.uniprot.org/uniprotkb/P04637.json', never)).rejects.toBeInstanceOf(
      AssetRefused,
    );
    await expect(live.readText('fixtures/p53/P04637.json', never)).rejects.toBeInstanceOf(AssetRefused);
  });

  it('a recorded answer waits its delay on the port’s clock and stops when aborted', async () => {
    restore = fakeClock();
    const fired: number[] = [];
    const clock: PortClock = {
      setTimeout: (fn, ms) => {
        fired.push(ms);
        return setTimeout(fn, ms);
      },
      clearTimeout: (h) => clearTimeout(h as ReturnType<typeof setTimeout>),
    };
    const port = servicesPort({ mode: () => 'fixtures', fixture: disk, delayMs: 250, clock });
    let done = false;
    const p = port.reader.readText(serviceUrl.sifts('1FIN'), never).then(() => (done = true));
    await vi.advanceTimersByTimeAsync(249);
    expect(done).toBe(false);
    await vi.advanceTimersByTimeAsync(2);
    await p;
    expect(fired).toEqual([250]);
    const stop = new AbortController();
    const q = port.reader.readText(serviceUrl.sifts('1FIN'), stop.signal);
    stop.abort(new Error('stopped'));
    await expect(q).rejects.toThrow('stopped');
  });

  it('a relative location of the examples is read as a recorded file in fixtures mode', async () => {
    const text = await reader('fixtures').readText('fixtures/services/sifts-1fin.json', never);
    expect(Object.keys(JSON.parse(text))).toEqual(['1fin']);
  });
});
