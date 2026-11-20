/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A unit of a complex read as a system by each consumer (EC-35): a section on the unit follows its sequence's section
// through the unit's step (EC-26), a hover on either lights the other, the selection of a unit's positions by the
// chain form is what an export of the same spec resolves and draws, and the export draws the unit's letters.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { locationOfResult, locationsEqual, resolveAddressIn } from '@molstar/seqx-core/query';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerConfig, viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type Hit, HoverModel } from '../models/hover.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { partsOf } from './fixture.js';
import { type Frames, installFrames } from './frames.js';

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN';
const Protein = '/data/p';
const Unit = '/data/c/A';
const ProteinSection = '/views/main/protein';
const UnitSection = '/views/main/unit';

/**
 * A protein of 30 residues and a complex whose unit A is its chain over [10, 30); the unit's section follows (or is
 * alone, `unitOnly`).
 */
function spec(unitOnly = false): Node {
  const b = createSeqSpec({ title: 'units' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  let unit = '';
  b.data.complex('c', { name: 'a complex' }, (c) => {
    unit = c.unit('A', { of: p, range: [10, 30], chain: 'A' });
    c.unit('B', { of: p, range: [10, 30], chain: 'B' });
  });
  b.views.view('main', (v) => {
    if (unitOnly) {
      v.section('unit', { of: unit as never })
        .track('letters')
        .layer('l', { type: 'letters' });
      return;
    }
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    v.section('unit', { of: unit as never, follow: protein })
      .track('letters')
      .layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('a unit read as a system (EC-35)', () => {
  it('follow: the unit’s section follows its sequence’s through the unit’s step, clipped to its range (EC-26)', () => {
    const view = new ViewModel(new SpecModel(spec()), viewerConfig());
    view.setWidth(ProteinSection, 20);
    view.setWidth(UnitSection, 20);
    view.setRange(ProteinSection, 12, 20);
    expect([view.get(UnitSection)?.start, view.get(UnitSection)?.end]).toEqual([2, 10]);
    view.setRange(ProteinSection, 12.5, 20);
    expect(view.get(UnitSection)?.start).toBeCloseTo(2.5, 10);
  });

  it('hover projection: a residue of the unit lights its sequence’s, and back; a copy is never reached (EC-25)', () => {
    const spec0 = new SpecModel(spec());
    const view = new ViewModel(spec0, viewerConfig());
    const hover = new HoverModel(spec0, view, viewerConfig());
    const onUnit: Hit = { layer: `${UnitSection}/letters/l`, system: Unit, position: 3 };
    hover.set(onUnit, 'p');
    expect(partsOf(hover.projected(ProteinSection))).toEqual([[Protein, [13, 14]]]);
    expect(partsOf(hover.projected(UnitSection))).toEqual([[Unit, [3, 4]]]);
    hover.set({ layer: `${ProteinSection}/letters/l`, system: Protein, position: 5 }, 'p');
    // a residue before the unit's range lights nothing on it
    expect(partsOf(hover.projected(UnitSection))).toEqual([]);
  });
});

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('navigator.form'), 'none'],
    ],
  });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

describe('a unit’s selection and export: two consumers of one address', () => {
  it('the chain form selects the unit’s positions; the viewer and an export resolve and draw it alike', async () => {
    const v = await make();
    const address = '/data/c/(chain :id A)/(range :from 2 :to 5)';
    v.addresses.select(address);
    frames!.flush();
    const root = spec();
    const model = resolveSpec(root, { runtimes: [SeqRuntime] });
    const r = resolveAddressIn(root, model, address);
    const exported = locationOfResult(root, model, r.result!).location;
    expect(partsOf(exported)).toEqual([[Unit, [1, 5]]]);
    expect(locationsEqual(exported, v.selection.value.value.shown)).toBe(true);

    const mine = await v.renderSvg();
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const theirs = await ctx.renderSvg(root, {
      size: { width: v.width, height: v.height },
      ...(v.viewport ? { viewport: v.viewport } : {}),
      emphasis: { selection: exported },
    });
    ctx.retire();
    expect(mine.svg).toBe(theirs.svg);
    // the selection is drawn on the unit's section and, carried up its step, on the protein's
    const selection = recordDrawLists(theirs.lists, ['selection']);
    expect(selection).toBe(recordDrawLists(v.frame!.lists, ['selection']));
    expect(selection.split('\n').filter((l) => l.includes('rect')).length).toBeGreaterThanOrEqual(2);
  });

  it('an export of the unit alone draws its own letters: its sequence’s over its range', async () => {
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const r = await ctx.renderSvg(spec(true), { size: { width: 800, height: 120 } });
    ctx.retire();
    expect(r.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    const glyphs = recordDrawLists(r.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('glyphs'))
      .join(' ');
    // W is the protein's 23rd residue, the unit's 13th; M (its first) is before the unit's range
    expect(glyphs).toContain('W ');
    expect(glyphs).not.toContain('M ');
  });
});
