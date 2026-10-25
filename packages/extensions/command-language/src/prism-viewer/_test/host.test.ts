/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterAll, afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { CommandLanguageExtension } from '../../prism-viewer.js';
import { prefixRows } from '../commands.js';
import type { ScriptRunValue } from '../operations.js';

// The commands on a headless Prism viewer (_spec §5, §6): what each leaves in the tree, one step per submission,
// all or nothing. No network: the structure is 1CBS from data/examples, given with the script as a file.

const cif = readFileSync(new URL('../../../../../../data/examples/1cbs_updated.cif', import.meta.url));
const file = () => new File([cif], '1cbs.cif');

// Measurement texts are drawn into a font atlas, which Node has no canvas for: a stand-in that draws nothing (as
// Prism's measurement tests).
beforeAll(() => {
  const context = {
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
// The module state this file changed, put back (files share their worker's modules: AGENTS.md).
afterAll(() => setCanvasModule(undefined as never));

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer() {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), CommandLanguageExtension] },
  });
  viewers.push(v);
  return v;
}

async function run(v: PrismViewer, text: string, dialect: 'pymol' | 'vmd' = 'pymol') {
  const r = await v.api.script!.run({ text, dialect, files: [file()] });
  return r;
}

async function ok(v: PrismViewer, text: string, dialect: 'pymol' | 'vmd' = 'pymol') {
  const r = await run(v, text, dialect);
  if (!r.ok) expect.fail(`${text}: ${r.code}: ${r.message}${r.hint ? ` (${r.hint})` : ''}`);
  return r as Extract<typeof r, { ok: true }>;
}

/** Each component of the (only) structure: its label, its representations' types, its atom count. */
function components(v: PrismViewer) {
  const s = v.plugin.managers.structure.hierarchy.current.structures[0];
  return (s?.components ?? []).map((c) => ({
    label: (c.cell.params?.values as { label?: string } | undefined)?.label || c.cell.obj?.label || '',
    types: c.representations.map((r) => (r.cell.transform.params as { type: { name: string } }).type.name),
    atoms: c.cell.obj?.data.elementCount ?? 0,
  }));
}

const typesShown = (v: PrismViewer) => [...new Set(components(v).flatMap((c) => c.types))].sort();

describe('PyMOL commands on the viewer', () => {
  it('loads a file given with the script as a named object', async () => {
    const v = await viewer();
    const r = await ok(v, 'load 1cbs.cif');
    expect(r.did).toContain('as 1cbs');
    const names = await v.api.script!.names({});
    expect(names.ok && names.did).toContain('0 1cbs');
  });

  it('shows and hides as PyMOL’s per-atom flags: one component per type, presets edited in place', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const before = components(v);
    expect(before.some((c) => c.types.includes('cartoon'))).toBe(true);
    await ok(v, 'show sticks, resn REA');
    // The ligand already shows as ball-and-stick (the preset): nothing is drawn twice.
    expect(components(v).filter((c) => c.label === 'sticks')).toEqual([]);
    await ok(v, 'show sticks, resi 10-20');
    const sticks = components(v).find((c) => c.label === 'sticks')!;
    expect(sticks.types).toEqual(['ball-and-stick']);
    await ok(v, 'hide cartoon, resi 1-50');
    const cartoon = components(v).find((c) => c.types.includes('cartoon'))!;
    const was = before.find((c) => c.types.includes('cartoon'))!;
    // The preset's polymer component, in place: the same name, fewer atoms.
    expect(cartoon.label).toBe(was.label);
    expect(cartoon.atoms).toBeLessThan(was.atoms);
  });

  it('as: only that representation for the selection', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'as spheres');
    expect(typesShown(v)).toEqual(['spacefill']);
    await ok(v, 'hide everything');
    expect(components(v)).toEqual([]);
  });

  it('shows a cartoon right after a load: on the polymer only, nothing for the waters and the ligand', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'show cartoon');
    await ok(v, 'as cartoon');
    expect(typesShown(v)).toEqual(['cartoon']);
    // Only polymer atoms take it: as many as `count_atoms polymer` reads.
    const cartoon = components(v).filter((c) => c.types.includes('cartoon'));
    const polymer = (await ok(v, 'count_atoms polymer')).value as { outputs: { text: string }[] };
    expect(String(cartoon.reduce((n, c) => n + c.atoms, 0))).toBe(polymer.outputs[0].text.match(/\d+/)?.[0]);
    await ok(v, 'show cartoon, resn REA');
  });

  it('is one step per submission: Undo takes it back whole', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const before = components(v);
    await ok(v, 'hide everything; show cartoon; show sticks, resn REA; color red, resn REA');
    expect(components(v)).not.toEqual(before);
    await v.api.history.undo({});
    expect(components(v)).toEqual(before);
  });

  it('runs nothing when a line does not read, and rolls back a failure', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const before = components(v);
    const bad = await run(v, 'as spheres\ncolor notacolor');
    expect(bad.ok).toBe(false);
    expect(!bad.ok && bad.message).toContain('line 2');
    expect(components(v)).toEqual(before);
    // A runtime failure (a distance between many atoms) after a change: the change is rolled back.
    const failed = await run(v, 'as spheres; distance resn REA, resn REA');
    expect(failed.ok).toBe(false);
    expect(!failed.ok && failed.code).toBe('failed');
    expect(components(v)).toEqual(before);
  });

  it('names selections, and reads counts', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'select lig, resn REA, enable=0');
    const r = await ok(v, 'count_atoms lig and not hydro');
    expect((r.value as ScriptRunValue).outputs.map((o) => o.text)).toContain('22');
    // Undo of the naming step forgets the name with it (custom state).
    await v.api.history.undo({});
    await v.api.history.undo({});
    const after = await run(v, 'count_atoms lig');
    expect(!after.ok && after.message).toContain("No object or selection is named 'lig'");
  });

  it('measures and labels single atoms', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'distance /1cbs//A/45/CA, /1cbs//A/60/CA');
    await ok(v, 'label /1cbs//A/45/CA, resn+resi');
    const list = await v.api.measure!.list({});
    expect(list.ok && list.did).toMatch(/distance/);
    expect(list.ok && list.did).toMatch(/label/);
    await ok(v, "label all, ''");
    const after = await v.api.measure!.list({});
    expect(after.ok && after.did).not.toMatch(/label/);
  });

  it('measures every pair within a cutoff (distance name, s1, s2, cutoff), and says when there is none', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    // The carboxylate's polar contacts: four O–N/O pairs within 4 Å in 1CBS.
    const r = await ok(v, 'distance polar, resn REA and elem O, polymer and (elem N or elem O), 4.0');
    const n = Number(r.did.match(/Added (\d+) distance/)?.[1]);
    expect(n, r.did).toBe(4);
    const list = await v.api.measure!.list({});
    expect(list.ok && (list.did.match(/distance/g) ?? []).length).toBeGreaterThanOrEqual(n);
    expect((await ok(v, 'distance none, resn REA, resi 1, 2.0')).did).toMatch(/No atoms/);
    const many = await run(v, 'distance all, polymer, polymer, 4');
    expect(many.ok).toBe(false);
  });
});

describe('VMD commands on the viewer', () => {
  it('loads with VMD’s look and keeps reps by index', async () => {
    const v = await viewer();
    await ok(v, 'mol new 1cbs.cif', 'vmd');
    expect(components(v).map((c) => c.types)).toEqual([['line']]);
    await ok(
      v,
      'mol delrep 0 top\nmol representation NewCartoon\nmol color Chain\nmol selection protein\nmol addrep top',
      'vmd',
    );
    expect(components(v).map((c) => c.types)).toEqual([['cartoon']]);
    const rep = v.plugin.managers.structure.hierarchy.current.structures[0].components[0].representations[0];
    expect((rep.cell.transform.params as { colorTheme: { name: string } }).colorTheme.name).toBe('chain-id');
    await ok(v, 'mol modstyle 0 top VDW', 'vmd');
    expect(components(v).map((c) => c.types)).toEqual([['spacefill']]);
    expect(
      (
        v.plugin.managers.structure.hierarchy.current.structures[0].components[0].representations[0].cell.transform
          .params as {
          colorTheme: { name: string };
        }
      ).colorTheme.name,
    ).toBe('chain-id');
    const bad = await run(v, 'mol delrep 3 top', 'vmd');
    expect(!bad.ok && bad.message).toContain('has no rep 3');
  });

  it('reads atomselect handles across submissions', async () => {
    const v = await viewer();
    await ok(v, 'mol new 1cbs.cif', 'vmd');
    await ok(v, 'set sel [atomselect top "resname REA and not hydrogen"]', 'vmd');
    const r = await ok(v, '$sel num', 'vmd');
    expect((r.value as ScriptRunValue).outputs.map((o) => o.text)).toEqual(['22']);
  });
});

describe('the menu’s rows', () => {
  it('say what a line will do, or what does not read, and complete', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const rows = prefixRows(v.api, 'pymol ', 'pymol', 'show sticks, resn REA');
    expect(rows[0].label).toBe("Show 'resn REA' as sticks");
    const bad = prefixRows(v.api, 'pymol ', 'pymol', 'shw sticks');
    expect(bad[0].label).toContain("PyMOL has no command 'shw'");
    const completing = prefixRows(v.api, 'pymol ', 'pymol', 'show sticks, chain ');
    expect(completing.map((r) => r.label)).toContain('A');
    expect(completing.find((r) => r.label === 'A')!.run()).toEqual({ text: 'pymol show sticks, chain A' });
  });
});

/** The style layers' lines (`representation.layers`). */
async function layers(v: PrismViewer) {
  const r = await v.api.representation.layers({});
  return r.ok ? r.did : '';
}

describe('what a first review found', () => {
  it('colors by element, whatever elements the selection has', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'util.cbag resn REA');
    await ok(v, 'color atomic, resn REA');
    await ok(v, 'util.cnc');
    const none = await run(v, 'util.cbag resn XYZ');
    expect(!none.ok && none.message).toContain('selects nothing');
  });

  it('shows a theme over everything after a color: the colors set before are cleared', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'color red');
    expect(await layers(v)).toContain('#ff0000');
    await ok(v, 'util.cbc');
    expect(await layers(v)).not.toContain('#ff0000');
    // On part of a component, the theme covers the whole component: said.
    const part = await ok(v, 'spectrum count, rainbow, resi 1-50');
    expect(part.warnings?.join(' ')).toContain('whole components');
  });

  it('does not throw from a preview when a substitution does not read', async () => {
    const v = await viewer();
    const r = v.api.script!.preview('set s [atomselect top {protein}x]', 'vmd');
    expect(r.diagnostics[0]).toMatchObject({ severity: 'error', code: 'syntax' });
  });

  it('keeps a rep’s own params when its material changes', async () => {
    const v = await viewer();
    await ok(v, 'mol new 1cbs.cif\nmol representation Licorice\nmol addrep top', 'vmd');
    await ok(v, 'mol modmaterial 1 top Transparent', 'vmd');
    const reps = v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap(
      (c) => c.representations,
    );
    const params = (reps[1].cell.transform.params as { type: { params: { sizeAspectRatio: number; alpha: number } } })
      .type.params;
    expect(params.sizeAspectRatio).toBe(1);
    expect(params.alpha).toBe(0.5);
  });

  it('loads a second copy for VMD’s mol new, and the names follow', async () => {
    const v = await viewer();
    await ok(v, 'mol new 1cbs.cif', 'vmd');
    await ok(v, 'mol new 1cbs.cif\nmol addrep top', 'vmd');
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(2);
  });

  it('fixes a named selection to the objects there were', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    await ok(v, 'select prot, polymer, enable=0');
    await ok(v, 'load 1cbs.cif, other');
    const r = await ok(v, 'count_atoms prot');
    const n = (r.value as ScriptRunValue).outputs.map((o) => o.text)[0];
    const all = await ok(v, 'count_atoms polymer');
    expect(Number(n)).toBeLessThan(Number((all.value as ScriptRunValue).outputs[0].text.split(' ')[0]));
  });

  it('refuses pair_fit pairs in another order than the structures hold them', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif\nload 1cbs.cif, other');
    const bad = await run(v, 'pair_fit /1cbs//A/10/CA, /other//A/20/CA, /1cbs//A/5/CA, /other//A/30/CA');
    expect(!bad.ok && bad.message).toContain('not in the same order');
  });

  it('says both size warnings', async () => {
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const r = await ok(v, 'set sphere_scale, 20');
    expect(r.warnings).toHaveLength(2);
  });
});

describe('fetch of a name another structure has', () => {
  it('names the new structure apart, and says so', async () => {
    // The PDB's file for 1CBS, without the network: any address naming it answers the example file.
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(input instanceof Request ? input.url : input);
        return /1cbs/i.test(url) ? new Response(cif) : new Response('', { status: 404 });
      }),
    );
    const v = await viewer();
    await ok(v, 'load 1cbs.cif');
    const r = await ok(v, 'fetch 1cbs');
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(2);
    expect(r.warnings?.join(' ')).toContain('this one is 1cbs_2');
    const names = await v.api.script!.names({});
    expect(names.ok && names.did).toMatch(/1cbs, 1 1cbs_2/);
  });
});

describe('the assistant’s run', () => {
  it('says what its loads reach, and refuses what would ask the person', async () => {
    const v = await viewer();
    const def = v.api.describe('script.run')!.definition;
    expect(def.approval).toBe('network');
    expect(def.addresses!({ dialect: 'pymol', text: 'fetch 1cbs; load https://example.org/x.cif', files: [] })).toEqual(
      ['https://example.org/x.cif'],
    );
    await ok(v, 'load 1cbs.cif');
    const agent = v.api.as({ kind: 'agent', turn: 't1' });
    const removed = await agent.script!.run({ dialect: 'pymol', text: 'delete 1cbs' });
    expect(!removed.ok && removed.code).toBe('blocked');
    const shown = await agent.script!.run({ dialect: 'pymol', text: 'show sticks, resi 10-20' });
    expect(shown.ok).toBe(true);
  });

  it('is a skill with one tool over script.run', async () => {
    const { CommandLanguageSkill } = await import('../skill.js');
    expect(CommandLanguageSkill.tools.map((t) => t.operation)).toEqual(['script.run']);
  });
});

describe('the assistant’s run and Tcl variables', () => {
  it('starts without variables each call, so none can name what its approval did not see', async () => {
    const v = await viewer();
    await ok(v, 'mol new 1cbs.cif', 'vmd');
    // The person's variables, and an earlier call's, are not the assistant's.
    await ok(v, 'set u https://evil.example/x.pdb; set m 0', 'vmd');
    const agent = v.api.as({ kind: 'agent', turn: 't1' });
    const load = await agent.script!.run({ dialect: 'vmd', text: 'mol new $u' });
    expect(!load.ok && load.message).toContain("no variable 'u'");
    const removed = await agent.script!.run({ dialect: 'vmd', text: 'mol delete $m' });
    expect(removed.ok).toBe(false);
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(1);
    // A variable set in the same text is seen by the approval too.
    const def = v.api.describe('script.run')!.definition;
    expect(def.addresses!({ dialect: 'vmd', text: 'set u https://a.example/x.pdb; mol new $u', files: [] })).toEqual([
      'https://a.example/x.pdb',
    ]);
  });
});
