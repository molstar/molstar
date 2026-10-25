/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { EmptyEnv, type ScriptEnv } from '../env.js';
import type { Command } from '../ir.js';
import { PymolCommands } from '../pymol/commands.js';
import { splitPymol } from '../pymol/tokenizer.js';
import { PymolUnsupported } from '../pymol/unsupported.js';
import { compileScript, hasErrors } from '../script.js';

// PyMOL's dialect (_spec §3.1, pymol.md): the tokenizer, every table example, the diagnostics, names and objects.

const env: ScriptEnv = {
  ...EmptyEnv,
  objects: [
    { name: '1cbs', molecule: 0 },
    { name: '1tqn', molecule: 1 },
  ],
  nextMolecule: 2,
  files: ['model.pdb'],
};

const compile = (text: string, e: ScriptEnv = env) => compileScript(text, 'pymol', e);

function ok(text: string, e: ScriptEnv = env): readonly Command[] {
  const r = compile(text, e);
  const errors = r.diagnostics.filter((d) => d.severity === 'error');
  if (errors.length) expect.fail(`${text}: ${errors.map((d) => `${d.code}: ${d.message}`).join('; ')}`);
  return r.commands;
}

function errorOf(text: string, e: ScriptEnv = env) {
  const r = compile(text, e);
  const d = r.diagnostics.find((x) => x.severity === 'error');
  if (!d) expect.fail(`${text}: no error`);
  return d;
}

describe('the PyMOL tokenizer', () => {
  it('splits commands on ; and newlines, arguments on commas at depth 0', () => {
    const { calls } = splitPymol('show sticks, (resn HEM, resi 1); color red\nhide lines # gone');
    expect(calls.map((c) => c.name)).toEqual(['show', 'color', 'hide']);
    expect(calls[0].args.map((a) => a.text)).toEqual(['sticks', '(resn HEM, resi 1)']);
    expect(calls[2].args.map((a) => a.text)).toEqual(['lines']);
  });

  it('keeps places, reads key=value, and joins continued lines', () => {
    const text = 'select near, lig around 5, \\\n  enable=0';
    const { calls } = splitPymol(text);
    expect(calls).toHaveLength(1);
    const [name, sel, en] = calls[0].args;
    expect(text.slice(name.at.from, name.at.to)).toBe('near');
    expect(sel.text).toBe('lig around 5');
    expect(en.key).toBe('enable');
    expect(en.text).toBe('0');
    // `b=1`-like selections stay positional when the key is not the command's (`PymolArgs` decides).
    expect(splitPymol('color red, b = 1').calls[0].args[1].key).toBe('b');
  });

  it('reads quoted arguments with commas and semicolons inside', () => {
    const { calls } = splitPymol(`label name CA, "a, b; c"`);
    expect(calls).toHaveLength(1);
    expect(calls[0].args[1].text).toBe('"a, b; c"');
  });

  it('says an unclosed quote is one', () => {
    expect(splitPymol(`label all, "oops`).diagnostics[0].code).toBe('syntax');
  });

  it('marks a line that ends after a comma or space as open (completion offers the next argument)', () => {
    expect(splitPymol('show sticks, ').calls[0].open).toBe(true);
    expect(splitPymol('show ').calls[0].open).toBe(true);
    expect(splitPymol('show sticks').calls[0].open).toBe(false);
  });
});

describe('the PyMOL table', () => {
  // Names the examples use: two objects, a named selection.
  const examples = {
    ...env,
    names: new Map(
      ok('select lig, resn HEM', env)
        .filter((c) => c.kind === 'name.set')
        .map((c) => [(c as { name: string }).name, (c as { sel: never }).sel]),
    ),
  };

  for (const c of PymolCommands) {
    for (const example of c.examples) {
      it(`compiles ${example}`, () => {
        ok(example, { ...examples, objects: [...env.objects, { name: 'heme', molecule: 2 }] });
      });
    }
  }

  it('fails every command Mol* does not do, saying so', () => {
    for (const u of PymolUnsupported) {
      const d = errorOf(u.name === 'python' ? 'python\nprint(1)\npython end' : `${u.name} x`);
      expect(d.code, u.name).toBe('unsupported');
    }
  });
});

describe('PyMOL commands to IR', () => {
  it('fetches into named objects, the names known to the next command', () => {
    const cmds = ok('fetch 4hhb; show sticks, /4hhb//A');
    expect(cmds[0]).toMatchObject({
      kind: 'load',
      name: '4hhb',
      molecule: 2,
      from: { entry: { db: 'pdb', id: '4hhb' } },
    });
    const display = cmds[1] as Extract<Command, { kind: 'display' }>;
    expect(display.sel.parts.map((p) => (p.scope === 'others' ? 'others' : p.scope.object))).toEqual([
      '4hhb',
      'others',
    ]);
  });

  it('fetches extended PDB ids: of the four-character era by those four, a later one as written', () => {
    expect(ok('fetch pdb_00001cbs')[0]).toMatchObject({ from: { entry: { db: 'pdb', id: '1cbs' } } });
    expect(ok('fetch PDB_00012ABC')[0]).toMatchObject({ from: { entry: { db: 'pdb', id: 'pdb_00012abc' } } });
    for (const bad of ['fetch 1cb', 'fetch pdb_1cbs', 'fetch pdb_0000abcd']) {
      expect(errorOf(bad)?.hint, bad).toContain('pdb_00001cbs');
    }
  });

  it('inlines a named selection, keeping what it meant when defined', () => {
    const cmds = ok('select lig, resn HEM; show sticks, lig and chain A; select lig, none');
    const shown = cmds.find((c) => c.kind === 'display') as Extract<Command, { kind: 'display' }>;
    expect(JSON.stringify(shown.sel.parts)).toContain('HEM');
  });

  it('reads `select expr` as the name sele, and selects it', () => {
    const cmds = ok('select resn HEM');
    expect(cmds.map((c) => c.kind)).toEqual(['name.set', 'select']);
    expect((cmds[0] as { name: string }).name).toBe('sele');
  });

  it('compiles as, show and hide with PyMOL’s defaults', () => {
    expect(ok('show')[0]).toMatchObject({ kind: 'display', mode: 'show', reps: ['lines'] });
    expect(ok('hide')[0]).toMatchObject({ kind: 'display', mode: 'hide', reps: 'all' });
    expect(ok('as cartoon')[0]).toMatchObject({ kind: 'display', mode: 'only', reps: ['cartoon'] });
    expect(ok('show stick, resn HEM')[0]).toMatchObject({ reps: ['sticks'] });
  });

  it('colors carbons and the rest by element for util.cbag', () => {
    const cmds = ok('util.cbag resn HEM');
    // The selection must select atoms; each element's color is one of the group (selecting none is fine).
    expect(cmds[0]).toMatchObject({ kind: 'require' });
    expect(cmds[1]).toMatchObject({ kind: 'color', color: { rgb: 0x33ff33 }, optional: true });
    expect(cmds.length).toBeGreaterThan(3);
  });

  it('reads labels’ expressions', () => {
    const cmds = ok('label name CA, "%s-%s" % (resn, resi)');
    expect(cmds[0]).toMatchObject({
      kind: 'label',
      text: [{ property: 'resn' }, { text: '-' }, { property: 'resi' }],
    });
    expect(ok('label name CA, resn+resi')[0]).toMatchObject({ text: [{ property: 'resn' }, { property: 'resi' }] });
    expect(ok("label all, ''")[0].kind).toBe('label.clear');
    expect(errorOf('label name CA, resn.lower()').code).toBe('unsupported');
  });

  it('reads distance with and without a name', () => {
    expect((ok('distance A/45/CA, A/60/CA')[0] as { atoms: unknown[] }).atoms).toHaveLength(2);
    expect((ok('distance d1, A/45/CA, A/60/CA')[0] as { atoms: unknown[] }).atoms).toHaveLength(2);
    expect(errorOf('distance A/45/CA, A/60/CA, mode=2').code).toBe('unsupported-argument');
  });

  it('warns about what it ignores, and still compiles', () => {
    const r = compile('png out.png, 800, 600, ray=1');
    expect(hasErrors(r.diagnostics)).toBe(false);
    expect(r.diagnostics.some((d) => d.severity === 'note' && d.message.includes('ray'))).toBe(true);
    expect(r.commands[0]).toMatchObject({ kind: 'export.image', name: 'out', width: 800, height: 600 });
    const a = compile('align 1tqn, 1cbs, cycles=5');
    expect(a.diagnostics.map((d) => d.severity)).toEqual(['warning']);
  });

  it('sets the settings Mol* has, and fails others listing them', () => {
    expect(ok('set bg_rgb, white')[0]).toMatchObject({ kind: 'canvas', background: 0xffffff });
    expect(ok('set orthoscopic, on')[0]).toMatchObject({ kind: 'canvas', projection: 'orthographic' });
    expect(ok('set stick_radius, 0.5')[0]).toMatchObject({ kind: 'size', factor: 2, reps: ['sticks'] });
    expect(ok('set ortho, 1')[0]).toMatchObject({ projection: 'orthographic' });
    const d = errorOf('set ray_trace_frames, 1');
    expect(d.code).toBe('unsupported');
    expect(d.hint).toContain('bg_rgb');
  });

  it('loads a file given with the script, and fails one that is not', () => {
    expect(ok('load model.pdb')[0]).toMatchObject({ kind: 'load', from: { file: 'model.pdb' }, name: 'model' });
    expect(errorOf('load other.pdb').code).toBe('no-name');
  });

  it('says what an unknown command is, and what it may have been', () => {
    const d = errorOf('shw sticks');
    expect(d.code).toBe('unknown-command');
    expect(d.hint).toContain('show');
    expect(errorOf('mol new 1cbs').hint).toContain('vmd mol');
  });

  it('places a selection’s error on its argument', () => {
    const text = 'show sticks, resn HEM and and';
    const d = errorOf(text);
    expect(d.code).toBe('selection');
    expect(text.slice(d.at.from, d.at.to)).toBe('resn HEM and and');
  });

  it('names what an unknown name could have been', () => {
    const d = errorOf('show sticks, ligand_x');
    expect(d.code).toBe('no-name');
    expect(d.hint).toContain('1cbs');
  });

  it('runs undo alone only', () => {
    expect(ok('undo')[0]).toMatchObject({ kind: 'history', action: 'undo' });
    expect(errorOf('show sticks; undo').message).toContain('runs alone');
  });

  it('reports all errors at once', () => {
    const r = compile('shw\ncolor notacolor\nshow sticks');
    expect(r.diagnostics.filter((d) => d.severity === 'error')).toHaveLength(2);
  });

  it('prints help from the table', () => {
    expect(ok('help show')[0]).toMatchObject({ kind: 'print' });
    expect((ok('help')[0] as { text: string }).text).toContain('fetch');
  });
});

describe('PyMOL scripts that run scripts', () => {
  it('runs @file.pml and load file.pml, and fails one not given', () => {
    const scripts = new Map([['style.pml', 'hide everything\nshow cartoon']]);
    const r = compileScript('@style.pml\nload style.pml', 'pymol', { ...env, scripts });
    expect(r.diagnostics).toEqual([]);
    expect(r.commands.map((c) => c.kind)).toEqual(['display', 'display', 'display', 'display']);
    expect(errorOf('@missing.pml').code).toBe('no-name');
  });
});

describe('selection names', () => {
  it('takes names that read as themselves, and refuses those PyMOL’s words take', () => {
    for (const name of ['lig', 'near', 'prot', 'site_1', 'pocket', 'my_sel', 'b2', 'both'])
      ok(`select ${name}, resn HEM`);
    for (const name of ['b', 'chain', 'all']) expect(errorOf(`select ${name}, resn HEM`).code, name).toBe('argument');
  });
});

describe('selections that do not read', () => {
  it('says where, for unfinished text too', () => {
    const d = errorOf('show sticks, chain ');
    expect(d.code).toBe('selection');
    expect(d.message).toContain('does not read as a PyMOL selection');
  });
});
