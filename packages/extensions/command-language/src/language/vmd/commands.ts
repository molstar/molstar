/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { vmdColor } from '../colors.js';
import { fail, note, type Span, warning } from '../diagnostics.js';
import { PdbIdHint, pdbId } from '../entries.js';
import type { ObjectInfo } from '../env.js';
import type { CommandBody, Sel } from '../ir.js';
import { vmdSelection } from '../selection.js';
import type { CommandHelp, Compiled, Compiler, ParamSpec } from '../table.js';
import { VmdColorMethodNames, VmdStyleNames, vmdColorMethod, vmdMaterial, vmdStyle } from './styles.js';

// VMD's commands (vmd.md §1): `mol …`, `atomselect` and its handles, `measure`, `label`, `color Display`, `display`,
// `rotate`/`scale`, `animate`, `render`, `puts`. Each entry is a command and its subcommand (`mol addrep`); the words
// after them are its arguments, already substituted (`dialect.ts`).

/** A command's argument words, after its name, with their places. */
export class VmdWords {
  constructor(
    readonly values: readonly string[],
    readonly ats: readonly Span[],
    readonly at: Span,
    readonly name: string,
  ) {}
  get length() {
    return this.values.length;
  }
  has(i: number) {
    return i < this.values.length;
  }
  get(i: number, what: string): string {
    if (i >= this.values.length) fail('argument', `${this.name} needs ${what}.`, this.at);
    return this.values[i];
  }
  opt(i: number): string | undefined {
    return this.values[i];
  }
  atOf(i: number): Span {
    return this.ats[i] ?? this.at;
  }
  /** The words from `i` on, as one text (VMD joins a selection's words). */
  rest(i: number): string {
    return this.values.slice(i).join(' ');
  }
  restAt(i: number): Span {
    if (i >= this.ats.length) return this.at;
    return { from: this.ats[i].from, to: this.ats[this.ats.length - 1].to };
  }
  int(i: number, what: string): number {
    const t = this.get(i, what);
    if (!/^-?\d+$/.test(t)) fail('argument', `${what} is a whole number, not '${t}'.`, this.atOf(i));
    return Number(t);
  }
  number(i: number, what: string): number {
    const t = this.get(i, what);
    const v = Number(t);
    if (!Number.isFinite(v)) fail('argument', `${what} is a number, not '${t}'.`, this.atOf(i));
    return v;
  }
  /** The words that are flags of a command (`type pdb`, `waitfor all`): name to value, from `i` on. */
  options(i: number): Map<string, { value: string; at: Span }> {
    const out = new Map<string, { value: string; at: Span }>();
    for (let k = i; k + 1 < this.values.length; k += 2)
      out.set(this.values[k].toLowerCase(), { value: this.values[k + 1], at: this.ats[k] });
    return out;
  }
}

export interface VmdCommand extends CommandHelp {
  readonly params: readonly ParamSpec[];
  compile(w: VmdWords, c: Compiler): Compiled;
}

const p = (name: string, kind: ParamSpec['kind'], optional = false, values?: readonly string[]): ParamSpec =>
  values ? { name, kind, optional, values } : { name, kind, optional };

/** A molecule by id or `top`. */
export function moleculeOf(c: Compiler, id: string, at: Span): ObjectInfo {
  const m = c.state.molecule(id);
  if (m) return m;
  const ids = c.state.objects.map((o) => `${o.molecule} (${o.name})`);
  return fail(
    'no-name',
    id === 'top' ? 'There is no top molecule: nothing is loaded.' : `There is no molecule ${id}.`,
    at,
    ids.length ? `Molecules: ${ids.join(', ')}.` : 'Load one with `mol new` or `mol pdbload`.',
  );
}

function repIndex(w: VmdWords, i: number): number {
  const n = w.int(i, 'a rep index');
  if (n < 0) fail('argument', 'Rep indices start at 0.', w.atOf(i));
  return n;
}

/** An atom by `molid/index` (VMD's labels) or `index` on a molecule: VMD's 0-based atom index. */
function atomSel(c: Compiler, spec: string, at: Span, molecule?: string): Sel {
  const m = spec.match(/^(?:(\w+)\/)?(\d+)$/);
  if (!m)
    fail(
      'argument',
      `'${spec}' is not an atom (molid/index).`,
      at,
      'Atoms are written 0/12: molecule 0, atom index 12.',
    );
  const object = moleculeOf(c, m[1] ?? molecule ?? 'top', at);
  const sel = vmdSelection(`index ${m[2]}`, object.name, at);
  return { ...sel, text: `${object.molecule}/${m[2]}` };
}

/** `{12 15}` or `12 15`: VMD atom indices. */
function indices(text: string, at: Span, n: number): string[] {
  const xs = text.split(/\s+/).filter(Boolean);
  if (xs.length !== n || xs.some((x) => !/^\d+$/.test(x))) {
    fail(
      'argument',
      `This takes ${n} atom indices, as {${Array.from({ length: n }, (_, i) => 10 + i).join(' ')}}.`,
      at,
    );
  }
  return xs;
}

const LabelKinds: Readonly<Record<string, 'atom' | 'distance' | 'angle' | 'dihedral'>> = {
  atoms: 'atom',
  bonds: 'distance',
  angles: 'angle',
  dihedrals: 'dihedral',
};

/** What `mol new` loads: a URL, a file given with the script, or an id (`type webpdb`, or a bare 4-character id). */
function newMolecule(w: VmdWords, c: Compiler, i: number, webpdb: boolean): CommandBody {
  const file = w.get(i, 'a file');
  const options = w.options(i + 1);
  const type = options.get('type')?.value.toLowerCase();
  for (const [k, v] of options) {
    if (k !== 'type' && k !== 'waitfor')
      c.report(warning('unsupported-argument', `${k} is ignored: the whole file loads.`, v.at));
  }
  const id = pdbId(file);
  const given = c.state.files.some((f) => f === file || f.split('/').pop() === file.split('/').pop());
  const url = /^(https?|ftp):\/\//i.test(file);
  const molecule = c.state.nextMolecule;
  if (webpdb || type === 'webpdb' || (id && !given && !url)) {
    if (!id) fail('argument', `'${file}' is not a PDB id.`, w.atOf(i), PdbIdHint);
    if (!webpdb && type !== 'webpdb') c.report(note(`${file} is read as a PDB id (there is no such file).`, w.atOf(i)));
    return {
      kind: 'load',
      from: { entry: { db: 'pdb', id } },
      name: c.state.freeName(file.toLowerCase()),
      molecule,
      look: 'vmd-lines',
      copy: true,
    };
  }
  if (!url && !given) {
    fail(
      'no-name',
      `There is no file '${file}' with this script.`,
      w.atOf(i),
      c.state.files.length
        ? `Files given: ${c.state.files.join(', ')}.`
        : 'Open the file with the script, or give a URL.',
    );
  }
  const stem = (file.split(/[\\/]/).pop() ?? file).replace(/\.(gz|zip)$/i, '').replace(/\.[^.]+$/, '');
  return {
    kind: 'load',
    from: url ? { url: file } : { file },
    name: c.state.freeName(stem),
    molecule,
    look: 'vmd-lines',
    copy: true,
    ...(type && type !== 'webpdb' ? { format: type } : {}),
  };
}

export const VmdCommands: readonly VmdCommand[] = [
  {
    name: 'mol new',
    summary: 'Load a file given with the script, a URL, or a PDB id (`type webpdb`), with VMD’s Lines look.',
    usage: 'mol new file [type t] [waitfor all]',
    params: [p('file', 'file'), p('type', 'enum', true, ['type']), p('format', 'enum', true, ['pdb', 'cif', 'webpdb'])],
    examples: ['mol new 1cbs type webpdb', 'mol new https://files.rcsb.org/download/1cbs.cif'],
    compile: (w, c) => [newMolecule(w, c, 0, false)],
  },
  {
    name: 'mol pdbload',
    summary: 'Load a PDB entry, with VMD’s Lines look.',
    usage: 'mol pdbload id',
    params: [p('id', 'text')],
    examples: ['mol pdbload 1cbs', 'mol pdbload pdb_00001cbs'],
    compile: (w, c) => [newMolecule(w, c, 0, true)],
  },
  {
    name: 'mol delete',
    summary: 'Remove a molecule (`all`: every one).',
    usage: 'mol delete molid|all',
    params: [p('molid', 'molecule')],
    examples: ['mol delete top', 'mol delete all'],
    compile(w, c) {
      const id = w.get(0, 'a molecule');
      if (id === 'all') return [{ kind: 'object.remove', objects: 'all' }];
      return [{ kind: 'object.remove', objects: [moleculeOf(c, id, w.atOf(0)).name] }];
    },
  },
  ...(['on', 'off'] as const).map(
    (state): VmdCommand => ({
      name: `mol ${state}`,
      summary: state === 'on' ? 'Show a molecule.' : 'Hide a molecule.',
      usage: `mol ${state} molid`,
      params: [p('molid', 'molecule')],
      examples: [`mol ${state} top`],
      compile: (w, c) => [
        {
          kind: 'object.visibility',
          objects: [moleculeOf(c, w.get(0, 'a molecule'), w.atOf(0)).name],
          visible: state === 'on',
        },
      ],
    }),
  ),
  {
    name: 'mol top',
    summary: 'Make a molecule the top one (what `top` means).',
    usage: 'mol top molid',
    params: [p('molid', 'molecule')],
    examples: ['mol top 0'],
    compile: (w, c) => [{ kind: 'object.top', object: moleculeOf(c, w.get(0, 'a molecule'), w.atOf(0)).name }],
  },
  {
    name: 'mol representation',
    aliases: ['mol rep'],
    summary: 'Set the style the next `mol addrep` takes.',
    usage: 'mol representation Style [params…]',
    params: [p('style', 'rep', false, VmdStyleNames)],
    examples: ['mol representation NewCartoon', 'mol representation Licorice 0.3 12 12'],
    compile(w, c) {
      const style = vmdStyle(w.get(0, 'a style'), w.atOf(0), c);
      if (w.length > 1)
        c.report(warning('unsupported-argument', 'The style’s numbers are ignored: Mol*’s sizes apply.', w.restAt(1)));
      return [{ kind: 'rep.defaults', patch: { style } }];
    },
  },
  {
    name: 'mol color',
    summary: 'Set the coloring the next `mol addrep` takes.',
    usage: 'mol color Method [args]',
    params: [p('method', 'color', false, VmdColorMethodNames)],
    examples: ['mol color Chain', 'mol color ColorID 1', 'mol color Structure'],
    compile: (w, c) => [
      {
        kind: 'rep.defaults',
        patch: { color: vmdColorMethod(w.values, w.restAt(0), c.state.molecule('top')?.molecule) },
      },
    ],
  },
  {
    name: 'mol selection',
    aliases: ['mol sel'],
    summary: 'Set the selection the next `mol addrep` takes.',
    usage: 'mol selection {selection}',
    params: [p('selection', 'selection')],
    examples: ['mol selection {protein and chain A}', 'mol selection "resname HEM"'],
    compile(w) {
      if (!w.length) fail('argument', 'mol selection needs a selection.', w.at);
      // Read now so a typo shows here, not at the next addrep (the selection is compiled again per molecule there).
      vmdSelection(w.rest(0), '', w.restAt(0));
      return [{ kind: 'rep.defaults', patch: { selection: w.rest(0) } }];
    },
  },
  {
    name: 'mol material',
    summary: 'Set the material the next `mol addrep` takes (Opaque, Transparent).',
    usage: 'mol material Material',
    params: [p('material', 'enum', false, ['Opaque', 'Transparent'])],
    examples: ['mol material Transparent'],
    compile: (w, c) => [
      { kind: 'rep.defaults', patch: { material: vmdMaterial(w.get(0, 'a material'), w.atOf(0), c) } },
    ],
  },
  {
    name: 'mol addrep',
    summary: 'Add a rep to a molecule with the current style, coloring, selection and material.',
    usage: 'mol addrep molid',
    params: [p('molid', 'molecule')],
    examples: ['mol addrep top'],
    compile(w, c) {
      const m = moleculeOf(c, w.get(0, 'a molecule'), w.atOf(0));
      vmdSelection(c.state.repDefaults.selection, m.name, w.at);
      return [{ kind: 'rep.add', object: m.name, rep: c.state.repDefaults }];
    },
  },
  {
    name: 'mol delrep',
    summary: 'Remove a rep of a molecule.',
    usage: 'mol delrep repid molid',
    params: [p('repid', 'int'), p('molid', 'molecule')],
    examples: ['mol delrep 0 top'],
    compile: (w, c) => [
      { kind: 'rep.remove', object: moleculeOf(c, w.get(1, 'a molecule'), w.atOf(1)).name, index: repIndex(w, 0) },
    ],
  },
  {
    name: 'mol modstyle',
    summary: 'Change a rep’s style.',
    usage: 'mol modstyle repid molid Style [params…]',
    params: [p('repid', 'int'), p('molid', 'molecule'), p('style', 'rep', false, VmdStyleNames)],
    examples: ['mol modstyle 0 top NewCartoon'],
    compile(w, c) {
      const style = vmdStyle(w.get(2, 'a style'), w.atOf(2), c);
      if (w.length > 3)
        c.report(warning('unsupported-argument', 'The style’s numbers are ignored: Mol*’s sizes apply.', w.restAt(3)));
      return [
        {
          kind: 'rep.modify',
          object: moleculeOf(c, w.get(1, 'a molecule'), w.atOf(1)).name,
          index: repIndex(w, 0),
          patch: { style },
        },
      ];
    },
  },
  {
    name: 'mol modcolor',
    summary: 'Change a rep’s coloring.',
    usage: 'mol modcolor repid molid Method [args]',
    params: [p('repid', 'int'), p('molid', 'molecule'), p('method', 'color', false, VmdColorMethodNames)],
    examples: ['mol modcolor 0 top Chain', 'mol modcolor 0 top ColorID 4'],
    compile(w, c) {
      const m = moleculeOf(c, w.get(1, 'a molecule'), w.atOf(1));
      const color = vmdColorMethod(w.values.slice(2), w.restAt(2), m.molecule);
      return [{ kind: 'rep.modify', object: m.name, index: repIndex(w, 0), patch: { color } }];
    },
  },
  {
    name: 'mol modselect',
    summary: 'Change a rep’s selection.',
    usage: 'mol modselect repid molid {selection}',
    params: [p('repid', 'int'), p('molid', 'molecule'), p('selection', 'selection')],
    examples: ['mol modselect 0 top {protein}'],
    compile(w, c) {
      const m = moleculeOf(c, w.get(1, 'a molecule'), w.atOf(1));
      if (w.length < 3) fail('argument', 'mol modselect needs a selection.', w.at);
      vmdSelection(w.rest(2), m.name, w.restAt(2));
      return [{ kind: 'rep.modify', object: m.name, index: repIndex(w, 0), patch: { selection: w.rest(2) } }];
    },
  },
  {
    name: 'mol modmaterial',
    summary: 'Change a rep’s material.',
    usage: 'mol modmaterial repid molid Material',
    params: [p('repid', 'int'), p('molid', 'molecule'), p('material', 'enum', false, ['Opaque', 'Transparent'])],
    examples: ['mol modmaterial 0 top Transparent'],
    compile: (w, c) => [
      {
        kind: 'rep.modify',
        object: moleculeOf(c, w.get(1, 'a molecule'), w.atOf(1)).name,
        index: repIndex(w, 0),
        patch: { material: vmdMaterial(w.get(2, 'a material'), w.atOf(2), c) },
      },
    ],
  },
  {
    name: 'mol showrep',
    summary: 'Show or hide a rep.',
    usage: 'mol showrep molid repid on|off',
    params: [p('molid', 'molecule'), p('repid', 'int'), p('state', 'enum', true, ['on', 'off'])],
    examples: ['mol showrep top 0 off'],
    compile(w, c) {
      const state = (w.opt(2) ?? 'on').toLowerCase();
      if (!['on', 'off', '1', '0'].includes(state)) fail('argument', `'${state}' is not on or off.`, w.atOf(2));
      return [
        {
          kind: 'rep.visibility',
          object: moleculeOf(c, w.get(0, 'a molecule'), w.atOf(0)).name,
          index: repIndex(w, 1),
          visible: state === 'on' || state === '1',
        },
      ];
    },
  },
  {
    name: 'mol list',
    summary: 'List the molecules.',
    usage: 'mol list',
    params: [],
    examples: ['mol list'],
    compile: () => [{ kind: 'read', read: { kind: 'objects' } }],
  },
  {
    name: 'measure bond',
    summary: 'Print the distance between two atoms (VMD indices).',
    usage: 'measure bond {i j} [molid m]',
    params: [p('atoms', 'text'), p('molid', 'enum', true, ['molid'])],
    examples: ['measure bond {0 1}'],
    compile: (w, c) => geometry(w, c, 'distance', 2),
  },
  {
    name: 'measure angle',
    summary: 'Print the angle of three atoms (VMD indices).',
    usage: 'measure angle {i j k} [molid m]',
    params: [p('atoms', 'text'), p('molid', 'enum', true, ['molid'])],
    examples: ['measure angle {0 1 2}'],
    compile: (w, c) => geometry(w, c, 'angle', 3),
  },
  {
    name: 'measure dihed',
    summary: 'Print the dihedral angle of four atoms (VMD indices).',
    usage: 'measure dihed {i j k l} [molid m]',
    params: [p('atoms', 'text'), p('molid', 'enum', true, ['molid'])],
    examples: ['measure dihed {0 1 2 3}'],
    compile: (w, c) => geometry(w, c, 'dihedral', 4),
  },
  ...(['center', 'rgyr', 'minmax'] as const).map(
    (what): VmdCommand => ({
      name: `measure ${what}`,
      summary:
        what === 'center'
          ? 'Print the center of an atom selection.'
          : what === 'rgyr'
            ? 'Print the radius of gyration of an atom selection.'
            : 'Print the bounding box of an atom selection.',
      usage: `measure ${what} $sel`,
      params: [p('sel', 'text')],
      examples: [`measure ${what} $sel`],
      compile(w, c) {
        const handle = c.state.session.handles.get(w.get(0, 'an atom selection'));
        if (!handle)
          fail(
            'no-name',
            `'${w.get(0, 'an atom selection')}' is not an atom selection.`,
            w.atOf(0),
            'Make one: set sel [atomselect top "protein"].',
          );
        return [{ kind: 'read', read: { kind: 'geometry', measure: what, atoms: [handle.sel] } }];
      },
    }),
  ),
  {
    name: 'label add',
    summary: 'Label an atom, or measure a bond, angle or dihedral (atoms as molid/index).',
    usage: 'label add Atoms|Bonds|Angles|Dihedrals molid/index …',
    params: [p('category', 'enum', false, ['Atoms', 'Bonds', 'Angles', 'Dihedrals']), p('atoms', 'text')],
    examples: ['label add Atoms 0/12', 'label add Bonds 0/12 0/15'],
    compile(w, c) {
      const category = w.get(0, 'a category');
      const kind = LabelKinds[category.toLowerCase()];
      if (!kind)
        fail('argument', `'${category}' is not a label category.`, w.atOf(0), 'Atoms, Bonds, Angles or Dihedrals.');
      const n = { atom: 1, distance: 2, angle: 3, dihedral: 4 }[kind];
      const specs = w.values
        .slice(1)
        .flatMap((v) => v.split(/\s+/))
        .filter(Boolean);
      if (specs.length !== n) fail('argument', `label add ${category} takes ${n} atom(s), as 0/12.`, w.restAt(1));
      const atoms = specs.map((s) => atomSel(c, s, w.restAt(1)));
      if (kind === 'atom') {
        return [
          {
            kind: 'label',
            sel: atoms[0],
            text: [{ property: 'resn' }, { text: ' ' }, { property: 'resi' }, { text: ':' }, { property: 'name' }],
          },
        ];
      }
      return [{ kind: 'measure', measure: kind, atoms }];
    },
  },
  {
    name: 'label delete',
    summary: 'Remove the labels or measurements of a category.',
    usage: 'label delete Atoms|Bonds|Angles|Dihedrals [all]',
    params: [
      p('category', 'enum', false, ['Atoms', 'Bonds', 'Angles', 'Dihedrals']),
      p('which', 'enum', true, ['all']),
    ],
    examples: ['label delete Bonds all', 'label delete Atoms'],
    compile(w) {
      const category = w.get(0, 'a category');
      const kind = LabelKinds[category.toLowerCase()];
      if (!kind)
        fail('argument', `'${category}' is not a label category.`, w.atOf(0), 'Atoms, Bonds, Angles or Dihedrals.');
      const which = w.opt(1);
      if (which !== undefined && which.toLowerCase() !== 'all') {
        fail(
          'unsupported-argument',
          'Removing one label by its index is not supported.',
          w.atOf(1),
          'Remove them all: label delete … all.',
        );
      }
      return kind === 'atom' ? [{ kind: 'label.clear' }] : [{ kind: 'measure.clear', measure: kind }];
    },
  },
  {
    name: 'color Display',
    summary: 'Set the background color (`color Display Background white`).',
    usage: 'color Display Background color',
    params: [p('item', 'enum', false, ['Background']), p('color', 'color', false, ['white', 'black', 'blue', 'gray'])],
    examples: ['color Display Background white'],
    compile(w) {
      const item = w.get(0, 'an item');
      if (item.toLowerCase() !== 'background') {
        fail('unsupported', `color Display ${item} is not supported.`, w.atOf(0), 'Only the Background is.');
      }
      const rgb = vmdColor(w.get(1, 'a color'));
      if (rgb === undefined)
        fail(
          'argument',
          `'${w.get(1, 'a color')}' is not a VMD color.`,
          w.atOf(1),
          'Give a name (white, black) or a ColorID.',
        );
      return [{ kind: 'canvas', background: rgb }];
    },
  },
  {
    name: 'display projection',
    summary: 'Orthographic or Perspective projection.',
    usage: 'display projection Orthographic|Perspective',
    params: [p('projection', 'enum', false, ['Orthographic', 'Perspective'])],
    examples: ['display projection Orthographic'],
    compile(w) {
      const t = w.get(0, 'a projection').toLowerCase();
      if (t !== 'orthographic' && t !== 'perspective')
        fail('argument', `'${w.get(0, '')}' is not a projection.`, w.atOf(0), 'Orthographic or Perspective.');
      return [{ kind: 'canvas', projection: t }];
    },
  },
  {
    name: 'display depthcue',
    summary: 'Fog on or off.',
    usage: 'display depthcue on|off',
    params: [p('state', 'enum', false, ['on', 'off'])],
    examples: ['display depthcue off'],
    compile(w) {
      const t = w.get(0, 'on or off').toLowerCase();
      if (t !== 'on' && t !== 'off') fail('argument', `'${t}' is not on or off.`, w.atOf(0));
      return [{ kind: 'canvas', fog: t === 'on' }];
    },
  },
  {
    name: 'display resetview',
    summary: 'Reset the camera.',
    usage: 'display resetview',
    params: [],
    examples: ['display resetview'],
    compile: () => [{ kind: 'camera', action: 'reset', animate: true }],
  },
  {
    name: 'display update',
    summary: 'Nothing to do: Mol* draws when the scene changes.',
    usage: 'display update [on|off]',
    params: [p('state', 'enum', true, ['on', 'off', 'ui'])],
    examples: ['display update'],
    compile(w, c) {
      c.report(note('display update is not needed: Mol* draws when the scene changes.', w.at));
      return [];
    },
  },
  {
    name: 'axes location',
    summary: 'Nothing to do with `off`: Mol* shows no axes in the scene.',
    usage: 'axes location off',
    params: [p('location', 'enum', false, ['off', 'origin', 'lowerleft', 'lowerright', 'upperleft', 'upperright'])],
    examples: ['axes location off'],
    compile(w, c) {
      if (w.get(0, 'a location').toLowerCase() !== 'off') {
        fail(
          'unsupported',
          'Axes in the scene are not supported.',
          w.atOf(0),
          'Mol* shows axes in exported images: png … with axes.',
        );
      }
      c.report(note('Mol* shows no axes in the scene: there are none to hide.', w.at));
      return [];
    },
  },
  {
    name: 'rotate',
    summary: 'Rotate the view about a screen axis, in degrees (`rotate y by 90`).',
    usage: 'rotate x|y|z by angle',
    params: [p('axis', 'enum', false, ['x', 'y', 'z']), p('by', 'enum', false, ['by']), p('angle', 'number')],
    examples: ['rotate y by 90'],
    compile(w) {
      const axis = w.get(0, 'an axis').toLowerCase();
      if (axis !== 'x' && axis !== 'y' && axis !== 'z')
        fail('argument', `'${axis}' is not an axis.`, w.atOf(0), 'x, y or z.');
      if (w.get(1, '`by`').toLowerCase() !== 'by')
        fail(
          'unsupported-argument',
          `rotate ${axis} ${w.get(1, '')} is not supported.`,
          w.atOf(1),
          'Write rotate y by 90.',
        );
      return [{ kind: 'camera.turn', axis, angle: w.number(2, 'the angle') }];
    },
  },
  {
    name: 'translate',
    summary: 'Move the view by screen units (`translate by x y z`).',
    usage: 'translate by x y z',
    params: [p('by', 'enum', false, ['by']), p('x', 'number'), p('y', 'number'), p('z', 'number')],
    examples: ['translate by 0.1 0 0'],
    compile(w) {
      if (w.get(0, '`by`').toLowerCase() !== 'by')
        fail('unsupported-argument', 'translate to is not supported.', w.atOf(0), 'Write translate by x y z.');
      const out: CommandBody[] = [];
      (['x', 'y', 'z'] as const).forEach((axis, k) => {
        const d = w.number(k + 1, axis);
        if (d) out.push({ kind: 'camera.move', axis, distance: d * 10 });
      });
      return out;
    },
  },
  {
    name: 'scale',
    summary: 'Zoom the view (`scale by 1.2`: closer).',
    usage: 'scale by factor',
    params: [p('by', 'enum', false, ['by']), p('factor', 'number')],
    examples: ['scale by 1.2'],
    compile(w) {
      if (w.get(0, '`by`').toLowerCase() !== 'by')
        fail('unsupported-argument', 'scale to is not supported.', w.atOf(0), 'Write scale by 1.2.');
      const f = w.number(1, 'the factor');
      if (f <= 0) fail('argument', 'The factor is positive.', w.atOf(1));
      return [{ kind: 'camera.zoom', factor: f }];
    },
  },
  {
    name: 'animate goto',
    summary: 'Show a frame of the trajectory (0-based; start, end).',
    usage: 'animate goto frame|start|end',
    params: [p('frame', 'text')],
    examples: ['animate goto 0', 'animate goto end'],
    compile(w) {
      const t = w.get(0, 'a frame').toLowerCase();
      if (t === 'start') return [{ kind: 'trajectory', frame: 0 }];
      if (t === 'end') return [{ kind: 'trajectory', frame: -1 }];
      return [{ kind: 'trajectory', frame: w.int(0, 'the frame') }];
    },
  },
  ...(['forward', 'stop'] as const).map(
    (play): VmdCommand => ({
      name: `animate ${play}`,
      summary: play === 'forward' ? 'Play the trajectory.' : 'Stop playing the trajectory.',
      usage: `animate ${play}`,
      params: [],
      examples: [`animate ${play}`],
      compile: () => [{ kind: 'trajectory', play }],
    }),
  ),
  {
    name: 'render',
    summary: 'Save an image (PNG; the renderer is Mol*’s).',
    usage: 'render snapshot|TachyonInternal|… file',
    params: [p('renderer', 'enum', false, ['snapshot', 'TachyonInternal']), p('file', 'file', true)],
    examples: ['render snapshot image.png'],
    compile(w, c) {
      const renderer = w.get(0, 'a renderer');
      if (renderer.toLowerCase() !== 'snapshot')
        c.report(note(`${renderer} is ignored: Mol* renders the image.`, w.atOf(0)));
      const file = w.opt(1) ?? 'vmdscene';
      return [{ kind: 'export.image', name: (file.split(/[\\/]/).pop() ?? file).replace(/\.[^.]+$/, '') }];
    },
  },
  {
    name: 'puts',
    summary: 'Print text.',
    usage: 'puts text',
    params: [p('text', 'text')],
    examples: ['puts "loaded"'],
    compile: (w) => [{ kind: 'print', text: w.rest(0) }],
  },
  {
    name: 'undo',
    summary: 'Undo the last step (alone on its line; Mol*’s, VMD has none).',
    usage: 'undo',
    params: [],
    examples: ['undo'],
    compile: () => [{ kind: 'history', action: 'undo' }],
  },
];

/** `measure bond|angle|dihed {i j …} [molid m]`: a read of the geometry of VMD atom indices. */
function geometry(w: VmdWords, c: Compiler, measure: 'distance' | 'angle' | 'dihedral', n: number): Compiled {
  const xs = indices(w.get(0, 'atom indices'), w.atOf(0), n);
  const options = w.options(1);
  const molecule = options.get('molid')?.value;
  for (const [k, v] of options) if (k !== 'molid') c.report(warning('unsupported-argument', `${k} is ignored.`, v.at));
  return [
    { kind: 'read', read: { kind: 'geometry', measure, atoms: xs.map((x) => atomSel(c, x, w.atOf(0), molecule)) } },
  ];
}

/** The commands with a subcommand (`mol addrep`): their first word. */
export const VmdGroups = new Set(VmdCommands.filter((x) => x.name.includes(' ')).map((x) => x.name.split(' ')[0]));

export const VmdCommandsByName: ReadonlyMap<string, VmdCommand> = new Map(
  VmdCommands.flatMap((x) => [x.name, ...(x.aliases ?? [])].map((n) => [n.toLowerCase(), x] as const)),
);
