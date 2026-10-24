/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ColorSpec, Command, RepType, Sel } from './ir.js';

// One line per command, before it runs: what the menu's Run row says (_spec §9), in Mol*'s words.

const sel = (s: Sel) => (s.text === 'all' ? 'everything' : `'${s.text}'`);

const reps = (r: readonly RepType[] | 'all') => (r === 'all' ? 'everything' : r.join(', '));

const hex = (rgb: number) => `#${rgb.toString(16).padStart(6, '0')}`;

const color = (c: ColorSpec) => ('rgb' in c ? hex(c.rgb) : `by ${c.theme}`);

export function describeCommand(c: Command): string {
  switch (c.kind) {
    case 'load':
      return 'entry' in c.from
        ? `Load ${c.from.entry.id.toUpperCase()} from ${c.from.entry.db === 'pdb' ? 'the PDB' : 'AlphaFold DB'} as ${c.name}`
        : 'url' in c.from
          ? `Load ${c.from.url} as ${c.name}`
          : `Load ${c.from.file} as ${c.name}`;
    case 'object.remove':
      return c.objects === 'all' ? 'Remove everything' : `Remove ${c.objects.join(', ')}`;
    case 'object.visibility':
      return `${c.visible ? 'Show' : 'Hide'} ${c.objects === 'all' ? 'every structure' : c.objects.join(', ')}`;
    case 'object.top':
      return `Make ${c.object} the top molecule`;
    case 'name.set':
      return `Name ${sel(c.sel)} ${c.name}`;
    case 'name.remove':
      return c.names === 'all' ? 'Forget every named selection' : `Forget ${c.names.join(', ')}`;
    case 'select':
      return `Select ${sel(c.sel)}${c.mode === 'set' ? '' : ` (${c.mode})`}`;
    case 'select.clear':
      return 'Clear the selection';
    case 'display':
      return c.mode === 'hide'
        ? `Hide ${reps(c.reps)} for ${sel(c.sel)}`
        : `Show ${sel(c.sel)} as ${c.mode === 'only' ? 'only ' : ''}${reps(c.reps)}`;
    case 'cartoon':
      return `Draw the cartoon of ${sel(c.sel)} as ${c.style}`;
    case 'rep.add':
      return `Add a ${c.rep.style} rep of '${c.rep.selection}' to ${c.object}, colored ${color(c.rep.color)}`;
    case 'rep.modify':
      return `Change rep ${c.index} of ${c.object}`;
    case 'rep.remove':
      return `Remove rep ${c.index} of ${c.object}`;
    case 'rep.visibility':
      return `${c.visible ? 'Show' : 'Hide'} rep ${c.index} of ${c.object}`;
    case 'rep.defaults':
      return `Set the next rep’s ${Object.keys(c.patch).join(', ')}`;
    case 'color':
      return `Color ${sel(c.sel)} ${color(c.color)}`;
    case 'transparency':
      return `Make ${sel(c.sel)} ${Math.round(c.alpha * 100)}% transparent`;
    case 'size':
      return `Size ${sel(c.sel)}${c.reps ? ` (${c.reps.join(', ')})` : ''} ×${c.factor}`;
    case 'hide-contents':
      return `Leave ${sel(c.sel)} out of the scene`;
    case 'hydrogens':
      return `${c.visible ? 'Show' : 'Hide'} hydrogens`;
    case 'camera':
      return c.action === 'reset'
        ? 'Reset the camera'
        : `${c.action === 'frame' ? 'Frame' : c.action === 'orient' ? 'Orient on' : 'Center on'} ${c.sel ? sel(c.sel) : 'everything'}`;
    case 'camera.zoom':
      return `Zoom by ${c.factor}`;
    case 'camera.turn':
      return `Turn the view ${c.angle}° about ${c.axis}`;
    case 'camera.move':
      return `Move the view ${c.distance} Å along ${c.axis}`;
    case 'measure':
      return `Measure the ${c.measure} of ${c.atoms.map(sel).join(', ')}`;
    case 'measure.clear':
      return `Remove the ${c.measure} measurements`;
    case 'label':
      return `Label ${sel(c.sel)}`;
    case 'label.clear':
      return c.sel ? `Remove the labels of ${sel(c.sel)}` : 'Remove the labels';
    case 'superpose':
      return `Superpose ${sel(c.mobile)} on ${sel(c.target)} by ${c.method}`;
    case 'canvas': {
      const parts = [
        c.background !== undefined ? `background ${hex(c.background)}` : '',
        c.projection ? `${c.projection} projection` : '',
        c.fov !== undefined ? `field of view ${c.fov}°` : '',
        c.fog !== undefined ? `fog ${c.fog ? 'on' : 'off'}` : '',
      ].filter(Boolean);
      return `Set the ${parts.join(', ')}`;
    }
    case 'trajectory':
      return c.play
        ? c.play === 'stop'
          ? 'Stop the trajectory'
          : 'Play the trajectory'
        : `Show frame ${c.frame === -1 ? 'last' : c.frame}`;
    case 'export.image':
      return `Save an image${c.name ? ` ${c.name}.png` : ''}`;
    case 'history':
      return c.action === 'undo' ? 'Undo' : 'Redo';
    case 'read':
      return c.read.kind === 'count'
        ? `Count the atoms of ${sel(c.read.sel)}`
        : c.read.kind === 'geometry'
          ? `Measure the ${c.read.measure}`
          : c.read.kind === 'values'
            ? `Read ${c.read.properties.join(', ')}`
            : `List the ${c.read.kind}`;
    case 'print':
      return c.text;
    case 'require':
      return '';
  }
}
