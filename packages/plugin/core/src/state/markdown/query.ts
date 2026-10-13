/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection } from '@molstar/core/state';
import { QueryContext, type QueryFn, StructureElement, StructureSelection } from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import type { MarkdownExtension } from '../manager/markdown-extensions.js';
import { PluginStateObject } from '../objects.js';
import { parseArray } from './helpers.js';

export const QueryMarkdownExtension: MarkdownExtension = {
  name: 'query',
  execute: ({ event, args, manager }) => {
    const expression = args['query'];
    if (!expression?.length) return;

    // supported languages: mol-script, pymol, vmd, jmol
    const language = args['lang'] || 'mol-script';
    // supported actions: highlight, focus
    const action = parseArray(args['action'] || 'highlight');
    const focusRadius = parseFloat(args['focus-radius'] || '3');

    if (event === 'mouse-leave') {
      if (action.includes('highlight')) {
        manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
      }
      return;
    }

    let query: QueryFn<StructureSelection>;
    try {
      query = Script.toQuery({
        language: language as Script.Language,
        expression,
      });
    } catch (e) {
      console.warn(`Failed to parse query '${expression}' (${language})`, e);
      return;
    }

    const structures = manager.plugin.state.data.selectQ((q) => q.rootsOfType(PluginStateObject.Molecule.Structure));

    if (event === 'mouse-enter') {
      if (!action.includes('focus')) {
        return;
      }
      manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
      for (const structure of structures) {
        if (!structure.obj?.data) continue;
        const selection = query(new QueryContext(structure.obj.data));
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        manager.plugin.managers.interactivity.lociHighlights.highlight(
          {
            loci,
          },
          false,
        );
      }
    }

    if (event === 'click') {
      if (!action.includes('focus')) {
        return;
      }
      const decorated = structures.map((s) =>
        StateSelection.getDecorated<PluginStateObject.Molecule.Structure>(manager.plugin.state.data, s.transform.ref),
      );
      const spheres = decorated
        .map((s) => {
          if (!s.obj?.data) return undefined;
          const selection = query(new QueryContext(s.obj.data));
          if (StructureSelection.isEmpty(selection)) return;

          const loci = StructureSelection.toLociWithSourceUnits(selection);
          return StructureElement.Loci.getBoundary(loci).sphere;
        })
        .filter((s) => !!s);

      if (spheres.length) {
        manager.plugin.managers.camera.focusSpheres(spheres, (s) => s, { extraRadius: focusRadius });
      }
    }
  },
};
