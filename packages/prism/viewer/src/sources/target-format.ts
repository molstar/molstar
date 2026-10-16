/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement, StructureProperties } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { CommandFormat, FormatCompletion } from '@molstar/prism-ui/core';
import type { ViewerCommandContext } from '../api/managers.js';
import type { QueryLanguages } from '../targets/languages.js';
import { KEYWORDS, PERSON_KEYWORDS, parseTarget } from '../targets/index.js';

// The `target` field of commands (targets.md §1; parity.md §3, the focus picker's entries): checked as it is typed
// (parsed, not evaluated), and completed with the keywords, the person's selection and focus, the loaded chains and
// ligands. What the command menu offers is what an operation takes: the same strings.

const KeywordHelp: Readonly<Record<string, string>> = {
  all: 'Everything',
  polymer: 'Proteins and nucleic acids',
  protein: 'Proteins',
  nucleic: 'Nucleic acids',
  branched: 'Carbohydrates',
  ligand: 'Ligands',
  ion: 'Ions',
  water: 'Waters',
  coarse: 'Coarse-grained parts',
  selection: 'What is selected',
  focus: 'What is focused',
};

export function createTargetFormat(
  plugin: PluginContext,
  languages: QueryLanguages,
): CommandFormat<ViewerCommandContext> {
  return {
    hint: () => 'ligand, A:120-130, resn HEM, selection, …',
    check(value) {
      const text = value.trim();
      if (!text)
        return { ok: false, message: 'Which part?', hint: 'A keyword (ligand), an address (A:120) or a query.' };
      const parsed = parseTarget(text, { languages });
      return parsed.ok ? { ok: true, value: text } : { ok: false, message: parsed.message, hint: parsed.hint };
    },
    complete(text) {
      const q = text.trim().toLowerCase();
      const out: FormatCompletion[] = [];
      const add = (value: string, label: string, detail?: string) => {
        if (!q || value.toLowerCase().startsWith(q) || label.toLowerCase().includes(q)) {
          out.push({ value, label, ...(detail ? { detail } : {}) });
        }
      };
      for (const k of [...KEYWORDS, ...PERSON_KEYWORDS]) add(k, k, KeywordHelp[k]);
      const seen = new Set<string>();
      for (const s of plugin.managers.structure.hierarchy.current.structures) {
        const data = s.cell.obj?.data;
        if (!data) continue;
        const loc = StructureElement.Location.create(data);
        for (const unit of data.units) {
          loc.unit = unit;
          loc.element = unit.elements[0];
          const chain = StructureProperties.chain.auth_asym_id(loc);
          if (!seen.has(`chain:${chain}`)) {
            seen.add(`chain:${chain}`);
            add(chain, `Chain ${chain}`, StructureProperties.entity.pdbx_description(loc).join(', '));
          }
          // A ligand: a non-polymer unit's first residue name.
          if (StructureProperties.entity.type(loc) === 'non-polymer') {
            const comp = StructureProperties.atom.label_comp_id(loc);
            if (comp !== 'HOH' && !seen.has(`resn:${comp}`)) {
              seen.add(`resn:${comp}`);
              add(`resn ${comp}`, comp, StructureProperties.entity.pdbx_description(loc).join(', ') || 'Ligand');
            }
          }
        }
      }
      return out.slice(0, 40);
    },
  };
}
