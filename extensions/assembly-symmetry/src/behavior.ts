/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { AssemblySymmetryProvider, AssemblySymmetryData } from '@molstar/assembly-symmetry-extension/prop';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { AssemblySymmetryClusterColorThemeProvider } from '@molstar/assembly-symmetry-extension/color';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import type { GenericRepresentationRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { AssemblySymmetryControls } from '@molstar/assembly-symmetry-extension/ui';
import { AssemblySymmetry3D, AssemblySymmetryPreset, InitAssemblySymmetry3D } from './transforms.js';

// The transforms, the preset and the config live in `transforms.ts` (no React), re-exported here as before.
export {
  AssemblySymmetry3D,
  AssemblySymmetryConfig,
  AssemblySymmetryPreset,
  getAssemblySymmetryConfig,
  InitAssemblySymmetry3D,
  tryCreateAssemblySymmetry,
} from './transforms.js';

const Tag = AssemblySymmetryData.Tag;

export const AssemblySymmetry = PluginBehavior.create<{ autoAttach: boolean }>({
  name: 'assembly-symmetry-prop',
  category: 'custom-props',
  display: {
    name: 'Assembly Symmetry',
    description: 'Assembly Symmetry data provided by RCSB PDB (calculated with BioJava) or by PDBe.',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
    private provider = AssemblySymmetryProvider;
    private unregisterEntry: (() => void) | undefined;

    register(): void {
      const entry: PluginRegistryEntry = {
        structure: {
          themes: { color: [AssemblySymmetryClusterColorThemeProvider] },
          presets: { representation: [AssemblySymmetryPreset] },
        },
        actions: [InitAssemblySymmetry3D],
      };
      this.unregisterEntry = this.ctx.register(entry);
      this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

      this.ctx.genericRepresentationControls.set(Tag.Representation, (selection) => {
        const refs: GenericRepresentationRef[] = [];
        selection.structures.forEach((structure) => {
          const symmRepr = structure.genericRepresentations?.filter(
            (r) => r.cell.transform.transformer.id === AssemblySymmetry3D.id,
          )[0];
          if (symmRepr) refs.push(symmRepr);
        });
        return [refs, 'Symmetries'];
      });
      this.ctx.customStructureControls.set(Tag.Representation, AssemblySymmetryControls as any);
    }

    update(p: { autoAttach: boolean }) {
      const updated = this.params.autoAttach !== p.autoAttach;
      this.params.autoAttach = p.autoAttach;
      this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
      return updated;
    }

    unregister() {
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
      this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

      this.ctx.genericRepresentationControls.delete(Tag.Representation);
      this.ctx.customStructureControls.delete(Tag.Representation);
    }
  },
  params: () => ({
    autoAttach: PD.Boolean(false),
    serverUrl: PD.Text(AssemblySymmetryData.DefaultServerUrl),
  }),
});

//
