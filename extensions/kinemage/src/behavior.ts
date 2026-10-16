/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { KinemageDataProvider, KinemageData } from '@molstar/kinemage-extension/prop';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { DefaultQueryRuntimeTable } from '@molstar/model/script/runtime/query/compiler';
import { KINFormatProvider, KinemageDragAndDropHandler } from '@molstar/kinemage-extension/model';
import { KinemageControls } from '@molstar/kinemage-extension/ui';

// The v5 plugin behavior: the data model and the .kin format are in model.ts (shared with the Prism extension), this
// adds the right-panel controls.

export * from '@molstar/kinemage-extension/model';

const Tag = KinemageData.Tag;

export const KinemageExtension = PluginBehavior.create<{ autoAttach: boolean }>({
  name: 'kinemage-data-prop',
  category: 'custom-props',
  display: {
    name: 'Kinemage data',
    description: 'Data loaded from Kinemage.',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
    private provider = KinemageDataProvider;
    private unregisterEntry: (() => void) | undefined;

    register(): void {
      DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

      this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

      // Register right-panel controls for Kinemage (show in the right-hand inspector)
      this.ctx.customStructureControls.set(Tag.Representation, KinemageControls as any);
      // Some app hosts expose a global customControls registry; register there too so the card is visible
      // even when no structure is loaded. Use `any` guards to avoid type errors if customControls isn't present.
      if ((this.ctx as any).customControls && typeof (this.ctx as any).customControls.set === 'function') {
        (this.ctx as any).customControls.set('kinemage', KinemageControls as any);
      }

      // Register the .kin file handler so opening/dropping .kin is supported via the data formats system
      const entry: PluginRegistryEntry = {
        formats: [KINFormatProvider],
        dragAndDrop: [KinemageDragAndDropHandler],
      };
      this.unregisterEntry = this.ctx.register(entry);
    }

    update(p: { autoAttach: boolean }) {
      const updated = this.params.autoAttach !== p.autoAttach;
      this.params.autoAttach = p.autoAttach;
      this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
      return updated;
    }

    unregister() {
      DefaultQueryRuntimeTable.removeCustomProp(this.provider.descriptor);

      this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

      this.ctx.genericRepresentationControls.delete(Tag.Representation);

      // Unregister the .kin data format provider and drag-and-drop handler
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;

      // Remove right-panel controls
      try {
        this.ctx.customStructureControls.delete(Tag.Representation);
      } catch {}
      if ((this.ctx as any).customControls && typeof (this.ctx as any).customControls.delete === 'function') {
        try {
          (this.ctx as any).customControls.delete('kinemage');
        } catch {}
      }
    }
  },
  params: () => ({
    autoAttach: PD.Boolean(false),
  }),
});
