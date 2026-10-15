/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginUIComponent } from '@molstar/plugin-ui/base';
import { CustomStructureControls } from '@molstar/plugin-ui/controls';
import { Icon, BuildSvg } from '@molstar/plugin-ui/controls/icons';
import { StructureComponentControls } from '@molstar/plugin-ui/structure/components';
import { StructureMeasurementsControls } from '@molstar/plugin-ui/structure/measurements';
import { StructureSourceControls } from '@molstar/plugin-ui/structure/source';
import { VolumeStreamingControls, VolumeSourceControls } from '@molstar/plugin-ui/structure/volume';
import { ParticleSourceControls } from '@molstar/plugin-ui/structure/particles';
import { PluginConfig } from '@molstar/plugin/config';
import { StructureSuperpositionControls } from '@molstar/plugin-ui/structure/superposition';
import { StructureQuickStylesControls } from '@molstar/plugin-ui/structure/quick-styles';
import { StructureProceduralAnimationControls } from '@molstar/plugin-ui/structure/procedural-animation';

export class DefaultStructureTools extends PluginUIComponent {
  render() {
    return (
      <>
        <div className="msp-section-header">
          <Icon svg={BuildSvg} />
          Structure Tools
        </div>

        <StructureSourceControls />
        <StructureMeasurementsControls />
        <StructureSuperpositionControls />
        <StructureQuickStylesControls />
        <StructureProceduralAnimationControls />
        <StructureComponentControls />
        {this.plugin.config.get(PluginConfig.VolumeStreaming.Enabled) && <VolumeStreamingControls />}
        <VolumeSourceControls />
        <ParticleSourceControls />

        <CustomStructureControls />
      </>
    );
  }
}
