/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginUIComponent } from 'mol-plugin/ui/base';
import * as React from 'react';
import { TransformUpdaterControl } from 'mol-plugin/ui/state/update-transform';

export class BasicWrapperControls extends PluginUIComponent {

    render() {
        return <div style={{ overflowY: 'auto', display: 'block', height: '100%' }}>
            <TransformUpdaterControl nodeRef='asm' />
            <TransformUpdaterControl nodeRef='seq-visual' header={{ name: 'Sequence Visual' }} />
            <TransformUpdaterControl nodeRef='het-visual' header={{ name: 'HET Visual' }} />
            <TransformUpdaterControl nodeRef='water-visual' header={{ name: 'Water Visual' }} initiallyCollapsed={true} />
            <TransformUpdaterControl nodeRef='ihm-visual' header={{ name: 'I/HM Visual' }} initiallyCollapsed={true} />
        </div>;
    }
}