/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginUIComponent } from '../../mol-plugin-ui/base';
import * as React from 'react';

export class CustomToastMessage extends PluginUIComponent {
    render() {
        return <>
            Custom <i>Toast</i> content. No timeout.
        </>;
    }
}