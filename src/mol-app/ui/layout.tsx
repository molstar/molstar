/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import * as React from 'react'
import { LayoutController, LayoutTarget, LayoutRegion, CollapsedControlsLayout } from '../controller/layout';
import { View } from './view';

export class Layout extends View<LayoutController, { }, { }> {

    private renderTarget(target: LayoutTarget) {
        const statics: any[] = [];
        const scrollable: any[] = [];

        for (let c of target.components) {
            if (c.isStatic) statics.push(<c.view key={c.key} controller={c.controller} />);
            else scrollable.push(<c.view key={c.key} controller={c.controller} />);
        }

        return <div key={target.cssClass} className={'molstar-layout-region molstar-layout-' + target.cssClass}>
            { statics.length ? <div className='molstar-layout-static'>{statics}</div> : void 0 }
            { scrollable.length ? <div className='molstar-layout-scrollable'>{scrollable}</div> : void 0 }
        </div>;
    }

    private updateTarget(name: string, regionType: LayoutRegion, layout: { regions: any[], layoutClass: string }) {
        const state = this.controller.latestState;
        const regionStates = state.regionStates;
        const region = this.controller.targets[regionType];
        let show: boolean;

        if (state.hideControls) {
            show = regionStates !== void 0 && regionStates[regionType] === 'Sticky' && region.components.length > 0;
        } else if (regionStates && regionStates[regionType] === 'Hidden') {
            show = false;
        } else {
            show = region.components.length > 0;
        }

        if (show) {
            layout.regions.push(this.renderTarget(region));
        } else {
            layout.layoutClass += ' molstar-layout-hide-' + name;
        }
    }

    render() {
        let layoutClass = '';

        const state = this.controller.latestState;
        let layoutType: string;

        if (state.isExpanded) {
            layoutType = 'molstar-layout-expanded';
        } else {
            layoutType = 'molstar-layout-standard ';
            switch (state.collapsedControlsLayout) {
                case CollapsedControlsLayout.Outside: layoutType += 'molstar-layout-standard-outside'; break;
                case CollapsedControlsLayout.Landscape: layoutType += 'molstar-layout-standard-landscape'; break;
                case CollapsedControlsLayout.Portrait: layoutType += 'molstar-layout-standard-portrait'; break;
                default: layoutType += 'molstar-layout-standard-outside'; break;
            }
        }

        const targets = this.controller.targets;
        const regions = [this.renderTarget(targets[LayoutRegion.Main])];

        const layout = { regions, layoutClass };
        this.updateTarget('top', LayoutRegion.Top, layout);
        this.updateTarget('right', LayoutRegion.Right, layout);
        this.updateTarget('bottom', LayoutRegion.Bottom, layout);
        this.updateTarget('left', LayoutRegion.Left, layout);
        layoutClass = layout.layoutClass;

        let root = targets[LayoutRegion.Root].components.map(c => <c.view key={c.key} controller={c.controller} />);

        return <div className='molstar-plugin'>
            <div className={'molstar-plugin-content ' + layoutType}>
                <div className={layoutClass}>
                    {regions}
                    {root}
                </div>
            </div>
        </div>;
    }
}