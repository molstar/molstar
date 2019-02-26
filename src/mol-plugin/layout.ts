/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginComponent } from './component';
import { PluginContext } from './context';
import { PluginCommands } from './command';

// TODO: support collapsed state control orientation
export const PluginLayoutStateParams = {
    isExpanded: PD.Boolean(false),
    showControls: PD.Boolean(true)
}

export type PluginLayoutStateProps = PD.Values<typeof PluginLayoutStateParams>

interface RootState {
    top: string | null,
    bottom: string | null,
    left: string | null,
    right: string | null,

    width: string | null,
    height: string | null,
    maxWidth: string | null,
    maxHeight: string | null,
    margin: string | null,
    marginLeft: string | null,
    marginRight: string | null,
    marginTop: string | null,
    marginBottom: string | null,

    scrollTop: number,
    scrollLeft: number,
    position: string | null,
    overflow: string | null,
    viewports: HTMLElement[],
    zindex: string | null
}

export class PluginLayout extends PluginComponent<PluginLayoutStateProps> {
    readonly events = {
        updated: this.ev()
    }

    private updateProps(state: Partial<PluginLayoutStateProps>) {
        let prevExpanded = !!this.state.isExpanded;
        this.updateState(state);
        if (this.root && typeof state.isExpanded === 'boolean' && state.isExpanded !== prevExpanded) this.handleExpand();

        this.events.updated.next();
    }

    private root: HTMLElement;
    private rootState: RootState | undefined = void 0;
    private expandedViewport: HTMLMetaElement;

    setProps(props: PluginLayoutStateProps) {
        this.updateState(props);
    }

    setRoot(root: HTMLElement) {
        this.root = root;
        if (this.state.isExpanded) this.handleExpand();
    }

    private getScrollElement() {
        if ((document as any).scrollingElement) return (document as any).scrollingElement;
        if (document.documentElement) return document.documentElement;
        return document.body;
    }

    private handleExpand() {
        try {
            let body = document.getElementsByTagName('body')[0];
            let head = document.getElementsByTagName('head')[0];

            if (!body || !head) return;

            if (this.state.isExpanded) {
                let children = head.children;
                let hasExp = false;
                let viewports: HTMLElement[] = [];
                for (let i = 0; i < children.length; i++) {
                    if (children[i] === this.expandedViewport) {
                        hasExp = true;
                    } else if (((children[i] as any).name || '').toLowerCase() === 'viewport') {
                        viewports.push(children[i] as any);
                    }
                }

                for (let v of viewports) {
                    head.removeChild(v);
                }

                if (!hasExp) head.appendChild(this.expandedViewport);


                let s = body.style;

                let doc = this.getScrollElement();
                let scrollLeft = doc.scrollLeft;
                let scrollTop = doc.scrollTop;

                this.rootState = {
                    top: s.top, bottom: s.bottom, right: s.right, left: s.left, scrollTop, scrollLeft, position: s.position, overflow: s.overflow, viewports, zindex: this.root.style.zIndex,
                    width: s.width, height: s.height,
                    maxWidth: s.maxWidth, maxHeight: s.maxHeight,
                    margin: s.margin, marginLeft: s.marginLeft, marginRight: s.marginRight, marginTop: s.marginTop, marginBottom: s.marginBottom
                };

                s.overflow = 'hidden';
                s.position = 'fixed';
                s.top = '0';
                s.bottom = '0';
                s.right = '0';
                s.left = '0';

                s.width = '100%';
                s.height = '100%';
                s.maxWidth = '100%';
                s.maxHeight = '100%';
                s.margin = '0';
                s.marginLeft = '0';
                s.marginRight = '0';
                s.marginTop = '0';
                s.marginBottom = '0';

                // TODO: setting this breaks viewport controls for some reason. Is there a fix?
                // this.root.style.zIndex = '100000';
            } else {
                let children = head.children;
                for (let i = 0; i < children.length; i++) {
                    if (children[i] === this.expandedViewport) {
                        head.removeChild(this.expandedViewport);
                        break;
                    }
                }

                if (this.rootState) {
                    let s = body.style, t = this.rootState;
                    for (let v of t.viewports) {
                        head.appendChild(v);
                    }
                    s.top = t.top;
                    s.bottom = t.bottom;
                    s.left = t.left;
                    s.right = t.right;

                    s.width = t.width;
                    s.height = t.height;
                    s.maxWidth = t.maxWidth;
                    s.maxHeight = t.maxHeight;
                    s.margin = t.margin;
                    s.marginLeft = t.marginLeft;
                    s.marginRight = t.marginRight;
                    s.marginTop = t.marginTop;
                    s.marginBottom = t.marginBottom;

                    s.position = t.position;
                    s.overflow = t.overflow;
                    let doc = this.getScrollElement();
                    doc.scrollTop = t.scrollTop;
                    doc.scrollLeft = t.scrollLeft;
                    this.rootState = void 0;
                    this.root.style.zIndex = t.zindex;
                }
            }
        } catch (e) {
            this.context.log.error('Layout change error, you might have to reload the page.');
            console.log('Layout change error, you might have to reload the page.', e);
        }
    }

    constructor(private context: PluginContext) {
        super({ ...PD.getDefaultValues(PluginLayoutStateParams), ...(context.spec.layout && context.spec.layout.initial) });

        PluginCommands.Layout.Update.subscribe(context, e => this.updateProps(e.state));

        // TODO how best make sure it runs on node.js as well as in the browser?
        if (typeof document !== 'undefined') {
            // <meta name='viewport' content='width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=0' />
            this.expandedViewport = document.createElement('meta') as any;
            this.expandedViewport.name = 'viewport';
            this.expandedViewport.content = 'width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=0';
        }
    }
}