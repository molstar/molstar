/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useContext } from 'react';
import ReactMarkdown, { Components } from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { PluginReactContext } from '../base';
import { PluginUIContext } from '../context';
import { PluginContext } from '../../mol-plugin/context';
import { MarkdownRenderer, parseMarkdownCommandArgs } from '../../mol-plugin/util/markdown-commands';
import { ColorLists } from '../../mol-util/color/lists';
import { getColorGradient } from '../../mol-util/color/utils';

export function Markdown({ children, components }: { children?: string, components?: Components }) {
    return <div className='msp-markdown'>
        <ReactMarkdown
            skipHtml
            components={{ a: MarkdownAnchor, img: MarkdownImg, ...components }}
            remarkPlugins={[remarkGfm]}
        >
            {children}
        </ReactMarkdown>
    </div>;
}

export function MarkdownImg({ src, children, element }: { src?: string, children?: any, element?: any }) {
    const plugin: PluginUIContext | undefined = useContext(PluginReactContext);

    if (!src) return element;

    if (src[0] === '!') {
        warnMissingPlugin(plugin);
        const args = parseMarkdownCommandArgs(src.substring(1));
        const result = plugin?.managers.markdownCommands.render(args, DefaultRenderers);
        return result ?? element;
    }

    return children;
}

export const DefaultRenderers: MarkdownRenderer[] = [
    {
        name: 'color-swatch',
        reactRenderFn: (args) => {
            const color = args['color-swatch'];
            if (!color) return null;
            return <span style={{ display: 'inline-block', width: '0.75em', height: '0.75em', backgroundColor: color, borderRadius: '25%' }}/>;
        }
    },
     {
        name: 'color-palette',
        reactRenderFn: (args) => {
            const name = args['color-palette'];
            const minWidth = args['color-palette-width'] ?? '150px';
            const height = args['color-palette-height'] ?? '0.5em';
            if (!name) return null;

            const list = ColorLists[name.toLowerCase() as keyof typeof ColorLists];
            if (!list) {
                console.warn(`Color palette '${name}' not found.`);
                return null;
            }

            return <span style={{
                display: 'inline-block',
                minWidth,
                height,
                background: getColorGradient(list.list),
                borderRadius: '2px'
            }} />;
        }
    }
];
export function MarkdownAnchor({ href, children, element }: { href?: string, children?: any, element?: any }) {
    const plugin: PluginUIContext | undefined = useContext(PluginReactContext);

    if (!href) return element;

    if (href[0] === '#') {
        warnMissingPlugin(plugin);
        return <a href='#' onClick={(e) => {
            e.preventDefault();
            plugin?.managers.snapshot.applyKey(href.substring(1));
        }}>{children}</a>;
    } else if (href[0] === '!') {
        const args = parseMarkdownCommandArgs(href.substring(1));
        warnMissingPlugin(plugin);
        return <a href='!command'
            onClick={(e) => {
                e.preventDefault();
                plugin?.managers.markdownCommands.execute('click', args);
            }}
            onMouseEnter={() => plugin?.managers.markdownCommands.execute('mouse-enter', args)}
            onMouseLeave={() => plugin?.managers.markdownCommands.execute('mouse-leave', args)}
        >
            {children}
        </a>;
    } else if (href) {
        return <a href={href} target='_blank' rel='noopener noreferrer'>{children}</a>;
    }


    return children;
}

function warnMissingPlugin(plugin: PluginContext | undefined) {
    if (plugin) return;
    console.warn('Markdown component requires a PluginReactContext to be set.');
}