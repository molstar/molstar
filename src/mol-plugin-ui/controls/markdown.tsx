/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useContext } from 'react';
import { PluginReactContext } from '../base';
import ReactMarkdown, { Components } from 'react-markdown';
import { PluginUIContext } from '../context';
import { PluginContext } from '../../mol-plugin/context';
import { parseMarkdownCommandArgs } from '../../mol-plugin/util/markdown-commands';

export function Markdown({ children, components }: { children?: string, components?: Components }) {
    return <ReactMarkdown skipHtml components={{ a: MarkdownAnchor, img: MarkdownImg, ...components }}>
        {children}
    </ReactMarkdown>;
}

export function MarkdownImg({ src, children, element }: { src?: string, children?: any, element?: any }) {
    // const plugin: PluginUIContext | undefined = useContext(PluginReactContext);

    if (!src) return element;

    if (src[0] === '!') {
        // TODO
        return children;
    }

    return children;
}

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
    console.warn('Markdown component requires a PluginReactContext to be set.');
}