/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useContext, useEffect, useRef, useState } from 'react';
import ReactMarkdown, { Components } from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { PluginReactContext } from '../base';
import { PluginUIContext } from '../context';
import { PluginContext } from '../../mol-plugin/context';
import { MarkdownExtension } from '../../mol-plugin-state/manager/markdown-extensions';
import { ColorLists } from '../../mol-util/color/lists';
import { getColorGradient, getColorGradientBanded, parseColorList } from '../../mol-util/color/utils';
import { useBehavior } from '../hooks/use-behavior';

export function Markdown({ children, components }: { children?: string, components?: Components }) {
    return <div className='msp-markdown'>
        <MarkdownAudioPlayer />
        <ReactMarkdown
            skipHtml
            components={{ a: MarkdownAnchor, img: MarkdownImg, ...components }}
            remarkPlugins={[remarkGfm]}
        >
            {children}
        </ReactMarkdown>
    </div>;
}

export function MarkdownAudioPlayer() {
    const parent = useRef<HTMLDivElement>(null);
    const plugin: PluginUIContext | undefined = useContext(PluginReactContext);
    const audio = useBehavior(plugin?.managers.markdownExtensions.state.audioPlayer);

    useEffect(() => {
        if (!parent.current) return;
        parent.current.appendChild(audio!);
        return () => { audio?.remove(); };
    }, [audio]);
    if (!audio) return null;

    return <div className='msp-markdown-audio-player' ref={parent} />;
}

export function MarkdownImg({ src, element, alt }: { src?: string, element?: any, alt?: string }) {
    const plugin: PluginUIContext | undefined = useContext(PluginReactContext);

    if (!src) return element;

    warnMissingPlugin(plugin);
    const args = plugin?.managers.markdownExtensions.parseArgs(src);
    if (args) {
        const result = plugin?.managers.markdownExtensions.tryRender(args, DefaultMarkdownExtensionRenderers);
        return result ?? element;
    } else {
        const data = plugin?.managers.markdownExtensions.tryResolveUri(src);
        if (typeof (data as Promise<string>)?.then === 'function') {
            return <LazyStaticImg alt={alt} data={data as Promise<string>} />;
        } else if (typeof data === 'string' && data) {
            return <img src={data} alt={alt} />;
        }
    }

    return <img src={src} alt={alt} />;
}

function LazyStaticImg({ alt, data }: { alt?: string, data: Promise<string> }) {
    const [src, setSrc] = useState<string | undefined>(undefined);
    useEffect(() => {
        let mounted = true;
        data.then(d => {
            if (mounted) setSrc(d);
        }).catch(e => {
            console.error('Failed to load static image', e);
            if (mounted) setSrc(undefined);
        });
        return () => { mounted = false; };
    }, [data]);
    if (!src) return null;
    return <img src={src} alt={alt} />;
}

export const DefaultMarkdownExtensionRenderers: MarkdownExtension[] = [
    {
        name: 'color-swatch',
        reactRenderFn: ({ args }) => {
            const color = args['color-swatch'];
            if (!color) return null;
            return <span style={{ display: 'inline-block', width: '0.75em', height: '0.75em', backgroundColor: color, borderRadius: '25%' }}/>;
        }
    },
     {
        name: 'color-palette',
        reactRenderFn: ({ args }) => {
            const name = args['color-palette-name'];
            const colors = args['color-palette-colors'];
            const minWidth = args['color-palette-width'] ?? '150px';
            const height = args['color-palette-height'] ?? '0.5em';
            const discrete = 'color-palette-discrete' in args;
            if (!name && !colors) return null;

            const list = colors
                ? parseColorList(colors)
                : ColorLists[name.toLowerCase() as keyof typeof ColorLists]?.list;

            if (!list?.length) {
                console.warn(`Color palette could not be resolved.`, args);
                return null;
            }

            return <span style={{
                display: 'inline-block',
                minWidth,
                height,
                background: (discrete ? getColorGradientBanded : getColorGradient)(list),
                borderRadius: '2px'
            }} />;
        }
    }
];

export function MarkdownAnchor({ href, children, element }: { href?: string, children?: any, element?: any }) {
    const plugin: PluginUIContext | undefined = useContext(PluginReactContext);

    if (!href) return element;

    warnMissingPlugin(plugin);
    const args = plugin?.managers.markdownExtensions.parseArgs(href);

    if (args) {
        return <a href='#'
            onClick={(e) => {
                e.preventDefault();
                plugin?.managers.markdownExtensions.tryExecute('click', args);
            }}
            onMouseEnter={() => plugin?.managers.markdownExtensions.tryExecute('mouse-enter', args)}
            onMouseLeave={() => plugin?.managers.markdownExtensions.tryExecute('mouse-leave', args)}
        >
            {children}
        </a>;
    } else if (href[0] === '#') {
        warnMissingPlugin(plugin);
        return <a href='#' onClick={(e) => {
            e.preventDefault();
            plugin?.managers.snapshot.applyKey(href.substring(1));
        }}>{children}</a>;
    } else if (href) {
        return <a href={href} target='_blank' rel='noopener noreferrer'>{children}⤴</a>;
    }

    return children;
}

function warnMissingPlugin(plugin: PluginContext | undefined) {
    if (plugin) return;
    console.warn('Markdown component requires a PluginReactContext to be set.');
}