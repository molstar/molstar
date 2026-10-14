/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The types of a Seq* story file. The stories are React (CSF-shaped: a default-exported meta and named story objects,
// each with a `render` returning a React element), shown in Mol*'s web-components Storybook: `src/preview.ts` is the
// decorator that mounts their element, and these types stand in for `@storybook/react-vite`'s `Meta` and `StoryObj`
// (the parts the stories use). No React runtime here.

import type { ReactElement } from 'react';

/** A control of the Controls panel (`argTypes` entry): a type, its options, a description. */
export interface StoryArgType {
  control?: string | { type: string; min?: number; max?: number; step?: number };
  options?: readonly unknown[];
  description?: string;
  table?: Record<string, unknown>;
}

/** The `argTypes` of args `A`: a control for some of its keys. */
export type StoryArgTypes<A> = { readonly [K in keyof A]?: StoryArgType };

/** What a story's `render` gets besides its args: the toolbar's globals, the story's id, title and name. */
export interface StoryContext<A = Record<string, unknown>> {
  args: A;
  globals: Record<string, any>;
  id: string;
  title: string;
  name: string;
  parameters: Record<string, any>;
}

/** The default export of a story file: the sidebar title, the parameters, and the args and controls its stories share. */
export interface Meta<A = {}> {
  title: string;
  parameters?: Record<string, any>;
  args?: Partial<A>;
  argTypes?: StoryArgTypes<A>;
  tags?: string[];
}

/** The args of a story: the ones its type names, or none when it is typed by its file's meta (`typeof meta`). */
type ArgsOf<T> = T extends { title: string } ? {} : T;

/** A named export of a story file. */
export interface StoryObj<T = {}> {
  name?: string;
  args?: Partial<ArgsOf<T>>;
  argTypes?: StoryArgTypes<ArgsOf<T>>;
  parameters?: Record<string, any>;
  tags?: string[];
  render: (args: ArgsOf<T>, context: StoryContext<ArgsOf<T>>) => ReactElement;
}
