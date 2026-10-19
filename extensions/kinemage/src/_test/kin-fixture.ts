/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small kinemage for the Prism tests: two masters, a group of dots, an animate group with two subgroups (dominant),
// a second animate group and one view.

export const Kin = `@kinemage 1
@caption Test kinemage for the Prism extension
@viewid {Overview}
@1viewid {Overview}
@1span 20
@1center 10 20 30
@master {bad overlap}
@master {clashes}
@group {Contacts}
@dotlist {Dots} color=yellow master={clashes}
{HN ALA 2}yellow 10.5,20.5,30.5
{HN ALA 3}yellow 11.5,21.5,31.5
@group {Structure} animate dominant
@subgroup {Backbone}
@vectorlist {CA trace} color=blue master={bad overlap}
{CA ALA 1}blue P 10.0,20.0,30.0 {CA ALA 2}blue 11.0,21.0,31.0
{"}blue 12.0,22.0,32.0 {CA ALA 3}blue 13.0,23.0,33.0
@subgroup {Sidechains}
@balllist {CB atoms} color=green radius=0.5
{CB ARG 1}green r=0.5 9.0,19.0,29.0
{CB ARG 2}green r=0.6 10.0,20.0,30.0
@group {Alternate} animate
@balllist {Waters} color=cyan radius=1.0
{HOH 101}cyan r=1.0 15.0,25.0,35.0
`;
