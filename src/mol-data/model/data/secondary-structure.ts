/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const enum SSF {
    None = 0x0,

    // category
    DoubleHelix = 0x1,
    Helix = 0x2,
    Beta = 0x4,
    Turn = 0x8,

    // category variant
    LeftHanded = 0x10,  // helix
    RightHanded = 0x20,

    ClassicTurn = 0x40,  // turn
    InverseTurn = 0x80,

    // sub-category
    HelixOther = 0x100,  // protein
    Helix27 = 0x200,
    Helix3Ten = 0x400,
    HelixAlpha = 0x800,
    HelixGamma = 0x1000,
    HelixOmega = 0x2000,
    HelixPi = 0x4000,
    HelixPolyproline = 0x8000,

    DoubleHelixOther = 0x10000,  // nucleic
    DoubleHelixZ = 0x20000,
    DoubleHelixA = 0x40000,
    DoubleHelixB = 0x80000,

    BetaOther = 0x100000,  // protein
    BetaStrand = 0x200000,  // single strand
    BetaSheet = 0x400000,  // multiple hydrogen bonded strands
    BetaBarell = 0x800000,  // closed series of sheets

    TurnOther = 0x1000000,  // protein
    Turn1 = 0x2000000,
    Turn2 = 0x4000000,
    Turn3 = 0x8000000,

    NA = 0x10000000,  // not applicable/available
}
export { SSF as SecondaryStructureFlag }

export const SecondaryStructureMmcif: { [value: string]: number } = {
    HELX_LH_27_P: SSF.Helix | SSF.LeftHanded | SSF.Helix27,  // left-handed 2-7 helix (protein)
    HELX_LH_3T_P: SSF.Helix | SSF.LeftHanded | SSF.Helix3Ten,  // left-handed 3-10 helix (protein)
    HELX_LH_AL_P: SSF.Helix | SSF.LeftHanded | SSF.HelixAlpha,  // left-handed alpha helix (protein)
    HELX_LH_A_N: SSF.DoubleHelix | SSF.LeftHanded | SSF.DoubleHelixA,  // left-handed A helix (nucleic acid)
    HELX_LH_B_N: SSF.DoubleHelix | SSF.LeftHanded | SSF.DoubleHelixB,  // left-handed B helix (nucleic acid)
    HELX_LH_GA_P: SSF.Helix | SSF.LeftHanded | SSF.HelixGamma,  // left-handed gamma helix (protein)
    HELX_LH_N: SSF.DoubleHelix | SSF.LeftHanded,  // left-handed helix with type not specified (nucleic acid)
    HELX_LH_OM_P: SSF.Helix | SSF.LeftHanded | SSF.HelixOmega,  // left-handed omega helix (protein)
    HELX_LH_OT_N: SSF.DoubleHelix | SSF.LeftHanded | SSF.DoubleHelixOther,  // left-handed helix with type that does not conform to an accepted category (nucleic acid)
    HELX_LH_OT_P: SSF.Helix | SSF.LeftHanded | SSF.HelixOther,  // left-handed helix with type that does not conform to an accepted category (protein)
    HELX_LH_P: SSF.Helix | SSF.LeftHanded,  // left-handed helix with type not specified (protein)
    HELX_LH_PI_P: SSF.Helix | SSF.LeftHanded | SSF.HelixPi,  // left-handed pi helix (protein)
    HELX_LH_PP_P: SSF.Helix | SSF.LeftHanded | SSF.HelixPolyproline,  // left-handed polyproline helix (protein)
    HELX_LH_Z_N: SSF.DoubleHelix | SSF.LeftHanded | SSF.DoubleHelixZ,  // left-handed Z helix (nucleic acid)
    HELX_N: SSF.DoubleHelix,  // helix with handedness and type not specified (nucleic acid)
    HELX_OT_N: SSF.DoubleHelix,  // helix with handedness and type that do not conform to an accepted category (nucleic acid)
    HELX_OT_P: SSF.Helix,  // helix with handedness and type that do not conform to an accepted category (protein)
    HELX_P: SSF.Helix,  // helix with handedness and type not specified (protein)
    HELX_RH_27_P: SSF.Helix | SSF.RightHanded | SSF.Helix27,  // right-handed 2-7 helix (protein)
    HELX_RH_3T_P: SSF.Helix | SSF.RightHanded | SSF.Helix3Ten,  // right-handed 3-10 helix (protein)
    HELX_RH_AL_P: SSF.Helix | SSF.RightHanded | SSF.HelixAlpha,  // right-handed alpha helix (protein)
    HELX_RH_A_N: SSF.DoubleHelix | SSF.RightHanded | SSF.DoubleHelixA,  // right-handed A helix (nucleic acid)
    HELX_RH_B_N: SSF.DoubleHelix | SSF.RightHanded | SSF.DoubleHelixB,  // right-handed B helix (nucleic acid)
    HELX_RH_GA_P: SSF.Helix | SSF.RightHanded | SSF.HelixGamma,  // right-handed gamma helix (protein)
    HELX_RH_N: SSF.DoubleHelix | SSF.RightHanded,  // right-handed helix with type not specified (nucleic acid)
    HELX_RH_OM_P: SSF.Helix | SSF.RightHanded | SSF.HelixOmega,  // right-handed omega helix (protein)
    HELX_RH_OT_N: SSF.DoubleHelix | SSF.RightHanded | SSF.DoubleHelixOther,  // right-handed helix with type that does not conform to an accepted category (nucleic acid)
    HELX_RH_OT_P: SSF.Helix | SSF.RightHanded | SSF.HelixOther,  // right-handed helix with type that does not conform to an accepted category (protein)
    HELX_RH_P: SSF.Helix | SSF.RightHanded,  // right-handed helix with type not specified (protein)
    HELX_RH_PI_P: SSF.Helix | SSF.RightHanded | SSF.HelixPi,  // right-handed pi helix (protein)
    HELX_RH_PP_P: SSF.Helix | SSF.RightHanded | SSF.HelixPolyproline,  // right-handed polyproline helix (protein)
    HELX_RH_Z_N: SSF.DoubleHelix | SSF.RightHanded | SSF.DoubleHelixZ,  // right-handed Z helix (nucleic acid)
    STRN: SSF.Beta | SSF.BetaStrand,  // beta strand (protein)
    TURN_OT_P: SSF.Turn | SSF.TurnOther,  // turn with type that does not conform to an accepted category (protein)
    TURN_P: SSF.Turn,  // turn with type not specified (protein)
    TURN_TY1P_P: SSF.Turn | SSF.InverseTurn | SSF.Turn1,  // type I prime turn (protein)
    TURN_TY1_P: SSF.Turn | SSF.ClassicTurn | SSF.Turn1,  // type I turn (protein)
    TURN_TY2P_P: SSF.Turn | SSF.InverseTurn | SSF.Turn2,  // type II prime turn (protein)
    TURN_TY2_P: SSF.Turn | SSF.ClassicTurn | SSF.Turn2,  // type II turn (protein)
    TURN_TY3P_P: SSF.Turn | SSF.InverseTurn | SSF.Turn3,  // type III prime turn (protein)
    TURN_TY3_P: SSF.Turn | SSF.ClassicTurn | SSF.Turn3,  // type III turn (protein)
}

export const SecondaryStructurePdb: { [value: string]: number } = {
    1: SSF.Helix | SSF.RightHanded | SSF.HelixAlpha,  // Right-handed alpha (default)
    2: SSF.Helix | SSF.RightHanded | SSF.HelixOmega,  // Right-handed omega
    3: SSF.Helix | SSF.RightHanded | SSF.HelixPi,  // Right-handed pi
    4: SSF.Helix | SSF.RightHanded | SSF.HelixGamma,  // Right-handed gamma
    5: SSF.Helix | SSF.RightHanded | SSF.Helix3Ten,  // Right-handed 310
    6: SSF.Helix | SSF.LeftHanded | SSF.HelixAlpha,  // Left-handed alpha
    7: SSF.Helix | SSF.LeftHanded | SSF.HelixOmega,  // Left-handed omega
    8: SSF.Helix | SSF.LeftHanded | SSF.HelixGamma,  // Left-handed gamma
    9: SSF.Helix | SSF.Helix27,  // 27 ribbon/helix
    10: SSF.Helix | SSF.HelixPolyproline,  // Polyproline
}

export const SecondaryStructureStride: { [value: string]: number } = {
    H: SSF.Helix | SSF.HelixAlpha,  // Alpha helix
    G: SSF.Helix | SSF.Helix3Ten,  // 3-10 helix
    I: SSF.Helix | SSF.HelixPi,  // PI-helix
    E: SSF.Beta | SSF.BetaSheet,  // Extended conformation
    B: SSF.Beta | SSF.BetaStrand,  // Isolated bridge
    b: SSF.Beta | SSF.BetaStrand,  // Isolated bridge
    T: SSF.Turn,  // Turn
    C: SSF.NA,  // Coil (none of the above)
}

export const SecondaryStructureDssp: { [value: string]: number } = {
    H: SSF.Helix | SSF.HelixAlpha,  // alpha-helix
    B: SSF.Beta | SSF.BetaStrand,  // residue in isolated beta-bridge
    E: SSF.Beta | SSF.BetaSheet,  // extended strand, participates in beta ladder
    G: SSF.Helix | SSF.Helix3Ten,  // 3-helix (310 helix)
    I: SSF.Helix | SSF.HelixPi,  // 5 helix (pi-helix)
    T: SSF.Turn,  // hydrogen bonded turn
    S: SSF.Turn,  // bend
}