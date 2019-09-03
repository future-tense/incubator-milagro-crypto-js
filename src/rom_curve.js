/*
    Licensed to the Apache Software Foundation (ASF) under one
    or more contributor license agreements.  See the NOTICE file
    distributed with this work for additional information
    regarding copyright ownership.  The ASF licenses this file
    to you under the Apache License, Version 2.0 (the
    "License"); you may not use this file except in compliance
    with the License.  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing,
    software distributed under the License is distributed on an
    "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
    KIND, either express or implied.  See the License for the
    specific language governing permissions and limitations
    under the License.
*/

/* Fixed Data in ROM - Field and Curve parameters */

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = {

        // BN254 Curve

        // Base Bits= 24

        Curve_Cof_I : 1,
        CURVE_A: 0,
        CURVE_B_I: 2,
        CURVE_B: [0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
        CURVE_Order: [0xD, 0x0, 0x10A100, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
        CURVE_Gx: [0x12, 0x0, 0x13A700, 0x0, 0x210000, 0x861, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
        CURVE_Gy: [0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],

        CURVE_Bnx: [0x1, 0x0, 0x4080, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
        CURVE_Cof: [0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
        CURVE_Cru: [0x7, 0x0, 0x6CD80, 0x0, 0x90000, 0x249, 0x400000, 0x49B362, 0x0, 0x0, 0x0],
        CURVE_Pxa: [0x3FB2B, 0x4224C8, 0xD91EE, 0x4898BF, 0x648BBB, 0xEDB6A4, 0x7E8C61, 0xEB8D8C, 0x9EB62F, 0x10BB51, 0x61A],
        CURVE_Pxb: [0xD54CF3, 0x34C1E7, 0xB70D8C, 0xAE3784, 0x4D746B, 0xAA5B1F, 0x8C5982, 0x310AA7, 0x737833, 0xAAF9BA, 0x516],
        CURVE_Pya: [0xCD2B9A, 0xE07891, 0xBD19F0, 0xBDBE09, 0xBD0AE6, 0x822329, 0x96698C, 0x9A90E0, 0xAF9343, 0x97A06B, 0x218],
        CURVE_Pyb: [0x3ACE9B, 0x1AEC6B, 0x578A2D, 0xD739C9, 0x9006FF, 0x8D37B0, 0x56F5F3, 0x8F6D44, 0x8B1526, 0x2B0E7C, 0xEBB],
        CURVE_W: [
            [0x3, 0x0, 0x20400, 0x0, 0x818000, 0x61, 0x0, 0x0, 0x0, 0x0, 0x0],
            [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
        ],
        CURVE_SB: [
            [
                [0x4, 0x0, 0x28500, 0x0, 0x818000, 0x61, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
            ],
            [
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0xA, 0x0, 0xE9D00, 0x0, 0x1E0000, 0x79E, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523]
            ]
        ],
        CURVE_WB: [
            [0x0, 0x0, 0x4080, 0x0, 0x808000, 0x20, 0x0, 0x0, 0x0, 0x0, 0x0],
            [0x5, 0x0, 0x54A80, 0x0, 0x70000, 0x1C7, 0x800000, 0x312241, 0x0, 0x0, 0x0],
            [0x3, 0x0, 0x2C580, 0x0, 0x838000, 0xE3, 0xC00000, 0x189120, 0x0, 0x0, 0x0],
            [0x1, 0x0, 0xC180, 0x0, 0x808000, 0x20, 0x0, 0x0, 0x0, 0x0, 0x0]
        ],
        CURVE_BB: [
            [
                [0xD, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0xC, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0xC, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0x2, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
            ],
            [
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0xC, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0xD, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0xC, 0x0, 0x106080, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523]
            ],
            [
                [0x2, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0x1, 0x0, 0x8100, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
            ],
            [
                [0x2, 0x0, 0x4080, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0x2, 0x0, 0x10200, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
                [0xA, 0x0, 0x102000, 0x0, 0x9F8000, 0x7FF, 0x800000, 0xBA344D, 0x1, 0x648240, 0x2523],
                [0x2, 0x0, 0x4080, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]
            ]
        ],

        USE_GLV: true,
        USE_GS_G2: true,
        USE_GS_GT: true,
        GT_STRONG: false,

        //debug: false,
    };
}
