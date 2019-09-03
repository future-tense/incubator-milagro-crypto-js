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

const BIG = require("./big").BIG;
const FP = require("./fp");
const ROM_FIELD = require("./rom_field");
const ROM_CURVE = require("./rom_curve");

/* Elliptic Curve Point class */

"use strict";

/**
 * Creates an instance of ECP
 *
 * @constructor
 * @this {ECP}
 */
const ECP = function(input) {
    if (input instanceof ECP) {
        // copy constructor
        this.x = new FP(input.x);
        this.y = new FP(input.y);
        this.z = new FP(input.z);
    } else {
        // default constructor (point at infinity)
        this.x = new FP(0);
        this.y = new FP(1);
        this.z = new FP(0);
    }
};

ECP.ATE_BITS = 66;

ECP.prototype = {

    /**
     * Tests for ECP point equal to infinity
     *
     * @this {ECP}
     * @param 1 if infinity, else returns 0
     */
    is_infinity: function() {

        this.x.reduce();
        this.z.reduce();
        this.y.reduce();
        return (this.x.iszilch() && this.z.iszilch());
    },

    /**
     * conditional swap of this and Q dependant on dCopy ECP point to another ECP point
     *
     * @this {ECP}
     */
    cswap: function(Q, d) {

        this.x.cswap(Q.x, d);
        this.y.cswap(Q.y, d);
        this.z.cswap(Q.z, d);

    },

    /**
     * conditional move of Q to P dependant on d
     *
     * @this {ECP}
     */
    cmove: function(Q, d) {

        this.x.cmove(Q.x, d);
        this.y.cmove(Q.y, d);
        this.z.cmove(Q.z, d);
    },

    /**
     * Constant time select from pre-computed table
     *
     * @this {ECP}
     */
    select: function(W, b) {
        const m = b >> 31;
        let babs = (b ^ m) - m;
        babs = (babs - 1) / 2;

        this.cmove(W[0], ECP.teq(babs, 0)); // conditional move
        this.cmove(W[1], ECP.teq(babs, 1));
        this.cmove(W[2], ECP.teq(babs, 2));
        this.cmove(W[3], ECP.teq(babs, 3));
        this.cmove(W[4], ECP.teq(babs, 4));
        this.cmove(W[5], ECP.teq(babs, 5));
        this.cmove(W[6], ECP.teq(babs, 6));
        this.cmove(W[7], ECP.teq(babs, 7));

        const MP = new ECP();
        MP.copy(this);
        MP.neg();
        this.cmove(MP, (m & 1));
    },

    /* Test P == Q */

    equals: function(Q) {

        const a = new FP(0);
        const b = new FP(0);
        a.copy(this.x);
        a.mul(Q.z);
        a.reduce();
        b.copy(Q.x);
        b.mul(this.z);
        b.reduce();

        if (!a.equals(b)) {
            return false;
        }

        a.copy(this.y);
        a.mul(Q.z);
        a.reduce();
        b.copy(Q.y);
        b.mul(this.z);
        b.reduce();
        return (a.equals(b));
    },

    /**
     * Copy ECP point to another ECP point
     *
     * @this {ECP}
     * @param P ECP instance
     */
    copy: function(P) {
        this.x.copy(P.x);
        this.y.copy(P.y);
        this.z.copy(P.z);
    },

    /**
     * set this=-this
     *
     * @this {ECP}
     */
    neg: function() {
        this.y.neg();
        this.y.norm();
    },

    /**
     * Set ECP to point-at-infinity
     *
     * @this {ECP}
     */
    inf: function() {
        this.x.zero();
        this.y.one();
        this.z.zero();
    },

    /**
     * set this=(x,y)
     *
     * @this {ECP}
     * @param ix x-value
     * @param iy y-value
     */
    setxy: function(ix, iy) {
        this.x = new FP(0);
        this.x.bcopy(ix);

        this.y = new FP(0);
        this.y.bcopy(iy);
        this.z = new FP(1);
        this.x.norm();
        const rhs = ECP.RHS(this.x);

        const y2 = new FP(0);
        y2.copy(this.y);
        y2.sqr();

        if (!y2.equals(rhs)) {
            this.inf();
        }
    },

    /**
     * set this=x, where x is BIG, y is derived from sign s
     *
     * @this {ECP}
     * @param ix x-value
     * @param s sign to derive y
     */
    setxi: function(ix, s) {
        this.x = new FP(0);
        this.x.bcopy(ix);
        this.x.norm();
        const rhs = ECP.RHS(this.x);
        this.z = new FP(1);

        if (rhs.jacobi() === 1) {
            const ny = rhs.sqrt();
            if (ny.redc().parity() !== s) {
                ny.neg();
            }
            this.y = ny;
        } else {
            this.inf();
        }
    },

    /**
     *
     * set this=x, y calculated from curve equation
     *
     * @this {ECP}
     * @param ix x-value
     */
    setx: function(ix) {
        this.x = new FP(0);
        this.x.bcopy(ix);
        this.x.norm();
        const rhs = ECP.RHS(this.x);
        this.z = new FP(1);

        if (rhs.jacobi() === 1) {
            this.y = rhs.sqrt();
        } else {
            this.inf();
        }
    },

    /**
     * convert this to affine, from (x,y,z) to (x,y)
     *
     * @this {ECP}
     */
    affine: function() {
        if (this.is_infinity()) {
            return;
        }

        const one = new FP(1);
        if (this.z.equals(one)) {
            return;
        }

        this.z.inverse();

        this.x.mul(this.z);
        this.x.reduce();
        this.y.mul(this.z);
        this.y.reduce();
        this.z = one;
    },

    /**
     * extract affine x as FP2
     *
     * @this {ECP}
     */
    getX: function() {
        const W = new ECP();
        W.copy(this);
        W.affine();
        return W.x.redc();
    },

    /**
     * extract affine y as FP2
     *
     * @this {ECP}
     */
    getY: function() {
        const W = new ECP();
        W.copy(this);
        W.affine();
        return W.y.redc();
    },

    /**
     * get sign of Y
     *
     * @this {ECP}
     */
    getS: function() {
        const y = this.getY();
        return y.parity();
    },

    /**
     * extract x as FP
     *
     * @this {ECP}
     */
    getx: function() {
        return this.x;
    },

    /**
     * extract y as FP
     *
     * @this {ECP}
     */
    gety: function() {
        return this.y;
    },

    /**
     * extract z as FP
     *
     * @this {ECP}
     */
    getz: function() {
        return this.z;
    },

    /**
     * convert this to byte arrayextract projective x
     *
     * @this {ECP}
     * @param b byte array output
     * @param compress boolean
     */
    toBytes: function(b, compress) {

        const t = [];

        const W = new ECP();
        W.copy(this);
        W.affine();
        W.x.redc().toBytes(t);

        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i + 1] = t[i];
        }

        if (compress) {
            b[0] = 0x02;
            if (W.y.redc().parity() === 1) {
                b[0] = 0x03;
            }
            return;
        }

        b[0] = 0x04;

        W.y.redc().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i + BIG.MODBYTES + 1] = t[i];
        }
    },

    /**
     * convert this to hex string
     *
     * @this {ECP}
     * @return hex string
     */
    toString: function() {
        const W = new ECP();
        W.copy(this);
        if (W.is_infinity()) {
            return "infinity";
        }

        W.affine();
        return "(" + W.x.redc().toString() + "," + W.y.redc().toString() + ")";
    },

    /**
     * this+=this
     *
     * @this {ECP}
     */
    dbl: function() {

        const t0 = new FP(0);
        t0.copy(this.y);
        t0.sqr();
        const t1 = new FP(0);
        t1.copy(this.y);
        t1.mul(this.z);
        const t2 = new FP(0);
        t2.copy(this.z);
        t2.sqr();

        this.z.copy(t0);
        this.z.add(t0);
        this.z.norm();
        this.z.add(this.z);
        this.z.add(this.z);
        this.z.norm();

        t2.imul(3 * ROM_CURVE.CURVE_B_I);

        const x3 = new FP(0);
        x3.copy(t2);
        x3.mul(this.z);
        const y3 = new FP(0);
        y3.copy(t0);
        y3.add(t2);
        y3.norm();
        this.z.mul(t1);
        t1.copy(t2);
        t1.add(t2);
        t2.add(t1);
        t0.sub(t2);
        t0.norm();
        y3.mul(t0);
        y3.add(x3);
        t1.copy(this.x);
        t1.mul(this.y);
        this.x.copy(t0);
        this.x.norm();
        this.x.mul(t1);
        this.x.add(this.x);

        this.x.norm();
        this.y.copy(y3);
        this.y.norm();
    },

    /**
     * Adds ECP instances
     *
     * param Q ECP instance
     * @this {ECP}
     */
    add: function(Q) {
        const b = 3 * ROM_CURVE.CURVE_B_I;
        const t0 = new FP(0);
        t0.copy(this.x);
        t0.mul(Q.x);
        const t1 = new FP(0);
        t1.copy(this.y);
        t1.mul(Q.y);
        const t2 = new FP(0);
        t2.copy(this.z);
        t2.mul(Q.z);
        const t3 = new FP(0);
        t3.copy(this.x);
        t3.add(this.y);
        t3.norm();
        const t4 = new FP(0);
        t4.copy(Q.x);
        t4.add(Q.y);
        t4.norm();
        t3.mul(t4);
        t4.copy(t0);
        t4.add(t1);

        t3.sub(t4);
        t3.norm();
        t4.copy(this.y);
        t4.add(this.z);
        t4.norm();
        const x3 = new FP(0);
        x3.copy(Q.y);
        x3.add(Q.z);
        x3.norm();

        t4.mul(x3);
        x3.copy(t1);
        x3.add(t2);

        t4.sub(x3);
        t4.norm();
        x3.copy(this.x);
        x3.add(this.z);
        x3.norm();
        const y3 = new FP(0);
        y3.copy(Q.x);
        y3.add(Q.z);
        y3.norm();
        x3.mul(y3);
        y3.copy(t0);
        y3.add(t2);
        y3.rsub(x3);
        y3.norm();
        x3.copy(t0);
        x3.add(t0);
        t0.add(x3);
        t0.norm();
        t2.imul(b);

        const z3 = new FP(0);
        z3.copy(t1);
        z3.add(t2);
        z3.norm();
        t1.sub(t2);
        t1.norm();
        y3.imul(b);

        x3.copy(y3);
        x3.mul(t4);
        t2.copy(t3);
        t2.mul(t1);
        x3.rsub(t2);
        y3.mul(t0);
        t1.mul(z3);
        y3.add(t1);
        t0.mul(t3);
        z3.mul(t4);
        z3.add(t0);

        this.x.copy(x3);
        this.x.norm();
        this.y.copy(y3);
        this.y.norm();
        this.z.copy(z3);
        this.z.norm();
    },

    /**
     * Subtracts ECP instance Q  from this
     *
     * @this {ECP}
     * @param Q ECP instance
     */
    sub: function(Q) {
        const NQ = new ECP();
        NQ.copy(Q);
        NQ.neg();
        this.add(NQ);
    },

    /**
     * constant time multiply by small integer of length bts - use ladder
     *
     * @this {ECP}
     * @param e small integer
     * @param bts e bit length
     */
    pinmul: function(e, bts) {
        const P = new ECP();
        const R0 = new ECP();
        const R1 = new ECP();
        R1.copy(this);

        for (let i = bts - 1; i >= 0; i--) {
            const b = (e >> i) & 1;
            P.copy(R1);
            P.add(R0);
            R0.cswap(R1, b);
            R1.copy(P);
            R0.dbl();
            R0.cswap(R1, b);
        }

        P.copy(R0);
        P.affine();

        return P;
    },

    /**
     * Multiplies an ECP instance P by a BIG, side-channel resistant
     *
     * @this {ECP}
     * @param e BIG number multiplier
     */
    mul: function(e) {
        if (e.iszilch() || this.is_infinity()) {
            return new ECP();
        }

        const P = new ECP();

        // fixed size windows
        const mt = new BIG();
        const t = new BIG();
        const Q = new ECP();
        const C = new ECP();
        const W = [];
        const w = [];

        // precompute table
        Q.copy(this);
        Q.dbl();
        W[0] = new ECP();
        W[0].copy(this);

        for (let i = 1; i < 8; i++) {
            W[i] = new ECP();
            W[i].copy(W[i - 1]);
            W[i].add(Q);
        }

        // make exponent odd - add 2P if even, P if odd
        t.copy(e);
        const s = t.parity();
        t.inc(1);
        t.norm();
        const ns = t.parity();
        mt.copy(t);
        mt.inc(1);
        mt.norm();
        t.cmove(mt, s);
        Q.cmove(this, ns);
        C.copy(Q);

        const nb = 1 + Math.floor((t.nbits() + 3) / 4);

        // convert exponent to signed 4-bit window
        for (let i = 0; i < nb; i++) {
            w[i] = (t.lastbits(5) - 16);
            t.dec(w[i]);
            t.norm();
            t.fshr(4);
        }
        w[nb] = t.lastbits(5);

        P.copy(W[Math.floor((w[nb] - 1) / 2)]);
        for (let i = nb - 1; i >= 0; i--) {
            Q.select(W, w[i]);
            P.dbl();
            P.dbl();
            P.dbl();
            P.dbl();
            P.add(Q);
        }
        P.sub(C);
        P.affine();

        return P;
    },

    /**
     * Return e.this+f.Q
     *
     * @this {ECP}
     * @param e BIG number multiplier
     * @param Q ECP instance
     * @param f BIG number multiplier
     */
    mul2: function(e, Q, f) {
        const te = new BIG();
        te.copy(e);

        const tf = new BIG();
        tf.copy(f);

        const S = new ECP();
        const T = new ECP();

        // precompute table
        const W = [];
        W[1] = new ECP();
        W[1].copy(this);
        W[1].sub(Q);
        W[2] = new ECP();
        W[2].copy(this);
        W[2].add(Q);
        S.copy(Q);
        S.dbl();
        W[0] = new ECP();
        W[0].copy(W[1]);
        W[0].sub(S);
        W[3] = new ECP();
        W[3].copy(W[2]);
        W[3].add(S);
        T.copy(this);
        T.dbl();
        W[5] = new ECP();
        W[5].copy(W[1]);
        W[5].add(T);
        W[6] = new ECP();
        W[6].copy(W[2]);
        W[6].add(T);
        W[4] = new ECP();
        W[4].copy(W[5]);
        W[4].sub(S);
        W[7] = new ECP();
        W[7].copy(W[6]);
        W[7].add(S);

        // if multiplier is odd, add 2, else add 1 to multiplier, and add 2P or P to correction

        const mt = new BIG();
        const C = new ECP();

        let s = te.parity();
        te.inc(1);
        te.norm();
        let ns = te.parity();
        mt.copy(te);
        mt.inc(1);
        mt.norm();
        te.cmove(mt, s);
        T.cmove(this, ns);
        C.copy(T);

        s = tf.parity();
        tf.inc(1);
        tf.norm();
        ns = tf.parity();
        mt.copy(tf);
        mt.inc(1);
        mt.norm();
        tf.cmove(mt, s);
        S.cmove(Q, ns);
        C.add(S);

        mt.copy(te);
        mt.add(tf);
        mt.norm();
        const nb = 1 + Math.floor((mt.nbits() + 1) / 2);

        // convert exponent to signed 2-bit window
        const w = [];
        for (let i = 0; i < nb; i++) {
            const a = (te.lastbits(3) - 4);
            te.dec(a);
            te.norm();
            te.fshr(2);
            const b = (tf.lastbits(3) - 4);
            tf.dec(b);
            tf.norm();
            tf.fshr(2);
            w[i] = (4 * a + b);
        }
        w[nb] = (4 * te.lastbits(3) + tf.lastbits(3));
        S.copy(W[Math.floor((w[nb] - 1) / 2)]);

        for (let i = nb - 1; i >= 0; i--) {
            T.select(W, w[i]);
            S.dbl();
            S.dbl();
            S.add(T);
        }
        S.sub(C); /* apply correction */
        S.affine();

        return S;
    }
};

/**
  * Set group generator
  *
  * @this {ECP}
  */
ECP.generator = function() {
    const gx = new BIG(0);
    gx.rcopy(ROM_CURVE.CURVE_Gx);

    const gy = new BIG(0);
    gy.rcopy(ROM_CURVE.CURVE_Gy);

    const G = new ECP();
    G.setxy(gx, gy);
    return G;
};

/* return 1 if b===c, no branching */
ECP.teq = function(b, c) {
    const x = (b ^ c) - 1;
    return ((x >> 31) & 1);
};

/**
  * convert from byte array to point
  *
  * @this {ECP}
  * @param b input byte array
  */
ECP.fromBytes = function(b) {

    const p = new BIG(0);
    p.rcopy(ROM_FIELD.Modulus);

    const t = [];
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = b[i + 1];
    }

    const P = new ECP();
    const px = BIG.fromBytes(t);
    if (BIG.comp(px, p) >= 0) {
        return P;
    }

    if (b[0] === 0x04) {
        for (let i = 0; i < BIG.MODBYTES; i++) {
            t[i] = b[i + BIG.MODBYTES + 1];
        }

        const py = BIG.fromBytes(t);

        if (BIG.comp(py, p) >= 0) {
            return P;
        }

        P.setxy(px, py);

        return P;
    }

    if (b[0] === 0x02 || b[0] === 0x03) {
        P.setxi(px, b[0] & 1);
        return P;
    }

    return P;
};

/**
  * Calculate RHS of the curve equation
  *
  * @this {ECP}
  * @param x x-value
  */
ECP.RHS = function(x) {
    const r = new FP(0);
    r.copy(x);
    r.sqr();
    r.mul(x);

    const b = new FP(ROM_CURVE.CURVE_B_I);
    r.add(b);

    r.reduce();
    return r;
};

// CommonJS module exports
if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = ECP;
}
