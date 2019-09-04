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
const FP2 = require("./fp2");
const ROM_CURVE = require("./rom_curve");

/* AMCL Weierstrass elliptic curve functions over FP2 */

"use strict";

/**
 * Creates an instance of ECP2
 *
 * @constructor
 * @this {ECP2}
 */
const ECP2 = function(input) {
    if (input instanceof ECP2) {
        // copy constructor
        this.x = new FP2(input.x);
        this.y = new FP2(input.y);
        this.z = new FP2(input.z);
    } else {
        // default constructor (point at infinity)
        this.x = new FP2(0);
        this.y = new FP2(1);
        this.z = new FP2(0);
    }
};

ECP2.prototype = {

    /**
     * Tests for ECP2 point equal to infinity
     *
     * @this {ECP2}
     * @param 1 if infinity, else returns 0
     */
    is_infinity: function() {

        this.x.reduce();
        this.y.reduce();
        this.z.reduce();
        return (this.x.iszilch() && this.z.iszilch());
    },

    /**
     * Copy ECP2 point to another ECP2 point
     *
     * @this {ECP2}
     * @param P ECP2 instance
     */
    copy: function(P) {
        this.x.copy(P.x);
        this.y.copy(P.y);
        this.z.copy(P.z);
        return this;
    },

    /**
     * Set ECP2 to point-at-infinity
     *
     * @this {ECP2}
     */
    inf: function() {
        this.x.zero();
        this.y.one();
        this.z.zero();
    },

    /**
     * conditional move of Q to P dependant on d
     *
     * @this {ECP2}
     */
    cmove: function(Q, d) {
        this.x.cmove(Q.x, d);
        this.y.cmove(Q.y, d);
        this.z.cmove(Q.z, d);
    },

    /**
     * Constant time select from pre-computed table
     *
     * @this {ECP2}
     */
    select: function(W, b) {

        const m = b >> 31;
        let babs = (b ^ m) - m;
        babs = (babs - 1) / 2;

        this.cmove(W[0], ECP2.teq(babs, 0)); // conditional move
        this.cmove(W[1], ECP2.teq(babs, 1));
        this.cmove(W[2], ECP2.teq(babs, 2));
        this.cmove(W[3], ECP2.teq(babs, 3));
        this.cmove(W[4], ECP2.teq(babs, 4));
        this.cmove(W[5], ECP2.teq(babs, 5));
        this.cmove(W[6], ECP2.teq(babs, 6));
        this.cmove(W[7], ECP2.teq(babs, 7));

        const MP = new ECP2();
        MP.copy(this);
        MP.neg();
        this.cmove(MP, (m & 1));
    },

    /**
     * Test P === Q
     *
     * @this {ECP2}
     * @param Q ECP2 instance
     */
    equals: function(Q) {

        const a = new FP2(0);
        a.copy(this.x);
        const b = new FP2(0);
        b.copy(Q.x);

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
     * set this=-this
     *
     * @this {ECP2}
     */
    neg: function() {
        this.y.norm();
        this.y.neg();
        this.y.norm();
        return this;
    },

    /**
     * convert this to affine, from (x,y,z) to (x,y)
     *
     * @this {ECP2}
     */
    affine: function() {
        if (this.is_infinity()) {
            return this;
        }

        const one = new FP2(1);

        if (this.z.equals(one)) {
            this.x.reduce();
            this.y.reduce();
            return this;
        }

        this.z.inverse();

        this.x.mul(this.z);
        this.x.reduce();
        this.y.mul(this.z);
        this.y.reduce();
        this.z.copy(one);

        return this;
    },

    /**
     * extract affine x as FP2
     *
     * @this {ECP2}
     */
    getX: function() {
        const W = new ECP2();
        W.copy(this);
        W.affine();
        return W.x;
    },


    /**
     * extract affine y as FP2
     *
     * @this {ECP2}
     */
    getY: function() {
        const W = new ECP2();
        W.copy(this);
        W.affine();
        return W.y;
    },

    /**
     * extract projective x
     *
     * @this {ECP2}
     */
    getx: function() {
        return this.x;
    },

    /**
     * extract projective y
     *
     * @this {ECP2}
     */
    gety: function() {
        return this.y;
    },

    /**
     * extract projective z
     *
     * @this {ECP2}
     */
    getz: function() {
        return this.z;
    },

    /**
     * convert this to byte arrayextract projective x
     *
     * @this {ECP2}
     * @param b byte array output
     */
    toBytes: function(b) {

        const W = new ECP2();
        W.copy(this);
        W.affine();

        const t = [];
        W.x.getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i] = t[i];
        }
        W.x.getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i + BIG.MODBYTES] = t[i];
        }

        W.y.getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i + 2 * BIG.MODBYTES] = t[i];
        }
        W.y.getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            b[i + 3 * BIG.MODBYTES] = t[i];
        }
    },

    /**
     * convert this to hex string
     *
     * @this {ECP2}
     * @return hex string
     */
    toString: function() {
        const W = new ECP2();
        W.copy(this);

        if (W.is_infinity()) {
            return "infinity";
        }

        W.affine();
        return "(" + W.x.toString() + "," + W.y.toString() + ")";
    },


    /**
     * set this=(x,y)
     *
     * @this {ECP2}
     * @param ix x-value
     * @param iy y-value
     */
    setxy: function(ix, iy) {
        this.x.copy(ix);
        this.y.copy(iy);
        this.z.one();
        this.x.norm();

        const rhs = ECP2.RHS(this.x);

        const y2 = new FP2(this.y);
        y2.sqr();

        if (!y2.equals(rhs)) {
            this.inf();
        }
    },

    /**
     * set this=(x,.)
     *
     * @this {ECP2}
     * @param ix x-value
     */
    setx: function(ix) {
        this.x.copy(ix);
        this.z.one();
        this.x.norm();

        const rhs = ECP2.RHS(this.x);

        if (rhs.sqrt()) {
            this.y.copy(rhs);
        } else {
            this.inf();
        }
    },

    /**
     * set this*=q, where q is Modulus, using Frobenius
     *
     * @this {ECP2}
     */
    frob: function(X) {
        const X2 = new FP2(X); //X2.copy(X);
        X2.sqr();
        this.x.conj();
        this.y.conj();
        this.z.conj();
        this.z.reduce();
        this.x.mul(X2);
        this.y.mul(X2);
        this.y.mul(X);
        return this;
    },

    /**
     * this+=this
     *
     * @this {ECP2}
     */
    dbl: function() {

        const iy = new FP2(0);
        iy.copy(this.y);
        iy.mul_ip();
        iy.norm();

        const t0 = new FP2(0);
        t0.copy(this.y);
        t0.sqr();
        t0.mul_ip();
        const t1 = new FP2(0);
        t1.copy(iy);
        t1.mul(this.z);
        const t2 = new FP2(0);
        t2.copy(this.z);
        t2.sqr();

        this.z.copy(t0);
        this.z.add(t0);
        this.z.norm();
        this.z.add(this.z);
        this.z.add(this.z);
        this.z.norm();

        t2.imul(3 * ROM_CURVE.CURVE_B_I);

        const x3 = new FP2(0);
        x3.copy(t2);
        x3.mul(this.z);

        const y3 = new FP2(0);
        y3.copy(t0);

        y3.add(t2);
        y3.norm();
        this.z.mul(t1);
        t1.copy(t2);
        t1.add(t2);
        t2.add(t1);
        t2.norm();
        t0.sub(t2);
        t0.norm(); //y^2-9bz^2
        y3.mul(t0);
        y3.add(x3); //(y^2+3z*2)(y^2-9z^2)+3b.z^2.8y^2
        t1.copy(this.x);
        t1.mul(iy); //
        this.x.copy(t0);
        this.x.norm();
        this.x.mul(t1);
        this.x.add(this.x); //(y^2-9bz^2)xy2

        this.x.norm();
        this.y.copy(y3);
        this.y.norm();

        return 1;
    },

    /**
     * Adds ECP2 instances
     *
     * param Q ECP2 instance
     * @this {ECP2}
     */
    add: function(Q) {

        const b = 3 * ROM_CURVE.CURVE_B_I;

        const t0 = new FP2(0);
        t0.copy(this.x);
        t0.mul(Q.x); // x.Q.x
        const t1 = new FP2(0);
        t1.copy(this.y);
        t1.mul(Q.y); // y.Q.y

        const t2 = new FP2(0);
        t2.copy(this.z);
        t2.mul(Q.z);
        const t3 = new FP2(0);
        t3.copy(this.x);
        t3.add(this.y);
        t3.norm(); //t3=X1+Y1
        const t4 = new FP2(0);
        t4.copy(Q.x);
        t4.add(Q.y);
        t4.norm(); //t4=X2+Y2
        t3.mul(t4); //t3=(X1+Y1)(X2+Y2)
        t4.copy(t0);
        t4.add(t1); //t4=X1.X2+Y1.Y2

        t3.sub(t4);
        t3.norm();
        t3.mul_ip();
        t3.norm(); //t3=(X1+Y1)(X2+Y2)-(X1.X2+Y1.Y2) = X1.Y2+X2.Y1

        t4.copy(this.y);
        t4.add(this.z);
        t4.norm(); //t4=Y1+Z1
        const x3 = new FP2(0);
        x3.copy(Q.y);
        x3.add(Q.z);
        x3.norm(); //x3=Y2+Z2

        t4.mul(x3); //t4=(Y1+Z1)(Y2+Z2)
        x3.copy(t1); //
        x3.add(t2); //X3=Y1.Y2+Z1.Z2

        t4.sub(x3);
        t4.norm();
        t4.mul_ip();
        t4.norm(); //t4=(Y1+Z1)(Y2+Z2) - (Y1.Y2+Z1.Z2) = Y1.Z2+Y2.Z1

        x3.copy(this.x);
        x3.add(this.z);
        x3.norm(); // x3=X1+Z1

        const y3 = new FP2(0);
        y3.copy(Q.x);
        y3.add(Q.z);
        y3.norm(); // y3=X2+Z2
        x3.mul(y3); // x3=(X1+Z1)(X2+Z2)
        y3.copy(t0);
        y3.add(t2); // y3=X1.X2+Z1+Z2
        y3.rsub(x3);
        y3.norm(); // y3=(X1+Z1)(X2+Z2) - (X1.X2+Z1.Z2) = X1.Z2+X2.Z1

        t0.mul_ip();
        t0.norm(); // x.Q.x
        t1.mul_ip();
        t1.norm(); // y.Q.y

        x3.copy(t0);
        x3.add(t0);
        t0.add(x3);
        t0.norm();
        t2.imul(b);

        const z3 = new FP2(0);
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

        return 0;
    },

    /**
     * Subtracts ECP instance Q  from this
     *
     * @this {ECP2}
     * @param Q ECP2 instance
     */
    sub: function(Q) {
        const NQ=new ECP2();
        NQ.copy(Q);
        NQ.neg();
        return this.add(NQ);
    },

    /**
     * Multiplies an ECP2 instance P by a BIG, side-channel resistant
     *
     * @this {ECP2}
     * @param e BIG number multiplier
     */
    mul: function(e) {
        /* fixed size windows */

        if (this.is_infinity()) {
            return new ECP2();
        }

        // precompute table
        const Q = new ECP2();
        Q.copy(this);
        Q.dbl();

        const W = [];
        W[0] = new ECP2();
        W[0].copy(this);

        for (let i = 1; i < 8; i++) {
            W[i] = new ECP2();
            W[i].copy(W[i - 1]);
            W[i].add(Q);
        }

        // make exponent odd - add 2P if even, P if odd
        const t = new BIG();
        t.copy(e);
        const s = t.parity();
        t.inc(1);
        t.norm();
        const ns = t.parity();
        const mt = new BIG();
        mt.copy(t);
        mt.inc(1);
        mt.norm();
        t.cmove(mt, s);
        Q.cmove(this, ns);
        const C = new ECP2();
        C.copy(Q);

        const nb = 1 + Math.floor((t.nbits() + 3) / 4);

        // convert exponent to signed 4-bit window
        const w = [];
        for (let i = 0; i < nb; i++) {
            w[i] = (t.lastbits(5) - 16);
            t.dec(w[i]);
            t.norm();
            t.fshr(4);
        }
        w[nb] = t.lastbits(5);

        const P = new ECP2();
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
    }
};

/**
  * Set group generator
  *
  * @this {ECP2}
  */
ECP2.generator = function() {
    const A = new BIG(0);
    const B = new BIG(0);
    A.rcopy(ROM_CURVE.CURVE_Pxa);
    B.rcopy(ROM_CURVE.CURVE_Pxb);
    const QX = new FP2(0);
    QX.bset(A, B);

    A.rcopy(ROM_CURVE.CURVE_Pya);
    B.rcopy(ROM_CURVE.CURVE_Pyb);
    const QY = new FP2(0);
    QY.bset(A, B);

    const G = new ECP2();
    G.setxy(QX, QY);
    return G;
};

/**
  * convert from byte array to point
  *
  * @this {ECP2}
  * @param b input byte array
  */
ECP2.fromBytes = function(b) {

    const rax = BIG.fromBytes(b.slice(0, 32));
    const rbx = BIG.fromBytes(b.slice(32, 64));
    const rx = new FP2(rax, rbx);

    const ray = BIG.fromBytes(b.slice(64, 96));
    const rby = BIG.fromBytes(b.slice(96, 128));
    const ry = new FP2(ray, rby);

    const P = new ECP2();
    P.setxy(rx, ry);
    return P;
};

/**
  * Calculate RHS of curve equation x^3+B
  *
  * @this {ECP2}
  * @param x x-value
  */
ECP2.RHS = function(x) {
    const r = new FP2(x);
    r.sqr();
    r.mul(x);

    const b = new FP2(ROM_CURVE.CURVE_B_I);
    b.div_ip();
    r.add(b);

    r.reduce();
    return r;
};

/* P=u0.Q0+u1*Q1+u2*Q2+u3*Q3 */
// Bos & Costello https://eprint.iacr.org/2013/458.pdf
// Faz-Hernandez & Longa & Sanchez  https://eprint.iacr.org/2013/158.pdf
// Side channel attack secure

/**
  * Calculate P=u0.Q0+u1*Q1+u2*Q2+u3*Q3
  *
  * @this {ECP2}
  */
ECP2.mul4 = function(Q, u) {

    const t = [];
    for (let i = 0; i < 4; i++) {
        t[i] = new BIG(u[i]);
        t[i].norm();
    }

    const T = [];
    T[0] = new ECP2(); T[0].copy(Q[0]); // Q[0]
    T[1] = new ECP2(); T[1].copy(T[0]); T[1].add(Q[1]); // Q[0]+Q[1]
    T[2] = new ECP2(); T[2].copy(T[0]); T[2].add(Q[2]); // Q[0]+Q[2]
    T[3] = new ECP2(); T[3].copy(T[1]); T[3].add(Q[2]); // Q[0]+Q[1]+Q[2]
    T[4] = new ECP2(); T[4].copy(T[0]); T[4].add(Q[3]); // Q[0]+Q[3]
    T[5] = new ECP2(); T[5].copy(T[1]); T[5].add(Q[3]); // Q[0]+Q[1]+Q[3]
    T[6] = new ECP2(); T[6].copy(T[2]); T[6].add(Q[3]); // Q[0]+Q[2]+Q[3]
    T[7] = new ECP2(); T[7].copy(T[3]); T[7].add(Q[3]); // Q[0]+Q[1]+Q[2]+Q[3]

    // Make it odd
    const pb = 1 - t[0].parity();
    t[0].inc(pb);
    t[0].norm();

    // Number of bits
    const mt = new BIG();
    mt.zero();
    for (let i = 0; i < 4; i++) {
        mt.or(t[i]);
    }

    const nb = 1 + mt.nbits();

    // Sign pivot
    const s = [];
    s[nb - 1] = 1;
    for (let i = 0; i < nb - 1; i++) {
        t[0].fshr(1);
        s[i] = 2 * t[0].parity() - 1;
    }

    // Recoded exponent
    const w = [];
    for (let i = 0; i < nb; i++) {
        w[i] = 0;
        let k = 1;
        for (let j = 1; j < 4; j++) {
            const bt = s[i] * t[j].parity();
            t[j].fshr(1);
            t[j].dec(bt >> 1);
            t[j].norm();
            w[i] += bt * k;
            k *= 2;
        }
    }

    // Main loop
    const P = new ECP2();
    const W = new ECP2();
    P.select(T, 2 * w[nb - 1] + 1);
    for (let i = nb - 2; i >= 0; i--) {
        P.dbl();
        W.select(T, 2 * w[i] + s[i]);
        P.add(W);
    }

    // apply correction
    W.copy(P);
    W.sub(Q[0]);
    P.cmove(W, pb);
    P.affine();
    return P;
};


/* return 1 if b===c, no branching */
ECP2.teq = function(b, c) {
    const x = (b ^ c) -1;
    return ((x >> 31) & 1);
};

// CommonJS module exports
if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = ECP2;
}
