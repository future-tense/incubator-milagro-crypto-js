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
const FP2 = require("./fp2");
const FP4 = require("./fp4");
const ROM_FIELD = require("./rom_field");

/* AMCL Fp^12 functions */

/* FP12 elements are of the form a+i.b+i^2.c */

"use strict";

/**
  * Creates an instance of FP12.
  *
  * @constructor
  * @this {FP12}
  */
const FP12 = function(d, e, f) {
    if (d instanceof FP12) {
        // ignore e, d, which are assumed be undefined in this case
        this.a = new FP4(d.a);
        this.b = new FP4(d.b);
        this.c = new FP4(d.c);
        this.stype = FP.DENSE;
    } else if (typeof d !== "undefined" && typeof e !== "undefined" && typeof f !== "undefined") {
        // all 3 components set to (can be anything that the FP4 constructor supports)
        this.a = new FP4(d);
        this.b = new FP4(e);
        this.c = new FP4(f);
        this.stype = FP.DENSE;
    } else if (typeof d === "number") {
        // first component is number
        this.a = new FP4(d);
        this.b = new FP4(0);
        this.c = new FP4(0);
        if (d === 1) {
            this.stype = FP.ONE;
        } else {
            this.stype = FP.SPARSER;
        }
    } else {
        // other cases, including `new FP12()` fall back to zero
        this.a = new FP4(0);
        this.b = new FP4(0);
        this.c = new FP4(0);
        this.stype = FP.ZERO;
    }
};

FP12.prototype = {

    /**
     * Reduces all components of possibly unreduced FP12 mod Modulus
     *
     * @this {FP12}
     */
    reduce: function() {
        this.a.reduce();
        this.b.reduce();
        this.c.reduce();
        return this;
    },

    /**
     * Normalises the components of an FP12
     *
     * @this {FP12}
     */
    norm: function() {
        this.a.norm();
        this.b.norm();
        this.c.norm();
        return this;
    },

    /**
     * Tests for FP12 equal to zero
     *
     * @this {FP12}
     */
    iszilch: function() {
        return (this.a.iszilch() && this.b.iszilch() && this.c.iszilch());
    },

    /**
     * Tests for FP12 equal to unity
     *
     * @this {FP12}
     */
    isunity: function() {
        const one = new FP4(1);
        return (this.a.equals(one) && this.b.iszilch() && this.c.iszilch());
    },

    /**
     * Conditional copy of FP12 number
     *
     * @this {FP12}
     * @param g FP12 instance
     * @param d copy depends on this value
     */
    cmove: function(g, d) {
        this.a.cmove(g.a, d);
        this.b.cmove(g.b, d);
        this.c.cmove(g.c, d);
        d = ~(d-1);
        this.stype ^= (this.stype ^ g.stype) & d;
    },

    /**
     * Constant time select from pre-computed table 
     *
     * @this {FP12}
     */
    select: function(g, b) {
        const invf = new FP12(0);
        const m = b >> 31;
        let babs = (b ^ m) - m;
        babs = (babs - 1) / 2;

        this.cmove(g[0], FP12.teq(babs, 0)); // conditional move
        this.cmove(g[1], FP12.teq(babs, 1));
        this.cmove(g[2], FP12.teq(babs, 2));
        this.cmove(g[3], FP12.teq(babs, 3));
        this.cmove(g[4], FP12.teq(babs, 4));
        this.cmove(g[5], FP12.teq(babs, 5));
        this.cmove(g[6], FP12.teq(babs, 6));
        this.cmove(g[7], FP12.teq(babs, 7));

        invf.copy(this);
        invf.conj();
        this.cmove(invf, (m & 1));
    },

    settype: function(w) {
        this.stype = w;
    },

    gettype: function() {
        return this.stype;
    },

    /**
     * extract a from this
     *
     * @this {FP12}
     */
    geta: function() {
        return this.a;
    },

    /**
     * extract b from this
     *
     * @this {FP12}
     */
    getb: function() {
        return this.b;
    },

    /**
     * extract c from this
     *
     * @this {FP12}
     */
    getc: function() {
        return this.c;
    },

    /**
     * Tests for equality of two FP12s
     *
     * @this {FP12}
     * @param x FP12 instance to compare
     */
    equals: function(x) {
        return (this.a.equals(x.a) && this.b.equals(x.b) && this.c.equals(x.c));
    },

    /**
     * Copy FP12 to another FP12
     *
     * @this {FP12}
     * @param x FP12 instance to be copied
     */
    copy: function(x) {
        this.a.copy(x.a);
        this.b.copy(x.b);
        this.c.copy(x.c);
        this.stype = x.stype;
        return this;
    },

    /**
     * Set FP12 to unity
     *
     * @this {FP12}
     * @param x FP12 instance to be set to one
     */
    one: function() {
        this.a.one();
        this.b.zero();
        this.c.zero();
        this.stype = FP.ONE;
    },

    /**
     * Set FP12 to zero
     *
     * @this {FP12}
     */
    zero: function() {
        this.a.zero();
        this.b.zero();
        this.c.zero();
        this.stype = FP.ZERO;
    },

    /**
     * Conjugation of FP12
     *
     * @this {FP12}
     */
    conj: function() {
        this.a.conj();
        this.b.nconj();
        this.c.conj();
        return this;
    },

    /**
     * Set FP12 from three FP4 values
     *
     * @this {FP12}
     * @param d FP4 instance
     * @param e FP4 instance
     * @param f FP4 instance
     */
    set: function(d, e, f) {
        this.a.copy(d);
        this.b.copy(e);
        this.c.copy(f);
        this.stype = FP.DENSE;
    },

    /**
     * Set FP12 from one FP4 value
     *
     * @this {FP12}
     * @param d FP4 instance
     */
    seta: function(d) {
        this.a.copy(d);
        this.b.zero();
        this.c.zero();
        this.stype = FP.SPARSER;
    },

    /**
     * Fast Squaring of an FP12 in "unitary" form
     *
     * @this {FP12}
     */
    usqr: function() {
        const A = new FP4(this.a).nconj().dbl();
        const B = new FP4(this.c).sqr().times_i();
        const C = new FP4(this.b).sqr();

        this.a.sqr();
        const D = new FP4(this.a).dbl();
        this.a.add(D).add(A);       //  3*a^2 + A

        D.copy(B).dbl().add(B);     //  D = 3*B
        A.copy(C).dbl().add(C);     //  A = 3*C

        this.b.conj().dbl().add(D);
        this.c.nconj().dbl().add(A);
        this.stype = FP.DENSE;
        return this.reduce();
    },

    /**
     * Fast Squaring of an FP12
     *
     * @this {FP12}
     */
    sqr: function() {
        if (this.stype === FP.ONE) {
            return this;
        }

        const A = new FP4(this.a).sqr();
        const B = new FP4(this.b).mul(this.c).dbl();
        const C = new FP4(this.c).sqr();
        const D = new FP4(this.a).mul(this.b).dbl();

        this.c.add(this.a).add(this.b).norm().sqr();

        this.a.copy(A);

        A.add(B).add(C).add(D).neg();
        B.times_i();
        C.times_i();

        this.a.add(B);
        this.b.copy(C).add(D);
        this.c.add(A);

        if (this.stype === FP.SPARSER) {
            this.stype = FP.SPARSE;
        } else {
            this.stype = FP.DENSE;
        }

        return this.norm();
    },

    /**
     * Full unconditional Multiplication of two FP12s
     *
     * @this {FP12}
     * @param y FP12 instance, the multiplier
     */
    mul: function(y) {
        const z0 = new FP4(this.a).mul(y.a);
        const z2 = new FP4(this.b).mul(y.b);

        const t0 = new FP4(this.a).add(this.b).norm();
        const t1 = new FP4(y.a).add(y.b).norm();

        const z1 = new FP4(t0).mul(t1);
        t0.copy(this.b).add(this.c).norm();
        t1.copy(y.b).add(y.c).norm();

        const z3 = new FP4(t0).mul(t1);
        t0.copy(z0).neg();
        t1.copy(z2).neg();

        z1.add(t0);
        this.b.copy(z1).add(t1);

        z3.add(t1);
        z2.add(t0);

        t1.copy(y.a).add(y.c).norm();
        t0.copy(this.a).add(this.c).norm().mul(t1);
        z2.add(t0);

        t0.copy(this.c).mul(y.c);
        t1.copy(t0).neg();

        this.c.copy(z2);
        this.c.add(t1);
        z3.add(t1);
        t0.times_i();
        z3.times_i();
        this.b.add(t0);
        this.a.copy(z0).add(z3);
        this.stype = FP.DENSE;
        return this.norm();
    },

    /* FP12 multiplication w=w*y */
    /* catering for special case that arises from special form of ATE pairing line function */
    /* w and y are both sparser line functions - cost = 6m */

    /**
     * Fast multiplication of two sparse FP12s that arises from ATE pairing line functions
     *
     * @this {FP12}
     * @param y FP12 instance, the multiplier
     */
    smul: function(y) {
        const w1 = new FP2(this.a.geta()).mul(y.a.geta());
        const w2 = new FP2(this.a.getb()).mul(y.a.getb());
        const w3 = new FP2(this.b.geta()).mul(y.b.geta());

        const ta = new FP2(this.a.geta()).add(this.a.getb()).norm();
        const tb = new FP2(y.a.geta()).add(y.a.getb()).norm();
        const t  = new FP2(w1).add(w2).neg();
        const tc = new FP2(ta).mul(tb).add(t);

        ta.copy(this.a.geta()).add(this.b.geta()).norm();
        tb.copy(y.a.geta()).add(y.b.geta()).norm();
        t.copy(w1).add(w3).neg();
        const td = new FP2(ta).mul(tb).add(t);

        ta.copy(this.a.getb()).add(this.b.geta()).norm();
        tb.copy(y.a.getb()).add(y.b.geta()).norm();
        t.copy(w2).add(w3).neg();
        const te = new FP2(ta).mul(tb).add(t);

        w2.mul_ip();
        w1.add(w2);

        this.a.geta().copy(w1);
        this.a.getb().copy(tc);
        this.b.geta().copy(td);
        this.b.getb().copy(te);
        this.c.geta().copy(w3);
        this.c.getb().zero();

        this.a.norm();
        this.b.norm();

        this.stype = FP.SPARSE;
    },

    /* FP12 full multiplication w=w*y */
    /* Supports sparse multiplicands */
    /* Usually w is denser than y */

    /**
     * Fast multiplication of what may be sparse multiplicands
     *
     * @this {FP12}
     * @param y FP12 instance, the multiplier
     */
    ssmul: function(y) {
        if (this.stype === FP.ONE) {
            this.copy(y);
            return;
        }
        if (y.stype === FP.ONE) {
            return;
        }

        if (y.stype >= FP.SPARSE) {
            const z0 = new FP4(this.a).mul(y.a);
            const z2 = new FP4(this.b).mul(y.b);
            const t0 = new FP4(this.a).add(this.b).norm();
            const t1 = new FP4(y.a).add(y.b).norm();
            const z1 = new FP4(t0).mul(t1);

            t0.copy(this.b).add(this.c).norm();
            t1.copy(y.b).add(y.c).norm();
            const z3 = new FP4(t0).mul(t1);

            t0.copy(z0).neg();
            t1.copy(z2).neg();
            z1.add(t0).add(t1);
            z2.add(t0);
            z3.add(t1);

            this.b.copy(z1);

            t1.copy(y.a).add(y.c).norm();
            t0.copy(this.a).add(this.c).norm().mul(t1);
            z2.add(t0);

            if (y.stype === FP.SPARSE || this.stype === FP.SPARSE) {
                t0.geta().copy(this.c.geta()).mul(y.c.geta());
                t0.getb().zero();
                if (y.stype !== FP.SPARSE) {
                    t0.getb().copy(this.c.geta()).mul(y.c.getb());
                }
                if (this.stype !== FP.SPARSE) {
                    t0.getb().copy(this.c.getb()).mul(y.c.geta());
                }
            } else {
                t0.copy(this.c).mul(y.c);
            }
            t1.copy(t0).neg();

            this.c.copy(z2).add(t1);
            t0.times_i();
            this.b.add(t0);

            z3.add(t1).norm().times_i();
            this.a.copy(z0).add(z3);
        } else {
            if (this.stype === FP.SPARSER) {
                this.smul(y);
                return;
            }

            // dense by sparser - 13m
            const z0 = new FP4(this.a).mul(y.a);
            const z2 = new FP4(this.b).pmul(y.b.real());
            const z3 = new FP4(this.b).add(this.c).norm().pmul(y.b.real());

            const t1 = new FP4(y.a);
            t1.real().add(y.b.real());
            t1.norm();

            this.b.add(this.a).norm().mul(t1).sub(z0).sub(z2);
            z3.sub(z2).norm().times_i().add(z0);

            const t0 = new FP4(this.a).add(this.c).norm().mul(y.a);
            this.a.copy(z3);
            this.c.copy(z2).sub(z0).add(t0);
        }

        this.stype = FP.DENSE;
        this.norm();
    },

    /**
     * Inverting an FP12
     *
     * @this {FP12}
     */
    inverse: function() {
        const f0 = new FP4(this.a);
        const f1 = new FP4(this.b);
        const f2 = new FP4(this.a);
        const f3 = new FP4(0);

        f0.sqr();
        f1.mul(this.c);
        f1.times_i();
        f0.sub(f1);
        f0.norm();

        f1.copy(this.c);
        f1.sqr();
        f1.times_i();
        f2.mul(this.b);
        f1.sub(f2);
        f1.norm();

        f2.copy(this.b);
        f2.sqr();
        f3.copy(this.a);
        f3.mul(this.c);
        f2.sub(f3);
        f2.norm();

        f3.copy(this.b);
        f3.mul(f2);
        f3.times_i();
        this.a.mul(f0);
        f3.add(this.a);
        this.c.mul(f1);
        this.c.times_i();

        f3.add(this.c);
        f3.norm();
        f3.inverse();
        this.a.copy(f0);
        this.a.mul(f3);
        this.b.copy(f1);
        this.b.mul(f3);
        this.c.copy(f2);
        this.c.mul(f3);
        this.stype = FP.DENSE;

        return this;
    },

    /**
     * Raises an FP12 to the power of the internal modulus p, using the Frobenius
     *
     * @this {FP12}
     * @param f Modulus
     */
    frob: function(f) {
        const f2 = new FP2(f).sqr();
        const f3 = new FP2(f).mul(f2);

        this.a.frob(f3);
        this.b.frob(f3).pmul(f);
        this.c.frob(f3).pmul(f2);

        this.stype = FP.DENSE;

        return this;
    },

    /**
     * Calculate the trace of an FP12
     *
     * @this {FP12}
     */
    trace: function() {
        return new FP4(this.a).imul(3).reduce();
    },

    /**
     * convert this to hex string
     *
     * @this {FP12}
     */
    toString: function() {
        return ("[" + this.a.toString() + "," + this.b.toString() + "," + this.c.toString() + "]");
    },

    /**
     * convert this to byte array
     *
     * @this {FP12}
     * @param w Byte array
     */
    toBytes: function(w) {
        const t = [];
        this.a.geta().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i] = t[i];
        }
        this.a.geta().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + BIG.MODBYTES] = t[i];
        }
        this.a.getb().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 2 * BIG.MODBYTES] = t[i];
        }
        this.a.getb().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 3 * BIG.MODBYTES] = t[i];
        }

        this.b.geta().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 4 * BIG.MODBYTES] = t[i];
        }
        this.b.geta().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 5 * BIG.MODBYTES] = t[i];
        }
        this.b.getb().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 6 * BIG.MODBYTES] = t[i];
        }
        this.b.getb().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 7 * BIG.MODBYTES] = t[i];
        }

        this.c.geta().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 8 * BIG.MODBYTES] = t[i];
        }
        this.c.geta().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 9 * BIG.MODBYTES] = t[i];
        }
        this.c.getb().getA().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 10 * BIG.MODBYTES] = t[i];
        }
        this.c.getb().getB().toBytes(t);
        for (let i = 0; i < BIG.MODBYTES; i++) {
            w[i + 11 * BIG.MODBYTES] = t[i];
        }
    },

    /**
     * Raises an FP12 to the power of a BIG
     *
     * @this {FP12}
     * @param e BIG instance exponent
     */
    pow: function(e) {
        const e1 = new BIG(e).norm();
        const e3 = new BIG(e1).imul(3).norm();
        const sf = new FP12(this).norm();
        const w = new FP12(sf);

        const nb = e3.nbits();
        for (let i = nb - 2; i >= 1; i--) {
            w.usqr();
            const bt = e3.bit(i) - e1.bit(i);

            if (bt === 1) {
                w.mul(sf);
            } else if (bt === -1) {
                sf.conj();
                w.mul(sf);
                sf.conj();
            }
        }

        return w.reduce();
    },

    /**
     * Raises an FP12 instance x to a small integer power, side-channel resistant
     *
     * @this {FP12}
     * @param e small integer exponent
     * @param bts maximum number of bits in exponent
     */
    pinpow: function(e, bts) {
        const R = [];
        R[0] = new FP12(1);
        R[1] = new FP12(this);

        for (let i = bts - 1; i >= 0; i--) {
            const b = (e >> i) & 1;
            R[1 - b].mul(R[b]);
            R[b].usqr();
        }

        this.copy(R[0]);
    },

    /**
     * Raises an FP12 instance to a BIG power, compressed to FP4
     *
     * @this {FP12}
     * @param e BIG exponent
     * @param r BIG group order
     */
    compow: function(e, r) {

        const fa = new BIG(0).rcopy(ROM_FIELD.Fra);
        const fb = new BIG(0).rcopy(ROM_FIELD.Frb);
        const f  = new FP2(fa, fb);

        const q = new BIG(0).rcopy(ROM_FIELD.Modulus);
        const m = new BIG(q).mod(r);

        const a = new BIG(e).mod(m);
        const b = new BIG(e).div(m);

        const g1 = new FP12(this);

        let c = g1.trace();
        if (b.iszilch()) {
            c = c.xtr_pow(e);
            return c;
        }

        const g2 = new FP12(g1).frob(f);

        const cp = g2.trace();
        g1.conj();
        g2.mul(g1);

        const cpm1 = g2.trace();
        g2.mul(g1);

        const cpm2 = g2.trace();

        c = c.xtr_pow2(cp, cpm1, cpm2, a, b);
        return c;
    }
};

/**
  * convert from byte array to FP12 
  *
  * @this {FP12}
  * @param w Byte array
  */    
FP12.fromBytes = function(w) {
    const t = [];
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i];
    }
    let a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + BIG.MODBYTES];
    }
    let b = BIG.fromBytes(t);
    let c = new FP2(a, b);

    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 2 * BIG.MODBYTES];
    }
    a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 3 * BIG.MODBYTES];
    }
    b = BIG.fromBytes(t);
    let d = new FP2(a, b);

    const e = new FP4(c, d);

    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 4 * BIG.MODBYTES];
    }
    a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 5 * BIG.MODBYTES];
    }
    b = BIG.fromBytes(t);
    c = new FP2(a, b); 

    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 6 * BIG.MODBYTES];
    }
    a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 7 * BIG.MODBYTES];
    }
    b = BIG.fromBytes(t);
    d = new FP2(a, b);

    const f = new FP4(c, d);

    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 8 * BIG.MODBYTES];
    }
    a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 9 * BIG.MODBYTES];
    }
    b = BIG.fromBytes(t);
    c = new FP2(a, b); 

    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 10 * BIG.MODBYTES];
    }
    a = BIG.fromBytes(t);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        t[i] = w[i + 11 * BIG.MODBYTES];
    }
    b = BIG.fromBytes(t);
    d = new FP2(a, b); 

    const g = new FP4(c, d);

    return new FP12(e, f, g);
};


/**
  * return 1 if b==c, no branching 
  *
  * @this {FP12}
  */
FP12.teq = function(b, c) {
    let x = b ^ c;
    x -= 1; // if x=0, x now -1
    return ((x >> 31) & 1);
};

/* p=q0^u0.q1^u1.q2^u2.q3^u3 */
// Bos & Costello https://eprint.iacr.org/2013/458.pdf
// Faz-Hernandez & Longa & Sanchez  https://eprint.iacr.org/2013/158.pdf
// Side channel attack secure

/**
  * p=q0^u0.q1^u1.q2^u2.q3^u3 
  *
  * @this {FP12}
  */    
FP12.pow4 = function(q, u) {

    const t = [];
    for (let i = 0; i < 4; i++) {
        t[i] = new BIG(u[i]).norm();
    }

    const g = [];
    g[0] = new FP12(q[0]);  // q[0]
    g[1] = new FP12(g[0]); g[1].mul(q[1]);  // q[0].q[1]
    g[2] = new FP12(g[0]); g[2].mul(q[2]);  // q[0].q[2]
    g[3] = new FP12(g[1]); g[3].mul(q[2]);  // q[0].q[1].q[2]
    g[4] = new FP12(q[0]); g[4].mul(q[3]);  // q[0].q[3]
    g[5] = new FP12(g[1]); g[5].mul(q[3]);  // q[0].q[1].q[3]
    g[6] = new FP12(g[2]); g[6].mul(q[3]);  // q[0].q[2].q[3]
    g[7] = new FP12(g[3]); g[7].mul(q[3]);  // q[0].q[1].q[2].q[3]

    // Make it odd
    const pb = 1 - t[0].parity();
    t[0].inc(pb).norm();

    // Number of bits
    const mt = new BIG(0);
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

    const p = new FP12(0);
    const r = new FP12(0);

    p.select(g, 2 * w[nb - 1] + 1);
    for (let i = nb - 2; i >= 0; i--) {
        p.usqr();
        r.select(g, 2 * w[i] + s[i]);
        p.mul(r);
    }

    // apply correction
    r.copy(q[0]);
    r.conj();
    r.mul(p);
    p.cmove(r, pb);

    p.reduce();
    return p;
};

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = FP12;
}
