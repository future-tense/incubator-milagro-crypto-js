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

/* Finite Field arithmetic  Fp^4 functions */

/* FP4 elements are of the form a+ib, where i is sqrt(-1+sqrt(-1))  */

"use strict";

/**
  * Creates an instance of FP4
  *
  * @constructor
  * @this {FP4}
  */
const FP4 = function(c, d) {
    if (c instanceof FP4) {
        this.a = new FP2(c.a);
        this.b = new FP2(c.b);
    } else {
        this.a = new FP2(c);
        this.b = new FP2(d);
    }
};

FP4.prototype = {

    /**
     * Reduces all components of possibly unreduced FP4 mod Modulus
     *
     * @this {FP4}
     */
    reduce: function() {
        this.a.reduce();
        this.b.reduce();
        return this;
    },

    /**
     * Normalises the components of an FP4
     *
     * @this {FP4}
     */
    norm: function() {
        this.a.norm();
        this.b.norm();
        return this;
    },

    /**
     * Tests for FP4 equal to zero
     *
     * @this {FP4}
     */
    iszilch: function() {
        return (this.a.iszilch() && this.b.iszilch());
    },

    /**
     * Tests for FP4 equal to unity
     *
     * @this {FP4}
     */
    isunity: function() {
        const one = new FP2(1);
        return (this.a.equals(one) && this.b.iszilch());
    },

    /**
     * Conditional copy of FP4 number
     *
     * @this {FP4}
     * @param g FP4 instance
     * @param d copy depends on this value
     */
    cmove: function(g, d) {
        this.a.cmove(g.a, d);
        this.b.cmove(g.b, d);
    },

    /**
     * test is w real? That is in a+ib test b is zero 
     *
     * @this {FP4}
     */
    isreal: function() {
        return this.b.iszilch();
    },

    /**
     * extract real part a
     *
     * @this {FP4}
     */
    real: function() {
        return this.a;
    },

    /**
     * extract a from this
     *
     * @this {FP4}
     */		
    geta: function() {
        return this.a;
    },

    /**
     * extract b from this
     *
     * @this {FP4}
     */	
    getb: function() {
        return this.b;
    },

    /**
     * Tests for equality of two FP4s
     *
     * @this {FP4}
     * @param x FP4 instance to compare
     */
    equals: function(x) {
        return (this.a.equals(x.a) && this.b.equals(x.b));
    },

    /**
     * Copy FP4 to another FP4
     *
     * @this {FP4}
     * @param x FP4 instance to be copied
     */
    copy: function(x) {
        this.a.copy(x.a);
        this.b.copy(x.b);
        return this;
    },

    /**
     * Set FP4 to zero
     *
     * @this {FP4}
     */
    zero: function() {
        this.a.zero();
        this.b.zero();
    },

    /**
     * Set FP4 to unity
     *
     * @this {FP4}
     * @param x FP4 instance to be set to one
     */
    one: function() {
        this.a.one();
        this.b.zero();
    },

    /**
     * Set FP4 from two FP2 values
     *
     * @this {FP4}
     * @param c FP2 instance
     * @param d FP2 instance
     */
    set: function(c, d) {
        this.a.copy(c);
        this.b.copy(d);
    },

    /**
     * Set FP4 from one FP2 value
     *
     * @this {FP4}
     * @param c FP2 instance
     */
    seta: function(c) {
        this.a.copy(c);
        this.b.zero();
    },

    /**
     * Negation of FP4
     *
     * @this {FP4}
     */
    neg: function() {
        this.norm();
        const m = new FP2(this.a).add(this.b).neg();
        const t = new FP2(m).add(this.b);
        this.b.copy(m).add(this.a);
        this.a.copy(t);
        return this.norm();
    },

    /**
     * Conjugation of FP4
     *
     * @this {FP4}
     */
    conj: function() {
        this.b.neg();
        return this.norm();
    },

    /**
     * Negative conjugation of FP4
     *
     * @this {FP4}
     */
    nconj: function() {
        this.a.neg();
        return this.norm();
    },

    dbl: function() {
        this.a.add(this.a);
        this.b.add(this.b);
        return this;
    },

    /**
     * addition of two FP4s
     *
     * @this {FP4}
     * @param x FP4 instance
     */
    add: function(x) {
        this.a.add(x.a);
        this.b.add(x.b);
        return this;
    },

    /**
     * subtraction of two FP4s
     *
     * @this {FP4}
     * @param x FP4 instance
     */
    sub: function(x) {
        const m = new FP4(x).neg();
        this.add(m);
        return this;
    },

    rsub: function(x) {
        this.neg();
        this.add(x);
    },

    /**
     * Multiplication of an FP4 by an FP
     *
     * @this {FP4}
     * @param s FP instance
     */
    pmul: function(s) {
        this.a.mul(s);
        this.b.mul(s);
        return this;
    },

    /**
     * Multiplication of an FP4 by an integer
     *
     * @this {FP4}
     * @param s integer multiplier
     */
    imul: function(c) {
        this.a.imul(c);
        this.b.imul(c);
    },

    /**
     * Fast Squaring of an FP4
     *
     * @this {FP4}
     */
    sqr: function() {
        const t1 = new FP2(this.a).add(this.b).norm();
        const t2 = new FP2(this.b).mul_ip().add(this.a).norm();
        const t3 = new FP2(this.a).mul(this.b);

        this.a.copy(t1).mul(t2);
        t2.copy(t3).mul_ip().add(t3).norm();

        this.a.sub(t2);
        this.b.copy(t3).add(t3);
        return this.norm();
    },

    /**
     * Full unconditional Multiplication of two FP4s
     *
     * @this {FP4}
     * @param y FP4 instance, the multiplier
     */
    mul: function(y) {
        const t1 = new FP2(this.a).mul(y.a);
        const t2 = new FP2(this.b).mul(y.b);
        const t3 = new FP2(y.b).add(y.a).norm();
        const t4 = new FP2(this.b).add(this.a).norm().mul(t3).sub(t1);

        t3.copy(t2);
        t2.mul_ip();

        this.a.copy(t2).add(t1);
        this.b.copy(t4).sub(t3);
        return this.norm();
    },

    /**
     * convert to hex string
     *
     * @this {FP4}
     */
    toString: function() {
        return ("[" + this.a.toString() + "," + this.b.toString() + "]");
    },

    /**
     * Inverting an FP4
     *
     * @this {FP4}
     */
    inverse: function() {
        this.norm();
        const t1 = new FP2(this.a).sqr();
        const t2 = new FP2(this.b).sqr().mul_ip().norm();
        t1.sub(t2).inverse();
        this.a.mul(t1);
        t1.neg().norm();
        this.b.mul(t1);
    },

    /**
     * multiplies an FP4 instance by irreducible polynomial sqrt(1+sqrt(-1))
     *
     * @this {FP4}
     */
    times_i: function() {
        const s = new FP2(this.b).times_i();
        const t = new FP2(this.b).add(s);
        this.b.copy(this.a);
        this.a.copy(t);
        return this.norm();
    },

    /**
     * Raises an FP4 to the power of the internal modulus p, using the Frobenius
     *
     * @this {FP4}
     * @param f Modulus
     */
    frob: function(f) {
        this.a.conj();
        this.b.conj();
        this.b.mul(f);
        return this;
    },

    /**
     * Raises an FP4 to the power of a BIG
     *
     * @this {FP4}
     * @param e BIG instance exponent
     */
    pow: function(e) {
        const w = new FP4(this).norm();
        const z = new BIG(e).norm();

        const r = new FP4(1);
        for (;;) {
            const bt = z.parity();
            z.fshr(1);

            if (bt === 1) {
                r.mul(w);
            }

            if (z.iszilch()) {
                break;
            }

            w.sqr();
        }

        return r.reduce();
    },

    /**
     * Calculates the XTR addition function r=w*x-conj(x)*y+z
     *
     * @this {FP4}
     * @param w FP4 instance
     * @param y FP4 instance
     * @param z FP4 instance
     */
    xtr_A: function(w, y, z) {
        const r = new FP4(w).sub(y).norm().pmul(this.a);
        const t = new FP4(w).add(y).norm().pmul(this.b).times_i();
        return this.copy(r).add(t).add(z).reduce();
    },

    /**
     * Calculates the XTR doubling function r=x^2-2*conj(x)
     *
     * @this {FP4}
     */
    xtr_D: function() {
        const w = new FP4(this).conj();
        w.add(w);
        this.sqr().sub(w).reduce();
    },

    /**
     * Calculates FP4 trace of an FP4 raised to the power of a BIG number
     *
     * @this {FP4}
     * @param n Big number
     */
    xtr_pow: function(n) {
        const sf = new FP4(this);
        sf.norm();
        const a = new FP4(3);
        const b = new FP4(sf);
        const c = new FP4(b);
        const t = new FP4(0);
        const r = new FP4(0);

        c.xtr_D();

        const par = n.parity();
        const v = new BIG(n);
        v.norm();
        v.fshr(1);

        if (par === 0) {
            v.dec(1);
            v.norm();
        }

        const nb = v.nbits();
        for (let i = nb - 1; i >= 0; i--) {
            if (v.bit(i) !== 1) {
                t.copy(b);
                sf.conj();
                c.conj();
                b.xtr_A(a, sf, c);
                sf.conj();
                c.copy(t);
                c.xtr_D();
                a.xtr_D();
            } else {
                t.copy(a);
                t.conj();
                a.copy(b);
                a.xtr_D();
                b.xtr_A(c, sf, t);
                c.xtr_D();
            }
        }

        if (par === 0) {
            r.copy(c);
        } else {
            r.copy(b);
        }
        r.reduce();

        return r;
    },

    /**
     * Calculates FP4 trace of c^a.d^b, where c and d are derived from FP4 traces of FP4s
     *
     * @this {FP4}
     */
    xtr_pow2: function(ck, ckml, ckm2l, a, b) {

        const e = new BIG(a);
        const d = new BIG(b);
        const w = new BIG(0);
        const cu = new FP4(ck);
        const cv = new FP4(this);
        const cumv = new FP4(ckml);
        const cum2v = new FP4(ckm2l);
        const r = new FP4(0);
        const t = new FP4(0);

        e.norm();
        d.norm();

        let f2 = 0;
        while (d.parity() === 0 && e.parity() === 0) {
            d.fshr(1);
            e.fshr(1);
            f2++;
        }

        while (BIG.comp(d, e) !== 0) {
            if (BIG.comp(d, e) > 0) {
                w.copy(e);
                w.imul(4);
                w.norm();

                if (BIG.comp(d, w) <= 0) {
                    w.copy(d);
                    d.copy(e);
                    e.rsub(w);
                    e.norm();

                    t.copy(cv);
                    t.xtr_A(cu, cumv, cum2v);
                    cum2v.copy(cumv);
                    cum2v.conj();
                    cumv.copy(cv);
                    cv.copy(cu);
                    cu.copy(t);

                } else if (d.parity() === 0) {
                    d.fshr(1);
                    r.copy(cum2v);
                    r.conj();
                    t.copy(cumv);
                    t.xtr_A(cu, cv, r);
                    cum2v.copy(cumv);
                    cum2v.xtr_D();
                    cumv.copy(t);
                    cu.xtr_D();
                } else if (e.parity() === 1) {
                    d.sub(e);
                    d.norm();
                    d.fshr(1);
                    t.copy(cv);
                    t.xtr_A(cu, cumv, cum2v);
                    cu.xtr_D();
                    cum2v.copy(cv);
                    cum2v.xtr_D();
                    cum2v.conj();
                    cv.copy(t);
                } else {
                    w.copy(d);
                    d.copy(e);
                    d.fshr(1);
                    e.copy(w);
                    t.copy(cumv);
                    t.xtr_D();
                    cumv.copy(cum2v);
                    cumv.conj();
                    cum2v.copy(t);
                    cum2v.conj();
                    t.copy(cv);
                    t.xtr_D();
                    cv.copy(cu);
                    cu.copy(t);
                }
            }
            if (BIG.comp(d, e) < 0) {
                w.copy(d);
                w.imul(4);
                w.norm();

                if (BIG.comp(e, w) <= 0) {
                    e.sub(d);
                    e.norm();
                    t.copy(cv);
                    t.xtr_A(cu, cumv, cum2v);
                    cum2v.copy(cumv);
                    cumv.copy(cu);
                    cu.copy(t);
                } else if (e.parity() === 0) {
                    w.copy(d);
                    d.copy(e);
                    d.fshr(1);
                    e.copy(w);
                    t.copy(cumv);
                    t.xtr_D();
                    cumv.copy(cum2v);
                    cumv.conj();
                    cum2v.copy(t);
                    cum2v.conj();
                    t.copy(cv);
                    t.xtr_D();
                    cv.copy(cu);
                    cu.copy(t);
                } else if (d.parity() === 1) {
                    w.copy(e);
                    e.copy(d);
                    w.sub(d);
                    w.norm();
                    d.copy(w);
                    d.fshr(1);
                    t.copy(cv);
                    t.xtr_A(cu, cumv, cum2v);
                    cumv.conj();
                    cum2v.copy(cu);
                    cum2v.xtr_D();
                    cum2v.conj();
                    cu.copy(cv);
                    cu.xtr_D();
                    cv.copy(t);
                } else {
                    d.fshr(1);
                    r.copy(cum2v);
                    r.conj();
                    t.copy(cumv);
                    t.xtr_A(cu, cv, r);
                    cum2v.copy(cumv);
                    cum2v.xtr_D();
                    cumv.copy(t);
                    cu.xtr_D();
                }
            }
        }
        r.copy(cv);
        r.xtr_A(cu, cumv, cum2v);
        for (let i = 0; i < f2; i++) {
            r.xtr_D();
        }
        r = r.xtr_pow(d);
        return r;
    },

    /* New stuff for ecp4.js */

    /**
     * Divide an FP4 by 2
     *
     * @this {FP4}
     */	
    div2: function() {
        this.a.div2();
        this.b.div2();
    },

    /**
     * Divide FP4 number by QNR
     *
     * @this {FP4}
     */		
    div_i: function() {
        const u = new FP2(this.a);
        u.div_ip();

        const v = new FP2(this.b);
        this.a.copy(v);
        this.b.copy(u);
    },

    /**
     * Divide an FP4 by QNR/2
     *
     * @this {FP4}
     */			
    div_2i: function() {
        const u = new FP2(this.a);
        u.div_ip2();

        const v = new FP2(this.b);
        v.add(v);
        v.norm();

        this.a.copy(v);
        this.b.copy(u);
    },

    /**
     * Calculate square root of an FP4
     *
     * @this {FP4}
     */					
    sqrt: function() {
        if (this.iszilch()) {
            return true;
        }
        const wa = new FP2(this.a);
        const ws = new FP2(this.b);
        const wt = new FP2(this.a);

        if (ws.iszilch()) {
            if (wt.sqrt()) {
                this.a.copy(wt);
                this.b.zero();
            } else {
                wt.div_ip();
                wt.sqrt();
                this.b.copy(wt);
                this.a.zero();
            }
            return true;
        }

        ws.sqr();
        wa.sqr();
        ws.mul_ip();
        ws.norm();
        wa.sub(ws);

        ws.copy(wa);
        if (!ws.sqrt()) {
            return false;
        }

        wa.copy(wt);
        wa.add(ws);
        wa.norm();
        wa.div2();

        if (!wa.sqrt()) {
            wa.copy(wt);
            wa.sub(ws);
            wa.norm();
            wa.div2();
            if (!wa.sqrt()) {
                return false;
            }
        }
        wt.copy(this.b);
        ws.copy(wa);
        ws.add(wa);
        ws.inverse();

        wt.mul(ws);
        this.a.copy(wa);
        this.b.copy(wt);

        return true;
    }
};

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = FP4;
}
