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
const DBIG = require("./big").DBIG;
const FP = require("./fp");
const ROM_FIELD = require("./rom_field");

/* Finite Field arithmetic  Fp^2 functions */

/* FP2 elements are of the form a+ib, where i is sqrt(-1) */

"use strict";

/**
  * Creates an instance of FP2.
  *
  * @constructor
  * @this {FP2}
  */
const FP2 = function(c, d) {
    if (c instanceof FP2) {
        this.a = new FP(c.a);
        this.b = new FP(c.b);
    } else {
        this.a = new FP(c);
        this.b = new FP(d);
    }
};

FP2.prototype = {

    /**
     * Reduces all components of possibly unreduced FP2 mod Modulus
     *
     * @this {FP2}
     */
    reduce: function() {
        this.a.reduce();
        this.b.reduce();
    },

    /**
     * Normalises the components of an FP2
     *
     * @this {FP2}
     */
    norm: function() {
        this.a.norm();
        this.b.norm();
    },

    /**
     * Tests for FP2 equal to zero
     *
     * @this {FP2}
     */
    iszilch: function() {
        return (this.a.iszilch() && this.b.iszilch());
    },

    /**
     * Tests for FP2 equal to unity
     *
     * @this {FP2}
     */
    isunity: function() {
        const one = new FP(1);
        return (this.a.equals(one) && this.b.iszilch());
    },

    /**
     * Tests for FP2 being odd
     */
    isodd: function() {
        return this.a.isodd();
    },

    /**
     * Conditional copy of FP2 number
     *
     * @this {FP2}
     * @param g FP2 instance
     * @param d copy depends on this value
     */
    cmove: function(g, d) {
        this.a.cmove(g.a, d);
        this.b.cmove(g.b, d);
    },

    /**
     * Tests for equality of two FP2 instances
     *
     * @this {FP2}
     * @param x FP2 instance to compare
     */
    equals: function(x) {
        return (this.a.equals(x.a) && this.b.equals(x.b));
    },

    /**
     * extract a from this
     *
     * @this {FP2}
     */
    getA: function() {
        return this.a.redc();
    },

    /**
     * extract b from this
     *
     * @this {FP2}
     */
    getB: function() {
        return this.b.redc();
    },

    /**
     * Set FP2 from two FP values
     *
     * @this {FP2}
     * @param c FP instance
     * @param d FP instance
     */
    set: function(c, d) {
        this.a.copy(c);
        this.b.copy(d);
    },

    /**
     * Set FP2 from one FP value
     *
     * @this {FP2}
     * @param c FP instance
     */
    seta: function(c) {
        this.a.copy(c);
        this.b.zero();
    },

    /**
     * Set FP2 from two BIG values
     *
     * @this {FP2}
     * @param c BIG instance
     * @param d BIG instance
     */
    bset: function(c, d) {
        this.a.bcopy(c);
        this.b.bcopy(d);
    },

    /**
     * Set FP2 from one BIG value
     *
     * @this {FP2}
     * @param c BIG instance
     */
    bseta: function(c) {
        this.a.bcopy(c);
        this.b.zero();
    },

    /**
     * Copy FP2 to another FP2
     *
     * @this {FP2}
     * @param x FP2 instance to be copied
     */
    copy: function(x) {
        this.a.copy(x.a);
        this.b.copy(x.b);
    },

    /**
     * Set FP2 to zero
     *
     * @this {FP2}
     */
    zero: function() {
        this.a.zero();
        this.b.zero();
    },

    /**
     * Set FP2 to unity
     *
     * @this {FP2}
     * @param x FP2 instance to be set to one
     */
    one: function() {
        this.a.one();
        this.b.zero();
    },

    /**
     * negate this
     *
     * @this {FP2}
     */
    neg: function() {
        const m = new FP(this.a);
        m.add(this.b);
        m.neg();
        const t = new FP(0);
        t.copy(m);
        t.add(this.b);
        this.b.copy(m);
        this.b.add(this.a);
        this.a.copy(t);
    },

    /**
     * Conjugation of FP2
     *
     * @this {FP2}
     */
    conj: function() {
        this.b.neg();
        this.b.norm();
    },

    /**
     * addition of two FP2s
     *
     * @this {FP2}
     * @param x FP2 instance
     */
    add: function(x) {
        this.a.add(x.a);
        this.b.add(x.b);
    },

    /**
     * subtraction of two FP2s
     *
     * @this {FP2}
     * @param x FP2 instance
     */
    sub: function(x) {
        const m = new FP2(x);
        m.neg();
        this.add(m);
    },

    rsub: function(x) {
        this.neg();
        this.add(x);
    },

    /**
     * Multiplication of an FP2 by an FP
     *
     * @this {FP2}
     * @param s FP instance
     */
    pmul: function(s) {
        this.a.mul(s);
        this.b.mul(s);
    },

    /**
     * Multiplication of an FP2 by a small integer
     *
     * @this {FP2}
     * @param c integer
     */
    imul: function(c) {
        this.a.imul(c);
        this.b.imul(c);
    },

    /**
     * Fast Squaring of an FP2
     *
     * @this {FP2}
     */
    sqr: function() {
        const w1 = new FP(this.a);
        const w3 = new FP(this.a);
        const mb = new FP(this.b);

        w1.add(this.b);
        w3.add(this.a);
        w3.norm();
        this.b.mul(w3);

        mb.neg();
        this.a.add(mb);

        this.a.norm();
        w1.norm();

        this.a.mul(w1);
    },

    /**
     * Full unconditional Multiplication of two FP2s
     *
     * @this {FP2}
     * @param y FP2 instance, the multiplier
     */
    mul: function(y) {
        const p = new BIG(0);
        p.rcopy(ROM_FIELD.Modulus);

        const pR = new DBIG(0);
        pR.ucopy(p);

        if ((this.a.XES + this.b.XES) * (y.a.XES + y.b.XES) > FP.FEXCESS) {
            if (this.a.XES > 1) {
                this.a.reduce();
            }

            if (this.b.XES > 1) {
                this.b.reduce();
            }
        }

        const A = BIG.mul(this.a.f, y.a.f);
        const B = BIG.mul(this.b.f, y.b.f);

        const C = new BIG(this.a.f);
        const D = new BIG(y.a.f);

        C.add(this.b.f);
        C.norm();
        D.add(y.b.f);
        D.norm();

        const E = BIG.mul(C, D);
        const F = new DBIG(0);
        F.copy(A);
        F.add(B);
        B.rsub(pR);

        A.add(B);
        A.norm();
        E.sub(F);
        E.norm();

        this.a.f.copy(FP.mod(A));
        this.a.XES = 3;
        this.b.f.copy(FP.mod(E));
        this.b.XES = 2;
    },

    /**
     * sqrt(a+ib) = sqrt(a+sqrt(a*a-n*b*b)/2)+ib/(2*sqrt(a+sqrt(a*a-n*b*b)/2))
     *
     * @this {FP2}
     * @return true if this is QR
     */
    sqrt: function() {

        if (this.iszilch()) {
            return true;
        }

        let w1 = new FP(this.b);
        let w2 = new FP(this.a);

        w1.sqr();
        w2.sqr();
        w1.add(w2);
        if (w1.jacobi() !== 1) {
            this.zero();
            return false;
        }
        w1 = w1.sqrt();
        w2.copy(this.a);
        w2.add(w1);
        w2.norm();
        w2.div2();
        if (w2.jacobi() !== 1) {
            w2.copy(this.a);
            w2.sub(w1);
            w2.norm();
            w2.div2();
            if (w2.jacobi() !== 1) {
                this.zero();
                return false;
            }
        }
        w2 = w2.sqrt();
        this.a.copy(w2);
        w2.add(w2);
        w2.inverse();
        this.b.mul(w2);

        return true;
    },

    /**
     * convert this to hex string
     *
     * @this {FP2}
     */
    toString: function() {
        return ("[" + this.a.toString() + "," + this.b.toString() + "]");
    },

    /**
     * Inverting an FP2
     *
     * @this {FP2}
     */
    inverse: function() {

        this.norm();
        const w1 = new FP(this.a);
        const w2 = new FP(this.b);

        w1.sqr();
        w2.sqr();
        w1.add(w2);
        w1.inverse();
        this.a.mul(w1);
        w1.neg();
        w1.norm();
        this.b.mul(w1);
    },

    /**
     * Divide an FP2 by 2
     *
     * @this {FP2}
     */
    div2: function() {
        this.a.div2();
        this.b.div2();
    },

    /**
     * Multiply an FP2 by sqrt(-1)
     *
     * @this {FP2}
     */
    times_i: function() {
        const z = new FP(this.a);
        this.a.copy(this.b);
        this.a.neg();
        this.b.copy(z);
    },

    /**
     * Multiply an FP2 by (1+sqrt(-1))
     *
     * @this {FP2}
     */
    mul_ip: function() {
        const t = new FP2(this);
        const z = new FP(this.a);

        this.a.copy(this.b);
        this.a.neg();
        this.b.copy(z);
        this.add(t);
    },

    /**
     * Divide an FP2 by (1+sqrt(-1))/2
     *
     * @this {FP2}
     */
    div_ip2: function() {
        this.norm();
        const t = new FP2(0);
        t.a.copy(this.a);
        t.a.add(this.b);
        t.b.copy(this.b);
        t.b.sub(this.a);
        this.copy(t);
        this.norm();
    },

    /**
     * Divide an FP2 by (1+sqrt(-1))
     *
     * @this {FP2}
     */
    div_ip: function() {
        this.norm();
        const t = new FP2(0);
        t.a.copy(this.a);
        t.a.add(this.b);
        t.b.copy(this.b);
        t.b.sub(this.a);
        this.copy(t);
        this.norm();
        this.div2();
    },

    /**
     * Raises an FP2 to the power of a BIG
     *
     * @this {FP2}
     * @param e BIG instance exponent
     */
    pow: function(e) {
        this.norm();
        e.norm();

        const x = new FP2(this);
        const r = new FP2(1);

        for (;;) {
            const bt = e.parity();
            e.fshr(1);

            if (bt === 1) {
                r.mul(x);
            }

            if (e.iszilch()) {
                break;
            }
            x.sqr();
        }

        r.reduce();
        return r;
    }
};

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = FP2;
}
