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

/* Finite Field arithmetic */
/* AMCL mod p functions */

var BIG = require("./big").BIG;
var DBIG = require("./big").DBIG;
var ROM_FIELD = require("./rom_field");

"use strict";

/**
  * Creates an instance of FP.
  *
  * @constructor
  * @this {FP}
  * @param x FP / BIG instance
  */
var FP = function(x) {
    if (x instanceof FP) {
        this.f = new BIG(x.f);
        this.XES = x.XES;
    } else {
        this.f = new BIG(x);
        this.XES = 1;
        if (!this.f.iszilch()) {
            this.nres();
        }
    }
};

FP.ZERO = 0;
FP.ONE = 1;
FP.SPARSER = 2;
FP.SPARSE = 3;
FP.DENSE= 4;

FP.MODBITS = 254;

FP.FEXCESS = ((1 << 10) - 1); // 2^(BASEBITS*NLEN-MODBITS)-1
FP.OMASK = (-1) << FP.TBITS;
FP.TBITS = FP.MODBITS % BIG.BASEBITS;
FP.TMASK = (1 << FP.TBITS) - 1;

FP.prototype = {

    /**
     * Set FP to zero
     *
     * @this {FP}
     */
    zero: function() {
        this.XES = 1;
        this.f.zero();
    },

    /**
     * copy from a BIG in ROM
     *
     * @this {FP}
     * @param y FP instance to be copied
     */
    rcopy: function(y) {
        this.f.rcopy(y);
        this.nres();
    },

    /**
     * copy from another BIG
     *
     * @this {FP}
     * @param y FP instance to be copied
     */
    bcopy: function(y) {
        this.f.copy(y);
        this.nres();
    },

    /**
     * Copy FP to another FP
     *
     * @this {FP}
     * @param y FP instance to be copied
     */
    copy: function(y) {
        this.XES = y.XES;
        this.f.copy(y.f);
    },

    /**
     * Conditional constant time swap of two FP numbers
     *
     * @this {BIG}
     * @parameter b FP number
     * @parameter d Integer
     */
    cswap: function(b, d) {
        this.f.cswap(b.f, d);
        var t, c = d;
        c = ~(c - 1);
        t = c & (this.XES ^ b.XES);
        this.XES ^= t;
        b.XES ^= t;
    },

    /**
     * Conditional copy of FP number
     *
     * @this {FP}
     * @param b FP instance
     * @param d copy depends on this value
     */
    cmove: function(b, d) {
        var c = d;

        c = ~(c - 1);

        this.f.cmove(b.f, d);
        this.XES ^= (this.XES ^ b.XES) & c;
    },

    /**
     * Converts from BIG integer to residue form mod Modulus
     *
     * @this {FP}
     */
    nres: function() {
        var r, d;

        r = new BIG();
        r.rcopy(ROM_FIELD.R2modp);

        d = BIG.mul(this.f, r);
        this.f.copy(FP.mod(d));
        this.XES = 2;

        return this;
    },

    /**
     * Converts from residue form back to BIG integer form
     *
     * @this {FP}
     */
    redc: function() {
        var r = new BIG(0),
            d, w;

        r.copy(this.f);

        d = new DBIG(0);
        d.hcopy(this.f);
        w = FP.mod(d);
        r.copy(w);
        return r;
    },

    /**
     * convert to hex string
     *
     * @this {FP}
     */
    toString: function() {
        var s = this.redc().toString();
        return s;
    },

    /**
     * Tests for FP equal to zero
     *
     * @this {FP}
     */
    iszilch: function() {
        var c=new FP(0); c.copy(this);
        c.reduce();
        return c.f.iszilch();
    },

    /**
     * Tests for FP being odd
     */
    isodd: function() {
        return !!this.f.parity();
    },

    /**
     * Reduces all components of possibly unreduced FP mod Modulus
     *
     * @this {FP}
     */
    reduce: function() {
        var q,carry,sr,sb,m = new BIG(0);
        m.rcopy(ROM_FIELD.Modulus);
        var r = new BIG(0);
        r.rcopy(ROM_FIELD.Modulus);
        this.f.norm();

        if (this.XES>16) {
            q=FP.quo(this.f,m);
            carry=r.pmul(q);
            r.w[BIG.NLEN-1]+=(carry<<BIG.BASEBITS); // correction - put any carry out back in again
            this.f.sub(r);
            this.f.norm();
            sb=2;
        } else {
            sb=FP.logb2(this.XES-1);
        }
        m.fshl(sb);

        while (sb>0) {
            // constant time...
            sr=BIG.ssn(r,this.f,m);  // optimized combined shift, subtract and norm
            this.f.cmove(r,1-sr);
            sb--;
        }

        this.XES = 1;
    },

    /**
     * Set FP to unity
     *
     * @this {FP}
     */
    one: function() {
        this.f.one();
        this.nres();
    },

    /**
     * Normalises the components of an FP
     *
     * @this {FP}
     */
    norm: function() {
        return this.f.norm();
    },

    /**
     * Fast Modular multiplication of two FPs, mod Modulus
     *
     * @this {FP}
     * @param b FP number, the multiplier
     */
    mul: function(b) {
        var d;

        if (this.XES * b.XES > FP.FEXCESS) {
            this.reduce();
        }

        d = BIG.mul(this.f, b.f);
        this.f.copy(FP.mod(d));
        this.XES = 2;

        return this;
    },

    /**
     * Multiplication of an FP by a small integer
     *
     * @this {FP}
     * @param s integer
     */
    imul: function(c) {
        var s = false,
            d, n;

        if (c < 0) {
            c = -c;
            s = true;
        }

        if (this.XES * c <= FP.FEXCESS) {
            this.f.pmul(c);
            this.XES *= c;
        } else {
            n = new FP(c);
            this.mul(n);
        }

        if (s) {
            this.neg();
            this.norm();
        }
        return this;
    },

    /**
     * Fast Squaring of an FP
     *
     * @this {FP}
     */
    sqr: function() {
        var d, t;

        if (this.XES * this.XES > FP.FEXCESS) {
            this.reduce();
        }

        d = BIG.sqr(this.f);
        t = FP.mod(d);
        this.f.copy(t);
        this.XES = 2;

        return this;
    },

    /* this+=b */
    add: function(b) {
        this.f.add(b.f);
        this.XES += b.XES;

        if (this.XES > FP.FEXCESS) {
            this.reduce();
        }

        return this;
    },

    /**
     * negate this
     *
     * @this {FP}
     * @param x FP instance to be set to one
     */
    neg: function() {
        var m = new BIG(0),
            sb;

        m.rcopy(ROM_FIELD.Modulus);

        sb = FP.logb2(this.XES - 1);

        m.fshl(sb);
        this.XES = (1 << sb)+1;
        this.f.rsub(m);

        if (this.XES > FP.FEXCESS) {
            this.reduce();
        }

        return this;
    },

    /**
     * subtraction of two FPs
     *
     * @this {FP}
     * @param x FP instance
     */
    sub: function(b) {
        var n = new FP(0);

        n.copy(b);
        n.neg();
        this.add(n);

        return this;
    },

    rsub: function(b) {
        var n = new FP(0);

        n.copy(this);
        n.neg();
        this.copy(b);
        this.add(n);
    },

    /**
     * Divide an FP by 2
     *
     * @this {FP}
     */
    div2: function() {
        var p;

        if (this.f.parity() === 0) {
            this.f.fshr(1);
        } else {
            p = new BIG(0);
            p.rcopy(ROM_FIELD.Modulus);

            this.f.add(p);
            this.f.norm();
            this.f.fshr(1);
        }

        return this;
    },

    // return this^(p-3)/4 or this^(p-5)/8
    // See https://eprint.iacr.org/2018/1038

    /**
     * return this^(p-3)/4 or this^(p-5)/8
     *
     * @this {FP}
     */

    fpow: function() {
        var i,j,k,bw,w,c,nw,lo,m,n;
        var xp=[];
        var ac=[1,2,3,6,12,15,30,60,120,240,255];
        // phase 1

        xp[0]=new FP(this);	// 1
        xp[1]=new FP(this); xp[1].sqr(); // 2
        xp[2]=new FP(xp[1]); xp[2].mul(this);  //3
        xp[3]=new FP(xp[2]); xp[3].sqr();  // 6
        xp[4]=new FP(xp[3]); xp[4].sqr();  // 12
        xp[5]=new FP(xp[4]); xp[5].mul(xp[2]);  // 15
        xp[6]=new FP(xp[5]); xp[6].sqr();  // 30
        xp[7]=new FP(xp[6]); xp[7].sqr();  // 60
        xp[8]=new FP(xp[7]); xp[8].sqr();  // 120
        xp[9]=new FP(xp[8]); xp[9].sqr();  // 240
        xp[10]=new FP(xp[9]); xp[10].mul(xp[5]);  // 255


        n=FP.MODBITS;
        n-=2;
        c=(ROM_FIELD.MConst+3)/4;

        bw=0; w=1;
        while (w<c) {
            w *= 2;
            bw += 1;
        }
        k=w-c;

        i=10;
        var key=new FP(0);
        if (k !== 0) {
            while (ac[i]>k) {
                i--;
            }
            key.copy(xp[i]);
            k-=ac[i];
        }
        while (k !== 0) {
            i--;
            if (ac[i]>k) {
                continue;
            }
            key.mul(xp[i]);
            k-=ac[i];
        }

        // phase 2
        xp[1].copy(xp[2]);
        xp[2].copy(xp[5]);
        xp[3].copy(xp[10]);

        j=3; m=8;
        nw=n-bw;
        var t=new FP(0);
        while (2*m<nw) {
            t.copy(xp[j++]);
            for (i=0;i<m;i++) {
                t.sqr();
            }
            xp[j].copy(xp[j-1]);
            xp[j].mul(t);
            m*=2;
        }
        lo=nw-m;
        var r=new FP(xp[j]);

        while (lo!==0) {
            m/=2; j--;
            if (lo<m) {
                continue;
            }
            lo-=m;
            t.copy(r);
            for (i=0;i<m;i++) {
                t.sqr();
            }
            r.copy(t);
            r.mul(xp[j]);
        }

        // phase 3
        if (bw!==0) {
            for (i=0;i<bw;i++ ) {
                r.sqr();
            }
            r.mul(key);
        }

        return r;
    },

    /**
     * Inverting an FP
     *
     * @this {FP}
     */
    inverse: function() {

        var m2=new BIG(0);
        m2.rcopy(ROM_FIELD.Modulus);
        m2.dec(2); m2.norm();
        this.copy(this.pow(m2));
        return this;
    },

    /**
     * Tests for equality of two FP instances
     *
     * @this {FP}
     * @param x FP instance to compare
     */
    equals: function(a) {
        var ft=new FP(0); ft.copy(this);
        var sd=new FP(0); sd.copy(a);
        ft.reduce();
        sd.reduce();

        if (BIG.comp(ft.f, sd.f) === 0) {
            return true;
        }

        return false;
    },

    /**
     * Raises an FP to the power of a BIG
     *
     * @this {FP}
     * @param e BIG instance exponent
     */
    pow: function(e) {
        var i,w=[],
            tb=[],
            t=new BIG(e),
            nb, lsbs, r;
        this.norm();
        t.norm();
        nb= 1 + Math.floor((t.nbits() + 3) / 4);

        for (i=0;i<nb;i++) {
            lsbs=t.lastbits(4);
            t.dec(lsbs);
            t.norm();
            w[i]=lsbs;
            t.fshr(4);
        }
        tb[0]=new FP(1);
        tb[1]=new FP(this);
        for (i=2;i<16;i++) {
            tb[i]=new FP(tb[i-1]);
            tb[i].mul(this);
        }
        r=new FP(tb[w[nb-1]]);
        for (i=nb-2;i>=0;i--) {
            r.sqr();
            r.sqr();
            r.sqr();
            r.sqr();
            r.mul(tb[w[i]]);
        }
        r.reduce();
        return r;
    },

    /**
     * return jacobi symbol (this/Modulus)
     *
     * @this {FP}
     */
    jacobi: function() {
        var p = new BIG(0),
            w = this.redc();

        p.rcopy(ROM_FIELD.Modulus);

        return w.jacobi(p);
    },

    /**
     * Fast Modular square root of a an FP, mod Modulus
     *
     * @this {FP}
     */
    sqrt: function() {
        this.reduce();
        var b = new BIG(0);
        b.rcopy(ROM_FIELD.Modulus);
        b.inc(1);
        b.norm();
        b.shr(2);
        return this.pow(b);
    }
};

FP.logb2 = function(v) {
    var r;

    v |= v >>> 1;
    v |= v >>> 2;
    v |= v >>> 4;
    v |= v >>> 8;
    v |= v >>> 16;

    v = v - ((v >>> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
    r = ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;

    return r;
};

FP.quo = function(n,m) {
    var num,den,hb=BIG.CHUNK>>1;
    if (FP.TBITS < hb) {
        var sh=hb-FP.TBITS;
        num=(n.w[BIG.NLEN-1]<<sh)|(n.w[BIG.NLEN-2]>>(BIG.BASEBITS-sh));
        den=(m.w[BIG.NLEN-1]<<sh)|(m.w[BIG.NLEN-2]>>(BIG.BASEBITS-sh));
    } else {
        num=n.w[BIG.NLEN-1];
        den=m.w[BIG.NLEN-1];
    }
    return Math.floor(num/(den+1));
};

/**
  * reduce a DBIG to a BIG using a "special" modulus
  *
  * @this {FP}
  */
FP.mod = function(d) {
    var b = new BIG(0),
        i, t, v, tw, tt, lo, carry, m, dd;

    m = new BIG(0);
    m.rcopy(ROM_FIELD.Modulus);

    b.copy(BIG.monty(m, ROM_FIELD.MConst, d));

    return b;
};

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = FP;
}
