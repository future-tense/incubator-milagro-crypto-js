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

"use strict";

/* AMCL BIG number class */

/**
 * General purpose Constructor
 *
 * @constructor
 * @this {BIG}
 */

const BIG = function(x) {
    this.w = new Array(BIG.NLEN);

    switch (typeof(x)) {
        case "object":
            this.copy(x);
            break;

        case "number":
            this.zero();
            this.w[0] = x;
            break;

        default:
            this.zero();
    }
};

BIG.CHUNK = 32;
BIG.MODBYTES = 32;
BIG.BASEBITS = 24;
BIG.NLEN = (1 + (Math.floor((8 * BIG.MODBYTES - 1) / BIG.BASEBITS)));
BIG.DNLEN = 2 * BIG.NLEN;
BIG.BMASK = (1 << BIG.BASEBITS) - 1;
BIG.BIGBITS = (8 * BIG.MODBYTES);
BIG.NEXCESS = (1 << (BIG.CHUNK - BIG.BASEBITS - 1));
BIG.MODINV = (Math.pow(2, -BIG.BASEBITS));

BIG.prototype = {

    /**
     * set to zero
     *
     * @this {BIG}
     * @return BIG number
     */
    zero: function() {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    /**
     * set to one
     *
     * @this {BIG}
     * @return BIG number
     */
    one: function() {

        this.w[0] = 1;
        for (let i = 1; i < BIG.NLEN; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    get: function(i) {
        return this.w[i];
    },

    set: function(i, x) {
        this.w[i] = x;
    },

    /**
     * test for zero
     *
     * @this {BIG}
     * @return boolean - True if zero
     */
    iszilch: function() {
        for (let i = 0; i < BIG.NLEN; i++) {
            if (this.w[i] !== 0) {
                return false;
            }
        }

        return true;
    },

    /**
     * test for unity
     *
     * @this {BIG}
     * @return boolean - True if one
     */
    isunity: function() {
        for (let i = 1; i < BIG.NLEN; i++) {
            if (this.w[i] !== 0) {
                return false;
            }
        }

        return (this.w[0] === 1);
    },

    /**
     * Conditional swap of two BIGs depending on d using XOR - no branches
     *
     * @this {BIG}
     * @parameter b BIG number
     * @parameter d BIG number
     */
    cswap: function(b, d) {
        const c = ~(d - 1);
        for (let i = 0; i < BIG.NLEN; i++) {
            const t = c & (this.w[i] ^ b.w[i]);
            this.w[i] ^= t;
            b.w[i] ^= t;
        }
    },

    /**
     * Conditional move of BIG depending on d using XOR - no branches
     *
     * @this {BIG}
     * @parameter b BIG number
     * @parameter d BIG number
     */
    cmove: function(b, d) {
        const c = ~(d - 1);
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] ^= (this.w[i] ^ b.w[i]) & c;
        }
    },

    /**
     * Copy from another BIG
     *
     * @this {BIG}
     * @parameter y BIG number
     * @return BIG - The BIG object
     */
    copy: function(y) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = y.w[i];
        }

        return this;
    },

    /**
     * copy from bottom half of DBIG
     *
     * @this {BIG}
     * @parameter y BIG number
     * @return BIG The new BIG object
     */
    hcopy: function(y) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = y.w[i];
        }

        return this;
    },

    /**
     *  copy from ROM
     *
     * @this {BIG}
     * @parameter y BIG number in ROM
     * @return BIG - The BIG object
     */
    rcopy: function(y) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = y[i];
        }

        return this;
    },

    xortop: function(x) {
        this.w[BIG.NLEN - 1] ^= x;
    },

    ortop: function(x) {
        this.w[BIG.NLEN - 1] |= x;
    },

    /**
     *  normalise BIG - force all digits < 2^BASEBITS
     *
     * @this {BIG}
     * @return BIG number
     */
    norm: function() {
        let carry = 0;
        for (let i = 0; i < BIG.NLEN - 1; i++) {
            const d = this.w[i] + carry;
            this.w[i] = d & BIG.BMASK;
            carry = d >> BIG.BASEBITS;
        }

        this.w[BIG.NLEN - 1] = (this.w[BIG.NLEN - 1] + carry);
//        return (this.w[BIG.NLEN - 1] >> ((8 * BIG.MODBYTES) % BIG.BASEBITS));
        return this;
    },

    /**
     * Quick Fast shifts a BIG right by a small number of bits - input must be normalised, output will be normalised
     *
     * @this {BIG}
     * @parameter k Number of bits to shift
     * @return r The shifted out part
     */
    fshr: function(k) {

        const r = this.w[0] & ((1 << k) - 1); /* shifted out part */

        for (let i = 0; i < BIG.NLEN - 1; i++) {
            this.w[i] = (this.w[i] >> k) | ((this.w[i + 1] << (BIG.BASEBITS - k)) & BIG.BMASK);
        }

        this.w[BIG.NLEN - 1] = this.w[BIG.NLEN - 1] >> k;

        return r;
    },

    /**
     * General shift right by k bits
     *
     * @this {BIG}
     * @parameter k Number of bits to shift
     * @return BIG number
     */
    shr: function(k) {
        const n = k % BIG.BASEBITS;
        const m = Math.floor(k / BIG.BASEBITS);

        for (let i = 0; i < BIG.NLEN - m - 1; i++) {
            this.w[i] = (this.w[m + i] >> n) | ((this.w[m + i + 1] << (BIG.BASEBITS - n)) & BIG.BMASK);
        }

        this.w[BIG.NLEN - m - 1] = this.w[BIG.NLEN - 1] >> n;

        for (let i = BIG.NLEN - m; i < BIG.NLEN; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    /**
     * Fast shifts a BIG left by a small number of bits - input must be normalised, output will be normalised
     *
     * @this {BIG}
     * @parameter k Number of bits to shift
     * @return r The shifted out part
     */
    fshl: function(k) {

        this.w[BIG.NLEN - 1] = ((this.w[BIG.NLEN - 1] << k)) | (this.w[BIG.NLEN - 2] >> (BIG.BASEBITS - k));

        for (let i = BIG.NLEN - 2; i > 0; i--) {
            this.w[i] = ((this.w[i] << k) & BIG.BMASK) | (this.w[i - 1] >> (BIG.BASEBITS - k));
        }

        this.w[0] = (this.w[0] << k) & BIG.BMASK;

        return (this.w[BIG.NLEN - 1] >> ((8 * BIG.MODBYTES) % BIG.BASEBITS)); /* return excess - only used in FF.java */
    },

    /**
     * General shift left by k bits
     *
     * @this {BIG}
     * @parameter k Number of bits to shift
     * @return BIG number
     */
    shl: function(k) {
        const n = k % BIG.BASEBITS;
        const m = Math.floor(k / BIG.BASEBITS);

        this.w[BIG.NLEN - 1] = (this.w[BIG.NLEN - 1 - m] << n);

        if (BIG.NLEN > m + 2) {
            this.w[BIG.NLEN - 1] |= (this.w[BIG.NLEN - m - 2] >> (BIG.BASEBITS - n));
        }

        for (let i = BIG.NLEN - 2; i > m; i--) {
            this.w[i] = ((this.w[i - m] << n) & BIG.BMASK) | (this.w[i - m - 1] >> (BIG.BASEBITS - n));
        }

        this.w[m] = (this.w[0] << n) & BIG.BMASK;

        for (let i = 0; i < m; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    /**
     * length in bits
     *
     * @this {BIG}
     * @return The number of bigs in BIG object
     */
    nbits: function() {

        const t = new BIG(0);
        t.copy(this);
        t.norm();

        let k = BIG.NLEN - 1;
        while (k >= 0 && t.w[k] === 0) {
            k--;
        }

        if (k < 0) {
            return 0;
        }

        let bts = BIG.BASEBITS * k;
        let c = t.w[k];

        while (c !== 0) {
            c = Math.floor(c / 2);
            bts++;
        }

        return bts;
    },

    /**
     * Convert to string
     *
     * @this {BIG}
     * @return string representation of a BIG number
     */
    toString: function() {

        let len = this.nbits();
        if (len % 4 === 0) {
            len = Math.floor(len / 4);
        } else {
            len = Math.floor(len / 4);
            len++;
        }

        if (len < BIG.MODBYTES * 2) {
            len = BIG.MODBYTES * 2;
        }

        let s = "";
        for (let i = len - 1; i >= 0; i--) {
            const b = new BIG(0);
            b.copy(this);
            b.shr(i * 4);
            s += (b.w[0] & 15).toString(16);
        }

        return s;
    },

    /**
     * Sum two BIG mumbers
     *
     * @this {BIG}
     * @parameter y BIG object
     * @return BIG - this+=y
     */
    add: function(y) {

        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] += y.w[i];
        }

        return this;
    },


    /**
     * OR two BIG mumbers
     *
     * @this {BIG}
     * @parameter y BIG object
     * @return BIG - this|=y
     */
    or: function(y) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] |= y.w[i];
        }

        return this;
    },

    /**
     * Sum two BIG mumbers
     *
     * @this {BIG}
     * @parameter x BIG object
     * @return BIG this+x
     */
    plus: function(x) {
        const s = new BIG(0);
        for (let i = 0; i < BIG.NLEN; i++) {
            s.w[i] = this.w[i] + x.w[i];
        }

        return s;
    },

    /**
     * Sum BIG and integer
     *
     * @this {BIG}
     * @parameter i Integer to add
     * @return BIG this+=i
     */
    inc: function(i) {
        this.norm();
        this.w[0] += i;
        return this;
    },

    /**
     * Subtract BIG from one another
     *
     * @this {BIG}
     * @parameter y BIG object
     * @return BIG this-=y
     */
    sub: function(y) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] -= y.w[i];
        }

        return this;
    },

    /**
     * Reverse subtract BIG from one another
     *
     * @this {BIG}
     * @parameter x BIG object
     * @return BIG this=x-this
     */
    rsub: function(x) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = x.w[i] - this.w[i];
        }

        return this;
    },

    /**
     * Subtract integer from BIG
     *
     * @this {BIG}
     * @parameter i Integer to subtract
     * @return BIG this-=i
     */
    dec: function(i) {
        this.norm();
        this.w[0] -= i;
        return this;
    },

    /**
     * Subtract BIG
     *
     * @this {BIG}
     * @parameter x BIG object
     * @return New BIG object
     */
    minus: function(x) {
        const d = new BIG(0);
        for (let i = 0; i < BIG.NLEN; i++) {
            d.w[i] = this.w[i] - x.w[i];
        }

        return d;
    },

    /**
     * Multiply by small integer
     *
     * @this {BIG}
     * @parameter c small integer
     * @return this*c
     */
    imul: function(c) {
        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] *= c;
        }

        return this;
    },

    /**
     * convert this BIG to byte array
     *
     * @this {BIG}
     */
    tobytearray: function(b, n) {
        const c = new BIG(0);
        c.copy(this);
        c.norm();

        for (let i = BIG.MODBYTES - 1; i >= 0; i--) {
            b[i + n] = c.w[0] & 0xff;
            c.fshr(8);
        }

        return this;
    },

    /**
     * convert this to byte array
     *
     * @this {BIG}
     */
    toBytes: function(b) {
        this.tobytearray(b, 0);
    },

    /**
     * this[i]+=x*y+c, and return high part
     *
     * @this {BIG}
     */
    muladd: function(x, y, c, i) {
        const prod = x * y + c + this.w[i];
        this.w[i] = prod & BIG.BMASK;
        return ((prod - this.w[i]) * BIG.MODINV);
    },

    /**
     *  multiply by larger int
     *
     * @this {BIG}
     * @parameter c large integer
     * @return number - carry value
     */
    pmul: function(c) {
        let carry = 0;
        for (let i = 0; i < BIG.NLEN; i++) {
            const ak = this.w[i];
            this.w[i] = 0;
            carry = this.muladd(ak, c, carry, i);
        }

        return carry;
    },

    /**
     *  multiply by still larger int - results requires a DBIG
     *
     * @this {BIG}
     * @parameter c large integer
     * @return DBIG object
     */
    pxmul: function(c) {
        const m = new DBIG(0);
        let carry = 0;
        for (let j = 0; j < BIG.NLEN; j++) {
            carry = m.muladd(this.w[j], c, carry, j);
        }

        m.w[BIG.NLEN] = carry;

        return m;
    },

    /**
     *  divide by 3
     *
     * @this {BIG}
     * @return number - carry value
     */
    div3: function() {
        this.norm();
        const base = (1 << BIG.BASEBITS);

        let carry = 0;
        for (let i = BIG.NLEN - 1; i >= 0; i--) {
            const ak = (carry * base + this.w[i]);
            this.w[i] = Math.floor(ak / 3);
            carry = ak % 3;
        }
        return carry;
    },

    /**
     * set x = x mod 2^m
     *
     * @this {BIG}
     * @parameter m Exponent
     * @return BIG object
     */
    mod2m: function(m) {

        const wd = Math.floor(m / BIG.BASEBITS);
        const bt = m % BIG.BASEBITS;
        const msk = (1 << bt) - 1;
        this.w[wd] &= msk;

        for (let i = wd + 1; i < BIG.NLEN; i++) {
            this.w[i] = 0;
        }
    },

    /**
     * a=1/a mod 2^256. This is very fast!
     *
     * @this {BIG}
     * @return BIG object
     */
    invmod2m: function() {
        const U = new BIG(0);
        let b = new BIG(0);
        const c = new BIG(0);

        U.inc(BIG.invmod256(this.lastbits(8)));

        for (let i = 8; i < BIG.BIGBITS; i <<= 1) {
            U.norm();
            b.copy(this);
            b.mod2m(i);
            const t1 = BIG.smul(U, b);
            t1.shr(i);
            c.copy(this);
            c.shr(i);
            c.mod2m(i);

            const t2 = BIG.smul(U, c);
            t2.mod2m(i);
            t1.add(t2);
            t1.norm();
            b = BIG.smul(t1, U);
            t1.copy(b);
            t1.mod2m(i);

            t2.one();
            t2.shl(i);
            t1.rsub(t2);
            t1.norm();
            t1.shl(i);
            U.add(t1);
        }

        U.mod2m(BIG.BIGBITS);
        this.copy(U);
        this.norm();
    },

    /**
     * reduce this mod m
     *
     * @this {BIG}
     */
    mod: function(m1) {
        this.norm();

        const m = new BIG(0);
        m.copy(m1);

        if (BIG.comp(this, m) < 0) {
            return this;
        }


        let k = 0;
        do {
            m.fshl(1);
            k++;
        } while (BIG.comp(this, m) >= 0);

        const r = new BIG(0);
        while (k > 0) {
            m.fshr(1);

            r.copy(this);
            r.sub(m);
            r.norm();
            this.cmove(r, (1 - ((r.w[BIG.NLEN - 1] >> (BIG.CHUNK - 1)) & 1)));
            k--;
        }
        return this;
    },

    /**
     * this/=m1
     *
     * @this {BIG}
     * @paramter m1 divisor
     * @return BIG number
     */
    div: function(m1) {

        const m = new BIG(0);
        m.copy(m1);

        this.norm();
        const b = new BIG(0);
        b.copy(this);
        this.zero();

        const e = new BIG(1);
        let k = 0;
        while (BIG.comp(b, m) >= 0) {
            e.fshl(1);
            m.fshl(1);
            k++;
        }

        const r = new BIG(0);
        while (k > 0) {
            m.fshr(1);
            e.fshr(1);

            r.copy(b);
            r.sub(m);
            r.norm();
            const d = (1 - ((r.w[BIG.NLEN - 1] >> (BIG.CHUNK - 1)) & 1));
            b.cmove(r, d);
            r.copy(this);
            r.add(e);
            r.norm();
            this.cmove(r, d);

            k--;
        }
    },

    /**
     * return parity of this
     *
     * @this {BIG}
     * @return number
     */
    parity: function() {
        return this.w[0] & 1;
    },

    /**
     * return n-th bit of this
     *
     * @this {BIG}
     * @parameter nth bit to return
     * @return number - bit value
     */
    bit: function(n) {
        if ((this.w[Math.floor(n / BIG.BASEBITS)] & (1 << (n % BIG.BASEBITS))) > 0) {
            return 1;
        } else {
            return 0;
        }
    },

    /**
     * return last n bits of this
     *
     * @this {BIG}
     * @parameter n bits to return
     * @return number - bit values
     */
    lastbits: function(n) {
        const msk = (1 << n) - 1;
        this.norm();
        return (this.w[0]) & msk;
    },

    isok: function() {

        for (let i = 0; i < BIG.NLEN; i++) {
            if ((this.w[i] >> BIG.BASEBITS) !== 0) {
                return false;
            }
        }

        return true;
    },

    /**
     * Jacobi Symbol (this/p)
     *
     * @this {BIG}
     * @parameter p BIG number
     * @return 0, 1 or -1
     */
    jacobi: function(p) {
        const zilch = new BIG(0);
        const one = new BIG(1);

        if (p.parity() === 0 || BIG.comp(this, zilch) === 0 || BIG.comp(p, one) <= 0) {
            return 0;
        }

        this.norm();
        const x = new BIG(0);
        x.copy(this);
        x.mod(p);

        const n = new BIG(0);
        n.copy(p);

        const t = new BIG(0);
        let m = 0;
        while (BIG.comp(n, one) > 0) {
            if (BIG.comp(x, zilch) === 0) {
                return 0;
            }

            const n8 = n.lastbits(3);
            let k = 0;

            while (x.parity() === 0) {
                k++;
                x.shr(1);
            }

            if (k % 2 === 1) {
                m += (n8 * n8 - 1) / 8;
            }

            m += (n8 - 1) * (x.lastbits(2) - 1) / 4;
            t.copy(n);
            t.mod(x);
            n.copy(x);
            x.copy(t);
            m &= 1;
        }

        if (m === 0) {
            return 1;
        } else {
            return -1;
        }
    },

    /**
     * this=1/this mod p. Binary method
     *
     * @this {BIG}
     * @parameter p The BIG Modulus
     * @return BIG object
     */
    invmodp: function(p) {
        this.mod(p);
        const u = new BIG(0);
        u.copy(this);

        const v = new BIG(0);
        v.copy(p);

        const x1 = new BIG(1);
        const x2 = new BIG(0);
        const t = new BIG(0);
        const one = new BIG(1);
        while (BIG.comp(u, one) !== 0 && BIG.comp(v, one) !== 0) {
            while (u.parity() === 0) {
                u.fshr(1);
                if (x1.parity() !== 0) {
                    x1.add(p);
                    x1.norm();
                }
                x1.fshr(1);
            }

            while (v.parity() === 0) {
                v.fshr(1);
                if (x2.parity() !== 0) {
                    x2.add(p);
                    x2.norm();
                }
                x2.fshr(1);
            }

            if (BIG.comp(u, v) >= 0) {
                u.sub(v);
                u.norm();
                if (BIG.comp(x1, x2) >= 0) {
                    x1.sub(x2);
                } else {
                    t.copy(p);
                    t.sub(x2);
                    x1.add(t);
                }
                x1.norm();
            } else {
                v.sub(u);
                v.norm();
                if (BIG.comp(x2, x1) >= 0) {
                    x2.sub(x1);
                } else {
                    t.copy(p);
                    t.sub(x1);
                    x2.add(t);
                }
                x2.norm();
            }
        }

        if (BIG.comp(u, one) === 0) {
            this.copy(x1);
        } else {
            this.copy(x2);
        }
    },

    /**
     * Exponentation modulo m
     *
     * @this {BIG}
     * @parameter e1 BIG number
     * @parameter m  The BIG Modulus
     * @return BIG - this^e mod m
     */
    powmod: function(e1, m) {
        this.norm();

        const e = new BIG(0);
        e.copy(e1);
        e.norm();

        const z = new BIG(0);
        z.copy(e);

        let s = new BIG(0);
        s.copy(this);

        let a = new BIG(1);
        for (;;) {
            const bt = z.parity();
            z.fshr(1);
            if (bt === 1) {
                a = BIG.modmul(a, s, m);
            }

            if (z.iszilch()) {
                break;
            }

            s = BIG.modsqr(s, m);
        }

        return a;
    }
};

BIG.ssn = function(r, a, m) {
    const n = BIG.NLEN-1;

    m.w[0] = (m.w[0] >> 1) | ((m.w[1] << (BIG.BASEBITS - 1)) & BIG.BMASK);
    r.w[0] = a.w[0] - m.w[0];

    let carry = r.w[0] >> BIG.BASEBITS;
    r.w[0] &= BIG.BMASK;
    for (let i = 1; i < n; i++) {
        m.w[i] = (m.w[i] >> 1) | ((m.w[i + 1] << (BIG.BASEBITS - 1)) & BIG.BMASK);
        r.w[i] = a.w[i] - m.w[i] + carry;
        carry = r.w[i] >> BIG.BASEBITS;
        r.w[i] &= BIG.BMASK;
    }

    m.w[n] >>= 1;
    r.w[n] = a.w[n] - m.w[n] + carry;
    return ((r.w[n] >> (BIG.CHUNK - 1)) & 1);
};

/**
 * convert from byte array to BIG
 *
 * @this {BIG}
 * @parameter b Bytearray
 * @return BIG object
 */
BIG.frombytearray = function(b, n) {
    const m = new BIG(0);
    for (let i = 0; i < BIG.MODBYTES; i++) {
        m.fshl(8);
        m.w[0] += b[i + n] & 0xff;
    }

    return m;
};

BIG.fromBytes = function(b) {
    return BIG.frombytearray(b, 0);
};

/**
 * return a*b where product fits a BIG
 *
 * @this {BIG}
 * @parameter a BIG number
 * @parameter b BIG number
 * @return BIG - a*b
 */
BIG.smul = function(a, b) {
    const c = new BIG(0);
    for (let i = 0; i < BIG.NLEN; i++) {
        let carry = 0;

        for (let j = 0; j < BIG.NLEN; j++) {
            if (i + j < BIG.NLEN) {
                carry = c.muladd(a.w[i], b.w[j], carry, i + j);
            }
        }
    }

    return c;
};

/**
 * Compare a and b
 *
 * @this {BIG}
 * @parameter a BIG number (normalised)
 * @parameter b BIG number (normalised
 * @return number - 0 if a==b, -1 if a<b, +1 if a>b
 */
BIG.comp = function(a, b) {

    for (let i = BIG.NLEN - 1; i >= 0; i--) {
        if (a.w[i] === b.w[i]) {
            continue;
        }

        if (a.w[i] > b.w[i]) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
};

/**
 * Multiple two BIG numbers
 *
 * @this {BIG}
 * @parameter a BIG number
 * @parameter b BIG number
 * @return DBIG - a*b as a DBIG number
 */
BIG.mul = function(a, b) {

    const d = new Array(BIG.NLEN);
    for (let i = 0; i < BIG.NLEN; i++) {
        d[i] = a.w[i] * b.w[i];
    }

    let s = d[0];
    let t = s;

    const c = new DBIG(0);
    c.w[0] = t;

    for (let k = 1; k < BIG.NLEN; k++) {
        s += d[k];
        t = s;
        for (let i = k; i >= 1 + Math.floor(k / 2); i--) {
            t += (a.w[i] - a.w[k - i]) * (b.w[k - i] - b.w[i]);
        }
        c.w[k] = t;
    }
    for (let k = BIG.NLEN; k < 2 * BIG.NLEN - 1; k++) {
        s -= d[k - BIG.NLEN];
        t = s;
        for (let i = BIG.NLEN - 1; i >= 1 + Math.floor(k / 2); i--) {
            t += (a.w[i] - a.w[k - i]) * (b.w[k - i] - b.w[i]);
        }
        c.w[k] = t;
    }

    let co = 0;
    for (let i = 0; i < BIG.DNLEN - 1; i++) {
        const n = c.w[i] + co;
        c.w[i] = n & BIG.BMASK;
        co = (n - c.w[i]) * BIG.MODINV;
    }
    c.w[BIG.DNLEN - 1] = co;

    return c;
};

/**
 * Square two BIG numbers
 *
 * @this {BIG}
 * @parameter a BIG number
 * @return a*2 as a DBIG number
 */
BIG.sqr = function(a) {

    const c = new DBIG(0);
    c.w[0] = a.w[0] * a.w[0];

    for (let j = 1; j < BIG.NLEN - 1;) {
        let t = a.w[j] * a.w[0];
        for (let i = 1; i < (j + 1) >> 1; i++) {
            t += a.w[j - i] * a.w[i];
        }
        t += t;
        c.w[j] = t;
        j++;
        t = a.w[j] * a.w[0];
        for (let i = 1; i < (j + 1) >> 1; i++) {
            t += a.w[j - i] * a.w[i];
        }
        t += t;
        t += a.w[j >> 1] * a.w[j >> 1];
        c.w[j] = t;
        j++;
    }

    for (let j = BIG.NLEN - 1 + BIG.NLEN % 2; j < BIG.DNLEN - 3;) {
        let t = a.w[BIG.NLEN - 1] * a.w[j - BIG.NLEN + 1];
        for (let i = j - BIG.NLEN + 2; i < (j + 1) >> 1; i++) {
            t += a.w[j - i] * a.w[i];
        }
        t += t;
        c.w[j] = t;
        j++;
        t = a.w[BIG.NLEN - 1] * a.w[j - BIG.NLEN + 1];
        for (let i = j - BIG.NLEN + 2; i < (j + 1) >> 1; i++) {
            t += a.w[j - i] * a.w[i];
        }
        t += t;
        t += a.w[j >> 1] * a.w[j >> 1];
        c.w[j] = t;
        j++;
    }

    let t = a.w[BIG.NLEN - 2] * a.w[BIG.NLEN - 1];
    t += t;
    c.w[BIG.DNLEN - 3] = t;

    t = a.w[BIG.NLEN - 1] * a.w[BIG.NLEN - 1];
    c.w[BIG.DNLEN - 2] = t;

    let co = 0;
    for (let i = 0; i < BIG.DNLEN - 1; i++) {
        const n = c.w[i] + co;
        c.w[i] = n & BIG.BMASK;
        co = (n - c.w[i]) * BIG.MODINV;
    }
    c.w[BIG.DNLEN - 1] = co;

    return c;
};

BIG.monty = function(m, nd, d) {

    let t = d.w[0];
    const v = new Array(BIG.NLEN);
    v[0] = ((t & BIG.BMASK) * nd) & BIG.BMASK;
    t += v[0] * m.w[0];

    const dd = new Array(BIG.NLEN);
    let c = d.w[1] + (t * BIG.MODINV);
    let s = 0;

    for (let k = 1; k < BIG.NLEN; k++) {
        let t = c + s + v[0] * m.w[k];
        for (let i = k - 1; i > Math.floor(k / 2); i--) {
            t += (v[k - i] - v[i]) * (m.w[i] - m.w[k - i]);
        }
        v[k] = ((t & BIG.BMASK) * nd) & BIG.BMASK;
        t += v[k] * m.w[0];
        c = (t * BIG.MODINV) + d.w[k + 1];

        dd[k] = v[k] * m.w[k];
        s += dd[k];
    }

    const b = new BIG(0);
    for (let k = BIG.NLEN; k < 2 * BIG.NLEN - 1; k++) {
        let t = c + s;
        for (let i = BIG.NLEN - 1; i >= 1 + Math.floor(k / 2); i--) {
            t += (v[k - i] - v[i]) * (m.w[i] - m.w[k - i]);
        }
        b.w[k - BIG.NLEN] = t & BIG.BMASK;
        c = ((t - b.w[k - BIG.NLEN]) * BIG.MODINV) + d.w[k + 1];

        s -= dd[k - BIG.NLEN + 1];
    }

    b.w[BIG.NLEN - 1] = c & BIG.BMASK;

    return b;
};

/**
 * Multiple two BIG numbers modulo m
 *
 * @this {BIG}
 * @parameter a1 BIG number
 * @parameter b1 BIG number
 * @parameter m The BIG Modulus
 * @return BIG - a1*b1 mod m  as a BIG number
 */
BIG.modmul = function(a1, b1, m) {
    const a = new BIG(0);
    a.copy(a1);
    a.mod(m);

    const b = new BIG(0);
    b.copy(b1);
    b.mod(m);

    const d = BIG.mul(a, b);
    return d.mod(m);
};

/**
 * Square a BIG number modulo m
 *
 * @this {BIG}
 * @parameter a1 BIG number
 * @parameter m The BIG Modulus
 * @return a*2 mod m  as a BIG number
 */
BIG.modsqr = function(a1, m) {
    const a = new BIG(0);
    a.copy(a1);
    a.mod(m);

    const d = BIG.sqr(a);
    return d.mod(m);
};

/**
 * Inversion
 *
 * @this {BIG}
 * @parameter a1 BIG number
 * @parameter m The BIG Modulus
 * @return -a1 mod m
 */
BIG.modneg = function(a1, m) {
    const a = new BIG(0);
    a.copy(a1);
    a.mod(m);

    return m.minus(a);
};

/**
 *  Arazi and Qi inversion mod 256
 *
 * @this {BIG}
 * @parameter a BIG number
 * @return BIG number
 */
BIG.invmod256 = function(a) {

    let t1 = 0;
    let c = (a >> 1) & 1;
    t1 += c;
    t1 &= 1;
    t1 = 2 - t1;
    t1 <<= 1;
    let U = t1 + 1;

    // i=2
    let b = a & 3;
    t1 = U * b;
    t1 >>= 2;
    c = (a >> 2) & 3;
    let t2 = (U * c) & 3;
    t1 += t2;
    t1 *= U;
    t1 &= 3;
    t1 = 4 - t1;
    t1 <<= 2;
    U += t1;

    // i=4
    b = a & 15;
    t1 = U * b;
    t1 >>= 4;
    c = (a >> 4) & 15;
    t2 = (U * c) & 15;
    t1 += t2;
    t1 *= U;
    t1 &= 15;
    t1 = 16 - t1;
    t1 <<= 4;
    U += t1;

    return U;
};

/* AMCL double length DBIG number class */

/**
 * General purpose Constructor
 *
 * @constructor
 * @this {DBIG}
 */
const DBIG = function(x) {
    this.w = new Array(BIG.DNLEN);
    this.zero();
    this.w[0] = x;
};

DBIG.prototype = {

    /**
     * set to zero
     *
     * @this {DBIG}
     * @return BIG number
     */
    zero: function() {
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] = 0;
        }
        return this;
    },

    /**
     * set to b
     *
     * @this {DBIG}
     * @parameter b DBIG number
     * @return DBIG number
     */
    copy: function(b) {
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] = b.w[i];
        }
        return this;
    },

    /**
     * copy from BIG
     *
     * @this {DBIG}
     * @parameter b BIG number
     * @return DBIG number
     */
    hcopy: function(b) {

        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = b.w[i];
        }

        for (let i = BIG.NLEN; i < BIG.DNLEN; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    ucopy: function(b) {

        for (let i = 0; i < BIG.NLEN; i++) {
            this.w[i] = 0;
        }

        for (let i = BIG.NLEN; i < BIG.DNLEN; i++) {
            this.w[i] = b.w[i - BIG.NLEN];
        }

        return this;
    },

    /**
     *  normalise DBIG - force all digits < 2^BASEBITS
     *
     * @this {DBIG}
     * @return DBIG number
     */
    norm: function() {

        let carry = 0;
        for (let i = 0; i < BIG.DNLEN - 1; i++) {
            const d = this.w[i] + carry;
            this.w[i] = d & BIG.BMASK;
            carry = d >> BIG.BASEBITS;
        }
        this.w[BIG.DNLEN - 1] = (this.w[BIG.DNLEN - 1] + carry);

        return this;
    },

    /**
     * this[i]+=x*y+c, and return high part
     *
     * @this {DBIG}
     */
    muladd: function(x, y, c, i) {
        const prod = x * y + c + this.w[i];
        this.w[i] = prod & BIG.BMASK;
        return ((prod - this.w[i]) * BIG.MODINV);
    },

    /**
     * General shift right by k bits
     *
     * @this {DBIG}
     * @parameter k Number of bits to shift
     * @return BIG number
     */
    shr: function(k) {
        const n = k % BIG.BASEBITS;
        const m = Math.floor(k / BIG.BASEBITS);

        for (let i = 0; i < BIG.DNLEN - m - 1; i++) {
            this.w[i] = (this.w[m + i] >> n) | ((this.w[m + i + 1] << (BIG.BASEBITS - n)) & BIG.BMASK);
        }

        this.w[BIG.DNLEN - m - 1] = this.w[BIG.DNLEN - 1] >> n;

        for (let i = BIG.DNLEN - m; i < BIG.DNLEN; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    /**
     * General shift left by k bits
     *
     * @this {BIG}
     * @parameter k Number of bits to shift
     * @return BIG number
     */
    shl: function(k) {
        const n = k % BIG.BASEBITS;
        const m = Math.floor(k / BIG.BASEBITS);

        this.w[BIG.DNLEN - 1] = ((this.w[BIG.DNLEN - 1 - m] << n)) | (this.w[BIG.DNLEN - m - 2] >> (BIG.BASEBITS - n));

        for (let i = BIG.DNLEN - 2; i > m; i--) {
            this.w[i] = ((this.w[i - m] << n) & BIG.BMASK) | (this.w[i - m - 1] >> (BIG.BASEBITS - n));
        }

        this.w[m] = (this.w[0] << n) & BIG.BMASK;

        for (let i = 0; i < m; i++) {
            this.w[i] = 0;
        }

        return this;
    },

    /**
     * Conditional move of BIG depending on d using XOR - no branches
     *
     * @this {DBIG}
     * @parameter b DBIG number
     * @parameter d DBIG number
     */
    cmove: function(b, d) {
        const c = ~(d - 1);
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] ^= (this.w[i] ^ b.w[i]) & c;
        }
    },

    /**
     * Sum two DBIG mumbers
     *
     * @this {DBIG}
     * @parameter x DBIG object
     * @return this+=x
     */
    add: function(x) {
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] += x.w[i];
        }
        return this;
    },

    /**
     * Subtract DBIG from one another
     *
     * @this {DBIG}
     * @parameter x BIG object
     * @return this-=x
     */
    sub: function(x) {
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] -= x.w[i];
        }
        return this;
    },

    rsub: function(x) {
        for (let i = 0; i < BIG.DNLEN; i++) {
            this.w[i] = x.w[i] - this.w[i];
        }
    },

    /**
     * length in bits
     *
     * @this {DBIG}
     * @return number - The number of bits in DBIG object
     */
    nbits: function() {

        const t = new DBIG(0);
        t.copy(this);
        t.norm();

        let k = BIG.DNLEN - 1;
        while (k >= 0 && t.w[k] === 0) {
            k--;
        }

        if (k < 0) {
            return 0;
        }

        let bts = BIG.BASEBITS * k;
        let c = t.w[k];

        while (c !== 0) {
            c = Math.floor(c / 2);
            bts++;
        }

        return bts;
    },

    /**
     * Convert to string
     *
     * @this {DBIG}
     * @return string representation of a BIG number
     */
    toString: function() {

        let len = this.nbits();
        if (len % 4 === 0) {
            len = Math.floor(len / 4);
        } else {
            len = Math.floor(len / 4);
            len++;
        }

        let s = "";
        for (let i = len - 1; i >= 0; i--) {
            const b = new DBIG(0);
            b.copy(this);
            b.shr(i * 4);
            s += (b.w[0] & 15).toString(16);
        }

        return s;
    },

    /**
     * reduces this DBIG mod a BIG, and returns the BIG
     *
     * @this {DBIG}
     * @return BIG object
     */
    mod: function(c) {
        this.norm();

        const m = new DBIG(0);
        m.hcopy(c);

        const r = new BIG(0);
        r.hcopy(this);

        if (DBIG.comp(this, m) < 0) {
            return r;
        }

        let k = 0;
        do {
            m.shl(1);
            k++;
        } while (DBIG.comp(this, m) >= 0);

        const dr = new DBIG(0);
        while (k > 0) {
            m.shr(1);

            dr.copy(this);
            dr.sub(m);
            dr.norm();
            this.cmove(dr, (1 - ((dr.w[BIG.DNLEN - 1] >> (BIG.CHUNK - 1)) & 1)));
            k--;
        }

        r.hcopy(this);

        return r;
    },

    /**
     * this/=c
     *
     * @this {DBIG}
     * @param c divisor
     * @return BIG number
     */
    div: function(c) {

        this.norm();

        const m = new DBIG(0);
        m.hcopy(c);

        const e = new BIG(1);
        let k = 0;
        while (DBIG.comp(this, m) >= 0) {
            e.fshl(1);
            m.shl(1);
            k++;
        }

        const a = new BIG(0);
        const r = new BIG(0);
        const dr = new DBIG(0);
        while (k > 0) {
            m.shr(1);
            e.shr(1);

            dr.copy(this);
            dr.sub(m);
            dr.norm();
            const d = (1 - ((dr.w[BIG.DNLEN - 1] >> (BIG.CHUNK - 1)) & 1));
            this.cmove(dr, d);
            r.copy(a);
            r.add(e);
            r.norm();
            a.cmove(r, d);

            k--;
        }

        return a;
    },

    /**
     * split this DBIG at position n, return higher half, keep lower half
     *
     * @this {DBIG}
     * @param n Position to splitdivisor
     * @return lower half BIG number
     */
    split: function(n) {
        const m = n % BIG.BASEBITS;
        const t = new BIG(0);
        let carry = this.w[BIG.DNLEN - 1] << (BIG.BASEBITS - m);
        for (let i = BIG.DNLEN - 2; i >= BIG.NLEN - 1; i--) {
            const nw = (this.w[i] >> m) | carry;
            carry = (this.w[i] << (BIG.BASEBITS - m)) & BIG.BMASK;
            t.w[i - BIG.NLEN + 1] = nw;
        }

        this.w[BIG.NLEN - 1] &= ((1 << m) - 1);

        return t;
    }

};

/**
 * Compare a and b
 *
 * @this {DBIG}
 * @parameter a DBIG number (normalised)
 * @parameter b DBIG number (normalised
 * @return number - 0 if a==b, -1 if a<b, +1 if a>b
 */
DBIG.comp = function(a, b) {

    for (let i = BIG.DNLEN - 1; i >= 0; i--) {
        if (a.w[i] === b.w[i]) {
            continue;
        }

        if (a.w[i] > b.w[i]) {
            return 1;
        } else {
            return -1;
        }
    }

    return 0;
};

if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = {
        BIG: BIG,
        DBIG: DBIG
    };
}
