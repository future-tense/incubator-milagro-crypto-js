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

var BIG = require("./big").BIG;
var FP = require("./fp");
var FP2 = require("./fp2");
var FP4 = require("./fp4");
var FP12 = require("./fp12");
var ECP = require("./ecp");
var ECP2 = require("./ecp2");
var ROM_FIELD = require("./rom_field");
var ROM_CURVE = require("./rom_curve");

"use strict";

/**
  * Creates an instance of PAIR
  *
  * @constructor
  * @this {PAIR}
  */
var PAIR = {

    /**
     * Line function
     *
     * @this {PAIR}
     */
    line: function(A, B, Qx, Qy) {
        var r = new FP12(1),
            c = new FP4(0),
            XX, YY, ZZ, YZ, sb,
            X1, Y1, T1, T2,
            a, b;

        if (A === B) { /* Doubling */
            XX = new FP2(A.getx());
            YY = new FP2(A.gety());
            ZZ = new FP2(A.getz());
            YZ = new FP2(YY);

            YZ.mul(ZZ); //YZ
            XX.sqr(); //X^2
            YY.sqr(); //Y^2
            ZZ.sqr(); //Z^2

            YZ.imul(4);
            YZ.neg();
            YZ.norm(); //-2YZ
            YZ.pmul(Qy); //-2YZ.Ys

            XX.imul(6); //3X^2
            XX.pmul(Qx); //3X^2.Xs

            sb = 3 * ROM_CURVE.CURVE_B_I;
            ZZ.imul(sb);
            ZZ.div_ip2();
            ZZ.norm(); // 3b.Z^2

            YY.add(YY);
            ZZ.sub(YY);
            ZZ.norm(); // 3b.Z^2-Y^2

            a = new FP4(YZ, ZZ); // -2YZ.Ys | 3b.Z^2-Y^2 | 3X^2.Xs
            b = new FP4(XX); // L(0,1) | L(0,0) | L(1,0)
            c = new FP4(0);
            A.dbl();
        } else { /* Addition */
            X1 = new FP2(A.getx()); // X1
            Y1 = new FP2(A.gety()); // Y1
            T1 = new FP2(A.getz()); // Z1
            T2 = new FP2(A.getz()); // Z1

            T1.mul(B.gety()); // T1=Z1.Y2
            T2.mul(B.getx()); // T2=Z1.X2

            X1.sub(T2);
            X1.norm(); // X1=X1-Z1.X2
            Y1.sub(T1);
            Y1.norm(); // Y1=Y1-Z1.Y2

            T1.copy(X1); // T1=X1-Z1.X2
            X1.pmul(Qy); // X1=(X1-Z1.X2).Ys
            T1.mul(B.gety()); // T1=(X1-Z1.X2).Y2

            T2.copy(Y1); // T2=Y1-Z1.Y2
            T2.mul(B.getx()); // T2=(Y1-Z1.Y2).X2
            T2.sub(T1);
            T2.norm(); // T2=(Y1-Z1.Y2).X2 - (X1-Z1.X2).Y2
            Y1.pmul(Qx);
            Y1.neg();
            Y1.norm(); // Y1=-(Y1-Z1.Y2).Xs

            a = new FP4(X1, T2); // (X1-Z1.X2).Ys  |  (Y1-Z1.Y2).X2 - (X1-Z1.X2).Y2  | - (Y1-Z1.Y2).Xs
            b = new FP4(Y1);
            c = new FP4(0);
            A.add(B);
        }

        r.set(a, b, c);
        r.settype(FP.SPARSER);

        return r;
    },

    /**
     * prepare for multi-pairing
     *
     * @this {PAIR}
     */
    initmp: function() {
        var r=[];
        for (var i=0;i<ECP.ATE_BITS;i++) {
            r[i] = new FP12(1);
        }
        return r;
    },

    /**
     * basic Miller loop
     *
     * @this {PAIR}
     * @param r FP12 precomputed array of accumulated line functions
     * @param res FP12 result
     */
    miller: function(r) {
        var res=new FP12(1);
        for (var i=ECP.ATE_BITS-1; i>=1; i--) {
            res.sqr();
            res.ssmul(r[i]);
        }

        res.conj();
        res.ssmul(r[0]);

        return res;
    },

    /**
     * Precompute line functions for n-pairing
     *
     * @this {PAIR}
     * @param r array of precomputed FP48 products of line functions
     * @param P1 An element of G2
     * @param Q1 An element of G1
     */
    another: function(r,P1,Q1) {

        var f;
        var n=new BIG(0);
        var n3=new BIG(0);
        var K=new ECP2();
        var lv,lv2;
        var bt;

        // P is needed in affine form for line function, Q for (Qx,Qy) extraction
        var P=new ECP2(); P.copy(P1); P.affine();
        var Q=new ECP(); Q.copy(Q1); Q.affine();

        P.affine();
        Q.affine();

        var fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        var fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);

        var Qx=new FP(Q.getx());
        var Qy=new FP(Q.gety());

        var A=new ECP2();
        A.copy(P);

        var MP=new ECP2();
        MP.copy(P); MP.neg();

        var nb=PAIR.lbits(n3,n);

        for (var i=nb-2;i>=1;i--) {
            lv=PAIR.line(A,A,Qx,Qy);

            bt=n3.bit(i)-n.bit(i);
            if (bt === 1) {
                lv2=PAIR.line(A,P,Qx,Qy);
                lv.smul(lv2);
            }
            if (bt === -1) {
                lv2=PAIR.line(A,MP,Qx,Qy);
                lv.smul(lv2);
            }
            r[i].ssmul(lv);
        }

        /* R-ate fixup required for BN curves */
        A.neg();

        K.copy(P);
        K.frob(f);
        lv=PAIR.line(A,K,Qx,Qy);
        K.frob(f);
        K.neg();
        lv2=PAIR.line(A,K,Qx,Qy);
        lv.smul(lv2);
        r[0].ssmul(lv);
    },

    /**
     * Calculate Miller loop for Optimal ATE pairing e(P,Q)
     *
     * @this {PAIR}
     * @param P1 An element of G2
     * @param Q1 An element of G1
     * @result r An element of GT i.e. result of the pairing calculation e(P,Q)
     */
    ate: function(P1, Q1) {
        var fa, fb, f, x, n, n3, K, lv, lv2,
            Qx, Qy, A, NP, r, nb, bt,
            i;

        n = new BIG(0);
        n3 = new BIG(0);
        K = new ECP2();

        fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);
        var P=new ECP2(); P.copy(P1); P.affine();
        var Q=new ECP(); Q.copy(Q1); Q.affine();

        Qx = new FP(Q.getx());
        Qy = new FP(Q.gety());

        A = new ECP2();
        r = new FP12(1);
        A.copy(P);

        NP = new ECP2();
        NP.copy(P);
        NP.neg();

        nb = PAIR.lbits(n3,n);

        for (i = nb - 2; i >= 1; i--) {
            r.sqr();
            lv = PAIR.line(A, A, Qx, Qy);
            bt=n3.bit(i)-n.bit(i);

            if (bt === 1) {
                lv2 = PAIR.line(A, P, Qx, Qy);
                lv.smul(lv2);
            }
            if (bt === -1) {
                lv2 = PAIR.line(A, NP, Qx, Qy);
                lv.smul(lv2);
            }
            r.ssmul(lv);
        }

        r.conj();

        /* R-ate fixup */
        A.neg();

        K.copy(P);
        K.frob(f);

        lv = PAIR.line(A, K, Qx, Qy);
        K.frob(f);
        K.neg();
        lv2 = PAIR.line(A, K, Qx, Qy);
        lv.smul(lv2);
        r.ssmul(lv);

        return r;
    },

    /**
     * Calculate Miller loop for Optimal ATE double-pairing e(P,Q).e(R,S)
     *
     * @this {PAIR}
     * @param P1 An element of G2
     * @param Q1 An element of G1
     * @param R1 An element of G2
     * @param S1 An element of G1
     * @result r An element of GT i.e. result of the double pairing calculation e(P,Q).e(R,S)
     */
    ate2: function(P1, Q1, R1, S1) {
        var fa, fb, f, x, n, n3, K, lv, lv2,
            Qx, Qy, Sx, Sy, A, B, NP,NR,r, nb, bt,
            i;

        n = new BIG(0);
        n3 = new BIG(0);
        K = new ECP2();

        fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);

        var P=new ECP2(); P.copy(P1); P.affine();
        var Q=new ECP(); Q.copy(Q1); Q.affine();
        var R=new ECP2(); R.copy(R1); R.affine();
        var S=new ECP(); S.copy(S1); S.affine();


        Qx = new FP(Q.getx());
        Qy = new FP(Q.gety());

        Sx = new FP(S.getx());
        Sy = new FP(S.gety());

        A = new ECP2();
        B = new ECP2();
        r = new FP12(1);

        A.copy(P);
        B.copy(R);

        NP = new ECP2();
        NP.copy(P);
        NP.neg();
        NR = new ECP2();
        NR.copy(R);
        NR.neg();

        nb = PAIR.lbits(n3,n);

        for (i = nb - 2; i >= 1; i--) {
            r.sqr();
            lv = PAIR.line(A, A, Qx, Qy);
            lv2 = PAIR.line(B, B, Sx, Sy);
            lv.smul(lv2);
            r.ssmul(lv);

            bt=n3.bit(i)-n.bit(i);

            if (bt === 1) {
                lv = PAIR.line(A, P, Qx, Qy);
                lv2 = PAIR.line(B, R, Sx, Sy);
                lv.smul(lv2);
                r.ssmul(lv);
            }
            if (bt === -1) {
                lv = PAIR.line(A, NP, Qx, Qy);
                lv2 = PAIR.line(B, NR, Sx, Sy);
                lv.smul(lv2);
                r.ssmul(lv);
            }
        }

        r.conj();

        // R-ate fixup required for BN curves
        A.neg();
        B.neg();
        K.copy(P);
        K.frob(f);

        lv = PAIR.line(A, K, Qx, Qy);
        K.frob(f);
        K.neg();
        lv2 = PAIR.line(A, K, Qx, Qy);
        lv.smul(lv2);
        r.ssmul(lv);

        K.copy(R);
        K.frob(f);

        lv = PAIR.line(B, K, Sx, Sy);
        K.frob(f);
        K.neg();
        lv2 = PAIR.line(B, K, Sx, Sy);
        lv.smul(lv2);
        r.ssmul(lv);

        return r;
    },

    /**
     * Final exponentiation of pairing, converts output of Miller loop to element in GT
     *
     * @this {PAIR}
     * @param m FP12 value
     * @result r m^((p^12-1)/r) where p is modulus and r is the group order
     */
    fexp: function(m) {
        var fa, fb, f, x, r, lv,
            x0, x1, x2, x3, x4, x5,
            y0, y1, y2, y3;

        fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);
        x = new BIG(0);
        x.rcopy(ROM_CURVE.CURVE_Bnx);

        r = new FP12(m);

        /* Easy part of final exp */
        lv = new FP12(r);
        lv.inverse();
        r.conj();
        r.mul(lv);
        lv.copy(r);
        r.frob(f);
        r.frob(f);
        r.mul(lv);

        /* Hard part of final exp */
        lv.copy(r);
        lv.frob(f);
        x0 = new FP12(lv); //x0.copy(lv);
        x0.frob(f);
        lv.mul(r);
        x0.mul(lv);
        x0.frob(f);
        x1 = new FP12(r); //x1.copy(r);
        x1.conj();

        x4 = r.pow(x);

        x3 = new FP12(x4); //x3.copy(x4);
        x3.frob(f);
        x2 = x4.pow(x);
        x5 = new FP12(x2); /*x5.copy(x2);*/
        x5.conj();
        lv = x2.pow(x);
        x2.frob(f);
        r.copy(x2);
        r.conj();

        x4.mul(r);
        x2.frob(f);

        r.copy(lv);
        r.frob(f);
        lv.mul(r);

        lv.usqr();
        lv.mul(x4);
        lv.mul(x5);
        r.copy(x3);
        r.mul(x5);
        r.mul(lv);
        lv.mul(x2);
        r.usqr();
        r.mul(lv);
        r.usqr();
        lv.copy(r);
        lv.mul(x1);
        r.mul(x0);
        lv.usqr();
        r.mul(lv);
        r.reduce();

        return r;
    }
};

/**
  * prepare ate parameter, n=6u+2 (BN) or n=u (BLS), n3=3*n
  *
  * @this {PAIR}
  */
PAIR.lbits = function(n3,n) {
    n.rcopy(ROM_CURVE.CURVE_Bnx);
    n.pmul(6);
    n.dec(2);

    n.norm();
    n3.copy(n);
    n3.pmul(3);
    n3.norm();
    return n3.nbits();
};

/**
  * GLV method
  *
  * @this {PAIR}
  */
PAIR.glv = function(e) {
    var u = [],
        t, q, v, d, x, x2, i, j;

    t = new BIG(0);
    q = new BIG(0);
    v = [];

    q.rcopy(ROM_CURVE.CURVE_Order);

    for (i = 0; i < 2; i++) {
        t.rcopy(ROM_CURVE.CURVE_W[i]);
        d = BIG.mul(t, e);
        v[i] = new BIG(d.div(q));
        u[i] = new BIG(0);
    }

    u[0].copy(e);

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            t.rcopy(ROM_CURVE.CURVE_SB[j][i]);
            t.copy(BIG.modmul(v[j], t, q));
            u[i].add(q);
            u[i].sub(t);
            u[i].mod(q);
        }
    }

    return u;
};

/**
  * Galbraith & Scott Method
  *
  * @this {PAIR}
  */
PAIR.gs = function(e) {
    var u = [],
        i, j, t, q, v, d, x, w;

    t = new BIG(0);
    q = new BIG(0);
    q.rcopy(ROM_CURVE.CURVE_Order);

    v = [];

    for (i = 0; i < 4; i++) {
        t.rcopy(ROM_CURVE.CURVE_WB[i]);
        d = BIG.mul(t, e);
        v[i] = new BIG(d.div(q));
        u[i] = new BIG(0);
    }

    u[0].copy(e);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            t.rcopy(ROM_CURVE.CURVE_BB[j][i]);
            t.copy(BIG.modmul(v[j], t, q));
            u[i].add(q);
            u[i].sub(t);
            u[i].mod(q);
        }
    }

    return u;
};

/**
  * Fast point multiplication of a member of the group G1 by a BIG number
  *
  * @this {PAIR}
  * @param P Member of G1
  * @param e BIG multiplier
  * @return R Member of G1 R=e.P
  */
PAIR.G1mul = function(P, e) {
    var R, Q, q, bcru, cru, t, u, np, nn;

    if (ROM_CURVE.USE_GLV) {
        R = new ECP();
        R.copy(P);
        Q = new ECP();
        Q.copy(P); Q.affine();
        q = new BIG(0);
        q.rcopy(ROM_CURVE.CURVE_Order);
        bcru = new BIG(0);
        bcru.rcopy(ROM_CURVE.CURVE_Cru);
        cru = new FP(bcru);
        t = new BIG(0);
        u = PAIR.glv(e);

        Q.getx().mul(cru);

        np = u[0].nbits();
        t.copy(BIG.modneg(u[0], q));
        nn = t.nbits();
        if (nn < np) {
            u[0].copy(t);
            R.neg();
        }

        np = u[1].nbits();
        t.copy(BIG.modneg(u[1], q));
        nn = t.nbits();
        if (nn < np) {
            u[1].copy(t);
            Q.neg();
        }
        u[0].norm();
        u[1].norm();
        R = R.mul2(u[0], Q, u[1]);
    } else {
        R = P.mul(e);
    }

    return R;
};

/**
  * Multiply P by e in group G2
  *
  * @this {PAIR}
  * @param P Member of G2
  * @param e BIG multiplier
  * @return R Member of G2 R=e.P
  */
PAIR.G2mul = function(P, e) {
    var R, Q, fa, fb, f, q, u, t, i, np, nn;

    if (ROM_CURVE.USE_GS_G2) {
        Q = [];
        fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);

        q = new BIG(0);
        q.rcopy(ROM_CURVE.CURVE_Order);

        u = PAIR.gs(e);
        t = new BIG(0);
        Q[0] = new ECP2();
        Q[0].copy(P);

        for (i = 1; i < 4; i++) {
            Q[i] = new ECP2();
            Q[i].copy(Q[i - 1]);
            Q[i].frob(f);
        }

        for (i = 0; i < 4; i++) {
            np = u[i].nbits();
            t.copy(BIG.modneg(u[i], q));
            nn = t.nbits();

            if (nn < np) {
                u[i].copy(t);
                Q[i].neg();
            }
            u[i].norm();
        }

        R = ECP2.mul4(Q, u);
    } else {
        R = P.mul(e);
    }
    return R;
};

/**
  * Fast raising of a member of GT to a BIG power
  *
  * @this {PAIR}
  * @param d Member of GT
  * @param e BIG exponent
  * @return r d^e
  */
PAIR.GTpow = function(d, e) {
    var r, g, fa, fb, f, q, t, u, i, np, nn;

    if (ROM_CURVE.USE_GS_GT) {
        g = [];
        fa = new BIG(0);
        fa.rcopy(ROM_FIELD.Fra);
        fb = new BIG(0);
        fb.rcopy(ROM_FIELD.Frb);
        f = new FP2(fa, fb);
        q = new BIG(0);
        q.rcopy(ROM_CURVE.CURVE_Order);
        t = new BIG(0);
        u = PAIR.gs(e);

        g[0] = new FP12(d);

        for (i = 1; i < 4; i++) {
            g[i] = new FP12(0);
            g[i].copy(g[i - 1]);
            g[i].frob(f);
        }

        for (i = 0; i < 4; i++) {
            np = u[i].nbits();
            t.copy(BIG.modneg(u[i], q));
            nn = t.nbits();

            if (nn < np) {
                u[i].copy(t);
                g[i].conj();
            }
            u[i].norm();
        }

        r = FP12.pow4(g, u);
    } else {
        r = d.pow(e);
    }

    return r;
};

// CommonJS module exports
if (typeof module !== "undefined" && typeof module.exports !== "undefined") {
    module.exports = PAIR;
}
