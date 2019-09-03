import {FP12} from "./fp12";
import {ECP} from "./ecp";
import {ECP2} from "./ecp2";

export class PAIR {
    static G1mul(ECP, BIG): ECP;
    static G2mul(ECP2, BIG): ECP2;
    static ate(p1: ECP2, q1: ECP): FP12;
    static fexp(p: FP12): FP12;
}
